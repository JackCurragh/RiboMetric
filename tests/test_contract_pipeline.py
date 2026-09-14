"""End-to-end contract tests on the synthetic transcriptome in tests/synthetic.py.

These run the real command line on a BAM whose geometry is known exactly, so
each assertion is a consequence of a processing decision in
docs/METRIC_CONTRACT.md: coordinates, offsets, frame assignment, read sets and
weights.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

from RiboMetric.bam_processing import (
    _recover_alignment_tags_from_pysam,
    ox_parse_entire_bam,
    process_reads,
)

from . import synthetic

RUN = "import sys; from RiboMetric.cli import main; sys.exit(main())"


def _run(outdir: Path, *args: str) -> dict:
    cmd = [sys.executable, "-c", RUN, "run", *args, "--threads", "2", "--json"]
    cmd += ["-o", str(outdir), "-n", "synth"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    return json.loads((outdir / "synth_RiboMetric.json").read_text())["results"]


@pytest.fixture(scope="module")
def fixture_files(tmp_path_factory):
    return synthetic.build_fixture(tmp_path_factory.mktemp("synthetic"))


@pytest.fixture(scope="module")
def runs(fixture_files, tmp_path_factory):
    """Run each configuration once; tests read the cached results."""
    bam, annotation = fixture_files
    cache = {}

    def get(key, *args):
        if key not in cache:
            cache[key] = _run(tmp_path_factory.mktemp(key.replace(" ", "_")), *args)
        return cache[key]

    get.bam = str(bam)
    get.annotation = str(annotation)
    return get


def _annotated(runs, method="tripsviz", target="a_site", *extra):
    return runs(
        f"{method} {target} {' '.join(extra)}",
        "-b",
        runs.bam,
        "-a",
        runs.annotation,
        "--offset-calculation-method",
        method,
        "--offset-target",
        target,
        *extra,
    )


# --- parse: coordinates and orientation -----------------------------------------


def _write_bam(path: Path, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "TX1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, start, seq, cigar, tags in reads:
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 255
            a.cigarstring = cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            for tag, value in tags:
                a.set_tag(tag, value)
            out.write(a)
    pysam.index(str(path))


def test_read_start_is_0_based_and_includes_5prime_soft_clip(tmp_path):
    bam = tmp_path / "c.bam"
    body = "ACGTACGTACGTACGTACGTACGTACGTAC"
    _write_bam(
        bam,
        [
            ("r0_x3", 100, body, "30M", [("NH", 1)]),
            ("r1_x3", 200, "GG" + body, "2S30M", [("NH", 1)]),
        ],
    )
    reads, _ = ox_parse_entire_bam(str(bam))
    assert reads["reference_start"].tolist() == [100.0, 198.0]
    assert reads["read_length"].astype(int).tolist() == [30, 30]


def test_last_dinucleotide_is_read_5prime_to_3prime(tmp_path):
    bam = tmp_path / "d.bam"
    _write_bam(bam, [("r0_x1", 100, "ACGTACGTACGTACGTACGTACGTACGTAC", "30M", [("NH", 1)])])
    reads, _ = ox_parse_entire_bam(str(bam))
    assert reads["first_dinucleotide"].astype(str).tolist() == ["AC"]
    assert reads["last_dinucleotide"].astype(str).tolist() == ["AC"]


def test_tag_recovery_joins_on_read_name_and_keeps_string_xa(tmp_path):
    bam = tmp_path / "t.bam"
    seq = "ACGTACGTACGTACGTACGTACGTACGTAC"
    _write_bam(
        bam,
        [
            ("r1", 100, seq, "30M", []),
            ("r2", 110, seq, "30M", [("XA", "chr2,+5,30M,0;chr3,-9,30M,1;")]),
            ("r3", 120, seq, "30M", []),
        ],
    )
    # oxbow-shaped frame in a different order, MAPQ lost as null (STAR 255)
    oxbow_df = pd.DataFrame(
        {
            "qname": ["r3", "r1", "r2"],
            "seq": [seq] * 3,
            "cigar": ["30M"] * 3,
            "rname": ["TX1"] * 3,
            "pos": [121, 101, 111],
            "mapq": [np.nan] * 3,
        }
    )
    _recover_alignment_tags_from_pysam(oxbow_df, str(bam))
    assert oxbow_df["mapq"].tolist() == [255.0, 255.0, 255.0]
    reads = process_reads(oxbow_df)
    assert reads["xa"].tolist()[2] == 2.0  # two alternative loci for r2
    assert bool(reads["mapq_available"].all())


# --- the synthetic transcriptome -----------------------------------------------


@pytest.mark.parametrize("method", ["tripsviz", "changepoint", "ribowaltz"])
def test_offsets_match_known_geometry(runs, method):
    a = _annotated(runs, method, "a_site")
    p = _annotated(runs, method, "p_site")
    assert a["computed_offsets"]["28"] == 12 + 3
    assert a["computed_offsets"]["29"] == 13 + 3
    assert p["computed_offsets"]["28"] == 12
    assert p["computed_offsets"]["29"] == 13
    # TXC is uniform, so no length-30 peak falls in bounds: the A-site default
    assert a["computed_offsets"]["30"] == 15


def test_frame_table_excludes_multimappers_and_secondaries(runs):
    r = _annotated(runs)
    frames = {
        int(k): {int(f): n for f, n in v.items()} for k, v in r["read_frame_distribution"].items()
    }
    # 187 codons (3..189) inside the 9 nt exclusion, weight 5, from TXA and TXF
    assert frames[28] == {0: 2 * 187 * 5, 1: 0, 2: 0}
    assert frames[29] == {0: 187 * 5, 1: 0, 2: 0}
    # TXC: A-sites 70..650, one read each, frames (a - 60) mod 3
    assert frames[30] == {0: 193, 1: 194, 2: 194}


def test_periodicity_dominance_on_known_frames(runs):
    d = _annotated(runs)["metrics"]["periodicity_dominance"]
    assert d["28"] == 1.0
    assert d["29"] == 1.0
    assert d["30"] == pytest.approx(194 / 581)
    total = 1870 + 935 + 581
    assert d["global"] == pytest.approx((1870 + 935 + 193) / total)
    assert d["global_by_read_length_max"] == pytest.approx((1870 + 935 + 194) / total)


def test_read_lengths_and_weights_come_from_primary_records(runs):
    r = _annotated(runs)
    expected = synthetic.read_length_distribution()
    assert {int(k): v for k, v in r["read_length_distribution"].items()} == expected
    total = sum(expected.values())
    rec = r["recommended_read_lengths"]
    assert rec["recommended_lengths"] == [28, 29]
    assert rec["recommended_read_proportion"] == pytest.approx(
        round((expected[28] + expected[29]) / total, 4)
    )
    prims = synthetic.primary_records()
    multi = sum(p.weight for p in prims if p.nh > 1)
    m = r["metrics"]
    assert m["rpf_multimapper_rate"] == pytest.approx(multi / total)
    assert m["duplicate_rate"] == pytest.approx(1 - len(prims) / total)


def test_mode_and_no_cds_transcript(runs):
    r = _annotated(runs)
    assert r["mode"] == "annotation"
    # TXD contributes leader/trailer reads only
    assert r["metrics"]["prop_reads_leader"]["global"] > 0
    assert r["metrics"]["prop_reads_trailer"]["global"] > 0


def test_annotation_free_run_completes_without_frame_metrics(runs):
    r = runs("annotation-free", "-b", runs.bam)
    assert r["mode"] == "annotation_free"
    assert "periodicity_dominance" not in r["metrics"]
    assert r["read_frame_distribution"] == {}


def test_removing_sequences_changes_only_sequence_metrics(runs):
    full = _annotated(runs)
    skipped = _annotated(runs, "tripsviz", "a_site", "--skip-sequence-metrics")
    sequence_keys = {k for k in full["metrics"] if k.startswith("terminal_bias")}
    assert sequence_keys, "the full run should report terminal bias"
    assert not sequence_keys & set(skipped["metrics"])
    for key, value in full["metrics"].items():
        if key in sequence_keys:
            continue
        assert skipped["metrics"][key] == value, key


def test_start_and_stop_window_ratios(runs):
    m = _annotated(runs)["metrics"]
    # Unique annotated reads, weighted, by distance from cds_start: -5..20
    # against 30..50. TXA and TXB give 200 + 5 * 5 near the start and 7 * 5 in
    # the body each, TXF 6 * 5 and 7 * 5, TXC 26 and 21; TXE is not unique.
    assert m["start_codon_enrichment_ratio"] == pytest.approx(round(506 / 126, 4))
    # By distance from cds_end: 1..30 after against -30..-1 before. Only TXC
    # runs past the stop; TXA, TXB and TXF each put one weight-5 read at -30.
    assert m["stop_codon_readthrough_ratio"] == pytest.approx(round(30 / 45, 4))


def test_mapping_hygiene_on_known_records(runs):
    m = _annotated(runs)["metrics"]
    prims = synthetic.primary_records()
    multi_rows = sum(1 for p in prims if p.nh > 1)
    assert m["alignment_multimapper_rate"] == pytest.approx(multi_rows / len(prims))
    assert m["soft_clip_rate_5prime"] == 0.0
    assert m["disome_proportion"] == 0.0


# --- golden output --------------------------------------------------------------

GOLDEN = Path(__file__).parent / "golden" / "synthetic_tripsviz_a_site.json"


def _assert_close(actual, expected, path="results"):
    if isinstance(expected, dict):
        assert isinstance(actual, dict), path
        assert set(actual) == set(expected), f"{path}: keys differ {set(actual) ^ set(expected)}"
        for key in expected:
            _assert_close(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for i, (a, e) in enumerate(zip(actual, expected)):
            _assert_close(a, e, f"{path}[{i}]")
    elif isinstance(expected, float):
        assert actual == pytest.approx(expected, rel=1e-9, abs=1e-12), path
    else:
        assert actual == expected, path


def test_output_matches_golden_file(runs):
    """Any change to RiboMetric's output on the synthetic transcriptome fails
    here. If the change is intended, run ``python -m tests.golden.regenerate``
    and say in the commit message which numbers moved and why."""
    from .golden.regenerate import golden_view

    expected = json.loads(GOLDEN.read_text())
    actual = json.loads(json.dumps(golden_view(_annotated(runs))))
    _assert_close(actual, expected)
