"""Regression tests for the audit fixes.

Each test pins a specific bug that previously passed silently:
  * S1 — 5' and 3' terminal-bias backgrounds were identical.
  * S2 — GFF exon lengths were off by one (frame-scrambling, wrong tx length).
  * S6 — --offset-global was never applied.
  * U1 — evaluate gate ignored metric directionality.
  * E1 — bare `RiboMetric` crashed instead of printing help.
  * E2 — --output-offsets crashed because Path was shadowed inside main().
"""

import os

import pandas as pd

from RiboMetric.bam_processing import process_sequences
from RiboMetric.file_parser import deterministic_subsample, gff_df_to_cds_df, parse_gff
from RiboMetric.modules import a_site_calculation, read_frame_distribution_annotated
from RiboMetric.results_output import evaluate_qc_status

TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")


# --------------------------------------------------------------------------- #
# S1 — terminal-bias 5' vs 3' background must differ
# --------------------------------------------------------------------------- #
def test_five_and_three_prime_backgrounds_differ():
    seqs = ["AACCGGTT", "TTGGCCAA", "ACGTACGT", "TACGTACG"]
    out = process_sequences(seqs, [1, 1, 1, 1], pattern_length=2)
    diff = sum(abs(out["5_prime_bg"][k] - out["3_prime_bg"][k]) for k in out["5_prime_bg"])
    assert diff > 0, "5' and 3' backgrounds are identical (S1 regression)"


def test_three_prime_background_excludes_terminal_pattern():
    # Reads whose only variation is at the 3' end should produce a 3' background
    # that differs from the 5' background (the 3' terminal pattern is excluded).
    seqs = ["ACGTACGT", "ACGTACGA"]
    out = process_sequences(seqs, [1, 1], pattern_length=2)
    diff = sum(abs(out["5_prime_bg"][k] - out["3_prime_bg"][k]) for k in out["5_prime_bg"])
    assert diff > 0


def test_terminal_background_uses_count_weights():
    sequences = ["AAAAAA", "CCCCCC"]
    weighted = process_sequences(sequences, [9, 1], pattern_length=2)
    unweighted = process_sequences(sequences, [1, 1], pattern_length=2)
    assert weighted["5_prime_bg"]["AA"] > unweighted["5_prime_bg"]["AA"]


def test_subsample_is_seeded_and_keeps_repeated_names_together():
    reads = pd.DataFrame({"read_name": ["a", "a", "b", "c", "d"], "value": range(5)})
    first, info = deterministic_subsample(reads, 2, seed=7)
    repeat, repeat_info = deterministic_subsample(reads, 2, seed=7)
    assert first.equals(repeat)
    assert info == repeat_info
    assert info["realised_count"] == len(first)
    assert first["read_name"].value_counts().max() <= 2


# --------------------------------------------------------------------------- #
# S2 — GFF coordinates are 1-based inclusive
# --------------------------------------------------------------------------- #
def test_gff_cds_lengths_are_whole_codons():
    gff_df, _ = parse_gff(os.path.join(TEST_DATA, "1000_entry.gff"), 1000)
    res = gff_df_to_cds_df(gff_df)
    cds_len = res["cds_end"] - res["cds_start"]
    # Real CDS spans are whole codons; the previous off-by-one broke this.
    assert (cds_len % 3 == 0).all()
    assert (res["cds_start"] >= 0).all()
    assert (res["cds_end"] <= res["transcript_length"]).all()


def test_gff_transcript_length_matches_inclusive_exon_sum():
    gff_df, _ = parse_gff(os.path.join(TEST_DATA, "1000_entry.gff"), 1000)
    res = gff_df_to_cds_df(gff_df).set_index("transcript_id")
    exons = gff_df[gff_df["type"] == "exon"]
    inclusive = exons.groupby("transcript_id").apply(
        lambda d: int((d["end"] - d["start"] + 1).sum()), include_groups=False
    )
    common = res.index.intersection(inclusive.index)
    assert (res.loc[common, "transcript_length"] == inclusive.loc[common]).all()


def test_zero_cds_start_is_retained_and_completeness_is_explicit():
    """A CDS beginning at transcript offset zero is not automatically partial."""
    gff_df = pd.DataFrame(
        [
            ["chr", "test", "exon", 1, 300, ".", "+", ".", "transcript_id=TX_COMPLETE"],
            ["chr", "test", "CDS", 1, 300, ".", "+", "0", "transcript_id=TX_COMPLETE"],
            ["chr", "test", "start_codon", 1, 3, ".", "+", "0", "transcript_id=TX_COMPLETE"],
            ["chr", "test", "exon", 1001, 1300, ".", "+", ".", "transcript_id=TX_UNKNOWN"],
            ["chr", "test", "CDS", 1001, 1300, ".", "+", "0", "transcript_id=TX_UNKNOWN"],
        ],
        columns=[
            "seq_id",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )
    gff_df["transcript_id"] = gff_df["attributes"].str.extract(r"transcript_id=([^;]+)")

    result = gff_df_to_cds_df(gff_df).set_index("transcript_id")

    assert result.loc["TX_COMPLETE", "cds_start"] == 0
    assert bool(result.loc["TX_COMPLETE", "cds_start_complete"])
    assert result.loc["TX_UNKNOWN", "cds_start"] == 0
    assert not bool(result.loc["TX_UNKNOWN", "cds_start_complete"])


def test_gff_stop_codon_is_normalised_only_when_included_in_cds():
    """GFF3 CDS spans including stop_codon end at the stop's first base."""
    rows = [
        ["chr", "gff3", "exon", 1, 300, ".", "+", ".", "transcript_id=TX_IN"],
        ["chr", "gff3", "CDS", 1, 300, ".", "+", "0", "transcript_id=TX_IN"],
        ["chr", "gff3", "stop_codon", 298, 300, ".", "+", "0", "transcript_id=TX_IN"],
        ["chr", "gff3", "exon", 1001, 1300, ".", "+", ".", "transcript_id=TX_OUT"],
        ["chr", "gff3", "CDS", 1001, 1297, ".", "+", "0", "transcript_id=TX_OUT"],
        ["chr", "gff3", "stop_codon", 1298, 1300, ".", "+", "0", "transcript_id=TX_OUT"],
    ]
    gff_df = pd.DataFrame(
        rows,
        columns=[
            "seq_id",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )
    gff_df["transcript_id"] = gff_df["attributes"].str.extract(r"transcript_id=([^;]+)")

    result = gff_df_to_cds_df(gff_df).set_index("transcript_id")

    assert result.loc["TX_IN", "cds_end"] == 297
    assert result.loc["TX_OUT", "cds_end"] == 297


def test_frame_distribution_keeps_complete_zero_start_transcripts():
    reads = pd.DataFrame(
        {
            "read_length": [28] * 3,
            "a_site": [30, 33, 36],
            "cds_start": [0] * 3,
            "cds_end": [300] * 3,
            "cds_start_complete": [True] * 3,
            "read_name": ["r1", "r2", "r3"],
            "mapq": [255] * 3,
            "count": [1] * 3,
        }
    )
    result = read_frame_distribution_annotated(reads, exclusion_length=0)
    assert result[28] == {0: 3, 1: 0, 2: 0}


def test_frame_distribution_keeps_all_read_lengths():
    reads = pd.DataFrame(
        {
            "read_length": [18, 28, 45],
            "a_site": [30, 30, 30],
            "cds_start": [0, 0, 0],
            "cds_end": [300, 300, 300],
            "cds_start_complete": [True, True, True],
            "read_name": ["r18", "r28", "r45"],
            "mapq": [255] * 3,
            "count": [1] * 3,
        }
    )
    result = read_frame_distribution_annotated(reads, exclusion_length=0)
    assert set(result) == {18, 28, 45}


# --------------------------------------------------------------------------- #
# S6 — global offset is applied
# --------------------------------------------------------------------------- #
def test_global_offset_is_applied():
    df = pd.DataFrame(
        {
            "read_name": ["r1", "r2"],
            "read_length": [28, 30],
            "reference_start": [100, 200],
        }
    )
    out = a_site_calculation(df, offset_type="global", global_offset=12)
    assert list(out["a_site"]) == [112, 212]


# --------------------------------------------------------------------------- #
# U1 — evaluate respects metric directionality
# --------------------------------------------------------------------------- #
def test_evaluate_lower_is_better_fails_on_high_value():
    thresholds = {"duplicate_rate": {"pass": 0.3, "warn": 0.5}}
    bad = evaluate_qc_status({"metrics": {"duplicate_rate": 0.9}}, "s", thresholds)
    good = evaluate_qc_status({"metrics": {"duplicate_rate": 0.1}}, "s", thresholds)
    assert bad["overall_status"] == "FAIL"
    assert good["overall_status"] == "PASS"


def test_evaluate_explicit_direction_override():
    thresholds = {"my_metric": {"pass": 0.3, "warn": 0.5, "direction": "lower"}}
    bad = evaluate_qc_status({"metrics": {"my_metric": 0.9}}, "s", thresholds)
    assert bad["overall_status"] == "FAIL"
    assert bad["checks"][0]["direction"] == "lower"


def test_evaluate_higher_is_better_unchanged():
    thresholds = {"periodicity_dominance": {"pass": 0.7, "warn": 0.5}}
    good = evaluate_qc_status(
        {"metrics": {"periodicity_dominance": {"global": 0.9}}}, "s", thresholds
    )
    assert good["overall_status"] == "PASS"


# --------------------------------------------------------------------------- #
# E1 — no subcommand prints help and exits cleanly
# --------------------------------------------------------------------------- #
def test_cli_no_command_prints_help(capsys, monkeypatch):
    from RiboMetric import cli

    monkeypatch.setattr("sys.argv", ["RiboMetric"])
    rc = cli.main()
    assert rc == 0
    captured = capsys.readouterr()
    assert "subcommands" in captured.out or "usage" in captured.out


# --------------------------------------------------------------------------- #
# E2 — --output-offsets must run through the CLI entry point
# --------------------------------------------------------------------------- #
def test_cli_output_offsets_writes_file(tmp_path, monkeypatch):
    from RiboMetric import cli

    output_offsets = tmp_path / "applied_offsets.tsv"
    monkeypatch.setattr(
        "sys.argv",
        [
            "RiboMetric",
            "run",
            "--bam",
            os.path.join(TEST_DATA, "test.bam"),
            "--annotation",
            os.path.join(TEST_DATA, "1000_entry_RiboMetric.tsv"),
            "--output",
            str(tmp_path),
            "--config",
            os.path.join(TEST_DATA, "../../config.yml"),
            "--offset-global",
            "15",
            "--json",
            "--output-offsets",
            str(output_offsets),
            "--subsample",
            "1000",
            "--threads",
            "1",
        ],
    )

    rc = cli.main()

    assert rc == 0
    assert output_offsets.exists()
    assert (
        output_offsets.read_text()
        .splitlines()[0]
        .startswith("sample\toffset_source\toffset_target")
    )
