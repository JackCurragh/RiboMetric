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
import pytest

from RiboMetric.bam_processing import process_sequences
from RiboMetric.file_parser import parse_gff, gff_df_to_cds_df
from RiboMetric.modules import a_site_calculation
from RiboMetric.results_output import evaluate_qc_status


TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")


# --------------------------------------------------------------------------- #
# S1 — terminal-bias 5' vs 3' background must differ
# --------------------------------------------------------------------------- #
def test_five_and_three_prime_backgrounds_differ():
    seqs = ["AACCGGTT", "TTGGCCAA", "ACGTACGT", "TACGTACG"]
    out = process_sequences(seqs, [1, 1, 1, 1], pattern_length=2)
    diff = sum(abs(out["5_prime_bg"][k] - out["3_prime_bg"][k])
               for k in out["5_prime_bg"])
    assert diff > 0, "5' and 3' backgrounds are identical (S1 regression)"


def test_three_prime_background_excludes_terminal_pattern():
    # Reads whose only variation is at the 3' end should produce a 3' background
    # that differs from the 5' background (the 3' terminal pattern is excluded).
    seqs = ["ACGTACGT", "ACGTACGA"]
    out = process_sequences(seqs, [1, 1], pattern_length=2)
    diff = sum(abs(out["5_prime_bg"][k] - out["3_prime_bg"][k])
               for k in out["5_prime_bg"])
    assert diff > 0


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


# --------------------------------------------------------------------------- #
# S6 — global offset is applied
# --------------------------------------------------------------------------- #
def test_global_offset_is_applied():
    df = pd.DataFrame({
        "read_name": ["r1", "r2"],
        "read_length": [28, 30],
        "reference_start": [100, 200],
    })
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
    assert output_offsets.read_text().splitlines()[0].startswith(
        "sample\toffset_source\toffset_target"
    )
