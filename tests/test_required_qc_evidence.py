"""
Regression tests for two fail-quiet bugs in QC gating and summary export.

Both are reachable from an ordinary cohort run: ``qc.py`` only populates the
terminal nucleotide bias metrics when a sequence background is available, so
the metric key set differs between a run with ``--fasta`` and one without.
"""

import csv
import math
from argparse import Namespace

import pytest

from RiboMetric.evaluate import evaluate, EXIT_PASS, EXIT_FAIL
from RiboMetric.results_output import evaluate_qc_status, generate_summary_tsv


def _results(metrics):
    return {
        "mode": "annotated",
        "read_length_distribution": {28: 100},
        "metrics": metrics,
    }


# ---------------------------------------------------------------------------
# Explicit threshold policies must evaluate every metric they name
# ---------------------------------------------------------------------------

def test_missing_required_metric_fails_rather_than_passing():
    """A policy naming a metric the run never produced must not report PASS."""
    status = evaluate_qc_status(
        _results({"periodicity_dominance": 0.12}),
        "no_fasta_run",
        {"terminal_bias_kl_5prime": {"pass": 0.8, "warn": 0.6}},
    )
    assert status["overall_status"] == "FAIL"
    assert status["summary"]["total_checks"] == 1
    check = status["checks"][0]
    assert check["status"] == "FAIL"
    assert check["value"] is None
    assert "not present" in check["reason"]


def test_misspelled_metric_in_policy_fails():
    status = evaluate_qc_status(
        _results({"periodicity_dominance": 0.02}),
        "terrible",
        {"periodicty_dominance": {"pass": 0.7, "warn": 0.5}},
    )
    assert status["overall_status"] == "FAIL"
    assert status["summary"]["failed"] == 1


@pytest.mark.parametrize("bad_value", [None, "n/a", float("nan"), float("inf")])
def test_unusable_metric_values_fail(bad_value):
    status = evaluate_qc_status(
        _results({"periodicity_dominance": bad_value}),
        "sample",
        {"periodicity_dominance": {"pass": 0.7, "warn": 0.5}},
    )
    assert status["overall_status"] == "FAIL"
    assert status["checks"][0]["value"] is None
    assert status["checks"][0]["reason"]


def test_per_read_length_metric_without_global_fails():
    status = evaluate_qc_status(
        _results({"periodicity_dominance": {"28": 0.9, "29": 0.8}}),
        "sample",
        {"periodicity_dominance": {"pass": 0.7, "warn": 0.5}},
    )
    assert status["overall_status"] == "FAIL"
    assert "global" in status["checks"][0]["reason"]


@pytest.mark.parametrize("policy", [
    {},
    {"periodicity_dominance": 0.7},
    {"periodicity_dominance": {"pass": 0.7}},
    {"periodicity_dominance": {"pass": "high", "warn": 0.5}},
    {"periodicity_dominance": {"pass": float("nan"), "warn": 0.5}},
    {"periodicity_dominance": {"pass": 0.7, "warn": 0.5, "direction": "up"}},
])
def test_malformed_policies_are_rejected(policy):
    with pytest.raises(ValueError):
        evaluate_qc_status(_results({"periodicity_dominance": 0.9}), "s", policy)


# ---------------------------------------------------------------------------
# Existing behaviour for usable values must be unchanged
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value,expected", [
    (0.90, "PASS"), (0.60, "WARNING"), (0.30, "FAIL"),
])
def test_higher_is_better_ordering_preserved(value, expected):
    status = evaluate_qc_status(
        _results({"periodicity_dominance": value}),
        "s",
        {"periodicity_dominance": {"pass": 0.7, "warn": 0.5}},
    )
    assert status["overall_status"] == expected
    assert status["checks"][0]["value"] == pytest.approx(value)


@pytest.mark.parametrize("value,expected", [
    (0.05, "PASS"), (0.25, "WARNING"), (0.80, "FAIL"),
])
def test_lower_is_better_direction_preserved(value, expected):
    status = evaluate_qc_status(
        _results({"duplicate_rate": value}),
        "s",
        {"duplicate_rate": {"pass": 0.1, "warn": 0.3}},
    )
    assert status["overall_status"] == expected


def test_global_key_still_resolved():
    status = evaluate_qc_status(
        _results({"periodicity_dominance": {"global": 0.9, "28": 0.4}}),
        "s",
        {"periodicity_dominance": {"pass": 0.7, "warn": 0.5}},
    )
    assert status["overall_status"] == "PASS"


def test_none_thresholds_still_delegates_to_scored_resolver():
    """The default (unscored) path must be untouched by the required-check rule."""
    status = evaluate_qc_status(
        _results({"periodicity_dominance": 0.9}), "s", None
    )
    assert "overall_status" in status
    assert status["overall_status"] in {"PASS", "WARNING", "FAIL"}


# ---------------------------------------------------------------------------
# CLI exit codes
# ---------------------------------------------------------------------------

def _write_json(tmp_path, metrics):
    import json
    path = tmp_path / "results.json"
    path.write_text(json.dumps({"results": _results(metrics)}))
    return path


def test_cli_fails_on_missing_required_metric(tmp_path):
    results = _write_json(tmp_path, {"periodicity_dominance": 0.9})
    expected = tmp_path / "expected.yml"
    expected.write_text("terminal_bias_kl_5prime:\n  pass: 0.8\n  warn: 0.6\n")
    code = evaluate(Namespace(
        input=str(results), expected=str(expected), name="s", output=None
    ))
    assert code == EXIT_FAIL


def test_cli_fails_on_malformed_policy(tmp_path):
    results = _write_json(tmp_path, {"periodicity_dominance": 0.9})
    expected = tmp_path / "expected.yml"
    expected.write_text("periodicity_dominance:\n  pass: 0.7\n")
    code = evaluate(Namespace(
        input=str(results), expected=str(expected), name="s", output=None
    ))
    assert code == EXIT_FAIL


def test_cli_passes_a_good_sample(tmp_path):
    results = _write_json(tmp_path, {"periodicity_dominance": 0.9})
    expected = tmp_path / "expected.yml"
    expected.write_text("periodicity_dominance:\n  pass: 0.7\n  warn: 0.5\n")
    code = evaluate(Namespace(
        input=str(results), expected=str(expected), name="s", output=None
    ))
    assert code == EXIT_PASS


# ---------------------------------------------------------------------------
# Summary TSV column identity
# ---------------------------------------------------------------------------

def _read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_append_keeps_values_under_their_own_columns(tmp_path):
    """Sample B lacks the terminal bias metrics; its values must not shift left."""
    with_fasta = _results({
        "periodicity_dominance": 0.91,
        "terminal_bias_kl_5prime": 0.80,
        "uniformity_entropy": 0.55,
    })
    without_fasta = _results({
        "periodicity_dominance": 0.42,
        "uniformity_entropy": 0.31,
    })
    generate_summary_tsv(with_fasta, {}, "sampleA", "s.tsv", str(tmp_path))
    generate_summary_tsv(without_fasta, {}, "sampleB", "s.tsv", str(tmp_path))

    rows = _read_rows(tmp_path / "s.tsv")
    assert rows[1]["sample"] == "sampleB"
    assert rows[1]["periodicity_dominance"] == "0.42"
    assert rows[1]["uniformity_entropy"] == "0.31"
    assert rows[1]["terminal_bias_kl_5prime"] == ""


def test_append_respects_header_order_not_row_order(tmp_path):
    first = _results({"periodicity_dominance": 0.91, "uniformity_entropy": 0.55})
    reordered = _results({"uniformity_entropy": 0.21, "periodicity_dominance": 0.74})
    generate_summary_tsv(first, {}, "sampleA", "s.tsv", str(tmp_path))
    generate_summary_tsv(reordered, {}, "sampleB", "s.tsv", str(tmp_path))

    rows = _read_rows(tmp_path / "s.tsv")
    assert rows[1]["periodicity_dominance"] == "0.74"
    assert rows[1]["uniformity_entropy"] == "0.21"


def test_new_metric_column_is_rejected(tmp_path):
    generate_summary_tsv(
        _results({"periodicity_dominance": 0.9}), {}, "sampleA", "s.tsv", str(tmp_path)
    )
    with pytest.raises(ValueError, match="fresh summary file"):
        generate_summary_tsv(
            _results({"periodicity_dominance": 0.5, "brand_new_metric": 0.1}),
            {}, "sampleB", "s.tsv", str(tmp_path),
        )


def test_malformed_header_is_rejected(tmp_path):
    path = tmp_path / "s.tsv"
    path.write_text("sample\t\tmode\n")
    with pytest.raises(ValueError, match="malformed header"):
        generate_summary_tsv(
            _results({"periodicity_dominance": 0.9}), {}, "sampleA", "s.tsv", str(tmp_path)
        )


def test_first_write_and_matching_append_still_work(tmp_path):
    for name, value in (("sampleA", 0.91), ("sampleB", 0.42)):
        generate_summary_tsv(
            _results({"periodicity_dominance": value}), {}, name, "s.tsv", str(tmp_path)
        )
    rows = _read_rows(tmp_path / "s.tsv")
    assert [r["sample"] for r in rows] == ["sampleA", "sampleB"]
    assert [r["periodicity_dominance"] for r in rows] == ["0.91", "0.42"]
    assert not math.isnan(float(rows[0]["periodicity_dominance"]))
