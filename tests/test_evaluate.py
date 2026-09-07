"""Tests for the `evaluate` subcommand and its QC scoring helpers."""
import json
from argparse import Namespace

import pytest

from RiboMetric.results_output import evaluate_qc_status, DEFAULT_QC_THRESHOLDS
from RiboMetric.evaluate import (
    evaluate,
    _load_results,
    _load_thresholds,
    EXIT_PASS,
    EXIT_WARNING,
    EXIT_FAIL,
)


THRESHOLDS = {
    "periodicity_dominance": {"pass": 0.7, "warn": 0.5},
    "prop_reads_CDS": {"pass": 0.7, "warn": 0.5},
}


def test_evaluate_qc_status_pass():
    results = {"metrics": {"periodicity_dominance": {"global": 0.8},
                           "prop_reads_CDS": 0.9}}
    status = evaluate_qc_status(results, "s", THRESHOLDS)
    assert status["overall_status"] == "PASS"
    assert status["summary"]["passed"] == 2


def test_evaluate_qc_status_warning():
    results = {"metrics": {"periodicity_dominance": 0.6, "prop_reads_CDS": 0.9}}
    status = evaluate_qc_status(results, "s", THRESHOLDS)
    assert status["overall_status"] == "WARNING"


def test_evaluate_qc_status_fail():
    results = {"metrics": {"periodicity_dominance": 0.3, "prop_reads_CDS": 0.9}}
    status = evaluate_qc_status(results, "s", THRESHOLDS)
    assert status["overall_status"] == "FAIL"


def test_evaluate_qc_status_fails_on_missing_metrics():
    """An explicit policy is a required-check contract.

    Previously every named metric was skipped when absent, so a results file
    with no metrics at all reported PASS with zero checks performed.
    """
    status = evaluate_qc_status({"metrics": {}}, "s", THRESHOLDS)
    assert status["overall_status"] == "FAIL"
    assert len(status["checks"]) == len(THRESHOLDS)
    assert all(c["status"] == "FAIL" for c in status["checks"])
    assert all(c["value"] is None and c["reason"] for c in status["checks"])


def test_load_results_json(tmp_path):
    p = tmp_path / "r.json"
    p.write_text(json.dumps(
        {"results": {"metrics": {"prop_reads_CDS": 0.5}}, "config": {}}))
    results = _load_results(p)
    assert results["metrics"]["prop_reads_CDS"] == 0.5


def test_load_results_json_bare(tmp_path):
    p = tmp_path / "r.json"
    p.write_text(json.dumps({"metrics": {"prop_reads_CDS": 0.5}}))
    assert _load_results(p)["metrics"]["prop_reads_CDS"] == 0.5


def test_load_results_json_invalid(tmp_path):
    p = tmp_path / "r.json"
    p.write_text(json.dumps({"foo": "bar"}))
    with pytest.raises(ValueError):
        _load_results(p)


def test_load_results_csv_collapses_global(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text(
        "sample,metric,read_length_or_region,value,description\n"
        "s,prop_reads_CDS,global,0.5,x\n"
        "s,periodicity_dominance,28,0.7,x\n"
        "s,periodicity_dominance,29,0.8,x\n"
    )
    results = _load_results(p)
    # single-global metric collapses to a scalar
    assert results["metrics"]["prop_reads_CDS"] == 0.5
    # per-region metric stays a dict
    assert results["metrics"]["periodicity_dominance"] == {"28": 0.7, "29": 0.8}


def test_load_results_unsupported(tmp_path):
    p = tmp_path / "r.txt"
    p.write_text("nope")
    with pytest.raises(ValueError):
        _load_results(p)


def test_load_thresholds_nested_and_flat(tmp_path):
    nested = tmp_path / "n.yml"
    nested.write_text("thresholds:\n  prop_reads_CDS:\n    pass: 0.7\n    warn: 0.5\n")
    flat = tmp_path / "f.yml"
    flat.write_text("prop_reads_CDS:\n  pass: 0.7\n  warn: 0.5\n")
    assert _load_thresholds(nested) == _load_thresholds(flat)


def test_evaluate_exit_codes(tmp_path):
    res = tmp_path / "r.json"
    res.write_text(json.dumps({"metrics": {"periodicity_dominance": 0.9,
                                           "prop_reads_CDS": 0.9}}))
    thr = tmp_path / "t.yml"
    thr.write_text("prop_reads_CDS:\n  pass: 0.7\n  warn: 0.5\n")

    args = Namespace(input=str(res), expected=str(thr), output=None, name=None)
    assert evaluate(args) == EXIT_PASS

    res.write_text(json.dumps({"metrics": {"prop_reads_CDS": 0.3}}))
    assert evaluate(args) == EXIT_FAIL


def test_evaluate_missing_file(tmp_path):
    args = Namespace(input=str(tmp_path / "nope.json"),
                     expected=None, output=None, name=None)
    assert evaluate(args) == EXIT_FAIL


def test_default_thresholds_used_when_none():
    status = evaluate_qc_status(
        {"metrics": {"periodicity_dominance": 0.9}}, "s", None)
    # periodicity_dominance is in the defaults
    assert any(c["metric"] == "periodicity_dominance" for c in status["checks"])
    assert "periodicity_dominance" in DEFAULT_QC_THRESHOLDS
