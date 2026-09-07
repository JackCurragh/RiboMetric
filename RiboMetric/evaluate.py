"""
The ``evaluate`` subcommand.

Takes a previously produced RiboMetric result (JSON or metrics-table CSV) and a
YAML file of expected metric thresholds, and reports whether the sample passes,
warns, or fails. Intended for pipeline gating: the process exit code reflects the
outcome (0 = PASS, 1 = WARNING, 2 = FAIL) so it can be branched on in a workflow.
"""
import csv
import json
from pathlib import Path
from argparse import Namespace
from typing import Any, Dict, cast

import yaml

from .results_output import evaluate_qc_status, DEFAULT_QC_THRESHOLDS


# Exit codes used to gate downstream pipeline steps.
EXIT_PASS = 0
EXIT_WARNING = 1
EXIT_FAIL = 2


def _load_results(path: Path) -> Dict[str, Any]:
    """Load a results file into a dict with a top-level "metrics" key.

    Supports RiboMetric JSON output ({"results": {"metrics": ...}} or a bare
    results dict) and a metrics-table CSV (columns: metric, read_length_or_region,
    value).
    """
    suffix = path.suffix.lower()
    if suffix == ".json":
        with open(path) as f:
            data = json.load(f)
        # JSON output is wrapped as {"results": ..., "config": ...}
        results = data.get("results", data)
        if "metrics" not in results:
            raise ValueError(
                f"{path} does not contain a 'metrics' section; "
                "is it a RiboMetric JSON output?"
            )
        return cast(Dict[str, Any], results)
    if suffix == ".csv":
        metrics: Dict[str, Any] = {}
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            required = {"metric", "value"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(
                    f"{path} must have at least 'metric' and 'value' columns "
                    "(a RiboMetric metrics-table CSV)."
                )
            for row in reader:
                name = row["metric"]
                region = row.get("read_length_or_region", "global") or "global"
                try:
                    value: Any = float(row["value"])
                except (TypeError, ValueError):
                    continue
                metrics.setdefault(name, {})[region] = value
        # Collapse single-global metrics to scalars for cleaner downstream use
        for name, by_region in list(metrics.items()):
            if set(by_region) == {"global"}:
                metrics[name] = by_region["global"]
        return {"metrics": metrics}
    raise ValueError(f"Unsupported results file type: {suffix} (expected .json or .csv)")


def _load_thresholds(path: Path) -> Dict[str, Dict[str, float]]:
    """Load a thresholds YAML.

    Accepts either a top-level mapping of {metric: {pass, warn}} or that mapping
    nested under a "thresholds:" key.
    """
    with open(path) as f:
        data = yaml.safe_load(f)
    if isinstance(data, dict) and "thresholds" in data:
        data = data["thresholds"]
    if not isinstance(data, dict):
        raise ValueError(
            f"{path} must contain a mapping of metric -> {{pass, warn}} thresholds."
        )
    return cast(Dict[str, Dict[str, float]], data)


def evaluate(args: Namespace) -> int:
    """Entry point for ``RiboMetric evaluate``. Returns a gating exit code."""
    results_path = Path(args.input)
    if not results_path.exists():
        print(f"Error: results file not found: {results_path}")
        return EXIT_FAIL

    if getattr(args, "expected", None):
        thresholds = _load_thresholds(Path(args.expected))
    else:
        print("No --expected thresholds provided; using built-in defaults.")
        thresholds = DEFAULT_QC_THRESHOLDS

    results = _load_results(results_path)
    sample_name = getattr(args, "name", None) or results_path.stem

    try:
        status = evaluate_qc_status(results, sample_name, thresholds)
    except ValueError as exc:
        # A policy that cannot be evaluated is a configuration error, not a
        # verdict on the sample. Exit non-zero so a pipeline stops rather than
        # treating an unusable policy as a pass.
        print(f"Error: {exc}")
        return EXIT_FAIL

    # Human-readable report
    print(f"\nQC evaluation for '{sample_name}': {status['overall_status']}")
    for check in status["checks"]:
        cmp = "<=" if check.get("direction") == "lower" else ">="
        if check["value"] is None:
            print(
                f"  [{check['status']:<7}] {check['metric']} = n/a "
                f"({check.get('reason', 'no value available')})"
            )
            continue
        print(
            f"  [{check['status']:<7}] {check['metric']} = "
            f"{check['value']:.4g} "
            f"(pass{cmp}{check['threshold_pass']}, warn{cmp}{check['threshold_warn']})"
        )
    if not status["checks"]:
        print("  (no thresholded metrics were found in the results file)")
    print(f"\n{status['recommendation']}")

    if getattr(args, "output", None):
        with open(args.output, "w") as f:
            json.dump(status, f, indent=2)
        print(f"\nEvaluation written to {args.output}")

    return {
        "PASS": EXIT_PASS,
        "WARNING": EXIT_WARNING,
        "FAIL": EXIT_FAIL,
    }[status["overall_status"]]
