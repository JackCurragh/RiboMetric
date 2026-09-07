"""
Comprehensive output generation for RiboMetric

This module provides multiple output formats optimized for different use cases:
- JSON: Complete data for reanalysis
- CSV: Legacy format (kept for backwards compatibility)
- Summary TSV: Single-line summaries for pipeline integration
- QC Status JSON: Pass/warn/fail status for automated QC
- Comparison CSV: Wide format for multi-sample comparison
- Metrics Table CSV: Detailed long-format metrics
- Offsets TSV: Applied offset audit table
"""

import csv
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, cast

import pandas as pd

# =============================================================================
# Legacy Functions (Backwards Compatibility)
# =============================================================================


def generate_json(
    results_dict: dict,
    config: dict,
    name: str = "RiboMetric_data.json",
    output_directory: str = "",
) -> None:
    """
    Generate a machine readable format of the RiboMetric results
    (Legacy function - kept for backwards compatibility)

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        name: Name of the output file
        output_directory: Directory to write the output file to

    Output:
        Writes to a json file
    """
    if "sequence_slice" in results_dict:
        del results_dict["sequence_slice"]

    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name + ".json"

    data = {"results": results_dict, "config": config}

    with open(output, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Results written in {output}")


def normalise_score(score: float, min_score: float, max_score: float) -> float:
    """
    Normalise the score of a metric

    Input:
        score: The score of the metric
        min_score: The minimum score of the metric
        max_score: The maximum score of the metric

    Output:
        The normalised score
    """
    if score > max_score:
        print(f"Score {score} is greater than the maximum score {max_score}")
        return 1
    elif score < min_score:
        return 0
    else:
        return (score - min_score) / (max_score - min_score)


def generate_csv(
    results_dict: dict,
    config: dict,
    name: str = "RiboMetric_data.json",
    output_directory: str = "",
) -> None:
    """
    Generate a csv file containing the different metrics and their
    corresponding score
    (Legacy function - kept for backwards compatibility)

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        name: Name of the output file
        output_directory: Directory to write the output file to

    Output:
        Writes to a csv file
    """
    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name + ".csv"

    def _range_for(metric_key: str) -> Optional[Sequence[float]]:
        """Look up [min,max] for a metric with robust fallbacks.

        Order of preference:
          1) exact key in config["max_mins"]
          2) strip a trailing "_metric"
          3) strip a trailing "_global"
          4) None (caller should treat value as already-normalised)
        """
        mm = cast(Dict[str, Sequence[float]], config.get("max_mins", {}))
        if metric_key in mm:
            return mm[metric_key]
        base = metric_key.replace("_metric", "")
        if base in mm:
            return mm[base]
        if metric_key.endswith("_global") and metric_key[:-7] in mm:
            return mm[metric_key[:-7]]
        return None

    columns = ["Metric", "Score", "MaxMinScore"]
    metrics_dict = []
    for key, value in results_dict.get("metrics", {}).items():
        if isinstance(value, (float, int)):
            rng = _range_for(key)
            max_min_score = normalise_score(value, rng[0], rng[1]) if rng else None
            metrics_dict.append(
                {
                    "Metric": key,
                    "Score": value,
                    "MaxMinScore": max_min_score,
                }
            )
        elif isinstance(value, dict):
            rng = _range_for(key)
            for k, v in value.items():
                max_min_score = normalise_score(v, rng[0], rng[1]) if rng else None
                metrics_dict.append(
                    {
                        "Metric": f"{key}_{k}",
                        "Score": v,
                        "MaxMinScore": max_min_score,
                    }
                )

    with open(output, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()
        for data in metrics_dict:
            writer.writerow(data)
    print(f"Metrics written in {output}")


# =============================================================================
# Improved Functions (v1.0 - Pipeline Integration & Sample Review)
# =============================================================================


def generate_summary_tsv(
    results_dict: dict,
    config: dict,
    sample_name: str,
    name: str = "RiboMetric_summary.tsv",
    output_directory: str = "",
) -> None:
    """
    Generate a single-line TSV summary perfect for pipeline integration.
    Can be concatenated across samples for easy comparison.

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        sample_name: Name of the sample
        name: Name of the output file
        output_directory: Directory to write the output file to

    Output:
        Writes a TSV file with one row per sample
    """
    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name

    # Extract key global metrics
    metrics = results_dict.get("metrics", {})

    summary_row = {
        "sample": sample_name,
        "timestamp": datetime.now().isoformat(),
        "mode": results_dict.get("mode", "unknown"),
        "total_reads": sum(results_dict.get("read_length_distribution", {}).values()),
    }

    # Add global metrics (skip per-read-length metrics)
    for metric_name, metric_value in metrics.items():
        if isinstance(metric_value, dict):
            # Only add the global value
            if "global" in metric_value:
                summary_row[metric_name] = metric_value["global"]
        elif isinstance(metric_value, (int, float)):
            summary_row[metric_name] = metric_value

    # Write with headers. When appending, the existing header is authoritative:
    # the metric key set is not constant across runs (e.g. the terminal bias
    # metrics are only produced when a sequence background is available), so
    # writing this row's own key order would silently place values under the
    # wrong column names.
    path = Path(output)
    file_exists = path.exists() and path.stat().st_size > 0

    if file_exists:
        with open(output, newline="") as f:
            header = next(csv.reader(f, delimiter="\t"), [])
        if not header or any(not col for col in header):
            raise ValueError(
                f"{output} has a missing or malformed header row; refusing to "
                "append. Remove the file or write to a new one."
            )
        if len(set(header)) != len(header):
            raise ValueError(
                f"{output} has duplicate column names; refusing to append. "
                "Remove the file or write to a new one."
            )
        new_columns = [col for col in summary_row if col not in header]
        if new_columns:
            raise ValueError(
                f"{output} has no column(s) for {', '.join(sorted(new_columns))}. "
                "Appending would produce a row inconsistent with the header. "
                "Write this sample to a fresh summary file, or use the "
                "comparison CSV writer, which handles a changing metric set."
            )
        # Columns present in the header but absent from this sample are left
        # blank rather than shifting every later value one column to the left.
        fieldnames = header
        row = {col: summary_row.get(col, "") for col in header}
    else:
        fieldnames = list(summary_row)
        row = summary_row

    with open(output, "a" if file_exists else "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

    print(f"Summary line written to {output}")


# Default pass/warn thresholds used when no override YAML is supplied.
DEFAULT_QC_THRESHOLDS: Dict[str, Dict[str, float]] = {
    "periodicity_dominance": {"pass": 0.7, "warn": 0.5},
    "uniformity_entropy": {"pass": 0.7, "warn": 0.5},
    "prop_reads_CDS": {"pass": 0.7, "warn": 0.5},
    "read_length_distribution_IQR_metric": {"pass": 0.3, "warn": 0.2},
    "recommended_read_proportion": {"pass": 0.5, "warn": 0.3},
}

# Metrics where a *lower* value is better, so a sample passes when its value is
# at or below the "pass" threshold. Everything else is treated as higher-is-
# better. A thresholds YAML can override this per metric with
# ``direction: lower`` / ``direction: higher``.
LOWER_IS_BETTER_METRICS = frozenset(
    {
        "duplicate_rate",
        "multimapper_rate",
        "rpf_multimapper_rate",
        "alignment_multimapper_rate",
        "soft_clip_rate_5prime",
        "disome_proportion",
        "terminal_bias_kl_5prime_raw",
        "terminal_bias_kl_3prime_raw",
        "stop_codon_readthrough_ratio",
        # High marginal discovery = library still un-saturated (under-sequenced).
        "marginal_position_discovery_rate",
        # High FLOSS heterogeneity = more transcripts with aberrant length profiles.
        "floss_median",
        "floss_aberrant_transcript_fraction",
    }
)


def _metric_direction(metric_name: str, threshold_dict: Dict) -> str:
    """Return 'higher' or 'lower' for a metric.

    Honours an explicit ``direction`` key in the threshold spec; otherwise
    falls back to LOWER_IS_BETTER_METRICS, defaulting to 'higher'.
    """
    explicit = threshold_dict.get("direction")
    if explicit == "lower":
        return "lower"
    if explicit == "higher":
        return "higher"
    return "lower" if metric_name in LOWER_IS_BETTER_METRICS else "higher"


def _evaluate_qc_status_scored(results_dict: dict, sample_name: str) -> dict:
    """QC status from the unified scoring resolver (default path).

    The overall verdict uses only the gated (Tier-1) metrics; non-gated
    metrics are reported as checks but do not fail the sample.
    """
    from .scoring import build_scored_metrics, overall_gate_status

    scored = build_scored_metrics(results_dict, None)
    overall_status = overall_gate_status(scored)

    qc_checks = [
        {
            "metric": m["key"],
            "value": m["raw"],
            "score": m["score"],
            "status": m["status"],
            "gate": m["gate"],
            "tier": m["tier"],
        }
        for m in scored
    ]
    gated_checks = [c for c in qc_checks if c["gate"]]
    return {
        "sample": sample_name,
        "timestamp": datetime.now().isoformat(),
        "overall_status": overall_status,
        "checks": qc_checks,
        "summary": {
            "total_checks": len(qc_checks),
            "gated_checks": len(gated_checks),
            "passed": sum(1 for c in qc_checks if c["status"] == "PASS"),
            "warnings": sum(1 for c in qc_checks if c["status"] == "WARNING"),
            "failed": sum(1 for c in qc_checks if c["status"] == "FAIL"),
        },
        "recommendation": _get_recommendation(overall_status, [c for c in gated_checks]),
    }


def _validate_explicit_thresholds(thresholds: Dict) -> None:
    """Reject a threshold policy that cannot be evaluated as written.

    An explicit policy is a contract: every metric it names must actually be
    checked. A policy that is empty or malformed cannot express that, so it is
    an error rather than something to work around silently.
    """
    if not isinstance(thresholds, dict) or not thresholds:
        raise ValueError(
            "Threshold policy is empty or not a mapping; expected " "{metric: {pass: x, warn: y}}."
        )
    for metric_name, spec in thresholds.items():
        if not isinstance(spec, dict):
            raise ValueError(
                f"Threshold policy for '{metric_name}' must be a mapping with "
                f"'pass' and 'warn' keys, got {type(spec).__name__}."
            )
        for bound in ("pass", "warn"):
            if bound not in spec:
                raise ValueError(f"Threshold policy for '{metric_name}' is missing '{bound}'.")
            bound_value = spec[bound]
            if isinstance(bound_value, bool) or not isinstance(bound_value, (int, float)):
                raise ValueError(
                    f"Threshold '{bound}' for '{metric_name}' must be a number, "
                    f"got {bound_value!r}."
                )
            if not math.isfinite(bound_value):
                raise ValueError(
                    f"Threshold '{bound}' for '{metric_name}' must be finite, "
                    f"got {bound_value!r}."
                )
        direction = spec.get("direction")
        if direction is not None and direction not in ("higher", "lower"):
            raise ValueError(
                f"Threshold direction for '{metric_name}' must be 'higher' or "
                f"'lower', got {direction!r}."
            )


def _resolve_metric_value(
    metric_value: object,
) -> Tuple[Optional[float], Optional[str]]:
    """Reduce a stored metric to a comparable number.

    Returns ``(value, reason)``. ``reason`` is None when the value is usable;
    otherwise it explains why the metric could not be compared, and ``value``
    is None.
    """
    if isinstance(metric_value, dict):
        if "global" in metric_value:
            return _resolve_metric_value(metric_value["global"])
        return None, (
            "reported per read length or region only; no 'global' value to "
            "compare against the threshold"
        )
    if isinstance(metric_value, bool) or not isinstance(metric_value, (int, float)):
        return None, (f"value is not numeric ({type(metric_value).__name__}: " f"{metric_value!r})")
    if not math.isfinite(metric_value):
        return None, f"value is not finite ({metric_value!r})"
    return float(metric_value), None


def evaluate_qc_status(
    results_dict: dict,
    sample_name: str,
    thresholds: Optional[Dict] = None,
) -> dict:
    """
    Score a results dict against pass/warn thresholds.

    Pure function: performs no I/O. Used by both ``generate_qc_status`` (which
    writes the result to disk during a run) and the ``evaluate`` subcommand.

    Input:
        results_dict: Dictionary containing the qc results (expects a "metrics" key)
        sample_name: Name of the sample
        thresholds: Optional dict of {metric: {"pass": x, "warn": y}}; falls back
            to DEFAULT_QC_THRESHOLDS when None

    Output:
        Dictionary with overall_status, per-check detail, summary counts and a
        recommendation.
    """
    # Default path: use the single scoring resolver (scoring.py) so the gate
    # agrees with the report cards/badges and respects per-metric gate
    # membership. An explicit thresholds dict (e.g. an external --expected YAML
    # for the `evaluate` subcommand) keeps the legacy raw-value comparison.
    if thresholds is None:
        return _evaluate_qc_status_scored(results_dict, sample_name)

    # An explicit policy names the metrics the caller requires. A metric that is
    # absent or uncomparable is missing evidence, so it fails the gate loudly
    # instead of being skipped, which previously let an empty or misspelled
    # policy report PASS with no checks performed at all.
    _validate_explicit_thresholds(thresholds)

    metrics = results_dict.get("metrics", {})
    if not isinstance(metrics, dict):
        raise ValueError("Results 'metrics' must be a mapping, got " f"{type(metrics).__name__}.")
    qc_checks = []
    overall_status = "PASS"

    for metric_name, threshold_dict in thresholds.items():
        direction = _metric_direction(metric_name, threshold_dict)
        base_check = {
            "metric": metric_name,
            "direction": direction,
            "threshold_pass": threshold_dict["pass"],
            "threshold_warn": threshold_dict["warn"],
        }

        if metric_name not in metrics:
            qc_checks.append(
                {
                    **base_check,
                    "value": None,
                    "status": "FAIL",
                    "reason": (
                        "required metric not present in results; QC evidence is "
                        "incomplete for this sample"
                    ),
                }
            )
            overall_status = "FAIL"
            continue

        value, reason = _resolve_metric_value(metrics[metric_name])
        if reason is not None:
            qc_checks.append(
                {
                    **base_check,
                    "value": None,
                    "status": "FAIL",
                    "reason": f"required metric could not be evaluated: {reason}",
                }
            )
            overall_status = "FAIL"
            continue

        # Determine status, respecting metric directionality. For lower-is-
        # better metrics (e.g. duplicate_rate) a value at or below "pass" is a
        # PASS; treating these as higher-is-better would let a bad sample pass.
        if direction == "lower":
            if value <= threshold_dict["pass"]:
                status = "PASS"
            elif value <= threshold_dict["warn"]:
                status = "WARNING"
                if overall_status == "PASS":
                    overall_status = "WARNING"
            else:
                status = "FAIL"
                overall_status = "FAIL"
        else:
            if value >= threshold_dict["pass"]:
                status = "PASS"
            elif value >= threshold_dict["warn"]:
                status = "WARNING"
                if overall_status == "PASS":
                    overall_status = "WARNING"
            else:
                status = "FAIL"
                overall_status = "FAIL"

        qc_checks.append(
            {
                **base_check,
                "value": value,
                "status": status,
                "reason": None,
            }
        )

    return {
        "sample": sample_name,
        "timestamp": datetime.now().isoformat(),
        "overall_status": overall_status,
        "checks": qc_checks,
        "summary": {
            "total_checks": len(qc_checks),
            "passed": sum(1 for c in qc_checks if c["status"] == "PASS"),
            "warnings": sum(1 for c in qc_checks if c["status"] == "WARNING"),
            "failed": sum(1 for c in qc_checks if c["status"] == "FAIL"),
        },
        "recommendation": _get_recommendation(overall_status, qc_checks),
    }


def generate_qc_status(
    results_dict: dict,
    config: dict,
    sample_name: str,
    thresholds: Optional[Dict] = None,
    name: str = "RiboMetric_qc_status.json",
    output_directory: str = "",
) -> None:
    """
    Generate a QC status file with pass/warning/fail indicators for pipeline decision-making.

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        sample_name: Name of the sample
        thresholds: Optional dict of thresholds for pass/warning/fail
        name: Name of the output file
        output_directory: Directory to write the output file to

    Output:
        Writes a JSON file with QC status and recommendations
    """
    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name

    qc_status = evaluate_qc_status(results_dict, sample_name, thresholds)

    with open(output, "w") as f:
        json.dump(qc_status, f, indent=2)

    print(f"QC status written to {output}")
    print(f"Overall QC Status: {qc_status['overall_status']}")


def _get_recommendation(status: str, checks: List[Dict]) -> str:
    """Generate a recommendation based on QC status"""
    if status == "INFO":
        return "No gated (Tier-1) metrics were available to score; provide an annotation for an automatic QC verdict."
    if status == "PASS":
        return "Sample passed all QC checks. Proceed with downstream analysis."
    elif status == "WARNING":
        warning_metrics = [c["metric"] for c in checks if c["status"] == "WARNING"]
        return f"Sample has warnings in: {', '.join(warning_metrics)}. Review these metrics before proceeding."
    else:
        failed_metrics = [c["metric"] for c in checks if c["status"] == "FAIL"]
        return f"Sample failed QC for: {', '.join(failed_metrics)}. Consider excluding from analysis or investigating issues."


def generate_comparison_ready_csv(
    results_dict: dict,
    config: dict,
    sample_name: str,
    name: str = "RiboMetric_comparison.csv",
    output_directory: str = "",
) -> None:
    """
    Generate a CSV optimized for comparing multiple samples.
    Wide format with one row per sample, all metrics as columns.

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        sample_name: Name of the sample
        name: Name of the output file
        output_directory: Directory to write the output file to

    Output:
        Writes a CSV file ready for multi-sample comparison
    """
    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name

    metrics = results_dict.get("metrics", {})

    # Flatten all metrics
    row = {
        "sample": sample_name,
        "timestamp": datetime.now().isoformat(),
        "mode": results_dict.get("mode", "unknown"),
    }

    for metric_name, metric_value in metrics.items():
        if isinstance(metric_value, dict):
            # Add global and potentially a few key read lengths
            if "global" in metric_value:
                row[f"{metric_name}_global"] = metric_value["global"]
            # Add most abundant read lengths (28-32)
            for rl in [28, 29, 30, 31, 32]:
                if rl in metric_value:
                    row[f"{metric_name}_rl{rl}"] = metric_value[rl]
        elif isinstance(metric_value, (int, float)):
            row[metric_name] = metric_value

    # Append to file or create new
    file_exists = Path(output).exists()

    # Read existing data if file exists
    if file_exists:
        df_existing = pd.read_csv(output)
        df_new = pd.DataFrame([row])
        df_combined = pd.concat([df_existing, df_new], ignore_index=True)
        df_combined.to_csv(output, index=False)
    else:
        df = pd.DataFrame([row])
        df.to_csv(output, index=False)

    print(f"Comparison-ready CSV written to {output}")


def generate_metrics_table_csv(
    results_dict: dict,
    config: dict,
    sample_name: str,
    name: str = "RiboMetric_metrics_table.csv",
    output_directory: str = "",
) -> None:
    """
    Generate a comprehensive CSV table with all metrics, including per-read-length.
    Better organized than the original CSV format.

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        sample_name: Name of the sample
        name: Name of the output file
        output_directory: Directory to write the output file to

    Output:
        Writes a CSV file with columns: sample, metric, read_length, value, description
    """
    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name

    metrics = results_dict.get("metrics", {})
    rows = []

    # Metric descriptions for context
    descriptions = {
        "periodicity_dominance": "Proportion of reads in dominant reading frame; global uses one shared dominant frame",
        "uniformity_entropy": "Shannon entropy of codon-binned start-codon metagene window",
        "read_length_distribution_IQR_metric": "Normalized IQR of read length distribution",
        "terminal_nucleotide_bias_distribution_5_prime_metric": "Normalized 5' terminal bias score from KL divergence",
        "terminal_nucleotide_bias_distribution_3_prime_metric": "Normalized 3' terminal bias score from KL divergence",
        "terminal_bias_kl_5prime": "Normalized 5' terminal bias score from KL divergence",
        "terminal_bias_kl_3prime": "Normalized 3' terminal bias score from KL divergence",
        "terminal_bias_kl_5prime_score": "Normalized 5' terminal bias score from KL divergence",
        "terminal_bias_kl_3prime_score": "Normalized 3' terminal bias score from KL divergence",
        "terminal_bias_kl_5prime_raw": "Raw 5' terminal KL divergence in bits",
        "terminal_bias_kl_3prime_raw": "Raw 3' terminal KL divergence in bits",
        "CDS_coverage_metric": "Proportion of CDS covered by reads",
        "prop_reads_CDS": "Proportion of reads mapping to CDS",
        "prop_reads_leader": "Proportion of reads mapping to 5' leader",
        "prop_reads_trailer": "Proportion of reads mapping to 3' trailer",
        "ratio_cds:leader": "Ratio of CDS to 5'UTR reads",
    }

    for metric_name, metric_value in metrics.items():
        desc = descriptions.get(metric_name, "")

        if isinstance(metric_value, dict):
            # Per-read-length or per-region metrics
            for key, value in metric_value.items():
                rows.append(
                    {
                        "sample": sample_name,
                        "metric": metric_name,
                        "read_length_or_region": str(key),
                        "value": value,
                        "description": desc,
                    }
                )
        elif isinstance(metric_value, (int, float)):
            # Global metrics
            rows.append(
                {
                    "sample": sample_name,
                    "metric": metric_name,
                    "read_length_or_region": "global",
                    "value": metric_value,
                    "description": desc,
                }
            )

    # Write CSV
    df = pd.DataFrame(rows)
    df.to_csv(output, index=False)

    print(f"Detailed metrics table written to {output}")


def generate_offsets_tsv(
    results_dict: dict,
    sample_name: str,
    name: str = "RiboMetric_offsets.tsv",
    output_directory: str = "",
) -> None:
    """Write a flat table of the offsets actually applied during the run."""
    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name

    offsets = results_dict.get("offsets", {})
    applied = offsets.get("applied_by_read_length", {}) or {}
    computed = offsets.get("computed_offsets", {}) or {}
    final = offsets.get("final_offsets", computed) or {}
    raw = offsets.get("raw_offsets", {}) or {}
    frame_adjustments = offsets.get("frame_adjustments", {}) or {}
    columns = [
        "sample",
        "offset_source",
        "offset_target",
        "offset_calculation_method",
        "read_length",
        "n_reads",
        "n_unique_offsets",
        "applied_offsets",
        "raw_offset",
        "final_offset",
        "min_offset",
        "max_offset",
        "computed_offset",
        "global_offset",
        "frame_adjusted",
        "old_offset",
        "new_offset",
        "dominant_frame",
        "dominant_fraction",
        "frame_adjustment_reads",
    ]
    rows = []

    if applied:
        for read_length, record in sorted(applied.items(), key=lambda item: int(item[0])):
            adjustment = frame_adjustments.get(str(read_length), {})
            offset_values = record.get("offsets", [])
            rows.append(
                {
                    "sample": sample_name,
                    "offset_source": offsets.get("source"),
                    "offset_target": offsets.get("target"),
                    "offset_calculation_method": offsets.get("offset_calculation_method"),
                    "read_length": read_length,
                    "n_reads": record.get("n_reads"),
                    "n_unique_offsets": record.get("n_unique_offsets"),
                    "applied_offsets": "|".join(str(v) for v in offset_values),
                    "min_offset": record.get("min_offset"),
                    "max_offset": record.get("max_offset"),
                    "raw_offset": raw.get(str(read_length)),
                    "final_offset": final.get(str(read_length)),
                    "computed_offset": final.get(str(read_length), computed.get(str(read_length))),
                    "global_offset": offsets.get("global_offset"),
                    "frame_adjusted": bool(adjustment),
                    "old_offset": adjustment.get("old_offset"),
                    "new_offset": adjustment.get("new_offset"),
                    "dominant_frame": adjustment.get("dominant_frame"),
                    "dominant_fraction": adjustment.get("dominant_fraction"),
                    "frame_adjustment_reads": adjustment.get("reads"),
                }
            )
    else:
        rows.append(
            {
                "sample": sample_name,
                "offset_source": offsets.get("source"),
                "offset_target": offsets.get("target"),
                "offset_calculation_method": offsets.get("offset_calculation_method"),
                "read_length": "",
                "n_reads": "",
                "n_unique_offsets": "",
                "applied_offsets": "",
                "min_offset": "",
                "max_offset": "",
                "computed_offset": "",
                "global_offset": offsets.get("global_offset"),
                "frame_adjusted": False,
                "old_offset": "",
                "new_offset": "",
                "dominant_frame": "",
                "dominant_fraction": "",
                "frame_adjustment_reads": "",
            }
        )

    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Offsets TSV written to {output}")


def generate_all_outputs(
    results_dict: dict,
    config: dict,
    sample_name: str,
    output_directory: str = "",
    thresholds: Optional[Dict] = None,
) -> None:
    """
    Convenience function to generate all improved output formats at once.

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        sample_name: Name of the sample
        output_directory: Directory to write the output files to
        thresholds: Optional dict of thresholds for QC status

    Output:
        Generates all output files
    """
    print("Generating improved output formats...")

    generate_summary_tsv(
        results_dict, config, sample_name, f"{sample_name}_summary.tsv", output_directory
    )

    generate_metrics_table_csv(
        results_dict, config, sample_name, f"{sample_name}_metrics_table.csv", output_directory
    )

    generate_qc_status(
        results_dict,
        config,
        sample_name,
        thresholds,
        f"{sample_name}_qc_status.json",
        output_directory,
    )

    generate_comparison_ready_csv(
        results_dict, config, sample_name, f"{sample_name}_comparison.csv", output_directory
    )

    generate_offsets_tsv(results_dict, sample_name, f"{sample_name}_offsets.tsv", output_directory)

    print("All improved outputs generated successfully!")
