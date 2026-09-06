"""
Comprehensive output generation for RiboMetric

This module provides multiple output formats optimized for different use cases:
- JSON: Complete data for reanalysis
- CSV: Legacy format (kept for backwards compatibility)
- Summary TSV: Single-line summaries for pipeline integration
- QC Status JSON: Pass/warn/fail status for automated QC
- Comparison CSV: Wide format for multi-sample comparison
- Metrics Table CSV: Detailed long-format metrics
"""

import json
import csv
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Sequence, cast


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
            max_min_score = (
                normalise_score(value, rng[0], rng[1]) if rng else None
            )
            metrics_dict.append({
                "Metric": key,
                "Score": value,
                "MaxMinScore": max_min_score,
            })
        elif isinstance(value, dict):
            rng = _range_for(key)
            for k, v in value.items():
                max_min_score = (
                    normalise_score(v, rng[0], rng[1]) if rng else None
                )
                metrics_dict.append({
                    "Metric": f"{key}_{k}",
                    "Score": v,
                    "MaxMinScore": max_min_score,
                })

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
    """Append metrics by saved column name, never by current dictionary order.

    Null metrics retain their column. Missing existing metrics are written as
    blank; new columns require a new file or the comparison-CSV writer instead
    of silently changing the schema of an existing TSV. This is a sequential
    writer, not a concurrent multi-process append API.
    """
    output = Path(output_directory) / name if output_directory else Path(name)
    summary_row = {
        "sample": sample_name,
        "timestamp": datetime.now().isoformat(),
        "mode": results_dict.get("mode", "unknown"),
        "total_reads": sum(results_dict.get("read_length_distribution", {}).values()),
    }
    for metric_name, metric_value in results_dict.get("metrics", {}).items():
        if metric_name in summary_row:
            raise ValueError(f"Metric collides with summary metadata: {metric_name}")
        if isinstance(metric_value, dict):
            if "global" in metric_value:
                summary_row[metric_name] = metric_value["global"]
        elif isinstance(metric_value, (int, float)) or metric_value is None:
            summary_row[metric_name] = metric_value

    file_exists = output.exists() and output.stat().st_size > 0
    if file_exists:
        with output.open(newline="") as f:
            header = next(csv.reader(f, delimiter="\t"), [])
        if (not header or len(header) != len(set(header))
                or not {"sample", "timestamp", "mode", "total_reads"}.issubset(header)):
            raise ValueError("Existing summary TSV has an invalid header")
        extra = set(summary_row) - set(header)
        if extra:
            raise ValueError(
                "Cannot append new metrics to existing summary TSV: "
                + ", ".join(sorted(extra)) + ". Use a new file or comparison CSV."
            )
    else:
        header = list(summary_row)
    # Header validation happens before the file is opened for writing.
    with output.open("a" if file_exists else "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        if not file_exists:
            writer.writeheader()
        writer.writerow(summary_row)
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
LOWER_IS_BETTER_METRICS = frozenset({
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
})


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
        "recommendation": _get_recommendation(
            overall_status, [c for c in gated_checks]
        ),
    }


def evaluate_qc_status(
    results_dict: dict,
    sample_name: str,
    thresholds: Optional[Dict] = None,
) -> dict:
    """Evaluate every explicitly requested check; unavailable is not PASS.

    ``thresholds=None`` retains the unified scoring resolver. An explicit
    threshold mapping is a required-check contract: missing/nonnumeric metrics
    fail the gate and invalid or empty policies raise ValueError. This failure
    indicates incomplete evidence, not necessarily a poor biological sample.
    Raw KL nulls with the documented zero-background status represent positive
    infinity; they fail a finite lower-is-better threshold without putting an
    Infinity token into the JSON report.
    """
    import math
    from numbers import Real

    if thresholds is None:
        return _evaluate_qc_status_scored(results_dict, sample_name)
    if not isinstance(thresholds, dict) or not thresholds:
        raise ValueError("Explicit QC thresholds must be a non-empty mapping")
    if not isinstance(results_dict, dict) or not isinstance(results_dict.get("metrics"), dict):
        raise ValueError("QC results must contain a metrics mapping")

    # Validate the complete policy before evaluating any values. A misspelled
    # direction or reversed thresholds must not silently weaken the gate.
    policy = []
    for name, spec in thresholds.items():
        if not isinstance(name, str) or not name or not isinstance(spec, dict):
            raise ValueError("Each QC threshold needs a metric name and mapping")
        if set(spec) - {"pass", "warn", "direction"}:
            raise ValueError(f"Unknown threshold option for {name}")
        if "direction" in spec and spec["direction"] not in ("lower", "higher"):
            raise ValueError(f"Invalid threshold direction for {name}")
        values = []
        for key in ("pass", "warn"):
            value = spec.get(key)
            if isinstance(value, bool) or not isinstance(value, Real):
                raise ValueError(f"Threshold {name}.{key} must be a finite number")
            try:
                numeric = float(value)
            except (ValueError, OverflowError) as exc:
                raise ValueError(f"Invalid threshold {name}.{key}") from exc
            if not math.isfinite(numeric):
                raise ValueError(f"Threshold {name}.{key} must be finite")
            values.append(numeric)
        passed, warned = values
        direction = _metric_direction(name, spec)
        if (direction == "lower" and passed > warned) or (direction == "higher" and passed < warned):
            raise ValueError(f"Inconsistent pass/warn ordering for {name}")
        policy.append((name, direction, passed, warned))

    metrics = results_dict["metrics"]
    checks = []
    infinite_kl_keys = {"terminal_bias_kl_5prime_raw", "terminal_bias_kl_3prime_raw"}
    for name, direction, passed, warned in policy:
        reason = None
        raw_status = None
        value = None
        raw_value = metrics.get(name)
        if name not in metrics:
            reason = "missing_metric"
        elif isinstance(raw_value, dict):
            if "global" not in raw_value:
                reason = "missing_global_value"
            else:
                raw_value = raw_value["global"]

        if reason is None:
            if (name in infinite_kl_keys and raw_value is None
                    and metrics.get(name + "_status") == "zero_background_support"):
                value = math.inf
                raw_status = "zero_background_support"
            elif raw_value is None:
                reason = "unavailable_metric"
            elif isinstance(raw_value, bool) or not isinstance(raw_value, Real):
                reason = "nonnumeric_metric"
            else:
                try:
                    value = float(raw_value)
                except (ValueError, OverflowError):
                    reason = "nonfinite_metric"
                if reason is None and not math.isfinite(value):
                    reason = "nonfinite_metric"

        if reason is not None:
            status = "FAIL"
            report_value = None
        else:
            pass_ok = value <= passed if direction == "lower" else value >= passed
            warn_ok = value <= warned if direction == "lower" else value >= warned
            status = "PASS" if pass_ok else "WARNING" if warn_ok else "FAIL"
            report_value = None if raw_status else value
        check = {
            "metric": name,
            "value": report_value,
            "status": status,
            "direction": direction,
            "threshold_pass": passed,
            "threshold_warn": warned,
        }
        if reason is not None:
            check["reason"] = reason
        if raw_status is not None:
            check["value_status"] = raw_status
        checks.append(check)

    overall = "FAIL" if any(c["status"] == "FAIL" for c in checks) else (
        "WARNING" if any(c["status"] == "WARNING" for c in checks) else "PASS"
    )
    unavailable = [c["metric"] for c in checks if "reason" in c]
    recommendation = _get_recommendation(overall, checks)
    if unavailable:
        recommendation = (
            "Required QC evidence is missing or invalid for: " + ", ".join(unavailable)
            + ". The sample has not passed the requested QC policy. "
            "Check metric names and regenerate the required measurements before proceeding."
        )
    return {
        "sample": sample_name,
        "timestamp": datetime.now().isoformat(),
        "overall_status": overall,
        "checks": checks,
        "summary": {
            "total_checks": len(checks),
            "passed": sum(c["status"] == "PASS" for c in checks),
            "warnings": sum(c["status"] == "WARNING" for c in checks),
            "failed": sum(c["status"] == "FAIL" for c in checks),
            "unavailable": len(unavailable),
        },
        "recommendation": recommendation,
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
                rows.append({
                    "sample": sample_name,
                    "metric": metric_name,
                    "read_length_or_region": str(key),
                    "value": value,
                    "description": desc
                })
        elif isinstance(metric_value, (int, float)):
            # Global metrics
            rows.append({
                "sample": sample_name,
                "metric": metric_name,
                "read_length_or_region": "global",
                "value": metric_value,
                "description": desc
            })

    # Write CSV
    df = pd.DataFrame(rows)
    df.to_csv(output, index=False)

    print(f"Detailed metrics table written to {output}")


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
        results_dict, config, sample_name,
        f"{sample_name}_summary.tsv", output_directory
    )

    generate_metrics_table_csv(
        results_dict, config, sample_name,
        f"{sample_name}_metrics_table.csv", output_directory
    )

    generate_qc_status(
        results_dict, config, sample_name, thresholds,
        f"{sample_name}_qc_status.json", output_directory
    )

    generate_comparison_ready_csv(
        results_dict, config, sample_name,
        f"{sample_name}_comparison.csv", output_directory
    )

    print("All improved outputs generated successfully!")
