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
from typing import Dict, List, Optional


# =============================================================================
# Legacy Functions (Backwards Compatibility)
# =============================================================================

def generate_json(
    results_dict: dict,
    config: dict,
    name: str = "RiboMetric_data.json",
    output_directory: str = "",
):
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
):
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

    columns = ["Metric", "Score", "MaxMinScore"]
    metrics_dict = []
    for key, value in results_dict["metrics"].items():
        if isinstance(value, float) or isinstance(value, int):
            if key not in config["max_mins"]:
                max_min_score = value
            else:
                max_min_score = normalise_score(
                    value,
                    config["max_mins"]['_'.join(key.split("_")[:-1])][0],
                    config["max_mins"]['_'.join(key.split("_")[:-1])][1]
                )
            metrics_dict.append(
                {"Metric": key,
                 "Score": value,
                 "MaxMinScore": max_min_score,
                 }
                 )
        elif isinstance(value, dict):
            for k, v in value.items():
                if key not in config["max_mins"]:
                    max_min_score = v
                else:
                    max_min_score = normalise_score(
                        v,
                        config["max_mins"][key][0],
                        config["max_mins"][key][1]
                    )
                metrics_dict.append(
                    {
                        "Metric": f"{key}_{k}",
                        "Score": v,
                        "MaxMinScore": max_min_score
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
):
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

    # Write with headers
    file_exists = Path(output).exists()

    with open(output, "a" if file_exists else "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_row.keys(), delimiter='\t')
        if not file_exists:
            writer.writeheader()
        writer.writerow(summary_row)

    print(f"Summary line written to {output}")


def generate_qc_status(
    results_dict: dict,
    config: dict,
    sample_name: str,
    thresholds: Optional[Dict] = None,
    name: str = "RiboMetric_qc_status.json",
    output_directory: str = "",
):
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

    # Default thresholds (can be overridden)
    if thresholds is None:
        thresholds = {
            "periodicity_dominance": {"pass": 0.7, "warn": 0.5},
            "uniformity_entropy": {"pass": 0.7, "warn": 0.5},
            "prop_reads_CDS": {"pass": 0.7, "warn": 0.5},
            "read_length_distribution_IQR_metric": {"pass": 0.3, "warn": 0.2},
        }

    metrics = results_dict.get("metrics", {})
    qc_checks = []
    overall_status = "PASS"

    for metric_name, threshold_dict in thresholds.items():
        if metric_name not in metrics:
            continue

        metric_value = metrics[metric_name]

        # Handle global metrics
        if isinstance(metric_value, dict) and "global" in metric_value:
            value = metric_value["global"]
        elif isinstance(metric_value, (int, float)):
            value = metric_value
        else:
            continue

        # Determine status
        if value >= threshold_dict["pass"]:
            status = "PASS"
        elif value >= threshold_dict["warn"]:
            status = "WARNING"
            if overall_status == "PASS":
                overall_status = "WARNING"
        else:
            status = "FAIL"
            overall_status = "FAIL"

        qc_checks.append({
            "metric": metric_name,
            "value": value,
            "status": status,
            "threshold_pass": threshold_dict["pass"],
            "threshold_warn": threshold_dict["warn"],
        })

    qc_status = {
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

    with open(output, "w") as f:
        json.dump(qc_status, f, indent=2)

    print(f"QC status written to {output}")
    print(f"Overall QC Status: {overall_status}")


def _get_recommendation(status: str, checks: List[Dict]) -> str:
    """Generate a recommendation based on QC status"""
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
):
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
):
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
        "periodicity_dominance": "Proportion of reads in dominant reading frame",
        "uniformity_entropy": "Shannon entropy of metagene distribution",
        "read_length_distribution_IQR_metric": "Normalized IQR of read length distribution",
        "terminal_nucleotide_bias_distribution_5_prime_metric": "KL divergence for 5' nucleotide bias",
        "terminal_nucleotide_bias_distribution_3_prime_metric": "KL divergence for 3' nucleotide bias",
        "CDS_coverage_metric": "Proportion of CDS covered by reads",
        "prop_reads_CDS": "Proportion of reads mapping to CDS",
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
):
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
