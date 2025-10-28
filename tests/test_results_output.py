"""
Tests for results_output.py - both legacy and improved functions
"""

import pytest
import json
import csv
import pandas as pd
from pathlib import Path
from RiboMetric.results_output import (
    generate_json,
    generate_csv,
    normalise_score,
    generate_summary_tsv,
    generate_qc_status,
    generate_comparison_ready_csv,
    generate_metrics_table_csv,
    generate_all_outputs,
)


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def sample_results_dict():
    """Sample results dictionary for testing"""
    return {
        "mode": "test",
        "read_length_distribution": {28: 100, 29: 200, 30: 150},
        "metrics": {
            "periodicity_dominance": {"global": 0.85, 28: 0.82, 29: 0.88, 30: 0.84},
            "uniformity_entropy": {"global": 0.72, 28: 0.70, 29: 0.75, 30: 0.71},
            "prop_reads_CDS": 0.78,
            "read_length_distribution_IQR_metric": 0.25,
            "terminal_nucleotide_bias_distribution_5_prime_metric": 0.15,
        },
    }


@pytest.fixture
def sample_config():
    """Sample config for testing"""
    return {
        "max_mins": {
            # The legacy CSV uses '_'.join(key.split("_")[:-1]) for simple metrics
            # So "prop_reads_CDS" becomes "prop_reads"
            "prop_reads": [0, 1],
            "read_length_distribution_IQR": [0, 1],
            "terminal_nucleotide_bias_distribution_5_prime": [0, 0.5],
            "periodicity_dominance": [0, 1],
            "uniformity_entropy": [0, 1],
        },
        "metrics": {
            "default": ["periodicity_dominance", "uniformity_entropy"],
            "optional": [],
        },
    }


@pytest.fixture
def qc_thresholds():
    """Sample QC thresholds for testing"""
    return {
        "periodicity_dominance": {"pass": 0.7, "warn": 0.5},
        "uniformity_entropy": {"pass": 0.7, "warn": 0.5},
        "prop_reads_CDS": {"pass": 0.7, "warn": 0.5},
    }


# =============================================================================
# Legacy Function Tests
# =============================================================================

class TestLegacyFunctions:
    """Test backwards-compatible legacy functions"""

    def test_normalise_score_within_range(self):
        """Test score normalization within valid range"""
        score = normalise_score(0.5, 0, 1)
        assert score == 0.5

    def test_normalise_score_above_max(self):
        """Test score normalization above maximum"""
        score = normalise_score(1.5, 0, 1)
        assert score == 1

    def test_normalise_score_below_min(self):
        """Test score normalization below minimum"""
        score = normalise_score(-0.5, 0, 1)
        assert score == 0

    def test_generate_json(self, sample_results_dict, sample_config, tmp_path):
        """Test JSON generation"""
        output_file = tmp_path / "test.json"
        generate_json(
            sample_results_dict,
            sample_config,
            "test",
            str(tmp_path),
        )

        assert output_file.exists()
        with open(output_file) as f:
            data = json.load(f)
            assert "results" in data
            assert "config" in data
            assert data["results"]["mode"] == "test"

    def test_generate_csv(self, sample_results_dict, sample_config, tmp_path):
        """Test CSV generation"""
        output_file = tmp_path / "test.csv"
        generate_csv(
            sample_results_dict,
            sample_config,
            "test",
            str(tmp_path),
        )

        assert output_file.exists()
        df = pd.read_csv(output_file)
        assert "Metric" in df.columns
        assert "Score" in df.columns
        assert "MaxMinScore" in df.columns
        assert len(df) > 0


# =============================================================================
# Improved Function Tests
# =============================================================================

class TestSummaryTSV:
    """Test generate_summary_tsv function"""

    def test_creates_file_with_headers(self, sample_results_dict, sample_config, tmp_path):
        """Test that TSV file is created with proper headers"""
        output_file = tmp_path / "summary.tsv"
        generate_summary_tsv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "summary.tsv",
            str(tmp_path),
        )

        assert output_file.exists()
        df = pd.read_csv(output_file, sep="\t")
        assert "sample" in df.columns
        assert "timestamp" in df.columns
        assert "mode" in df.columns
        assert "total_reads" in df.columns

    def test_single_row_per_sample(self, sample_results_dict, sample_config, tmp_path):
        """Test that one row is written per sample"""
        output_file = tmp_path / "summary.tsv"
        generate_summary_tsv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "summary.tsv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 1
        assert df["sample"].iloc[0] == "Sample1"

    def test_appends_to_existing_file(self, sample_results_dict, sample_config, tmp_path):
        """Test that subsequent samples are appended"""
        output_file = tmp_path / "summary.tsv"

        # First sample
        generate_summary_tsv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "summary.tsv",
            str(tmp_path),
        )

        # Second sample
        generate_summary_tsv(
            sample_results_dict,
            sample_config,
            "Sample2",
            "summary.tsv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 2
        assert "Sample1" in df["sample"].values
        assert "Sample2" in df["sample"].values

    def test_includes_global_metrics(self, sample_results_dict, sample_config, tmp_path):
        """Test that global metrics are included"""
        output_file = tmp_path / "summary.tsv"
        generate_summary_tsv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "summary.tsv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file, sep="\t")
        assert "periodicity_dominance" in df.columns
        assert df["periodicity_dominance"].iloc[0] == 0.85

    def test_excludes_per_read_length_metrics(self, sample_results_dict, sample_config, tmp_path):
        """Test that per-read-length metrics are not added as separate columns"""
        output_file = tmp_path / "summary.tsv"
        generate_summary_tsv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "summary.tsv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file, sep="\t")
        # Should have global metric but not "periodicity_dominance_28" etc
        assert "periodicity_dominance" in df.columns
        assert "periodicity_dominance_28" not in df.columns


class TestQCStatus:
    """Test generate_qc_status function"""

    def test_creates_json_file(self, sample_results_dict, sample_config, qc_thresholds, tmp_path):
        """Test that QC status JSON is created"""
        output_file = tmp_path / "qc_status.json"
        generate_qc_status(
            sample_results_dict,
            sample_config,
            "Sample1",
            qc_thresholds,
            "qc_status.json",
            str(tmp_path),
        )

        assert output_file.exists()
        with open(output_file) as f:
            data = json.load(f)
            assert "sample" in data
            assert "overall_status" in data
            assert "checks" in data
            assert "summary" in data
            assert "recommendation" in data

    def test_pass_status_for_good_metrics(self, sample_results_dict, sample_config, qc_thresholds, tmp_path):
        """Test that PASS status is assigned when metrics exceed thresholds"""
        output_file = tmp_path / "qc_status.json"
        generate_qc_status(
            sample_results_dict,
            sample_config,
            "Sample1",
            qc_thresholds,
            "qc_status.json",
            str(tmp_path),
        )

        with open(output_file) as f:
            data = json.load(f)
            assert data["overall_status"] == "PASS"
            # Check individual metrics passed
            periodicity_check = [c for c in data["checks"] if c["metric"] == "periodicity_dominance"][0]
            assert periodicity_check["status"] == "PASS"

    def test_warning_status(self, sample_results_dict, sample_config, qc_thresholds, tmp_path):
        """Test WARNING status when metrics are between warn and pass"""
        # Modify results to trigger warning
        sample_results_dict["metrics"]["periodicity_dominance"]["global"] = 0.6

        output_file = tmp_path / "qc_status.json"
        generate_qc_status(
            sample_results_dict,
            sample_config,
            "Sample1",
            qc_thresholds,
            "qc_status.json",
            str(tmp_path),
        )

        with open(output_file) as f:
            data = json.load(f)
            periodicity_check = [c for c in data["checks"] if c["metric"] == "periodicity_dominance"][0]
            assert periodicity_check["status"] == "WARNING"

    def test_fail_status(self, sample_results_dict, sample_config, qc_thresholds, tmp_path):
        """Test FAIL status when metrics are below warn threshold"""
        # Modify results to trigger failure
        sample_results_dict["metrics"]["periodicity_dominance"]["global"] = 0.3

        output_file = tmp_path / "qc_status.json"
        generate_qc_status(
            sample_results_dict,
            sample_config,
            "Sample1",
            qc_thresholds,
            "qc_status.json",
            str(tmp_path),
        )

        with open(output_file) as f:
            data = json.load(f)
            assert data["overall_status"] == "FAIL"
            periodicity_check = [c for c in data["checks"] if c["metric"] == "periodicity_dominance"][0]
            assert periodicity_check["status"] == "FAIL"

    def test_summary_counts(self, sample_results_dict, sample_config, qc_thresholds, tmp_path):
        """Test that summary counts are correct"""
        output_file = tmp_path / "qc_status.json"
        generate_qc_status(
            sample_results_dict,
            sample_config,
            "Sample1",
            qc_thresholds,
            "qc_status.json",
            str(tmp_path),
        )

        with open(output_file) as f:
            data = json.load(f)
            summary = data["summary"]
            assert summary["total_checks"] == len(data["checks"])
            assert summary["passed"] + summary["warnings"] + summary["failed"] == summary["total_checks"]

    def test_default_thresholds(self, sample_results_dict, sample_config, tmp_path):
        """Test that default thresholds are used when none provided"""
        output_file = tmp_path / "qc_status.json"
        generate_qc_status(
            sample_results_dict,
            sample_config,
            "Sample1",
            None,  # No thresholds provided
            "qc_status.json",
            str(tmp_path),
        )

        assert output_file.exists()
        with open(output_file) as f:
            data = json.load(f)
            assert "checks" in data
            assert len(data["checks"]) > 0


class TestComparisonCSV:
    """Test generate_comparison_ready_csv function"""

    def test_creates_csv_file(self, sample_results_dict, sample_config, tmp_path):
        """Test that comparison CSV is created"""
        output_file = tmp_path / "comparison.csv"
        generate_comparison_ready_csv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "comparison.csv",
            str(tmp_path),
        )

        assert output_file.exists()
        df = pd.read_csv(output_file)
        assert "sample" in df.columns
        assert len(df) == 1

    def test_wide_format_with_global_metrics(self, sample_results_dict, sample_config, tmp_path):
        """Test that metrics are in wide format with global suffix"""
        output_file = tmp_path / "comparison.csv"
        generate_comparison_ready_csv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "comparison.csv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file)
        assert "periodicity_dominance_global" in df.columns
        assert df["periodicity_dominance_global"].iloc[0] == 0.85

    def test_includes_key_read_lengths(self, sample_results_dict, sample_config, tmp_path):
        """Test that key read lengths (28-32) are included"""
        output_file = tmp_path / "comparison.csv"
        generate_comparison_ready_csv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "comparison.csv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file)
        assert "periodicity_dominance_rl28" in df.columns
        assert "periodicity_dominance_rl29" in df.columns
        assert "periodicity_dominance_rl30" in df.columns

    def test_appends_multiple_samples(self, sample_results_dict, sample_config, tmp_path):
        """Test that multiple samples are appended correctly"""
        output_file = tmp_path / "comparison.csv"

        # First sample
        generate_comparison_ready_csv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "comparison.csv",
            str(tmp_path),
        )

        # Second sample
        generate_comparison_ready_csv(
            sample_results_dict,
            sample_config,
            "Sample2",
            "comparison.csv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file)
        assert len(df) == 2
        assert "Sample1" in df["sample"].values
        assert "Sample2" in df["sample"].values


class TestMetricsTableCSV:
    """Test generate_metrics_table_csv function"""

    def test_creates_csv_file(self, sample_results_dict, sample_config, tmp_path):
        """Test that metrics table CSV is created"""
        output_file = tmp_path / "metrics_table.csv"
        generate_metrics_table_csv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "metrics_table.csv",
            str(tmp_path),
        )

        assert output_file.exists()
        df = pd.read_csv(output_file)
        assert "sample" in df.columns
        assert "metric" in df.columns
        assert "read_length_or_region" in df.columns
        assert "value" in df.columns
        assert "description" in df.columns

    def test_long_format_structure(self, sample_results_dict, sample_config, tmp_path):
        """Test that data is in long format (one row per metric/read_length)"""
        output_file = tmp_path / "metrics_table.csv"
        generate_metrics_table_csv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "metrics_table.csv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file)
        # Should have rows for each metric and read length combination
        periodicity_rows = df[df["metric"] == "periodicity_dominance"]
        assert len(periodicity_rows) == 4  # global + 28, 29, 30

    def test_includes_descriptions(self, sample_results_dict, sample_config, tmp_path):
        """Test that metric descriptions are included"""
        output_file = tmp_path / "metrics_table.csv"
        generate_metrics_table_csv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "metrics_table.csv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file)
        periodicity_row = df[df["metric"] == "periodicity_dominance"].iloc[0]
        assert periodicity_row["description"] == "Proportion of reads in dominant reading frame"

    def test_global_metrics_labeled_correctly(self, sample_results_dict, sample_config, tmp_path):
        """Test that global metrics are labeled as 'global'"""
        output_file = tmp_path / "metrics_table.csv"
        generate_metrics_table_csv(
            sample_results_dict,
            sample_config,
            "Sample1",
            "metrics_table.csv",
            str(tmp_path),
        )

        df = pd.read_csv(output_file)
        prop_reads_row = df[df["metric"] == "prop_reads_CDS"].iloc[0]
        assert prop_reads_row["read_length_or_region"] == "global"
        assert prop_reads_row["value"] == 0.78


class TestGenerateAllOutputs:
    """Test generate_all_outputs convenience function"""

    def test_creates_all_files(self, sample_results_dict, sample_config, qc_thresholds, tmp_path):
        """Test that all output files are created"""
        generate_all_outputs(
            sample_results_dict,
            sample_config,
            "Sample1",
            str(tmp_path),
            qc_thresholds,
        )

        # Check all expected files exist
        assert (tmp_path / "Sample1_summary.tsv").exists()
        assert (tmp_path / "Sample1_metrics_table.csv").exists()
        assert (tmp_path / "Sample1_qc_status.json").exists()
        assert (tmp_path / "Sample1_comparison.csv").exists()

    def test_files_have_correct_content(self, sample_results_dict, sample_config, qc_thresholds, tmp_path):
        """Test that generated files have correct content"""
        generate_all_outputs(
            sample_results_dict,
            sample_config,
            "Sample1",
            str(tmp_path),
            qc_thresholds,
        )

        # Verify TSV
        tsv_df = pd.read_csv(tmp_path / "Sample1_summary.tsv", sep="\t")
        assert tsv_df["sample"].iloc[0] == "Sample1"

        # Verify QC JSON
        with open(tmp_path / "Sample1_qc_status.json") as f:
            qc_data = json.load(f)
            assert qc_data["sample"] == "Sample1"

        # Verify comparison CSV
        comp_df = pd.read_csv(tmp_path / "Sample1_comparison.csv")
        assert comp_df["sample"].iloc[0] == "Sample1"

        # Verify metrics table
        metrics_df = pd.read_csv(tmp_path / "Sample1_metrics_table.csv")
        assert "Sample1" in metrics_df["sample"].values
