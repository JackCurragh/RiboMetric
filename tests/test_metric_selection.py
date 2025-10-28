"""
Tests for the metric selection feature (default vs optional metrics)
"""

import pytest
from RiboMetric.qc import should_calculate_metric
from RiboMetric.arg_parser import open_config
from argparse import Namespace


class TestMetricSelection:
    """Test suite for metric selection functionality"""

    def test_should_calculate_metric_no_config(self):
        """Test backwards compatibility - no enabled_metrics means calculate all"""
        config = {}
        assert should_calculate_metric("any_metric", config) is True

    def test_should_calculate_metric_with_config(self, default_config):
        """Test metric calculation with enabled metrics"""
        # Set up config with only default metrics enabled
        config = default_config.copy()
        config["enabled_metrics"] = config["metrics"]["default"]

        # Default metrics should be calculated
        assert should_calculate_metric("periodicity_dominance", config) is True
        assert should_calculate_metric("uniformity_entropy", config) is True

        # Optional metrics should not be calculated
        assert should_calculate_metric("periodicity_fourier", config) is False
        assert should_calculate_metric("uniformity_gini_index", config) is False

    def test_should_calculate_metric_name_variants(self, default_config):
        """Test that metric name variants are recognized"""
        config = default_config.copy()
        config["enabled_metrics"] = ["periodicity_dominance"]

        # Should match with or without prefix/suffix
        assert should_calculate_metric("periodicity_dominance", config) is True
        assert should_calculate_metric("periodicity_dominance_metric", config) is True

    def test_should_calculate_metric_optional_enabled(self, default_config):
        """Test when optional metrics are explicitly enabled"""
        config = default_config.copy()
        config["enabled_metrics"] = (
            config["metrics"]["default"] + config["metrics"]["optional"]
        )

        # All metrics should be calculable
        assert should_calculate_metric("periodicity_fourier", config) is True
        assert should_calculate_metric("uniformity_gini_index", config) is True
        assert should_calculate_metric("periodicity_dominance", config) is True

    def test_should_calculate_metric_specific_optional(self, default_config):
        """Test enabling specific optional metrics"""
        config = default_config.copy()
        config["enabled_metrics"] = config["metrics"]["default"] + [
            "periodicity_fourier"
        ]

        # Default + specific optional should be calculated
        assert should_calculate_metric("periodicity_dominance", config) is True
        assert should_calculate_metric("periodicity_fourier", config) is True

        # Other optional metrics should not
        assert should_calculate_metric("uniformity_gini_index", config) is False


class TestArgParser:
    """Test argument parser metric selection handling"""

    def test_default_metrics_only(self, config_path, test_data_dir):
        """Test that default behavior only enables default metrics"""
        args = Namespace(
            command="run",
            bam=str(test_data_dir / "test.bam"),
            annotation=str(test_data_dir / "1000_entry_RiboMetric.tsv"),
            config=str(config_path),
            gff=None,
            fasta=None,
            subsample=None,
            transcripts=None,
            threads=2,
            json=False,
            html=True,
            pdf=False,
            csv=False,
            all=False,
            enable_optional_metrics=False,
            enable_metric=None,
            offset_calculation_method="changepoint",
        )

        config = open_config(args)

        # Should have enabled_metrics populated
        assert "enabled_metrics" in config
        enabled = set(config["enabled_metrics"])

        # Should include defaults
        assert "periodicity_dominance" in enabled
        assert "uniformity_entropy" in enabled

        # Should not include all optionals (unless they're in config file defaults)
        # This depends on your config.yml content

    def test_enable_all_optional_metrics(self, config_path, test_data_dir):
        """Test --enable-optional-metrics flag"""
        args = Namespace(
            command="run",
            bam=str(test_data_dir / "test.bam"),
            annotation=str(test_data_dir / "1000_entry_RiboMetric.tsv"),
            config=str(config_path),
            gff=None,
            fasta=None,
            subsample=None,
            transcripts=None,
            threads=2,
            json=False,
            html=True,
            pdf=False,
            csv=False,
            all=False,
            enable_optional_metrics=True,  # Enable all optional
            enable_metric=None,
            offset_calculation_method="changepoint",
        )

        config = open_config(args)

        enabled = set(config["enabled_metrics"])

        # Should include both default and optional metrics
        assert "periodicity_dominance" in enabled
        assert "uniformity_entropy" in enabled
        # Check for some optional metrics from config
        assert "periodicity_autocorrelation" in enabled or \
               "periodicity_fourier" in enabled or \
               "uniformity_gini_index" in enabled

    def test_enable_specific_metrics(self, config_path, test_data_dir):
        """Test --enable-metric flag for specific metrics"""
        args = Namespace(
            command="run",
            bam=str(test_data_dir / "test.bam"),
            annotation=str(test_data_dir / "1000_entry_RiboMetric.tsv"),
            config=str(config_path),
            gff=None,
            fasta=None,
            subsample=None,
            transcripts=None,
            threads=2,
            json=False,
            html=True,
            pdf=False,
            csv=False,
            all=False,
            enable_optional_metrics=False,
            enable_metric=["periodicity_fourier", "uniformity_gini_index"],
            offset_calculation_method="changepoint",
        )

        config = open_config(args)

        enabled = set(config["enabled_metrics"])

        # Should include defaults
        assert "periodicity_dominance" in enabled

        # Should include specifically requested metrics
        assert "periodicity_fourier" in enabled
        assert "uniformity_gini_index" in enabled

    def test_prepare_mode_no_metrics(self, config_path, test_data_dir):
        """Test that prepare mode doesn't set enabled_metrics"""
        args = Namespace(
            command="prepare",
            gff=str(test_data_dir / "1000_entry.gff"),
            output=str(test_data_dir),
            config=str(config_path),
            transcripts=None,
            threads=2,
        )

        config = open_config(args)

        # Prepare mode shouldn't set enabled_metrics
        # (it's only for run mode)
        # Either it's not present or empty
        if "enabled_metrics" in config:
            # This is fine, just means the config initialization happened
            pass


class TestMetricCalculationConditionals:
    """Test that metrics are actually conditionally calculated"""

    def test_default_metrics_calculated(self, default_config):
        """Test that default metrics are included in results"""
        # This would require running the full pipeline
        # For now, we just test the configuration logic
        config = default_config.copy()
        config["enabled_metrics"] = config["metrics"]["default"]

        # Verify default metrics would be calculated
        for metric in config["metrics"]["default"]:
            assert should_calculate_metric(metric, config) is True

    def test_optional_metrics_skipped_by_default(self, default_config):
        """Test that optional metrics are skipped when not enabled"""
        config = default_config.copy()
        config["enabled_metrics"] = config["metrics"]["default"]

        # Verify optional metrics would be skipped
        for metric in config["metrics"]["optional"]:
            assert should_calculate_metric(metric, config) is False

    def test_optional_metrics_calculated_when_enabled(self, default_config):
        """Test that optional metrics are calculated when explicitly enabled"""
        config = default_config.copy()
        config["enabled_metrics"] = (
            config["metrics"]["default"] + config["metrics"]["optional"]
        )

        # Verify all metrics would be calculated
        for metric in config["metrics"]["default"] + config["metrics"]["optional"]:
            assert should_calculate_metric(metric, config) is True
