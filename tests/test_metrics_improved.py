"""
Improved and expanded tests for RiboMetric metrics
"""

import pytest
import numpy as np
from RiboMetric.metrics import (
    read_length_distribution_IQR_normalised_metric,
    read_length_distribution_coefficient_of_variation_metric,
    read_length_distribution_max_prop_metric,
    read_length_distribution_bimodality,
    read_length_distribution_normality_metric,
    periodicity_dominance,
    periodicity_autocorrelation,
    uniformity_entropy,
    uniformity_gini_index,
    uniformity_theil_index,
    fourier_transform,
    read_frame_information_content,
    information_metric_cutoff,
    read_frame_information_weighted_score,
    terminal_nucleotide_bias_KL_metric,
    terminal_nucleotide_bias_max_absolute_metric,
    cds_coverage_metric,
    region_region_ratio_metric,
    proportion_of_reads_in_region,
)


class TestReadLengthMetrics:
    """Tests for read length distribution metrics"""

    def test_iqr_metric(self, sample_read_length_dict):
        """Test IQR normalized metric"""
        metric = read_length_distribution_IQR_normalised_metric(
            sample_read_length_dict
        )
        assert 0 <= metric <= 1
        assert isinstance(metric, (int, float))

    def test_iqr_metric_single_length(self):
        """Test IQR metric with single read length"""
        single_length = {28: 1000}
        metric = read_length_distribution_IQR_normalised_metric(single_length)
        # Should handle edge case gracefully
        assert isinstance(metric, (int, float))

    def test_coefficient_of_variation(self, sample_read_length_dict):
        """Test coefficient of variation metric"""
        metric = read_length_distribution_coefficient_of_variation_metric(
            sample_read_length_dict
        )
        assert 0 < metric <= 1
        assert isinstance(metric, float)

    def test_max_prop_metric(self, sample_read_length_dict):
        """Test maximum proportion metric"""
        metric = read_length_distribution_max_prop_metric(
            sample_read_length_dict, num_top_readlens=1
        )
        assert 0 < metric <= 1
        # Should be the proportion of the most frequent read length
        max_count = max(sample_read_length_dict.values())
        total_count = sum(sample_read_length_dict.values())
        expected = max_count / total_count
        assert abs(metric - expected) < 1e-10

    def test_max_prop_metric_top3(self, sample_read_length_dict):
        """Test maximum proportion metric with top 3 read lengths"""
        metric = read_length_distribution_max_prop_metric(
            sample_read_length_dict, num_top_readlens=3
        )
        assert 0 < metric <= 1
        # Should be higher than single length
        metric_single = read_length_distribution_max_prop_metric(
            sample_read_length_dict, num_top_readlens=1
        )
        assert metric >= metric_single

    def test_bimodality(self, sample_read_length_dict):
        """Test bimodality coefficient"""
        metric = read_length_distribution_bimodality(sample_read_length_dict)
        assert isinstance(metric, float)
        assert 0 < metric <= 1

    def test_normality_metric(self, sample_read_length_dict):
        """Test normality metric (p-value from normality test)"""
        metric = read_length_distribution_normality_metric(
            sample_read_length_dict
        )
        assert 0 <= metric <= 1
        assert isinstance(metric, float)


class TestPeriodicityMetrics:
    """Tests for periodicity metrics"""

    def test_dominance_strong_periodicity(self):
        """Test dominance metric with strong periodicity"""
        strong_periodicity = {
            28: {0: 1000, 1: 10, 2: 10},
            29: {0: 900, 1: 5, 2: 5},
            30: {0: 1100, 1: 15, 2: 15},
        }
        metric = periodicity_dominance(strong_periodicity)

        # Global metric should be high (strong dominance)
        assert metric["global"] > 0.8
        # Per-read-length metrics should also be high
        assert all(metric[rl] > 0.8 for rl in [28, 29, 30])

    def test_dominance_weak_periodicity(self):
        """Test dominance metric with weak periodicity"""
        weak_periodicity = {
            28: {0: 100, 1: 90, 2: 80},
            29: {0: 110, 1: 100, 2: 95},
            30: {0: 105, 1: 100, 2: 98},
        }
        metric = periodicity_dominance(weak_periodicity)

        # Global metric should be low (weak dominance)
        assert metric["global"] < 0.5

    def test_dominance_no_reads(self):
        """Test dominance metric with zero reads"""
        no_reads = {
            28: {0: 0, 1: 0, 2: 0},
        }
        metric = periodicity_dominance(no_reads)
        # Should handle gracefully
        assert metric[28] == 0
        assert metric["global"] == 0

    def test_information_content(self, sample_read_frame_dict):
        """Test information content metric"""
        metric = read_frame_information_content(sample_read_frame_dict)

        # Should return scores and counts for each read length
        for read_len in sample_read_frame_dict.keys():
            assert read_len in metric
            score, count = metric[read_len]
            assert 0 <= score <= 1  # Information content normalized
            assert count > 0

    def test_information_cutoff(self, sample_read_frame_dict):
        """Test information metric with cutoff"""
        frame_info = read_frame_information_content(sample_read_frame_dict)
        metric = information_metric_cutoff(frame_info, min_count_threshold=0.05)

        # Should have global score
        assert "global" in metric
        assert 0 <= metric["global"] <= 1

    def test_weighted_score(self, sample_read_frame_dict):
        """Test weighted information score"""
        frame_info = read_frame_information_content(sample_read_frame_dict)
        metric = read_frame_information_weighted_score(frame_info)

        assert 0 <= metric <= 1
        assert isinstance(metric, float)

    def test_autocorrelation_periodicity(self, sample_metagene_profile):
        """Test autocorrelation-based periodicity"""
        metric = periodicity_autocorrelation(sample_metagene_profile, lag=3)

        # Should have global and per-read-length scores
        assert "global" in metric
        for read_len in sample_metagene_profile["start"].keys():
            assert read_len in metric

    def test_fourier_transform(self, sample_metagene_profile):
        """Test Fourier transform periodicity"""
        metric = fourier_transform(sample_metagene_profile)

        # Should have global and per-read-length scores
        assert "global" in metric
        assert 0 <= metric["global"] <= 1


class TestUniformityMetrics:
    """Tests for uniformity metrics"""

    def test_entropy_uniform(self):
        """Test entropy with uniform distribution"""
        # Create perfectly uniform metagene profile
        uniform_profile = {"start": {}}
        for read_len in [28, 29, 30]:
            uniform_profile["start"][read_len] = {i: 100 for i in range(90)}

        metric = uniformity_entropy(uniform_profile)

        # High uniformity should give high entropy
        assert metric["global"] > 0.9
        assert all(metric[rl] > 0.9 for rl in [28, 29, 30])

    def test_entropy_non_uniform(self):
        """Test entropy with non-uniform distribution"""
        # Create peaked distribution
        peaked_profile = {"start": {}}
        for read_len in [28, 29, 30]:
            peaked_profile["start"][read_len] = {
                i: (1000 if i == 0 else 1) for i in range(90)
            }

        metric = uniformity_entropy(peaked_profile)

        # Low uniformity should give low entropy
        assert metric["global"] < 0.5

    def test_gini_index(self, sample_metagene_profile):
        """Test Gini index uniformity metric"""
        metric = uniformity_gini_index(sample_metagene_profile)

        # Should have global score
        assert "global" in metric
        # Gini coefficient between 0 and 1
        assert 0 <= metric["global"] <= 1

    def test_theil_index(self, sample_metagene_profile):
        """Test Theil index uniformity metric"""
        metric = uniformity_theil_index(sample_metagene_profile)

        # Should have global score
        assert "global" in metric
        assert 0 <= metric["global"] <= 1


class TestTerminalNucleotideBias:
    """Tests for terminal nucleotide bias metrics"""

    def test_kl_divergence_no_bias(self, sample_sequence_background):
        """Test KL divergence with no bias (observed = expected)"""
        observed = {
            "five_prime": sample_sequence_background["5_prime_bg"].copy(),
            "three_prime": sample_sequence_background["3_prime_bg"].copy()
        }

        metric_5 = terminal_nucleotide_bias_KL_metric(
            observed, sample_sequence_background["5_prime_bg"], "five_prime"
        )
        metric_3 = terminal_nucleotide_bias_KL_metric(
            observed, sample_sequence_background["3_prime_bg"], "three_prime"
        )

        # No divergence from expected
        assert metric_5 == pytest.approx(1.0)
        assert metric_3 == pytest.approx(1.0)

    def test_kl_divergence_with_bias(self, sample_sequence_background):
        """Test KL divergence with bias"""
        # Create biased distribution
        observed = {
            "five_prime": sample_sequence_background["5_prime_bg"].copy(),
            "three_prime": sample_sequence_background["3_prime_bg"].copy()
        }
        # Heavily bias towards AA
        observed["five_prime"]["AA"] = 0.5
        # Reduce others proportionally
        remaining = 0.5 / 15
        for nt in observed["five_prime"]:
            if nt != "AA":
                observed["five_prime"][nt] = remaining

        metric = terminal_nucleotide_bias_KL_metric(
            observed, sample_sequence_background["5_prime_bg"], "five_prime"
        )

        # Should have significant divergence
        assert metric < 0.9

    def test_max_absolute_metric(self, sample_sequence_background):
        """Test maximum absolute deviation metric"""
        observed = {
            "five_prime": sample_sequence_background["5_prime_bg"].copy(),
            "three_prime": sample_sequence_background["3_prime_bg"].copy()
        }
        # Bias one dinucleotide
        observed["five_prime"]["AA"] = 0.5
        remaining = 0.5 / 15
        for nt in observed["five_prime"]:
            if nt != "AA":
                observed["five_prime"][nt] = remaining

        metric = terminal_nucleotide_bias_max_absolute_metric(
            observed, sample_sequence_background["5_prime_bg"], "five_prime"
        )

        # Max deviation should be approximately 0.5 - 1/16
        expected_max = abs(0.5 - 1.0/16)
        expected_score = 1 - expected_max
        assert abs(metric - expected_score) < 0.01


class TestRegionalMetrics:
    """Tests for regional distribution metrics"""

    def test_region_ratio(self):
        """Test region-to-region ratio metric"""
        mRNA_dist = {
            28: {"leader": 100, "CDS": 1000, "trailer": 50},
            29: {"leader": 90, "CDS": 900, "trailer": 45},
            30: {"leader": 110, "CDS": 1100, "trailer": 55},
        }

        ratio = region_region_ratio_metric(
            mRNA_dist, region1="CDS", region2="leader"
        )

        # Should have global ratio
        assert "global" in ratio
        # CDS should be ~10x leader
        assert 9 < ratio["global"] < 11

    def test_region_ratio_zero_denominator(self):
        """Test region ratio with zero denominator"""
        mRNA_dist = {
            28: {"leader": 0, "CDS": 1000, "trailer": 50},
        }

        ratio = region_region_ratio_metric(
            mRNA_dist, region1="CDS", region2="leader"
        )

        # Should handle zero gracefully
        assert ratio["global"] == 0 or ratio[28] == 0

    def test_proportion_in_region(self):
        """Test proportion of reads in region"""
        mRNA_dist = {
            28: {"five_leader": 100, "CDS": 800, "three_trailer": 100},
            29: {"five_leader": 90, "CDS": 810, "three_trailer": 100},
        }

        prop = proportion_of_reads_in_region(mRNA_dist, region="CDS")

        # Should have global proportion
        assert "global" in prop
        # CDS should be ~80%
        assert 0.79 < prop["global"] < 0.81
