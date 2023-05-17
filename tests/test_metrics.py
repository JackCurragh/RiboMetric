

import pandas as pd
from RiboMetric.modules import (
    read_length_distribution,
    ligation_bias_distribution,
    calculate_expected_dinucleotide_freqs,
)

from RiboMetric.metrics import (
    read_length_distribution_metric as rld_metric,
    ligation_bias_distribution_metric as lbd_metric,
    read_frame_distribution_information_content_metric as rfd_metric,
    triplet_periodicity_weighted_score,
    triplet_periodicity_best_read_length_score as tpbrl_metric,
    information_metric_cutoff,
)


def test_read_length_distribution_metric():
    """
    Test the read length distribution metric
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    read_length_dict = read_length_distribution(read_df)
    read_length_metric = rld_metric(read_length_dict)
    assert read_length_metric == 6.0


def test_ligation_bias_distribution_metric():
    """
    Test the ligation bias distribution metric
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    ligation_bias_dict = ligation_bias_distribution(read_df)
    expected_freqs = calculate_expected_dinucleotide_freqs(read_df)
    ligation_bias_metric = lbd_metric(ligation_bias_dict, expected_freqs)

    assert round(ligation_bias_metric, 2) == 1.17


def test_read_frame_distribution_metric():
    """
    Test the information content metric for read frame distribution
    """
    read_frame_dict = {
        27: {0: 100, 1: 5, 2: 5},
        28: {0: 1000, 1: 50000, 2: 50},
        29: {0: 100, 1: 2000, 2: 5},
        30: {0: 100000, 1: 50000, 2: 50000},
        31: {0: 100, 1: 50, 2: 500},
        32: {0: 1000, 1: 5000, 2: 5000},
    }
    pre_scores = rfd_metric(read_frame_dict)
    read_frame_metric = information_metric_cutoff(pre_scores)
    assert round(read_frame_metric[30], 2) == 0.23


def test_read_frame_distribution_metric_best_read_length():
    """
    Test the information content metric for read frame distribution
    using the best_read length score
    """
    read_frame_dict = {
        27: {0: 100, 1: 5, 2: 5},
        28: {0: 1000, 1: 50000, 2: 50},
        29: {0: 100, 1: 2000, 2: 5},
        30: {0: 100000, 1: 50000, 2: 50000},
        31: {0: 100, 1: 50, 2: 500},
        32: {0: 1000, 1: 5000, 2: 5000},
    }
    pre_scores = rfd_metric(read_frame_dict)
    read_frame_metric = information_metric_cutoff(pre_scores)
    read_frame_metric_best_read_length = tpbrl_metric(
        read_frame_metric
    )
    assert round(read_frame_metric_best_read_length, 2) == 0.95


def test_triplet_periodicity_weighted_score():
    """
    Test the triplet periodicity weighted score
    """
    read_frame_dict = {
        27: {0: 100, 1: 5, 2: 5},
        28: {0: 1000, 1: 50000, 2: 50},
        29: {0: 100, 1: 2000, 2: 5},
        30: {0: 100000, 1: 50000, 2: 50000},
        31: {0: 100, 1: 50, 2: 500},
        32: {0: 1000, 1: 5000, 2: 5000},
    }
    pre_scores = rfd_metric(read_frame_dict)
    read_frame_metric = information_metric_cutoff(pre_scores)

    weighted_score = triplet_periodicity_weighted_score(
        read_frame_metric,
        read_frame_dict,
    )
    assert round(weighted_score, 2) == 0.36
