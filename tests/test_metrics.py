

import pandas as pd
from RiboMetric.modules import (
    read_length_distribution,
    ligation_bias_distribution,
)

from RiboMetric.metrics import (
    read_length_distribution_metric as rld_metric,
    ligation_bias_distribution_metric as lbd_metric,
    read_frame_information_content as rfd_metric,
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
    assert round(read_length_metric, 3) == 0.333


def test_ligation_bias_distribution_metric():
    """
    Test the ligation bias distribution metric
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    categories = ["first_dinucleotide", "last_dinucleotide"]
    read_df[categories] = read_df[categories].astype("category")
    ligation_bias_dict = ligation_bias_distribution(read_df)
    sequence_background = {2:
                           {"5_prime_bg": {
                                'AA': 0.24390243902439024,
                                'AC': 0.0,
                                'AG': 0.0975609756097561,
                                'AT': 0.0,
                                'CA': 0.04878048780487805,
                                'CC': 0.04317073170731707,
                                'CG': 0.0,
                                'CT': 0.03,
                                'GA': 0.0,
                                'GC': 0.0,
                                'GG': 0.17073170731707318,
                                'GT': 0.12195121951219512,
                                'TA': 0.024390243902439025,
                                'TC': 0.12195121951219512,
                                'TG': 0.0,
                                'TT': 0.0975609756097561}}}
    ligation_bias_metric = lbd_metric(
        ligation_bias_dict, sequence_background[2]["5_prime_bg"]
        )

    assert round(ligation_bias_metric, 2) == 1.18


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

    weighted_score = triplet_periodicity_weighted_score(
        pre_scores,
    )
    assert round(weighted_score, 2) == 0.38
