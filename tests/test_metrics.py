import pandas as pd
import pytest

from RiboMetric.metrics import (
    information_metric_cutoff,
    read_frame_information_weighted_score,
    read_length_iqr_fraction,
    terminal_nucleotide_bias_KL_divergence,
)
from RiboMetric.metrics import (
    read_frame_information_content as rfd_metric,
)
from RiboMetric.modules import (
    read_length_distribution,
    terminal_nucleotide_bias_distribution,
)


def test_read_length_distribution_metric():
    """
    Test the read length distribution metric
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[read_df_pre.index.repeat(read_df_pre["count"])].reset_index(drop=True)
    read_length_dict = read_length_distribution(read_df)
    # Reported as the spread itself now (lower is better), not 1 - spread.
    read_length_metric = read_length_iqr_fraction(read_length_dict)
    assert round(read_length_metric, 3) == 0.667


def test_terminal_nucleotide_bias_KL_divergence():
    """
    Test the ligation bias distribution metric
    """
    # Collapsed rows: each row stands for `count` reads and the module weights
    # by it. This used to expand the rows while keeping `count`, which only
    # gave the right answer because repeated read names switched the weights
    # off library-wide (O-037).
    read_df = pd.read_csv("tests/test_data/test.csv")
    categories = ["first_dinucleotide", "last_dinucleotide"]
    read_df[categories] = read_df[categories].astype("category")
    terminal_nucleotide_bias_dict = terminal_nucleotide_bias_distribution(read_df)
    sequence_background = {
        2: {
            "5_prime_bg": {
                "AA": 0.24390243902439024,
                "AC": 0.0,
                "AG": 0.0975609756097561,
                "AT": 0.0,
                "CA": 0.04878048780487805,
                "CC": 0.04317073170731707,
                "CG": 0.0,
                "CT": 0.03,
                "GA": 0.0,
                "GC": 0.0,
                "GG": 0.17073170731707318,
                "GT": 0.12195121951219512,
                "TA": 0.024390243902439025,
                "TC": 0.12195121951219512,
                "TG": 0.0,
                "TT": 0.0975609756097561,
            }
        }
    }
    # Reported in bits now (higher = more bias). The old 1/(1 + KL) goodness
    # value of 0.4581 corresponds to KL = 1.1829 bits.
    terminal_nucleotide_bias_metric = terminal_nucleotide_bias_KL_divergence(
        terminal_nucleotide_bias_dict, sequence_background[2]["5_prime_bg"]
    )

    assert terminal_nucleotide_bias_metric == pytest.approx(1.1829, rel=1e-3)
    assert 1 / (1 + terminal_nucleotide_bias_metric) == pytest.approx(0.4581, rel=1e-3)


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
    assert round(read_frame_metric[30], 2) == 0.05


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

    weighted_score = read_frame_information_weighted_score(
        pre_scores,
    )
    assert round(weighted_score, 2) == 0.23
