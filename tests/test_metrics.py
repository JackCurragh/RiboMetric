

import pandas as pd
from RiboMetric.modules import (
    read_length_distribution,
    ligation_bias_distribution,
    calculate_expected_dinucleotide_freqs,
)

from RiboMetric.metrics import (
    read_length_distribution_metric as rld_metric,
    ligation_bias_distribution_metric as lbd_metric,
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

    assert round(ligation_bias_metric, 2) == 1.26
