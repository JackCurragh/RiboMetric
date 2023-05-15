'''
This script contains the functions used to calculate individual metrics
for different aspects of ribosome profiling data. The functions are called
are called from qc and the input data comes from the output of their
respective modules

'''

import pandas as pd
import math


def read_length_distribution_metric(
        rld_dict: dict,
        ) -> pd.DataFrame:
    """
    Calculate the read length distribution metric from the output of
    the read_length_distribution module.

    This metric is the IQR of the read length distribution and is
    calculated as the difference between the 75th and 25th percentile

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module

    Outputs:
        rld_df: Dataframe containing the read length distribution metric
    """
    rld_df = pd.DataFrame.from_dict(rld_dict, orient="index")
    rld_df = rld_df.reset_index()
    rld_df.columns = ["read_length", "read_count"]

    Q3 = rld_df["read_length"].quantile(0.75)
    Q1 = rld_df["read_length"].quantile(0.25)

    return Q3 - Q1


def ligation_bias_distribution_metric(
        observed_freq: dict,
        expected_freq: dict,
        ) -> pd.DataFrame:
    """
    Calculate the ligation bias metric from the output of
    the ligation_bias_distribution module.

    This metric is the K-L divergence of the ligation bias distribution
    of the observed frequencies from the expected frequencies. The
    expected frequencies are calculated from the nucleotide composition
    of the genome.

    Inputs:
        observed_freq: Dictionary containing the output of the
                ligation_bias_distribution module
        expected_freq: Dictionary containing the expected frequencies

    Outputs:
        lbd_df: Dataframe containing the ligation bias metric in bits
    """
    kl_divergence = 0.0

    for dinucleotide, observed_prob in observed_freq.items():
        expected_prob = expected_freq[dinucleotide]
        kl_divergence += observed_prob * math.log2(
                                            observed_prob / expected_prob
                                            )

    return kl_divergence

