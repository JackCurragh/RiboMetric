'''
This script contains the functions used to calculate individual metrics
for different aspects of ribosome profiling data. The functions are called
are called from qc and the input data comes from the output of their
respective modules

'''

import pandas as pd


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
