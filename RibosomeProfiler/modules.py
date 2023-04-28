'''
This script contains the functions required to run individual modules of the RibosomeProfiler pipeline

'''

import pandas as pd
import numpy as np

def read_length_distribution(read_df: pd.DataFrame, config: dict) -> dict:
    '''
    Calculate the read length distribution for the full dataset

    Inputs:
        read_df: Dataframe containing the read information
        config: Dictionary containing the configuration information

    Outputs:
        dict: Dictionary containing the read length distribution
    '''
    read_lengths, read_counts = np.unique(read_df['read_length'], return_counts=True)
    return dict(zip(read_lengths, read_counts))
