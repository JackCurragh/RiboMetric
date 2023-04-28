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

def ligation_bias_distribution(read_df: pd.DataFrame, config: dict) -> dict:
    '''
    Calculate the proportion of the occurence in the first 2 nucleotides of the reads to check for ligation bias

    Inputs:
        read_df: Dataframe containing the read information
        config: Dictionary containing the configuration information
    
    Outputs:
        read_start_df: Dictionary containing the distribution of the first two nucleotides in the reads
    '''
    two_sequence_dict = dict(read_df['sequence'].str.slice(stop=2).value_counts(normalize=True).sort_index())
    ligation_bias_dict = {k:v for k,v in two_sequence_dict.items() if not k.startswith("N")}
    ligation_bias_dict.update({k:v for k,v in two_sequence_dict.items() if k.startswith("N")})
    return ligation_bias_dict