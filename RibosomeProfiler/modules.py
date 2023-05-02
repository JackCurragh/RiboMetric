"""
This script contains the functions required to run individual modules
of the RibosomeProfiler pipeline

"""

import pandas as pd
import numpy as np


def read_length_distribution(read_df: pd.DataFrame) -> dict:
    """
    Calculate the read length distribution for the full dataset

    Inputs:
        read_df: Dataframe containing the read information

    Outputs:
        dict: Dictionary containing the read length distribution
    """
    read_lengths, read_counts = np.unique(read_df["read_length"],
                                          return_counts=True)
    return dict(zip(read_lengths, read_counts))


def ligation_bias_distribution(
    read_df: pd.DataFrame, num_bases: int = 2, five_prime: bool = True
) -> dict:
    """
    Calculate the proportion of the occurence in the first or last n
    nucleotides of the reads to check for ligation bias

    Inputs:
        read_df: Dataframe containing the read information
        num_bases: Number of bases to be read (Default = 2)
        five_prime: Start at 5' end (True) or 3' end (False) of read
        (Default = True)

    Outputs:
        read_start_df: Dictionary containing the distribution of the
        first two nucleotides in the reads
    """
    if five_prime:
        two_sequence_dict = dict(
            read_df["sequence"]
            .str.slice(stop=num_bases)
            .value_counts(normalize=True)
            .sort_index()
        )
    else:
        two_sequence_dict = dict(
            read_df["sequence"]
            .str.slice(start=-num_bases)
            .value_counts(normalize=True)
            .sort_index()
        )
    ligation_bias_dict = {
        k: v for k, v in two_sequence_dict.items() if "N" not in k
        }
    ligation_bias_dict.update(
        {k: v for k, v in two_sequence_dict.items() if "N" in k}
        )
    return ligation_bias_dict


def nucleotide_composition(read_df: pd.DataFrame) -> dict:
    """
    Calculate the nucleotide composition

    Inputs:
        read_df: Dataframe containing the read information

    Outputs:
        dict: Dictionary containing the nucleotide distribution for every read position.
    """
    readlen = read_df['sequence'].str.len().max()
    nucleotide_composition_dict = {"A": list(), "C": list(), "G": list(), "T": list()}
    base_nts = pd.Series([0,0,0,0], index =['A','C','G','T'])
    for i in range(readlen):
        nucleotide_counts = read_df.sequence.str.slice(i, i + 1).value_counts()
        nucleotide_counts.drop("", errors="ignore", inplace=True)
        nucleotide_counts = base_nts.add(nucleotide_counts,fill_value=0)
        nucleotide_sum = nucleotide_counts.sum()
        nucleotide_composition_dict["A"].append(nucleotide_counts["A"] / nucleotide_sum)
        nucleotide_composition_dict["C"].append(nucleotide_counts["C"] / nucleotide_sum)
        nucleotide_composition_dict["G"].append(nucleotide_counts["G"] / nucleotide_sum)
        nucleotide_composition_dict["T"].append(nucleotide_counts["T"] / nucleotide_sum)
    return nucleotide_composition_dict
