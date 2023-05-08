"""
This script contains the functions required to run individual modules
of the RibosomeProfiler pipeline

"""

import pandas as pd
import numpy as np
from xhtml2pdf import pisa


def read_df_to_cds_read_df(
    a_site_df: pd.DataFrame, annotation_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Convert the a_site_df to a cds_read_df by removing reads that do not
    map to the CDS

    Inputs:
        a_site_df: Dataframe containing the read information
        cds_df: Dataframe containing the coordinates of the CDS per tx

    Outputs:
        cds_read_df: Dataframe containing the read information for reads
                    that map to the CDS
    """
    cds_read_df = pd.DataFrame()
    for tx in annotation_df["transcript_id"]:
        tx_df = a_site_df[a_site_df["reference_name"].str.contains(str(tx))]
        idx = annotation_df[annotation_df["transcript_id"] == tx].index[0]
        tx_df = tx_df[
            tx_df["a_site"].between(
                annotation_df.loc[idx, "cds_start"],
                annotation_df.loc[idx, "cds_end"]
            )
        ]
        cds_read_df = pd.concat([cds_read_df, tx_df])

    return cds_read_df


def a_site_calculation(read_df: pd.DataFrame, offset=15) -> pd.DataFrame:
    """
    Adds a column to the read_df containing the A-site for the reads

    Inputs:
        read_df: Dataframe containing the read information
        offset: Offset from the start of the read to the A-site (Default = 15)

    Outputs:
        asite_df: Dataframe containing the read information with an added
                    column for the A-site
    """
    a_site_df = read_df.assign(a_site=read_df.reference_start.add(offset))
    return a_site_df


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
        sequence_dict = dict(
            read_df["sequence"]
            .str.slice(stop=num_bases)
            .value_counts(normalize=True)
            .sort_index()
        )
    else:
        sequence_dict = dict(
            read_df["sequence"]
            .str.slice(start=-num_bases)
            .value_counts(normalize=True)
            .sort_index()
        )
    ligation_bias_dict = {k: v for k, v in sequence_dict.items()
                          if "N" not in k}
    ligation_bias_dict.update({k: v for k, v in sequence_dict.items()
                               if "N" in k})
    return ligation_bias_dict


def nucleotide_composition(
    read_df: pd.DataFrame, nucleotides=["A", "C", "G", "T"]
) -> dict:
    """
    Calculate the nucleotide composition

    Inputs:
        read_df: Dataframe containing the read information

    Outputs:
        dict: Dictionary containing the nucleotide distribution for every
            read position.
    """
    readlen = read_df["sequence"].str.len().max()
    nucleotide_composition_dict = {nt: [] for nt in nucleotides}
    base_nts = pd.Series([0, 0, 0, 0], index=nucleotides)
    for i in range(readlen):
        nucleotide_counts = read_df.sequence.str.slice(i, i + 1).value_counts()
        nucleotide_counts.drop("", errors="ignore", inplace=True)
        nucleotide_counts = base_nts.add(nucleotide_counts, fill_value=0)
        nucleotide_sum = nucleotide_counts.sum()
        for nt in nucleotides:
            nt_proportion = nucleotide_counts[nt] / nucleotide_sum
            nucleotide_composition_dict[nt].append(nt_proportion)
    return nucleotide_composition_dict


def read_frame_distribution(a_site_df: pd.DataFrame) -> dict:
    """
    Calculate the distribution of the reading frame over the dataset

    Inputs:
        a_site_df: Dataframe ontaining the read information with an added 
        column for the a-site location

    Outputs:
        read_frame_dict: Nested dictionary containing counts for every reading
        frame at the different read lengths
    """
    frame_df = (
        a_site_df.assign(read_frame=a_site_df.a_site.mod(3))
        .groupby(["read_length", "read_frame"])
        .size()
    )
    read_frame_dict = {}
    for index, value in frame_df.items():
        read_length, read_frame = index
        if read_length not in read_frame_dict:
            read_frame_dict[read_length] = {0:0,1:0,2:0}
        read_frame_dict[read_length][read_frame] = value
    return read_frame_dict

def read_frame_cull(read_frame_dict: dict, config: dict) -> dict:
    """
    Culls the read_frame_dict according to config so only read lengths of interest are kept
    
    Inputs:
    read_frame_dict: 
    config: 
    
    Outputs:
    culled_read_frame_dict
    """
    culled_read_frame_dict = read_frame_dict.copy()
    cull_list = list(culled_read_frame_dict.keys())
    for k in cull_list:
        if (
            k > config["plots"]["read_frame_distribution"]["upper_limit"]
            or k < config["plots"]["read_frame_distribution"]["lower_limit"]
        ):
            del culled_read_frame_dict[k]
    
    return culled_read_frame_dict

def read_frame_score(read_frame_dict:dict) -> dict:
    """
    Generates scores for each read_length seperately as well as a global score
    Can be used after read_frame_cull to calculate the global score of the region of interest
    The calculation for this score is: 1 - sum(2nd highest peak count)/sum(highest peak count)
    A score close to 1 has good periodicity, while a score closer to 0 has a random spread
    
    Inputs:
    read_frame_dict: dictionary containing the distribution of the reading frames over the different read lengths
    
    Outputs:
    scored_read_frame_dict: dictionary containing read frame distribution scores for each read length and a global score
    """
    scored_read_frame_dict = {}
    highest_peak_sum, second_peak_sum = 0, 0
    for k, inner_dict in read_frame_dict.items():
        top_two_values = sorted(inner_dict.values(), reverse=True)[:2]
        highest_peak_sum += top_two_values[0]
        second_peak_sum += top_two_values[1]
        scored_read_frame_dict[k] = 1-top_two_values[1]/top_two_values[0]
    scored_read_frame_dict["global"] = (1-second_peak_sum/highest_peak_sum)
    return scored_read_frame_dict


def convert_html_to_pdf(source_html, output_filename):
    result_file = open(output_filename, "w+b")

    pisa_status = pisa.CreatePDF(source_html, dest=result_file)
    result_file.close()
    return pisa_status.err
