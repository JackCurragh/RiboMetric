"""
Main script for running qc analysis

Three main modes:
    annoation free: no gff file provided just use the bam file
    annotation based: gff file provided and use the bam file
    sequence based: gff file and transcriptome fasta file
                    provided and use the bam file

"""

import pandas as pd
from .modules import (
    read_length_distribution,
    a_site_calculation,
    read_df_to_cds_read_df,
)
from .file_parser import gff_df_to_cds_df


def annotation_free_mode(read_df: pd.DataFrame, config: str) -> dict:
    """
    Run the annotation free mode of the qc analysis

    Inputs:
        read_df: dataframe containing the read information
                (keys are the read names)
        config: Path to the config file

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    """
    results_dict = {}

    results_dict['read_len_dist'] = read_length_distribution(read_df)
    return results_dict


def annotation_mode(
    read_df: pd.DataFrame,
    gff_df: pd.DataFrame,
    transcript_list: list,
    config: str
) -> dict:
    """
    Run the annotation mode of the qc analysis

    Inputs:
        read_df: Dataframe containing the read information
                (keys are the read names)
        gffdf: Dataframe containing the gff information
        transcript_list: List of the top N transcripts

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    """
    print("Extracting CDS information from gff file")
    cds_df = gff_df_to_cds_df(gff_df, transcript_list)
    print("Calculating A site information")
    a_site_df = a_site_calculation(read_df)
    print("Subsetting to CDS reads")
    cds_read_df = read_df_to_cds_read_df(a_site_df, cds_df)
    read_length_distribution(cds_read_df)
    results_dict = {}
    results_dict['read_len_dist'] = read_length_distribution(a_site_df)

    return results_dict


def sequence_mode(
    read_df: pd.DataFrame,
    gff_path: str,
    transcript_list: list,
    fasta_path: str
) -> dict:
    """
    Run the sequence mode of the qc analysis

    Inputs:
        read_df: dataframe containing the read information
                (keys are the read names)
        gff_path: Path to the gff file
        transcript_list: List of the top N transcripts
        fasta_path: Path to the transcriptome fasta file

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    """
    results_dict = {}

    return results_dict
