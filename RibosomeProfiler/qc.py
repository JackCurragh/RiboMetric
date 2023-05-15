"""
Main script for running qc analysis

Three main modes:
    annotation free: no gff file provided just use the bam file
    annotation based: gff file provided and use the bam file
    sequence based: gff file and transcriptome fasta file
                    provided and use the bam file

"""

import pandas as pd
from .modules import (
    read_length_distribution,
    read_df_to_cds_read_df,
    ligation_bias_distribution,
    nucleotide_composition,
    read_frame_distribution,
    mRNA_distribution,
    annotate_reads,
)


def annotation_free_mode(read_df: pd.DataFrame, config: dict) -> dict:
    """
    Run the annotation free mode of the qc analysis

    Inputs:
        read_df: dataframe containing the read information
                (keys are the read names)
        config:  Dictionary containing the configuration information

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    """

    print("Running modules")
    results_dict = {
        "mode": "annotation_free_mode",
        "read_length_distribution": read_length_distribution(read_df),
        "ligation_bias_distribution": ligation_bias_distribution(read_df),
        "nucleotide_composition": nucleotide_composition(read_df),
        "read_frame_distribution": read_frame_distribution(read_df),
    }

    return results_dict


def annotation_mode(
    read_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    config: dict
) -> dict:
    """
    Run the annotation mode of the qc analysis

    Inputs:
        read_df: Dataframe containing the read information
                (keys are the read names)
        annotation_df: Dataframe containing the annotation information
        transcript_list: List of the top N transcripts
        config: Dictionary containing the configuration information

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    """
    print("Subsetting to CDS reads")
    annotation_df["transcript_id"] = annotation_df["transcript_id"].replace(
        "\"", "", regex=True
        )

    cds_read_df = read_df_to_cds_read_df(read_df, annotation_df)
    print("Merging annotation and reads")
    annotated_read_df = annotate_reads(read_df, annotation_df)
    print("Running modules")
    results_dict = {
        "mode": "annotation_mode",
        "read_length_distribution": read_length_distribution(read_df),
        "ligation_bias_distribution": ligation_bias_distribution(read_df),
        "nucleotide_composition": nucleotide_composition(read_df),
    }
    results_dict["read_frame_distribution"] = read_frame_distribution(
                                                    cds_read_df)\
        if config["qc"]["use_cds_subset"]["read_frame_distribution"]\
        else read_frame_distribution(read_df)
    results_dict["mRNA_distribution"] = mRNA_distribution(annotated_read_df)

    return results_dict


def sequence_mode(
    read_df: pd.DataFrame,
    gff_path: str,
    transcript_list: list,
    fasta_path: str,
    config: dict
) -> dict:
    """
    Run the sequence mode of the qc analysis

    Inputs:
        read_df: dataframe containing the read information
                (keys are the read names)
        gff_path: Path to the gff file
        transcript_list: List of the top N transcripts
        fasta_path: Path to the transcriptome fasta file
        config: Dictionary containing the configuration information

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    """
    results_dict = {
        "mode": "sequence_mode",
        "read_length_distribution": read_length_distribution(read_df),
        "ligation_bias_distribution": ligation_bias_distribution(read_df),
        "nucleotide_composition": nucleotide_composition(read_df),
        "read_frame_distribution": read_frame_distribution(read_df)
    }
    # results_dict["read_frame_distribution"] = read_frame_distribution(
    #   cds_read_df)\
    #     if config["qc"]["use_cds_subset"]["read_frame_distribution"]\
    #     else read_frame_distribution(read_df)

    return results_dict
