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
    calculate_expected_dinucleotide_freqs,
    normalise_ligation_bias,
    nucleotide_composition,
    read_frame_distribution,
    mRNA_distribution,
    annotate_reads,
    sequence_slice,
    metagene_profile,
)

from .metrics import (
    read_length_distribution_metric as rld_metric,
    ligation_bias_distribution_metric as lbd_metric,
    read_frame_distribution_metric as rfd_metric,
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
    results_dict = {}
    results_dict["mode"] = "annotation_free_mode"

    print("> read_length_distribution")
    results_dict["read_length_distribution"] = read_length_distribution(
        read_df
        )
    results_dict["read_length_distribution_metric"] = rld_metric(
        results_dict["read_length_distribution"]
        )

    print("> ligation_bias_distribution")
    results_dict["ligation_bias_distribution"] = ligation_bias_distribution(
        read_df
        )
    results_dict["ligation_bias_distribution_metric"] = lbd_metric(
        results_dict["ligation_bias_distribution"],
        calculate_expected_dinucleotide_freqs(
            read_df,
        ),
        )

    print("> nucleotide_composition")
    results_dict["nucleotide_composition"] = nucleotide_composition(
        read_df
        )

    print("> read_frame_distribution")
    results_dict["read_frame_distribution"] = read_frame_distribution(
        read_df
        )
    results_dict["read_frame_distribution_metric"] = rfd_metric(
        results_dict["read_frame_distribution"]
    )

    print("> sequence_slice")
    results_dict["sequence_slice"] = sequence_slice(
        read_df,
        nt_start=config["plots"]["nucleotide_proportion"]["nucleotide_start"],
        nt_count=config["plots"]["nucleotide_proportion"]["nucleotide_count"],
    )

    print("> summary_metrics")
    results_dict["summary_metrics"] = {}

    return results_dict


def annotation_mode(
    read_df: pd.DataFrame, annotation_df: pd.DataFrame, config: dict
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
    print("Merging annotation and reads")
    annotated_read_df = annotate_reads(read_df, annotation_df)
    print("Subsetting to CDS reads")
    cds_read_df = read_df_to_cds_read_df(annotated_read_df)
    print("Running modules")

    results_dict = {}
    results_dict["mode"] = "annotation_mode"
    print("> read_length_distribution")
    results_dict["read_length_distribution"] = read_length_distribution(
        read_df
        )
    results_dict["read_length_distribution_metric"] = rld_metric(
        results_dict["read_length_distribution"]
        )

    print("> ligation_bias_distribution")
    results_dict["ligation_bias_distribution"] = ligation_bias_distribution(
        read_df,
        num_bases=config["plots"]["ligation_bias_distribution"]["nucleotide_count"],
        five_prime=config["plots"]["ligation_bias_distribution"]["five_prime"],
        )
    results_dict["ligation_bias_distribution_metric"] = lbd_metric(
        results_dict["ligation_bias_distribution"],
        calculate_expected_dinucleotide_freqs(read_df,),
        )
    if config["plots"]["ligation_bias_distribution"]["background_freq"]:
        results_dict["ligation_bias_distribution"] = normalise_ligation_bias(
            read_df,
            results_dict["ligation_bias_distribution"],
            num_bases=config["plots"]["ligation_bias_distribution"]["nucleotide_count"],
            five_prime=config["plots"]["ligation_bias_distribution"]["five_prime"],
        )

    print("> nucleotide_composition")
    results_dict["nucleotide_composition"] = nucleotide_composition(
        read_df
        )

    print("> sequence_slice")
    results_dict["sequence_slice"] = sequence_slice(
        read_df,
        nt_start=config["plots"]["nucleotide_proportion"]["nucleotide_start"],
        nt_count=config["plots"]["nucleotide_proportion"]["nucleotide_count"],
    )

    print("> read_frame_distribution")
    results_dict["read_frame_distribution"] = (
        read_frame_distribution(cds_read_df)
        if config["qc"]["use_cds_subset"]["read_frame_distribution"]
        else read_frame_distribution(read_df)
    )
    results_dict["read_frame_distribution_metric"] = rfd_metric(
        results_dict["read_frame_distribution"]
    )

    print("> mRNA_distribution")
    results_dict["mRNA_distribution"] = mRNA_distribution(
        annotated_read_df
        )

    print("> metagene_profile")
    results_dict["metagene_profile"] = metagene_profile(
        annotated_read_df,
        config["plots"]["metagene_profile"]["distance_target"],
        config["plots"]["metagene_profile"]["distance_range"],
    )
    print("> summary_metrics")
    results_dict["summary_metrics"] = {}

    return results_dict


def sequence_mode(
    read_df: pd.DataFrame,
    gff_path: str,
    transcript_list: list,
    fasta_path: str,
    config: dict,
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
        "read_frame_distribution": read_frame_distribution(read_df),
    }
    # results_dict["read_frame_distribution"] = read_frame_distribution(
    #   cds_read_df)\
    #     if config["qc"]["use_cds_subset"]["read_frame_distribution"]\
    #     else read_frame_distribution(read_df)

    return results_dict
