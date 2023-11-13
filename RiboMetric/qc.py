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
    chunked_annotate_reads,
    assign_mRNA_category,
    read_length_distribution,
    read_df_to_cds_read_df,
    ligation_bias_distribution,
    normalise_ligation_bias,
    nucleotide_composition,
    read_frame_distribution,
    read_frame_distribution_annotated,
    mRNA_distribution,
    metagene_profile,
    reading_frame_triangle,
    read_frame_score,
    read_frame_cull,
)

from .metrics import (
    read_length_distribution_metric as rld_metric,
    ligation_bias_distribution_metric as lbd_metric,
    read_frame_information_content as rfd_metric,
    triplet_periodicity_weighted_score,
    triplet_periodicity_best_read_length_score as tpbrl_metric,
    information_metric_cutoff,
    triplet_periodicity_weighted_score_best_3_read_lengths as tpw3rl_metric,
    cds_coverage_metric
)


def annotation_mode(
    read_df: pd.DataFrame,
    sequence_data: dict,
    sequence_background: dict,
    annotation_df: pd.DataFrame = pd.DataFrame(),
    config: dict = {}
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
    if len(annotation_df) > 0:
        annotation = True
        print("Merging annotation and reads")
        annotated_read_df = chunked_annotate_reads(read_df, annotation_df)

        print("assigning mRNA categories")
        annotated_read_df = assign_mRNA_category(annotated_read_df)
        print("Subsetting to CDS reads")
        cds_read_df = read_df_to_cds_read_df(annotated_read_df)
    else:
        annotation = False
    print("Running modules")

    results_dict = {}
    results_dict["mode"] = ("annotation"
                            if annotation
                            else "annotation_free")
    results_dict["metrics"] = {}
    print("> read_length_distribution")
    results_dict["read_length_distribution"] = read_length_distribution(
        read_df
    )
    results_dict["metrics"]["read_length_distribution_metric"] = rld_metric(
        results_dict["read_length_distribution"]
    )

    if sequence_background:
        print("> ligation_bias_distribution")
        results_dict["ligation_bias_distribution"] = ligation_bias_distribution(
            read_df,
            pattern_length=config["plots"]["ligation_bias_distribution"][
                "nucleotide_count"
            ]
        )
        results_dict["metrics"]["ligation_bias_distribution_metric"] = lbd_metric(
            results_dict["ligation_bias_distribution"],
            sequence_background["5_prime_bg"],
        )
        if config["plots"]["ligation_bias_distribution"]["background_freq"]:
            results_dict["ligation_bias_distribution"] = normalise_ligation_bias(
                results_dict["ligation_bias_distribution"],
                sequence_background=sequence_background,
                pattern_length=config["plots"]["ligation_bias_distribution"][
                    "nucleotide_count"
                ],
            )

        print("> nucleotide_composition")
        results_dict["nucleotide_composition"] = nucleotide_composition(
            sequence_data)

    print("> read_frame_distribution")
    if annotation:
        read_frame_dist = (
            read_frame_distribution_annotated(cds_read_df)
            if config["qc"]["use_cds_subset"]["read_frame_distribution"]
            and annotation
            else read_frame_distribution_annotated(annotated_read_df)
            )
        print("> reading_frame_proportions")
        results_dict["reading_frame_triangle"] = reading_frame_triangle(
                annotated_read_df
            )

    else:
        read_frame_dist = (read_frame_distribution(read_df))

    frame_info_content_dict = rfd_metric(read_frame_dist)
    results_dict["read_frame_distribution"] = read_frame_dist
    results_dict["metrics"]["read_frame_information_metric"] =\
        information_metric_cutoff(
            frame_info_content_dict,
            config['qc']['read_frame_distribution']['3nt_count_cutoff']
        )

    culled_read_frame_dict = read_frame_cull(read_frame_dist, config)
    results_dict["metrics"]["read_frame_bias"] = read_frame_score(
        culled_read_frame_dict)

    results_dict["metrics"]["3nt_weighted_score"] = \
        triplet_periodicity_weighted_score(
            frame_info_content_dict,
        )
    results_dict["metrics"]["3nt_weighted_score_best_3_read_lengths"] = \
        tpw3rl_metric(
            frame_info_content_dict,
    )
    results_dict["metrics"]["3nt_best_read_length_score"] = tpbrl_metric(
        frame_info_content_dict,
    )

    if annotation:
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
        print("> cds_coverage_metric")
        results_dict["metrics"]["CDS_coverage_metric"] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=config["qc"]["cds_coverage"]["in_frame_coverage"]
            )
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
