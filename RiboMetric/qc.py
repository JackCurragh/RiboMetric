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
    read_frame_score_trips_viz,
    read_frame_cull,
    asite_calculation_per_readlength,
    a_site_calculation_variable_offset
)

from .metrics import (
    read_length_distribution_spread_metric as rld_metric,
    read_length_distribution_variation_metric as rldv_metric,
    ligation_bias_distribution_metric as lbd_metric,
    ligation_bias_max_proportion_metric as lbmp_metric,
    read_frame_information_content as rf_info_metric,
    read_frame_information_weighted_score,
    read_frame_information_best_read_length_score as tpbrl_metric,
    information_metric_cutoff,
    read_frame_information_weighted_score_best_3_read_lengths as tpw3rl_metric,
    cds_coverage_metric,
    leader_cds_ratio_metric,
    read_length_distribution_prop_at_peak_metric as rldpp_metric,
    autocorrelation,
    uniformity,
    theil_index,
    gini_index,
    read_frame_dominance,
    fourier_transform,
    multitaper,
    wavelet_transform
)
from typing import Any, Dict


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
        offsets = asite_calculation_per_readlength(annotated_read_df)
        annotated_read_df = a_site_calculation_variable_offset(
            annotated_read_df, offsets
            )
        print("Subsetting to CDS reads")
        cds_read_df = read_df_to_cds_read_df(annotated_read_df)
        if cds_read_df.empty:
            print("No CDS reads found")
            annotation = False

    else:
        annotation = False
    print("Running modules")

    results_dict: Dict[str, Any] = {
        "mode": ("annotation" if annotation else "annotation_free"),
        "metrics": {}
    }

    print("> read_length_distribution")
    results_dict["read_length_distribution"] = read_length_distribution(
        read_df
    )
    results_dict["metrics"]["read_length_distribution_metric"] = rld_metric(
        results_dict["read_length_distribution"]
    )
    results_dict["metrics"]["read_length_distribution_variation_metric"] =\
        rldv_metric(
            results_dict["read_length_distribution"]
        )
    results_dict["metrics"]["rld_prop_top_1_metric"] = rldpp_metric(
        results_dict["read_length_distribution"],
        num_top_readlens=1
    )
    results_dict["metrics"]["rld_prop_top_3_metric"] = rldpp_metric(
        results_dict["read_length_distribution"],
        num_top_readlens=3
    )
    results_dict["metrics"]["rld_prop_top_5_metric"] = rldpp_metric(
        results_dict["read_length_distribution"],
        num_top_readlens=5
    )

    if sequence_background:
        print("> ligation_bias_distribution")
        results_dict[
            "ligation_bias_distribution"
            ] = ligation_bias_distribution(
            read_df,
            pattern_length=config["plots"]["ligation_bias_distribution"][
                "nucleotide_count"
            ],
        )
        results_dict["metrics"][
            "ligation_bias_distribution_5_prime_metric"
            ] = lbd_metric(
            results_dict["ligation_bias_distribution"],
            sequence_background["5_prime_bg"],
            prime="five_prime",
        )
        results_dict["metrics"][
            "ligation_bias_distribution_3_prime_metric"
            ] = lbd_metric(
            results_dict["ligation_bias_distribution"],
            sequence_background["3_prime_bg"],
            prime="three_prime",
        )
        results_dict["metrics"][
            "ligation_bias_max_proportion_metric_5_prime"
            ] = lbmp_metric(
            results_dict["ligation_bias_distribution"],
            sequence_background["5_prime_bg"],
            prime="five_prime",
        )
        results_dict["metrics"][
            "ligation_bias_max_proportion_metric_3_prime"
            ] = lbmp_metric(
            results_dict["ligation_bias_distribution"],
            sequence_background["3_prime_bg"],
            prime="three_prime",
        )
        if config["plots"]["ligation_bias_distribution"]["background_freq"]:
            results_dict[
                "ligation_bias_distribution"
                ] = normalise_ligation_bias(
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
        read_frame_dist_exclude_15 = (
            read_frame_distribution_annotated(cds_read_df, exclusion_length=15)
            if config["qc"]["use_cds_subset"]["read_frame_distribution"]
            and annotation
            else read_frame_distribution_annotated(annotated_read_df)
            )
        read_frame_dist_exclude_60 = (
            read_frame_distribution_annotated(cds_read_df, exclusion_length=60)
            if config["qc"]["use_cds_subset"]["read_frame_distribution"]
            and annotation
            else read_frame_distribution_annotated(annotated_read_df)
            )
        print("> reading_frame_proportions")
        results_dict["reading_frame_triangle"] = reading_frame_triangle(
                annotated_read_df
            )
        read_frame_dist_28_to_32 = (
            read_frame_distribution_annotated(
                cds_read_df,
                read_length_range=(28, 32)
                )
            if config["qc"]["use_cds_subset"]["read_frame_distribution"]
            and annotation
            else read_frame_distribution_annotated(annotated_read_df)
            )
        frame_info_content_dict = rf_info_metric(read_frame_dist)
        results_dict["read_frame_distribution"] = read_frame_dist
        results_dict["metrics"]["read_frame_information_metric"] =\
            information_metric_cutoff(
                frame_info_content_dict,
                config['qc']['read_frame_distribution']['3nt_count_cutoff']
            )
        frame_info_content_dict_exclude_15 = rf_info_metric(read_frame_dist_exclude_15)
        results_dict["read_frame_distribution"] = read_frame_dist_exclude_15
        results_dict["metrics"]["read_frame_information_metric_exclude_15"] =\
            information_metric_cutoff(
                frame_info_content_dict_exclude_15,
                config['qc']['read_frame_distribution']['3nt_count_cutoff']
            )
        frame_info_content_dict_exclude_60 = rf_info_metric(read_frame_dist_exclude_60)
        results_dict["read_frame_distribution"] = read_frame_dist_exclude_60
        results_dict["metrics"]["read_frame_information_metric_exclude_60"] =\
            information_metric_cutoff(
                frame_info_content_dict_exclude_60,
                config['qc']['read_frame_distribution']['3nt_count_cutoff']
            )

        frame_info_content_dict_28_to_32 = rf_info_metric(read_frame_dist_28_to_32)
        results_dict["read_frame_distribution_28_to_32"] = read_frame_dist_28_to_32
        results_dict["metrics"]["read_frame_information_metric_28_to_32"] =\
            information_metric_cutoff(
                frame_info_content_dict_28_to_32,
                config['qc']['read_frame_distribution']['3nt_count_cutoff']
            )

        results_dict["metrics"]["3nt_weighted_score"] = \
            read_frame_information_weighted_score(
                frame_info_content_dict,
            )
        results_dict["metrics"]["3nt_weighted_score_best_3_read_lengths"] = \
            tpw3rl_metric(
                frame_info_content_dict,
        )
        results_dict["metrics"]["3nt_best_read_length_score"] = tpbrl_metric(
            frame_info_content_dict,
        )
        read_frame_dist = read_frame_distribution(read_df)
        results_dict[
            "read_frame_distribution_best_frame_per_tx"] = read_frame_dist
    else:
        read_frame_dist = read_frame_distribution(read_df)
        results_dict["read_frame_distribution"] = read_frame_dist

    culled_read_frame_dict = read_frame_cull(read_frame_dist, config)
    results_dict["metrics"]["read_frame_score_trips-viz"] = read_frame_score_trips_viz(
        culled_read_frame_dict)["global"]
    results_dict["metrics"]["read_frame_dominance"] = read_frame_dominance(
        culled_read_frame_dict
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

        coding_metagene = metagene_profile(
                annotated_read_df,
                target="start",
                distance_range=[15, 100],
            )
        results_dict["metrics"]["autocorrelation"] = autocorrelation(
            coding_metagene
        )
        results_dict["metrics"]["uniformity"] = uniformity(
            coding_metagene
        )
        results_dict["metrics"]["thiel_index"] = theil_index(
            coding_metagene
        )
        results_dict["metrics"]["gini_index"] = gini_index(
            coding_metagene
        )
        results_dict["metrics"]["fourier"] = fourier_transform(
            coding_metagene
        )
        print()
        print()
        print(type(fourier_transform(coding_metagene)))
        results_dict["metrics"]["multitaper"] = multitaper(
            coding_metagene
        )
        results_dict["metrics"]["wavelet"] = wavelet_transform(
            coding_metagene
        )
        print("> cds_coverage_metric")
        results_dict["metrics"]["CDS_coverage_metric"] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=config["qc"]["cds_coverage"]["in_frame_coverage"]
            )
        results_dict["metrics"][
            "CDS_coverage_metric_not_inframe"
            ] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=False
            )
        results_dict["metrics"][
            "CDS_coverage_metric_top_10"
            ] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=True,
            num_transcripts=10
            )
        results_dict["metrics"][
            "CDS_coverage_metric_top_100"
            ] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=True,
            num_transcripts=100
            )

        results_dict["metrics"][
            "CDS_coverage_metric_top_500"
            ] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=True,
            num_transcripts=500
            )

        results_dict["metrics"]["ratio_cds:leader"] = leader_cds_ratio_metric(
            mRNA_distribution=results_dict["mRNA_distribution"]
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
