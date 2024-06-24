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
    terminal_nucleotide_bias_distribution,
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
    read_length_distribution_IQR_normalised_metric as rld_metric,
    read_length_distribution_coefficient_of_variation_metric as rldv_metric,
    read_length_distribution_normality_metric as rldn_metric,
    terminal_nucleotide_bias_KL_metric as lbd_metric,
    terminal_nucleotide_bias_max_absolute_metric as lbmp_metric,
    read_frame_information_content as rf_info_metric,
    read_frame_information_weighted_score,
    information_metric_cutoff,
    cds_coverage_metric,
    region_region_ratio_metric,
    read_length_distribution_max_prop_metric as rldpp_metric,
    periodicity_autocorrelation,
    uniformity_autocorrelation,
    uniformity_entropy,
    uniformity_theil_index,
    uniformity_gini_index,
    periodicity_dominance,
    fourier_transform,
    read_length_distribution_bimodality,
    proportion_of_reads_in_region
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

    #######################################################################
    # READ LENGTH DISTRIBUTION
    #######################################################################
    print("> read_length_distribution")
    results_dict["read_length_distribution"] = read_length_distribution(
        read_df
    )
    results_dict["metrics"][
        "read_length_distribution_IQR_metric"
        ] = rld_metric(
        results_dict["read_length_distribution"]
    )
    results_dict["metrics"][
        "read_length_distribution_bimodality_metric"
        ] = read_length_distribution_bimodality(
            results_dict["read_length_distribution"]
        )
    results_dict["metrics"][
        "read_length_distribution_normality_metric"
        ] = rldn_metric(
            results_dict["read_length_distribution"]
        )

    results_dict["metrics"][
        "read_length_distribution_coefficient_of_variation_metric"
        ] = rldv_metric(
            results_dict["read_length_distribution"]
        )
    results_dict["metrics"][
        "read_length_distribution_maxprop_metric"] = rldpp_metric(
        results_dict["read_length_distribution"],
        num_top_readlens=1
    )

    #######################################################################
    # TERMINAL NUCLEOTIDE BIAS
    #######################################################################
    if sequence_background:
        print("> terminal_nucleotide_bias_distribution")
        results_dict[
            "terminal_nucleotide_bias_distribution"
            ] = terminal_nucleotide_bias_distribution(
            read_df,
            pattern_length=config[
                "plots"][
                "terminal_nucleotide_bias_distribution"][
                "nucleotide_count"
            ],
        )
        results_dict["metrics"][
            "terminal_nucleotide_bias_distribution_5_prime_metric"
            ] = lbd_metric(
            results_dict["terminal_nucleotide_bias_distribution"],
            sequence_background["5_prime_bg"],
            prime="five_prime",
        )
        results_dict["metrics"][
            "terminal_nucleotide_bias_distribution_3_prime_metric"
            ] = lbd_metric(
            results_dict["terminal_nucleotide_bias_distribution"],
            sequence_background["3_prime_bg"],
            prime="three_prime",
        )
        results_dict["metrics"][
            "terminal_nucleotide_bias_max_absolute_metric_5_prime_metric"
            ] = lbmp_metric(
            results_dict["terminal_nucleotide_bias_distribution"],
            sequence_background["5_prime_bg"],
            prime="five_prime",
        )
        results_dict["metrics"][
            "terminal_nucleotide_bias_max_absolute_metric_3_prime_metric"
            ] = lbmp_metric(
            results_dict["terminal_nucleotide_bias_distribution"],
            sequence_background["3_prime_bg"],
            prime="three_prime",
        )
        if config["plots"][
                "terminal_nucleotide_bias_distribution"][
                "background_freq"]:
            results_dict[
                "terminal_nucleotide_bias_distribution"
                ] = normalise_ligation_bias(
                results_dict["terminal_nucleotide_bias_distribution"],
                sequence_background=sequence_background,
                pattern_length=config[
                    "plots"][
                    "terminal_nucleotide_bias_distribution"][
                    "nucleotide_count"
                ],
            )

        print("> nucleotide_composition")
        results_dict["nucleotide_composition"] = nucleotide_composition(
            sequence_data)

    print("> read_frame_distribution")
    if annotation:
        coding_metagene = metagene_profile(
                annotated_read_df,
                target="start",
                distance_range=[30, 117],
            )

        #######################################################################
        # Periodicity
        #######################################################################
        results_dict["metrics"][
            "periodicity_autocorrelation"
            ] = periodicity_autocorrelation(
            coding_metagene.copy()
        )
        results_dict["metrics"]["periodicity_fourier"] = fourier_transform(
            coding_metagene.copy()
        )

        results_dict["reading_frame_triangle"] = reading_frame_triangle(
                annotated_read_df
            )
        read_frame_dist = (
            read_frame_distribution_annotated(cds_read_df)
            if config["qc"]["use_cds_subset"]["read_frame_distribution"]
            and annotation
            else read_frame_distribution_annotated(annotated_read_df)
            )
        frame_info_content_dict = rf_info_metric(read_frame_dist)
        results_dict["read_frame_distribution"] = read_frame_dist
        results_dict["metrics"]["periodicity_information"] =\
            information_metric_cutoff(
                frame_info_content_dict,
                config['qc']['read_frame_distribution']['3nt_count_cutoff']
            )

        results_dict["metrics"]["periodicity_information_weighted_score"] = \
            read_frame_information_weighted_score(
                frame_info_content_dict,
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

        #######################################################################
        # UNIFORMITY
        #######################################################################
        results_dict["metrics"][
            "uniformity_autocorrelation"
            ] = uniformity_autocorrelation(
            coding_metagene.copy()
        )
        results_dict["metrics"]["uniformity_entropy"] = uniformity_entropy(
            coding_metagene.copy()
        )
        results_dict["metrics"][
            "uniformity_theil_index"
            ] = uniformity_theil_index(
            coding_metagene.copy()
        )
        results_dict["metrics"][
            "uniformity_gini_index"
            ] = uniformity_gini_index(
            coding_metagene.copy()
        )

        #######################################################################
        # COVERAGE
        #######################################################################
        print("> cds_coverage_metric")
        results_dict["metrics"]["CDS_coverage_metric"] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=config["qc"]["cds_coverage"]["in_frame_coverage"]
            )
        results_dict["metrics"][
            "CDS_coverage_metric_not_inframe_1read_1000tx"
            ] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=False,
            num_transcripts=1000
            )
        results_dict["metrics"][
            "CDS_coverage_metric_not_inframe_100read_100tx"
            ] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=100,
            in_frame_coverage=False,
            num_transcripts=100
            )
        results_dict["metrics"][
            "CDS_coverage_metric_inframe_1read_1000tx"
            ] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=True,
            num_transcripts=1000
            )
        results_dict["metrics"][
            "CDS_coverage_metric_inframe_100read_100tx"
            ] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=100,
            in_frame_coverage=True,
            num_transcripts=100
            )

        #######################################################################
        # RNA REGIONAL SUPPORT
        #######################################################################
        results_dict["metrics"]["ratio_cds:leader"] =\
            region_region_ratio_metric(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region1="five_leader",
                region2="CDS",
            )
        results_dict["metrics"]["ratio_cds:trailer"] =\
            region_region_ratio_metric(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region1="CDS",
                region2="three_trailer",
            )
        results_dict["metrics"]["ratio_leader:trailer"] =\
            region_region_ratio_metric(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region1="five_leader",
                region2="three_trailer",
            )

        results_dict["metrics"]["prop_reads_CDS"] =\
            proportion_of_reads_in_region(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region="CDS",
            )
        results_dict["metrics"]["prop_reads_leader"] =\
            proportion_of_reads_in_region(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region="five_leader",
            )
        results_dict["metrics"]["prop_reads_trailer"] =\
            proportion_of_reads_in_region(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region="three_trailer",
            )
    else:
        read_frame_dist = read_frame_distribution(read_df)
        results_dict["read_frame_distribution"] = read_frame_dist

    culled_read_frame_dict = read_frame_cull(read_frame_dist, config)
    results_dict["metrics"][
        "periodicity_trips-viz"
        ] = read_frame_score_trips_viz(
        culled_read_frame_dict)

    results_dict["metrics"]["periodicity_dominance"] = periodicity_dominance(
        culled_read_frame_dict
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
        "terminal_nucleotide_bias_distribution": terminal_nucleotide_bias_distribution(read_df),
        "nucleotide_composition": nucleotide_composition(read_df),
        "read_frame_distribution": read_frame_distribution(read_df),
    }
    # results_dict["read_frame_distribution"] = read_frame_distribution(
    #   cds_read_df)\
    #     if config["qc"]["use_cds_subset"]["read_frame_distribution"]\
    #     else read_frame_distribution(read_df)

    return results_dict
