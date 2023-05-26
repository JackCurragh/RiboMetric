'''
This script contains the functions used to calculate individual metrics
for different aspects of ribosome profiling data. The functions are called
are called from qc and the input data comes from the output of their
respective modules

'''

import pandas as pd
import math


def read_length_distribution_metric(
        rld_dict: dict,
        ) -> pd.DataFrame:
    """
    Calculate the read length distribution metric from the output of
    the read_length_distribution module.

    This metric is the IQR of the read length distribution and is
    calculated as the difference between the 75th and 25th percentile

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module

    Outputs:
        rld_df: Dataframe containing the read length distribution metric
    """
    rld_df = pd.DataFrame.from_dict(rld_dict, orient="index")
    rld_df = rld_df.reset_index()
    rld_df.columns = ["read_length", "read_count"]

    Q3 = rld_df["read_length"].quantile(0.75)
    Q1 = rld_df["read_length"].quantile(0.25)

    return Q3 - Q1


def ligation_bias_distribution_metric(
        observed_freq: dict,
        expected_freq: dict,
        ) -> float:
    """
    Calculate the ligation bias metric from the output of
    the ligation_bias_distribution module.

    This metric is the K-L divergence of the ligation bias distribution
    of the observed frequencies from the expected frequencies. The
    expected frequencies are calculated from the nucleotide composition
    of the genome.

    Inputs:
        observed_freq: Dictionary containing the output of the
                ligation_bias_distribution module
        expected_freq: Dictionary containing the expected frequencies

    Outputs:
        lbd_df: Dataframe containing the ligation bias metric in bits
    """
    kl_divergence = 0.0

    for dinucleotide, observed_prob in observed_freq.items():
        expected_prob = expected_freq[dinucleotide]
        kl_divergence += observed_prob * math.log2(
                                            observed_prob / expected_prob
                                            )

    return kl_divergence


def cds_coverage_metric(
        cds_read_df: pd.DataFrame,
        minimum_reads: int = 1,
        in_frame_coverage: bool = True
        ) -> float:
    """
    Calculates the proportion of CDS covered by ribosomal protected fragments

    Inputs:
        annotated_read_df: Dataframe containing the reads that have
        a transcript available in the provided annotation
        minimum_reads: The minimum amount of reads that should cover
        a specific nucleotide to be counted for the proportion
        in_frame_count: If set to True, only controls the coverage in frame

    Outputs:
        cds_coverage: A proportion of the amount of individual nucleotides
        represented by the A-sites over the total number of nucleotides in
        the CDS of transcripts present in the reads
    """
    # Create the cds_coverage_df that contains only the required columns and
    # the "name_pos" column, combining the transcript_id and a_site
    cds_coverage_df = cds_read_df[["transcript_id",
                                   "a_site",
                                   "cds_start",
                                   "cds_end"]].copy()
    cds_coverage_df["name_pos"] = (cds_coverage_df["transcript_id"]
                                   .astype("object")
                                   + cds_coverage_df["a_site"]
                                   .astype(str)
                                   ).astype("category")

    # Calculate the total combined length of the CDS of transcripts that have
    # reads aligned to them
    cds_transcripts = cds_coverage_df[~cds_coverage_df["transcript_id"]
                                      .duplicated()].copy()
    cds_transcripts["cds_length"] = (cds_transcripts
                                     .apply(lambda x: x['cds_end']
                                            - x['cds_start'],
                                            axis=1))
    cds_length_total = cds_transcripts["cds_length"].sum()
    del cds_transcripts

    # If in_frame_coverage is true, take only reads that are in frame for
    # their transcript and divide the combined CDS length by 3
    if in_frame_coverage:
        cds_coverage_df = cds_coverage_df[
            (cds_coverage_df["a_site"]
             - cds_coverage_df["cds_start"]
             ) % 3 == 0]
        cds_length_total = cds_length_total/3

    # Calculate the count of nucleotides covered by the reads after filtering
    cds_reads_count = sum(cds_coverage_df.value_counts("name_pos")
                          > minimum_reads)

    return cds_reads_count/cds_length_total


def calculate_score(probabilities):
    '''
    Calculate the triplet periodicity score for a given probability of a read
    being in frame. The score is the square root of the bits of information in
    the triplet distribution.

    Numerator is the Maximum Entropy of the triplet distribution minus the
    entropy of the triplet distribution.
    Denominator is the Maximum Entropy of the triplet distribution.

    Inputs:
        probability (float): The probability of a read being in frame.

    Returns:
        result (float): The triplet periodicity score.
    '''
    maximum_entropy = math.log2(3)
    entropy = 0
    for probability in probabilities:
        entropy += -(probability * math.log2(probability))

    result = math.sqrt((maximum_entropy - entropy) / maximum_entropy)
    return result


def read_frame_distribution_information_content_metric(
    read_frame_distribution: dict,
        ) -> float:
    """
    Calculate the read frame distribution metric from the output of
    the read_frame_distribution module.

    This metric is the Shannon entropy of the read frame distribution

    Inputs:
        read_frame_distribution: Dictionary containing the output of the
                read_frame_distribution module

    Outputs:
        read_frame_distribution_metric: Shannon entropy of the read frame
                distribution
    """
    pseudocount = 1e-100
    pre_scores = {}
    for read_length in read_frame_distribution:
        total_count = sum(read_frame_distribution[read_length].values())

        probabilities = []
        for frame, count in read_frame_distribution[read_length].items():
            prob = (count + pseudocount) / total_count
            probabilities.append(prob)

        score = calculate_score(probabilities)

        pre_scores[read_length] = score, total_count

    return pre_scores


def information_metric_cutoff(
    pre_scores: dict,
    min_count_threshold: float = 0.05,
        ) -> dict:
    """
    Apply the cut off to the information content metric

    Inputs:
        pre_scores: Dictionary containing the output of the
                read_frame_distribution_information_content_metric module
        min_count_threshold: Minimum count threshold for a read length to be
                included in the metric

    Outputs:
        information_content_metric: Dictionary containing the information
                content metric for each read length
    """
    information_content_metric = {}
    total_reads = sum(
        pre_scores[key][1]
        for key in pre_scores
        )
    for read_length in pre_scores:
        score, count = pre_scores[read_length]
        if count > total_reads * min_count_threshold:
            information_content_metric[read_length] = score
    return information_content_metric


def triplet_periodicity_best_read_length_score(information_content_metric):
    '''
    Produce a single metric for the triplet periodicity by taking the maximum
    score across all read lengths.

    Inputs:
        information_content_metric (dict): The information content metric
            for each read length.

    Returns:
        result (float): The triplet periodicity score.
    '''
    return max(information_content_metric.values())


def triplet_periodicity_weighted_score(
    pre_scores: dict,
        ):
    '''
    Produce a single metric for the triplet periodicity by taking the weighted
    average of the scores for each read length.

    Inputs:
        pre_scores (dict): Dictionary containing the information

    Returns:
        result (float): The triplet periodicity score.
    '''
    total_reads = sum(
        pre_scores[key][1]
        for key in pre_scores
        )
    weighted_scores = []
    for _, score in pre_scores.items():
        weighted_score = score[0] * score[1]
        weighted_scores.append(weighted_score)

    return sum(weighted_scores) / total_reads


def triplet_periodicity_weighted_score_best_3_read_lengths(
        pre_scores: dict,) -> float:
    """
    Produce a single metric for the triplet periodicity by taking the weighted
    average of the scores for the best 3 read lengths.

    Inputs:
        pre_scores: Dictionary containing the information
                content metric and total counts for each read length

    Returns:
        result: The triplet periodicity score
    """
    total_reads = sum(
        pre_scores[key][1]
        for key in pre_scores
        )
    sorted_counts = sorted(
        pre_scores.items(),
        key=lambda x: x[1][1],
        reverse=True,
        )[:3]
    weighted_scores = []
    for _, score in sorted_counts:
        weighted_score = score[0] * score[1]
        weighted_scores.append(weighted_score)

    return sum(weighted_scores) / total_reads
