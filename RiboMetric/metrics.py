'''
This script contains the functions used to calculate individual metrics
for different aspects of ribosome profiling data. The functions are called
are called from qc and the input data comes from the output of their
respective modules

'''

import pandas as pd
import math


def find_category_by_cumulative_percentage(df, percentage):
    """
    Calculate the read_length with cumulative percentages
    """
    df['cumulative_percentage'] = (df['read_count'].cumsum()
                                   / df['read_count'].sum())
    read_length = df.loc[df['cumulative_percentage'] >= percentage,
                         'read_length'].iloc[0]
    return read_length


def read_length_distribution_metric(
        rld_dict: dict,
        ) -> pd.DataFrame:
    """
    Calculate the read length distribution metric from the output of
    the read_length_distribution module.

    This metric is the IQR of the read length distribution and is
    calculated as the difference between the 75th and 25th percentile
    The metric is then normalised by dividing by the range of the
    read length distribution between the 10th and 90th percentile

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module

    Outputs:
        rld_df: Dataframe containing the read length distribution metric
    """
    rld_df = pd.DataFrame.from_dict(rld_dict, orient="index")
    rld_df = rld_df.reset_index()
    rld_df.columns = ["read_length", "read_count"]

    Q3 = find_category_by_cumulative_percentage(rld_df, 0.75)
    Q1 = find_category_by_cumulative_percentage(rld_df, 0.25)
    inter_quartile_range = Q3 - Q1

    max_range = find_category_by_cumulative_percentage(rld_df, 0.9)\
        - find_category_by_cumulative_percentage(rld_df, 0.1)

    return 1 - (inter_quartile_range / max_range)


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
    # Needs possible rewrite using normalised ligation bias.
    # Current iteration only accounts for five_prime
    # division by 0 if background is non-existent, Only patterns that occur
    # at least once are used (needs to be changed in ligation bias)
    kl_divergence = 0.0

    for dinucleotide, observed_prob in observed_freq["five_prime"].items():
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


def calculate_3nt_periodicity_score(probabilities):
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


def read_frame_information_content(
    read_frame_distribution: dict,
        ) -> dict:
    """
    Calculate the read frame distribution metric from the output of
    the read_frame_distribution module.

    This metric is the Shannon entropy of the read frame distribution

    Inputs:
        read_frame_distribution: Dictionary containing the output of the
                read_frame_distribution module

    Outputs:
        frame_info_content_dict: Shannon entropy of the read frame
                distribution where keys are read length and values are tuples
                containing information content in bits and number of reads in
                frame
    """
    pseudocount = 1e-100
    frame_info_content_dict = {}
    for read_length in read_frame_distribution:
        total_count = sum(read_frame_distribution[read_length].values())

        probabilities = []
        for frame, count in read_frame_distribution[read_length].items():
            prob = (count + pseudocount) / (total_count + pseudocount)
            probabilities.append(prob)

        score = calculate_3nt_periodicity_score(probabilities)

        frame_info_content_dict[read_length] = score, total_count

    return frame_info_content_dict


def information_metric_cutoff(
    frame_info_content_dict: dict,
    min_count_threshold: float = 0.05,
        ) -> dict:
    """
    Apply the cut off to the information content metric

    Inputs:
        frame_info_content_dict: Dictionary containing the output of the
                information_metric_cutoff module
        min_count_threshold: Minimum count threshold for a read length to be
                included in the metric

    Outputs:
        information_content_metric: Dictionary containing the information
                content metric for each read length
    """
    information_content_metric = {}
    total_reads = sum(
        frame_info_content_dict[key][1]
        for key in frame_info_content_dict
        )
    for read_length in frame_info_content_dict:
        score, count = frame_info_content_dict[read_length]
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
    frame_info_content_dict: dict,
        ):
    '''
    Produce a single metric for the triplet periodicity by taking the weighted
    average of the scores for each read length.

    Inputs:
        frame_info_content_dict (dict): Dictionary containing the information

    Returns:
        result (float): The triplet periodicity score.
    '''
    total_reads = sum(
        frame_info_content_dict[key][1]
        for key in frame_info_content_dict
        )
    weighted_scores = []
    for _, score in frame_info_content_dict.items():
        weighted_score = score[0] * score[1]
        weighted_scores.append(weighted_score)

    return sum(weighted_scores) / total_reads


def triplet_periodicity_weighted_score_best_3_read_lengths(
        frame_info_content_dict: dict,) -> float:
    """
    Produce a single metric for the triplet periodicity by taking the weighted
    average of the scores for the best 3 read lengths.

    Inputs:
        frame_info_content_dict: Dictionary containing the information
                content metric and total counts for each read length

    Returns:
        result: The triplet periodicity score
    """
    total_reads = sum(
        frame_info_content_dict[key][1]
        for key in frame_info_content_dict
        )
    sorted_counts = sorted(
        frame_info_content_dict.items(),
        key=lambda x: x[1][1],
        reverse=True,
        )[:3]
    weighted_scores = []
    for _, score in sorted_counts:
        weighted_score = score[0] * score[1]
        weighted_scores.append(weighted_score)

    return sum(weighted_scores) / total_reads


def frame_bias_metric(
        read_df: pd.DataFrame,
        ) -> dict:
    """
    Calculate the frame bias metric. This metric is the proportion 
    of reads in each frame

    Inputs:
        read_df: Dataframe containing the reads

    Outputs:
        frame_bias_dict: Dictionary containing the frame bias metric
    """
    frame_bias_dict = {}
    total_reads = len(read_df)

    for frame in range(3):
        frame_reads = len(read_df[read_df["a_site"] % 3 == frame])
        frame_bias_dict[frame] = frame_reads / total_reads

    return frame_bias_dict
