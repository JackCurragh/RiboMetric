'''
This script contains the functions used to calculate individual metrics
for different aspects of ribosome profiling data. The functions are called
are called from qc and the input data comes from the output of their
respective modules

'''

import pandas as pd
import math
import numpy as np

from typing import Dict
import numpy as np
from scipy.stats import skew, kurtosis, normaltest
import scipy.signal as signal


def find_category_by_cumulative_percentage(df, percentage):
    """
    Calculate the read_length with cumulative percentages
    """
    df['cumulative_percentage'] = (df['read_count'].cumsum()
                                   / df['read_count'].sum())
    read_length = df.loc[df['cumulative_percentage'] >= percentage,
                         'read_length'].iloc[0]
    return read_length


def read_length_distribution_IQR_normalised_metric(
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
    rld_df.columns = pd.Index(["read_length", "read_count"])

    Q3 = find_category_by_cumulative_percentage(rld_df, 0.75)
    Q1 = find_category_by_cumulative_percentage(rld_df, 0.25)
    inter_quartile_range = Q3 - Q1

    max_range = find_category_by_cumulative_percentage(rld_df, 0.9)\
        - find_category_by_cumulative_percentage(rld_df, 0.1)

    return 1 - (inter_quartile_range / max_range)


def read_length_distribution_coefficient_of_variation_metric(
        rld_dict: dict,
        ) -> float:
    """
    Calculate the read length distribution metric from the output of
    the read_length_distribution module.

    This metric is the coefficient of variation of the read length
    distribution and is calculated as the standard deviation of the
    read length distribution divided by the mean of the read length
    distribution

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module

    Outputs:
        variation_metric (float): The coefficient of variation of the read
    """
    rld_df = pd.DataFrame.from_dict(rld_dict, orient="index")
    rld_df = rld_df.reset_index()
    rld_df.columns = pd.Index(["read_length", "read_count"])

    mean = (rld_df["read_length"] * rld_df["read_count"]).sum()\
        / rld_df["read_count"].sum()
    variance = (
        (rld_df["read_length"] - mean)**2 * rld_df["read_count"]
        ).sum()\
        / rld_df["read_count"].sum()
    return math.sqrt(variance) / mean


def read_length_distribution_bimodality(data):
    """
    Calculate the bimodality coefficient for a given dataset.

    Args:
        data (dict): A dictionary containing the read length distribution.

    Returns:
        float: The bimodality coefficient.
    """
    read_lens = np.array(list(data.keys()))
    counts = np.array(list(data.values()))

    data = np.repeat(read_lens, counts)
    n = len(data)
    skew_value = skew(data)
    kurt_value = kurtosis(data)

    numerator = (skew_value ** 2) + 1
    denominator = kurt_value + (3 * ((n - 1) ** 2 / ((n - 2) * (n - 3))))
    bimodality_coeff = numerator / denominator

    return bimodality_coeff


def read_length_distribution_normality_metric(
        rld_dict: dict,
        ) -> float:
    """
    Calculate the read length distribution (non) normality metric from the output of
    the read_length_distribution module.

    This metric is normaltest statistic of the read length distribution

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module

    Outputs:
        non_normality_metric (float): The normaltest statistic of the read
    """
    read_lens = np.array(list(rld_dict.keys()))
    counts = np.array(list(rld_dict.values()))

    data = np.repeat(read_lens, counts)
    res = normaltest(data)
    return res.pvalue


def read_length_distribution_max_prop_metric(
        rld_dict: dict,
        num_top_readlens: int = 1,
        ) -> float:
    """
    Calculate the proportion of reads in the most frequent read length

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module
        num_top_readlens: The number of top read lengths to consider

    Outputs:
        prop_at_peak (float): The proportion of reads in the most frequent
    """
    max_count = sum(sorted(rld_dict.values(), reverse=True)[:num_top_readlens])
    total_count = sum(rld_dict.values())

    return max_count / total_count


def terminal_nucleotide_bias_KL_metric(
        observed_freq: dict,
        expected_freq: dict,
        prime: str = "five_prime",
        ) -> float:
    """
    Calculate the ligation bias metric from the output of
    the terminal_nucleotide_bias_distribution module.

    This metric is the K-L divergence of the ligation bias distribution
    of the observed frequencies from the expected frequencies. The
    expected frequencies are calculated from the nucleotide composition
    of the genome.

    Inputs:
        observed_freq: Dictionary containing the output of the
                terminal_nucleotide_bias_distribution module
        expected_freq: Dictionary containing the expected frequencies
        prime: The prime end to consider

    Outputs:
        lbd_df: Dataframe containing the ligation bias metric in bits
    """
    # Needs possible rewrite using normalised ligation bias.
    # Current iteration only accounts for five_prime
    # division by 0 if background is non-existent, Only patterns that occur
    # at least once are used (needs to be changed in ligation bias)
    kl_divergence = 0.0

    for dinucleotide, observed_prob in observed_freq[prime].items():
        expected_prob = expected_freq[dinucleotide]
        kl_divergence += observed_prob * math.log2(
                                            observed_prob / expected_prob
                                            )
    return -kl_divergence


def terminal_nucleotide_bias_max_absolute_metric(
        observed_freq: dict,
        expected_freq: dict,
        prime: str = "five_prime",
        ) -> float:
    """
    Calculate the ligation bias metric from the output of
    the terminal_nucleotide_bias_distribution module.

    This metric is the maximum difference in observed and expected
    frequencies of dinucleotides

    Inputs:
        observed_freq: Dictionary containing the output of the
                terminal_nucleotide_bias_distribution module
        expected_freq: Dictionary containing the expected frequencies

    Outputs:
        lbd_df: Dataframe containing the ligation bias metric in bits
    """
    scores = {}
    for dinucleotide, observed_prob in observed_freq[prime].items():
        expected_prob = expected_freq[dinucleotide]
        scores[dinucleotide] = abs(
            observed_prob - expected_prob)

    return max(scores.values())


def cds_coverage_metric(
        cds_read_df: pd.DataFrame,
        minimum_reads: int = 1,
        in_frame_coverage: bool = True,
        num_transcripts: int = 100,
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

    top_transcripts = cds_coverage_df[
        "transcript_id"
        ].value_counts().index[:num_transcripts]
    cds_coverage_df = cds_coverage_df[
        cds_coverage_df["transcript_id"].isin(top_transcripts)]

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
    Apply the cut off to the information content metric and calculate a global score

    Inputs:
        frame_info_content_dict: Dictionary containing the output of the
                information_metric_cutoff module
        min_count_threshold: Minimum count threshold for a read length to be
                included in the metric

    Outputs:
        information_content_metric: Dictionary containing the information
                content metric for each read length and a global score
    """
    information_content_metric = {}
    total_reads = sum(
        frame_info_content_dict[key][1]
        for key in frame_info_content_dict
    )
    total_weighted_score = 0
    total_count_above_threshold = 0

    for read_length in frame_info_content_dict:
        score, count = frame_info_content_dict[read_length]
        if count > total_reads * min_count_threshold:
            information_content_metric[read_length] = score
            total_weighted_score += score * count
            total_count_above_threshold += count

    # Calculate global score
    if total_count_above_threshold > 0:
        global_score = total_weighted_score / total_count_above_threshold
    else:
        global_score = 0

    # Add global score to the output dictionary
    information_content_metric['global'] = global_score

    return information_content_metric


def read_frame_information_weighted_score(
    frame_info_content_dict: dict,
        ):
    '''
    Produce a single metric for the triplet periodicity by taking the weighted
    average of the scores for each read length.

    Inputs:
        frame_info_content_dict (dict): Dictionary containing the information
            content metric and total counts for each read length

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


def region_region_ratio_metric(
    mRNA_distribution: dict,
    region1: str = "leader",
    region2: str = "CDS",
    read_length_range: tuple = (20, 40),
    ) -> dict:
    """
    Calculate the region-region ratio metric. This metric is the ratio of
    reads in region1 relative to region2.

    Inputs:
    mRNA_distribution: Dictionary containing the output of the
        mRNA_distribution module
    region1: String specifying the first region
    region2: String specifying the second region
    read_length_range: Tuple containing the minimum and maximum read
        length to consider for the metric

    Outputs:
    region_region_ratio: Dictionary containing the region-region ratio metric
    """
    region_region_ratio: Dict[str, float] = {}
    region1_total, region2_total = 0, 0
    read_lengths = [i for i in range(
        read_length_range[0], read_length_range[1]
    )]
    for read_len in mRNA_distribution:
        if read_len in read_lengths:
            region1_total += mRNA_distribution[read_len][region1]
            region2_total += mRNA_distribution[read_len][region2]
            if mRNA_distribution[read_len][region2] == 0:
                region_region_ratio[read_len] = 0
            else:
                region_region_ratio[read_len] = mRNA_distribution[
                    read_len][region1] / mRNA_distribution[read_len][region2]

    if region2_total == 0:
        region_region_ratio["global"] = 0
    else:
        region_region_ratio["global"] = region1_total / region2_total
    return region_region_ratio


def proportion_of_reads_in_region(
        mRNA_distribution: dict,
        region: str = "CDS",
        ):
    """
    Calculate the proportion of reads in a specific region

    Inputs:
    mRNA_distribution: Dictionary containing the output of the
        mRNA_distribution module
    region: String specifying the region

    Outputs:
    proportion: Dictionary containing the proportion of reads in the region
    """
    proportion = {}
    total = 0
    read_len_total = {}
    for read_len in mRNA_distribution:
        read_len_total[read_len] = sum(mRNA_distribution[read_len].values())
        total += read_len_total[read_len]

    for read_len in mRNA_distribution:
        proportion[read_len] = mRNA_distribution[read_len][region] / read_len_total[read_len]
    proportion["global"] = sum(
        mRNA_distribution[read_len][region] for read_len in mRNA_distribution
    ) / total
    return proportion


def autocorrelate(signal: np.array) -> np.ndarray:
    """
    Computes the autocorrelation of a signal

    Inputs:
        signal: np.array
            The signal to compute the autocorrelation of.

    Returns:
        correlation_score: float
            The autocorrelation scores for all lags
    """
    np.seterr(divide='ignore', invalid='ignore')  # ignore divide by zero here
    autocorr = np.correlate(signal, signal, mode='full')
    autocorr = autocorr[len(signal)-1:].astype(float)
    autocorr /= autocorr[0]
    np.seterr(divide='warn', invalid='warn')  # reset to default
    return autocorr


def autocorrelate_counts(
        metagene_profile: dict,
        mode: str = "uniformity",
        lag: int = 0
        ) -> dict:
    """
    Computes the autocorrelation of the ribosome counts at a given lag.

    Parameters:
    -----------
    metagene_profile: dict
        The metagene profile to compute the autocorrelation of.

    lag: int
        The lag to compute the autocorrelation at.

    Returns:
    --------
    read_length_scores: dict
        The autocorrelation scores at the given lag.
    """
    read_length_scores = {}
    global_counts = []

    for read_length in metagene_profile:
        if not global_counts:
            global_counts = list(metagene_profile[read_length].values())
        else:
            global_counts = [
                i + j for i, j in zip(
                    global_counts,
                    list(metagene_profile[read_length].values())
                    )
                    ]
        counts = list(metagene_profile[read_length].values())

        if counts[0] is not None:
            if mode == "uniformity":
                triplet_counts = [
                    sum(counts[i:i+3]) for i in range(0, len(counts), 3)
                    ]
                count_list = np.array(triplet_counts)
                auto_correlation = autocorrelate(count_list)
                read_length_scores[read_length] = float(
                    auto_correlation[:5].mean()
                )
            elif mode == "periodicity":
                count_list = np.array(counts)
                auto_correlation = autocorrelate(count_list)
                read_length_scores[read_length] = auto_correlation[lag]
        else:
            read_length_scores[read_length] = 0

    if mode == "uniformity":
        global_counts = [
            sum(global_counts[i:i+3]) for i in range(0, len(global_counts), 3)
            ]
        global_count_list = np.array(global_counts)
        global_auto_correlation = autocorrelate(global_count_list)
        read_length_scores['global'] = global_auto_correlation[:5].mean()
    elif mode == "periodicity":
        global_auto_correlation = autocorrelate(np.array(global_counts))
        read_length_scores['global'] = (
            global_auto_correlation[lag] - global_auto_correlation.mean()
            ) / global_auto_correlation.mean()
    else:
        read_length_scores['global'] = 0
    return read_length_scores


def periodicity_autocorrelation(metagene_profile: dict, lag: int = 3) -> dict:
    """
    Computes the autocorrelation of the ribosome counts at a given lag.

    Inputs:
        metagene_profile: dict
            The metagene profile to compute the autocorrelation of.

        lag: int
            The lag to compute the autocorrelation at.

    Returns:
        read_length_scores: dict
            The autocorrelation scores at the given lag.
    """
    for read_len in metagene_profile['start']:
        for i in range(0, int(max(metagene_profile['start'][read_len]))):
            if i not in metagene_profile['start'][read_len]:
                metagene_profile['start'][read_len][i] = 0
    return autocorrelate_counts(metagene_profile['start'], mode="periodicity", lag=lag)


def uniformity_autocorrelation(metagene_profile: dict, lag: int = 3) -> dict:
    """
    Computes the autocorrelation of the ribosome counts at a given lag.

    Inputs:
        metagene_profile: dict
            The metagene profile to compute the autocorrelation of.

        lag: int
            The lag to compute the autocorrelation at.

    Returns:
        read_length_scores: dict
            The autocorrelation scores at the given lag.
    """
    for read_len in metagene_profile['start']:
        for i in range(0, int(max(metagene_profile['start'][read_len]))):
            if i not in metagene_profile['start'][read_len]:
                metagene_profile['start'][read_len][i] = 0
    return autocorrelate_counts(metagene_profile['start'], mode="uniformity")


def uniformity_entropy(metagene_profile: dict) -> dict:
    """
    Computes the uniformity of the metagene profile. Inspired by ORQAS

    Inputs:
        metagene_profile: dict
            The metagene profile to compute the uniformity of.

    Returns:
        read_length_scores: dict
            The uniformity scores for each read length.
    """
    read_len_uniformity = {}

    global_counts = []
    for read_len in metagene_profile['start']:
        if not global_counts:
            global_counts = list(metagene_profile['start'][read_len].values())
        else:
            global_counts = [
                i + j for i, j in zip(
                    global_counts,
                    list(metagene_profile['start'][read_len].values())
                    )
                    ]
        counts = list(metagene_profile['start'][read_len].values())
        total_counts = sum(counts)
        entropy = 0.0

        triplet_counts = [sum(counts[i:i+3]) for i in range(0, len(counts), 3)]

        for count in triplet_counts:
            if count > 0:
                probability = count / total_counts
                entropy -= probability * math.log(probability, 2)
        max_entropy = math.log(len(metagene_profile['start'][read_len]), 2)
        uniformity = entropy / max_entropy
        read_len_uniformity[read_len] = uniformity

    global_total_counts = sum(global_counts)
    global_entropy = 0.0
    global_triplet_counts = [
        sum(global_counts[i:i+3]) for i in range(0, len(global_counts), 3)
        ]
    for count in global_triplet_counts:
        if count > 0:
            probability = count / global_total_counts
            global_entropy -= probability * math.log(probability, 2)
    global_max_entropy = math.log(len(global_counts), 2)
    global_uniformity = global_entropy / global_max_entropy
    read_len_uniformity["global"] = global_uniformity
    return read_len_uniformity


def uniformity_theil_index(
        profile,
        read_lengths=[25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
        ):
    """
    Calculates the Theil index for a Ribo-Seq profile.

    Inputs:
        profile (dict): A dictionary where keys represent positions,
        and values represent counts.

    Returns:
        dict: The Theil index for the given profile.
    """
    theils = {}
    global_counts = []
    global_sum = 0
    for read_len in profile['start']:
        if read_len in read_lengths:
            if not global_counts:
                global_counts = list(profile['start'][read_len].values())
            else:
                global_counts = [
                    i + j for i, j in zip(
                        global_counts,
                        list(profile['start'][read_len].values())
                        )
                        ]
        total_sum = sum(profile['start'][read_len].values())
        global_sum += total_sum if read_len in read_lengths else 0

        theil_sum = 0

        for i in range(0, len(profile['start'][read_len]), 3):
            triplet_counts = list(profile['start'][read_len].values())[i:i+3]
            total_triplet_sum = sum(triplet_counts)

            if total_triplet_sum > 0:
                triplet_proportions = [
                    count / total_triplet_sum
                    for count in triplet_counts
                    ]
                theil_sum += sum(
                    [
                        proportion * math.log(1 / proportion)
                        for proportion in triplet_proportions
                        if proportion > 0
                        ]
                        )

        theils[read_len] = theil_sum

    global_theil_sum = 0
    for i in range(0, len(global_counts), 3):
        triplet_counts = global_counts[i:i+3]
        total_triplet_sum = sum(triplet_counts)

        if total_triplet_sum > 0:
            triplet_proportions = [
                count / total_triplet_sum
                for count in triplet_counts
                ]
            global_theil_sum += sum(
                [
                    proportion * math.log(1 / proportion)
                    for proportion in triplet_proportions
                    if proportion > 0
                    ]
                    )

    theils["global"] = 1 / (global_theil_sum + 0.0001)

    return theils


def uniformity_gini_index(profile):
    """

    Calculates the Gini index for a Ribo-Seq profile.

    See: https://kimberlyfessel.com/mathematics/applications/gini-use-cases/

    Inputs:
        profile (dict): A dictionary where keys represent positions,
        and values represent counts.

    Returns:
        dict: The Gini index for the given profile.
    """
    ginis = {}
    global_raw_counts = []

    for read_len in profile['start']:
        if not global_raw_counts:
            global_raw_counts = list(profile['start'][read_len].values())
        else:
            global_raw_counts = [
                i + j for i, j in zip(
                    global_raw_counts,
                    list(profile['start'][read_len].values())
                    )
                    ]
        counts = list(profile['start'][read_len].values())
        total_sum = sum(counts)
        if total_sum == 0:
            ginis[read_len] = 0
            continue
        triplet_counts = [sum(counts[i:i+3]) for i in range(0, len(counts), 3)]
        triplet_counts.sort()

        gini_sum = 0
        for i, count in enumerate(triplet_counts):
            gini_sum += count * ((2 * i) - len(triplet_counts) - 1)

        denominator = len(triplet_counts) * sum(triplet_counts)
        ginis[read_len] = abs(gini_sum / denominator)

    global_counts = [
        sum(counts[i:i+3]) for i in range(0, len(global_raw_counts), 3)
        ]
    global_counts.sort()

    global_gini_sum = 0
    for i, count in enumerate(global_counts):
        global_gini_sum += count * (2 * i - len(global_counts) - 1)

    global_denominator = len(global_counts) * sum(global_counts)
    if global_denominator != 0:
        ginis["global"] = abs(global_gini_sum / global_denominator)
    else:
        ginis["global"] = 0
    return ginis


def periodicity_dominance(read_frame_dict):
    """
    Calculate the read frame dominance metric from the output of
    the read_frame_distribution module.

    This metric is the proportion of reads in the dominant frame

    Inputs:
        read_frame_dict: Dictionary containing the output of the
                read_frame_distribution module

    Outputs:
        read_frame_dominance: Dictionary containing the read frame dominance
    """
    read_frame_dominance = {}
    global_total = 0
    global_max_frame = 0
    for read_length in read_frame_dict:
        total_count = sum(read_frame_dict[read_length].values())
        if total_count == 0:
            read_frame_dominance[read_length] = 0
            continue
        max_frame = max(
            read_frame_dict[read_length], key=read_frame_dict[read_length].get
            )
        read_frame_dominance[read_length] = read_frame_dict[
            read_length][max_frame] / total_count

        global_total += total_count
        global_max_frame += read_frame_dict[read_length][max_frame]

    read_frame_dominance["global"] = global_max_frame / global_total
    return read_frame_dominance


def counts_to_codon_proportions(counts: list) -> list:
    """
    Convert a list of counts to proportions of codon.
    Codons are windows of 3 nucleotides and there is no overlap between windows

    Inputs:
        counts: list
            A list of counts for each position

    Returns:
        dict: A dictionary where keys represent positions,
        and values represent codon proportions.
    """
    codon_proportions = []
    for i in range(0, len(counts), 3):
        codon_counts = counts[i:i+3]
        total_count = sum(codon_counts)
        for count in codon_counts:
            if total_count != 0:
                codon_proportions.append(count / total_count)
            else:
                codon_proportions.append(0)
    return codon_proportions


def fourier_transform(
        metagene_profile,
        read_lengths=[25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
        ):
    """
    Calculate the Fourier transform of the metagene profile and extract the
    amplitude at the expected 3nt periodicity frequency.

    Inputs:
        metagene_profile: dict
            The metagene profile to compute the Fourier transform of.
        read_lengths: list (optional)
            The list of read lengths to consider.

    Returns:
        fourier_scores: dict
            The Fourier transform scores for each read length, including the
            amplitude at the expected 3nt periodicity frequency.
    """
    fourier_scores = {}
    global_counts = []
    for read_len in read_lengths:
        if not global_counts:
            global_counts = list(metagene_profile['start'][read_len].values())
        else:
            global_counts = [i + j for i, j in zip(
                global_counts, list(metagene_profile['start'][read_len].values()))]
        counts = list(metagene_profile['start'][read_len].values())
        if len(counts) < 2:
            fourier_scores[read_len] = 0
        else:
            fourier_transform = np.fft.fft(counts)
            frequencies = np.fft.fftfreq(len(counts), 1/len(counts))

            expected_frequency = 1/3
            idx_3nt = np.argmin(np.abs(frequencies - expected_frequency))

            amplitudes = np.abs(fourier_transform) ** 2
            total_power = np.sum(amplitudes)
            triplet_power = amplitudes[idx_3nt]
            fourier_scores[read_len] = 1 - (triplet_power / total_power)

    if len(global_counts) < 2:
        fourier_scores["global"] = 0
    else:
        global_fourier_transform = np.fft.fft(global_counts)
        frequencies = np.fft.fftfreq(len(global_counts), 1/len(global_counts))

        expected_frequency = 1/3
        idx_3nt = np.argmin(np.abs(frequencies - expected_frequency))

        amplitudes = np.abs(global_fourier_transform) ** 2
        total_power = np.sum(amplitudes)
        triplet_power = amplitudes[idx_3nt]

        fourier_scores["global"] = 1 - (triplet_power / total_power)

    return fourier_scores


def multitaper(
        metagene_profile,
        read_lengths=[25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
        nperseg=8,
        noverlap=4):
    """
    Calculate the multitaper transform of the metagene profile.

    Inputs:
        metagene_profile: dict
            The metagene profile to compute the multitaper transform of.

        read_lengths: list, optional
            The list of read lengths to calculate the multitaper transform for.
            Default is [28, 29, 30, 31, 32].
        nperseg: int, optional
            The length of the segments to use for the multitaper transform.
            Default is 8.
        noverlap: int, optional
            The number of points of overlap between segments.
            Default is 4.

    Returns:
        multitaper_scores: dict
            The multitaper transform scores for each read length.

    Explanation:
        The multitaper transform is a spectral analysis technique that estimates
        the power spectrum of a signal. In this case, the metagene profile is
        transformed using the multitaper method for each specified read length.
        The resulting scores represent the maximum power in the multitaper spectrum
        for each read length, indicating the strength of periodicity in the metagene
        profile at different read lengths.
    """
    multitaper_scores = {}
    global_counts = []
    for read_len in read_lengths:
        if not global_counts:
            global_counts = list(metagene_profile['start'][read_len].values())
        else:
            global_counts = [
                i + j for i, j in zip(
                    global_counts,
                    list(metagene_profile['start'][read_len].values())
                    )
                    ]
        counts = list(metagene_profile['start'][read_len].values())
        multitaper_transform = signal.spectrogram(
                                        np.array(counts),
                                        window='hann',
                                        nperseg=nperseg,
                                        noverlap=noverlap
                                        )
        multitaper_scores[read_len] = np.max(multitaper_transform[2])

    global_multitaper_transform = signal.spectrogram(
                                        np.array(global_counts),
                                        window='hann',
                                        nperseg=nperseg,
                                        noverlap=noverlap
                                        )
    multitaper_scores["global"] = np.max(global_multitaper_transform[2])
    return multitaper_scores


def wavelet_transform(
        metagene_profile,
        read_lengths=[25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
        ):
    """
    Calculate the wavelet transform of the metagene profile.

    Inputs:
        metagene_profile: dict
            The metagene profile to compute the wavelet transform of.

    Returns:
        wavelet_scores: dict
            The wavelet transform scores for each read length.
    """
    wavelet_scores = {}
    global_counts = []
    for read_len in read_lengths:
        if not global_counts:
            global_counts = list(metagene_profile['start'][read_len].values())
        else:
            global_counts = [
                i + j for i, j in zip(
                    global_counts,
                    list(metagene_profile['start'][read_len].values())
                    )
                    ]
        counts = list(metagene_profile['start'][read_len].values())
        wavelet_transform = signal.cwt(np.array(counts), signal.ricker, [1])
        wavelet_scores[read_len] = np.max(wavelet_transform) / np.sum(np.abs(wavelet_transform))

    global_wavelet_transform = signal.cwt(np.array(global_counts), signal.ricker, [1])
    wavelet_scores["global"] = np.max(global_wavelet_transform) / np.sum(np.abs(global_wavelet_transform))
    return wavelet_scores
