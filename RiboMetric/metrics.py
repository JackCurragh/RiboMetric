"""
This script contains the functions used to calculate individual metrics
for different aspects of ribosome profiling data. The functions are called
are called from qc and the input data comes from the output of their
respective modules

"""

import math
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import kurtosis, normaltest, skew


def find_category_by_cumulative_percentage(df: pd.DataFrame, percentage: float) -> int:
    """
    Calculate the read_length with cumulative percentages
    """
    df["cumulative_percentage"] = df["read_count"].cumsum() / df["read_count"].sum()
    read_length = df.loc[df["cumulative_percentage"] >= percentage, "read_length"].iloc[0]
    return int(read_length)


def read_length_distribution_IQR_normalised_metric(
    rld_dict: Dict[int, int],
) -> float:
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

    max_range = find_category_by_cumulative_percentage(
        rld_df, 0.9
    ) - find_category_by_cumulative_percentage(rld_df, 0.1)

    if max_range == 0:
        return 1.0  # Perfect concentration, IQR = 0
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

    mean = (rld_df["read_length"] * rld_df["read_count"]).sum() / rld_df["read_count"].sum()
    variance = ((rld_df["read_length"] - mean) ** 2 * rld_df["read_count"]).sum() / rld_df[
        "read_count"
    ].sum()
    coefficient_of_variation = math.sqrt(variance) / mean if mean != 0 else 0
    # Transform so that higher scores indicate tighter (better) distributions.
    return 1 / (1 + coefficient_of_variation)


def read_length_distribution_bimodality(data: Dict[int, int]) -> float:
    """
    Calculate the bimodality coefficient for a given dataset.

    Args:
        data (dict): A dictionary containing the read length distribution.

    Returns:
        float: The bimodality coefficient.
    """
    read_lens = np.array(list(data.keys()))
    counts = np.array(list(data.values()))

    expanded = np.repeat(read_lens, counts)
    n = len(expanded)
    skew_value = skew(expanded)
    kurt_value = kurtosis(expanded)

    numerator = (skew_value**2) + 1
    denominator = kurt_value + (3 * ((n - 1) ** 2 / ((n - 2) * (n - 3))))
    bimodality_coeff = numerator / denominator if denominator != 0 else 0

    # Lower bimodality is preferred; map to (0,1] where 1 is best (unimodal).
    bimodality_coeff = max(bimodality_coeff, 0)
    return 1 / (1 + bimodality_coeff)


def read_length_distribution_normality_metric(
    rld_dict: Dict[int, int],
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

    expanded = np.repeat(read_lens, counts)
    try:
        # SciPy >=1.9 returns a Result object with .pvalue
        pvalue = float(normaltest(expanded).pvalue)
    except Exception:
        # Older SciPy returns tuple (stat, pvalue)
        _, pvalue = normaltest(expanded)
        pvalue = float(pvalue)
    # Non-normal distributions are expected; invert p-value so higher is better.
    return max(0.0, min(1.0, 1 - pvalue))


def read_length_distribution_max_prop_metric(
    rld_dict: Dict[int, int],
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


def terminal_nucleotide_bias_KL_divergence(
    observed_freq: Dict[str, Dict[str, float]],
    expected_freq: Dict[str, float] | Dict[str, Dict[str, float]],
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
        kl_divergence: Raw Kullback-Leibler divergence in bits.
    """
    # Needs possible rewrite using normalised ligation bias.
    # Current iteration only accounts for five_prime
    # division by 0 if background is non-existent, Only patterns that occur
    # at least once are used (needs to be changed in ligation bias)
    kl_divergence = 0.0

    exp_map: Dict[str, float]
    ref = expected_freq.get(prime) if isinstance(expected_freq, dict) else None
    if isinstance(ref, dict):
        exp_map = ref
    else:
        exp_map = expected_freq  # type: ignore[assignment]

    for dinucleotide, observed_prob in observed_freq[prime].items():
        expected_prob = exp_map[dinucleotide]
        if observed_prob <= 0:
            continue
        if expected_prob <= 0:
            continue
        kl_divergence += observed_prob * math.log2(observed_prob / expected_prob)
    return max(0.0, kl_divergence)


def terminal_nucleotide_bias_KL_metric(
    observed_freq: Dict[str, Dict[str, float]],
    expected_freq: Dict[str, float] | Dict[str, Dict[str, float]],
    prime: str = "five_prime",
) -> float:
    """
    Calculate a normalized terminal nucleotide bias score from raw KL divergence.

    The returned value is a goodness score, not the raw divergence: 1.0 means
    observed terminal dinucleotide frequencies match the expected background,
    and lower values indicate stronger bias.
    """
    kl_divergence = terminal_nucleotide_bias_KL_divergence(
        observed_freq,
        expected_freq,
        prime=prime,
    )
    # Higher values signal less bias; map divergence to (0,1].
    return 1 / (1 + kl_divergence)


def terminal_nucleotide_bias_max_absolute_metric(
    observed_freq: Dict[str, Dict[str, float]],
    expected_freq: Dict[str, float] | Dict[str, Dict[str, float]],
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
    exp_map2: Dict[str, float]
    ref2 = expected_freq.get(prime) if isinstance(expected_freq, dict) else None
    if isinstance(ref2, dict):
        exp_map2 = ref2
    else:
        exp_map2 = expected_freq  # type: ignore[assignment]

    for dinucleotide, observed_prob in observed_freq[prime].items():
        expected_prob = exp_map2[dinucleotide]
        scores[dinucleotide] = abs(observed_prob - expected_prob)

    max_diff = max(scores.values()) if scores else 0
    # Perfect agreement should score 1, larger deviations trend towards 0.
    return max(0.0, 1 - max_diff)


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
    required_columns = ["transcript_id", "a_site", "cds_start", "cds_end"]
    optional_columns = ["count"] if "count" in cds_read_df.columns else []
    cds_coverage_df = cds_read_df[required_columns + optional_columns].copy()
    if "count" not in cds_coverage_df.columns:
        cds_coverage_df["count"] = 1
    cds_coverage_df["count"] = (
        pd.to_numeric(cds_coverage_df["count"], errors="coerce").fillna(0).astype(int)
    )
    # Build a stable string key "transcript_id_aSite" using Python strings to
    # avoid Arrow-backed string arithmetic issues in some pandas builds.
    a_site_num = (
        pd.to_numeric(cds_coverage_df["a_site"], errors="coerce").fillna(-1).astype("int64")
    )
    cds_coverage_df["name_pos"] = (
        cds_coverage_df["transcript_id"].astype("object").astype(str) + "_" + a_site_num.astype(str)
    ).astype("category")

    top_transcripts = (
        cds_coverage_df.groupby("transcript_id", observed=True)["count"]
        .sum()
        .sort_values(ascending=False)
        .index[:num_transcripts]
    )
    cds_coverage_df = cds_coverage_df[cds_coverage_df["transcript_id"].isin(top_transcripts)]

    # Calculate the total combined length of the CDS of transcripts that have
    # reads aligned to them
    cds_transcripts = cds_coverage_df.drop_duplicates("transcript_id").copy()
    interior_cds_length = (cds_transcripts["cds_end"] - cds_transcripts["cds_start"] - 1).clip(
        lower=0
    )
    cds_length_total = int(interior_cds_length.sum())
    if in_frame_coverage:
        cds_length_total = int((interior_cds_length // 3).sum())
    del cds_transcripts

    # If in_frame_coverage is true, take only reads that are in frame for
    # their transcript.
    if in_frame_coverage:
        cds_coverage_df = cds_coverage_df[
            (cds_coverage_df["a_site"] - cds_coverage_df["cds_start"]) % 3 == 0
        ]

    # Calculate the count of nucleotides covered by the reads after filtering
    # Weighted “coverage”: sum positions whose weighted count exceeds threshold
    # Ensure numeric weights for aggregation
    pos_counts = cds_coverage_df.groupby("name_pos", observed=True)["count"].sum()
    cds_reads_count = int((pos_counts >= minimum_reads).sum())
    denom = float(cds_length_total) if float(cds_length_total) != 0 else 1.0
    return float(cds_reads_count) / denom


def calculate_3nt_periodicity_score(probabilities: List[float]) -> float:
    """
    Calculate the triplet periodicity score for a given probability of a read
    being in frame. The score is the normalised entropy reduction of the
    triplet distribution: 0 = random (maximum entropy), 1 = perfect
    single-frame occupancy (zero entropy).

    Numerator is the Maximum Entropy of the triplet distribution minus the
    entropy of the triplet distribution.
    Denominator is the Maximum Entropy of the triplet distribution.

    The previous implementation took the square root of this ratio, which
    inflated middling values without adding interpretability (see
    docs/METRICS_DESIGN.md, periodicity_information). The ratio is already
    naturally anchored, so it is returned directly.

    Inputs:
        probability (float): The probability of a read being in frame.

    Returns:
        result (float): The triplet periodicity score.
    """
    maximum_entropy = math.log2(3)
    entropy = 0.0
    for probability in probabilities:
        entropy += -(probability * math.log2(probability))

    result = (maximum_entropy - entropy) / maximum_entropy
    return float(result)


def read_frame_information_content(
    read_frame_distribution: Dict[int, Dict[int, int]],
) -> Dict[int, Tuple[float, int]]:
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
    frame_info_content_dict: Dict[int, Tuple[float, int]] = {}
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
    frame_info_content_dict: Dict[int, Tuple[float, int]],
    min_count_threshold: float = 0.05,
) -> Dict[int | str, float]:
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
    information_content_metric: Dict[int | str, float] = {}
    total_reads = sum(frame_info_content_dict[key][1] for key in frame_info_content_dict)
    total_weighted_score: float = 0.0
    total_count_above_threshold: float = 0.0

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
        global_score = 0.0

    # Add global score to the output dictionary
    information_content_metric["global"] = global_score

    return information_content_metric


def read_frame_information_weighted_score(
    frame_info_content_dict: Dict[int, Tuple[float, int]],
) -> float:
    """
    Produce a single metric for the triplet periodicity by taking the weighted
    average of the scores for each read length.

    Inputs:
        frame_info_content_dict (dict): Dictionary containing the information
            content metric and total counts for each read length

    Returns:
        result (float): The triplet periodicity score.
    """
    total_reads = sum(frame_info_content_dict[key][1] for key in frame_info_content_dict)
    weighted_scores = []
    for _, score in frame_info_content_dict.items():
        weighted_score = score[0] * score[1]
        weighted_scores.append(weighted_score)

    return sum(weighted_scores) / total_reads if total_reads > 0 else 0


def region_region_ratio_metric(
    mRNA_distribution: Dict[int, Dict[str, int]],
    region1: str = "leader",
    region2: str = "CDS",
    read_length_range: tuple = (20, 40),
) -> Dict[int | str, float]:
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
    region_region_ratio: Dict[int | str, float] = {}
    region1_total, region2_total = 0, 0
    read_lengths = [i for i in range(read_length_range[0], read_length_range[1])]
    for read_len in mRNA_distribution:
        if read_len in read_lengths:
            region1_total += mRNA_distribution[read_len][region1]
            region2_total += mRNA_distribution[read_len][region2]
            if mRNA_distribution[read_len][region2] == 0:
                region_region_ratio[read_len] = 0
            else:
                region_region_ratio[read_len] = (
                    mRNA_distribution[read_len][region1] / mRNA_distribution[read_len][region2]
                )

    if region2_total == 0:
        region_region_ratio["global"] = 0
    else:
        region_region_ratio["global"] = region1_total / region2_total
    return region_region_ratio


def proportion_of_reads_in_region(
    mRNA_distribution: Dict[int, Dict[str, int]],
    region: str = "CDS",
) -> Dict[int | str, float]:
    """
    Calculate the proportion of reads in a specific region

    Inputs:
    mRNA_distribution: Dictionary containing the output of the
        mRNA_distribution module
    region: String specifying the region

    Outputs:
    proportion: Dictionary containing the proportion of reads in the region
    """
    proportion: Dict[int | str, float] = {}
    total = 0
    read_len_total: Dict[int, int] = {}
    for read_len in mRNA_distribution:
        read_len_total[read_len] = sum(mRNA_distribution[read_len].values())
        total += read_len_total[read_len]

    for read_len in mRNA_distribution:
        proportion[read_len] = (
            mRNA_distribution[read_len][region] / read_len_total[read_len]
            if read_len_total[read_len] > 0
            else 0.0
        )
    total_reads = float(total)
    proportion["global"] = (
        float(sum(mRNA_distribution[rl][region] for rl in mRNA_distribution)) / total_reads
        if total_reads > 0
        else 0.0
    )
    return proportion


def autocorrelate(signal: npt.NDArray[np.floating[Any]]) -> npt.NDArray[np.float64]:
    """
    Computes the autocorrelation of a signal

    Inputs:
        signal: np.array
            The signal to compute the autocorrelation of.

    Returns:
        correlation_score: float
            The autocorrelation scores for all lags
    """
    np.seterr(divide="ignore", invalid="ignore")  # ignore divide by zero here
    autocorr = np.correlate(signal, signal, mode="full")
    autocorr = autocorr[len(signal) - 1 :].astype(float)
    autocorr /= autocorr[0]
    np.seterr(divide="warn", invalid="warn")  # reset to default
    return autocorr.astype(np.float64)


def autocorrelate_counts(
    metagene_profile: Dict[int, Dict[int, int]], mode: str = "uniformity", lag: int = 0
) -> Dict[int | str, float]:
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
    read_length_scores: Dict[int | str, float] = {}
    global_counts: npt.NDArray[np.float64] | None = None

    for read_length in metagene_profile:
        counts = list(metagene_profile[read_length].values())
        counts_arr = np.array(counts, dtype=float)
        if global_counts is None:
            global_counts = counts_arr.copy()
        else:
            global_counts = global_counts + counts_arr

        if counts[0] is not None and sum(counts) > 0:
            if mode == "uniformity":
                triplet_counts = [sum(counts[i : i + 3]) for i in range(0, len(counts), 3)]
                count_list = np.array(triplet_counts, dtype=float)
                auto_correlation = autocorrelate(count_list)
                read_length_scores[read_length] = float(auto_correlation[:5].mean())
            elif mode == "periodicity":
                count_list = np.array(counts, dtype=float)
                auto_correlation = autocorrelate(count_list)
                read_length_scores[read_length] = auto_correlation[lag]
        else:
            read_length_scores[read_length] = 0

    if global_counts is None:
        read_length_scores["global"] = 0
        return read_length_scores

    if mode == "uniformity" and global_counts.sum() > 0:
        triplet_arr = np.array(
            [global_counts[i : i + 3].sum() for i in range(0, len(global_counts), 3)], dtype=float
        )
        global_auto_correlation = autocorrelate(triplet_arr)
        read_length_scores["global"] = float(global_auto_correlation[:5].mean())
    elif mode == "periodicity" and global_counts.sum() > 0:
        global_auto_correlation = autocorrelate(global_counts)
        read_length_scores["global"] = (
            global_auto_correlation[lag] - global_auto_correlation.mean()
        ) / global_auto_correlation.mean()
    else:
        read_length_scores["global"] = 0
    return read_length_scores


def periodicity_autocorrelation(
    metagene_profile: Dict[str, Dict[int, Dict[int, int]]], lag: int = 3
) -> Dict[int | str, float]:
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
    # Build arrays over actual observed distance ranges per read length
    return autocorrelate_counts(metagene_profile["start"], mode="periodicity", lag=lag)


def uniformity_autocorrelation(
    metagene_profile: Dict[str, Dict[int, Dict[int, int]]], lag: int = 3
) -> Dict[int | str, float]:
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
    # Build arrays over actual observed distance ranges per read length
    return autocorrelate_counts(metagene_profile["start"], mode="uniformity")


def uniformity_entropy(
    metagene_profile: Dict[str, Dict[int, Dict[int, int]]],
) -> Dict[int | str, float]:
    """
    Computes the uniformity of the metagene profile. Inspired by ORQAS

    Inputs:
        metagene_profile: dict
            The metagene profile to compute the uniformity of.

    Returns:
        read_length_scores: dict
            The uniformity scores for each read length.
    """
    read_len_uniformity: Dict[int | str, float] = {}

    global_counts: List[int] = []
    for read_len in metagene_profile["start"]:
        if not global_counts:
            global_counts = list(metagene_profile["start"][read_len].values())
        else:
            global_counts = [
                i + j
                for i, j in zip(global_counts, list(metagene_profile["start"][read_len].values()))
            ]
        counts = list(metagene_profile["start"][read_len].values())
        total_counts = sum(counts)
        entropy = 0.0

        triplet_counts = [sum(counts[i : i + 3]) for i in range(0, len(counts), 3)]

        for count in triplet_counts:
            if count > 0:
                probability = count / total_counts
                entropy -= probability * math.log(probability, 2)
        max_entropy = math.log(len(triplet_counts), 2)
        uniformity = entropy / max_entropy if max_entropy > 0 else 0
        read_len_uniformity[read_len] = uniformity

    global_total_counts = sum(global_counts)
    global_entropy = 0.0
    global_triplet_counts = [sum(global_counts[i : i + 3]) for i in range(0, len(global_counts), 3)]
    for count in global_triplet_counts:
        if count > 0:
            probability = count / global_total_counts
            global_entropy -= probability * math.log(probability, 2)
    global_max_entropy = math.log(len(global_triplet_counts), 2)
    global_uniformity = global_entropy / global_max_entropy if global_max_entropy > 0 else 0
    read_len_uniformity["global"] = global_uniformity
    return read_len_uniformity


def uniformity_theil_index(
    profile: Mapping[str, Dict[int, Dict[int, int]]],
    read_lengths: Optional[List[int]] = None,
) -> Dict[int | str, float]:
    """
    Calculates the Theil index for a Ribo-Seq profile.

    Inputs:
        profile (dict): A dictionary where keys represent positions,
        and values represent counts.

    Returns:
        dict: The Theil index for the given profile.
    """
    theils: Dict[int | str, float] = {}
    global_counts: List[int] = []
    global_sum = 0
    # Default to all observed read lengths rather than the human 25-35 nt window.
    if read_lengths is None:
        read_lengths = list(profile["start"].keys())
    for read_len in profile["start"]:
        if read_len in read_lengths:
            if not global_counts:
                global_counts = list(profile["start"][read_len].values())
            else:
                global_counts = [
                    i + j for i, j in zip(global_counts, list(profile["start"][read_len].values()))
                ]
        read_len_counts = list(profile["start"][read_len].values())
        total_sum = sum(read_len_counts)
        global_sum += total_sum if read_len in read_lengths else 0

        theil_sum = 0.0

        for i in range(0, len(read_len_counts), 3):
            triplet_counts = read_len_counts[i : i + 3]
            total_triplet_sum = sum(triplet_counts)

            if total_triplet_sum > 0:
                triplet_proportions = [count / total_triplet_sum for count in triplet_counts]
                theil_sum += sum(
                    [
                        proportion * math.log(1 / proportion)
                        for proportion in triplet_proportions
                        if proportion > 0
                    ]
                )

        theils[read_len] = 1 / (1 + theil_sum)

    global_theil_sum: float = 0.0
    for i in range(0, len(global_counts), 3):
        triplet_counts = global_counts[i : i + 3]
        total_triplet_sum = sum(triplet_counts)

        if total_triplet_sum > 0:
            triplet_proportions = [count / total_triplet_sum for count in triplet_counts]
            global_theil_sum += sum(
                [
                    proportion * math.log(1 / proportion)
                    for proportion in triplet_proportions
                    if proportion > 0
                ]
            )

    theils["global"] = 1 / (1 + global_theil_sum)

    return theils


def uniformity_gini_index(
    profile: Mapping[str, Dict[int, Dict[int, int]]],
) -> Dict[int | str, float]:
    """

    Calculates the Gini index for a Ribo-Seq profile.

    See: https://kimberlyfessel.com/mathematics/applications/gini-use-cases/

    Inputs:
        profile (dict): A dictionary where keys represent positions,
        and values represent counts.

    Returns:
        dict: The Gini index for the given profile.
    """
    ginis: Dict[int | str, float] = {}
    global_raw_counts: List[int] = []

    for read_len in profile["start"]:
        if not global_raw_counts:
            global_raw_counts = list(profile["start"][read_len].values())
        else:
            global_raw_counts = [
                i + j for i, j in zip(global_raw_counts, list(profile["start"][read_len].values()))
            ]
        counts = list(profile["start"][read_len].values())
        total_sum = sum(counts)
        if total_sum == 0:
            ginis[read_len] = 0
            continue
        triplet_counts = [sum(counts[i : i + 3]) for i in range(0, len(counts), 3)]
        triplet_counts.sort()

        gini_sum = 0
        for i, count in enumerate(triplet_counts):
            gini_sum += count * ((2 * i) - len(triplet_counts) - 1)

        denominator = len(triplet_counts) * sum(triplet_counts)
        gini_value = abs(gini_sum / denominator) if denominator != 0 else 0
        gini_value = max(0, min(1, gini_value))
        ginis[read_len] = 1 - gini_value

    global_counts: List[int] = [
        sum(global_raw_counts[i : i + 3]) for i in range(0, len(global_raw_counts), 3)
    ]
    global_counts.sort()

    global_gini_sum = 0
    for i, count in enumerate(global_counts):
        global_gini_sum += count * (2 * i - len(global_counts) - 1)

    global_denominator = len(global_counts) * sum(global_counts)
    if global_denominator != 0:
        global_gini_value = abs(global_gini_sum / global_denominator)
    else:
        global_gini_value = 0
    global_gini_value = max(0, min(1, global_gini_value))
    ginis["global"] = 1 - global_gini_value
    return ginis


def periodicity_dominance(read_frame_dict: Dict[int, Dict[int, int]]) -> Dict[int | str, float]:
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
    read_frame_dominance: Dict[int | str, float] = {}
    global_total: int = 0
    global_by_read_length_max: int = 0
    global_frame_counts: Dict[int, int] = {0: 0, 1: 0, 2: 0}
    for read_length in read_frame_dict:
        total_count = sum(read_frame_dict[read_length].values())
        if total_count == 0:
            read_frame_dominance[read_length] = 0
            continue
        max_frame = max(read_frame_dict[read_length], key=lambda k: read_frame_dict[read_length][k])
        read_frame_dominance[read_length] = (
            read_frame_dict[read_length][max_frame] / total_count if total_count > 0 else 0
        )

        global_total += total_count
        global_by_read_length_max += read_frame_dict[read_length][max_frame]
        for frame, count in read_frame_dict[read_length].items():
            global_frame_counts[int(frame)] = global_frame_counts.get(int(frame), 0) + int(count)

    read_frame_dominance["global"] = (
        max(global_frame_counts.values()) / global_total if global_total > 0 else 0
    )
    read_frame_dominance["global_by_read_length_max"] = (
        global_by_read_length_max / global_total if global_total > 0 else 0
    )
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
        codon_counts = counts[i : i + 3]
        total_count = sum(codon_counts)
        for count in codon_counts:
            if total_count != 0:
                codon_proportions.append(count / total_count)
            else:
                codon_proportions.append(0)
    return codon_proportions


def fourier_transform(
    metagene_profile: Dict[str, Dict[int, Dict[int, int]]],
    read_lengths: Optional[List[int]] = None,
) -> Dict[int | str, float]:
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
    fourier_scores: Dict[int | str, float] = {}
    global_counts: List[int] = []
    present = list(metagene_profile["start"].keys())
    # Default to every observed read length rather than a hardcoded human
    # monosome window (25-35 nt), which silently excluded sub-codon footprints
    # and non-human / alternative-nuclease libraries.
    read_lengths = present if read_lengths is None else [i for i in read_lengths if i in present]
    if not read_lengths:
        read_lengths = present

    for read_len in read_lengths:
        series = metagene_profile["start"][read_len]
        positions = sorted(series.keys())
        counts = [series[p] for p in positions]
        if not global_counts:
            global_counts = counts.copy()
        else:
            # pad to same length if needed
            if len(global_counts) < len(counts):
                global_counts += [0] * (len(counts) - len(global_counts))
            elif len(counts) < len(global_counts):
                counts += [0] * (len(global_counts) - len(counts))
            global_counts = [i + j for i, j in zip(global_counts, counts)]
        if len(counts) < 2:
            fourier_scores[read_len] = 0
        else:
            fourier_transform = np.fft.fft(counts)
            frequencies = np.fft.fftfreq(len(counts), 1 / len(counts))

            expected_frequency = 1 / 3
            idx_3nt = np.argmin(np.abs(frequencies - expected_frequency))

            amplitudes = np.abs(fourier_transform) ** 2
            total_power = np.sum(amplitudes)
            triplet_power = amplitudes[idx_3nt]
            fourier_scores[read_len] = triplet_power / total_power if total_power != 0 else 0

    if len(global_counts) < 2:
        fourier_scores["global"] = 0
    else:
        global_fourier_transform = np.fft.fft(global_counts)
        frequencies = np.fft.fftfreq(len(global_counts), 1 / len(global_counts))

        expected_frequency = 1 / 3
        idx_3nt = np.argmin(np.abs(frequencies - expected_frequency))

        amplitudes = np.abs(global_fourier_transform) ** 2
        total_power = np.sum(amplitudes)
        triplet_power = amplitudes[idx_3nt]

        fourier_scores["global"] = triplet_power / total_power if total_power > 0 else 0

    return fourier_scores


def recommend_read_lengths(
    read_frame_distribution: Dict[int, Dict[int, int]],
    read_length_distribution: Dict[int, int],
    offsets: Optional[Dict] = None,
    min_periodicity: float = 0.5,
    min_read_proportion: float = 0.05,
) -> Dict[str, Any]:
    """Recommend the read lengths carrying clean 3-nt periodicity.

    Synthesises the per-read-length frame distribution and (when available)
    computed P-site offsets into a single actionable recommendation: which read
    lengths to keep for downstream P-site assignment / ORF calling, and the
    fraction of the library they represent. A read length is recommended when
    its dominant-frame fraction is at least ``min_periodicity`` and it carries at
    least ``min_read_proportion`` of all reads (so a length with great
    periodicity but only a handful of reads is not recommended).

    Returns a dict with per-read-length detail plus library-level summaries:
        recommended_lengths            sorted list of recommended read lengths
        n_recommended                  count of recommended read lengths
        recommended_read_proportion    fraction of all reads in those lengths
    """
    total_reads = sum(read_length_distribution.values()) or 1
    by_read_length: Dict[int, Dict[str, Any]] = {}
    for read_length, frames in read_frame_distribution.items():
        rl = int(read_length)
        frame_total = sum(frames.values())
        periodicity = max(frames.values()) / frame_total if frame_total > 0 else 0.0
        proportion = read_length_distribution.get(rl, 0) / total_reads
        recommended = periodicity >= min_periodicity and proportion >= min_read_proportion
        entry: Dict[str, Any] = {
            "periodicity": round(periodicity, 4),
            "read_proportion": round(proportion, 4),
            "recommended": bool(recommended),
        }
        if offsets is not None and rl in offsets:
            entry["offset"] = int(offsets[rl])
        elif offsets is not None and str(rl) in offsets:
            entry["offset"] = int(offsets[str(rl)])
        by_read_length[rl] = entry

    recommended_lengths = sorted(rl for rl, e in by_read_length.items() if e["recommended"])
    recommended_read_proportion = (
        sum(read_length_distribution.get(rl, 0) for rl in recommended_lengths) / total_reads
    )
    return {
        "by_read_length": by_read_length,
        "recommended_lengths": recommended_lengths,
        "n_recommended": len(recommended_lengths),
        "recommended_read_proportion": round(recommended_read_proportion, 4),
    }


def classify_library_type(
    periodicity: float,
    prop_reads_cds: float,
    start_codon_enrichment_ratio: Optional[float],
    min_periodicity: float = 0.4,
    min_cds_proportion: float = 0.5,
    initiation_ratio_cutoff: float = 3.0,
) -> Dict[str, Any]:
    """Heuristically classify the Ribo-Seq library type from summary metrics.

    Distinguishes three broad regimes that change how a dataset should be
    interpreted downstream:

    * ``low_quality``  — weak periodicity and/or little CDS enrichment; the
      footprints do not look ribosome-protected (degraded RNA / RNA-seq-like).
    * ``initiation``   — a sharp start-codon peak relative to the CDS body,
      consistent with initiation-inhibitor treatment (harringtonine / LTM) or
      heavy stalling at initiation.
    * ``elongation``   — periodic, CDS-enriched, flat start/body ratio; the
      standard elongating profile (cycloheximide / untreated).

    Returns the label plus the evidence used, so the call is auditable rather
    than a black box.
    """
    if periodicity < min_periodicity or prop_reads_cds < min_cds_proportion:
        label = "low_quality"
    elif (
        start_codon_enrichment_ratio is not None
        and start_codon_enrichment_ratio >= initiation_ratio_cutoff
    ):
        label = "initiation"
    else:
        label = "elongation"
    return {
        "label": label,
        "evidence": {
            "periodicity": round(float(periodicity), 4),
            "prop_reads_CDS": round(float(prop_reads_cds), 4),
            "start_codon_enrichment_ratio": (
                round(float(start_codon_enrichment_ratio), 4)
                if start_codon_enrichment_ratio is not None
                else None
            ),
        },
    }


def cds_enrichment_ratio(annotated_read_df: pd.DataFrame) -> Optional[float]:
    """CDS enrichment ratio E = observed_CDS_fraction / expected_CDS_fraction.

    Per METRICS_DESIGN.md §5 (R-O1). Returns None if required columns are
    absent, no reads are present, or the expected fraction cannot be computed.
    Self-consistency: uniform-per-nucleotide reads ⇒ E ≈ 1.
    """
    required = {"transcript_id", "mRNA_category", "cds_start", "cds_end", "transcript_length"}
    if not required.issubset(annotated_read_df.columns):
        return None
    if len(annotated_read_df) == 0:
        return None

    total_reads = len(annotated_read_df)
    cds_reads = int((annotated_read_df["mRNA_category"] == "CDS").sum())
    observed_fraction = cds_reads / total_reads

    # Per-transcript annotation: one row per transcript (cds_start/end/length)
    tx_ann = (
        annotated_read_df[["transcript_id", "cds_start", "cds_end", "transcript_length"]]
        .drop_duplicates("transcript_id")
        .set_index("transcript_id")
    )
    tx_reads = annotated_read_df.groupby("transcript_id", observed=True).size()

    # Eligibility: transcript_length > 0 and CDS body length > 0
    tx_ann = tx_ann[tx_ann["transcript_length"] > 0].copy()
    tx_ann["cds_body_len"] = (tx_ann["cds_end"] - tx_ann["cds_start"] - 1).clip(lower=0)
    tx_ann = tx_ann[tx_ann["cds_body_len"] > 0]

    common = tx_ann.index.intersection(tx_reads.index)
    if len(common) == 0:
        return None
    tx_ann = tx_ann.loc[common]
    tx_reads_common = tx_reads.loc[common]

    total_weighted = float(tx_reads_common.sum())
    if total_weighted == 0:
        return None

    expected_fraction = float(
        (tx_reads_common * (tx_ann["cds_body_len"] / tx_ann["transcript_length"])).sum()
        / total_weighted
    )
    if expected_fraction <= 0:
        return None

    return float(observed_fraction / expected_fraction)
