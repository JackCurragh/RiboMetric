"""
This script contains the functions used to calculate individual metrics
for different aspects of ribosome profiling data. The functions are called
are called from qc and the input data comes from the output of their
respective modules

"""

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import kurtosis, normaltest, skew

# Minimum frame-assigned reads before a per-read-length dominant-frame fraction
# is reported. At n=100 the sampling standard error on a fraction near the 1/3
# random floor is ~0.047; below that the estimate is not separable from noise.
DOMINANCE_MIN_READS = 100


def find_category_by_cumulative_percentage(df: pd.DataFrame, percentage: float) -> int:
    """
    Calculate the read_length with cumulative percentages
    """
    df["cumulative_percentage"] = df["read_count"].cumsum() / df["read_count"].sum()
    read_length = df.loc[df["cumulative_percentage"] >= percentage, "read_length"].iloc[0]
    return int(read_length)


def read_length_iqr_fraction(
    rld_dict: Dict[int, int],
) -> Optional[float]:
    """
    Interquartile range of the read length distribution as a fraction of its
    10th-90th percentile range.

    Direction: LOWER IS BETTER. A tight footprint length distribution is the
    good case, so this returns the spread itself rather than ``1 - spread``.
    The score derived from it (``read_length_concentration_score``) applies the
    flip; see docs/METRIC_NAMING.md.

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module

    Outputs:
        iqr_fraction (float): IQR / (P90 - P10), or None when unavailable
    """
    if not rld_dict or sum(rld_dict.values()) <= 0:
        return None
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
        return None
    return float(inter_quartile_range / max_range)


def read_length_cv(
    rld_dict: dict,
) -> Optional[float]:
    """
    Coefficient of variation of the read length distribution: the standard
    deviation divided by the mean.

    Direction: LOWER IS BETTER, and reported as the CV itself. The previous
    ``1/(1 + CV)`` transform was monotonic but not anchored to anything, and it
    hid the direction behind a badness-shaped name.

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module

    Outputs:
        variation_metric (float): The coefficient of variation of the read
    """
    if not rld_dict or sum(rld_dict.values()) <= 0:
        return None
    rld_df = pd.DataFrame.from_dict(rld_dict, orient="index")
    rld_df = rld_df.reset_index()
    rld_df.columns = pd.Index(["read_length", "read_count"])

    mean = (rld_df["read_length"] * rld_df["read_count"]).sum() / rld_df["read_count"].sum()
    variance = ((rld_df["read_length"] - mean) ** 2 * rld_df["read_count"]).sum() / rld_df[
        "read_count"
    ].sum()
    if mean == 0:
        return None
    return float(math.sqrt(variance) / mean)


def read_length_bimodality_coefficient(data: Dict[int, int]) -> Optional[float]:
    """
    Sarle's bimodality coefficient of the read length distribution.

    Direction: LOWER IS BETTER (a unimodal footprint distribution is the good
    case), and reported as the coefficient itself rather than ``1/(1 + BC)``.

    Args:
        data (dict): A dictionary containing the read length distribution.

    Returns:
        float: The bimodality coefficient.
    """
    read_lens = np.array(list(data.keys()))
    counts = np.array(list(data.values()))

    expanded = np.repeat(read_lens, counts)
    n = len(expanded)
    # The small-sample correction divides by (n - 2)(n - 3): undefined for
    # n <= 3, and skewness is undefined for a single repeated length.
    if n < 4 or np.count_nonzero(counts) < 2:
        return None
    skew_value = skew(expanded)
    kurt_value = kurtosis(expanded)

    numerator = (skew_value**2) + 1
    denominator = kurt_value + (3 * ((n - 1) ** 2 / ((n - 2) * (n - 3))))
    if denominator == 0 or not np.isfinite(numerator) or not np.isfinite(denominator):
        return None
    return float(max(numerator / denominator, 0))


def read_length_normality_pvalue(
    rld_dict: Dict[int, int],
) -> Optional[float]:
    """
    p-value of D'Agostino-Pearson normaltest on the read length distribution.

    Reported as the p-value itself. It has no "good" direction: a Ribo-seq
    footprint distribution is not expected to be normal, so this is a shape
    diagnostic, not a quality score. The previous ``1 - p`` transform under a
    "normality" name implied both a direction and a quality judgement it does
    not carry.

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module

    Outputs:
        non_normality_metric (float): The normaltest statistic of the read
    """
    read_lens = np.array(list(rld_dict.keys()))
    counts = np.array(list(rld_dict.values()))

    expanded = np.repeat(read_lens, counts)
    # D'Agostino-Pearson needs at least 8 observations.
    if len(expanded) < 8:
        return None
    try:
        # SciPy >=1.9 returns a Result object with .pvalue
        pvalue = float(normaltest(expanded).pvalue)
    except Exception:
        # Older SciPy returns tuple (stat, pvalue)
        _, pvalue = normaltest(expanded)
        pvalue = float(pvalue)
    if not np.isfinite(pvalue):
        return None
    return max(0.0, min(1.0, pvalue))


def read_length_max_proportion(
    rld_dict: Dict[int, int],
    num_top_readlens: int = 1,
) -> Optional[float]:
    """
    Proportion of reads carried by the most frequent read length(s).

    Direction: CONTEXT. High concentration is normal for a clean monosome
    library and abnormal for a mixed one, so this is reported as a diagnostic
    with no score attached.

    Inputs:
        rld_dict: Dictionary containing the output of the
                read_length_distribution module
        num_top_readlens: The number of top read lengths to consider

    Outputs:
        prop_at_peak (float): The proportion of reads in the most frequent
    """
    max_count = sum(sorted(rld_dict.values(), reverse=True)[:num_top_readlens])
    total_count = sum(rld_dict.values())
    if total_count == 0:
        return None
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


def terminal_nucleotide_bias_max_deviation(
    observed_freq: Dict[str, Dict[str, float]],
    expected_freq: Dict[str, float] | Dict[str, Dict[str, float]],
    prime: str = "five_prime",
) -> float:
    """
    Largest absolute difference between an observed terminal dinucleotide
    frequency and its expected background frequency.

    Direction: LOWER IS BETTER, and reported as the deviation itself. This
    previously returned ``1 - max_diff`` under a key named for bias, so the
    name promised a badness and the value delivered a goodness. Scoring it as
    a rate then inverted it: a perfectly unbiased library scored 0.00 FAIL.
    See docs/METRIC_NAMING.md section 1.

    Inputs:
        observed_freq: Dictionary containing the output of the
                terminal_nucleotide_bias_distribution module
        expected_freq: Dictionary containing the expected frequencies

    Outputs:
        max_deviation (float): max ``|observed - expected|`` over dinucleotides
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

    return float(max(scores.values())) if scores else 0.0


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
    if cds_length_total <= 0:
        return None
    return float(cds_reads_count) / float(cds_length_total)


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
        global_score = None

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

    return sum(weighted_scores) / total_reads if total_reads > 0 else None


def region_region_ratio_metric(
    mRNA_distribution: Dict[int, Dict[str, int]],
    region1: str = "leader",
    region2: str = "CDS",
    read_length_range: Optional[tuple] = None,
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
    read_lengths = (
        set(mRNA_distribution)
        if read_length_range is None
        else set(range(read_length_range[0], read_length_range[1]))
    )
    for read_len in mRNA_distribution:
        if read_len in read_lengths:
            region1_total += mRNA_distribution[read_len][region1]
            region2_total += mRNA_distribution[read_len][region2]
            if mRNA_distribution[read_len][region2] == 0:
                region_region_ratio[read_len] = None
            else:
                region_region_ratio[read_len] = (
                    mRNA_distribution[read_len][region1] / mRNA_distribution[read_len][region2]
                )

    if region2_total == 0:
        region_region_ratio["global"] = None
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
            else None
        )
    total_reads = float(total)
    proportion["global"] = (
        float(sum(mRNA_distribution[rl][region] for rl in mRNA_distribution)) / total_reads
        if total_reads > 0
        else None
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
        metagene_profile: The metagene profile to compute the
            autocorrelation of.
        lag: The lag to compute the autocorrelation at.

    Returns:
        read_length_scores: The autocorrelation scores at the given lag.
    """

    def _score(values: List[float]) -> float:
        if not values or sum(values) <= 0:
            return None
        if mode == "uniformity":
            codons = np.array(
                [sum(values[i : i + 3]) for i in range(0, len(values), 3)], dtype=float
            )
            auto_correlation = autocorrelate(codons)
            # Lags 1-4. Lag 0 is 1 by construction and only inflated the mean.
            return float(auto_correlation[1:5].mean()) if len(auto_correlation) > 1 else None
        auto_correlation = autocorrelate(np.array(values, dtype=float))
        return float(auto_correlation[lag]) if lag < len(auto_correlation) else None

    # Series are read in position order and the global series sums read
    # lengths position by position. The global periodicity value is now the
    # same statistic as the per-length one; it used to be
    # (r_lag - mean r) / mean r, a different quantity on a different scale.
    read_length_scores: Dict[int | str, float] = {}
    global_series: Dict[int, float] = {}
    for read_length, profile in metagene_profile.items():
        positions = sorted(profile)
        values = [float(profile[p]) for p in positions]
        read_length_scores[read_length] = _score(values)
        for p, v in zip(positions, values):
            global_series[p] = global_series.get(p, 0.0) + v
    read_length_scores["global"] = _score([global_series[p] for p in sorted(global_series)])
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

    def _normalised_codon_entropy(values: List[float]) -> float:
        codons = [sum(values[i : i + 3]) for i in range(0, len(values), 3)]
        total = sum(codons)
        if total <= 0 or len(codons) < 2:
            return None
        entropy = -sum((c / total) * math.log2(c / total) for c in codons if c > 0)
        return entropy / math.log2(len(codons))

    # Positions are read in position order and the global profile sums read
    # lengths position by position. Both used to rely on dict insertion
    # order, which put zero-count positions last (O-024).
    read_len_uniformity: Dict[int | str, float] = {}
    global_series: Dict[int, float] = {}
    for read_len, profile in metagene_profile["start"].items():
        positions = sorted(profile)
        values = [float(profile[p]) for p in positions]
        read_len_uniformity[read_len] = _normalised_codon_entropy(values)
        for p, v in zip(positions, values):
            global_series[p] = global_series.get(p, 0.0) + v
    read_len_uniformity["global"] = _normalised_codon_entropy(
        [global_series[p] for p in sorted(global_series)]
    )
    return read_len_uniformity


def _codon_bins(series: Mapping[int, float]) -> List[float]:
    """Sum a position -> count series into consecutive 3-nt bins, in position order."""
    values = [float(series[p]) for p in sorted(series)]
    return [sum(values[i : i + 3]) for i in range(0, len(values), 3)]


def _sum_by_position(series_list: Sequence[Mapping[int, float]]) -> Dict[int, float]:
    total: Dict[int, float] = {}
    for series in series_list:
        for p, v in series.items():
            total[p] = total.get(p, 0.0) + float(v)
    return total


def _theil_t(values: List[float]) -> float:
    x = np.asarray(values, dtype=float)
    if x.size == 0 or x.sum() <= 0:
        return None
    ratio = x[x > 0] / x.mean()
    return float(np.sum(ratio * np.log(ratio)) / x.size)


def _gini(values: List[float]) -> float:
    x = sorted(float(v) for v in values)
    n, total = len(x), sum(x)
    if n == 0 or total <= 0:
        return None
    # Rank i is 1-based: G = sum_i (2i - n - 1) x_(i) / (n sum x).
    return float(sum((2 * (i + 1) - n - 1) * v for i, v in enumerate(x)) / (n * total))


def uniformity_theil_index(
    profile: Mapping[str, Dict[int, Dict[int, int]]],
    read_lengths: Optional[List[int]] = None,
) -> Dict[int | str, float]:
    """
    Theil T index of codon-binned coverage in the start-codon window.

    ``T = (1/K) sum_k (x_k / mu) ln(x_k / mu)`` over the K codon bins, in
    [0, ln K]: 0 when coverage is even across codons, ln K when every read
    falls in one codon. Lower is more uniform.

    This used to return ``1 / (1 + sum_codons H(within-codon frame
    proportions))``, which measures lack of triplet periodicity, not
    inequality of coverage, under a name and registry entry that promised a
    lower-is-better Theil index (O-032).

    Inputs:
        profile: metagene profile ({"start": {read_length: {position: count}}})
        read_lengths: read lengths to include (default: all observed)

    Returns:
        dict: Theil T per read length and ``global`` over the included lengths.
    """
    if read_lengths is None:
        read_lengths = list(profile["start"].keys())
    theils: Dict[int | str, float] = {}
    included = []
    for read_len, series in profile["start"].items():
        theils[read_len] = _theil_t(_codon_bins(series))
        if read_len in read_lengths:
            included.append(series)
    theils["global"] = _theil_t(_codon_bins(_sum_by_position(included)))
    return theils


def uniformity_gini_index(
    profile: Mapping[str, Dict[int, Dict[int, int]]],
) -> Dict[int | str, float]:
    """
    Gini coefficient of codon-binned coverage in the start-codon window.

    0 when coverage is even across the K codon bins, (K - 1)/K when every read
    falls in one codon. Lower is more uniform, as the registry states.

    This used to return ``1 - G``, computed with a 0-based rank that shifted
    the coefficient by 2/n before an ``abs()`` hid the sign (O-032).

    See: https://kimberlyfessel.com/mathematics/applications/gini-use-cases/

    Inputs:
        profile: metagene profile ({"start": {read_length: {position: count}}})

    Returns:
        dict: Gini coefficient per read length and ``global``.
    """
    ginis: Dict[int | str, float] = {}
    for read_len, series in profile["start"].items():
        ginis[read_len] = _gini(_codon_bins(series))
    ginis["global"] = _gini(_codon_bins(_sum_by_position(list(profile["start"].values()))))
    return ginis


def periodicity_dominance(
    read_frame_dict: Dict[int, Dict[int, int]],
    min_reads: int = DOMINANCE_MIN_READS,
) -> Dict[int | str, float]:
    """
    Calculate the read frame dominance metric from the output of
    the read_frame_distribution module.

    This metric is the proportion of reads in the dominant frame

    Read lengths carrying fewer than ``min_reads`` frame-assigned reads are
    omitted from the per-read-length map rather than reported. A dominant-frame
    fraction estimated from a handful of reads is noise, and it is published as
    a confident number: a read length with a single read reports a dominance of
    exactly 1.000, which then appears in cohort tables and as a full-height bar
    in the recommended-read-lengths plot. ``information_metric_cutoff`` already
    applies the same convention to the sibling periodicity metric.

    The global values are unaffected: every read still contributes to
    ``global`` and ``global_by_read_length_max``, so the gated Tier-1 metric
    keeps its existing meaning.

    Inputs:
        read_frame_dict: Dictionary containing the output of the
                read_frame_distribution module
        min_reads: Minimum frame-assigned reads for a per-read-length
                dominance value to be reported

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
            continue
        max_frame = max(read_frame_dict[read_length], key=lambda k: read_frame_dict[read_length][k])
        if total_count >= min_reads:
            read_frame_dominance[read_length] = (
                read_frame_dict[read_length][max_frame] / total_count
            )

        global_total += total_count
        global_by_read_length_max += read_frame_dict[read_length][max_frame]
        for frame, count in read_frame_dict[read_length].items():
            global_frame_counts[int(frame)] = global_frame_counts.get(int(frame), 0) + int(count)

    read_frame_dominance["global"] = (
        max(global_frame_counts.values()) / global_total if global_total > 0 else None
    )
    read_frame_dominance["global_by_read_length_max"] = (
        global_by_read_length_max / global_total if global_total > 0 else None
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


def _triplet_power_fraction(counts: List[float]) -> float:
    """Share of a series' variance carried at a period of 3 positions.

    The series is mean-centred, so the DC term carries no power, and the power
    at the frequency bins nearest +1/3 and -1/3 cycle per position is divided
    by the total. A series that repeats every 3 positions scores 1; a flat one
    has no variance and scores 0.
    """
    x = np.asarray(counts, dtype=float)
    if x.size < 3:
        return 0.0
    x = x - x.mean()
    power = np.abs(np.fft.fft(x)) ** 2
    total = float(power.sum())
    if total <= 0:
        return 0.0
    freqs = np.fft.fftfreq(x.size)  # cycles per position
    pos = int(np.argmin(np.abs(freqs - 1 / 3)))
    neg = int(np.argmin(np.abs(freqs + 1 / 3)))
    triplet = float(power[pos]) + (float(power[neg]) if neg != pos else 0.0)
    return triplet / total


def fourier_transform(
    metagene_profile: Dict[str, Dict[int, Dict[int, int]]],
    read_lengths: Optional[List[int]] = None,
) -> Dict[int | str, float]:
    """
    Fraction of each metagene series' variance at a period of 3 nt.

    For every read length, and for the position-wise sum over read lengths
    (``global``), the start-codon metagene is read in position order and
    scored by :func:`_triplet_power_fraction`: 1 means all variation is
    triplet, 0 means none.

    The frequency grid used to be ``np.fft.fftfreq(n, 1/n)``, whose bins are
    whole numbers of cycles per series, so the bin "nearest 1/3" was always
    the DC term and the score was DC power over total power: a perfectly
    periodic profile scored 0.40 and a flat one 1.00 (O-001).

    Inputs:
        metagene_profile: dict
            The metagene profile ({"start": {read_length: {position: count}}}).
        read_lengths: list (optional)
            Read lengths to consider (default: every observed read length).

    Returns:
        fourier_scores: dict
            Triplet power fraction per read length and ``global``.
    """
    fourier_scores: Dict[int | str, float] = {}
    profiles = metagene_profile["start"]
    present = list(profiles.keys())
    # Default to every observed read length rather than a hardcoded human
    # monosome window (25-35 nt), which silently excluded sub-codon footprints
    # and non-human / alternative-nuclease libraries.
    read_lengths = present if read_lengths is None else [i for i in read_lengths if i in present]
    if not read_lengths:
        read_lengths = present

    for read_len in read_lengths:
        series = profiles[read_len]
        fourier_scores[read_len] = _triplet_power_fraction(
            [float(series[p]) for p in sorted(series)]
        )
    global_series = _sum_by_position([profiles[read_len] for read_len in read_lengths])
    fourier_scores["global"] = _triplet_power_fraction(
        [global_series[p] for p in sorted(global_series)]
    )
    return fourier_scores


def recommend_read_lengths(
    read_frame_distribution: Dict[int, Dict[int, int]],
    read_length_distribution: Dict[int, int],
    offsets: Optional[Dict] = None,
    min_periodicity: float = 0.5,
    min_read_proportion: float = 0.05,
    min_frame_reads: int = DOMINANCE_MIN_READS,
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
        # Skip read lengths with too few frame-assigned reads to estimate a
        # dominant-frame fraction. Without this a length backed by one read
        # reports periodicity 1.0 and is drawn as a full-height bar.
        if frame_total < min_frame_reads:
            continue
        periodicity = max(frames.values()) / frame_total
        proportion = read_length_distribution.get(rl, 0) / total_reads
        recommended = periodicity >= min_periodicity and proportion >= min_read_proportion
        entry: Dict[str, Any] = {
            "periodicity": round(periodicity, 4),
            "read_proportion": round(proportion, 4),
            "n_frame_reads": int(frame_total),
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
    periodicity: Optional[float],
    prop_reads_cds: Optional[float],
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
    if (
        periodicity is None
        or prop_reads_cds is None
        or periodicity < min_periodicity
        or prop_reads_cds < min_cds_proportion
    ):
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
            "periodicity": round(float(periodicity), 4) if periodicity is not None else None,
            "prop_reads_CDS": (
                round(float(prop_reads_cds), 4) if prop_reads_cds is not None else None
            ),
            "start_codon_enrichment_ratio": (
                round(float(start_codon_enrichment_ratio), 4)
                if start_codon_enrichment_ratio is not None
                else None
            ),
        },
    }


def cds_enrichment_ratio(annotated_read_df: pd.DataFrame) -> Optional[float]:
    """CDS enrichment ratio E = observed_CDS_fraction / expected_CDS_fraction.

    Per docs/METRIC_CONTRACT.md section 3.2 (METRICS_DESIGN.md R-O1). Both
    fractions are taken over the same reads -- those on *eligible* transcripts
    (transcript_length > 0 and a CDS body of at least one nucleotide) -- and
    both are weighted by the collapse ``count``:

        observed = sum_t C_t / sum_t W_t
        expected = sum_t W_t * (body_t / length_t) / sum_t W_t

    where W_t is the weight of reads on transcript t and C_t the weight of its
    CDS-body reads. Uniform-per-nucleotide coverage gives E = 1.

    Previously ``observed`` counted rows rather than weights and was taken over
    every annotated read, including transcripts left out of ``expected``, so
    the ratio compared two different read sets (O-031).

    Returns None if required columns are absent, or no eligible transcript
    has reads.
    """
    required = {"transcript_id", "mRNA_category", "cds_start", "cds_end", "transcript_length"}
    if not required.issubset(annotated_read_df.columns):
        return None
    if len(annotated_read_df) == 0:
        return None

    df = annotated_read_df
    weights = (
        pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(float)
        if "count" in df.columns
        else pd.Series(1.0, index=df.index)
    )

    tx_ann = (
        df[["transcript_id", "cds_start", "cds_end", "transcript_length"]]
        .drop_duplicates("transcript_id")
        .set_index("transcript_id")
    )
    tx_ann = tx_ann[tx_ann["transcript_length"] > 0].copy()
    tx_ann["cds_body_len"] = (tx_ann["cds_end"] - tx_ann["cds_start"] - 1).clip(lower=0)
    tx_ann = tx_ann[tx_ann["cds_body_len"] > 0]

    eligible = df["transcript_id"].isin(tx_ann.index).to_numpy()
    if not eligible.any():
        return None
    el_weights = weights[eligible]
    total_weight = float(el_weights.sum())
    if total_weight <= 0:
        return None

    is_cds = (df.loc[eligible, "mRNA_category"] == "CDS").to_numpy()
    observed_fraction = float(el_weights[is_cds].sum()) / total_weight

    tx_weight = el_weights.groupby(df.loc[eligible, "transcript_id"].astype(str).to_numpy()).sum()
    tx_weight = tx_weight[tx_weight > 0]
    tx_ann.index = tx_ann.index.astype(str)
    body_share = (
        tx_ann.loc[tx_weight.index, "cds_body_len"]
        / tx_ann.loc[tx_weight.index, "transcript_length"]
    )
    expected_fraction = float((tx_weight * body_share).sum() / tx_weight.sum())
    if expected_fraction <= 0:
        return None

    return float(observed_fraction / expected_fraction)
