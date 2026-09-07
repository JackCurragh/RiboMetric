"""
This script contains the functions required to run individual modules
of the RibosomeProfiler pipeline

"""

import numpy as np
import pandas as pd

try:
    from xhtml2pdf import pisa

    HAS_PDF = True
except ImportError:
    HAS_PDF = False

from typing import Dict, List, Optional, Tuple, cast

DEFAULT_OFFSET_BOUNDS: Tuple[int, int] = (8, 20)
DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION = 2 / 3
OFFSET_TARGET_SHIFTS: Dict[str, int] = {"p_site": 0, "a_site": 3}


def offset_shift_for_target(offset_target: str) -> int:
    """Return the nt shift from a P-site offset to the requested target."""
    try:
        return OFFSET_TARGET_SHIFTS[str(offset_target)]
    except KeyError as exc:
        raise ValueError(
            "offset_target must be one of: " f"{', '.join(sorted(OFFSET_TARGET_SHIFTS))}"
        ) from exc


def _get_weights(df: pd.DataFrame) -> Optional[pd.Series]:
    """Return integer weights if unexpanded; otherwise None.

    Heuristic: if a 'count' column exists but rows are already expanded
    (duplicate read_name values), skip weighting to avoid double counting.
    """
    if "count" not in df.columns:
        return None
    if "read_name" in df.columns and df["read_name"].duplicated().any():
        return None
    return df["count"].astype(int)


def filter_unique_mappers(df: pd.DataFrame, enabled: bool = True) -> pd.DataFrame:
    """Return only uniquely-mapped reads when reliable metadata is present.

    Prefer the SAM NH tag when available (NH=1), because it is aligner-neutral.
    For STAR transcriptome BAMs the historical frame-sensitive behaviour was
    MAPQ=255, which means exactly one genomic locus even if multiple transcript
    rows are reported for that fragment. Keep that path for offset and
    periodicity calculations so transcript-level alternate-hit metadata does
    not discard reads that still belong to one genomic locus.

    XA/NH can still be used elsewhere for multimapper reporting, but this
    filter intentionally preserves the STAR MAPQ convention for frame-sensitive
    analyses.
    """
    if not enabled:
        return df

    if "nh" in df.columns and df["nh"].notna().any():
        return df[df["nh"].astype(float) == 1]

    if "mapq" not in df.columns:
        return df.iloc[0:0].copy()
    mapq_available = (
        df["mapq_available"].astype(bool)
        if "mapq_available" in df.columns
        else pd.Series(True, index=df.index)
    )
    return df[mapq_available & (df["mapq"] == 255)]


def filter_unique_mappers_strict(df: pd.DataFrame, enabled: bool = True) -> pd.DataFrame:
    """Return a conservative unique-mapper set for frame-sensitive analyses.

    MAPQ/NH remains the primary uniqueness signal. If a BAM still contains
    multiple primary rows with the same read name after that filter, exclude
    those names from offset and periodicity calculations rather than allowing
    ambiguous loci to drive frame-sensitive metrics.
    """
    unique_df = filter_unique_mappers(df, enabled=enabled)
    if not enabled or "read_name" not in unique_df.columns:
        return unique_df
    return unique_df[~unique_df["read_name"].duplicated(keep=False)]


def _reference_metadata(reference_names: pd.Series) -> pd.DataFrame:
    """Parse transcriptome reference metadata from STAR transcript IDs."""
    ref_str = reference_names.astype(str)
    parts = ref_str.str.split("|")
    meta = pd.DataFrame(index=reference_names.index)
    meta["transcript_id"] = parts.str[0]
    meta["gene_id"] = parts.str[1].where(parts.str.len() > 1, parts.str[0])
    meta["transcript_name"] = parts.str[4].where(parts.str.len() > 4, parts.str[0])
    meta["gene_name"] = parts.str[5].where(parts.str.len() > 5, meta["gene_id"])
    return meta


def representative_transcripts_for_offset_metagene(
    annotated_read_df: pd.DataFrame,
) -> pd.DataFrame:
    """Keep one representative transcript per gene and read length.

    Representative transcripts are chosen from strict unique fragments using the
    transcript with the highest weighted support for each ``gene_id`` and
    ``read_length`` pair. This avoids multi-isoform genes contributing several
    transcript-level start profiles to the offset metagene.
    """
    if annotated_read_df.empty:
        return annotated_read_df
    if (
        "gene_id" not in annotated_read_df.columns
        or "transcript_id" not in annotated_read_df.columns
    ):
        return annotated_read_df

    weights = (
        annotated_read_df["count"].astype(int)
        if "count" in annotated_read_df.columns
        else pd.Series(1, index=annotated_read_df.index, dtype=int)
    )
    support = (
        annotated_read_df.assign(_w=weights)
        .groupby(["gene_id", "read_length", "transcript_id"], observed=True)["_w"]
        .sum()
        .reset_index()
        .sort_values(
            ["gene_id", "read_length", "_w", "transcript_id"],
            ascending=[True, True, False, True],
        )
        .drop_duplicates(["gene_id", "read_length"], keep="first")
    )
    chosen = set(
        zip(
            support["gene_id"].astype(str),
            support["read_length"].astype(str),
            support["transcript_id"].astype(str),
        )
    )
    keys = list(
        zip(
            annotated_read_df["gene_id"].astype(str),
            annotated_read_df["read_length"].astype(str),
            annotated_read_df["transcript_id"].astype(str),
        )
    )
    mask = pd.Series([key in chosen for key in keys], index=annotated_read_df.index)
    return annotated_read_df[mask]


def unique_fragments_single_gene_for_offset_metagene(
    annotated_read_df: pd.DataFrame,
) -> pd.DataFrame:
    """Keep fragments whose unique-mapper rows resolve to exactly one gene.

    STAR transcriptome BAMs can emit several transcript rows for one fragment.
    For offset calling, keep fragments only when all surviving unique-mapper
    rows agree on one ``gene_id``. Ambiguous cross-gene fragments are excluded.
    """
    if annotated_read_df.empty:
        return annotated_read_df
    if "read_name" not in annotated_read_df.columns or "gene_id" not in annotated_read_df.columns:
        return annotated_read_df

    gene_counts = annotated_read_df.groupby("read_name", observed=True)["gene_id"].nunique(
        dropna=True
    )
    keep_names = gene_counts[gene_counts == 1].index
    return annotated_read_df[annotated_read_df["read_name"].isin(keep_names)]


def frame_safe_unique_fragments(
    annotated_read_df: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse annotated reads to frame-safe fragments for periodicity.

    Transcriptome BAMs can produce multiple annotated rows per fragment. For
    frame-sensitive metrics, keep a fragment only when all surviving rows agree
    on one gene, one transcript, and one CDS frame. This is the conservative
    correctness-first path: transcript-ambiguous fragments are excluded even if
    the alternate transcript rows happen to imply the same frame.
    """
    if annotated_read_df.empty:
        return annotated_read_df
    required = {"read_name", "gene_id", "transcript_id", "read_frame"}
    if required - set(annotated_read_df.columns):
        return annotated_read_df

    group_cols = ["read_name"]
    if "read_length" in annotated_read_df.columns:
        group_cols.append("read_length")

    # Vectorised equivalent of "keep groups that agree on one gene, one
    # transcript and one frame, then collapse to a single row (count = max)".
    # A per-group Python loop here is O(n_fragments) with a DataFrame slice and
    # concat per group, which dominates runtime on deep libraries.
    grouped = annotated_read_df.groupby(group_cols, observed=True, sort=False)
    single_valued = (
        (grouped["gene_id"].transform("nunique") == 1)
        & (grouped["transcript_id"].transform("nunique") == 1)
        & (grouped["read_frame"].transform("nunique") == 1)
    )
    base = annotated_read_df[single_valued]
    if base.empty:
        return annotated_read_df.iloc[0:0].copy()

    if "count" in base.columns:
        # ``count`` arrives as a (non-ordered) categorical from BAM parsing;
        # pandas cannot aggregate ``max`` over an unordered Categorical, so cast
        # to int *before* the grouped reduction rather than after.
        max_count = (
            base.assign(count=base["count"].astype(int))
            .groupby(group_cols, observed=True, sort=False)["count"]
            .transform("max")
            .astype(int)
        )
        base = base.assign(count=max_count)

    return base.drop_duplicates(group_cols, keep="first").reset_index(drop=True)


def is_valid_offset(
    read_length: int,
    offset: Optional[int],
    offset_bounds: Tuple[int, int] = DEFAULT_OFFSET_BOUNDS,
    max_read_length_fraction: Optional[float] = DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION,
) -> bool:
    """Return True when an offset is numeric and plausible for a read length."""
    if offset is None:
        return False
    try:
        read_length_i = int(read_length)
        offset_i = int(offset)
    except (TypeError, ValueError):
        return False

    min_offset, max_offset = offset_bounds
    if max_read_length_fraction is not None:
        max_offset = min(
            max_offset,
            int(np.floor(read_length_i * float(max_read_length_fraction))),
        )

    return offset_i > 0 and offset_i < read_length_i and min_offset <= offset_i <= max_offset


def sanitise_offset(
    read_length: int,
    offset: Optional[int],
    default_offset: int = 15,
    offset_bounds: Tuple[int, int] = DEFAULT_OFFSET_BOUNDS,
    max_read_length_fraction: Optional[float] = DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION,
) -> int:
    """Return a usable offset, falling back when a caller produced noise."""
    if is_valid_offset(
        read_length,
        offset,
        offset_bounds,
        max_read_length_fraction,
    ):
        return int(cast(int, offset))
    if is_valid_offset(
        read_length,
        default_offset,
        offset_bounds,
        max_read_length_fraction,
    ):
        return int(default_offset)

    min_offset, max_offset = offset_bounds
    if max_read_length_fraction is not None:
        max_offset = min(
            max_offset,
            int(np.floor(int(read_length) * float(max_read_length_fraction))),
        )
    bounded_default = min(max(int(default_offset), min_offset), max_offset)
    return min(bounded_default, int(read_length) - 1)


def sanitise_offset_dict(
    offset_dict: Dict[int, int],
    default_offset: int = 15,
    offset_bounds: Tuple[int, int] = DEFAULT_OFFSET_BOUNDS,
    max_read_length_fraction: Optional[float] = DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION,
) -> Dict[int, int]:
    """Validate per-read-length offsets and replace implausible calls."""
    return {
        int(read_length): sanitise_offset(
            int(read_length),
            int(offset),
            default_offset=default_offset,
            offset_bounds=offset_bounds,
            max_read_length_fraction=max_read_length_fraction,
        )
        for read_length, offset in offset_dict.items()
    }


def read_df_to_cds_read_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the a_site_df to a cds_read_df by removing reads that do not
    map to the CDS

    Inputs:
        df: Dataframe containing the read information and annotation

    Outputs:
        cds_read_df: Dataframe containing the read information for reads
                    that map to the CDS
    """
    cds_read_df = df[(df["cds_start"] < df["a_site"]) & (df["a_site"] < df["cds_end"])]
    return cds_read_df


def a_site_calculation(
    read_df: pd.DataFrame,
    offset_file: str = "None",
    offset_type: str = "calculate",
    global_offset: int = 15,
    offset_bounds: Tuple[int, int] = DEFAULT_OFFSET_BOUNDS,
    max_read_length_fraction: Optional[float] = DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION,
) -> pd.DataFrame:
    """
    Adds a column to the read_df containing the A-site for the reads

    Inputs:
        read_df: Dataframe containing the read information
        offset_file: Path to a file containing offsets for each read length
        offset_type: Method to calculate offsets
                     Options: 'calculate', 'variable', 'file'

    Outputs:
        asite_df: Dataframe containing the read information with an added
                    column for the A-site
    """
    # Ensure reference_start is numeric before any arithmetic
    if "reference_start" in read_df.columns:
        read_df = read_df.copy()
        read_df["reference_start"] = pd.to_numeric(read_df["reference_start"], errors="coerce")

    if offset_type == "calculate":
        print("Calculating offsets")
        a_site_df = a_site_calculation_variable_offset(
            read_df,
            default_offset=global_offset,
            offset_bounds=offset_bounds,
            max_read_length_fraction=max_read_length_fraction,
        )
    elif offset_type == "read_length":
        # TSV with two columns: read_length<tab>offset.
        # Accepts an optional header row (e.g., "read_len\toffset").
        rl_table = pd.read_csv(
            offset_file,
            sep="\t",
            header=None,
            names=["read_length", "offset"],
            dtype={"read_length": str, "offset": str},
        )
        # Coerce to numeric; drop any non-numeric header-like rows
        rl_table["read_length_num"] = pd.to_numeric(rl_table["read_length"], errors="coerce")
        rl_table["offset_num"] = pd.to_numeric(rl_table["offset"], errors="coerce")
        rl_table = rl_table.dropna(subset=["read_length_num", "offset_num"]).astype(
            {
                "read_length_num": int,
                "offset_num": int,
            }
        )
        offset_dict = dict(zip(rl_table["read_length_num"], rl_table["offset_num"]))
        a_site_df = a_site_calculation_variable_offset(
            read_df,
            offset_dict,
            default_offset=global_offset,
            offset_bounds=offset_bounds,
            max_read_length_fraction=max_read_length_fraction,
        )
    elif offset_type == "global":
        df = read_df.copy()
        # Coerce again for safety under different pandas versions
        df["reference_start"] = pd.to_numeric(df["reference_start"], errors="coerce")
        df["offset"] = int(global_offset)
        df["a_site"] = df["reference_start"] + df["offset"]
        a_site_df = df
    elif offset_type == "read_specific":
        read_offsets = pd.read_csv(offset_file, sep="\t", names=["read_name", "offset"])

        merged_df = read_df.merge(read_offsets, on="read_name", how="left")
        # Robust numeric coercion to avoid object/str arithmetic issues on CI
        merged_df["offset"] = pd.to_numeric(merged_df["offset"], errors="coerce")
        merged_df["offset"] = merged_df["offset"].fillna(global_offset).astype(int)
        merged_df["reference_start"] = pd.to_numeric(merged_df["reference_start"], errors="coerce")
        merged_df["a_site"] = merged_df["reference_start"] + merged_df["offset"]

        a_site_df = merged_df[read_df.columns.tolist() + ["a_site", "offset"]]
    else:
        # Fallback to global behavior with robust numeric handling
        df = read_df.copy()
        df["reference_start"] = pd.to_numeric(df["reference_start"], errors="coerce")
        df["offset"] = int(global_offset)
        df["a_site"] = df["reference_start"] + df["offset"]
        a_site_df = df
    return a_site_df


def a_site_calculation_variable_offset(
    read_df: pd.DataFrame,
    offset_dict: Optional[dict] = None,
    default_offset: int = 15,
    offset_bounds: Tuple[int, int] = DEFAULT_OFFSET_BOUNDS,
    max_read_length_fraction: Optional[float] = DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION,
    validate_offsets: bool = False,
) -> pd.DataFrame:
    """
    Adds a column to the read_df containing the A-site for the reads

    Inputs:
        read_df: Dataframe containing the read information
        offset_dict: Dictionary containing offsets for each read length
                     Keys: read_length, Values: offset
                     If offset_dict is None, a default offset of 15 is
                     used for all read lengths.

    Outputs:
        asite_df: Dataframe containing the read information with an added
                    column for the A-site
    """
    default_offset = sanitise_offset(
        10**9,
        default_offset,
        default_offset,
        offset_bounds,
        max_read_length_fraction,
    )

    # If offset_dict is not provided, use default offset of
    # 15 for all read lengths
    if offset_dict is None:
        offset = default_offset
    else:
        if validate_offsets:
            offset_dict = sanitise_offset_dict(
                offset_dict,
                default_offset=default_offset,
                offset_bounds=offset_bounds,
                max_read_length_fraction=max_read_length_fraction,
            )
        # Map offsets to corresponding read lengths (cast to built-in int to avoid numpy int hash mismatch)
        read_len_int = read_df["read_length"].astype(int)
        read_df["offset"] = read_len_int.map(lambda l: int(offset_dict.get(int(l), default_offset)))
        read_df["offset"] = read_df["offset"].astype("int64")
        offset = read_df["offset"]

    read_df["reference_start"] = read_df["reference_start"].astype(int)

    # Calculate A-site based on offset for each read
    a_site_df = read_df.assign(a_site=read_df["reference_start"] + offset)
    return a_site_df


def read_length_distribution(read_df: pd.DataFrame) -> dict:
    """
    Calculate the read length distribution for the full dataset

    Inputs:
        read_df: Dataframe containing the read information

    Outputs:
        dict: Dictionary containing the read length distribution
    """
    weights = _get_weights(read_df)
    if weights is not None:
        counts = read_df.assign(_w=weights).groupby("read_length", observed=True)["_w"].sum()
        return {int(k): int(v) for k, v in counts.to_dict().items()}
    else:
        read_lengths, read_counts = np.unique(read_df["read_length"], return_counts=True)
        return dict(zip(read_lengths.tolist(), read_counts.tolist()))


def terminal_nucleotide_bias_distribution(
    read_df: pd.DataFrame,
    pattern_length: int = 2,
    keep_N: bool = False,
    target: str = "both",
) -> dict:
    """
    Calculate the proportion of the occurrence in the first or last n
    nucleotides of the reads to check for ligation bias

    Inputs:
        read_df: Dataframe containing read information
        # pattern_length: Length of nucleotide pattern
        keep_N: Keep nucleotide patterns with 'N', or discard if False
        target: Calculate ligation bias for 5', 3' or both

    Outputs:
        terminal_nucleotide_bias_dict: Dictionary containing the distribution
        of the first pattern of nucleotides in the reads
    """
    terminal_nucleotide_bias_dict: dict = (
        {target: {}} if target != "both" else {"five_prime": {}, "three_prime": {}}
    )

    weights = _get_weights(read_df)
    total_counts = weights.sum() if weights is not None else len(read_df)
    if weights is not None:
        tmp = read_df.assign(_w=weights)
        prime_counts = {
            "five_prime": tmp.groupby("first_dinucleotide", observed=True)["_w"].sum(),
            "three_prime": tmp.groupby("last_dinucleotide", observed=True)["_w"].sum(),
        }
    else:
        prime_counts = {
            "five_prime": read_df["first_dinucleotide"].value_counts(),
            "three_prime": read_df["last_dinucleotide"].value_counts(),
        }

    categories = {
        "five_prime": read_df["first_dinucleotide"].cat.categories.to_list(),
        "three_prime": read_df["last_dinucleotide"].cat.categories.to_list(),
    }

    pattern_list = read_df["first_dinucleotide"].cat.categories.to_list()
    pattern_list += read_df["last_dinucleotide"].cat.categories.to_list()
    pattern_list = sorted(list(set(categories["five_prime"]) | set(categories["three_prime"])))

    if keep_N:
        pattern_list = sorted(pattern_list, key=lambda x: ("N" in x, x))
    else:
        pattern_list = [pattern for pattern in pattern_list if "N" not in pattern]

    for pattern in pattern_list:
        for prime in terminal_nucleotide_bias_dict:
            if pattern in categories[prime]:
                # prime_counts[prime] is a pandas Series indexed by pattern;
                # .get returns 0 for patterns absent from this batch.
                val = prime_counts[prime].get(pattern, 0)
                terminal_nucleotide_bias_dict[prime][pattern] = (
                    float(val) / float(total_counts) if total_counts else 0.0
                )

    return terminal_nucleotide_bias_dict


def normalise_ligation_bias(
    terminal_nucleotide_bias_dict: dict,
    sequence_background: dict,
    pattern_length: int = 2,
) -> dict:
    """
    Calculate the difference between the observed and expected nucleotide
    pattern at the start and end of the sequences.

    Inputs:
        terminal_nucleotide_bias_dict: Dictionary containing observed
                            proportions for 5' and 3' ends of the sequences
        sequence_background: Dictionary containing expected proportions for 5'
                            and 3' directions of sequences
        # pattern_length: Length of nucleotide pattern

    Outputs:
        terminal_nucleotide_bias_dict_norm: Modified
                                terminal_nucleotide_bias_dict to show the
                                difference between observed and expected
                                distributions
    """
    # Work on a copy to avoid mutating the input
    terminal_nucleotide_bias_dict_norm = {
        prime: inner.copy() for prime, inner in terminal_nucleotide_bias_dict.items()
    }
    expected_distribution = {
        "five_prime": sequence_background["5_prime_bg"],
        "three_prime": sequence_background["3_prime_bg"],
    }

    for prime in terminal_nucleotide_bias_dict_norm:
        for pattern in terminal_nucleotide_bias_dict_norm[prime]:

            if pattern in expected_distribution[prime]:
                terminal_nucleotide_bias_dict_norm[prime][pattern] -= expected_distribution[prime][
                    pattern
                ]

    return terminal_nucleotide_bias_dict_norm


def slicer_vectorized(array: np.ndarray, start: int, end: int) -> np.ndarray:
    """
    String slicer for numpy arrays

    Note: https://stackoverflow.com/a/39045337

    Inputs:
        array: A numpy array of strings
        start: The start position of the slice
        end: The end position of the slice

    Outputs:
        sliced_array: An array consisting of only the selected characters
        from the input string array
    """
    sliced_array = array.view(str).reshape(len(array), -1)[:, start:end]
    return np.frombuffer(sliced_array.tobytes(), dtype=(str, end - start))


def nucleotide_composition(sequence_data_single: dict) -> dict:
    """
    Calculate the proportions of nucleotides for each read position

    Inputs:
        sequence_data_single: A dictionary containing the counts for single
        nucleotides on each read position

    Outputs:
        nucleotide_composition_dict: A dictionary containing the proportion
        for single nucleotides on each read position
    """
    read_length = len(sequence_data_single["A"])
    nucleotide_composition_dict: dict = {nt: [] for nt in ["A", "C", "G", "T"]}
    for position in range(read_length):
        position_count = 0
        for nt in sequence_data_single:
            position_count += sequence_data_single[nt][position]
        for nt in sequence_data_single:
            nucleotide_composition_dict[nt].append(
                sequence_data_single[nt][position] / position_count
            )

    return nucleotide_composition_dict


def read_frame_cull(read_frame_dict: dict, config: dict) -> dict:
    """
    Culls the read_frame_dict according to config so only read lengths of
    interest are kept

    Inputs:
    read_frame_dict:
    config:

    Outputs:
    culled_read_frame_dict
    """
    culled_read_frame_dict = read_frame_dict.copy()
    cull_list = list(culled_read_frame_dict.keys())
    for k in cull_list:
        if (
            k > config["plots"]["read_frame_distribution"]["upper_limit"]
            or k < config["plots"]["read_frame_distribution"]["lower_limit"]
        ):
            del culled_read_frame_dict[k]

    return culled_read_frame_dict


def read_frame_score_trips_viz(read_frame_dict: dict) -> dict:
    """
    Generates scores for each read_length separately as well as a global score
    Can be used after read_frame_cull to calculate the global score of the
    region of interest. The calculation for this score is: 1 - sum(2nd highest
    peak count)/sum(highest peak count). A score close to 1 has good
    periodicity, while a score closer to 0 has a random spread

    Inputs:
    read_frame_dict: dictionary containing the distribution of the reading
                    frames over the different read lengths

    Outputs:
    scored_read_frame_dict: dictionary containing read frame distribution
                            scores for each read length and a global score
    """
    scored_read_frame_dict: Dict[str, float] = {}
    highest_peak_sum, second_peak_sum = 0, 0
    for k, inner_dict in read_frame_dict.items():
        top_two_values = sorted(inner_dict.values(), reverse=True)[:2]
        if top_two_values[0] == 0:
            scored_read_frame_dict[k] = 0
            continue
        elif top_two_values[1] == 0:
            scored_read_frame_dict[k] = 1
        else:
            highest_peak_sum += top_two_values[0]
            second_peak_sum += top_two_values[1]
            scored_read_frame_dict[k] = 1 - top_two_values[1] / top_two_values[0]

    if highest_peak_sum == 0:
        scored_read_frame_dict["global"] = 0.0
    else:
        scored_read_frame_dict["global"] = 1 - second_peak_sum / highest_peak_sum
    return scored_read_frame_dict


def read_frame_distribution(a_site_df: pd.DataFrame) -> dict:
    """
    Calculate the distribution of the reading frame over the dataset

    Inputs:
    a_site_df: Dataframe containing the read information with an
                added column for the a-site location

    Outputs:
    read_frame_dict: Nested dictionary containing counts for every
                    reading frame at the different read lengths
    """
    read_frame_dict = {}

    # Weighted counts (use heuristic to avoid double-counting expanded rows)
    weights = _get_weights(a_site_df)

    # Compute read_frame per row
    a_site_df = a_site_df.assign(read_frame=(a_site_df["a_site"] % 3).astype(int))

    # Group by read_length and frame; sum weights or sizes
    if weights is not None:
        grouped = (
            a_site_df.assign(_w=weights)
            .groupby(["read_length", "read_frame"], observed=True)["_w"]
            .sum()
        )
    else:
        grouped = a_site_df.groupby(["read_length", "read_frame"], observed=True).size()

    for (read_length, frame), value in grouped.items():
        read_length = int(read_length)
        if read_length not in read_frame_dict:
            read_frame_dict[read_length] = {0: 0, 1: 0, 2: 0}
        read_frame_dict[read_length][int(frame)] = int(value)

    return read_frame_dict


def read_frame_distribution_annotated(
    annotated_read_df: pd.DataFrame,
    exclusion_length: int = 0,
    read_length_range: tuple = (20, 40),
    unique_only: bool = True,
) -> dict:
    """
    Calculate the distribution of the reading frame over the dataset

    Inputs:
        a_site_df: Dataframe containing the read information with an added
        column for the a-site location

    Outputs:
        read_frame_dict: Nested dictionary containing counts for every reading
        frame at the different read lengths
    """
    read_lengths: List[int] = [i for i in range(read_length_range[0], read_length_range[1])]

    df_slice = filter_unique_mappers(annotated_read_df, enabled=unique_only)
    df_slice = df_slice[df_slice["cds_start"] != 0]
    df_slice = df_slice[
        (df_slice["a_site"] > df_slice["cds_start"] + exclusion_length)
        & (df_slice["a_site"] < df_slice["cds_end"] - exclusion_length)
    ]
    base = df_slice.assign(read_frame=(df_slice.a_site - df_slice.cds_start).mod(3))
    base = frame_safe_unique_fragments(base)
    weights = _get_weights(base)
    if weights is not None:
        frame_df = (
            base.assign(_w=weights)
            .groupby(["read_length", "read_frame"], observed=True)["_w"]
            .sum()
        )
    else:
        frame_df = base.groupby(["read_length", "read_frame"], observed=True).size()
    read_frame_dict: Dict[int, Dict[int, int]] = {}
    for index, value in frame_df.items():
        read_length: int
        read_frame: int
        read_length, read_frame = index
        if read_length in read_lengths:
            if read_length not in read_frame_dict:
                read_frame_dict[read_length] = {0: 0, 1: 0, 2: 0}
            read_frame_dict[read_length][read_frame] = value
    return read_frame_dict


def annotate_reads(a_site_df: pd.DataFrame, annotation_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merges the annotation dataframe with the read dataframe

    Inputs:
        a_site_df: Dataframe containing the read information with an added
        column for the a-site location
        annotation_df: Dataframe containing the CDS start/stop
        and transcript id from a gff file.

    Outputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along
        with the columns from the gff file
    """
    metadata = _reference_metadata(a_site_df.reference_name)
    annotated_read_df = a_site_df.assign(
        transcript_id=metadata["transcript_id"],
        gene_id=metadata["gene_id"],
        transcript_name=metadata["transcript_name"],
        gene_name=metadata["gene_name"],
    ).merge(annotation_df, on="transcript_id")
    annotated_read_df["transcript_id"] = annotated_read_df["transcript_id"].astype("category")
    annotated_read_df["gene_id"] = annotated_read_df["gene_id"].astype("category")
    annotated_read_df["transcript_name"] = annotated_read_df["transcript_name"].astype("category")
    annotated_read_df["gene_name"] = annotated_read_df["gene_name"].astype("category")
    return annotated_read_df.drop(["reference_name"], axis=1)


def chunked_annotate_reads(
    a_site_df: pd.DataFrame, annotation_df: pd.DataFrame, chunk_size: int = 10000000
) -> pd.DataFrame:
    """
    Merges the annotation dataframe with the read dataframe in smaller chunks.

    Inputs:
        a_site_df: DataFrame containing the read information with an added
        column for the a-site location.
        annotation_df: DataFrame containing the CDS start/stop
        and transcript id from a gff file.
        chunk_size: Size of each processing chunk.

    Outputs:
        annotated_read_df: DataFrame containing the read information
        with an added column for the a-site location along
        with the columns from the gff file.
    """
    # Initialize an empty list to store processed chunks
    processed_chunks = []

    # Split a_site_df into chunks
    num_chunks = len(a_site_df) // chunk_size + 1
    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(a_site_df))
        chunk = a_site_df.iloc[start_idx:end_idx]

        # Process the chunk
        metadata = _reference_metadata(chunk.reference_name)
        chunk = chunk.assign(
            transcript_id=metadata["transcript_id"],
            gene_id=metadata["gene_id"],
            transcript_name=metadata["transcript_name"],
            gene_name=metadata["gene_name"],
        )

        chunk = chunk.drop(["reference_name"], axis=1)
        chunk = chunk.merge(annotation_df, on="transcript_id")
        chunk["transcript_id"] = chunk["transcript_id"].astype("category")
        chunk["gene_id"] = chunk["gene_id"].astype("category")
        chunk["transcript_name"] = chunk["transcript_name"].astype("category")
        chunk["gene_name"] = chunk["gene_name"].astype("category")

        # Append the processed chunk to the list
        processed_chunks.append(chunk)

    # Concatenate the processed chunks
    annotated_read_df = pd.concat(processed_chunks)

    return annotated_read_df


def assign_mRNA_category(annotated_read_df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds the mRNA category column to the annotated_read_df, labelling the read
    according to the position of the A-site
    Assign an mRNA category based on the A-site of the read
    and the CDS start/stop, used through df.apply()

    Inputs:
        annotated_read_df: Dataframe with read data, added a-site positions
        and joined with annotation_df.

    Outputs:
        mRNA category: string with the category for the read
        ["five_leader", "start_codon", "CDS", "stop_codon", "three_trailer"]
    """
    # Calculate mRNA category based on conditions
    conditions = [
        annotated_read_df["a_site"] < annotated_read_df["cds_start"],
        annotated_read_df["a_site"] == annotated_read_df["cds_start"],
        (annotated_read_df["cds_start"] < annotated_read_df["a_site"])
        & (annotated_read_df["a_site"] < annotated_read_df["cds_end"]),
        annotated_read_df["a_site"] == annotated_read_df["cds_end"],
        annotated_read_df["a_site"] > annotated_read_df["cds_end"],
    ]
    choices = ["five_leader", "start_codon", "CDS", "stop_codon", "three_trailer"]
    annotated_read_df["mRNA_category"] = np.select(conditions, choices, "unknown")
    annotated_read_df["mRNA_category"] = annotated_read_df["mRNA_category"].astype("category")
    return annotated_read_df


def mRNA_distribution(annotated_read_df: pd.DataFrame) -> Dict[int, Dict[str, int]]:
    """
    Calculate the distribution of the mRNA categories over the read length

    Inputs:
        annotated_read_df: Dataframe containing the read information
                           with an added column for the a-site location along
                           with the columns from the gff file
    Outputs:
        mRNA_distribution_dict: Nested dictionary containing counts for every
                                mRNA category at the different read lengths
    """
    # Creating MultiIndex for reindexing
    categories = [
        "five_leader",
        "start_codon",
        "CDS",
        "stop_codon",
        "three_trailer",
    ]
    classes = annotated_read_df["read_length"].unique()
    idx = pd.MultiIndex.from_product([classes, categories], names=["class", "category"])
    # Group annotated_read_df
    weights = _get_weights(annotated_read_df)
    if weights is not None:
        grp = (
            annotated_read_df.assign(_w=weights)
            .groupby(["read_length", "mRNA_category"], observed=True)["_w"]
            .sum()
            .reindex(idx, fill_value=0)
            .sort_index()
            .to_frame(name=0)
            .reset_index()
        )
    else:
        grp = (
            annotated_read_df.groupby(["read_length", "mRNA_category"], observed=True)
            .size()
            .reindex(idx, fill_value=0)
            .sort_index()
            .to_frame()  # value column named 0
            .reset_index()
        )
    annotated_read_df = grp

    # Creating mRNA_distribution_dict from annotated_read_df
    mRNA_distribution_dict: dict = {"global": {}}
    for read_length, mRNA_category, value in annotated_read_df.itertuples(index=False, name=None):
        if read_length not in mRNA_distribution_dict:
            mRNA_distribution_dict[read_length] = {}
        mRNA_distribution_dict[read_length][mRNA_category] = value
        if mRNA_category in mRNA_distribution_dict["global"]:
            mRNA_distribution_dict["global"][mRNA_category] += value
        else:
            mRNA_distribution_dict["global"][mRNA_category] = value
    return mRNA_distribution_dict


def sum_mRNA_distribution(mRNA_distribution_dict: dict, config: dict) -> dict:
    """
    Calculate the sum of mRNA categories

    Inputs:
        annotated_read_dict: Dataframe containing the read information
        with an added column for the a-site location along
        with the columns from the gff file

    Outputs:
        read_frame_dict: Nested dictionary containing counts for every reading
        frame at the different read lengths
    """
    sum_mRNA_dict: dict = {}
    for inner_dict in mRNA_distribution_dict.values():
        for k, v in inner_dict.items():
            if k in sum_mRNA_dict:
                sum_mRNA_dict[k] += v
            else:
                sum_mRNA_dict[k] = v
    if not config["plots"]["mRNA_distribution"]["absolute_counts"]:
        sum_mRNA_dict = {k: (v / sum(sum_mRNA_dict.values())) for k, v in sum_mRNA_dict.items()}

    return sum_mRNA_dict


def metagene_distance(
    annotated_read_df: pd.DataFrame,
    target: str = "start",
    position: str = "a_site",
) -> pd.Series:
    """
    Calculate distance from a read position to start or stop codon

    Inputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along with data from
        the annotation file
        target: Target from which the distance is calculated
        position: Column to use as the read position. Use "a_site" (default)
            for the standard metagene profile, or "reference_start" to build
            a 5'-end metagene (required for accurate offset detection).

    Outputs:
    pd.Series
    """
    if target == "start":
        return annotated_read_df[position] - annotated_read_df["cds_start"]
    elif target == "stop":
        return annotated_read_df[position] - annotated_read_df["cds_end"]
    else:
        raise ValueError("Target must be start or stop")


def metagene_profile(
    annotated_read_df: pd.DataFrame,
    target: str = "both",
    distance_range: list = [-50, 50],
    position: str = "a_site",
    extend: bool = False,
) -> dict:
    """
    Groups the reads by read_length and distance to a target and counts them

    Inputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along with data from
        the annotation file
        target: Target from which the distance is calculated
        distance_range: The range of the plot
        position: The position to plot
        extend: Extend the range to include all read lengths

    Outputs:
        metagene_profile_dict: dictionary containing the read_length of
        the read and distance to the target as keys and the counts as values
    """
    target_loop = [target] if target != "both" else ["start", "stop"]
    metagene_profile_dict: Dict[str, Dict[str, dict]] = {"start": {}, "stop": {}}
    for current_target in target_loop:
        annotated_read_df = annotated_read_df.assign(
            metagene_info=metagene_distance(annotated_read_df, current_target, position)
        )
        filtered = annotated_read_df[
            (annotated_read_df["metagene_info"] > distance_range[0] - 1)
            & (annotated_read_df["metagene_info"] < distance_range[1] + 1)
        ]
        wts = _get_weights(filtered)
        if wts is not None:
            pre_series = (
                filtered.assign(_w=wts)
                .groupby(["read_length", "metagene_info"], observed=True)["_w"]
                .sum()
            )
        else:
            pre_series = filtered.groupby(["read_length", "metagene_info"], observed=True).size()
        pre_metaprofile_dict = pre_series.to_dict()
        if pre_metaprofile_dict == {}:
            if extend:  # If no reads in range
                w_all = _get_weights(annotated_read_df)
                if w_all is not None:
                    pre_metaprofile_dict = (
                        annotated_read_df.assign(_w=w_all)
                        .groupby(["read_length", "metagene_info"], observed=True)["_w"]
                        .sum()
                        .to_dict()
                    )
                else:
                    pre_metaprofile_dict = (
                        annotated_read_df.groupby(["read_length", "metagene_info"], observed=True)
                        .size()
                        .to_dict()
                    )
            else:
                return {
                    "start": {0: {1: 0}},
                    "stop": {0: {1: 0}},
                }

        # Fill empty read lengths with 0
        min_length = int(min([x[0] for x in list(pre_metaprofile_dict.keys())]))
        max_length = int(max([x[0] for x in list(pre_metaprofile_dict.keys())]))

        for y in range(min_length, max_length + 1):
            if y not in [x[0] for x in list(pre_metaprofile_dict.keys())]:
                pre_metaprofile_dict[(y, 0)] = 0

        neg_distance = int(min([x[1] for x in list(pre_metaprofile_dict.keys())]))
        pos_distance = int(max([x[1] for x in list(pre_metaprofile_dict.keys())]))
        position_range = range(neg_distance, pos_distance + 1)

        for key, value in pre_metaprofile_dict.items():
            if key[0] not in metagene_profile_dict[current_target]:
                metagene_profile_dict[current_target][key[0]] = {}
            metagene_profile_dict[current_target][key[0]][int(key[1])] = value

        # Fill empty distances with 0
        for position_dict in metagene_profile_dict[current_target].values():
            for pos in position_range:
                position_dict.setdefault(int(pos), 0)

    return metagene_profile_dict


def reading_frame_triangle(
    annotated_read_df: pd.DataFrame,
) -> Dict[str, List[int]]:
    """
    Get the per-transcript reading-frame counts used by the triangle plot.

    Inputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along with data from
        the annotation file

    Outputs:
        triangle_dict: Dictionary mapping each transcript_id to its raw
        per-frame A-site counts ``[frame0, frame1, frame2]``. The conversion
        to cartesian (ternary) coordinates happens in ``plots.py``.
    """
    if annotated_read_df.empty or "transcript_id" not in annotated_read_df.columns:
        return {}

    # Vectorised per-transcript frame counts (a_site % 3), weighted by 'count'
    # when present. Replaces a per-transcript Python groupby loop that scaled
    # linearly in the number of transcripts (~hundreds of thousands).
    frame = (annotated_read_df["a_site"].to_numpy() % 3).astype(int)
    counts = pd.DataFrame(
        {
            "transcript_id": annotated_read_df["transcript_id"].to_numpy(),
            "frame": frame,
        }
    )
    if "count" in annotated_read_df.columns:
        counts["w"] = annotated_read_df["count"].astype(int).to_numpy()
        grouped = counts.groupby(["transcript_id", "frame"], observed=True)["w"].sum()
    else:
        grouped = counts.groupby(["transcript_id", "frame"], observed=True).size()

    matrix = grouped.unstack("frame").reindex(columns=[0, 1, 2], fill_value=0).fillna(0).astype(int)
    return {
        tid: [int(row[0]), int(row[1]), int(row[2])]
        for tid, row in zip(matrix.index, matrix.to_numpy())
    }


def sequence_slice(read_df: pd.DataFrame, nt_start: int = 0, nt_count: int = 15) -> Dict[str, str]:
    sequence_slice_dict = {
        k: v[nt_start : nt_start + nt_count] for k, v in read_df["sequence"].to_dict().items()
    }
    return sequence_slice_dict


def convert_html_to_pdf(source_html: str, output_filename: str) -> int:
    if not HAS_PDF:
        raise ImportError(
            "PDF export requires xhtml2pdf. Install with: pip install RiboMetric[pdf]"
        )
    result_file = open(output_filename, "w+b")

    pisa_status = pisa.CreatePDF(source_html, dest=result_file)
    result_file.close()
    # xhtml2pdf does not provide type hints; coerce to int explicitly
    try:
        return int(getattr(pisa_status, "err"))
    except Exception:
        return 1


def ribowaltz_psite_prediction(
    read_counts: Dict[int, Dict[int, int]],
    flanking_length: int = 9,
    offset_range: Tuple[int, int] = (10, 18),
    min_prominence: Optional[float] = None,
) -> Dict[int, Optional[int]]:
    """
    Predict P-site offsets for each read length from a 5'-end metagene.

    The metagene must be built using reference_start (5' end of the read) as
    the position, so the P-site signal appears as a peak at a NEGATIVE position
    upstream of the start codon.  The returned offsets are positive integers
    representing the distance from the 5' end of the read to the P-site.

    Args:
        read_counts: {read_length: {metagene_position: count}}
            Positions should be reference_start - cds_start (negative = upstream).
        flanking_length: Exclude positions within this distance of the start codon
            to avoid initiation peak noise. Default: 9.
        offset_range: Allowed offsets (inclusive) in nucleotides. Default: (10, 18).
        min_prominence: If provided, require (peak / median) >= min_prominence
            within the search window; otherwise treat as missing and fall back
            to consensus. Example useful range: 1.5–2.0.

    Returns:
        {read_length: psite_offset (positive int or None)}
    """
    psite_offsets: Dict[int, Optional[int]] = {}
    read_lengths = sorted(read_counts.keys())

    lo, hi = offset_range
    # Negative positions to search (e.g., -18..-10); inclusive bounds
    window_positions = set(range(-hi, -lo + 1))

    candidate_offsets: Dict[int, Optional[int]] = {}
    for read_length in read_lengths:
        counts = read_counts[read_length]
        # Strictly upstream positions past flanking, constrained to window
        upstream_pos = [
            p for p in counts.keys() if (p < -flanking_length and p in window_positions)
        ]

        if not upstream_pos:
            # Fallback: allow any upstream position past flanking if window empty
            upstream_pos = [p for p in counts.keys() if p < -flanking_length]

        if not upstream_pos:
            candidate_offsets[read_length] = None
            continue

        upstream_pos.sort()
        upstream_vals = np.array([counts.get(p, 0) for p in upstream_pos], dtype=float)
        peak_idx = int(np.argmax(upstream_vals))
        peak_pos = upstream_pos[peak_idx]
        peak_cnt = upstream_vals[peak_idx]

        # Optional prominence check (peak vs. median of window)
        if min_prominence is not None:
            med = float(np.median(upstream_vals)) if upstream_vals.size else 0.0
            med = med if med > 0 else 1.0
            if (peak_cnt / med) < float(min_prominence):
                candidate_offsets[read_length] = None
                continue

        # P-site offset = |peak position|
        candidate_offsets[read_length] = abs(int(peak_pos))

    # Consensus fallback across read lengths
    votes: Dict[int, int] = {}
    for off in candidate_offsets.values():
        if off is not None:
            votes[off] = votes.get(off, 0) + 1
    consensus: Optional[int] = (
        max(votes, key=lambda off: cast(int, votes.get(off))) if votes else None
    )

    for rl in read_lengths:
        psite_offsets[rl] = (
            candidate_offsets[rl] if candidate_offsets[rl] is not None else consensus
        )

    return psite_offsets


def trips_psite_prediction(
    read_counts: Dict[int, Dict[int, int]],
    offset_range: Tuple[int, int] = (10, 18),
) -> Dict[int, Optional[int]]:
    """
    Predict P-site offsets for each read length using a Trips-Viz-style peak.

    The input is a 5' read-end metagene around the start codon. To keep all
    automatic modes on the same target, this returns P-site offsets; callers
    convert to A-site or P-site according to ``offset_target``.

    Args:
        read_counts (dict): Dictionary of read counts per position for each read length.
                            Format: {read_length: {position: count}}

    Returns:
        psite_offsets (dict): Dictionary of P-site offsets for each read length.
                              Format: {read_length: offset}
    """

    psite_offsets: Dict[int, Optional[int]] = {}
    lo, hi = offset_range
    window_positions = set(range(-hi, -lo + 1))
    for read_length, counts in read_counts.items():
        counts = {
            int(pos): count
            for pos, count in counts.items()
            if int(pos) < 0 and int(pos) in window_positions
        }

        if counts:
            max_pos = max(counts, key=lambda k: counts[k])
            offset = abs(int(max_pos))
        else:
            offset = None
        psite_offsets[read_length] = offset

    return psite_offsets


def trips_asite_prediction(
    read_counts: Dict[int, Dict[int, int]],
    offset_range: Tuple[int, int] = (10, 18),
) -> Dict[int, Optional[int]]:
    """Compatibility wrapper for the old Trips helper name."""
    return trips_psite_prediction(read_counts, offset_range=offset_range)


def asite_calculation_per_readlength(
    annotated_read_df: pd.DataFrame,
    method: str = "ribowaltz",
    offset_range: Tuple[int, int] = (10, 18),
    default_offset: int = 15,
    offset_bounds: Tuple[int, int] = DEFAULT_OFFSET_BOUNDS,
    max_read_length_fraction: Optional[float] = DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION,
    min_prominence: Optional[float] = None,
    unique_only: bool = True,
    offset_target: str = "a_site",
) -> Dict[int, int]:
    """
    Calculate offset values per read length for the A-site
    using an improved change point detection method.

    Input:
        annotated_read_df: DataFrame with read counts and CDS info
        offset_range: Range of allowed offsets
        default_offset: Default offset to use if no significant
                        change point is found

    Output:
        offset_dict: Dictionary containing the offset values for
                    each read length
    """
    offset_dict: Dict[int, int] = {}
    print(
        "Running offset calculation per read length "
        f"with {method} method targeting {offset_target}"
    )
    target_shift = offset_shift_for_target(offset_target)
    # Build the offset-calling metagene from fragment-level unique reads that
    # resolve to one gene, then choose one representative transcript per gene.
    unique_df = filter_unique_mappers(annotated_read_df, enabled=unique_only)
    unique_df = unique_fragments_single_gene_for_offset_metagene(unique_df)
    unique_df = representative_transcripts_for_offset_metagene(unique_df)
    for read_length in unique_df["read_length"].unique():
        # Build metagene using the 5' end of reads (reference_start) rather than
        # the preliminary a_site.  For typical A-site offsets of 12–18 nt the
        # start-codon pile-up sits at positions −12 to −17 in the 5'-end metagene,
        # which falls squarely in the detection range of both changepoint and
        # ribowaltz.  Using the preliminary a_site metagene placed that pile-up
        # near position 0, inside the exclusion zone of both methods, so they
        # were detecting noise rather than signal.
        read_length_metagene = metagene_profile(
            unique_df[unique_df["read_length"] == read_length],
            target="start",
            distance_range=[-50, 20],
            position="reference_start",
            extend=False,
        )
        if read_length not in read_length_metagene["start"]:
            offset_dict[read_length] = sanitise_offset(
                int(read_length),
                default_offset,
                default_offset=default_offset,
                offset_bounds=offset_bounds,
                max_read_length_fraction=max_read_length_fraction,
            )
            continue

        if method == "changepoint":
            # Find the position with the highest read count in the upstream
            # region of the 5'-end metagene. The pile-up from ribosomes
            # stalled at the start codon appears at position -(P-site offset)
            # (negative = upstream). The caller-selected target controls
            # whether that P-site offset is used directly or shifted to A-site.
            counts = read_length_metagene["start"][read_length]
            # Search positions corresponding to P-site offsets in offset_range
            candidate_positions = {
                pos: counts.get(pos, 0) for pos in range(-offset_range[1], -offset_range[0] + 1)
            }
            peak_count = max(candidate_positions.values(), default=0)
            if peak_count == 0:
                offset = default_offset
            else:
                peak_pos = max(candidate_positions, key=lambda k: candidate_positions[k])
                psite_offset = abs(int(peak_pos))
                offset = int(psite_offset + target_shift)

        elif method == "ribowaltz":
            # ribowaltz_psite_prediction returns a positive P-site offset
            # (distance from 5' end to P-site). Shift according to target.
            rw_offset: Optional[int] = ribowaltz_psite_prediction(
                {read_length: read_length_metagene["start"][read_length]},
                flanking_length=9,
                offset_range=offset_range,
                min_prominence=min_prominence,
            )[read_length]
            if rw_offset is None:
                offset = default_offset
            else:
                offset = int(rw_offset + target_shift)

        elif method == "tripsviz":
            t = trips_psite_prediction(
                {read_length: read_length_metagene["start"][read_length]},
                offset_range=offset_range,
            )[read_length]

            if t is None:
                offset = default_offset
            else:
                offset = int(t + target_shift)

        else:
            raise ValueError(f"Invalid method: {method}")

        raw_offset = abs(int(offset))
        offset_dict[read_length] = sanitise_offset(
            int(read_length),
            raw_offset,
            default_offset=default_offset,
            offset_bounds=offset_bounds,
            max_read_length_fraction=max_read_length_fraction,
        )
        if offset_dict[read_length] != raw_offset:
            print(
                "Warning: computed offset "
                f"{raw_offset} for read length {int(read_length)} is outside "
                f"bounds {offset_bounds[0]}-{offset_bounds[1]} or not below "
                f"read length; using {offset_dict[read_length]}"
            )

    return offset_dict


def gene_body_coverage_ramp(
    cds_read_df: pd.DataFrame,
    n_bins: int = 100,
) -> Dict[str, object]:
    """Metagene of A-site density across *relative* CDS position (0-100%).

    Unlike the start/stop metagenes (which are anchored at the termini in
    nucleotides), this profiles the whole CDS body in relative coordinates, so
    transcripts of different lengths are comparable. It exposes:

    * the 5' translation "ramp" — elevated density just after the start codon,
      reflecting slower early elongation (and run-off artefacts);
    * the 3' drop-off — declining density toward the stop codon.

    Returns:
        profile                 list[float] of length n_bins, mean-normalised so
                                a flat profile is ~1.0 everywhere
        five_prime_ramp_ratio   mean density in the first 10% / middle (40-60%)
        three_prime_drop_ratio  mean density in the last 10% / middle (40-60%)
    """
    empty: Dict[str, object] = {
        "profile": [0.0] * n_bins,
        "five_prime_ramp_ratio": None,
        "three_prime_drop_ratio": None,
    }
    if cds_read_df.empty:
        return empty

    df = cds_read_df
    cds_len = (df["cds_end"] - df["cds_start"]).astype(float)
    valid = cds_len > 0
    if not valid.any():
        return empty
    df = df[valid]
    cds_len = cds_len[valid]
    rel = (df["a_site"].astype(float) - df["cds_start"].astype(float)) / cds_len
    in_range = (rel >= 0) & (rel < 1)
    if not in_range.any():
        return empty
    rel = rel[in_range]
    weights = (
        df.loc[in_range, "count"].astype(float)
        if "count" in df.columns
        else pd.Series(1.0, index=rel.index)
    )

    bin_idx = np.minimum((rel * n_bins).astype(int), n_bins - 1)
    profile = np.zeros(n_bins, dtype=float)
    np.add.at(profile, bin_idx.to_numpy(), weights.to_numpy())

    mean_density = profile.mean()
    if mean_density <= 0:
        return empty
    norm_profile = profile / mean_density

    lo = max(1, n_bins // 10)
    mid_lo, mid_hi = int(n_bins * 0.4), int(n_bins * 0.6)
    middle = norm_profile[mid_lo:mid_hi].mean()
    if middle <= 0:
        ramp = drop = None
    else:
        ramp = round(float(norm_profile[:lo].mean() / middle), 4)
        drop = round(float(norm_profile[-lo:].mean() / middle), 4)

    return {
        "profile": [round(float(x), 4) for x in norm_profile],
        "five_prime_ramp_ratio": ramp,
        "three_prime_drop_ratio": drop,
    }


def library_complexity_curve(
    annotated_read_df: pd.DataFrame,
    n_points: int = 10,
) -> Dict[str, object]:
    """Analytic rarefaction (saturation) curve of distinct A-site positions.

    Answers "was this library sequenced deeply enough?" without random
    subsampling: given the per-position read counts c_i, the expected number of
    distinct (transcript, A-site) positions recovered when sampling a fraction f
    of reads is sum_i [1 - (1 - f)^c_i]. A curve that has flattened by f=1 means
    extra sequencing buys little new signal (low complexity / saturated); a still-
    rising curve means the library is under-sequenced.

    Returns:
        fractions                       list of sampling fractions
        distinct_positions              expected distinct positions at each f
        total_distinct_positions        distinct positions at full depth
        marginal_discovery_rate         fraction of reads at the margin (top 5%
                                        of depth) landing on a new position;
                                        high => still discovering (undersaturated)
    """
    empty: Dict[str, object] = {
        "fractions": [],
        "distinct_positions": [],
        "total_distinct_positions": 0,
        "marginal_discovery_rate": None,
    }
    if annotated_read_df.empty or "a_site" not in annotated_read_df.columns:
        return empty
    tx_col = (
        "transcript_id"
        if "transcript_id" in annotated_read_df.columns
        else "reference_name" if "reference_name" in annotated_read_df.columns else None
    )
    if tx_col is None:
        return empty

    weights = _get_weights(annotated_read_df)
    grp = (
        annotated_read_df.assign(_w=(weights if weights is not None else 1))
        .groupby([tx_col, "a_site"], observed=True)["_w"]
        .sum()
    )
    counts = grp.to_numpy(dtype=float)
    counts = counts[counts > 0]
    if counts.size == 0:
        return empty

    total_reads = float(counts.sum())

    def expected_distinct(frac: float) -> float:
        # sum_i [1 - (1 - f)^c_i]; clamp base to [0,1] for safety.
        base = max(0.0, min(1.0, 1.0 - frac))
        return float(np.sum(1.0 - np.power(base, counts)))

    fractions = [round((i + 1) / n_points, 4) for i in range(n_points)]
    distinct = [round(expected_distinct(f), 2) for f in fractions]
    total_distinct = distinct[-1]

    d_full = expected_distinct(1.0)
    d_margin = expected_distinct(0.95)
    margin_reads = 0.05 * total_reads
    marginal_rate = (
        round(float((d_full - d_margin) / margin_reads), 4) if margin_reads > 0 else None
    )

    return {
        "fractions": fractions,
        "distinct_positions": distinct,
        "total_distinct_positions": int(round(total_distinct)),
        "marginal_discovery_rate": marginal_rate,
    }


def floss_library_heterogeneity(
    annotated_read_df: pd.DataFrame,
    min_reads_per_transcript: int = 20,
    floss_cutoff: float = 0.3,
) -> Dict[str, object]:
    """Library-level FLOSS-style read-length heterogeneity.

    FLOSS (Fragment Length Organization Similarity Score; Ingolia et al. 2014)
    measures how far a feature's footprint-length histogram departs from a
    reference distribution:

        FLOSS = 0.5 * Σ_l | F_transcript(l) - F_reference(l) |

    Here the reference is the library's own aggregate CDS read-length
    distribution. Computed per transcript and summarised across the library, it
    becomes a *sample-level* QC of footprint-length homogeneity: a high fraction
    of transcripts with aberrant length profiles indicates a heterogeneous or
    contaminated library (non-ribosomal fragments, mixed protocols), even when
    the aggregate length distribution looks fine. This is a QC summary, not an
    ORF/translation classifier.

    Returns:
        floss_median                       median per-transcript FLOSS
        floss_aberrant_transcript_fraction fraction of transcripts with
                                           FLOSS > floss_cutoff
        n_transcripts_scored               transcripts meeting min_reads
    """
    empty: Dict[str, object] = {
        "floss_median": None,
        "floss_aberrant_transcript_fraction": None,
        "n_transcripts_scored": 0,
    }
    if "transcript_id" not in annotated_read_df.columns:
        return empty
    if annotated_read_df.empty:
        return empty

    weights = _get_weights(annotated_read_df)
    df = annotated_read_df.assign(_w=(weights if weights is not None else 1))
    # Per (transcript, read_length) weighted counts -> length histograms.
    grp = (
        df.groupby(["transcript_id", "read_length"], observed=True)["_w"]
        .sum()
        .unstack(fill_value=0)
    )
    if grp.empty:
        return empty

    read_lengths = list(grp.columns)
    totals = grp.sum(axis=1)
    # Reference = aggregate distribution over all reads, normalised.
    ref_counts = grp.sum(axis=0).to_numpy(dtype=float)
    ref_total = ref_counts.sum()
    if ref_total <= 0:
        return empty
    ref_dist = ref_counts / ref_total

    floss_scores: List[float] = []
    for transcript in grp.index:
        n = float(totals.loc[transcript])
        if n < min_reads_per_transcript:
            continue
        t_dist = grp.loc[transcript].to_numpy(dtype=float) / n
        floss = 0.5 * float(np.abs(t_dist - ref_dist).sum())
        floss_scores.append(floss)

    if not floss_scores:
        return empty

    scores = np.array(floss_scores, dtype=float)
    aberrant = float((scores > floss_cutoff).mean())
    return {
        "floss_median": round(float(np.median(scores)), 4),
        "floss_aberrant_transcript_fraction": round(aberrant, 4),
        "n_transcripts_scored": len(floss_scores),
        "read_lengths": [int(x) for x in read_lengths],
        "floss_scores": [round(float(s), 4) for s in floss_scores],
    }
