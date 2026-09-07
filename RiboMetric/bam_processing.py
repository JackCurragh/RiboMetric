"""
This script contains processing steps used to parse bam files.
"""
import io
import itertools
import os
from typing import Any, Dict, Iterable, List

import numpy as np
import oxbow as ox
import pandas as pd
import pyarrow.ipc
import pysam

from .file_splitting import format_progress, split_bam

# Process 1-in-N read chunks when accumulating sequence composition and
# terminal-bias backgrounds. 1 = use every read (accurate, the default);
# increase to sample for speed on very large inputs.
SEQUENCE_CHUNK_STRIDE = int(os.environ.get("RIBOMETRIC_SEQUENCE_CHUNK_STRIDE", "1"))


def validate_bam(bam_file: str) -> None:
    """
    Validate a bam file by attempting to read it with oxbow

    Inputs:
        bam_file: Path to the BAM file

    Outputs:
        None
    """
    try:
        ox.read_bam(bam_file)
    except Exception as e:
        if "InvalidReferenceSequenceName" in str(e):
            raise Exception("InvalidReferenceSequenceName - \
                            Likely an invalid character in sequence name \
                            eg. ), or ( ")
        else:
            raise Exception("Invalid bam file")


def ox_parse_reads(bam_file: str,
                   split_num: int,
                   reference_df: pd.DataFrame,
                   tempdir: str
                   ) -> tuple:
    """
    Splits a bam files using generated bed files, uses oxbow to process these
    batches of reads directly into a data frame and then processes the data
    from this generated dataframe into read and sequence data

    Inputs:
        bam_file: Path to the BAM file
        split_num: Number of the split
        reference_df: Reference dataframe, generated from samtools idxstats
        tempdir: Path to the temporary directory

    Outputs:
        tuple: A tuple containing:
            batch_df: Dataframe containing a processed batch of reads
            sequence_data: Dictionary containing the processed sequence data
    """
    formatted_num = f"{split_num+1:02d}"
    try:
        print_columns = os.get_terminal_size().columns // 25
    except Exception:
        print_columns = 4

    print("\n"*(split_num // print_columns),
          "\033[25C"*(split_num % print_columns),
          f"thread {formatted_num}: splitting.. | ",
          "\033[1A"*(split_num // print_columns),
          end="\r", flush=False, sep="")

    tmp_bam = split_bam(bam_file,
                        split_num,
                        reference_df,
                        tempdir)

    validate_bam(tmp_bam)

    print("\n"*(split_num // print_columns),
          "\033[25C"*(split_num % print_columns),
          f"thread {formatted_num}: parsing..   | ",
          "\033[1A"*(split_num // print_columns),
          end="\r", flush=False, sep="")

    oxbow_df = read_oxbow_df(tmp_bam)

    print("\n"*(split_num // print_columns),
          "\033[25C"*(split_num % print_columns),
          f"thread {formatted_num}: to pandas.. | ",
          "\033[1A"*(split_num // print_columns),
          end="\r", flush=False, sep="")

    return process_oxbow_batch(oxbow_df, split_num, formatted_num, print_columns)


def ox_parse_entire_bam(bam_file: str) -> tuple:
    """Parse a BAM directly without first creating a split BAM."""
    validate_bam(bam_file)
    oxbow_df = read_oxbow_df(bam_file)
    return process_oxbow_batch(oxbow_df)


def read_oxbow_df(bam_file: str) -> pd.DataFrame:
    try:
        arrow_ipc = ox.read_bam(bam_file)
    except Exception as e:
        if "InvalidReferenceSequenceName" in str(e):
            raise Exception("InvalidReferenceSequenceName - \
                            Likely an invalid character in sequence name \
                            eg. ), or ( ")
        else:
            raise Exception(f"Invalid bam file {e}")

    oxbow_df = pyarrow.ipc.open_file(io.BytesIO(arrow_ipc)).read_pandas()
    del arrow_ipc
    oxbow_df.attrs["bam_file"] = bam_file
    _recover_alignment_tags_from_pysam(oxbow_df, bam_file)
    return oxbow_df


def _recover_alignment_tags_from_pysam(oxbow_df: pd.DataFrame, bam_file: str) -> None:
    """Recover alignment tags that oxbow does not expose.

    Oxbow/Arrow treats BAM's 0xff sentinel as null, but STAR uses 255 as the
    unique-mapper code.  Recover the original integer values from pysam so
    unique-only filtering is based on real BAM MAPQ, not a fill value.  STAR
    transcriptome BAMs can also carry alternate-hit evidence in XA:i while all
    reported transcript rows have MAPQ=255, so recover XA for multimapper
    filtering and reporting.
    """
    if oxbow_df.empty:
        oxbow_df.attrs["mapq_recovered_from_pysam"] = False
        return

    recovered_mapq = []
    recovered_nh = []
    recovered_xa = []
    recovered_qname = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary:
                continue
            recovered_mapq.append(read.mapping_quality)
            recovered_nh.append(read.get_tag("NH") if read.has_tag("NH") else np.nan)
            recovered_xa.append(read.get_tag("XA") if read.has_tag("XA") else np.nan)
            recovered_qname.append(read.query_name)
            if len(recovered_mapq) >= len(oxbow_df):
                break

    if len(recovered_mapq) != len(oxbow_df):
        oxbow_df.attrs["mapq_recovered_from_pysam"] = False
        return

    # Row order is only trustworthy when oxbow and pysam enumerate primary reads
    # identically. That holds for STAR transcriptome BAMs (no secondaries), but
    # not necessarily for genome BAMs or BAMs with secondary alignments, where
    # oxbow's row order can diverge from ``fetch(until_eof=True)`` minus
    # secondaries. Verify per-read identity by qname before trusting the
    # positional join; bail out (leaving oxbow's values untouched) on any drift.
    name_column = next(
        (c for c in _READ_NAME_COLUMNS if c in oxbow_df.columns), None
    )
    if name_column is not None:
        oxbow_names = oxbow_df[name_column].astype(str).to_numpy()
        if not np.array_equal(oxbow_names, np.asarray(recovered_qname, dtype=str)):
            oxbow_df.attrs["mapq_recovered_from_pysam"] = False
            return

    if "mapq" in oxbow_df.columns:
        missing = oxbow_df["mapq"].isna()
        if missing.any():
            recovered_series = pd.Series(recovered_mapq, index=oxbow_df.index, dtype="float")
            oxbow_df.loc[missing, "mapq"] = recovered_series[missing]
            oxbow_df.attrs["mapq_recovered_from_pysam"] = True
        else:
            oxbow_df.attrs["mapq_recovered_from_pysam"] = False
    else:
        oxbow_df.attrs["mapq_recovered_from_pysam"] = False
    if "nh" not in oxbow_df.columns and any(not pd.isna(x) for x in recovered_nh):
        oxbow_df["nh"] = pd.Series(recovered_nh, index=oxbow_df.index, dtype="float")
    if "xa" not in oxbow_df.columns and any(not pd.isna(x) for x in recovered_xa):
        oxbow_df["xa"] = pd.Series(recovered_xa, index=oxbow_df.index, dtype="float")


def read_pysam_df(bam_file: str) -> pd.DataFrame:
    """Fallback BAM reader for oxbow empty/no-schema outputs."""
    records = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary:
                continue
            records.append({
                "qname": read.query_name,
                "seq": read.query_sequence or "",
                "cigar": read.cigarstring or "",
                "rname": bam.get_reference_name(read.reference_id),
                "pos": (
                    float(read.reference_start + 1)
                    if read.reference_start is not None else np.nan
                ),
                "mapq": read.mapping_quality,
                "nh": read.get_tag("NH") if read.has_tag("NH") else np.nan,
                "xa": read.get_tag("XA") if read.has_tag("XA") else np.nan,
            })
    return pd.DataFrame.from_records(
        records,
        columns=["qname", "seq", "cigar", "rname", "pos", "mapq", "nh", "xa"],
    )


def process_oxbow_batch(
        oxbow_df: pd.DataFrame,
        split_num: int = 0,
        formatted_num: str = "01",
        print_columns: int = 4,
        ) -> tuple:
    if oxbow_df.empty and not set(_REQUIRED_OXBOW_COLUMNS).issubset(oxbow_df.columns):
        oxbow_df = read_pysam_df(oxbow_df.attrs["bam_file"])
        if oxbow_df.empty:
            return (_empty_read_batch(), {1: [], 2: []})

    batch_df = process_reads(oxbow_df)

    print("\n"*(split_num // print_columns),
          "\033[25C"*(split_num % print_columns),
          f"thread {formatted_num}: sequencing..| ",
          "\033[1A"*(split_num // print_columns),
          end="\r", flush=False, sep="")

    sequence_data: Dict[int, list] = {1: [], 2: []}
    sequence_list = oxbow_df["seq"].tolist()
    count_list = batch_df["count"].tolist()

    # sequence_list batch size
    size = 10000
    list_length = len(sequence_list)

    if list_length < size and list_length != 0:
        size = list_length

    # Stride over the chunks of `size` reads used to accumulate sequence
    # composition / terminal-bias backgrounds. SEQUENCE_CHUNK_STRIDE = 1 uses
    # every chunk (all reads); a larger value samples 1-in-N chunks to trade
    # accuracy for speed on very large inputs. Reads are already capped by
    # --subsample, so the default processes all of them rather than the silent
    # 1-in-10 sample used previously.
    for pattern_length in sequence_data:
        count = -1
        progress = 0

        for i in range(0, len(sequence_list), size):
            count += 1
            if count % SEQUENCE_CHUNK_STRIDE != 0:
                continue
            section = sequence_list[i:i+size]
            counts = count_list[i:i+size]
            sequence_data[pattern_length].append(
                process_sequences(section,
                                  counts,
                                  pattern_length))

            progress += size
            formatted_progress = (format_progress((progress/list_length)*1000)
                                  if (progress/list_length)*1000 < 100
                                  else format_progress(100))
            print(
                "\n"*(split_num // print_columns),
                "\033[25C"*(split_num % print_columns),
                f"thread {formatted_num}: {pattern_length}: {formatted_progress}  | ",
                "\033[1A"*(split_num // print_columns),
                end="\r", flush=False, sep="")

    print("\n"*(split_num // print_columns),
          "\033[25C"*(split_num % print_columns),
          f"thread {formatted_num}: Parsed!     | ",
          "\033[1A"*(split_num // print_columns),
          end="\r", flush=False, sep="")

    return (batch_df, sequence_data)


_READ_NAME_COLUMNS = ("qname", "query_name", "read_name", "name")
_REQUIRED_OXBOW_COLUMNS = ("seq", "cigar", "rname", "pos", "mapq")


def _parse_nh_tag_value(value: Any) -> float:
    """Return NH as a float, or NaN when the tag is absent/unparseable."""
    if value is None:
        return np.nan
    if isinstance(value, float) and np.isnan(value):
        return np.nan
    if isinstance(value, (int, np.integer)):
        return float(value)
    if isinstance(value, (float, np.floating)):
        return float(value)
    if isinstance(value, dict):
        return _parse_nh_tag_value(value.get("NH") or value.get("nh"))
    if isinstance(value, (list, tuple)):
        for item in value:
            parsed = _parse_nh_tag_value(item)
            if not np.isnan(parsed):
                return parsed
        return np.nan

    text = str(value)
    if text.isdigit():
        return float(text)
    for marker in ("NH:i:", "NH:Z:", "NH="):
        if marker in text:
            tail = text.split(marker, 1)[1]
            digits = []
            for char in tail:
                if char.isdigit():
                    digits.append(char)
                elif digits:
                    break
            if digits:
                return float("".join(digits))
    return np.nan


def _parse_xa_value(value: Any) -> float:
    """Return the number of *alternative* loci encoded by an XA tag, or NaN.

    XA appears in two incompatible forms:

    * STAR transcriptome BAMs write ``XA:i:N`` — an integer count of alternative
      alignments, surfaced here directly as ``float(N)``.
    * BWA writes ``XA:Z:`` — a semicolon-terminated list of alternative loci
      (``chr,pos,CIGAR,NM;...``). The alt-locus count is the number of list
      entries.

    The downstream multimapper test is ``alt_loci > 0``, so both forms must
    reduce to a numeric count. Previously the string form was coerced to NaN,
    silently disabling XA-based multimapper detection on BWA BAMs.
    """
    if value is None:
        return np.nan
    if isinstance(value, float) and np.isnan(value):
        return np.nan
    if isinstance(value, (int, np.integer)):
        return float(value)
    if isinstance(value, (float, np.floating)):
        return float(value)
    if isinstance(value, (list, tuple)):
        # Already a list of alternative loci.
        return float(len([item for item in value if item not in (None, "")]))

    text = str(value).strip()
    # Strip a leading SAM tag prefix if present (``XA:Z:`` / ``XA:i:``).
    for marker in ("XA:Z:", "XA:i:", "XA="):
        if text.startswith(marker):
            text = text[len(marker):]
            break
    if text == "":
        return np.nan
    if text.lstrip("-").isdigit():
        # Integer count (STAR ``XA:i:``).
        return float(text)
    # BWA ``XA:Z:`` semicolon-terminated list of loci.
    return float(len([locus for locus in text.split(";") if locus]))


def _extract_nh_column(oxbow_df: pd.DataFrame) -> pd.Series:
    """Extract the SAM NH tag from common oxbow/pysam tag representations."""
    for column in ("nh", "NH", "tag_NH", "tags.NH", "aux_NH"):
        if column in oxbow_df.columns:
            return pd.to_numeric(oxbow_df[column], errors="coerce")
    for column in ("tags", "aux", "optional_fields"):
        if column in oxbow_df.columns:
            return oxbow_df[column].map(_parse_nh_tag_value)
    return pd.Series(np.nan, index=oxbow_df.index, dtype="float")


def _empty_read_batch() -> pd.DataFrame:
    return pd.DataFrame({
        "read_name": pd.Series(dtype="category"),
        "read_length": pd.Series(dtype="category"),
        "reference_name": pd.Series(dtype="category"),
        "reference_start": pd.Series(dtype="float"),
        "soft_clip_5": pd.Series(dtype="uint8"),
        "first_dinucleotide": pd.Series(dtype="category"),
        "last_dinucleotide": pd.Series(dtype="category"),
        "count": pd.Series(dtype="category"),
        "mapq": pd.Series(dtype="uint8"),
        "mapq_available": pd.Series(dtype="bool"),
        "mapq_recovered_from_pysam": pd.Series(dtype="bool"),
        "nh": pd.Series(dtype="float"),
        "xa": pd.Series(dtype="float"),
    })


def _get_oxbow_column(
    oxbow_df: pd.DataFrame,
    candidates: Iterable[str],
    label: str,
) -> pd.Series:
    for column in candidates:
        if column in oxbow_df.columns:
            return oxbow_df[column]
    if oxbow_df.empty:
        return pd.Series(dtype="object")
    raise KeyError(
        f"Missing {label} column in oxbow BAM output. "
        f"Available columns: {list(oxbow_df.columns)}"
    )


def process_reads(oxbow_df: pd.DataFrame) -> pd.DataFrame:
    """
    Process batches of reads from parse_bam, retrieving the data of interest
    and putting it in a dataframe.
    Ensure category columns are set to category type for memory efficiency.

    Inputs:
        oxbow_df: List of read contents from bam files, returned by pysam

    Outputs:
        batch_df: Dataframe containing a processed batch of reads
    """
    if oxbow_df.empty:
        return _empty_read_batch()

    batch_df = pd.DataFrame()
    read_names = _get_oxbow_column(oxbow_df, _READ_NAME_COLUMNS, "read name")
    batch_df['read_name'] = read_names.astype("category")

    # Compute read length from CIGAR as the number of read bases consumed by the
    # alignment (M, =, X, I). This is robust to reference-side gaps (D) and avoids
    # counting soft-clipped bases as part of the footprint.
    # Falls back to SEQ length if CIGAR is missing/malformed.
    def _read_len_from_cigar(cigar: str, seqlen: int) -> int:
        try:
            parts = [] if not isinstance(cigar, str) else cigar
            # Find all (length, op) pairs
            pairs = [] if not parts else [
                (int(n), op) for n, op in __import__("re").findall(r"(\d+)([MIDNSHP=XB])", parts)
            ]
            if not pairs:
                return seqlen
            return int(sum(n for n, op in pairs if op in ("M", "=", "X", "I")))
        except Exception:
            return seqlen

    seq_len_series = oxbow_df["seq"].str.len().fillna(0).astype(int)
    cigar_series = oxbow_df["cigar"].astype(str)
    rl_series = [
        _read_len_from_cigar(cig, seqlen)
        for cig, seqlen in zip(cigar_series.tolist(), seq_len_series.tolist())
    ]
    batch_df["read_length"] = pd.Series(rl_series, dtype="category")

    batch_df["reference_name"] = oxbow_df["rname"].astype("category")

    # Correct reference_start for 5' soft clips so it points to the true 5' end
    # of the read (first base of the protected fragment) rather than the first
    # aligned base. This is required for accurate A-site/P-site offset detection.
    # CIGAR pattern: optional leading hard-clip (\d+H), then optional soft-clip (\d+S)
    soft_clip_5 = (
        oxbow_df["cigar"]
        .str.extract(r"^(?:\d+H)?(\d+)S", expand=False)
        .fillna("0")
        .astype(int)
    )
    # Use float to preserve NaN for reads with no mapped position; NaN - int = NaN
    batch_df["reference_start"] = oxbow_df["pos"].astype("float") - soft_clip_5
    # Store the 5' soft-clip count so downstream code can compute soft-clip rates.
    batch_df["soft_clip_5"] = soft_clip_5.astype("uint8")
    batch_df["first_dinucleotide"] = (oxbow_df["seq"].str.slice(stop=2)
                                      .astype("category"))
    batch_df["last_dinucleotide"] = (oxbow_df["seq"].str.slice(stop=-3,
                                                               step=-1)
                                     .astype("category"))
    batch_df["count"] = pd.Series([int(query.split("_x")[-1]) if "_x" in query
                                   else 1 for query in read_names],
                                  dtype="category")
    mapq_values = pd.to_numeric(oxbow_df["mapq"], errors="coerce")
    mapq_available = mapq_values.notna()
    batch_df["mapq"] = mapq_values.fillna(0).astype("uint8")
    batch_df["mapq_available"] = mapq_available.astype("bool")
    batch_df["mapq_recovered_from_pysam"] = bool(
        oxbow_df.attrs.get("mapq_recovered_from_pysam", False)
    )
    batch_df["nh"] = _extract_nh_column(oxbow_df).astype("float")
    batch_df["xa"] = (
        oxbow_df["xa"].map(_parse_xa_value).astype("float")
        if "xa" in oxbow_df.columns
        else pd.Series(np.nan, index=oxbow_df.index, dtype="float")
    )
    return batch_df


def process_sequences(sequences: list,
                      counts: list,
                      pattern_length: int = 1,
                      max_sequence_length: int = -1,
                      ) -> dict:
    """
    Calculate the occurence of nucleotides patterns in the sequences from
    the reads. The nucleotides patterns are stored in lexicographic order
    (see pattern to index)

    Inputs:
        sequences_counts: List of tuples containing read_name and sequence
        pattern_length: Length of the nucleotide pattern
        (e.g. 1: [A,C,G,T], 2: [AA,AC,AG,..,GT,TT])
        max_sequence_length: Manually set the max sequence length, sequences
        will be cut to this length. If None, takes the max found sequence
        length in the list of sequences

    Outputs:
        condensed_arrays: Dictionary containing raw pattern counts, 5' and 3'
        background frequencies and number of sequences in the batch (used
        later for joining of background frequencies)
    """
    # Create the counts array
    counts_array = np.array(counts)

    # Set sequences and calculate array dimensions
    num_sequences = len(sequences)
    if max_sequence_length == -1:
        max_sequence_length = max(len(seq) for seq in sequences)

    if max_sequence_length < pattern_length:
        return {}

    # Create the 3D numpy array with zeros
    sequence_array = np.zeros((num_sequences,
                               max_sequence_length - pattern_length + 1,
                               4 ** pattern_length),
                              dtype=int)

    # Populate the sequence array with counts for the corresponding
    # nucleotide patterns
    for i, sequence in enumerate(sequences):
        for j in range(len(sequence) - pattern_length + 1):
            pattern = sequence[j:j + pattern_length]
            index = pattern_to_index(pattern)
            if index != -1:
                sequence_array[i, j, index] = 1
    if pattern_length == 2:
        # Calculate background frequencies
        three_prime_bg = calculate_background(sequence_array,
                                              sequences,
                                              pattern_length,
                                              five_prime=False)
        five_prime_bg = calculate_background(sequence_array,
                                             sequences,
                                             pattern_length,
                                             five_prime=True)

    condensed_arrays = {}

    if pattern_length == 1:
        # Perform element-wise multiplication of sequence array
        # and counts array
        result_array = sequence_array * counts_array[:, None, None]
        # Create the condensed 2D arrays for each nucleotide

        nucleotides = ["".join(nt) for nt in
                       itertools.product('ACGT', repeat=pattern_length)]
        for nucleotide in nucleotides:
            nucleotide_counts = np.sum(
                result_array[:, :, pattern_to_index(nucleotide)],
                axis=0)
            condensed_arrays[nucleotide] = nucleotide_counts

    # Add backgrounds and sequence_number to output dictionary
    if pattern_length == 2:
        condensed_arrays["3_prime_bg"] = three_prime_bg
        condensed_arrays["5_prime_bg"] = five_prime_bg
        condensed_arrays["sequence_number"] = num_sequences

    return condensed_arrays


def pattern_to_index(pattern: str) -> int:
    """
    Converts a nucleotide pattern to its corresponding index in
    the counts array. Ensure A,C,G,T ordered array.
    (i.e. AA, AC, AG, AT, CA... TG, TT)
    """
    index = 0
    base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for nucleotide in pattern:
        if nucleotide in base_to_index:
            index = index * 4 + base_to_index[nucleotide]
        else:
            return -1
    return index


def calculate_background(sequence_array: np.ndarray,
                         sequences: List[str],
                         pattern_length: int,
                         five_prime: bool
                         ) -> Dict[str, float]:
    """
    Calculate the background frequency for a list of sequences. The background
    frequency is the proportion of nucleotide patterns without the first or
    last pattern in the read, for five prime and three prime respectively.

    Inputs:
        sequence_array: 3D array of a batch of sequences
        sequences: list of sequences from a batch
        pattern_length: The length of nucleotide patterns being processed
        five_prime: If set to True, returns the 'five_prime_bg' background,
        else returns the 'three_prime_bg' background

    Outputs:
        sequence_bg: A dictionary with the nucleotide pattern as keys and
        their background proportion as values
    """
    condensed_arrays = {}
    sequence_bg = np.copy(sequence_array)

    # The background is the pattern composition of the read *excluding* the
    # terminal pattern whose enrichment we are testing. For the 5' end that is
    # the first position (index 0); for the 3' end it is the last populated
    # position of each read. Previously this always zeroed index 0, which made
    # the 5' and 3' backgrounds identical and gave the 3' ligation-bias metric
    # the wrong reference distribution.
    n_positions = sequence_bg.shape[1]
    for i, sequence in enumerate(sequences):
        if five_prime:
            sequence_bg[i, 0, :] = 0
        else:
            # Last valid pattern start for this read (patterns are length
            # pattern_length; a read of length L has L-pattern_length+1 starts).
            last_pos = len(sequence) - pattern_length
            if 0 <= last_pos < n_positions:
                sequence_bg[i, last_pos, :] = 0

    nucleotides = ["".join(nt) for nt in
                   itertools.product('ACGT', repeat=pattern_length)]
    for nucleotide in nucleotides:
        nucleotide_counts = np.sum(sequence_bg[:, :,
                                               pattern_to_index(nucleotide)])
        condensed_arrays[nucleotide] = nucleotide_counts
    total_bg_counts = sum(condensed_arrays.values())
    return {k: v/total_bg_counts for k, v in condensed_arrays.items()}


def join_batches(bam_batches: list) -> tuple:
    """
    Get and join the data returned from multiprocessed_batches

    Inputs:
        read_batches: List of dataframes containing read information returned
                    from multiprocessed batches
        full_sequence_batches: Dictionary containing sequence data (counts per
                    position and background) returned from multiprocessed
                    batches

    Outputs:
        tuple containing:
            read_df_pre: The read dataframe containing read information before
                        further modifications to the dataframe
            sequence_data: Dictionary containing the total counts of
                        nucleotide patterns per nucleotide position
            sequence_background: Dictionary containing the background
                        frequency of nucleotide patterns for five and
                        three prime
    """
    print("\nGetting data from async objects..")
    read_batches, background_batches, sequence_batches = \
        get_batch_data(bam_batches)

    print("Joining batch files..")
    # Joining reads
    read_df_pre = pd.concat(read_batches, ignore_index=True)
    category_columns = ["read_length",
                        "reference_name",
                        "first_dinucleotide",
                        "last_dinucleotide",
                        "count"]
    read_df_pre[category_columns] = (read_df_pre[category_columns]
                                     .astype("category"))
    # Joining sequence data
    sequence_data = {}

    for pattern in sequence_batches:
        # Determine the maximum length among the arrays
        max_length = max(len(arr) for arr in
                         sequence_batches[pattern])

        # Pad the arrays with zeros to match the maximum length
        padded_arrays = [np.pad(arr, (0, max_length - len(arr)),
                                mode='constant') for arr in
                         sequence_batches[pattern]]

        sequence_data[pattern] = np.sum(padded_arrays,
                                        axis=0)
    # Joining sequence backgrounds
    sequence_background: Dict = {}

    for background in background_batches.keys():
        if background == "sequence_number":
            continue

        sequence_background[background] = {}
        iterable = background_batches[
            background][0].keys()
        for pattern in iterable:
            total_weighted_sum = 0
            total_count = 0

        # Calculate the weighted sum for the current pattern
            sum_iter = background_batches[background]
            for i, dictionary in enumerate(sum_iter):
                proportion = dictionary[pattern]
                count = background_batches[
                    "sequence_number"][i]
                weighted_sum = proportion * count
                total_weighted_sum += weighted_sum
                total_count += count

            # No sequence content (e.g. a BAM with no stored sequences):
            # skip rather than dividing by zero. See issue #122.
            if total_count == 0:
                continue
            # Calculate the weighted average for the current key
            sequence_background[background][pattern] = \
                total_weighted_sum / total_count

    # If no background could be computed (sequence-less BAM), return empty
    # dicts so sequence-based metrics are cleanly skipped downstream.
    if not any(sequence_background.get(bg) for bg in sequence_background):
        sequence_data, sequence_background = {}, {}

    return (read_df_pre, sequence_data, sequence_background)


def get_batch_data(
        bam_batches: list
        ) -> tuple:
    """
    Return readable data from the multiprocessed pools, separating the
    full sequence data into backgrounds data and sequence data.
    Called in the join_batches function

    Inputs:
        read_batches: List of dataframes containing read information returned
                    from multiprocessed batches
        full_sequence_batches: Dictionary containing sequence data (counts per
                    position and background) returned from multiprocessed

    Outputs:
        tuple containing:
            read_batches: List of dataframes containing read information
            background_batches: Dictionary containing background data
            sequence_batches: Dictionary containing sequence data
    """
    if type(bam_batches[0]) is pd.DataFrame:
        read_batches = [bam_batches[0]]
        sequence_data: Dict = {}
        full_sequence_batches = [sequence_data]
        for pattern_length in bam_batches[1].keys():
            sequence_data[pattern_length] = [result.get() for result
                                             in bam_batches[1][pattern_length]
                                             ]

    elif isinstance(bam_batches[0], tuple):
        bam_tuples = bam_batches

        read_batches = [data[0] for data in bam_tuples]
        full_sequence_batches = [data[1] for data in bam_tuples]

    else:
        bam_tuples = [result.get() for result in bam_batches]

        read_batches = [data[0] for data in bam_tuples]
        full_sequence_batches = [data[1] for data in bam_tuples]

    background_batches, sequence_batches = {}, {}
    for pattern_length in full_sequence_batches[0].keys():

        for full_batch in full_sequence_batches:

            for result in full_batch[pattern_length]:
                result_dict = result

                for pattern, array in result_dict.items():
                    if "bg" in pattern or "sequence" in pattern:
                        if pattern not in background_batches:
                            background_batches[pattern] = [array]
                        else:
                            (background_batches[pattern]
                             .append(array))
                    else:
                        if pattern not in sequence_batches:
                            sequence_batches[pattern] = [array]
                        else:
                            (sequence_batches[pattern]
                             .append(array))

    return read_batches, background_batches, sequence_batches
