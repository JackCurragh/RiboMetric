"""
This script contains processing steps used to parse bam files.
"""
import pandas as pd
import numpy as np
import itertools
import oxbow as ox
import io
import os
import pyarrow.ipc
from .file_splitting import split_bam, format_progress
from multiprocessing import Pool
from typing import Dict


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

    try:
        arrow_ipc = ox.read_bam(tmp_bam)
    except Exception as e:
        if "InvalidReferenceSequenceName" in str(e):
            raise Exception("InvalidReferenceSequenceName - \
                            Likely an invalid character in sequence name \
                            eg. ), or ( ")
        else:
            raise Exception(f"Invalid bam file {e}")

    oxbow_df = pyarrow.ipc.open_file(io.BytesIO(arrow_ipc)).read_pandas()
    del arrow_ipc

    print("\n"*(split_num // print_columns),
          "\033[25C"*(split_num % print_columns),
          f"thread {formatted_num}: to pandas.. | ",
          "\033[1A"*(split_num // print_columns),
          end="\r", flush=False, sep="")

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

    for pattern_length in sequence_data:
        count = -1
        progress = 0

        for i in range(0, len(sequence_list), size):
            count += 1
            if count % 10 != 0:
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


def ox_server_parse_reads(bam_file: str,
                          num_processes: int = 4
                          ) -> tuple:
    """
    Functionally the same as ox_parse_reads,
    but without splitting the bam file.

    Inputs:
        bam_file: Path to the BAM file
        pool: Multiprocessing pool

    Outputs:
        tuple: A tuple containing:
            batch_df: Dataframe containing a processed batch of reads
            sequence_data: Dictionary containing the processed sequence data
    """
    pool = Pool(processes=num_processes)
    print("Running in server mode")
    print("Generating pyarrow object")
    arrow_ipc = ox.read_bam(bam_file)
    print("Transforming to pandas df")
    oxbow_df = pyarrow.ipc.open_file(io.BytesIO(arrow_ipc)).read_pandas()
    del arrow_ipc
    read_df = process_reads(oxbow_df)
    print("retrieving sequence data")
    sequence_data: Dict[int, list] = {1: [], 2: []}
    sequence_list = oxbow_df["seq"].tolist()
    count_list = read_df["count"].tolist()
    del oxbow_df
    # sequence_list batch size
    size = 10000
    if len(sequence_list) < size and len(sequence_list) != 0:
        size = len(sequence_list)
    for pattern_length in sequence_data:
        count = -1
        for i in range(0, len(sequence_list), size):
            count += 1
            if count % 10 != 0:
                continue
            section = sequence_list[i:i+size]
            counts = count_list[i:i+size]
            sequence_data[pattern_length].append(
                pool.apply_async(
                    process_sequences,
                    [section, counts, pattern_length]
                    )
                )

    pool.close()
    pool.join()

    return (read_df, sequence_data)


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
    batch_df = pd.DataFrame()
    batch_df["read_length"] = pd.Series(oxbow_df["end"] - oxbow_df["pos"] + 1,
                                        dtype="category")
    batch_df["reference_name"] = oxbow_df["rname"].astype("category")
    batch_df["reference_start"] = oxbow_df["pos"].dropna().astype("int")
    batch_df["first_dinucleotide"] = (oxbow_df["seq"].str.slice(stop=2)
                                      .astype("category"))
    batch_df["last_dinucleotide"] = (oxbow_df["seq"].str.slice(stop=-3,
                                                               step=-1)
                                     .astype("category"))
    batch_df["count"] = pd.Series([int(query.split("_x")[-1]) if "_x" in query
                                   else 1 for query in oxbow_df["qname"]],
                                  dtype="category")
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
            return 0
    return index


def calculate_background(sequence_array: np.array,
                         sequences,
                         pattern_length,
                         five_prime: bool
                         ) -> dict:
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

    for i, sequence in enumerate(sequences):
        sequence_bg[i, 0, :] = 0

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

            # Calculate the weighted average for the current key
            sequence_background[background][pattern] = \
                total_weighted_sum / total_count

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
