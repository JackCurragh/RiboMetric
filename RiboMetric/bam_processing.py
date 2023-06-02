"""
This script contains processing steps used to parse bam files.
"""
import pandas as pd
import numpy as np
import itertools


def process_reads(reads):
    """
    Process batches of reads from parse_bam, retrieving the data of interest
    and putting it in a dataframe.
    Ensure category columns are set to category type for memory efficiency.

    Inputs:
        reads: List of read contents from bam files, returned by pysam

    Outputs:
        batch_df: Dataframe containing a processed batch of reads
    """
    read_list = []
    for read in reads:
        if "_x" in read[0]:
            count = int(read[0].split("_x")[-1])
        else:
            count = 1
        read_list.append(
            [
                len(read[9]),      # read_length
                read[2],           # reference_name
                int(read[3]),      # reference_start
                read[9][0:2],      # first dinucleotide
                read[9][-2:],      # last dinucleotide
                count,             # count
            ]
        )
    batch_df = pd.DataFrame(read_list, columns=['read_length',
                                                'reference_name',
                                                'reference_start',
                                                'first_dinucleotide',
                                                'last_dinucleotide',
                                                'count'])

    batch_df["first_dinucleotide"] = (batch_df["first_dinucleotide"]
                                      .astype("category"))
    batch_df["last_dinucleotide"] = (batch_df["last_dinucleotide"]
                                     .astype("category"))
    batch_df["reference_name"] = (batch_df["reference_name"]
                                  .astype("category"))

    return batch_df


def process_sequences(sequences_counts: list,
                      pattern_length: int = 1,
                      max_sequence_length: int = None,
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
    read_names = [sequences[0] for sequences in sequences_counts]
    counts_array = []
    for read in read_names:
        if "_x" in read:
            counts_array.append(int(read.split("_x")[-1]))
        else:
            counts_array.append(1)
    counts_array = np.array(counts_array)

    # Set sequences and calculate array dimensions
    sequences = [sequences[1] for sequences in sequences_counts]
    num_sequences = len(sequences)
    if max_sequence_length is None:
        max_sequence_length = max(len(seq) for seq in sequences)

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

    # Calculate background frequencies
    three_prime_bg = calculate_background(sequence_array,
                                          sequences,
                                          pattern_length,
                                          five_prime=False)
    five_prime_bg = calculate_background(sequence_array,
                                         sequences,
                                         pattern_length,
                                         five_prime=True)

    # Perform element-wise multiplication of sequence array and counts array
    result_array = sequence_array * counts_array[:, None, None]
    # Create the condensed 2D arrays for each nucleotide
    condensed_arrays = {}
    nucleotides = ["".join(nt) for nt in
                   itertools.product('ACGT', repeat=pattern_length)]
    for nucleotide in nucleotides:
        nucleotide_counts = np.sum(result_array[:, :,
                                                pattern_to_index(nucleotide)],
                                   axis=0)
        condensed_arrays[nucleotide] = nucleotide_counts

    # Add backgrounds and sequence_number to output dictionary
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
    if not five_prime:
        # Flip array so 3' is at the start
        flipped_arr = np.flip(sequence_bg, axis=1)

        # Move rows containing only zeros to the end of each matrix
        nonzero_mask = np.any(flipped_arr != 0, axis=2)
        sequence_bg = np.zeros_like(flipped_arr)
        for i in range(flipped_arr.shape[0]):
            nonzeros = flipped_arr[i, nonzero_mask[i]]
            zeros = flipped_arr[i, ~nonzero_mask[i]]
            sequence_bg[i] = np.concatenate((nonzeros, zeros), axis=0)

    for i, sequence in enumerate(sequences):
        sequence_bg[i, 0, :] = 0

    nucleotides = ["".join(nt) for nt in
                   itertools.product('ACGT', repeat=pattern_length)]
    for nucleotide in nucleotides:
        nucleotide_counts = np.sum(sequence_bg[:, :,
                                               pattern_to_index(nucleotide)])
        condensed_arrays[nucleotide] = nucleotide_counts
    total_bg_counts = sum(condensed_arrays.values())
    sequence_bg = {k: v/total_bg_counts for k, v in condensed_arrays.items()}
    return sequence_bg


def join_batches(read_batches: list, full_sequence_batches: dict) -> tuple:
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
        get_batch_data(read_batches, full_sequence_batches)

    print("Joining batch files..")
    # Joining reads
    read_df_pre = pd.concat(read_batches, ignore_index=True)
    category_columns = ["reference_name",
                        "first_dinucleotide",
                        "last_dinucleotide"]
    read_df_pre[category_columns] = (read_df_pre[category_columns]
                                     .astype("category"))
    # Joining sequence data
    sequence_data = {}
    for pattern_length in sequence_batches:
        sequence_data[pattern_length] = {}
        for pattern in sequence_batches[pattern_length]:
            # Determine the maximum length among the arrays
            max_length = max(len(arr) for arr in
                             sequence_batches[pattern_length][pattern])

            # Pad the arrays with zeros to match the maximum length
            padded_arrays = [np.pad(arr, (0, max_length - len(arr)),
                                    mode='constant') for arr in
                             sequence_batches[pattern_length][pattern]]

            sequence_data[pattern_length][pattern] = np.sum(padded_arrays,
                                                            axis=0)
    # Joining sequence backgrounds
    sequence_background = {}
    for pattern_length in background_batches:
        sequence_background[pattern_length] = {}
        # Iterate over the patterns in the dictionaries
        for background in background_batches[pattern_length].keys():
            if background == "sequence_number":
                continue

            sequence_background[pattern_length][background] = {}
            iterable = background_batches[pattern_length][
                background][0].keys()
            for pattern in iterable:
                total_weighted_sum = 0
                total_count = 0

            # Calculate the weighted sum for the current pattern
                sum_iter = background_batches[pattern_length][background]
                for i, dictionary in enumerate(sum_iter):
                    proportion = dictionary[pattern]
                    count = background_batches[pattern_length][
                        "sequence_number"][i]
                    weighted_sum = proportion * count
                    total_weighted_sum += weighted_sum
                    total_count += count

                # Calculate the weighted average for the current key
                sequence_background[pattern_length][background][pattern] = \
                    total_weighted_sum / total_count

    return (read_df_pre, sequence_data, sequence_background)


def get_batch_data(
        read_batches: list,
        full_sequence_batches: dict
        ) -> tuple(
                list, dict, dict
        ):
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
    read_batches = [result.get() for result in read_batches]

    background_batches, sequence_batches = {}, {}
    for pattern_length in full_sequence_batches:
        background_batches[pattern_length] = {}
        sequence_batches[pattern_length] = {}

        for result in full_sequence_batches[pattern_length]:
            result_dict = result.get()

            for pattern, array in result_dict.items():
                if "bg" in pattern or "sequence" in pattern:
                    if pattern not in background_batches[pattern_length]:
                        background_batches[pattern_length][pattern] = [array]
                    else:
                        (background_batches[pattern_length][pattern]
                         .append(array))
                else:
                    if pattern not in sequence_batches[pattern_length]:
                        sequence_batches[pattern_length][pattern] = [array]
                    else:
                        (sequence_batches[pattern_length][pattern]
                         .append(array))

    return read_batches, background_batches, sequence_batches
