"""
This script contains processing steps used to parse bam files.
"""
import pandas as pd
import numpy as np
import itertools
#temp test imports
import time

def process_reads(reads):
    """
    Process batches of reads from parse_bam, retrieving the data of interest
    and putting it in a dataframe.

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
                                                'count'])
    batch_df["reference_name"] = batch_df["reference_name"].astype("category")

    return batch_df


def process_sequences(sequences_counts: list,
                      pattern_length: int = 1,
                      max_sequence_length: int = None):
    """
    Calculate the occurence of nucleotides or groups of nucleotides in the
    sequences from the reads. The nucleotides or groups are stored in
    lexicographic order.
    """
    t0 = time.time() 

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
    sequence_array = np.zeros((num_sequences, max_sequence_length - pattern_length + 1, 4 ** pattern_length), dtype=int)
    
    # Populate the sequence array with counts for the corresponding nucleotide patterns
    for i, sequence in enumerate(sequences):
        for j in range(len(sequence) - pattern_length + 1):  # Adjusted range to include last position
            pattern = sequence[j:j + pattern_length]
            index = pattern_to_index(pattern)
            if index != -1:
                sequence_array[i, j, index] = 1

    """Following 2 blocks (background calculation) could be in its own function"""
    # Calculate the background frequency for three prime patterns
    condensed_arrays = {}
    three_prime_bg = np.copy(sequence_array)

    for i, sequence in enumerate(sequences):
        last_pattern_index = len(sequence) - pattern_length # changed because of flipped array
        three_prime_bg[i, last_pattern_index, :] = 0
        
    nucleotides = ["".join(nt) for nt in itertools.product('ACGT', repeat=pattern_length)]
    for nucleotide in nucleotides:
        nucleotide_counts = np.sum(three_prime_bg[:, :, pattern_to_index(nucleotide)])
        condensed_arrays[nucleotide] = nucleotide_counts
    total_bg_counts = sum(condensed_arrays.values())
    three_prime_bg = {k: v/total_bg_counts for k,v in condensed_arrays.items()}

    # Calculate the background frequency for five prime patterns
    condensed_arrays = {}
    five_prime_bg = np.copy(sequence_array)
    for i, sequence in enumerate(sequences):
        five_prime_bg[i, 0, :] = 0

    for nucleotide in nucleotides:
        nucleotide_counts = np.sum(five_prime_bg[:, :, pattern_to_index(nucleotide)])
        condensed_arrays[nucleotide] = nucleotide_counts
    total_bg_counts = sum(condensed_arrays.values())
    five_prime_bg = {k: v/total_bg_counts for k,v in condensed_arrays.items()}

    # Perform element-wise multiplication of sequence array and counts array
    result_array = sequence_array * counts_array[:, None, None]

    # Create the condensed 2D arrays for each nucleotide
    condensed_arrays = {}
    nucleotides = ["".join(nt) for nt in itertools.product('ACGT', repeat=pattern_length)]
    for nucleotide in nucleotides:
        nucleotide_counts = np.sum(result_array[:, :, pattern_to_index(nucleotide)], axis=0)
        condensed_arrays[nucleotide] = nucleotide_counts
    condensed_arrays["3_prime_bg"] = three_prime_bg
    condensed_arrays["5_prime_bg"] = five_prime_bg
    condensed_arrays["sequence_number"] = num_sequences
    print(f"completion time: {time.time() - t0}")

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


def calculate_background(sequence_array: np.array, sequences, pattern_length, five_prime: bool) -> dict:
    """
    
    """
    condensed_arrays = {}
    three_prime_bg = np.copy(sequence_array)
    for i, sequence in enumerate(sequences):
        last_pattern_index = len(sequence) - pattern_length
        three_prime_bg[i, last_pattern_index, :] = 0
        
    nucleotides = ["".join(nt) for nt in itertools.product('ACGT', repeat=pattern_length)]
    for nucleotide in nucleotides:
        nucleotide_counts = np.sum(three_prime_bg[:, :, pattern_to_index(nucleotide)])
        condensed_arrays[nucleotide] = nucleotide_counts
    total_bg_counts = sum(condensed_arrays.values())
    three_prime_bg = {k: v/total_bg_counts for k,v in condensed_arrays.items()}
    return


def join_batches(read_batches: list, full_sequence_batches: dict) -> tuple:
    """

    """
    print("\nGetting data from async objects..")

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
                        background_batches[pattern_length][pattern].append(array)
                else:
                    if pattern not in sequence_batches[pattern_length]:
                        sequence_batches[pattern_length][pattern] = [array]
                    else:
                        sequence_batches[pattern_length][pattern].append(array)
    
    print("Joining batch files..")

    read_df_pre = pd.concat(read_batches, ignore_index=True)
    read_df_pre["reference_name"] = (read_df_pre["reference_name"]
                                         .astype("category"))
    sequence_data = {}
    for pattern_length in sequence_batches:
        sequence_data[pattern_length] = {}
        for pattern in sequence_batches[pattern_length]:
            # Determine the maximum length among the arrays
            max_length = max(len(arr) for arr in sequence_batches[pattern_length][pattern])

            # Pad the arrays with zeros to match the maximum length
            padded_arrays = [np.pad(arr, (0, max_length - len(arr)), mode='constant') for arr in sequence_batches[pattern_length][pattern]]

            sequence_data[pattern_length][pattern] = np.sum(padded_arrays, axis=0)

    sequence_background = {}
    for pattern_length in background_batches:
        sequence_background[pattern_length] = {}
        for background in background_batches[pattern_length].keys():
        # Iterate over the patterns in the dictionaries
            if "sequence" in background:
                continue
            sequence_background[pattern_length][background] = {}

            for pattern in background_batches[pattern_length][background][0].keys():
                total_weighted_sum = 0
                total_count = 0

                # Calculate the weighted sum for the current pattern
                for i, dictionary in enumerate(background_batches[pattern_length][background]):
                    proportion = dictionary[pattern]
                    count = background_batches[pattern_length]["sequence_number"][i]
                    weighted_sum = proportion * count
                    total_weighted_sum += weighted_sum
                    total_count += count

                # Calculate the weighted average for the current key
                sequence_background[pattern_length][background][pattern] = total_weighted_sum / total_count

    return (read_df_pre, sequence_data, sequence_background)
