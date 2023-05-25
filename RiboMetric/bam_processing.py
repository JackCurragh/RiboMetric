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
                count,             # count
            ]
        )
    batch_df = pd.DataFrame(read_list, columns=['read_length',
                                                'reference_name',
                                                'reference_start',
                                                'count'])
    batch_df["reference_name"] = batch_df["reference_name"].astype("category")

    return batch_df


def process_sequences(sequences_counts,
                      pattern_length=1):
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

    """Following 2 blocks could be in its own function"""
    # Calculate the background frequency for three prime patterns
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
