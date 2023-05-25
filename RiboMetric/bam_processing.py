"""
This script contains processing steps used to parse bam files.
"""
import pandas as pd
import numpy as np
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
    t0 = time.time()
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
    print(f"read process batch completion time: {time.time() - t0}")
    return batch_df


def process_sequences(sequences,
                      pattern_length=1,
                      sequence_length=50):
    """
    Calculate the occurence of nucleotides or groups of nucleotides in the
    sequences from the reads. The nucleotides or groups are stored in
    lexicographic order.
    """
    t0 = time.time()
    # Create an empty 2D array to store the counts
    counts_array = np.zeros((4 ** pattern_length,
                             sequence_length - pattern_length + 1
                             ), dtype=int)

    # Iterate over positions in the sequences
    for position in range(sequence_length - pattern_length + 1):
        # Get the nucleotides at the current position in all reads
        patterns = [sequence[position:position+pattern_length]
                    if position + pattern_length <= len(sequence)
                    else 0 for sequence in sequences]

        # Count the occurrences of each nucleotide pattern at
        # the current position
        counts = np.unique(patterns, return_counts=True)

        # Update the counts array
        for pattern, count in zip(counts[0], counts[1]):
            if pattern != 0:
                index = pattern_to_index(pattern)
                counts_array[index, position] = count
    print(f"sequence process batch completion time: {time.time() - t0}")
    return counts_array


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
