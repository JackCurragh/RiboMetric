"""
This script contains processing steps used to parse bam files.
"""
import pandas as pd
import numpy as np


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
                      pattern_length=1,
                      sequence_length=50):
    """
    Calculate the occurence of nucleotides or groups of nucleotides in the
    sequences from the reads. The nucleotides or groups are stored in
    lexicographic order.
    """
    read_names = [sequences[0] for sequences in sequences_counts]
    duplicate_counts = []
    for read in read_names:
        if "_x" in read:
            duplicate_counts.append(int(read.split("_x")[-1]))
        else:
            duplicate_counts.append(1)

    sequences = {}
    sequences["single"] = [sequences[1] for sequences in sequences_counts]
    sequences["duplicates"] = np.repeat(sequences["single"], np.array(duplicate_counts)-1)

    counts_array = {}
    for key in sequences:
        # Create an empty 2D array to store the counts
        counts_array[key] = np.zeros((4 ** pattern_length,
                                sequence_length - pattern_length + 1
                                ), dtype=int)

        # Iterate over positions in the sequences
        for position in range(sequence_length - pattern_length + 1):
            # Get the nucleotides at the current position in all reads
            patterns = [sequence[position:position+pattern_length]
                        if position + pattern_length <= len(sequence)
                        else 0 for sequence in sequences[key]]

            # Count the occurrences of each nucleotide pattern at
            # the current position
            counts = np.unique(patterns, return_counts=True)

            # Update the counts array
            for pattern, count in zip(counts[0], counts[1]):
                if pattern != 0:
                    index = pattern_to_index(pattern)
                    counts_array[key][index, position] = count
    
    background_frequency = np.delete(b,0,1)
    background_frequency = np.sum(background_frequency,axis=1)/np.sum(background_frequency)

    counts_array = np.add(counts_array["single"], counts_array["duplicates"])
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
