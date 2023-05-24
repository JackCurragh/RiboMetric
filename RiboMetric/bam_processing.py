def process_reads(reads):
    """
    Process batches of reads from parse_bam, retrieving the data of interest and putting it in a dataframe.

    Inputs:
        reads:

    Outputs:
        batch_df:
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
                read[9],           # sequence
                count,             # count
            ]
        )
    batch_df = pd.DataFrame(read_list, columns=['read_length',
                                    'reference_name', 'reference_start',
                                    'count'])  # Convert to DataFrame
    batch_df["reference_name"] = batch_df["reference_name"].astype("category")
    return batch_df


def process_sequences(sequences, pattern_length=1, sequence_length = 50):
    """
    Calculate the occurence of nucleotides or groups of nucleotides in the sequences from the reads.
    The nucleotides or groups are stored in lexicographic order, (i.e. AA, AC, AG, AT, CA... TG, TT)
    """
    # Create an empty 2D array to store the counts
    counts_array = np.zeros((4 ** pattern_length, sequence_length - pattern_length + 1), dtype=int)

    # Iterate over each position in the sequences
    for i in range(sequence_length - pattern_length + 1):
        # Get the nucleotides at the current position
        patterns = [sequence[i:i+pattern_length] if i + pattern_length <= len(sequence) else 0 for sequence in sequences]

        # Count the occurrences of each nucleotide pattern at the current position
        counts = np.unique(patterns, return_counts=True)

        # Update the counts array
        for pattern, count in zip(counts[0], counts[1]):
            if pattern is not None:
                index = pattern_to_index(pattern)
                counts_array[index, i] = count

    return counts_array

def pattern_to_index(pattern):
    """
    Converts a nucleotide pattern to its corresponding index in the counts array.
    """
    index = 0
    base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for nucleotide in pattern:
        if nucleotide in base_to_index:
            index = index * 4 + base_to_index[nucleotide]
        else:
            return 0
    return index