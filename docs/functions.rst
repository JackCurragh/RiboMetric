Table of Contents
=================

.. contents::
    :local:
    :backlinks: top
    :depth: 2
    
..
    Lukas: Doesn't seem to be working (or at least not on github), use github's sections at the top right of this viewer.
    
File Parsing (file_parser.py)
=================

parse_bam
---------------

Reads in a BAM file and performs processing on the reads. It uses multiprocessing to improve performance by parallelizing the processing. The function takes several parameters, including the path to the BAM file, the number of reads to parse, the batch size for processing reads, and the number of processes to create. The function returns a tuple containing three elements: `read_df_pre`, `sequence_data`, and `sequence_background`.

Parameters:
    - ``bam_file`` (str): Path to the BAM file.
    - ``num_reads`` (int, optional): Number of reads to parse (default: 1000000).
    - ``batch_size`` (int, optional): The number of reads processed at a time (default: 10000).
    - ``num_processes`` (int, optional): The maximum number of processes to create (default: 4).

Returns:
    - ``parsed_bam`` (tuple): Tuple containing:
    
        - ``read_df_pre`` (pd.DataFrame): The read dataframe containing read information before further modifications to the dataframe.
        - ``sequence_data`` (dict): Dictionary containing the total counts of nucleotide patterns per nucleotide position.
        - ``sequence_background`` (dict): Dictionary containing the background frequency of nucleotide patterns for five and three prime.

The function starts by opening the BAM file using `pysam.AlignmentFile` and creates a multiprocessing pool with the specified number of processes. It then iterates through the reads in the BAM file, appending each read to a list. If the batch size is reached or the specified number of reads is exceeded, the current batch of reads is processed using the `process_reads` and `process_sequences` functions, and the results are stored in the respective batches. The read list is then cleared for the next batch.

After processing all the reads, the multiprocessing pool is closed and joined. The `join_batches` function is called to combine the read batches and sequence batches into the `parsed_bam` tuple, which is returned by the function.

----

process_reads
---------------

Process batches of reads from the `parse_bam` function, retrieving the data of interest and putting it in a dataframe.

Parameters:
    - ``reads`` (list): List of read contents from BAM files, returned by `pysam`.

Returns:
    - ``batch_df`` (pd.DataFrame): Dataframe containing a processed batch of reads.

The function iterates through the provided list of reads and extracts the relevant information. For each read, it checks if the read name contains "_x" to determine the count of the read. It then creates a list with the following columns:

    - ``read_length``: Length of the read.
    - ``reference_name``: Name of the reference sequence where the read is aligned.
    - ``reference_start``: Start position of the read on the reference sequence.
    - ``first_dinucleotide``: First two nucleotides of the read sequence.
    - ``last_dinucleotide``: Last two nucleotides of the read sequence.
    - ``count``: Count of the read (either from the read name or default to 1).

The list of processed reads is used to create a Pandas DataFrame called ``batch_df`` with the specified column names. The columns ``first_dinucleotide``, ``last_dinucleotide``, and ``reference_name`` are cast to the categorical data type in the DataFrame.

The resulting DataFrame ``batch_df`` is returned by the function.

----

process_sequences
---------------

Calculate the occurrence of nucleotide patterns in the sequences from the reads. The nucleotide patterns are stored in lexicographic order. The function takes a list of tuples containing read names and sequences, the length of the nucleotide pattern, and the maximum sequence length. It returns a dictionary containing the raw pattern counts, 5' and 3' background frequencies, and the number of sequences in the batch.

Parameters:
    - ``sequences_counts`` (list): List of tuples containing read names and sequences.
    - ``pattern_length`` (int, optional): Length of the nucleotide pattern (default: 1).
    - ``max_sequence_length`` (int, optional): Manually set the maximum sequence length (default: None).

Returns:
    - ``condensed_arrays`` (dict): Dictionary containing raw pattern counts, 5' and 3' background frequencies, and the number of sequences in the batch.

The function starts by creating a counts array based on the read names. It checks if the read name contains "_x" to determine the count of the read. It then creates a 3D numpy array called ``sequence_array`` with zeros to store the counts for each nucleotide pattern. The array dimensions depend on the number of sequences, the maximum sequence length, and the pattern length.

Next, the function populates the ``sequence_array`` by iterating through the sequences and nucleotide positions. It extracts the nucleotide pattern, converts it to an index using the ``pattern_to_index`` helper function, and updates the count in the corresponding position of the array.

The function calculates the 5' and 3' background frequencies using the ``calculate_background`` helper function. The resulting ``sequence_array`` is multiplied element-wise by the counts array to obtain the final result array. The raw pattern counts for each nucleotide are condensed into 2D arrays, and the background frequencies and sequence count are added to the output dictionary.

----

pattern_to_index
~~~~~~~~~~~~~~~~~~

Converts a nucleotide pattern to its corresponding index in the counts array. The function ensures an A, C, G, T-ordered array.

Parameters:
    - ``pattern`` (str): Nucleotide pattern.

Returns:
    - ``index`` (int): Index of the nucleotide pattern in the counts array.

The function iterates through each nucleotide in the pattern, converts it to the corresponding index, and calculates the final index value based on the nucleotide order. If a nucleotide is not found in the base-to-index dictionary, the index is set to 0.

----

calculate_background
~~~~~~~~~~~~~~~~~~

Calculate the background frequency for a list of sequences. The background frequency is the proportion of nucleotide patterns without the first or last pattern in the read, for the 5' and 3' ends, respectively.

Parameters:
    - ``sequence_array`` (np.array): 3D array of a batch of sequences.
    - ``sequences`` (list): List of sequences from a batch.
    - ``pattern_length`` (int): Length of nucleotide patterns being processed.
    - ``five_prime`` (bool): If set to True, returns the 5' background; otherwise, returns the 3' background.

Returns:
    - ``sequence_bg`` (dict): A dictionary with the nucleotide pattern as keys and their background proportion as values.

The function initializes a copy of the ``sequence_array`` and flips it if calculating the 3' background. It then moves the rows containing only zeros to the end of each matrix. The function updates the first position of each sequence in the array to 0, representing the absence of the first nucleotide pattern.

Next, the function calculates the background proportions for each nucleotide pattern by summing the counts in the appropriate positions. The total background counts are calculated, and the proportions are obtained by dividing each count by the total.

The function returns a dictionary with the nucleotide patterns as keys and their background proportions as values.

----

join_batches
---------------


Get and join the data returned from multiprocessed batches. This function takes a list of dataframes containing read information and a dictionary containing sequence data returned from multiprocessed batches. It returns a tuple containing the joined read dataframe, the total counts of nucleotide patterns per nucleotide position, and the background frequency of nucleotide patterns for the 5' and 3' ends.

Parameters:
    - ``read_batches`` (list): List of dataframes containing read information returned from multiprocessed batches.
    - ``full_sequence_batches`` (dict): Dictionary containing sequence data (counts per position and background) returned from multiprocessed batches.

Returns:
    - ``tuple``: A tuple containing the joined read dataframe, the sequence data dictionary, and the sequence background dictionary.

The function starts by calling the `get_batch_data` helper function to retrieve the data from the asynchronous objects. It then proceeds to join the read data by concatenating the dataframes in the `read_batches` list. The categorical columns in the joined dataframe are converted to the "category" dtype.

Next, the function joins the sequence data by iterating over the pattern lengths and patterns in the `sequence_batches` dictionary. It determines the maximum length among the arrays and pads them with zeros to match the maximum length. The padded arrays are summed along the 0-axis to obtain the total counts of nucleotide patterns per position.

Similarly, the function joins the sequence backgrounds by iterating over the pattern lengths and backgrounds in the `background_batches` dictionary. For each pattern, it calculates the weighted sum of the background proportions using the counts from the "sequence_number" key in the dictionary. It then calculates the weighted average for each key and stores the results in the `sequence_background` dictionary.

Finally, the function returns a tuple containing the joined read dataframe, the sequence data dictionary, and the sequence background dictionary.

----

get_batch_data
~~~~~~~~~~~~~~~~~~

Return readable data from the multiprocessed pools, separating the full sequence data into background data and sequence data. This function is called in the `join_batches` function.

Parameters:
    - ``read_batches`` (list): List of dataframes containing read information returned from multiprocessed batches.
    - ``full_sequence_batches`` (dict): Dictionary containing sequence data (counts per position and background) returned from multiprocessed batches.

Returns:
    - ``tuple``: A tuple containing the read batches, background batches, and sequence batches.

The function iterates over the pattern lengths and results in the `full_sequence_batches` dictionary. For each result, it checks if the pattern contains "bg" or "sequence" to determine whether it belongs to the background or sequence data. It appends the array to the corresponding dictionary entry based on the pattern and pattern length.

Finally, the function returns a tuple containing the read batches, background batches, and sequence batches.


Read Data Frame Modifications
=============

a_site_calculation
---------------

Adds a column to the ``read_df`` containing the A-site for the reads.

Parameters:
    - ``read_df`` (pd.DataFrame): Dataframe containing the read information.
    - ``offset`` (int, optional): Offset from the start of the read to the A-site (default: 15).

Returns:
    - ``asite_df`` (pd.DataFrame): Dataframe containing the read information with an added column for the A-site.
    
----

annotate_reads
---------------

Merges the annotation dataframe with the read dataframe. This function takes a dataframe containing read information with an added column for the A-site location (`a_site_df`) and a dataframe containing CDS start/stop and transcript ID from a GFF file (`annotation_df`). It returns a dataframe containing the read information with an added column for the A-site location and columns from the GFF file.

Parameters:
    - ``a_site_df`` (pd.DataFrame): Dataframe containing the read information with an added column for the A-site location.
    - ``annotation_df`` (pd.DataFrame): Dataframe containing the CDS start/stop and transcript ID from a GFF file.

Returns:
    - ``pd.DataFrame``: A dataframe containing the read information with an added column for the A-site location and columns from the GFF file.

The function assigns the transcript ID to a new column in `a_site_df` by splitting the `reference_name` column on the "|" character and taking the first element. It then merges `a_site_df` with `annotation_df` based on the transcript ID, using the `merge` method.

Finally, the function converts the "transcript_id" column to the "category" dtype and drops the "reference_name" column from the resulting dataframe.

----

assign_mRNA_category
---------------

Adds the mRNA category column to the annotated read dataframe, labeling the read according to the position of the A-site. This function takes the annotated read dataframe (`annotated_read_df`) and returns a string with the category for each read.

Parameters:
    - ``annotated_read_df``: Dataframe with read data, added A-site positions, and joined with `annotation_df`.

Returns:
    - ``str``: The mRNA category for the read (one of ["five_leader", "start_codon", "CDS", "stop_codon", "three_trailer"]).

The function calculates the mRNA category based on the A-site position and CDS start/stop values using conditional statements. The conditions compare the A-site position with the CDS start and end positions, and the choices represent the corresponding mRNA categories.

The `np.select` function is used to apply the conditions and choices to each row of the dataframe, assigning the appropriate mRNA category to the "mRNA_category" column.

Finally, the function converts the "mRNA_category" column to the "category" dtype and returns the annotated read dataframe.


