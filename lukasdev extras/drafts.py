# Human Readable dictionary, not useful for plots
def nucleotide_composition(read_df: pd.DataFrame) -> dict:
    """
    Calculate the nucleotide composition

    Inputs:
        read_df: Dataframe containing the read information

    Outputs:
        dict: Dictionary containing the nucleotide distribution for every read position.
    """
    readlen = read_df.sequence.str.len().max()
    nucleotide_composition_dict = dict()
    for i in range(readlen):
        nucleotide_counts = read_df.sequence.str.slice(i,i+1).value_counts()
        nucleotide_counts.drop("", errors='ignore', inplace=True)
        nucleotide_sum = nucleotide_counts.sum()
        nucleotide_composition_dict[i] = {"A" : nucleotide_counts["A"]/nucleotide_sum,"C" : nucleotide_counts["C"]/nucleotide_sum,"G" : nucleotide_counts["G"]/nucleotide_sum,"T" : nucleotide_counts["T"]/nucleotide_sum}
    return nucleotide_composition_dict

# Original plotting function heatmap (using tuple keys instead of nested dictionary)
def plot_metagene_heatmap(metagene_heatmap_dict: dict,config: dict) -> dict:
    """
    Generate a heatmap of the reads depending on their distance
    to a target, read length and count

    Inputs:
        metagene_heatmap_dict: Dictionary containing the counts as values and distance
        from target as keys
        config: Dictionary containing the configuration information

    Outputs:
        plot_metagene_heatmap: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    fig = go.Figure(
        data=go.Heatmap(
            x = [x[1] for x in list(pre_heatmap_dict.keys())],
            y = [x[0] for x in list(pre_heatmap_dict.keys())],
            z = list(pre_heatmap_dict.values()),
        colorscale='electric',
        zmin = 0,
        zmax = config["plots"]["metagene_heatmap"]["max_colorscale"]
    ))
    fig.update_xaxes(range=config["plots"]["metagene_heatmap"]["distance_range"])
    fig.update_layout(
            title="Metagene Heatmap",
            xaxis_title="Read length",
            yaxis_title="Relative position",
            font=dict(
                family=config["plots"]["font_family"],
                size=18,
                color=config["plots"]["base_color"],
            ),
            legend={'traceorder':'normal'},
    )
    plot_metagene_heatmap = {
        "name": "Metagene Heatmap",
        "description": "Metagene heatmap showing the distance between the A-site and \
a target per read length and the counts in colorscale.",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig, config)
    }
    return plot_metagene_heatmap

# multiprocess.process
def parse_bam(bam_file, batch_size, num_processes, num_reads=1000000):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    # pool = mp.Pool(processes=num_processes)
    results = []
    read_list = []
    output_queue = mp.Queue()
    processes = []

    for i, read in enumerate(samfile.fetch()):
        if num_reads and i > num_reads:
            break
        read_list.append(read.to_string().split(sep="\t"))

        if len(read_list) == batch_size:
            # Use the output_queue to store the results asynchronously
            process = mp.Process(target=process_reads, args=(read_list, output_queue))
            process.start()
            processes.append(process)
            read_list = []

        read_percentage = round((i + 1) / num_reads * 100, 3)
        print(f"Processed {i + 1}/{num_reads} ({read_percentage}%)", end="\r")

    if read_list:
        process = mp.Process(target=process_reads, args=(read_list, output_queue))
        process.start()
        processes.append(process)

    print(processes)
    print("finished parsing")

    # Retrieve the results from the output_queue
    while not output_queue.empty():
        results.append(output_queue.get())
    print("appended results from output_queue")

    # # Wait for all processes to complete
    # for process in processes:
    #     process.join()
    # print("joined processes")

    # Yield the collected batch DataFrames
    for batch_df in results:
        yield batch_df
    print("yielded results... parse_bam complete")

# Memory intensive
def cds_coverage_metric(
        cds_read_df: pd.DataFrame,
        minimum_reads: int = 1,
        in_frame_coverage: bool = True
        ) -> float:
    """
    Calculates the proportion of CDS covered by ribosomal protected fragments

    Inputs:
        annotated_read_df: Dataframe containing the reads that have
        a transcript available in the provided annotation
        minimum_reads: The minimum amount of reads that should cover
        a specific nucleotide to be counted for the proportion
        in_frame_count: If set to True, only controls the coverage in frame

    Outputs:
        cds_coverage: A proportion of the amount of individual nucleotides
        represented by the A-sites over the total number of nucleotides in
        the CDS of transcripts present in the reads
    """
    # Calculate the total combined length of the CDS of transcripts that have
    # reads aligned to them
    cds_transcripts = cds_read_df[~cds_read_df["transcript_id"]
                                  .duplicated()
                                  ][["transcript_id", "cds_start", "cds_end"]]
    cds_transcripts["cds_length"] = (cds_transcripts
                                     .apply(lambda x: x['cds_end']
                                            - x['cds_start'],
                                            axis=1))
    cds_length_total = cds_transcripts["cds_length"].sum()
    del cds_transcripts

    # Subset the dataframe so that only reads where the count covering a
    # position is greater than minimum_reads
    cds_reads_count = cds_read_df[cds_read_df
                                  .groupby(["transcript_id", "a_site"])
                                  .transform('size') > minimum_reads
                                  ][["transcript_id", "a_site", "cds_start"]]
    
    # If in_frame_coverage is true, take only reads that are in frame for
    # their transcript and divide the combined CDS length by 3
    if in_frame_coverage:
        cds_reads_count[
        (cds_reads_count["a_site"] - cds_reads_count["cds_start"])%3 == 0
        ]
        cds_length_total = cds_length_total/3
    
    # Calculate the count of nucleotides covered by the reads after filtering
    cds_reads_count = (cds_reads_count
                     .groupby(["transcript_id", "a_site"])
                     .size()
                     .groupby('transcript_id')
                     )
    cds_reads_count = cds_reads_count.size().sum()

    return cds_reads_count/cds_length_total

# deprecated
def nucleotide_composition(
    read_df: pd.DataFrame, nucleotides=["A", "C", "G", "T"]
) -> dict:
    """
    Calculate the nucleotide composition, the proportion of nucleotides
    in each nucleotide position

    Inputs:
        read_df: Dataframe containing the read information
        nucleotides: list of counted nucleotides

    Outputs:
        dict: Dictionary containing the nucleotide distribution for every
            read position.
    """
    sequence_array = read_df["sequence"].to_numpy().astype(str)
    readlen = len(max(sequence_array, key=len))
    nucleotide_composition_dict = {nt: [] for nt in nucleotides}
    for i in range(readlen):
        nucleotide_array = slicer_vectorized(sequence_array, i, i+1)
        nucleotide_array = nucleotide_array[nucleotide_array != '']

        nucleotide, counts = np.unique(nucleotide_array, return_counts=True)

        nucleotide_counts = dict(zip(nucleotide, counts))
        nucleotide_sum = sum(nucleotide_counts.values())

        for nt in nucleotide_composition_dict.keys():
            if nt in nucleotide_counts.keys():
                nucleotide_composition_dict[nt].append(
                    nucleotide_counts[nt]/nucleotide_sum
                )
            else:
                nucleotide_composition_dict[nt].append(
                    0
                )
    return nucleotide_composition_dict

# deprecated
def global_nucleotide_proportion(
    read_df: pd.DataFrame, num_bases: int = 2, five_prime: bool = True
) -> dict:
    """
    Calculate the global proportion nucleotide groups in the reads,
    used as a background for ligation bias distribution

    Inputs:
        read_df: Dataframe containing the read information
        num_bases: Number of bases to be read (Default = 2)
        five_prime: Start at 5' end (True) or 3' end (False) of read
        (Default = True)

    Outputs:
        dinucleotide_counts: Dictionary containing the distribution of the
        nucleotide groups in the reads
    """
    # Remove first n characters and last character from sequences
    # with odd lengths
    if five_prime:
        series = (
            read_df["sequence"]
            .drop_duplicates()
            .apply(lambda x: x[num_bases: -(len(x) % num_bases)])
        )
    # If five_prime is false, remove last n characters and first if odd length
    else:
        series = (
            read_df["sequence"]
            .drop_duplicates()
            .apply(lambda x: x[len(x) % num_bases: -num_bases])
        )
    # Concatenate all strings in the Series
    concatenated = "".join(series.tolist())
    # Calculate dinucleotide occurrences
    expected_nucleotide_proportion = Counter(
        concatenated[i : i + num_bases]
        for i in range(0, len(concatenated) - 1, num_bases)
    )
    sumcounts = sum(expected_nucleotide_proportion.values())
    for k in expected_nucleotide_proportion:
        expected_nucleotide_proportion[k] = (
            expected_nucleotide_proportion[k] / sumcounts
        )
    return dict(expected_nucleotide_proportion)

# deprecated
def ligation_bias_distribution(
    read_df: pd.DataFrame,
    num_bases: int = 2,
    five_prime: bool = True,
) -> dict:
    """
    Calculate the proportion of the occurrence in the first or last n
    nucleotides of the reads to check for ligation bias

    Inputs:
        read_df: Dataframe containing the read information
        num_bases: Number of bases to be read (Default = 2)
        five_prime: Start at 5' end (True) or 3' end (False) of read
        (Default = True)

    Outputs:
        ligation_bias_dict: Dictionary containing the distribution of the
        first two nucleotides in the reads
    """
    if five_prime:
        sequence_dict = dict(
            read_df["sequence"]
            .str.slice(stop=num_bases)
            .value_counts(normalize=True)
            .sort_index()
        )
    else:
        sequence_dict = dict(
            read_df["sequence"]
            .str.slice(start=-num_bases)
            .value_counts(normalize=True)
            .sort_index()
        )
    ligation_bias_dict = {
        k: v for k, v in sequence_dict.items() if "N" not in k
    }
    ligation_bias_dict.update(
        {k: v for k, v in sequence_dict.items() if "N" in k}
    )
    return ligation_bias_dict
