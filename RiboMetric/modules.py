"""
This script contains the functions required to run individual modules
of the RibosomeProfiler pipeline

"""

import pandas as pd
import numpy as np
from xhtml2pdf import pisa
from collections import Counter

from typing import List, Dict


def read_df_to_cds_read_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the a_site_df to a cds_read_df by removing reads that do not
    map to the CDS

    Inputs:
        df: Dataframe containing the read information and annotation

    Outputs:
        cds_read_df: Dataframe containing the read information for reads
                    that map to the CDS
    """
    cds_read_df = df[
        (df["cds_start"] < df["a_site"]) & (df["a_site"] < df["cds_end"])
    ]
    return cds_read_df


def a_site_calculation(read_df: pd.DataFrame, offset=12) -> pd.DataFrame:
    """
    Adds a column to the read_df containing the A-site for the reads

    Inputs:
        read_df: Dataframe containing the read information
        offset: Offset from the start of the read to the A-site (Default = 15)

    Outputs:
        asite_df: Dataframe containing the read information with an added
                    column for the A-site
    """
    a_site_df = read_df.assign(a_site=read_df.reference_start.add(offset))
    return a_site_df


def a_site_calculation_variable_offset(
        read_df: pd.DataFrame,
        offset_dict: dict = {},
        ) -> pd.DataFrame:
    """
    Adds a column to the read_df containing the A-site for the reads

    Inputs:
        read_df: Dataframe containing the read information
        offset_dict: Dictionary containing offsets for each read length
                     Keys: read_length, Values: offset
                     If offset_dict is None, a default offset of 15 is
                     used for all read lengths.

    Outputs:
        asite_df: Dataframe containing the read information with an added
                    column for the A-site
    """
    # If offset_dict is not provided, use default offset of
    # 15 for all read lengths
    if offset_dict is None:
        offset = 15
    else:
        # Create a mapping from read_length to offset
        offset_mapping = {
            length: offset_dict.get(length, 15)
            for length in read_df['read_length'].unique()
            }
        # Map offsets to corresponding read lengths
        read_df['offset'] = read_df['read_length'].map(offset_mapping)
        read_df['offset'] = read_df['offset'].astype('int64')
        offset = read_df['offset']

    read_df['reference_start'] = read_df['reference_start'].astype(int)

    # Calculate A-site based on offset for each read
    a_site_df = read_df.assign(a_site=read_df['reference_start'] + offset)
    a_site_df.to_csv('a_site_df_offsets.csv')
    return a_site_df


def read_length_distribution(read_df: pd.DataFrame) -> dict:
    """
    Calculate the read length distribution for the full dataset

    Inputs:
        read_df: Dataframe containing the read information

    Outputs:
        dict: Dictionary containing the read length distribution
    """
    read_lengths, read_counts = np.unique(
        read_df["read_length"], return_counts=True
    )
    return dict(zip(read_lengths.tolist(), read_counts.tolist()))


def terminal_nucleotide_bias_distribution(
    read_df: pd.DataFrame,
    pattern_length: int = 2,
    keep_N: bool = False,
    target: str = "both",
) -> dict:
    """
    Calculate the proportion of the occurrence in the first or last n
    nucleotides of the reads to check for ligation bias

    Inputs:
        read_df: Dataframe containing read information
        # pattern_length: Length of nucleotide pattern
        keep_N: Keep nucleotide patterns with 'N', or discard if False
        target: Calculate ligation bias for 5', 3' or both

    Outputs:
        terminal_nucleotide_bias_dict: Dictionary containing the distribution of the
        first pattern of nucleotides in the reads
    """
    terminal_nucleotide_bias_dict: dict = (
        {target: {}} if target != "both" else {
            "five_prime": {}, "three_prime": {}
            }
                          )

    total_counts = len(read_df)
    prime_counts = {
        "five_prime": read_df["first_dinucleotide"].value_counts(),
        "three_prime": read_df["last_dinucleotide"].value_counts(),
    }

    categories = {
        "five_prime": read_df["first_dinucleotide"].cat.categories.to_list(),
        "three_prime": read_df["last_dinucleotide"].cat.categories.to_list(),
    }

    pattern_list = read_df["first_dinucleotide"].cat.categories.to_list()
    pattern_list += read_df["last_dinucleotide"].cat.categories.to_list()
    pattern_list = sorted(
        list(set(categories["five_prime"]) | set(categories["three_prime"]))
        )

    if keep_N:
        pattern_list = sorted(pattern_list, key=lambda x: ('N' in x, x))
    else:
        pattern_list = [
            pattern for pattern in pattern_list if 'N' not in pattern
            ]

    for pattern in pattern_list:
        for prime in terminal_nucleotide_bias_dict:
            if pattern in categories[prime]:
                terminal_nucleotide_bias_dict[prime][pattern] = \
                    prime_counts[prime][pattern]/total_counts

    return terminal_nucleotide_bias_dict


def normalise_ligation_bias(
    terminal_nucleotide_bias_dict: dict,
    sequence_background: dict,
    pattern_length: int = 2,
) -> dict:
    """
    Calculate the difference between the observed and expected nucleotide
    pattern at the start and end of the sequences.

    Inputs:
        terminal_nucleotide_bias_dict: Dictionary containing observed proportions for 5'
                            and 3' ends of the sequences
        sequence_background: Dictionary containing expected proportions for 5'
                            and 3' directions of sequences
        # pattern_length: Length of nucleotide pattern

    Outputs:
        terminal_nucleotide_bias_dict_norm: Modified terminal_nucleotide_bias_dict to show the
                                difference between observed and expected
                                distributions
    """
    terminal_nucleotide_bias_dict_norm = terminal_nucleotide_bias_dict
    expected_distribution = {
        "five_prime": sequence_background["5_prime_bg"],
        "three_prime": sequence_background["3_prime_bg"],
        }

    for prime in terminal_nucleotide_bias_dict_norm:
        for pattern in terminal_nucleotide_bias_dict_norm[prime]:

            if pattern in expected_distribution[prime]:
                terminal_nucleotide_bias_dict_norm[prime][pattern]\
                    -= expected_distribution[prime][pattern]

    return terminal_nucleotide_bias_dict_norm


def slicer_vectorized(array: np.ndarray, start: int, end: int):
    """
    String slicer for numpy arrays

    Note: https://stackoverflow.com/a/39045337

    Inputs:
        array: A numpy array of strings
        start: The start position of the slice
        end: The end position of the slice

    Outputs:
        sliced_array: An array consisting of only the selected characters
        from the input string array
    """
    sliced_array = array.view(str).reshape(len(array), -1)[:, start:end]
    return np.frombuffer(sliced_array.tobytes(), dtype=(str, end-start))


def nucleotide_composition(sequence_data_single: dict) -> dict:
    """
    Calculate the proportions of nucleotides for each read position

    Inputs:
        sequence_data_single: A dictionary containing the counts for single
        nucleotides on each read position

    Outputs:
        nucleotide_composition_dict: A dictionary containing the proportion
        for single nucleotides on each read position
    """
    read_length = len(sequence_data_single["A"])
    nucleotide_composition_dict: dict = {nt: [] for nt in ["A", "C", "G", "T"]}
    for position in range(read_length):
        position_count = 0
        for nt in sequence_data_single:
            position_count += sequence_data_single[nt][position]
        for nt in sequence_data_single:
            nucleotide_composition_dict[nt].append(
                sequence_data_single[nt][position] / position_count
                )

    return nucleotide_composition_dict


def read_frame_cull(read_frame_dict: dict, config: dict) -> dict:
    """
    Culls the read_frame_dict according to config so only read lengths of
    interest are kept

    Inputs:
    read_frame_dict:
    config:

    Outputs:
    culled_read_frame_dict
    """
    culled_read_frame_dict = read_frame_dict.copy()
    cull_list = list(culled_read_frame_dict.keys())
    for k in cull_list:
        if (
            k > config["plots"]["read_frame_distribution"]["upper_limit"]
            or k < config["plots"]["read_frame_distribution"]["lower_limit"]
        ):
            del culled_read_frame_dict[k]

    return culled_read_frame_dict


def read_frame_score_trips_viz(read_frame_dict: dict) -> dict:
    """
    Generates scores for each read_length separately as well as a global score
    Can be used after read_frame_cull to calculate the global score of the
    region of interest. The calculation for this score is: 1 - sum(2nd highest
    peak count)/sum(highest peak count). A score close to 1 has good
    periodicity, while a score closer to 0 has a random spread

    Inputs:
    read_frame_dict: dictionary containing the distribution of the reading
                    frames over the different read lengths

    Outputs:
    scored_read_frame_dict: dictionary containing read frame distribution
                            scores for each read length and a global score
    """
    scored_read_frame_dict: Dict[str, float] = {}
    highest_peak_sum, second_peak_sum = 0, 0
    for k, inner_dict in read_frame_dict.items():
        top_two_values = sorted(inner_dict.values(), reverse=True)[:2]
        if top_two_values[0] == 0:
            scored_read_frame_dict[k] = 0
            continue
        elif top_two_values[1] == 0:
            scored_read_frame_dict[k] = 1
        else:
            highest_peak_sum += top_two_values[0]
            second_peak_sum += top_two_values[1]
            scored_read_frame_dict[k] = 1 -\
                top_two_values[1] / top_two_values[0]

    if highest_peak_sum == 0:
        scored_read_frame_dict["global"] = 0.0
    else:
        scored_read_frame_dict["global"] = 1 -\
            second_peak_sum / highest_peak_sum
    return scored_read_frame_dict


def read_frame_distribution(a_site_df: pd.DataFrame) -> dict:
    """
    Calculate the distribution of the reading frame over the dataset

    Inputs:
    a_site_df: Dataframe containing the read information with an added column for the a-site location

    Outputs:
    read_frame_dict: Nested dictionary containing counts for every reading frame at the different read lengths
    """
    read_frame_dict = {}

    # Iterate over unique combinations of transcript_id and read_length
    for (transcript_id, read_length), group in a_site_df.groupby(['reference_name', 'read_length']):
        group['read_frame'] = group['a_site'] % 3
        frame_counts = group['read_frame'].value_counts().sort_values(ascending=False).to_dict()

        # Assign frame numbers based on sorted order per transcript
        # The most frequent frame is assigned 0, the second most frequent is assigned 1, and the least frequent is assigned 2
        frame_count_dict = {idx: count for idx, count in enumerate(sorted(frame_counts.values(), reverse=True))}

        if read_length not in read_frame_dict:
            read_frame_dict[int(read_length)] = {0: 0, 1: 0, 2: 0}

        for frame, count in frame_count_dict.items():
            read_frame_dict[int(read_length)][frame] += count

    return read_frame_dict


def read_frame_distribution_annotated(
        annotated_read_df: pd.DataFrame,
        exclusion_length: int = 0,
        read_length_range: tuple = (20, 40)
        ) -> dict:
    """
    Calculate the distribution of the reading frame over the dataset

    Inputs:
        a_site_df: Dataframe containing the read information with an added
        column for the a-site location

    Outputs:
        read_frame_dict: Nested dictionary containing counts for every reading
        frame at the different read lengths
    """
    read_lengths: List[int] = [
        i for i in range(read_length_range[0], read_length_range[1])
        ]

    df_slice = annotated_read_df[annotated_read_df["cds_start"] != 0]
    df_slice = df_slice[
        (df_slice["a_site"] > df_slice["cds_start"] + exclusion_length) &
        (df_slice["a_site"] < df_slice["cds_end"] - exclusion_length)
    ]
    frame_df = (
        df_slice.assign(read_frame=(df_slice.a_site-df_slice.cds_start).mod(3))
        .groupby(["read_length", "read_frame"])
        .size()
    )
    read_frame_dict: Dict[int, Dict[int, int]] = {}
    for index, value in frame_df.items():
        read_length: int
        read_frame: int
        read_length, read_frame = index
        if read_length in read_lengths:
            if read_length not in read_frame_dict:
                read_frame_dict[read_length] = {0: 0, 1: 0, 2: 0}
            read_frame_dict[read_length][read_frame] = value
    return read_frame_dict


def annotate_reads(
    a_site_df: pd.DataFrame, annotation_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merges the annotation dataframe with the read dataframe

    Inputs:
        a_site_df: Dataframe containing the read information with an added
        column for the a-site location
        annotation_df: Dataframe containing the CDS start/stop
        and transcript id from a gff file.

    Outputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along
        with the columns from the gff file
    """
    annotated_read_df = a_site_df.assign(
        transcript_id=a_site_df.reference_name.str.split("|").str[0]
    ).merge(annotation_df, on="transcript_id")
    annotated_read_df["transcript_id"] = (annotated_read_df["transcript_id"]
                                          .astype("category"))
    return annotated_read_df.drop(["reference_name"], axis=1)


def chunked_annotate_reads(a_site_df: pd.DataFrame,
                           annotation_df: pd.DataFrame,
                           chunk_size: int = 10000000) -> pd.DataFrame:
    """
    Merges the annotation dataframe with the read dataframe in smaller chunks.

    Inputs:
        a_site_df: DataFrame containing the read information with an added
        column for the a-site location.
        annotation_df: DataFrame containing the CDS start/stop
        and transcript id from a gff file.
        chunk_size: Size of each processing chunk.

    Outputs:
        annotated_read_df: DataFrame containing the read information
        with an added column for the a-site location along
        with the columns from the gff file.
    """
    # Initialize an empty list to store processed chunks
    processed_chunks = []

    # Split a_site_df into chunks
    num_chunks = len(a_site_df) // chunk_size + 1
    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(a_site_df))
        chunk = a_site_df.iloc[start_idx:end_idx]

        # Process the chunk
        chunk = chunk.assign(
            transcript_id=chunk.reference_name.str.split("|").str[0]
        )

        chunk = chunk.drop(["reference_name"], axis=1)
        chunk = chunk.merge(annotation_df, on="transcript_id")
        chunk["transcript_id"] = chunk["transcript_id"].astype("category")

        # Append the processed chunk to the list
        processed_chunks.append(chunk)

    # Concatenate the processed chunks
    annotated_read_df = pd.concat(processed_chunks)

    return annotated_read_df


def assign_mRNA_category(annotated_read_df) -> pd.DataFrame:

    """
    Adds the mRNA category column to the annotated_read_df, labelling the read
    according to the position of the A-site
    Assign an mRNA category based on the A-site of the read
    and the CDS start/stop, used through df.apply()

    Inputs:
        annotated_read_df: Dataframe with read data, added a-site positions
        and joined with annotation_df.

    Outputs:
        mRNA category: string with the category for the read
        ["five_leader", "start_codon", "CDS", "stop_codon", "three_trailer"]
    """
    # Calculate mRNA category based on conditions
    conditions = [
        annotated_read_df["a_site"] < annotated_read_df["cds_start"],
        annotated_read_df["a_site"] == annotated_read_df["cds_start"],
        (annotated_read_df["cds_start"] < annotated_read_df["a_site"]) &
        (annotated_read_df["a_site"] < annotated_read_df["cds_end"]),
        annotated_read_df["a_site"] == annotated_read_df["cds_end"],
        annotated_read_df["a_site"] > annotated_read_df["cds_end"]
    ]
    choices = [
        "five_leader", "start_codon", "CDS", "stop_codon", "three_trailer"
    ]
    annotated_read_df["mRNA_category"] = np.select(
        conditions,
        choices,
        "unknown"
        )
    annotated_read_df["mRNA_category"] = \
        annotated_read_df["mRNA_category"].astype("category")
    return annotated_read_df


def mRNA_distribution(annotated_read_df: pd.DataFrame) -> dict:
    """
    Calculate the distribution of the mRNA categories over the read length

    Inputs:
        annotated_read_df: Dataframe containing the read information
                           with an added column for the a-site location along
                           with the columns from the gff file
    Outputs:
        mRNA_distribution_dict: Nested dictionary containing counts for every
                                mRNA category at the different read lengths
    """
    # Creating MultiIndex for reindexing
    categories = [
        "five_leader",
        "start_codon",
        "CDS",
        "stop_codon",
        "three_trailer",
    ]
    classes = annotated_read_df["read_length"].unique()
    idx = pd.MultiIndex.from_product(
        [classes, categories], names=["class", "category"]
    )
    # Group annotated_read_df
    annotated_read_df = (
        annotated_read_df.groupby(["read_length", "mRNA_category"])
        .size()
        .reindex(idx, fill_value=0)
        .sort_index()
        .to_frame()  # Convert the resulting series to a DataFrame
        .reset_index()  # Reset the index to a regular column
    )

    # Creating mRNA_distribution_dict from annotated_read_df
    mRNA_distribution_dict: dict = {"global": {}}
    for read_length, mRNA_category, value in annotated_read_df.itertuples(
                                                                index=False,
                                                                name=None):
        if read_length not in mRNA_distribution_dict:
            mRNA_distribution_dict[read_length] = {}
        mRNA_distribution_dict[read_length][mRNA_category] = value
        if mRNA_category in mRNA_distribution_dict["global"]:
            mRNA_distribution_dict["global"][mRNA_category] += value
        else:
            mRNA_distribution_dict["global"][mRNA_category] = value
    return mRNA_distribution_dict


def sum_mRNA_distribution(mRNA_distribution_dict: dict, config: dict) -> dict:
    """
    Calculate the sum of mRNA categories

    Inputs:
        annotated_read_dict: Dataframe containing the read information
        with an added column for the a-site location along
        with the columns from the gff file

    Outputs:
        read_frame_dict: Nested dictionary containing counts for every reading
        frame at the different read lengths
    """
    sum_mRNA_dict: dict = {}
    for inner_dict in mRNA_distribution_dict.values():
        for k, v in inner_dict.items():
            if k in sum_mRNA_dict:
                sum_mRNA_dict[k] += v
            else:
                sum_mRNA_dict[k] = v
    if not config["plots"]["mRNA_distribution"]["absolute_counts"]:
        sum_mRNA_dict = {
            k: (v / sum(sum_mRNA_dict.values()))
            for k, v in sum_mRNA_dict.items()
        }

    return sum_mRNA_dict


def metagene_distance(
    annotated_read_df: pd.DataFrame, target: str = "start"
) -> pd.Series:
    """
    Calculate distance from A-site to start or stop codon

    Inputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along with data from
        the annotation file
        target: Target from which the distance is calculated

    Outputs:
    pd.Series
    """
    if target == "start":
        return annotated_read_df["a_site"] - annotated_read_df["cds_start"]
    elif target == "stop":
        return annotated_read_df["a_site"] - annotated_read_df["cds_end"]
    else:
        raise ValueError("Target must be start or stop")


def metagene_profile(
    annotated_read_df: pd.DataFrame,
    target: str = "both",
    distance_range: list = [-50, 50],
    position: str = "a_site",
) -> dict:
    """
    Groups the reads by read_length and distance to a target and counts them

    Inputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along with data from
        the annotation file
        target: Target from which the distance is calculated
        distance_range: The range of the plot

    Outputs:
        metagene_profile_dict: dictionary containing the read_length of
        the read and distance to the target as keys and the counts as values
    """
    target_loop = [target] if target != "both" else ["start", "stop"]
    metagene_profile_dict: Dict[str, Dict[str, dict]] = {
        "start": {}, "stop": {}
        }
    for current_target in target_loop:
        annotated_read_df = annotated_read_df.assign(
            metagene_info=metagene_distance(annotated_read_df, current_target)
        )
        pre_metaprofile_dict = (
            annotated_read_df[
                (annotated_read_df["metagene_info"] > distance_range[0] - 1)
                & (annotated_read_df["metagene_info"] < distance_range[1] + 1)
            ]
            .groupby(["read_length", "metagene_info"])
            .size()
            .to_dict()
        )

        if pre_metaprofile_dict == {}:  # If no reads in range
            pre_metaprofile_dict = (
                annotated_read_df.groupby(["read_length", "metagene_info"])
                .size()
                .to_dict()
            )
        # Fill empty read lengths with 0
        min_length = int(min([x[0] for x
                              in list(pre_metaprofile_dict.keys())]))
        max_length = int(max([x[0] for x
                              in list(pre_metaprofile_dict.keys())]))
        for y in range(min_length, max_length):
            if y not in [x[0] for x in list(pre_metaprofile_dict.keys())]:
                pre_metaprofile_dict[(y, 0)] = 0

        neg_distance = int(min([x[1] for x
                                in list(pre_metaprofile_dict.keys())]))
        pos_distance = int(max([x[1] for x
                                in list(pre_metaprofile_dict.keys())]))
        position_range = range(neg_distance, pos_distance+1)

        for key, value in pre_metaprofile_dict.items():
            if key[0] not in metagene_profile_dict[current_target]:
                metagene_profile_dict[current_target][key[0]] = {}
            metagene_profile_dict[current_target][key[0]][int(key[1])] = value

        # Fill empty distances with 0
        for position_dict in metagene_profile_dict[current_target].values():
            for position in position_range:
                position_dict.setdefault(int(position), 0)

    return metagene_profile_dict


def proportion_of_kmer(
    annotated_read_df: pd.DataFrame,
) -> list:
    '''
    get proportion of reads with predicted a-sites in each frame for
    sequence of length k

    Inputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along with data from
        the annotation file

    '''
    f1 = annotated_read_df[annotated_read_df['a_site'] % 3 == 0]
    f2 = annotated_read_df[annotated_read_df['a_site'] % 3 == 1]
    f3 = annotated_read_df[annotated_read_df['a_site'] % 3 == 2]

    return [len(f1), len(f2), len(f3)]


def get_cart_point(ternary_point, vertices=[[0.0, 0.0], [1.0, 0.0], [0.5, 1]]):
    '''
    Get the cartesian coordinates of a point in ternary space

    Inputs:
        ternary_point: Point in ternary space
        vertices: Vertices of the triangle

    Outputs:
        cartesian_point: Point in cartesian space
    '''
    point = np.array(ternary_point)
    verts = np.array(vertices)

    return np.dot(point, verts)


def reading_frame_triangle(
        annotated_read_df: pd.DataFrame,
) -> dict:
    '''
    Get the cartesian coordinates of the triangle plot for the reading frame

    Inputs:
        annotated_read_df: Dataframe containing the read information
        with an added column for the a-site location along with data from
        the annotation file

    Outputs:
        triangle_dict: Dictionary containing the cartesian coordinates
        of the triangle plot for the reading frame
    '''
    triangle_dict = {}
    for transcript, df in annotated_read_df.groupby('transcript_id'):
        # Get the proportion of reads with predicted a-sites in each frame
        proportion = proportion_of_kmer(df)
        if len(proportion) < 3:
            continue

        # triangle_dict[transcript] = get_cart_point(proportion[0])
        triangle_dict[transcript] = proportion

    return triangle_dict


def sequence_slice(
    read_df: pd.DataFrame, nt_start: int = 0, nt_count: int = 15
) -> dict:
    sequence_slice_dict = {
        k: v[nt_start: nt_start + nt_count]
        for k, v in read_df["sequence"].to_dict().items()
    }
    return sequence_slice_dict


def convert_html_to_pdf(source_html, output_filename):
    result_file = open(output_filename, "w+b")

    pisa_status = pisa.CreatePDF(source_html, dest=result_file)
    result_file.close()
    return pisa_status.err


# Deprecated
def calculate_expected_dinucleotide_freqs(read_df: pd.DataFrame) -> dict:
    """
    Calculate the expected dinucleotide frequencies based on the
    nucleotide frequencies in the aligned reads

    Inputs:
        read_df: Dataframe containing the read information

    Outputs:
        expected_dinucleotide_freqs: Dictionary containing the expected
        dinucleotide frequencies
    """
    dinucleotides = []
    for read in read_df["sequence"].drop_duplicates():
        for i in range(len(read) - 1):
            dinucleotides.append(read[i: i + 2])

    observed_freq = Counter(dinucleotides)
    total_count = sum(observed_freq.values())

    expected_dinucleotide_freqs = {}
    for dinucleotide, count in observed_freq.items():
        expected_dinucleotide_freqs[dinucleotide] = count / total_count

    return expected_dinucleotide_freqs


def change_point_analysis(
        read_counts: dict,
        surrounding_range: tuple = (-30, 10),
        ) -> int:
    """
    Calculate the change point for the metagene profile
    This should reflect where the cds starts and as a result the
    offset to apply to get a-site

    Inputs:
        read_counts: Dictionary containing the read counts for each position
        surrounding_range: tuple of start and stop for change point analysis

    Outputs:
        max_shift_position: The position of the change point
    """
    max_shift = 0
    max_shift_position = None

    for i in range(surrounding_range[0], surrounding_range[1]):
        left = []
        for i in range(i-3, i+1):
            if i in read_counts:
                left.append(read_counts[i])
            else:
                left.append(0)
        right = []
        for i in range(i+1, i+5):
            if i in read_counts:
                right.append(read_counts[i])
            else:
                right.append(0)

        mean_left = sum(left) / 4
        mean_right = sum(right) / 4
        shift = abs(mean_right - mean_left)
        if shift > max_shift:
            max_shift = shift
            max_shift_position = i
    if not max_shift_position:
        return 15
    return max_shift_position


def asite_calculation_per_readlength(
        annotated_read_df: pd.DataFrame,
        offset_range: tuple = (10, 15),
) -> dict:
    """
    Calculate offset values per read length for the A-site
    Shoelaces based method using metagene counts

    Input:read_counts.get(i, 0) for i in range(i+1, i+5)
        with an cds info added
        offset_range: Range of offsets to test

    Output:
        offset_dict: Dictionary containing the offset values for each read len
    """
    offset_dict: dict = {}

    for read_length in annotated_read_df["read_length"].unique():
        offset_dict[read_length] = {}
        read_length_metagene = metagene_profile(
            annotated_read_df[annotated_read_df["read_length"] == read_length],
            target="start",
            distance_range=[-30, 10],
        )

        offset_dict[read_length] = change_point_analysis(
            read_length_metagene["start"][read_length],
            surrounding_range=(-26, 5)
        )
    return offset_dict
