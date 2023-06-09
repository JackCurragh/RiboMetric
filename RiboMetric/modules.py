"""
This script contains the functions required to run individual modules
of the RibosomeProfiler pipeline

"""

import pandas as pd
import numpy as np
from xhtml2pdf import pisa
from collections import Counter


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


def a_site_calculation(read_df: pd.DataFrame, offset=15) -> pd.DataFrame:
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


def ligation_bias_distribution(
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
        ligation_bias_dict: Dictionary containing the distribution of the
        first pattern of nucleotides in the reads
    """
    ligation_bias_dict = ({target: {}} if target != "both"
                          else {"five_prime": {}, "three_prime": {}})

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
        for prime in ligation_bias_dict:
            if pattern in categories[prime]:
                ligation_bias_dict[prime][pattern] = \
                    prime_counts[prime][pattern]/total_counts

    return ligation_bias_dict


def normalise_ligation_bias(
    ligation_bias_dict: dict,
    sequence_background: dict,
    pattern_length: int = 2,
) -> dict:
    """
    Calculate the difference between the observed and expected nucleotide
    pattern at the start and end of the sequences.

    Inputs:
        ligation_bias_dict: Dictionary containing observed proportions for 5'
                            and 3' ends of the sequences
        sequence_background: Dictionary containing expected proportions for 5'
                            and 3' directions of sequences
        # pattern_length: Length of nucleotide pattern

    Outputs:
        ligation_bias_dict_norm: Modified ligation_bias_dict to show the
                                difference between observed and expected
                                distributions
    """
    ligation_bias_dict_norm = ligation_bias_dict
    expected_distribution = {
        "five_prime": sequence_background["5_prime_bg"],
        "three_prime": sequence_background["3_prime_bg"],
        }

    for prime in ligation_bias_dict_norm:
        for pattern in ligation_bias_dict_norm[prime]:

            if pattern in expected_distribution[prime]:
                ligation_bias_dict_norm[prime][pattern]\
                    -= expected_distribution[prime][pattern]

    return ligation_bias_dict_norm


def slicer_vectorized(array: np.array, start: int, end: int):
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
    sliced_array = array.view((str, 1)).reshape(len(array), -1)[:, start:end]
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
    nucleotide_composition_dict = {nt: [] for nt in ["A", "C", "G", "T"]}
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


def read_frame_score(read_frame_dict: dict) -> dict:
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
    scored_read_frame_dict = {}
    highest_peak_sum, second_peak_sum = 0, 0

    for k, inner_dict in read_frame_dict.items():
        top_two_values = sorted(inner_dict.values(), reverse=True)[:2]
        highest_peak_sum += top_two_values[0]
        second_peak_sum += top_two_values[1]
        scored_read_frame_dict[k] = 1 - top_two_values[1] / top_two_values[0]

    scored_read_frame_dict["global"] = 1 - second_peak_sum / highest_peak_sum
    return scored_read_frame_dict


def read_frame_distribution(a_site_df: pd.DataFrame) -> dict:
    """
    Calculate the distribution of the reading frame over the dataset

    Inputs:
        a_site_df: Dataframe containing the read information with an added
        column for the a-site location

    Outputs:
        read_frame_dict: Nested dictionary containing counts for every reading
        frame at the different read lengths
    """
    frame_df = (
        a_site_df.assign(read_frame=a_site_df.a_site.mod(3))
        .groupby(["read_length", "read_frame"])
        .size()
    )
    read_frame_dict = {}
    for index, value in frame_df.items():
        read_length, read_frame = index
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


def assign_mRNA_category(annotated_read_df) -> str:

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
    )
    # Creating mRNA_distribution_dict from annotated_read_df
    mRNA_distribution_dict = {}
    for index, value in annotated_read_df.items():
        read_length, mRNA_category = index
        if read_length not in mRNA_distribution_dict:
            mRNA_distribution_dict[read_length] = {}
        mRNA_distribution_dict[read_length][mRNA_category] = value

    # Setting order of categories 5' to 3'
    for i in mRNA_distribution_dict:
        mRNA_distribution_dict[i] = {
            k: mRNA_distribution_dict[i][k]
            for k in categories
            if k in mRNA_distribution_dict[i]
        }

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
    sum_mRNA_dict = {}
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


def metagene_profile(
    annotated_read_df: pd.DataFrame,
    target: str = "both",
    distance_range: list = [-50, 50],
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
    metagene_profile_dict = {"start": {}, "stop": {}}
    for current_target in target_loop:
        annotated_read_df["metagene_info"] = metagene_distance(
            annotated_read_df, current_target
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
        if pre_metaprofile_dict == {}:
            print(
                "ERR - Metagene Profile: No reads found in specified range, \
removing boundaries..."
            )
            pre_metaprofile_dict = (
                annotated_read_df.groupby(["read_length", "metagene_info"])
                .size()
                .to_dict()
            )
        # Fill empty read length and distance keys with None
        min_length = min([x[0] for x in list(pre_metaprofile_dict.keys())])
        max_length = max([x[0] for x in list(pre_metaprofile_dict.keys())])
        for y in range(min_length, max_length):
            if y not in [x[0] for x in list(pre_metaprofile_dict.keys())]:
                pre_metaprofile_dict[(y, 0)] = None
        min_distance = min([x[1] for x in list(pre_metaprofile_dict.keys())])
        max_distance = max([x[1] for x in list(pre_metaprofile_dict.keys())])
        for z in range(min_distance, max_distance):
            if z not in [x[1] for x in list(pre_metaprofile_dict.keys())]:
                pre_metaprofile_dict[(min_length, z)] = None
        for key, value in pre_metaprofile_dict.items():
            if key[0] not in metagene_profile_dict[current_target]:
                metagene_profile_dict[current_target][key[0]] = {}
            metagene_profile_dict[current_target][key[0]][key[1]] = value
    return metagene_profile_dict


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
def calculate_expected_dinucleotide_freqs(read_df: pd.DataFrame) -> dict():
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
