"""
This script contains the functions used to load the required files
in the RiboMetric pipeline

The functions are called by the main script RiboMetric.py
"""
from Bio import SeqIO
from pysam import AlignmentFile
import gffpandas.gffpandas as gffpd

import pandas as pd
import numpy as np
import tempfile
import gzip
import os

from multiprocessing import Pool
from tempfile import TemporaryDirectory


from .bam_processing import join_batches, ox_parse_reads, ox_server_parse_reads
from .file_splitting import split_gff_file
from .file_splitting import run_samtools_idxstats, split_idxstats_df


def parse_annotation(annotation_path: str) -> pd.DataFrame:
    """
    Read in the annotation file at the provided path and return a dataframe

    Inputs:
        annotation_path: Path to the annotation file with columns:
                                    "transcript_id","cds_start",
                                    "cds_end","transcript_length",
                                    "genomic_cds_starts","genomic_cds_ends"

    Outputs:
        annotation_df: Dataframe containing the annotation information
    """
    return pd.read_csv(
        annotation_path,
        sep="\t",
        dtype={
            "transcript_id": str,
            "cds_start": int,
            "cds_end": int,
            "transcript_length": int,
            "genomic_cds_starts": str,
            "genomic_cds_ends": str,
        },
    )


def parse_fasta(fasta_path: str) -> dict:
    """
    Read in the transcriptome fasta file at the provided path
    and return a dictionary

    Inputs:
        fasta_path: Path to the transcriptome fasta file

    Outputs:
        transcript_dict: Dictionary containing the
                         transcript information
    """
    transcript_dict = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        transcript_dict[record.id] = record

    return transcript_dict


def check_bam(bam_path: str) -> bool:
    """
    Check whether the bam file and its index exists at the provided path
    Return True if both files exist, False otherwise

    Inputs:
        bam_path: Path to the bam file

    Outputs:
        bool: True if the bam file and its index exist, False otherwise
    """
    if os.path.exists(bam_path) and os.path.exists(bam_path + ".bai"):
        return True
    else:
        return False


def flagstat_bam(bam_path: str) -> dict:
    """
    Run samtools flagstat on the bam file at the provided path
    and return a dictionary

    Inputs:
        bam_path: Path to the bam file

    Outputs:
        flagstat_dict: Dictionary containing the flagstat information

    """
    flagstat_dict = {}
    with AlignmentFile(bam_path, "rb") as bamfile:
        flagstat_dict["total_reads"] = bamfile.mapped + bamfile.unmapped
        flagstat_dict["mapped_reads"] = bamfile.mapped
        flagstat_dict["unmapped_reads"] = bamfile.unmapped
        flagstat_dict["duplicates"] = bamfile.mapped + bamfile.unmapped
    return flagstat_dict


def parse_bam(bam_file: str,
              num_reads: int,
              batch_size: int = 10000000,
              num_processes: int = 4,
              server_mode: bool = False,
              ) -> tuple:
    """
    Read in the bam file at the provided path and return parsed read and
    sequence data

    Inputs:
        bam_file: Path to the bam file
        num_reads: Maximum number of reads to parse
        batch_size: The number of reads that are processed at a time
        num_processes: The maximum number of processes that this function can
                       create

    Outputs:
        parsed_bam: Tuple containing:
            read_df_pre: The read dataframe containing read information before
                         further modifications to the dataframe
            sequence_data: Dictionary containing the total counts of
                           nucleotide patterns per nucleotide position
            sequence_background: Dictionary containing the background
                                frequency of nucleotide patterns for five and
                                three prime
    """
    batch_size = int((num_reads/num_processes)*1.02)
    # Small percentage increase to ensure remaining reads aren't
    # in separate batch

    print(f"Splitting BAM into {batch_size} reads")
    if server_mode is False:
        pool = Pool(processes=num_processes)
        bam_batches = []
        with TemporaryDirectory() as tempdir:
            idxstats_df = run_samtools_idxstats(bam_file)
            reference_dfs = split_idxstats_df(idxstats_df,
                                              batch_size,
                                              num_reads)
            for split_num, reference_df in enumerate(reference_dfs):
                bam_batches.append(pool.apply_async(ox_parse_reads,
                                                    [bam_file,
                                                     split_num,
                                                     reference_df,
                                                     tempdir]))

            pool.close()
            pool.join()

            parsed_bam = join_batches(bam_batches)

    else:
        print("Warning: The server option is not working as intended.",
              "Regular runs are recommended.")
        bam_batches = ox_server_parse_reads(bam_file)
        parsed_bam = join_batches(bam_batches)

    return parsed_bam


def get_top_transcripts(read_df: dict, num_transcripts: int) -> list:
    """
    Get the top N transcripts with the most reads

    Inputs:
        read_df: DataFrame containing the read information
        num_transcripts: Number of transcripts to return

    Outputs:
        top_transcripts: List of the top N transcripts
    """
    count_sorted_df = (
        read_df.groupby("reference_name")
               .sum()
               .sort_values("count", ascending=False)
    )

    return count_sorted_df.index[:num_transcripts].tolist()


def subset_gff(gff_df: pd.DataFrame) -> pd.DataFrame:
    """
    Subset the gff dataframe to only include the CDS features

    Inputs:
        gff_df: Dataframe containing the gff information

    Outputs:
        gff_df: Dataframe containing the gff information
    """
    return gff_df[gff_df["feature"] == "CDS"]


def extract_transcript_id(attr_str):
    for attr in attr_str.split(";"):
        # Ensembl GFF3 support
        if attr.startswith("Parent=transcript:") \
                or attr.startswith("ID=transcript:"):
            return attr.split(":")[1]
        # Gencode GFF3 support
        elif attr.startswith("transcript_id="):
            return attr.split("=")[1]
        # Ensembl GTF support
        elif attr.startswith(" transcript_id "):
            return attr.split(" ")[2].replace('"', "")
    return np.nan


def check_annotation(file_path: str) -> bool:
    """
    Checks whether an annotation file exists and is in the right format

    Inputs:
        file_path: Path to the annotation or gff file

    Outputs:
        bool: True if an annotation exists, False otherwise
    """
    if os.path.exists(file_path):
        with open(file_path) as f:
            if "transcript_id" in str(f.readline()):
                return True
            else:
                return False
    else:
        return False



def prepare_annotation(
        gff_path: str,
        outdir: str,
        num_transcripts: int,
        num_processes: int = 4,
        ) -> pd.DataFrame:
    """
    Given a path to a gff file, produce a tsv file containing the
    transcript_id, tx_cds_start, tx_cds_end, tx_length,
    genomic_cds_starts, genomic_cds_ends for each transcript

    Inputs:
        gff_path: Path to the gff file
        outdir: Path to the output directory
        num_transcripts: Number of transcripts to include in the annotation
        num_processes: The maximum number of processes that this function can
                       create

    Outputs:
        annotation_df: Dataframe containing the annotation information
    """
    pool = Pool(processes=num_processes)

    print("Parsing gff..")
    gff_df, coding_tx_ids = parse_gff(gff_path, num_transcripts)
    split_length = (len(gff_df) // num_processes) + 1

    split_df_list = []
    for split in range(num_processes):
        
        split_df_list.append(gff_df.iloc[split_length * split : split_length * (split + 1)])


    print(split_df_list[0].iloc[-10:])
    print("="*20)
    print(split_df_list[1].iloc[0:10])
    exit()
    annotation_batches = []
    print("Subsetting CDS regions, Progress:")
    for split_num, split_df in enumerate(split_df_list):
        annotation_batches.append(gff_df_to_cds_df(split_df, coding_tx_ids, split_num))
        # annotation_batches.append(pool.apply_async(gff_df_to_cds_df,[split_df, coding_tx_ids, split_num]))

    pool.close()
    pool.join()

    print("\n"*(split_num // 4))
    results = [batch.get() for batch in annotation_batches]

    print([len(df) for df in results])

    total_len = 0
    for result in results:
        print(len(result))
        total_len += len(result)
    print(total_len)

    exit("WIP")

    basename = '.'.join(os.path.basename(gff_path).split(".")[:-1])
    output_name = f"{basename}_RiboMetric.tsv"
    annotation_df.to_csv(
        os.path.join(outdir, output_name),
        sep="\t",
        index=False
        )
    return annotation_df


def is_gzipped(file_path: str) -> bool:
    """
    Checks whether the file is gzipped or not

    Inputs:
        file_path: Path to the file to be checked

    Outputs:
        True if gzipped, otherwise False
    """
    try:
        with open(file_path, 'rb') as f:
            # Read the first two bytes of the file
            header = f.read(2)

        # Check if the file starts with the gzip magic number (0x1f 0x8b)
        return header == b'\x1f\x8b'

    except IOError:
        # File not found or unable to open
        return False


def parse_gff(gff_path: str, num_transcripts: int) -> pd.DataFrame:
    """
    Read in the gff file at the provided path and return a dataframe

    Inputs:
        gff_path: Path to the gff file

    Outputs:
        gff_df: Dataframe containing the gff information
    """
    if is_gzipped(gff_path):
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_filepath = temp_file.name
            with gzip.open(gff_path, 'rt') as f:
                for line in f:
                    temp_file.write(line.encode())

        gff_df = gffpd.read_gff3(temp_filepath).df

        os.remove(temp_filepath)

    else:
        gff_df = gffpd.read_gff3(gff_path).df

    gff_df.loc[:, "transcript_id"] = gff_df["attributes"].apply(
        extract_transcript_id
        )

    cds_df = gff_df[gff_df["type"] == "CDS"]
    coding_tx_ids = cds_df["transcript_id"].unique()[:num_transcripts]

    # subset GFF DataFrame to only include transcripts in the transcript_list
    gff_df = gff_df[gff_df["transcript_id"].isin(coding_tx_ids)]

    # Sort GFF DataFrame by transcript ID
    gff_df = gff_df.sort_values("transcript_id")

    return gff_df, coding_tx_ids


def gff_df_to_cds_df(
        gff_df: pd.DataFrame,
        transcript_list: list,
        split_num: int,
        ) -> pd.DataFrame:
    """
    Subset the gff dataframe to only include the CDS features
    with tx coordinates for a specific list of transcripts.

    Inputs:
        gff_df: Dataframe containing the gff information
        transcript_list: List of transcripts to subset

    Outputs:
        cds_df: Dataframe containing the CDS information
                columns: transcript_id, cds_start, cds_end
    """
    # Extract transcript ID from "attributes" column using regular expression
    rows = {
        "transcript_id": [],
        "cds_start": [],
        "cds_end": [],
        "transcript_length": [],
        "genomic_cds_starts": [],
        "genomic_cds_ends": [],
    }

    counter = 0
    for group_name, group_df in gff_df.groupby("transcript_id"):
        counter += 1
        if counter % 100 == 0:
            progress = format_progress((counter / len(transcript_list))*100)
            formatted_num = f"{split_num+1:02d}"
            # print("\n"*(split_num // 4),
            #       "\033[20C"*(split_num % 4),
            #       f"thread {formatted_num}: {progress} | ",
            #       "\033[1A"*(split_num // 4),
            #       end="\r", flush=False, sep="")
        
        if group_name in transcript_list:
            print(group_df)
            transcript_start = group_df["start"].min()

            cds_df_tx = group_df[group_df["type"] == "CDS"]
            cds_start_end_tuple_list = sorted(
                zip(cds_df_tx["start"], cds_df_tx["end"])
                )
            cds_tx_start = cds_start_end_tuple_list[0][0] - transcript_start
            cds_tx_end = cds_start_end_tuple_list[-1][1] - transcript_start

            for cds in cds_start_end_tuple_list:
                cds_length = cds[1] - cds[0]
                cds_tx_end += cds_length

            genomic_cds_starts = ",".join(
                [str(x[0]) for x in cds_start_end_tuple_list]
                )

            genomic_cds_ends = ",".join(
                [str(x[1]) for x in cds_start_end_tuple_list]
                )

            rows["transcript_id"].append(group_name)
            rows["cds_start"].append(cds_tx_start)
            rows["cds_end"].append(cds_tx_end)
            rows["transcript_length"].append(cds_tx_end - cds_tx_start)
            rows["genomic_cds_starts"].append(genomic_cds_starts)
            rows["genomic_cds_ends"].append(genomic_cds_ends)

    return pd.DataFrame(rows)


def format_progress(percentage):
    percentage = round(percentage,3)
    formatted_percentage = "{:.3f}%".format(percentage)
    if len(formatted_percentage) > 7:
        formatted_percentage = "{:.1f}%".format(percentage)
    elif len(formatted_percentage) > 6:
        formatted_percentage = "{:.2f}%".format(percentage)

    return formatted_percentage
