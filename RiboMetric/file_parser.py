"""
This script contains the functions used to load the required files
in the RiboMetric pipeline

The functions are called by the main script RiboMetric.py
"""
from Bio import SeqIO
import pysam
import pandas as pd
from multiprocessing import Pool
import gffpandas.gffpandas as gffpd
import os
import numpy as np

from .bam_processing import process_reads, process_sequences


def parse_gff(gff_path: str) -> gffpd.Gff3DataFrame:
    """
    Read in the gff file at the provided path and return a dataframe

    Inputs:
        gff_path: Path to the gff file

    Outputs:
        gff_df: Dataframe containing the gff information
    """
    return gffpd.read_gff3(gff_path)


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
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        flagstat_dict["total_reads"] = bamfile.mapped + bamfile.unmapped
        flagstat_dict["mapped_reads"] = bamfile.mapped
        flagstat_dict["unmapped_reads"] = bamfile.unmapped
        flagstat_dict["duplicates"] = bamfile.mapped + bamfile.unmapped
    return flagstat_dict


def parse_bam(bam_file,
              num_reads=1000000,
              batch_size=100000,
              num_processes=1
              ) -> list:
    """
    Read in the bam file at the provided path and return a list of dataframes

    Inputs:
        bam_file: Path to the bam file
        num_reads: Number of reads to parse
        batch_size: The number of reads that are processed at a time
        num_processes: The maximum number of processes that this function can
        create

    Outputs:
        batch_results: List containing dataframes for the parsed reads which
        will be grouped together in following steps
    """
    samfile = pysam.AlignmentFile(bam_file, "rb")
    pool = Pool(processes=num_processes)
    read_list, batch_results = [], []
    for idx, read in enumerate(samfile.fetch()):
        read_list.append(read.to_string().split(sep="\t"))
        if idx >= num_reads - 1:
            break

        if len(read_list) == batch_size:
            batch_results.append(pool.apply_async(process_reads, [read_list]))
            read_list = []
        read_percentage = round((idx) / num_reads * 100, 3)
        print(f"Processed {idx}/{num_reads} \
({read_percentage}%)", end="\r")

    if read_list:
        batch_results.append(pool.apply_async(process_reads, [read_list]))

    pool.close()
    pool.join()

    return [result.get() for result in batch_results]


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


def gff_df_to_cds_df(
        gff_df: pd.DataFrame,
        transcript_list: list
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

    # Sort GFF DataFrame by transcript ID
    gff_df = gff_df.sort_values("transcript_id")

    # subset GFF DataFrame to only include transcripts in the transcript_list
    gff_df = gff_df[gff_df["transcript_id"].isin(transcript_list)]

    counter = 0
    for group_name, group_df in gff_df.groupby("transcript_id"):
        counter += 1
        if counter % 100 == 0:
            prop = counter / len(transcript_list)
            print(f"Processing transcript ({prop*100}%)", end="\r")
        if group_name in transcript_list:
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


def prepare_annotation(
    gff_path: str, outdir: str, num_transcripts: int, config: str
) -> pd.DataFrame:
    """
    Given a path to a gff file, produce a tsv file containing the
    transcript_id, tx_cds_start, tx_cds_end, tx_length,
    genomic_cds_starts, genomic_cds_ends for each transcript

    Inputs:
        gff_path: Path to the gff file
        outdir: Path to the output directory
        num_transcripts: Number of transcripts to include in the annotation
        config: Path to the config file

    Outputs:
        annotation_df: Dataframe containing the annotation information
    """
    print("Parsing gff")
    gffdf = parse_gff(gff_path).df

    # transcript_id_regex = r"transcript_id=([^;]+)"
    # gffdf.loc[:, "transcript_id"] = gffdf["attributes"].str.extract(
    # transcript_id_regex
    # )

    gffdf.loc[:, "transcript_id"] = gffdf["attributes"].apply(
        extract_transcript_id
        )

    cds_df = gffdf[gffdf["type"] == "CDS"]

    coding_tx_ids = cds_df["transcript_id"].unique()[:num_transcripts]

    annotation_df = gff_df_to_cds_df(gffdf, coding_tx_ids)
    basename = '.'.join(os.path.basename(gff_path).split(".")[:-1])
    output_name = f"{basename}_RiboMetric.tsv"
    annotation_df.to_csv(
        os.path.join(outdir, output_name),
        sep="\t",
        index=False
        )
    return annotation_df
