"""
This script contains the functions used to load the required files
in the RibosomeProfiler pipeline

The functions are called by the main script RibosomeProfiler.py
"""
from Bio import SeqIO
import pysam
import pandas as pd
import subprocess
import gffpandas.gffpandas as gffpd


def parse_gff(gff_path: str) -> gffpd.Gff3DataFrame:
    """
    Read in the gff file at the provided path and return a dataframe

    Inputs:
        gff_path: Path to the gff file

    Outputs:
        gff_df: Dataframe containing the gff information
    """
    return gffpd.read_gff3(gff_path)


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


def parse_bam(bam_file: str, num_reads: int) -> pd.DataFrame:
    """
    Read in the bam file at the provided path and return a dictionary

    Inputs:
        bam_file: Path to the bam file
        Num_reads: Number of reads to parse

    Outputs:
        read_dict: Dictionary containing the read information
        (keys are the read names)
    """
    # Convert the BAM file to SAM format and read the output in chunks
    cmd = f"samtools view {bam_file}"
    print(f"Running {cmd}")
    process = subprocess.Popen(cmd,
                               shell=True,
                               stdout=subprocess.PIPE,
                               text=True)
    print("Processing reads...")
    read_list = []
    for line in iter(process.stdout.readline, ""):
        if line.startswith("@"):
            continue
        fields = line.strip().split("\t")

        if "_x" in fields[0]:
            count = int(fields[0].split("_x")[-1])
        else:
            count = 1

        read_list.append(
            {
                "read_name": fields[0],
                "read_length": len(fields[9]),
                "reference_name": fields[2],
                "reference_start": int(fields[3]) - 1,
                "sequence": fields[9],
                "sequence_qualities": fields[10],
                "tags": fields[11:],
                "count": count,
            }
        )

        if len(read_list) > num_reads:
            process.kill()  # kill the process if we've read enough data
            break
        else:
            read_percentage = round(len(read_list)/num_reads, 3) * 100
            print(
                f"Processed {len(read_list)}/{num_reads}({read_percentage}%)",
                end="\r",
            )

    process.kill()
    return pd.DataFrame(read_list)


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
        read_df.groupby("reference_name").sum().sort_values(
            "count",
            ascending=False
            )
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
    transcript_id_regex = r"transcript_id=([^;]+)"
    gff_df.loc[:, "transcript_id"] = gff_df["attributes"].str.extract(
        transcript_id_regex
        )
    rows = {"transcript_id": [], "cds_start": [], "cds_end": []}

    # Sort GFF DataFrame by transcript ID
    gff_df = gff_df.sort_values("transcript_id")

    # subset GFF DataFrame to only include transcripts in the transcript_list
    gff_df = gff_df[gff_df["transcript_id"].isin(transcript_list)]

    counter = 0
    for group_name, group_df in gff_df.groupby("transcript_id"):
        counter += 1
        if counter % 100 == 0:
            prop = counter/len(transcript_list)
            print(f"Processing transcript ({prop*100}%)", end="\r")
        if group_name in transcript_list:
            transcript_start = group_df["start"].min()

            cds_df_tx = group_df[group_df["type"] == "CDS"]
            cds_start_end_tuple_list = sorted(zip(
                cds_df_tx["start"],
                cds_df_tx["end"]
                ))
            cds_tx_start = cds_start_end_tuple_list[0][0] - transcript_start
            cds_tx_end = cds_start_end_tuple_list[-1][1] - transcript_start
            for cds in cds_start_end_tuple_list:
                cds_length = cds[1] - cds[0]
                cds_tx_end += cds_length

            rows["transcript_id"].append(group_name)
            rows["cds_start"].append(cds_tx_start)
            rows["cds_end"].append(cds_tx_end)

    return pd.DataFrame(rows)
