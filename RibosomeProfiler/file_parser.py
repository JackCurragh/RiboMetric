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


def parse_gff(gff_path: str) -> pd.DataFrame:
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


def subset_gff(gff_path: str, transcript_list: list, output_dir: str) -> str:
    """
    Subset the GFF file to only include the transcripts in the provided list

    Inputs:
        gff_path: Path to the annotation file
        transcript_list: List of transcripts to include in the
        subsetted GFF file

    Outputs:
        filtered_gff_path: GFF file containing only the
        transcripts in the provided list
    """
    # read in with pandas
    gff = pd.read_csv(gff_path, sep="\t", header=None, comment="#")
    subsetted_gff = gff[gff[8].str.contains("|".join(transcript_list))]
    subsetted_gff.to_csv(
        f"{output_dir}/subsetted.gff", sep="\t", header=None, index=False
    )

    return f"{output_dir}/subsetted.gff"
