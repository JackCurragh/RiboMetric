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


from .bam_processing import (join_batches,
                             ox_parse_reads,
                             ox_server_parse_reads,
                             validate_bam,
                             )
from .file_splitting import (split_gff_df,
                             run_samtools_idxstats,
                             split_idxstats_df,
                             )


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
    Check whether the bam file exists. Index the bam file if it exists
    Return True if both files exist, False otherwise

    Inputs:
        bam_path: Path to the bam file

    Outputs:
        bool: True if the bam file and its index exist, False otherwise
    """
    if os.path.exists(bam_path):
        if os.path.exists(bam_path + ".bai"):
            return True
        else:
            os.system(f"samtools index {bam_path}")
            if not os.path.exists(bam_path + ".bai"):
                raise Exception("""
                        Indexing failed - Ensure bam is sorted by coordinate
                            samtools sort -o {bam_path} {bam_path}
                                """)
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

            print("\n"*(split_num // 4))

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
            if len(str(f.readline()).split()) == 4:
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
    split_df_list = split_gff_df(gff_df, num_processes)
    del gff_df
    basename = '.'.join(os.path.basename(gff_path).split(".")[:-1])
    output_name = f"{basename}_RiboMetric.tsv"
    open(f"{outdir}/{output_name}", "w").write(
        "transcript_id\tcds_start\tcds_end\ttranscript_length\n")
    annotation_batches = []
    print("Subsetting CDS regions, Progress:")
    for split_df in split_df_list:
        annotation_batches.append(pool.apply_async(
                gff_df_to_cds_df,
                [split_df, f"{outdir}/{output_name}"]
            ))

    pool.close()
    pool.join()

    results = [batch.get() for batch in annotation_batches]

    if all(isinstance(result, pd.DataFrame) for result in results):
        print("Combining results..")
        annotation_df = pd.concat(results, ignore_index=True)
        print("Done")

        return annotation_df
    return pd.DataFrame()


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


def gff_df_to_cds_df(gff_df, outpath=None):
    rows = {
        "transcript_id": [],
        "cds_start": [],
        "cds_end": [],
        "transcript_length": [],
    }

    for group_name, group_df in gff_df.groupby("transcript_id"):
        if group_df.iloc[0]["strand"] == "+":
            cds_start = group_df[group_df["type"] == "CDS"]["start"].min()
            cds_end = group_df[group_df["type"] == "CDS"]["end"].max()
        else:
            cds_start = group_df[group_df["type"] == "CDS"]["end"].max()
            cds_end = group_df[group_df["type"] == "CDS"]["start"].min()

        leader_length, trailer_length, transcript_length = 0, 0, 0

        for idx, exon in group_df[group_df["type"] == "exon"].sort_values(
            "start", ascending=False
                ).iterrows():
            if group_df.iloc[0]["strand"] == "+":
                if exon['end'] <= cds_start:
                    leader_length += exon['end'] - exon['start']
                elif exon['start'] <= cds_start:
                    leader_length += cds_start - exon['start']
                elif exon['start'] >= cds_end and exon['end'] >= cds_end:
                    trailer_length += exon['end'] - exon['start']
                elif exon['end'] >= cds_end and exon['start'] <= cds_end:
                    trailer_length += exon['end'] - cds_end

            elif group_df.iloc[0]["strand"] == "-":
                if exon['start'] >= cds_start:
                    leader_length += exon['end'] - exon['start']
                elif exon['end'] >= cds_start:
                    leader_length += exon['end'] - cds_start
                elif exon['end'] <= cds_end:
                    trailer_length += exon['end'] - exon['start']
                elif exon['start'] <= cds_end:
                    trailer_length += cds_end - exon['start']
            transcript_length += exon['end'] - exon['start']

        cds_start_transcript = leader_length
        cds_end_transcript = transcript_length - trailer_length

        if outpath is None:
            rows["transcript_id"].append(group_name)
            rows["cds_start"].append(cds_start_transcript)
            rows["cds_end"].append(cds_end_transcript)
            rows["transcript_length"].append(transcript_length)
        else:
            with open(outpath, "a") as f:
                f.write(f"{group_name}\t{cds_start_transcript}\t"
                        f"{cds_end_transcript}\t{transcript_length}\n")

    if outpath is None:
        return pd.DataFrame(rows)
