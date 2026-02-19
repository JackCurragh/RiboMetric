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
import subprocess

from multiprocessing import Pool
from tempfile import TemporaryDirectory


from .bam_processing import (join_batches,
                             ox_parse_reads,
                             validate_bam,
                             )
from .file_splitting import (split_gff_df,
                             run_samtools_idxstats,
                             split_idxstats_df,
                             )


def parse_annotation(annotation_path: str) -> pd.DataFrame:
    """
    Read an annotation TSV and return a DataFrame.

    Accepts both minimal and extended schemas. Required columns:
    transcript_id, cds_start, cds_end, transcript_length.
    Optional columns (ignored if missing): genomic_cds_starts, genomic_cds_ends.
    """
    df = pd.read_csv(annotation_path, sep="\t")

    required = {"transcript_id", "cds_start", "cds_end", "transcript_length"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Annotation file missing required columns: {sorted(missing)}"
        )

    # Coerce dtypes for required columns
    df["transcript_id"] = df["transcript_id"].astype(str)
    for c in ["cds_start", "cds_end", "transcript_length"]:
        df[c] = df[c].astype(int)

    # Ensure optional cols exist as strings
    for c in ["genomic_cds_starts", "genomic_cds_ends"]:
        if c not in df.columns:
            df[c] = ""
        df[c] = df[c].astype(str)

    return df


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
            subprocess.run(["samtools", "index", bam_path], check=True)
            if not os.path.exists(bam_path + ".bai"):
                raise Exception(
                    "Indexing failed - Ensure bam is sorted by coordinate:\n"
                    f"    samtools sort -o {bam_path} {bam_path}"
                )
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
    return flagstat_dict


def parse_bam(bam_file: str,
              num_reads: int,
              batch_size: int = 10000000,
              num_processes: int = 4,
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
    # Standard (chunked) parsing path
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
    """Validate that an annotation TSV exists and has required columns."""
    if not os.path.exists(file_path):
        return False
    try:
        head = pd.read_csv(file_path, sep="\t", nrows=1)
    except Exception:
        return False
    required = {"transcript_id", "cds_start", "cds_end", "transcript_length"}
    return required.issubset(set(head.columns))


def prepare_annotation(
        gff_path: str,
        outdir: str,
        num_transcripts: int,
        num_processes: int = 4,
        ) -> pd.DataFrame:
    """
    Given a path to a gff file, produce a tsv file containing the
    transcript_id, tx_cds_start, tx_cds_end, tx_length for each transcript.

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

    annotation_batches = []
    print("Subsetting CDS regions..")
    for split_df in split_df_list:
        annotation_batches.append(pool.apply_async(
                gff_df_to_cds_df,
                [split_df]
            ))

    pool.close()
    pool.join()

    results = [batch.get() for batch in annotation_batches]

    if all(isinstance(result, pd.DataFrame) for result in results):
        print("Combining results..")
        annotation_df = pd.concat(results, ignore_index=True)

        if outdir is not None:
            basename = '.'.join(os.path.basename(gff_path).split(".")[:-1])
            output_name = f"{basename}_RiboMetric.tsv"
            output_path = os.path.join(outdir, output_name)
            annotation_df.to_csv(output_path, sep="\t", index=False)
            print(f"Annotation written to {output_path}")

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

    # Extract transcript_id from the attributes string using a single vectorized
    # regex that handles the common GFF3 and GTF attribute formats:
    #   Ensembl GFF3:  Parent=transcript:ENST...  or  ID=transcript:ENST...
    #   Gencode GFF3:  transcript_id=ENST...
    #   GTF:           transcript_id "ENST..."
    _TRANSCRIPT_ID_RE = (
        r'(?:Parent=transcript:|ID=transcript:|transcript_id=|transcript_id ")'
        r'([^;"]+)'
    )
    gff_df["transcript_id"] = gff_df["attributes"].str.extract(
        _TRANSCRIPT_ID_RE, expand=False
    )

    cds_df = gff_df[gff_df["type"] == "CDS"]
    coding_tx_ids = cds_df["transcript_id"].unique()[:num_transcripts]

    # subset GFF DataFrame to only include transcripts in the transcript_list
    gff_df = gff_df[gff_df["transcript_id"].isin(coding_tx_ids)]

    # Sort GFF DataFrame by transcript ID
    gff_df = gff_df.sort_values("transcript_id")

    return gff_df, coding_tx_ids



def gff_df_to_cds_df(gff_df, outpath=None):
    """
    Vectorized conversion of a GFF dataframe to a CDS annotation dataframe.

    For each transcript, computes transcript-space CDS start/end and total
    transcript length from exon and CDS features.  Handles both + and - strand.

    The outpath parameter is kept for API compatibility but is no longer used
    for incremental writing; the caller (prepare_annotation) handles file I/O.
    """
    exon_df = gff_df[gff_df["type"] == "exon"].copy()
    cds_df = gff_df[gff_df["type"] == "CDS"]

    if exon_df.empty or cds_df.empty:
        return pd.DataFrame(
            columns=["transcript_id", "cds_start", "cds_end", "transcript_length"]
        )

    # Transcript length = sum of exon lengths per transcript
    exon_df["_exon_len"] = exon_df["end"] - exon_df["start"]
    transcript_length = exon_df.groupby("transcript_id")["_exon_len"].sum()

    # Strand per transcript
    strand = gff_df.groupby("transcript_id")["strand"].first()

    # CDS genomic boundaries — strand-aware to match original coordinate convention:
    #   + strand: cgs = min CDS start  (5' end in genomic coords)
    #             cge = max CDS end    (3' end in genomic coords)
    #   - strand: cgs = max CDS end    (5' end — highest genomic coord)
    #             cge = min CDS start  (3' end — lowest genomic coord)
    cds_min_start = cds_df.groupby("transcript_id")["start"].min()
    cds_max_end = cds_df.groupby("transcript_id")["end"].max()

    is_plus = strand == "+"
    cgs = pd.Series(index=strand.index, dtype=float)  # genomic CDS start
    cge = pd.Series(index=strand.index, dtype=float)  # genomic CDS end
    cgs[is_plus] = cds_min_start.reindex(strand.index[is_plus])
    cgs[~is_plus] = cds_max_end.reindex(strand.index[~is_plus])
    cge[is_plus] = cds_max_end.reindex(strand.index[is_plus])
    cge[~is_plus] = cds_min_start.reindex(strand.index[~is_plus])

    # Join boundaries and strand onto exon rows
    exon_df = exon_df.join(cgs.rename("_cgs"), on="transcript_id")
    exon_df = exon_df.join(cge.rename("_cge"), on="transcript_id")
    exon_df = exon_df.join(strand.rename("_strand"), on="transcript_id")

    plus = exon_df["_strand"] == "+"
    s, e = exon_df["start"], exon_df["end"]

    # Leader contribution per exon (exon nucleotides before the CDS start):
    #   + strand: max(0, min(exon_end, cgs) - exon_start)
    #   - strand: max(0, exon_end - max(exon_start, cgs))
    leader = pd.Series(0.0, index=exon_df.index)
    leader[plus] = np.maximum(
        0, np.minimum(e[plus], exon_df.loc[plus, "_cgs"]) - s[plus]
    )
    leader[~plus] = np.maximum(
        0, e[~plus] - np.maximum(s[~plus], exon_df.loc[~plus, "_cgs"])
    )

    # Trailer contribution per exon (exon nucleotides after the CDS end):
    #   + strand: max(0, exon_end - max(exon_start, cge))
    #   - strand: max(0, min(exon_end, cge) - exon_start)
    trailer = pd.Series(0.0, index=exon_df.index)
    trailer[plus] = np.maximum(
        0, e[plus] - np.maximum(s[plus], exon_df.loc[plus, "_cge"])
    )
    trailer[~plus] = np.maximum(
        0, np.minimum(e[~plus], exon_df.loc[~plus, "_cge"]) - s[~plus]
    )

    exon_df["_leader"] = leader
    exon_df["_trailer"] = trailer

    leader_sum = exon_df.groupby("transcript_id")["_leader"].sum()
    trailer_sum = exon_df.groupby("transcript_id")["_trailer"].sum()

    tx_ids = transcript_length.index
    result = pd.DataFrame({
        "transcript_id": tx_ids,
        "cds_start": leader_sum.reindex(tx_ids).fillna(0).astype(int).values,
        "cds_end": (
            transcript_length - trailer_sum.reindex(tx_ids).fillna(0)
        ).astype(int).values,
        "transcript_length": transcript_length.astype(int).values,
    })

    return result
