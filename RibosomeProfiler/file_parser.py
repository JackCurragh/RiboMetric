'''
This script contains the functions used to load the required files in the RibosomeProfiler pipeline

The functions are called by the main script RibosomeProfiler.py
'''
from Bio import SeqIO
from gffutils import create_db, FeatureDB
import pysam
import pandas as pd
import random

def parse_fasta(fasta_path: str) -> dict:
    '''
    Read in the transcriptome fasta file at the provided path and return a dictionary

    Inputs:
        fasta_path: Path to the transcriptome fasta file

    Outputs:
        transcript_dict: Dictionary containing the transcript information
    '''
    #read in with biopython
    transcript_dict = {}
    for record in SeqIO.parse(fasta_path, 'fasta'):
        transcript_dict[record.id] = record

    return transcript_dict


def flagstat_bam(bam_path: str) -> dict:
    '''
    Run samtools flagstat on the bam file at the provided path and return a dictionary

    Inputs:
        bam_path: Path to the bam file

    Outputs:
        flagstat_dict: Dictionary containing the flagstat information
    
    '''
    flagstat_dict = {}
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        flagstat_dict['total_reads'] = bamfile.mapped + bamfile.unmapped
        flagstat_dict['mapped_reads'] = bamfile.mapped
        flagstat_dict['unmapped_reads'] = bamfile.unmapped
        flagstat_dict['duplicates'] = bamfile.mapped + bamfile.unmapped
    return flagstat_dict


def parse_bam(bam_path: str, num_reads: int) -> dict:
    '''
    Read in the bam file at the provided path and return a dictionary

    Inputs:
        bam_path: Path to the bam file
        num_reads: Number of reads to consider

    Outputs:
        read_dict: Dictionary containing the read information (keys are the read names)
    '''
    read_dict = {}
    count = pysam.view("-c", f"{bam_path}")
    fraction = num_reads / int(count)

    pysam.view("-s", f"123.{str(fraction).split('.')[1]}", "-bo",  f"{bam_path}.subsampled.bam", f"{bam_path}")
    pysam.index(f"{bam_path}.subsampled.bam")

    with pysam.AlignmentFile(f"{bam_path}.subsampled.bam", "rb") as bamfile:
        fraction = 0.1

        # Use the view() function to create a generator that yields a fraction of the reads

        # Loop over the subsampled reads and do something with each read
        for read in bamfile.fetch():
            # Do something with the read here
            if read.query_name not in read_dict:
                read_dict[read.query_name] = read

    return read_dict

def get_top_transcripts(read_dict: dict, num_transcripts: int) -> list:
    '''
    Get the top N transcripts based on the number of reads mapping to them

    Inputs:
        read_dict: Dictionary containing the read information (keys are the read names)
        num_transcripts: Number of transcripts to consider

    Outputs:
        top_transcripts: List of the top N transcripts
    '''
    #get the top N transcripts based on the number of reads mapping to them
    transcript_dict = {}
    for read in read_dict.values():
        if read.reference_name in transcript_dict:
            transcript_dict[read.reference_name] += 1
        else:
            transcript_dict[read.reference_name] = 1

    top_transcripts = sorted(transcript_dict, key=transcript_dict.get, reverse=True)[:num_transcripts]

    return top_transcripts


def subset_gff(gff_path: str, transcript_list: list, output_dir: str) -> str:
    '''
    Subset the GFF file to only include the transcripts in the provided list

    Inputs:
        gff_path: Path to the annotation file
        transcript_list: List of transcripts to include in the subsetted GFF file

    Outputs:
        filtered_gff_path: GFF file containing only the transcripts in the provided list
    '''
    #read in with pandas
    gff = pd.read_csv(gff_path, sep='\t', header=None, comment='#')
    subsetted_gff = gff[gff[8].str.contains('|'.join(transcript_list))]
    subsetted_gff.to_csv(f'{output_dir}/subsetted.gff', sep='\t', header=None, index=False)

    return f'{output_dir}/subsetted.gff'


def parse_gff(gff_path: str) -> dict:
    '''
    Read in the annotation file at the provided path and return a dictionary

    Inputs:
        gff_path: Path to the annotation file

    Outputs:
        gene_dict: Dictionary containing the gene information
    '''
    #read in with gffutils
    db = create_db(gff_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    gene_dict = {}
    for gene in db.features_of_type('gene'):
        gene_dict[gene.id] = gene

    return gene_dict
