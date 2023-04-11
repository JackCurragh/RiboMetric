'''
Main script for running qc analysis

Three main modes:
    annoation free: no gff file provided just use the bam file
    annotation based: gff file provided and use the bam file
    sequence based: gff file and transcriptome fasta file provided and use the bam file

'''

def annotation_free_mode(read_dict: dict) -> dict:
    '''
    Run the annotation free mode of the qc analysis

    Inputs:
        read_dict: Dictionary containing the read information (keys are the read names)

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    '''
    results_dict = {}
    
    return results_dict


def annotation_mode(read_dict: dict, gff_path: str, transcript_list: list) -> dict:
    '''
    Run the annotation mode of the qc analysis

    Inputs:
        read_dict: Dictionary containing the read information (keys are the read names)
        gff_path: Path to the gff file
        transcript_list: List of the top N transcripts

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    '''
    results_dict = {}
    
    return results_dict


def sequence_mode(read_dict: dict, gff_path: str, transcript_list: list, fasta_path: str) -> dict:
    '''
    Run the sequence mode of the qc analysis

    Inputs:
        read_dict: Dictionary containing the read information (keys are the read names)
        gff_path: Path to the gff file
        transcript_list: List of the top N transcripts
        fasta_path: Path to the transcriptome fasta file

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    '''
    results_dict = {}
    
    return results_dict



