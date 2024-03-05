"""

"""

import subprocess
import pandas as pd
import numpy as np


def run_samtools_idxstats(bam_file: str) -> pd.DataFrame:
    """
    Run 'samtools idxstats' for the bam file and return a dataframe with
    the results

    Inputs:
        bam_file: Path to the bam file

    Outputs:
        idxstats_df: dataframe containing idxstats for the bam file
    """
    # Run samtools idxstats command and capture the output
    process = subprocess.Popen(['samtools', 'idxstats', bam_file],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # Convert the output to a pandas DataFrame
    lines = stdout.decode().strip().split('\n')
    data = [line.split('\t') for line in lines]
    df = pd.DataFrame(data, columns=['Reference',
                                     'Length',
                                     'Mapped_Reads',
                                     'Unmapped_Reads'])

    return df


def split_idxstats_df(idxstats_df: pd.DataFrame,
                      batch_size: int,
                      num_reads: int) -> list:
    """
    Split the idxstats data frame into a list of data frames limited to
    the max read count while also preparing the dataframe for conversion to
    bed files

    Inputs:
        idxstats_df: Results from samtools idxstats
        batch_size: Maximum number of reads in a batch
        num_reads: Number of reads to parse

    Outputs:
        split_dfs
    """
    # Convert the 'Mapped_Reads' and 'Unmapped_Reads' columns to numpy arrays
    mapped_reads = idxstats_df['Mapped_Reads'].astype(int).values
    unmapped_reads = idxstats_df['Unmapped_Reads'].astype(int).values

    split_dfs = []
    current_sum = 0
    current_df = pd.DataFrame()
    last_index = 0
    split_num = 0

    for i in range(len(mapped_reads)):
        reads = mapped_reads[i] + unmapped_reads[i]
        if current_sum + reads <= batch_size:
            current_sum += reads
            if (current_sum + (split_num * batch_size)) > num_reads:
                break
        else:
            current_df = idxstats_df.iloc[last_index:i, [0, 1]].copy()
            current_df['Start'] = np.zeros(i - last_index, dtype=np.int8)
            current_df = current_df[['Reference',
                                     'Start',
                                     'Length']]
            split_dfs.append(current_df)
            split_num += 1
            current_df = pd.DataFrame()
            last_index = i
            current_sum = reads

    # Add the last remaining DataFrame
    current_df = idxstats_df.iloc[last_index:i, [0, 1]].copy()
    current_df['Start'] = np.zeros(i - last_index, dtype=np.int8)
    current_df = current_df[['Reference',
                             'Start',
                             'Length']]
    split_dfs.append(current_df)

    return split_dfs


def split_bam(bam_file: str,
              split_num: int,
              reference_df: pd.DataFrame,
              tempdir: str
              ) -> str:
    """
    Splits the bam files with the bed files generated from the idxstats

    Inputs:
        bam_file: Path to the bam file
        split_num: Count of splits
        reference_dfs: Dataframes containing the reference names and start and
        stop of the transcripts
        tempdir: Path to the temp directory

    Outputs:
        outfile: Path to the split bam file
    """
    bedfile = f"{tempdir}/bed_{split_num}.bed"
    reference_df.to_csv(bedfile, sep="\t", header=False, index=False)
    outfile = f"{tempdir}/split_sorted_{split_num}.bam"
    samview = subprocess.Popen(('samtools',
                                'view',
                                '-F 256',
                                '-h',
                                '-L',
                                bedfile,
                                bam_file),
                               stdout=subprocess.PIPE)
    with open("/dev/null", "w") as stderr:
        subprocess.run(('samtools',
                        'sort',
                        '-O',
                        'bam',
                        '-o',
                        outfile),
                       stdin=samview.stdout,
                       stderr=stderr)
    samview.wait()
    subprocess.run(('samtools',
                    'index',
                    outfile))

    return outfile


def split_gff_file(input_file, outdir, num_files) -> list:
    with open(input_file, 'r') as f:
        gff_lines = f.readlines()

    num_lines = len(gff_lines)
    lines_per_split = num_lines // num_files

    split_indices = [0]
    line_number = 0
    for i in range(num_lines):
        line_number += 1
        line = gff_lines[i]
        if line_number >= lines_per_split:
            if line.startswith("#"):
                continue
            if line.split('\t')[2] == "gene":
                split_indices.append(i)
                line_number = 0

    split_indices.append(num_lines)
    split_gffs = []
    for i in range(len(split_indices)-1):
        start_index = split_indices[i]
        end_index = split_indices[i + 1]
        file_lines = gff_lines[start_index:end_index]

        output_file = f"{outdir}/temp_gff_{i + 1}.gff"
        with open(output_file, 'w') as f_out:
            f_out.writelines(file_lines)

        split_gffs.append(output_file)

    return split_gffs


def split_gff_df(gff_df: pd.DataFrame, split_num: int) -> list:
    """
    Split the gff dataframe into a list of dataframes limited to the max

    Inputs:
        gff_df: Dataframe containing the gff data
        split_num: Number of splits

    Outputs:
        split_df_list: List of dataframes
    """
    df_length = len(gff_df)

    if df_length < split_num:
        split_num = 1

    split_length = (df_length // split_num) + 1
    prev_offset, offset = 0, 0
    split_df_list = []
    for split in range(split_num):

        lower_limit = (split_length * split)
        upper_limit = (split_length * (split + 1))

        if upper_limit + offset >= df_length:
            split_df = gff_df.iloc[lower_limit + prev_offset:
                                   df_length]
            split_df_list.append(split_df)
            return split_df_list

        last_transcript_id = gff_df["transcript_id"].iloc[upper_limit
                                                          - 1
                                                          + offset]
        next_transcript_id = gff_df["transcript_id"].iloc[upper_limit
                                                          + offset]

        while last_transcript_id == next_transcript_id:

            if upper_limit + offset + 1 > df_length:
                split_df = gff_df.iloc[lower_limit + prev_offset:
                                       df_length]
                split_df_list.append(split_df)
                return split_df_list

            offset += 1
            last_transcript_id = gff_df["transcript_id"].iloc[upper_limit
                                                              - 1
                                                              + offset]
            next_transcript_id = gff_df["transcript_id"].iloc[upper_limit
                                                              + offset]

        # print("="*50)
        # print(last_transcript_id, next_transcript_id, offset)

        split_df = gff_df.iloc[lower_limit + prev_offset:
                               upper_limit + offset]
        prev_offset = offset
        split_df_list.append(split_df)

    return split_df_list


def format_progress(percentage):
    round_percentage = round(percentage, 3)
    formatted_percentage = "{:.3f}%".format(round_percentage)

    if len(formatted_percentage) > 7:
        round_percentage = round(percentage, 1)
        formatted_percentage = "{:.1f}%".format(round_percentage)

    elif len(formatted_percentage) > 6:
        round_percentage = round(percentage, 2)
        formatted_percentage = "{:.2f}%".format(round_percentage)

    return formatted_percentage
