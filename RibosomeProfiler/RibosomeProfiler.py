'''
Main module for RibosomeProfiler
Handles the command line interface and calls the appropriate functions

Many different input combinations are possible. 

Minimal Set:
    -b, --bam <path> : Path to the bam file

With this set the calculations will potentially less reliable and no gene feature information will be included in the output

Standard Set:
    -b, --bam <path> : Path to the bam file
    -g, --gff <path> : Path to the gff file

with this set the calculations will be more reliable and gene feature information will be included in the output

Full Set:
    -b, --bam <path> : Path to the bam file
    -g, --gff <path> : Path to the gff file
    -t, --transcriptome <path> : Path to the transcriptome fasta file

with this set the calculations will contain the post information in its output but will take longest to run


Optional Arguments:
    -n, --name <str> : Name of the sample being analysed (default: filename of bam file)
    -S, --subsample <int> : Number of reads to subsample from the bam file (default: 10000000)
    -T, --transcripts <int> : Number of transcripts to consider (default: 100000)
    -c, --config <path> : Path to the config file (default: config.yaml)

Output:
    --json : Output the results as a json file
    --html : Output the results as an html file (default)
    --csv : Output the results as a csv file
    --all : Output the results as all of the above
'''

import argparse


def argumnet_parser():
    '''
    Parse the command line arguments and return the parser object

    Inputs:
        None

    Outputs:
        parser: ArgumentParser object containing the parsed arguments
    '''
    parser = argparse.ArgumentParser(
        description='''A python command-line utility for the generation of comprehensive reports on the quality of ribosome profiling (Ribo-Seq) datasets''',
    )
    parser.add_argument('-b', '--bam', type=str, required=True, help='Path to the bam file')
    parser.add_argument('-g', '--gff', type=str, required=False, help='Path to the gff file')
    parser.add_argument('-t', '--transcriptome', type=str, required=False, help='Path to the transcriptome fasta file')

    parser.add_argument('-n', '--name', type=str, required=False, help='Name of the sample being analysed (default: filename of bam file)')
    parser.add_argument('-S', '--subsample', type=int, required=False, default=10000000, help='Number of reads to subsample from the bam file (default: 10000000)')
    parser.add_argument('-T', '--transcripts', type=int, required=False, default=100000, help='Number of transcripts to consider (default: 100000)')
    parser.add_argument('-c', '--config', type=str, required=False, default='config.yaml', help='Path to the config file (default: config.yaml)')

    parser.add_argument('--json', action='store_true', default=False, help='Output the results as a json file')
    parser.add_argument('--html', action='store_true', help='Output the results as an html file (default)')
    parser.add_argument('--csv', action='store_true', default=False, help='Output the results as a csv file')
    parser.add_argument('--all', action='store_true', default=False, help='Output the results as all of the above')

    return parser

def main():
    '''
    Main function for the RibosomeProfiler command line interface

    Inputs:
        None

    Outputs:
        None
    '''
    parser = argumnet_parser()
    args = parser.parse_args()

    if args.all:
        args.json = True
        args.html = True
        args.csv = True

    if args.gff is None:
        print("Running annotation free mode")

    elif args.transcriptome is None:
        print("Running annotation mode")
    else:
        print("Running sequence mode")


if __name__ == '__main__':
    main()

