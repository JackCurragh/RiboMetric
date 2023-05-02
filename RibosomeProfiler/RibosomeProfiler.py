"""
Main module for RibosomeProfiler
Handles the command line interface and calls the appropriate functions

Many different input combinations are possible.

Minimal Set:
    -b, --bam <path> : Path to the bam file

With this set the calculations will potentially less reliable and no gene
feature information will be included in the output

Standard Set:
    -b, --bam <path> : Path to the bam file
    -g, --gff <path> : Path to the gff file

with this set the calculations will be more reliable and gene feature
information will be included in the output

Full Set:
    -b, --bam <path> : Path to the bam file
    -g, --gff <path> : Path to the gff file
    -t, --transcriptome <path> : Path to the transcriptome fasta file

with this set the calculations will contain the post information in its
output but will take longest to run


Optional Arguments:
    -n, --name <str> : Name of the sample being analysed
                        (default: filename of bam file)
    -S, --subsample <int> : Number of reads to subsample from the bam file
                        (default: 10000000)
    -T, --transcripts <int> : Number of transcripts to consider
                        (default: 100000)
    -c, --config <path> : Path to the config file
                        (default: config.yaml)

Output:
    --json : Output the results as a json file
    --html : Output the results as an html file (default)
    --csv : Output the results as a csv file
    --all : Output the results as all of the above
"""

import argparse
from rich.console import Console
from rich.text import Text
from rich.table import Table
from rich.emoji import Emoji

from file_parser import (
    parse_bam,
    get_top_transcripts,
    subset_gff,
    parse_fasta,
    parse_gff,
)
from qc import annotation_free_mode, annotation_mode, sequence_mode


def print_logo(console):
    """
    print the logo to the console
    """
    logo = Text(
        """
    ██████╗  ██╗ ██████╗  ██████╗ ███████╗ ██████╗ ███╗   ███╗███████╗
    ██╔══██╗ ██║ ██╔══██╗██╔═══██╗██╔════╝██╔═══██╗████╗ ████║██╔════╝
    ██████╔╝ ██║ ██████╔╝██║   ██║███████╗██║   ██║██╔████╔██║█████╗
    ██╔══██╗ ██║ ██╔══██╗██║   ██║╚════██║██║   ██║██║╚██╔╝██║██╔══╝
    ██║  ██║ ██║ ██████╔╝╚██████╔╝███████║╚██████╔╝██║ ╚═╝ ██║███████╗
    ╚═╝  ╚═╝ ╚═╝ ══════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝
    """,
        style="bold blue",
    )
    logo += Text(
        """
    ██████╗  ██████╗  ██████╗ ███████╗██╗██╗    ███████╗██████╗
    ██╔══██╗ ██╔══██╗██╔═══██╗██╔════╝██║██║    ██╔════╝██╔══██╗
    ██████╔╝ ██████╔╝██║   ██║█████╗  ██║██║    █████╗  ██████╔╝
    ██╔═══╝  ██╔══██╗██║   ██║██╔══╝  ██║██║    ██╔══╝  ██╔══██╗
    ██║      ██║  ██║╚██████╔╝██║     ██║██████╗███████╗██║  ██║
    ╚═╝      ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚═════╝╚══════╝╚═╝  ╚═╝
    """,
        style="bold red",
    )
    console.print(logo)


def print_table(args, console, mode):
    console = Console()

    Inputs = Table(show_header=True, header_style="bold magenta")
    Inputs.add_column("Parameters", style="dim", width=20)
    Inputs.add_column("Values")
    Inputs.add_row("Bam File:", args.bam)
    Inputs.add_row("Gff File:", args.gff)
    Inputs.add_row("Transcriptome File:", args.fasta)

    Configs = Table(show_header=True, header_style="bold yellow")
    Configs.add_column("Options", style="dim", width=20)
    Configs.add_column("Values")
    Configs.add_row("Mode:", mode)
    Configs.add_row("# of reads:", str(args.subsample))
    Configs.add_row("# of transcripts:", str(args.transcripts))
    Configs.add_row("Config file:", args.config)

    Output = Table(show_header=True, header_style="bold blue")
    Output.add_column("Output Options", style="dim", width=20)
    Output.add_column("Values")
    Output.add_row("JSON:", str(args.json))
    Output.add_row("HTML:", str(args.html))
    Output.add_row("CSV:", str(args.csv))
    Output.add_row("All:", str(args.all))

    # Print tables side by side
    console.print(Inputs, Configs, Output, justify="inline", style="bold")


def argumnet_parser():
    """
    Parse the command line arguments and return the parser object

    Inputs:
        None

    Outputs:
        parser: ArgumentParser object containing the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="""A python command-line utility for the generation of
                    comprehensive reports on the quality of ribosome profiling
                    (Ribo-Seq) datasets""",
        epilog=f"""

        Made with {Emoji('heart')} in LAPTI lab at University College Cork.
        For more information, please visit:
        https://ribosomeprofiler.readthedocs.io/en/latest/
        """,
    )
    parser.add_argument(
        "-b", "--bam", type=str, required=True, help="Path to the bam file"
    )
    parser.add_argument(
        "-g", "--gff", type=str, required=False, help="Path to the gff file"
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=False,
        help="Path to the transcriptome fasta file",
    )

    parser.add_argument(
        "-n",
        "--name",
        type=str,
        required=False,
        help="""Name of the sample being analysed
        (default: filename of bam file)""",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default=".",
        help="""Path to the output directory
        (default: current directory)""",
    )
    parser.add_argument(
        "-S",
        "--subsample",
        type=int,
        required=False,
        default=1000000,
        help="""Number of reads to subsample from the bam file
        (default: 10000000)""",
    )
    parser.add_argument(
        "-T",
        "--transcripts",
        type=int,
        required=False,
        default=100000,
        help="Number of transcripts to consider (default: 100000)",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=False,
        default="config.yaml",
        help="Path to the config file (default: config.yaml)",
    )

    parser.add_argument(
        "--json",
        action="store_true",
        default=False,
        help="Output the results as a json file",
    )
    parser.add_argument(
        "--html",
        action="store_true",
        default=True,
        help="Output the results as an html file (default)",
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        default=False,
        help="Output the results as a csv file",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        default=False,
        help="Output the results as all of the above",
    )

    return parser


def main(args):
    """
    Main function for the RibosomeProfiler command line interface

    Inputs:
        args: ArgumentParser object containing the parsed arguments

    Outputs:
        None
    """

    console = Console()
    print_logo(console)
    if args.all:
        args.json = True
        args.html = True
        args.csv = True

    if args.gff is None:
        print_table(args, console, "Annotation Free Mode")

    elif args.fasta is None:
        print_table(args, console, "Annotation Mode")

    else:
        print_table(args, console, "Sequence Mode")

    read_df = parse_bam(args.bam, args.subsample)
    if args.gff is None:
        results_dict = annotation_free_mode(read_df, args.config)

    elif args.fasta is None:
        top_transcripts = get_top_transcripts(read_df, args.transcripts)
        gff_path = subset_gff(args.gff, top_transcripts)
        gffdf = parse_gff(gff_path)
        results_dict = annotation_mode(read_df,
                                       gffdf.df,
                                       top_transcripts,
                                       args.config)

    else:
        fasta_dict = parse_fasta(args.fasta)

        top_transcripts = get_top_transcripts(read_df, args.transcripts)
        gff_path = subset_gff(args.gff, top_transcripts, args.output)
        gffdf = parse_gff(gff_path)
        print(gffdf.df)

        results_dict = annotation_mode(read_df,
                                       gffdf.df,
                                       top_transcripts,
                                       args.config)

        results_dict = sequence_mode(results_dict,
                                     read_df,
                                     fasta_dict,
                                     args.config)


if __name__ == "__main__":
    parser = argumnet_parser()
    args = parser.parse_args()

    main(args)
