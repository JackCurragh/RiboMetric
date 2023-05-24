"""
Main module for RiboMetric
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
    --pdf : Output the results as a pdf file
    --csv : Output the results as a csv file
    --all : Output the results as all of the above
"""

import argparse
from rich.console import Console
from rich.text import Text
from rich.table import Table
from rich.emoji import Emoji

import yaml
import os
import numpy as np
import pandas as pd

from .file_parser import (
    parse_bam,
    parse_fasta,
    parse_annotation,
    prepare_annotation,
    flagstat_bam,
    check_bam,
)
from .qc import annotation_mode, sequence_mode
from .plots import generate_plots
from .modules import a_site_calculation
from .html_report import generate_report
from .json_output import generate_json


def print_logo(console):
    """
    print the logo to the console
    """
    logo = Text(
        """
                   ██████╗  ██╗ ██████╗  ██████╗
                   ██╔══██╗ ██║ ██╔══██╗██╔═══██╗
                   ██████╔╝ ██║ ██████╔╝██║   ██║
                   ██╔══██╗ ██║ ██╔══██╗██║   ██║
                   ██║  ██║ ██║ ██████╔╝╚██████╔╝
                   ╚═╝  ╚═╝ ╚═╝ ══════╝  ╚═════╝
    """,
        style="bold blue",
    )
    logo += Text(
        """
    ███╗   ███╗███████╗█████████╗██████╗  ██╗  ██████╗  ███████╗
    ████╗ ████║██╔════╝╚══██╔═══╝██╔══██╗ ██║ ██╔════╝ ██╔═════╝
    ██╔████╔██║█████╗     ██║    ██████╔╝ ██║ ██║      ████████╗
    ██║╚██╔╝██║██╔══╝     ██║    ██╔══██╗ ██║ ██║      ╚═════██║
    ██║ ╚═╝ ██║███████╗   ██║    ██║  ██║ ██║ ╚██████╗ ████████║
    ╚═╝     ╚═╝╚══════╝   ╚═╝    ╚═╝  ╚═╝ ╚═╝  ╚═════╝ ╚═══════╝
    """,
        style="bold red",
    )
    console.print(logo)


def print_table_run(args, console, mode):
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
    Output.add_row("PDF:", str(args.pdf))
    Output.add_row("CSV:", str(args.csv))
    Output.add_row("All:", str(args.all))

    # Print tables side by side
    console.print(Inputs, Configs, Output, justify="inline", style="bold")


def print_table_prepare(args, console, mode):
    console = Console()

    Inputs = Table(show_header=True, header_style="bold magenta")
    Inputs.add_column("Parameters", style="dim", width=20)
    Inputs.add_column("Values")
    Inputs.add_row("Gff File:", args.gff)

    Configs = Table(show_header=True, header_style="bold yellow")
    Configs.add_column("Options", style="dim", width=20)
    Configs.add_column("Values")
    Configs.add_row("Mode:", mode)
    Configs.add_row("# of transcripts:", str(args.transcripts))
    Configs.add_row("Config file:", args.config)

    # Print tables side by side
    console.print(Inputs, Configs, justify="inline", style="bold")


def argument_parser():
    """
    Parse the command line arguments and return the parser object

    Inputs:
        None

    Outputs:
        parser: ArgumentParser object containing the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="""A python command-line utility for the generation
                        of comprehensive reports on the quality of ribosome
                        profiling (Ribo-Seq) datasets""",
        epilog=f"""

            Made with {Emoji('heart')} in LAPTI lab at University College Cork.
            For more information, please visit:
            https://RiboMetric.readthedocs.io/en/latest/
            """,
    )

    subparsers = parser.add_subparsers(dest="command", title="subcommands")

    # create the parser for the "run" command
    run_parser = subparsers.add_parser(
        "run", help="run RiboMetric in normal mode"
    )
    run_parser.add_argument(
        "-b", "--bam", type=str, required=True, help="Path to bam file"
    )
    run_parser.add_argument(
        "-a",
        "--annotation",
        type=str,
        required=False,
        help="Path to RiboMetric annotation file",
    )
    run_parser.add_argument(
        "-g", "--gff", type=str, required=False, help="Path to gff file"
    )
    run_parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=False,
        help="Path to the transcriptome fasta file",
    )
    run_parser.add_argument(
        "-n",
        "--name",
        type=str,
        required=False,
        help="""Name of the sample being analysed
            (default: filename of bam file)""",
    )
    run_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default=".",
        help="""Path to the output directory
            (default: current directory)""",
    )
    run_parser.add_argument(
        "-S",
        "--subsample",
        type=int,
        required=False,
        default=1000000,
        help="""Number of reads to subsample from the bam file
            (default: 10000000)""",
    )
    run_parser.add_argument(
        "-T",
        "--transcripts",
        type=int,
        required=False,
        default=100000,
        help="Number of transcripts to consider (default: 100000)",
    )
    run_parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=False,
        default="config.yml",
        help="Path to the config file (default: config.yml)",
    )
    run_parser.add_argument(
        "--json",
        action="store_true",
        default=False,
        help="Output the results as a json file",
    )
    run_parser.add_argument(
        "--html",
        action="store_true",
        default=True,
        help="Output the results as an html file (default)",
    )
    run_parser.add_argument(
        "--no-html",
        dest='html',
        action='store_false',
        help="Disable html output",
    )
    run_parser.add_argument(
        "--pdf",
        action="store_true",
        default=False,
        help="Output the results as a pdf file (default)",
    )
    run_parser.add_argument(
        "--csv",
        action="store_true",
        default=False,
        help="Output the results as a csv file",
    )
    run_parser.add_argument(
        "--all",
        action="store_true",
        default=False,
        help="Output the results as all of the above",
    )

    # create the parser for the "prepare" command
    prepare_parser = subparsers.add_parser(
        "prepare", help="run RiboMetric in preparation mode"
    )
    prepare_parser.add_argument(
        "-g", "--gff", type=str, required=True, help="Path to gff file"
    )
    prepare_parser.add_argument(
        "-T",
        "--transcripts",
        type=int,
        required=False,
        default=10000000000,
        help="Number of transcripts to consider (default: 100000)",
    )
    prepare_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default=".",
        help="""Path to the output directory
            (default: current directory)""",
    )
    prepare_parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=False,
        default="config.yml",
        help="Path to the config file (default: config.yml)",
    )
    return parser


def main(args):
    """
    Main function for the RiboMetric command line interface

    Inputs:
        args: Namespace object containing the parsed arguments

    Outputs:
        None
    """
    console = Console()
    print_logo(console)

    if os.path.exists(args.config):
        with open(args.config, "r") as yml:
            config = yaml.load(yml, Loader=yaml.Loader)
    else:
        # load default config file
        project_dir = os.path.dirname(os.path.abspath(__file__))
        config_file_path = os.path.join(project_dir, 'config.yml')

        with open(config_file_path, "r") as yml:
            config = yaml.load(yml, Loader=yaml.Loader)

    if args.command == "prepare":
        print_table_prepare(args, console, "Prepare Mode")
        prepare_annotation(args.gff, args.output, args.transcripts, config)

    else:
        if args.all:
            args.json = True
            args.html = True
            args.pdf = True
            args.csv = True
        print_table_run(args, console, "Run Mode")

        if args.html:
            if args.pdf:
                report_export = "both"
            else:
                report_export = "html"
        elif args.pdf:
            report_export = "pdf"

        if not check_bam(args.bam):
            raise Exception("""
            Either BAM file or it's index does not exist at given path

            To create an index for a BAM file, run:
            samtools index <bam_file>
            """)
        flagstat = flagstat_bam(args.bam)
        if flagstat['mapped_reads'] < args.subsample:
            read_limit = flagstat['mapped_reads']
        else:
            read_limit = args.subsample

        read_df_pre = pd.concat(parse_bam(
            args.bam,
            read_limit
            ), ignore_index=True)

        read_df_pre["reference_name"] = (read_df_pre["reference_name"]
                                         .astype("category"))

        # sequence_data = {}
        # for item in sequence_list:
        #     for key, value in item.items():
        #         if key in sequence_data:
        #             sequence_data[key] += value
        #         else:
        #             sequence_data[key] = value

        # del sequence_list
        # print(sequence_data) # temp
        # print(read_df_pre.head()) # temp
        print("Reads parsed")

        # Expand the dataframe to have one row per read
        if "count" not in read_df_pre.columns:
            read_df_pre["count"] = 1
            read_df = read_df_pre
        else:
            print("Expanding dataframe")
            repeat_indices = np.repeat(read_df_pre.index, read_df_pre["count"])
            read_df = read_df_pre.iloc[repeat_indices].reset_index(drop=True)
            print("Dataframe expanded")

        del read_df_pre
        print("Calculating A site information")
        read_df = a_site_calculation(read_df)

        if args.gff is None and args.annotation is None:
            results_dict = annotation_mode(read_df,
                                           config=config)

        else:
            if args.annotation is not None and args.gff is not None:
                print("Both annotation and gff provided, using annotation")
                annotation_df = parse_annotation(args.annotation)
            elif args.annotation is None and args.gff is not None:
                print("Gff provided, preparing annotation")
                annotation_df = prepare_annotation(
                    args.gff, args.output, args.transcripts, config
                )
                print("Annotation prepared")

            elif args.annotation is not None and args.gff is None:
                print("Annotation provided, parsing")
                annotation_df = parse_annotation(args.annotation)
                print("Annotation parsed")

                print("Running annotation mode")
                results_dict = annotation_mode(read_df, annotation_df, config)

            if args.fasta is not None:
                fasta_dict = parse_fasta(args.fasta)
                results_dict = sequence_mode(
                    results_dict, read_df, fasta_dict, config
                )

        filename = args.bam.split('/')[-1]
        if "." in filename:
            filename = filename.split('.')[:-1]
        report_prefix = f"{''.join(filename)}_RiboMetric"

        if args.html or args.pdf:
            plots_list = generate_plots(results_dict, config)
            generate_report(plots_list,
                            config,
                            report_export,
                            report_prefix,
                            args.output)

        if args.json:
            generate_json(results_dict,
                          config,
                          report_prefix,
                          args.output)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    main(args)
