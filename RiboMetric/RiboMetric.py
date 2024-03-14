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

from rich.console import Console
from rich.text import Text
from rich.table import Table

import numpy as np

from .file_parser import (
    parse_bam,
    parse_fasta,
    parse_annotation,
    prepare_annotation,
    flagstat_bam,
    check_bam,
    check_annotation
)
from .arg_parser import argument_parser, open_config
from .qc import annotation_mode, sequence_mode
from .plots import generate_plots
from .modules import a_site_calculation
from .html_report import generate_report, parse_json_input
from .results_output import generate_json, generate_csv


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
              ╚═╝  ╚═╝ ╚═╝ ╚═════╝  ╚═════╝
    """,
        style="bold blue",
    )
    logo += Text(
        """
    ███╗   ███╗███████╗█████████╗██████╗  ██╗  ██████╗
    ████╗ ████║██╔════╝╚══██╔═══╝██╔══██╗ ██║ ██╔════╝
    ██╔████╔██║█████╗     ██║    ██████╔╝ ██║ ██║
    ██║╚██╔╝██║██╔══╝     ██║    ██╔══██╗ ██║ ██║
    ██║ ╚═╝ ██║███████╗   ██║    ██║  ██║ ██║ ╚██████╗
    ╚═╝     ╚═╝╚══════╝   ╚═╝    ╚═╝  ╚═╝ ╚═╝  ╚═════╝
    """,
        style="bold red",
    )
    console.print(logo)


def print_table_run(args, config: dict, console, mode):
    console = Console()

    Inputs = Table(show_header=True, header_style="bold magenta")
    Inputs.add_column("Parameters", style="dim", width=20)
    Inputs.add_column("Values")
    if config["argument"]["bam"]:
        Inputs.add_row("Bam File:", config["argument"]["bam"])
        if config["argument"]["annotation"]:
            Inputs.add_row(
                "Annotation File:", config["argument"]["annotation"]
                )
        else:
            Inputs.add_row("GFF File:", config["argument"]["gff"])
        Inputs.add_row("Transcriptome File:", config["argument"]["fasta"])
    elif config["argument"]["json_in"]:
        Inputs.add_row("JSON File:", config["argument"]["json_in"])

    Configs = Table(show_header=True, header_style="bold yellow")
    Configs.add_column("Options", style="dim", width=20)
    Configs.add_column("Values")
    Configs.add_row("Mode:", mode)
    Configs.add_row("# of reads:", str(config["argument"]["subsample"]
                                       if not None else "Full file"))
    Configs.add_row("# of transcripts:", str(config["argument"]["transcripts"]
                                             if not None else "Full file"))
    Configs.add_row("# of threads:", str(config["argument"]["threads"]))
    Configs.add_row("Config file:", args.config)

    Output = Table(show_header=True, header_style="bold blue")
    Output.add_column("Output Options", style="dim", width=20)
    Output.add_column("Values")
    Output.add_row("JSON:", str(config["argument"]["json"]))
    Output.add_row("HTML:", str(config["argument"]["html"]))
    Output.add_row("PDF:", str(config["argument"]["pdf"]))
    Output.add_row("CSV:", str(config["argument"]["csv"]))

    # Print tables side by side
    console.print(Inputs, Configs, Output, justify="inline", style="bold")


def print_table_prepare(args, config, console, mode):
    console = Console()

    Inputs = Table(show_header=True, header_style="bold magenta")
    Inputs.add_column("Parameters", style="dim", width=20)
    Inputs.add_column("Values")
    Inputs.add_row("Gff File:", config["argument"]["gff"])

    Configs = Table(show_header=True, header_style="bold yellow")
    Configs.add_column("Options", style="dim", width=20)
    Configs.add_column("Values")
    Configs.add_row("Mode:", mode)
    Configs.add_row("# of transcripts:", str(config["argument"]["transcripts"]
                                             if not None else "Full file"))
    Configs.add_row("# of threads:", str(config["argument"]["threads"]))
    Configs.add_row("Config file:", args.config)

    # Print tables side by side
    console.print(Inputs, Configs, justify="inline", style="bold")


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

    config = open_config(args)
    export = config["argument"].copy()

    # Handle inputs and run modes appropriately
    if args.command == "prepare":
        print_table_prepare(args, config, console, "Prepare Mode")
        prepare_annotation(config["argument"]["gff"],
                           config["argument"]["output"],
                           config["argument"]["transcripts"],
                           config["argument"]["threads"],
                           )

    else:
        print_table_run(args, config, console, "Run Mode")

        if config["argument"]["bam"]:
            if not check_bam(config["argument"]["bam"]):
                raise Exception("""
                Either BAM file or it's index does not exist at given path

                To create an index for a BAM file, run:
                samtools index <bam_file>
                """)

            if config["argument"]["annotation"] is not None:
                if not check_annotation(config["argument"]["annotation"]):
                    raise Exception("""
                    Annotation file not found or not in the correct format.

                    To create an annotation file, run:
                    RiboMetric prepare -g <gff_file>
                    """)

            # To Do:
            # Ensure Flagstat is written to the report
            flagstat = flagstat_bam(config["argument"]["bam"])
            if (config["argument"]["subsample"] is None
                    or flagstat['mapped_reads'] < config[
                        "argument"]["subsample"]):
                read_limit = flagstat['mapped_reads']
            else:
                read_limit = config["argument"]["subsample"]

            # Parse the bam file
            read_df_pre, sequence_data, sequence_background = parse_bam(
                bam_file=config["argument"]["bam"],
                num_reads=read_limit,
                num_processes=config["argument"]["threads"],
                server_mode=config["argument"]["server"])
            if read_df_pre.empty:
                raise Exception("""
                No reads found in the given bam file.

                Please check the file and try again.
                """)
            print("Reads parsed")

            # Expand the dataframe to have one row per read
            # This makes calculations on reads easier but uses more memory
            if "count" not in read_df_pre.columns:
                read_df_pre["count"] = 1
                read_df = read_df_pre
            else:
                print("Expanding dataframe")
                repeat_indices = np.repeat(read_df_pre.index,
                                           read_df_pre["count"])
                read_df = (read_df_pre.iloc[repeat_indices]
                           .reset_index(drop=True))
                print("Dataframe expanded")

            del read_df_pre
            print("Calculating A site information...")

            read_df = a_site_calculation(read_df)
            if (config["argument"]["gff"] is None and
                    config["argument"]["annotation"] is None):
                results_dict = annotation_mode(read_df,
                                               sequence_data,
                                               sequence_background,
                                               config=config)

            else:
                if (config["argument"]["annotation"] is not None and
                        config["argument"]["gff"] is not None):
                    print("Running annotation mode")
                    annotation_df = parse_annotation(
                        config["argument"]["annotation"]
                        )
                elif (config["argument"]["annotation"] is None and
                        config["argument"]["gff"] is not None):
                    print("Gff provided, preparing annotation")
                    annotation_df = prepare_annotation(
                        config["argument"]["gff"],
                        config["argument"]["output"],
                        config["argument"]["transcripts"],
                        config
                    )
                    print("Annotation prepared")

                elif (config["argument"]["annotation"] is not None and
                        config["argument"]["gff"] is None):
                    print("Annotation provided, parsing")
                    annotation_df = parse_annotation(
                        config["argument"]["annotation"]
                        )
                    print("Annotation parsed")

                    print("Running annotation mode")
                    results_dict = annotation_mode(read_df,
                                                   sequence_data,
                                                   sequence_background,
                                                   annotation_df,
                                                   config)
                if config["argument"]["fasta"] is not None:
                    fasta_dict = parse_fasta(config["argument"]["fasta"])
                    results_dict = sequence_mode(
                        results_dict, read_df, fasta_dict, config
                    )

            filename = config["argument"]["bam"].split('/')[-1]
            if "." in filename:
                filename = filename.split('.')[:-1]

        elif config["argument"]["json_in"]:
            print("JSON input provided")
            filename = config["argument"]["json_in"].split('/')[-1]
            if "." in filename:
                filename = filename.split('.')[:-1]

            json_dicts = parse_json_input(config["argument"]["json_in"])
            results_dict = json_dicts[0]
            json_config = json_dicts[1]
            config["argument"] = json_config["argument"]

            if config["argument"]["json_config"]:
                config["plots"] = json_config["plots"]

        # Indentify output requirements
        if export["name"] is not None:
            filename = export["name"]

        report_prefix = f"{''.join(filename)}_RiboMetric"

        if export["html"]:
            if export["pdf"]:
                report_export = "both"
            else:
                report_export = "html"
        elif export["pdf"]:
            report_export = "pdf"
        else:
            report_export = None
        # Write out the specified output files
        if report_export is not None:
            plots_list = generate_plots(results_dict, config)
            generate_report(plots_list,
                            config,
                            report_export,
                            report_prefix,
                            export["output"])

        if export["json"]:
            generate_json(results_dict,
                          config,
                          report_prefix,
                          export["output"])

        if export["csv"]:
            generate_csv(results_dict,
                         config,
                         report_prefix,
                         export["output"])


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    if not vars(args):
        parser.print_help()

    if args.bam and args.json_in:
        parser.error(
            "Only one of -b/--bam or -j/--json-in should be specified.")
    elif not args.bam and not args.json_in:
        parser.error("Either -b/--bam or -j/--json-in must be specified.")

    main(args)
