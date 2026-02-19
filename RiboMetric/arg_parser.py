
import os
import yaml
import argparse
from rich.emoji import Emoji


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
        "run",
        help="run RiboMetric in normal mode"
    )
    group_run_parser = run_parser.add_mutually_exclusive_group(required=True)
    group_run_parser.add_argument(
        "-b",
        "--bam",
        type=str,
        help="Path to bam file"
    )
    group_run_parser.add_argument(
        "-j",
        "--json-in",
        type=str,
        help="Path to json input file"
    )
    run_parser.add_argument(
        "-a",
        "--annotation",
        type=str,
        required=False,
        help="Path to RiboMetric annotation file",
    )
    run_parser.add_argument(
        "-g",
        "--gff",
        type=str,
        required=False,
        help="Path to gff file"
    )
    run_parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=False,
        help="Path to the transcriptome fasta file",
    )
    run_parser.add_argument(
        "--offset-read-length",
        type=str,
        required=False,
        help="Path to the tsv file of read length specific offsets format: read_length <tab> offset",
    )
    run_parser.add_argument(
        "--offset-read-specific",
        type=str,
        required=False,
        help="Path to the tsv file of read specific offsets format: read_name <tab> offset",
    )
    run_parser.add_argument(
        "--offset-global",
        type=int,
        required=False,
        default=15,
        help="Global offset to be used for all read lengths (default: 15)",
    )
    run_parser.add_argument(
        "--offset-calculation-method",
        type=str,
        required=False,
        default="changepoint",
        help="Method to calculate offsets (default: changepoint) [changepoint, ribowaltz]",
    )
    run_parser.add_argument(
        "--enable-optional-metrics",
        action="store_true",
        default=False,
        help="Enable all optional (theoretical) metrics in addition to defaults",
    )
    run_parser.add_argument(
        "--enable-metric",
        type=str,
        action="append",
        help="Enable specific optional metric(s) (can be used multiple times)",
    )
    run_parser.add_argument(
        "--json-config",
        action="store_true",
        default=False,
        help="Use JSON config instead of active config for generating plots",
    )
    run_parser.add_argument(
        "-n",
        "--name",
        type=str,
        required=False,
        help="""Name of the sample being analysed for output files
            (default: filename of bam file)""",
    )
    run_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        help="Path to the output directory",
    )
    run_parser.add_argument(
        "-S",
        "--subsample",
        type=int,
        required=False,
        help="Number of reads to subsample from the bam file",
    )
    run_parser.add_argument(
        "-T",
        "--transcripts",
        type=int,
        required=False,
        help="Number of transcripts to consider",
    )
    run_parser.add_argument(
        "-p",
        "--threads",
        type=int,
        required=False,
        help="Number of threads used by RiboMetric",
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
        default=False,
        help="Output the results as an html file (default)",
    )
    run_parser.add_argument(
        "--pdf",
        action="store_true",
        default=False,
        help="Output the results as a pdf file",
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
    run_parser.add_argument(
        "--output-offsets",
        type=str,
        required=False,
        help="Path to write per-read-length offsets TSV (read_len<tab>offset). "
             "Only written when offsets are calculated internally (not from --offset-read-length).",
    )

    # create the parser for the "prepare" command
    prepare_parser = subparsers.add_parser(
        "prepare", help="run RiboMetric in preparation mode"
    )
    prepare_parser.add_argument(
        "-g", "--gff", type=str, required=True, help="Path to gff file"
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
        "-T",
        "--transcripts",
        type=int,
        required=False,
        help="Number of transcripts to consider",
    )
    prepare_parser.add_argument(
        "-p",
        "--threads",
        type=int,
        required=False,
        help="Number of threads used by RiboMetric",
    )
    prepare_parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=False,
        default="config.yml",
        help="Path to the config file (default: config.yml)",
    )

    # create the parser for the "view" command
    view_parser = subparsers.add_parser(
        "view",
        help="view RiboMetric results interactively in the terminal"
    )
    view_parser.add_argument(
        "json_file",
        type=str,
        help="Path to RiboMetric JSON output file (*_RiboMetric_data.json)"
    )

    return parser


def open_config(args) -> dict:
    """
    Opens config and overrides config dictionary with commandline arguments

    Inputs:
        args: Arguments passed from the commandline through argparse
        config: Config read from yaml file

    Outputs:
        config: Modified config with arguments
    """
    if os.path.exists(args.config):
        with open(args.config, "r") as yml:
            config = yaml.load(yml, Loader=yaml.Loader)
    else:
        # load default config file
        project_dir = os.path.dirname(os.path.abspath(__file__))
        config_file_path = os.path.join(project_dir, 'config.yml')

        with open(config_file_path, "r") as yml:
            config = yaml.load(yml, Loader=yaml.Loader)

    if args.command == "run" and args.all:
        args.json = True
        args.html = True
        args.pdf = True
        args.csv = True

    for arg in vars(args):
        if getattr(args, arg) is not False and getattr(args, arg) is not None:
            config["argument"][arg] = getattr(args, arg)

    # CI/Env override: allow RIBOMETRIC_THREADS to set threads
    try:
        env_threads = os.environ.get("RIBOMETRIC_THREADS")
        if env_threads:
            config["argument"]["threads"] = int(env_threads)
    except Exception:
        pass

    # Handle metric selection
    if args.command == "run":
        enabled_metrics = set(config.get("metrics", {}).get("default", []))

        if hasattr(args, 'enable_optional_metrics') and args.enable_optional_metrics:
            # Add all optional metrics
            enabled_metrics.update(config.get("metrics", {}).get("optional", []))

        if hasattr(args, 'enable_metric') and args.enable_metric:
            # Add specific metrics
            enabled_metrics.update(args.enable_metric)

        config["enabled_metrics"] = list(enabled_metrics)

    return config
