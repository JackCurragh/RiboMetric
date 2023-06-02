
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
        default=False,
        help="Output the results as an html file (default)",
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

    for arg in vars(args):
        print(f"{arg}: {getattr(args, arg)}")

    if args.command == "run" and args.all:
        args.json = True
        args.html = True
        args.pdf = True
        args.csv = True

    for arg in vars(args):
        if getattr(args, arg) is not False:
            config["argument"][arg] = getattr(args, arg)

    return config
