
import argparse
import os
from typing import Any, Dict, cast

import yaml
from rich.emoji import Emoji


def argument_parser() -> argparse.ArgumentParser:
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
        default=None,
        dest="global_offset",
        help="Apply this fixed offset to every read length instead of "
             "calculating per-read-length offsets. If omitted, offsets are "
             "calculated automatically.",
    )
    run_parser.add_argument(
        "--offset-min",
        type=int,
        required=False,
        default=None,
        help="Minimum plausible per-read-length offset (default: config.yml, 8)",
    )
    run_parser.add_argument(
        "--offset-max",
        type=int,
        required=False,
        default=None,
        help="Maximum plausible per-read-length offset (default: config.yml, 20)",
    )
    run_parser.add_argument(
        "--offset-max-read-length-fraction",
        type=float,
        required=False,
        default=None,
        help=(
            "Maximum plausible offset as a fraction of read length "
            "(default: config.yml, 0.6667)"
        ),
    )
    run_parser.add_argument(
        "--offset-calculation-method",
        type=str,
        required=False,
        default=None,  # None so config.yml controls the default
        choices=["changepoint", "ribowaltz", "tripsviz"],
        help=(
            "Method to calculate offsets. If omitted, the value in the active "
            "config.yml is used."
        ),
    )
    run_parser.add_argument(
        "--offset-target",
        type=str,
        required=False,
        default=None,
        choices=["a_site", "p_site"],
        help=(
            "Site represented by calculated offsets. Automatic methods infer "
            "a P-site peak and convert to this target. Default: config.yml."
        ),
    )
    run_parser.add_argument(
        "--min-periodicity",
        type=float,
        required=False,
        default=None,  # None so config.yml controls the default
        help=(
            "Minimum dominant-frame fraction for a read length to be "
            "recommended for downstream P-site/ORF analysis "
            "(default: config.yml, 0.5)."
        ),
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
    # Improved outputs
    run_parser.add_argument(
        "--improved-outputs",
        action="store_true",
        default=False,
        help=(
            "Emit all improved outputs: summary TSV, QC status JSON, "
            "comparison CSV, and detailed metrics table"
        ),
    )
    run_parser.add_argument(
        "--summary-tsv",
        action="store_true",
        default=False,
        help="Write {sample}_summary.tsv (one-line summary)",
    )
    run_parser.add_argument(
        "--qc-status",
        action="store_true",
        default=False,
        help="Write {sample}_qc_status.json (pass/warn/fail)",
    )
    run_parser.add_argument(
        "--comparison-csv",
        action="store_true",
        default=False,
        help="Append {sample}_comparison.csv (wide sample comparison)",
    )
    run_parser.add_argument(
        "--metrics-table",
        action="store_true",
        default=False,
        help="Write {sample}_metrics_table.csv (long-form metrics)",
    )
    run_parser.add_argument(
        "--output-offsets",
        type=str,
        required=False,
        help="Path to write per-read-length offsets TSV (read_len<tab>offset). "
             "Only written when offsets are calculated internally (not from --offset-read-length).",
    )
    run_parser.add_argument(
        "--skip-sequence-metrics",
        action="store_true",
        default=False,
        dest="skip_sequence_metrics",
        help=(
            "Skip all sequence-based metrics (terminal nucleotide bias, "
            "nucleotide composition). Use for BAMs with no stored sequences, "
            "or to speed up runs where sequence metrics are not needed."
        ),
    )
    run_parser.add_argument(
        "--multimap-filter",
        type=str,
        required=False,
        default=None,  # None so config.yml controls the default
        choices=["unique_only", "none"],
        dest="multimap_filter",
        help=(
            "How to handle multimapping reads for frame-sensitive calculations "
            "(offset detection and periodicity). "
            "'unique_only' restricts those calculations to MAPQ=255 reads, which "
            "in STAR output means exactly one genomic mapping locus. "
            "'none' uses all primary reads (pre-0.2.0 behaviour). "
            "If omitted, the value in config.yml is used (default: unique_only)."
        ),
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

    # create the parser for the "evaluate" command
    evaluate_parser = subparsers.add_parser(
        "evaluate",
        help="assess whether a result passes QC against expected thresholds",
    )
    evaluate_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to a RiboMetric result (JSON output or metrics-table CSV)",
    )
    evaluate_parser.add_argument(
        "-e",
        "--expected",
        type=str,
        required=False,
        help=(
            "Path to a YAML of expected thresholds "
            "(metric -> {pass, warn}). If omitted, built-in defaults are used."
        ),
    )
    evaluate_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        help="Optional path to write the evaluation result as JSON",
    )
    evaluate_parser.add_argument(
        "-n",
        "--name",
        type=str,
        required=False,
        help="Sample name for the report (default: input filename stem)",
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


def open_config(args: argparse.Namespace) -> Dict[str, Any]:
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
            config = cast(Dict[str, Any], yaml.load(yml, Loader=yaml.Loader))
    else:
        # load default config file
        project_dir = os.path.dirname(os.path.abspath(__file__))
        config_file_path = os.path.join(project_dir, 'config.yml')

        with open(config_file_path, "r") as yml:
            config = cast(Dict[str, Any], yaml.load(yml, Loader=yaml.Loader))

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
