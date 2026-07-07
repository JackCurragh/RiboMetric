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
from typing import Dict, Any
import argparse
import hashlib
import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

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
from .qc import annotation_mode
from .plots import generate_plots
from .html_report import generate_report, parse_json_input
from .results_output import (
    generate_json,
    generate_csv,
    generate_summary_tsv,
    generate_qc_status,
    generate_comparison_ready_csv,
    generate_metrics_table_csv,
    generate_offsets_tsv,
    generate_all_outputs,
)


_HASH_SAMPLE_BYTES = 1024 * 1024
_FULL_HASH_LIMIT_BYTES = 50 * 1024 * 1024


def _config_path_used(args: argparse.Namespace) -> str:
    if os.path.exists(args.config):
        return args.config
    return str(Path(__file__).with_name("config.yml"))


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _file_fingerprint(path_value: Any) -> Dict[str, Any]:
    """Return a reproducibility fingerprint without forcing huge full-file hashes."""
    if not path_value:
        return {"path": path_value, "exists": False}

    path = Path(str(path_value))
    record: Dict[str, Any] = {
        "path": str(path),
        "exists": path.exists(),
    }
    if not path.exists() or not path.is_file():
        return record

    stat = path.stat()
    record.update({
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    })

    full_hash = (
        os.environ.get("RIBOMETRIC_FULL_INPUT_HASH", "").lower()
        in {"1", "true", "yes"}
        or stat.st_size <= _FULL_HASH_LIMIT_BYTES
    )

    h = hashlib.sha256()
    if full_hash:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(chunk)
        record["sha256"] = h.hexdigest()
        record["hash_method"] = "full_sha256"
    else:
        with path.open("rb") as handle:
            first = handle.read(_HASH_SAMPLE_BYTES)
            if stat.st_size > _HASH_SAMPLE_BYTES:
                handle.seek(max(0, stat.st_size - _HASH_SAMPLE_BYTES))
                last = handle.read(_HASH_SAMPLE_BYTES)
            else:
                last = b""
        h.update(first)
        h.update(last)
        record["sha256_sampled"] = h.hexdigest()
        record["hash_method"] = (
            f"first_last_{_HASH_SAMPLE_BYTES}_bytes_sha256"
        )
    return record


def _build_run_provenance(args: argparse.Namespace, config: Dict[str, Any]) -> Dict[str, Any]:
    arg_cfg = config.get("argument", {})
    input_keys = [
        "bam",
        "annotation",
        "gff",
        "fasta",
        "json_in",
        "offset_read_length",
        "offset_read_specific",
    ]
    config_path = _config_path_used(args)
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "command": getattr(args, "command", None),
        "package_version": _package_version(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "config_file": _file_fingerprint(config_path),
        "effective_config_sha256": _sha256_bytes(
            json.dumps(config, sort_keys=True, default=str).encode("utf-8")
        ),
        "inputs": {
            key: _file_fingerprint(arg_cfg.get(key))
            for key in input_keys
            if arg_cfg.get(key) is not None
        },
    }


def _package_version() -> str:
    try:
        from . import __version__
        return str(__version__)
    except Exception:
        return "unknown"


def print_logo(console: Console) -> None:
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


def print_table_run(args: argparse.Namespace, config: Dict[str, Any], console: Console, mode: str) -> None:
    console = Console()

    Inputs = Table(show_header=True, header_style="bold magenta")
    Inputs.add_column("Parameters", style="dim", width=20)
    Inputs.add_column("Values")
    if config["argument"]["bam"]:
        Inputs.add_row("Bam File:", config["argument"]["bam"])
        if config["argument"].get("annotation"):
            Inputs.add_row("Annotation File:", config["argument"]["annotation"])
        elif config["argument"].get("gff"):
            Inputs.add_row("GFF File:", config["argument"]["gff"])
        if config["argument"].get("fasta"):
            Inputs.add_row("Transcriptome File:", config["argument"]["fasta"])
    elif config["argument"]["json_in"]:
        Inputs.add_row("JSON File:", config["argument"]["json_in"])

    Configs = Table(show_header=True, header_style="bold yellow")
    Configs.add_column("Options", style="dim", width=20)
    Configs.add_column("Values")
    Configs.add_row("Mode:", mode)
    subs = config["argument"].get("subsample")
    trans = config["argument"].get("transcripts")
    Configs.add_row("# of reads:", str(subs) if subs is not None else "Full file")
    Configs.add_row("# of transcripts:", str(trans) if trans is not None else "Full file")
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
    console.print(Inputs, Configs, Output, justify=None, style="bold")


def print_table_prepare(args: argparse.Namespace, config: Dict[str, Any], console: Console, mode: str) -> None:
    console = Console()

    Inputs = Table(show_header=True, header_style="bold magenta")
    Inputs.add_column("Parameters", style="dim", width=20)
    Inputs.add_column("Values")
    Inputs.add_row("Gff File:", config["argument"]["gff"])

    Configs = Table(show_header=True, header_style="bold yellow")
    Configs.add_column("Options", style="dim", width=20)
    Configs.add_column("Values")
    Configs.add_row("Mode:", mode)
    trans = config["argument"].get("transcripts")
    Configs.add_row("# of transcripts:", str(trans) if trans is not None else "Full file")
    Configs.add_row("# of threads:", str(config["argument"]["threads"]))
    Configs.add_row("Config file:", args.config)

    # Print tables side by side
    console.print(Inputs, Configs, justify=None, style="bold")


def main(args: argparse.Namespace) -> int:
    """
    Main function for the RiboMetric command line interface

    Inputs:
        args: Namespace object containing the parsed arguments

    Outputs:
        None
    """
    # Handle evaluate command separately (no logo or config needed)
    if args.command == "evaluate":
        from .evaluate import evaluate as run_evaluate
        return run_evaluate(args)

    # Handle view command separately (no logo or config needed)
    if args.command == "view":
        from .tui import run_tui
        import json
        from pathlib import Path

        # Validate file exists and is JSON
        file_path = Path(args.json_file)
        if not file_path.exists():
            print(f"Error: File not found: {args.json_file}")
            return 1

        if not file_path.suffix == '.json':
            print("Error: File must be a JSON file")
            return 1

        # Validate it's a RiboMetric JSON file
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                if "results" not in data or "config" not in data:
                    print("Warning: File may not be a valid RiboMetric JSON file.")
                    print("Expected structure: {'results': {...}, 'config': {...}}")
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON file: {e}")
            return 1
        except Exception as e:
            print(f"Error reading file: {e}")
            return 1

        # Launch TUI
        run_tui(str(file_path))
        return 0

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
            )
            if read_df_pre.empty:
                raise Exception("""
                No reads found in the given bam file.

                Please check the file and try again.
                """)
            print("Reads parsed")

            # Optionally drop sequence data so sequence-based metrics are
            # skipped (explicit flag, or BAMs with no stored sequences). #122
            if config["argument"].get("skip_sequence_metrics"):
                print("Skipping sequence-based metrics (--skip-sequence-metrics)")
                sequence_data, sequence_background = {}, {}
            elif not sequence_background:
                print("No stored sequences detected; "
                      "sequence-based metrics will be skipped.")

            # Use weighted computations downstream instead of expanding rows
            if "count" not in read_df_pre.columns:
                read_df_pre["count"] = 1
            read_df = read_df_pre
            del read_df_pre

            # Parse FASTA up-front so it can be passed into annotation_mode
            # for RUST and other sequence-level metrics.
            fasta_dict = None
            if config["argument"]["fasta"] is not None:
                print("Parsing FASTA for sequence-level metrics...")
                fasta_dict = parse_fasta(config["argument"]["fasta"])

            if (config["argument"]["gff"] is None and
                    config["argument"]["annotation"] is None):
                results_dict = annotation_mode(
                    read_df,
                    sequence_data,
                    sequence_background,
                    config=config,
                    fasta_dict=fasta_dict,
                )

            else:
                if (config["argument"]["annotation"] is not None and
                        config["argument"]["gff"] is not None):
                    print("Running annotation mode")
                    annotation_df = parse_annotation(
                        config["argument"]["annotation"]
                        )
                    # Ensure annotation mode actually runs for this branch
                    results_dict = annotation_mode(
                        read_df,
                        sequence_data,
                        sequence_background,
                        annotation_df,
                        config,
                        fasta_dict=fasta_dict,
                    )
                elif (config["argument"]["annotation"] is None and
                        config["argument"]["gff"] is not None):
                    print("Gff provided, preparing annotation")
                    annotation_df = prepare_annotation(
                        config["argument"]["gff"],
                        config["argument"]["output"],
                        config["argument"]["transcripts"],
                        config["argument"]["threads"]
                    )
                    print("Annotation prepared")
                    # Run annotation mode after preparing annotation
                    results_dict = annotation_mode(
                        read_df,
                        sequence_data,
                        sequence_background,
                        annotation_df,
                        config,
                        fasta_dict=fasta_dict,
                    )

                elif (config["argument"]["annotation"] is not None and
                        config["argument"]["gff"] is None):
                    print("Annotation provided, parsing")
                    annotation_df = parse_annotation(
                        config["argument"]["annotation"]
                        )
                    print("Annotation parsed")

                    print("Running annotation mode")
                    results_dict = annotation_mode(
                        read_df,
                        sequence_data,
                        sequence_background,
                        annotation_df,
                        config,
                        fasta_dict=fasta_dict,
                    )

            # Merge samtools flagstat data into alignment_stats so the report
            # can display total_reads, mapping_rate, etc.
            results_dict.setdefault("alignment_stats", {}).update(flagstat)
            results_dict["provenance"] = _build_run_provenance(args, config)
    
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
            results_dict["provenance"] = _build_run_provenance(args, config)

        # Indentify output requirements
        if export["name"] is not None:
            filename = export["name"]

        report_prefix = f"{''.join(filename)}_RiboMetric"

        if export.get("html"):
            if export["pdf"]:
                report_export = "both"
            else:
                report_export = "html"
        elif export.get("pdf"):
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

        if export.get("json"):
            generate_json(results_dict,
                          config,
                          report_prefix,
                          export["output"])

        if export.get("csv"):
            generate_csv(results_dict,
                         config,
                         report_prefix,
                         export["output"])

        # Improved outputs (low-hanging fruit)
        sample_name = "".join(filename)
        if export.get("improved_outputs"):
            generate_all_outputs(
                results_dict,
                config,
                sample_name,
                export.get("output", ""),
            )
        else:
            if export.get("summary_tsv"):
                generate_summary_tsv(
                    results_dict, config, sample_name,
                    f"{sample_name}_summary.tsv", export.get("output", "")
                )
            if export.get("metrics_table"):
                generate_metrics_table_csv(
                    results_dict, config, sample_name,
                    f"{sample_name}_metrics_table.csv", export.get("output", "")
                )
            if export.get("qc_status"):
                generate_qc_status(
                    results_dict, config, sample_name, None,
                    f"{sample_name}_qc_status.json", export.get("output", "")
                )
            if export.get("comparison_csv"):
                generate_comparison_ready_csv(
                    results_dict, config, sample_name,
                    f"{sample_name}_comparison.csv", export.get("output", "")
                )
            if export.get("offsets_tsv"):
                generate_offsets_tsv(
                    results_dict,
                    sample_name,
                    f"{sample_name}_offsets.tsv",
                    export.get("output", ""),
                )

        if export.get("output_offsets"):
            output_path = Path(export["output_offsets"])
            generate_offsets_tsv(
                results_dict,
                sample_name,
                output_path.name,
                str(output_path.parent) if str(output_path.parent) != "." else "",
            )


    return 0
if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    # -b/--bam vs -j/--json-in exclusivity is enforced by argparse's
    # mutually-exclusive group on the `run` subparser; no manual check needed
    # here (the previous one assumed those attributes existed for every
    # subcommand and crashed for prepare/evaluate/view).
    if getattr(args, "command", None) is None:
        parser.print_help()
        raise SystemExit(0)

    raise SystemExit(main(args) or 0)
