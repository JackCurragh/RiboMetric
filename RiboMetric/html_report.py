"""
Code in this script is used to generate the HTML and pdf output report
The functions are called by the main script RiboMetric.py
if the user specifies the --html flag
"""

from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from .modules import convert_html_to_pdf
import base64
import json

import os


from typing import List, Dict, Any, Tuple


LOWER_IS_BETTER = {
    "duplicate_rate",
    "multimapper_rate",
    "rpf_multimapper_rate",
    "alignment_multimapper_rate",
    "soft_clip_rate_5prime",
    "disome_proportion",
    "terminal_bias_kl_5prime_raw",
    "terminal_bias_kl_3prime_raw",
    "stop_codon_readthrough_ratio",
    "marginal_position_discovery_rate",
    "floss_median",
    "floss_aberrant_transcript_fraction",
}

METRIC_GROUPS = {
    "Mapping": (
        "duplicate_rate",
        "multimapper_rate",
        "rpf_multimapper",
        "unique_rpf",
        "alignment_multimapper",
        "soft_clip_rate",
    ),
    "Footprint Lengths": (
        "read_length_distribution",
        "recommended_read",
        "disome",
    ),
    "Periodicity": (
        "periodicity",
        "uniformity",
    ),
    "Annotation / CDS": (
        "prop_reads",
        "cds_coverage",
        "CDS_coverage",
        "region",
        "ratio_",
        "start_codon",
        "stop_codon",
    ),
    "Codon / Sequence": (
        "terminal_",
        "terminal_bias",
        "rust",
        "codon",
    ),
    "Saturation": (
        "floss",
        "library",
        "marginal",
    ),
}

CARD_METRICS = (
    "periodicity_dominance",
    "prop_reads_CDS",
    "recommended_read_proportion",
    "duplicate_rate",
    "multimapper_rate",
)


METRIC_DESCRIPTIONS = {
    "duplicate_rate": "Estimated fraction of weighted reads explained by duplicate collapsed sequences.",
    "multimapper_rate": "Backwards-compatible weighted RPF-level fraction with evidence of multiple genomic loci.",
    "rpf_multimapper_rate": "Weighted fraction of ribosome protected fragments with evidence of multiple genomic loci.",
    "unique_rpf_rate": "Weighted fraction of ribosome protected fragments considered uniquely mapped.",
    "alignment_multimapper_rate": "Fraction of reported alignment rows whose fragment has evidence of another genomic alignment.",
    "soft_clip_rate_5prime": "Fraction of weighted reads with 5 prime soft clipping, often indicating trimming or adapter issues.",
    "read_length_distribution_IQR_metric": "Scores how tightly footprint lengths concentrate around the central distribution.",
    "read_length_distribution_coefficient_of_variation_metric": "Scores footprint length consistency using the coefficient of variation.",
    "read_length_distribution_maxprop_metric": "Fraction of reads in the most abundant footprint length.",
    "read_length_distribution_bimodality_metric": "Scores whether the footprint length distribution is dominated by one mode rather than multiple modes.",
    "read_length_distribution_normality_metric": "Scores departure from a normal length distribution, as sharp RPF peaks are expected.",
    "disome_proportion": "Fraction of reads in the configured di-ribosome footprint length window.",
    "terminal_nucleotide_bias_distribution_5_prime_metric": "Scores nucleotide balance at the 5 prime end of reads.",
    "terminal_nucleotide_bias_distribution_3_prime_metric": "Scores nucleotide balance at the 3 prime end of reads.",
    "terminal_bias_kl_5prime": "KL-divergence based score for 5 prime terminal nucleotide bias.",
    "terminal_bias_kl_3prime": "KL-divergence based score for 3 prime terminal nucleotide bias.",
    "terminal_bias_maxabs_5prime": "Largest absolute 5 prime terminal nucleotide deviation from expectation.",
    "terminal_bias_maxabs_3prime": "Largest absolute 3 prime terminal nucleotide deviation from expectation.",
    "periodicity_dominance": "Fraction of coding reads in the dominant reading frame after offset assignment.",
    "periodicity_autocorrelation": "Triplet-periodicity score from autocorrelation of frame-specific signal.",
    "periodicity_fourier": "Triplet-periodicity score from Fourier power at the codon frequency.",
    "periodicity_information": "Information-content score for separation between reading frames.",
    "periodicity_information_weighted_score": "Read-depth weighted information-content score for reading-frame separation.",
    "periodicity_trips-viz": "Trips-Viz style triplet-periodicity score.",
    "uniformity_entropy": "Entropy-based score for how evenly reads cover coding regions.",
    "uniformity_autocorrelation": "Autocorrelation-based score for smoothness of coding-region coverage.",
    "uniformity_theil_index": "Inequality-based score for coding-region coverage uniformity.",
    "uniformity_gini_index": "Gini-based score for coding-region coverage uniformity.",
    "CDS_coverage_metric": "Score summarising how much annotated CDS sequence is covered by reads.",
    "cds_coverage": "Observed proportion of annotated CDS positions covered by reads.",
    "prop_reads_CDS": "Fraction of annotated reads assigned to CDS regions.",
    "prop_reads_leader": "Fraction of annotated reads assigned to leader or 5 prime UTR regions.",
    "prop_reads_trailer": "Fraction of annotated reads assigned to trailer or 3 prime UTR regions.",
    "stop_codon_readthrough_ratio": "Ratio of downstream to upstream signal around stop codons.",
    "start_codon_enrichment_ratio": "Read enrichment near start codons relative to downstream coding background.",
    "recommended_read_proportion": "Fraction of reads in lengths recommended for downstream frame-sensitive analysis.",
    "n_recommended_read_lengths": "Number of read lengths passing the recommendation criteria.",
    "five_prime_ramp_ratio": "Relative 5 prime CDS ramp signal compared with downstream coding signal.",
    "three_prime_drop_ratio": "Relative 3 prime CDS signal compared with upstream coding signal.",
    "marginal_position_discovery_rate": "Estimated rate at which additional reads discover new covered positions.",
    "complexity_distinct_positions": "Number of distinct positions covered by the library.",
    "floss_median": "Median FLOSS-style footprint length heterogeneity across transcripts.",
    "floss_aberrant_transcript_fraction": "Fraction of transcripts with aberrant footprint length distributions.",
    "rust_mean_kl_divergence": "Mean codon-level RUST divergence across the metagene window.",
    "codon_dwell_cv": "Coefficient of variation in codon dwell signal.",
    "codon_dwell_p90_p10": "Spread between high and low codon dwell signals.",
    "proline_dwell": "Relative dwell signal on proline codons.",
    "cga_dwell": "Relative dwell signal on CGA codons.",
}


def _metric_key(metric: Dict[str, Any]) -> str:
    return str(metric.get("key") or metric.get("name", "").lower().replace(" ", "_"))


def _metric_label(metric_key: str) -> str:
    labels = {
        "periodicity_dominance": "Periodicity",
        "prop_reads_CDS": "CDS enrichment",
        "recommended_read_proportion": "Recommended reads",
        "duplicate_rate": "Duplicate rate",
        "multimapper_rate": "Multimapper rate",
        "rpf_multimapper_rate": "RPF multimapper rate",
        "unique_rpf_rate": "Unique RPF rate",
        "alignment_multimapper_rate": "Alignment multimapper rate",
        "soft_clip_rate_5prime": "5' soft clips",
    }
    return labels.get(metric_key, metric_key.replace("_", " ").capitalize())


def _metric_description(metric_key: str) -> str:
    if metric_key in METRIC_DESCRIPTIONS:
        return METRIC_DESCRIPTIONS[metric_key]
    for suffix in (
        "_global",
        "_rl28",
        "_rl29",
        "_rl30",
        "_rl31",
        "_rl32",
        "_metric",
    ):
        if metric_key.endswith(suffix):
            base_key = metric_key[: -len(suffix)]
            if base_key in METRIC_DESCRIPTIONS:
                return METRIC_DESCRIPTIONS[base_key]
    return "Quality-control metric included in the configured RiboMetric summary."


def _format_score(score: Any) -> str:
    if not isinstance(score, (float, int)):
        return str(score)
    if 0 <= score <= 1:
        return f"{score:.1%}"
    return f"{score:.3g}"


_RAW_BITS_KEYS = {"terminal_bias_kl_5prime_raw", "terminal_bias_kl_3prime_raw"}
_RAW_RATIO_KEYS = {"cds_enrichment_ratio", "start_codon_enrichment_ratio",
                   "stop_codon_readthrough_ratio"}


def _format_raw(key: str, raw: Any) -> str:
    """Format a raw value in its natural units for display beside the score."""
    if raw is None:
        return "n/a"
    if not isinstance(raw, (int, float)):
        return str(raw)
    if key in _RAW_BITS_KEYS:
        return f"{raw:.3f} bits"
    if key in _RAW_RATIO_KEYS:
        return f"E = {raw:.2f}"
    if 0 <= raw <= 1:
        return f"{raw:.1%}"
    return f"{raw:.3g}"


def _metric_status(metric_key: str, score: Any) -> str:
    if not isinstance(score, (float, int)):
        return "info"
    if metric_key in LOWER_IS_BETTER:
        if score <= 0.2:
            return "pass"
        if score <= 0.5:
            return "warn"
        return "fail"
    if score >= 0.7:
        return "pass"
    if score >= 0.5:
        return "warn"
    return "fail"


def _group_for_metric(metric_key: str) -> str:
    for group_name, prefixes in METRIC_GROUPS.items():
        if any(metric_key.startswith(prefix) for prefix in prefixes):
            return group_name
    return "Other"


# Map the resolver's status vocabulary to the short forms used as CSS classes.
_STATUS_SHORT = {
    "PASS": "pass",
    "WARNING": "warn",
    "FAIL": "fail",
    "INFO": "info",
}


def build_report_context(summary: Dict[str, Any]) -> Dict[str, Any]:
    """Build presentation-only report structures from the summary plot payload.

    Status and score come pre-resolved from the single scoring resolver
    (scoring.py) via the summary payload; this function no longer re-derives
    them, so badges, cards, and the QC gate cannot disagree. Older payloads
    without a precomputed status fall back to the local heuristic.
    """
    rows = []
    for metric in summary.get("metrics", []):
        key = _metric_key(metric)
        score = metric.get("score")
        raw = metric.get("raw")
        resolver_status = metric.get("status")
        if resolver_status in _STATUS_SHORT:
            status = _STATUS_SHORT[resolver_status]
        else:
            status = _metric_status(key, score)
        rows.append({
            "key": key,
            "label": _metric_label(key),
            "score": score,
            "score_label": _format_score(score),
            "raw": raw,
            "raw_label": _format_raw(key, raw),
            "status": status,
            "gate": bool(metric.get("gate", False)),
            "tier": int(metric["tier"]) if metric.get("tier") is not None else None,
            "direction": "Higher is better",
            "description": metric.get("decision") or _metric_description(key),
            "group": _group_for_metric(key),
        })

    # --- Tier-based grouping (new primary layout) ---
    def _tier_rows(tier_num: int) -> list:
        return [r for r in rows if r.get("tier") == tier_num]

    tier_sections = [
        {
            "tier": 1,
            "heading": "Tier 1 — Is this Ribo-seq?",
            "subtitle": "These metrics decide whether frame-dependent analysis is defensible. "
                        "Failures here mean P-site assignment and ORF calling should not proceed.",
            "metrics": _tier_rows(1),
        },
        {
            "tier": 2,
            "heading": "Tier 2 — Is it usable for my analysis?",
            "subtitle": "These metrics describe whether enough usable signal remains after filtering. "
                        "Weak scores are caveats, not identity failures.",
            "metrics": _tier_rows(2),
        },
        {
            "tier": 3,
            "heading": "Tier 3 — Technical caveats",
            "subtitle": "These describe technical distortions and loss of usable depth. "
                        "They colour interpretation but do not automatically fail the sample.",
            "metrics": _tier_rows(3),
        },
    ]
    tier_sections = [s for s in tier_sections if s["metrics"]]

    # --- Legacy grouped view (kept for backwards-compat; not shown in new layout) ---
    grouped_metrics = []
    for group_name in list(METRIC_GROUPS.keys()) + ["Other"]:
        group_rows = [row for row in rows if row["group"] == group_name]
        if group_rows:
            grouped_metrics.append({"name": group_name, "metrics": group_rows})

    metric_map = {row["key"]: row for row in rows}
    top_cards = [
        metric_map[key] for key in CARD_METRICS if key in metric_map
    ][:5]

    status_rank = {"fail": 3, "warn": 2, "pass": 1, "info": 0}
    gated_rows = [row for row in rows if row.get("gate")]
    verdict_rows = gated_rows if gated_rows else rows
    overall_status = "info"
    if verdict_rows:
        overall_status = max(
            verdict_rows, key=lambda row: status_rank[row["status"]]
        )["status"]

    # Three-way framing per METRICS_DESIGN.md §4
    periodicity = metric_map.get("periodicity_dominance")
    cds = metric_map.get("cds_enrichment_ratio") or metric_map.get("prop_reads_CDS")
    recommended = metric_map.get("recommended_read_proportion")
    if overall_status == "pass":
        interpretation = "Passed Ribo-seq identity checks."
    elif overall_status == "warn":
        interpretation = "Borderline on one or more Ribo-seq identity checks — review before frame-dependent analysis."
    elif overall_status == "fail":
        interpretation = "Fails one or more Ribo-seq identity checks — frame-dependent analysis not advised without investigation."
    else:
        interpretation = "Gated identity metrics were not available."
    details = []
    if periodicity and periodicity.get("raw") is not None:
        details.append(f"periodicity {periodicity['raw_label']}")
    if cds and cds.get("raw") is not None:
        details.append(f"CDS enrichment E={cds['raw_label']}")
    if recommended and recommended.get("raw") is not None:
        details.append(f"recommended-read proportion {recommended['raw_label']}")
    if details:
        interpretation = f"{interpretation} Top-line: {', '.join(details)}."

    # --- Context strip --------------------------------------------------
    raw_context = summary.get("context", {})
    lib_label = raw_context.get("library_type") or "unknown"
    dom_rl = raw_context.get("dominant_read_length")
    total_reads = raw_context.get("total_reads")
    annotation = raw_context.get("annotation")
    bam = raw_context.get("bam", "")
    context_strip = {
        "library_type": lib_label,
        "dominant_read_length": f"{dom_rl} nt" if dom_rl else "n/a",
        "total_reads": f"{total_reads:,}" if total_reads else "n/a",
        "annotation": annotation.split("/")[-1] if annotation else "none",
        "bam": bam.split("/")[-1] if bam else "n/a",
    }

    # --- Diagnostics list -----------------------------------------------
    diagnostics = summary.get("diagnostics", [])

    return {
        "top_cards": top_cards,
        "tier_sections": tier_sections,
        "grouped_metrics": grouped_metrics,   # kept for backwards-compat
        "overall_status": overall_status,
        "interpretation": interpretation,
        "context_strip": context_strip,
        "diagnostics": diagnostics,
    }


def generate_report(
    plots: List[Dict[str, Any]],
    config: Dict[str, Any],
    export_mode: str = "html",
    name: str = "RiboMetric_report",
    outdir: str = "",
) -> None:
    """
    Generates a report of the RiboMetric results with plots

    Inputs:
        plots: A list containing the plots and metrics for the report
        export_mode: A string defining the mode of export: 'html', 'pdf' or
        'both' (Default: 'html')
        name: A string for the file name (Default: 'RiboMetric_report')
        outdir: A string for the output directory (Default: '')

    Outputs:
        No variables will be output
    """
    project_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env = Environment(
        loader=FileSystemLoader(
            ["templates", f"{project_path}/RiboMetric/templates"]
        ),
        autoescape=False,
    )

    file_names = {"bam": config["argument"]["bam"].split("/")[-1]}
    if config["argument"]["annotation"] is not None:
        file_names["annotation"] = (config["argument"]["annotation"]
                                    .split("/")[-1])
    elif config["argument"]["gff"] is not None:
        file_names["annotation"] = config["argument"]["gff"].split("/")[-1]

    completion_time = datetime.now().strftime("%H:%M:%S %d/%m/%Y")

    binary_logo = open(
        f"{project_path}/RiboMetric/templates/RiboMetric_logo.png",
        "rb"
        ).read()
    base64_logo = base64.b64encode(binary_logo).decode("utf-8")

    binary_icon = open(
        f"{project_path}/RiboMetric/templates/RiboMetric_favicon.png",
        "rb"
        ).read()
    base64_icon = base64.b64encode(binary_icon).decode("utf-8")

    if outdir == "":
        output = name
    else:
        if outdir.endswith("/") and outdir != "":
            outdir = outdir[:-1]
        output = outdir + "/" + name

    if export_mode == "both":
        export_mode_list: List[str] = ["html", "pdf"]
    else:
        export_mode_list = [export_mode]

    template = env.get_template("base.html")
    summary = plots.pop(0)
    context = {
        "summary": summary,
        "report_context": build_report_context(summary),
        "plots": plots,
        "file_names": file_names,
        "completion_time": completion_time,
        "logo": base64_logo,
        "favicon": base64_icon,
    }

    for filetype in export_mode_list:
        if filetype == "html":
            context["filetype"] = filetype
            jinja_render = template.render(context)
            out = output + ".html"
            with open(out, mode="w", encoding="utf-8") as f:
                f.write(jinja_render)
            print(f"Your {filetype} report can be found in {out}")
        else:
            context["filetype"] = filetype
            jinja_render = template.render(context)
            out = output + ".pdf"
            convert_html_to_pdf(jinja_render, out)
            print(f"Your {filetype} report can be found in {out}")


def int_keys_hook(data: Dict[Any, Any]) -> Dict[Any, Any]:
    """
    Custom object_hook for JSON parsing that converts number strings into
    integers
    """
    for key in list(data.keys()):
        if isinstance(key, str) and key.isdigit():
            data[int(key)] = data.pop(key)
    return data


def parse_json_input(json_path: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Parse json input from a previous RiboMetric run for use in generating plots

    Inputs:
        json_path: Path to json file

    Outputs:
        results_dict: Dictionary containing results from a RiboMetric analysis
        json_config: Config from the RiboMetric analysis
    """
    with open(json_path, 'r') as json_file:
        json_dict = json.load(json_file, object_hook=int_keys_hook)
    result_dict = json_dict["results"]
    json_config = json_dict["config"]
    return (result_dict, json_config)
