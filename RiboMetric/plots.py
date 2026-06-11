"""
This script contains the code for generating the plots for
RiboMetric reports
"""

import os
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from .modules import read_frame_cull, read_frame_score_trips_viz, sum_mRNA_distribution
from .scoring import build_scored_metrics
import plotly.io as pio
import base64

import plotly.express as px
import pandas as pd


def _prepare_report_figure(fig: go.Figure) -> go.Figure:
    """Apply report-wide sizing that leaves room for axis labels."""
    fig.update_layout(
        autosize=True,
        height=460,
        margin=dict(l=78, r=36, t=76, b=92),
    )
    fig.update_xaxes(automargin=True, title_standoff=18)
    fig.update_yaxes(automargin=True, title_standoff=18)
    return fig


def plotly_to_html(fig: go.Figure) -> str:
    """Return a report plot fragment without embedding Plotly.js repeatedly."""
    return pio.to_html(
        _prepare_report_figure(fig),
        full_html=False,
        include_plotlyjs=False,
        config={"responsive": True},
    )


def generate_plots(results_dict: dict, config: dict) -> list:
    """
    Wrapper function generating plots based on the results_dict from qc.py

    Input:
        results_dict: Dictionary containing result from modules after running
        through qc.py
        config: Dictionary containing the configuration information

    Output:

    """
    print("Generating plots")
    plots_list = [plot_metrics_summary(results_dict, config)]
    plots_list.append(
            plot_read_frame_distribution(
                results_dict["read_frame_distribution"], config
            )
    )
    if results_dict["mode"] == "annotation":
        plots_list.extend(
            [
                plot_metagene_profile(
                    results_dict["metagene_profile"],
                    config,
                    results_dict.get("metagene_profile_stats"),
                ),
                plot_metagene_heatmap(
                    results_dict["metagene_profile"],
                    config,
                    results_dict.get("metagene_profile_stats"),
                ),
                plot_mRNA_distribution(
                    results_dict["mRNA_distribution"],
                    config,
                ),
                plot_mRNA_read_breakdown(
                    results_dict["mRNA_distribution"],
                    config,
                ),
                plot_read_frame_triangle(
                    results_dict["reading_frame_triangle"],
                    config
                ),
            ]
        )
    if "nucleotide_composition" in results_dict:
        plots_list.extend([
                plot_terminal_nucleotide_bias_distribution(
                    results_dict["terminal_nucleotide_bias_distribution"], config
                ),
                plot_nucleotide_composition(
                    results_dict["nucleotide_composition"], config
                ),
        ])
    plots_list.extend([
            plot_read_length_distribution(
                results_dict["read_length_distribution"], config
            ),
            ])

    if results_dict.get("recommended_read_lengths"):
        plots_list.append(
            plot_recommended_read_lengths(
                results_dict["recommended_read_lengths"], config
            )
        )
    if results_dict.get("gene_body_coverage", {}).get("profile"):
        plots_list.append(
            plot_gene_body_coverage(results_dict["gene_body_coverage"], config)
        )
    if results_dict.get("library_complexity", {}).get("fractions"):
        plots_list.append(
            plot_library_complexity(results_dict["library_complexity"], config)
        )
    if results_dict.get("codon_dwell_times", {}).get("dwell_times"):
        plots_list.append(
            plot_codon_dwell_times(results_dict["codon_dwell_times"], config)
        )
    if results_dict.get("floss", {}).get("floss_scores"):
        plots_list.append(
            plot_floss_heterogeneity(results_dict["floss"], config)
        )

    return plots_list


def plot_codon_dwell_times(dwell: dict, config: dict) -> dict:
    """Bar chart of A-site codon dwell-times, sorted from slowest to fastest."""
    dwell_times = dwell.get("dwell_times", {})
    items = sorted(dwell_times.items(), key=lambda kv: kv[1], reverse=True)
    codons = [c for c, _ in items]
    values = [v for _, v in items]
    pro = {"CCT", "CCC", "CCA", "CCG"}
    colors = [
        "#c93434" if c == "CGA" else "#f0a81f" if c in pro else "#2e85db"
        for c in codons
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=codons, y=values, marker_color=colors, name="",
            hovertemplate="<b>Codon</b>: %{x}"
                          "<br><b>Dwell-time</b>: %{y:.2f}<extra></extra>",
        )
    )
    fig.add_hline(y=1.0, line_dash="dash", line_color=config["plots"]["base_color"])
    fig.update_layout(
        title="A-site Codon Dwell-times",
        xaxis_title="Codon (slowest -> fastest)",
        yaxis_title="Dwell-time (observed / expected)",
        font=dict(
            family=config["plots"]["font_family"],
            size=14,
            color=config["plots"]["base_color"],
        ),
    )
    return {
        "name": "Codon Dwell-times",
        "description": (
            "A-site occupancy per codon relative to its abundance. Values >1 "
            "(above the dashed line) are slow/pausing codons; proline codons are "
            "orange and the inhibitory CGA codon is red."
        ),
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(
            fig, config["plots"]["image_size"][0],
            config["plots"]["image_size"][1]
        ),
    }


def plot_floss_heterogeneity(floss: dict, config: dict) -> dict:
    """Histogram of per-transcript FLOSS read-length heterogeneity scores."""
    scores = floss.get("floss_scores", [])
    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=scores, nbinsx=30, name="",
            hovertemplate="<b>FLOSS</b>: %{x:.2f}"
                          "<br><b>Transcripts</b>: %{y}<extra></extra>",
        )
    )
    fig.update_layout(
        title=(
            "FLOSS Read-length Heterogeneity "
            f"(median {floss.get('floss_median')}, "
            f"{floss.get('floss_aberrant_transcript_fraction', 0):.0%} aberrant)"
        ),
        xaxis_title="Per-transcript FLOSS (departure from library length profile)",
        yaxis_title="Number of transcripts",
        font=dict(
            family=config["plots"]["font_family"],
            size=16,
            color=config["plots"]["base_color"],
        ),
    )
    return {
        "name": "FLOSS Heterogeneity",
        "description": (
            "Distribution of per-transcript FLOSS scores (footprint-length "
            "departure from the library aggregate). A long right tail means many "
            "transcripts have anomalous length profiles — a heterogeneous or "
            "contaminated library."
        ),
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(
            fig, config["plots"]["image_size"][0],
            config["plots"]["image_size"][1]
        ),
    }


def plot_recommended_read_lengths(recommended: dict, config: dict) -> dict:
    """Bar chart of per-read-length periodicity, highlighting recommended lengths."""
    by_rl = recommended.get("by_read_length", {})
    read_lengths = sorted(by_rl.keys())
    periodicity = [by_rl[rl]["periodicity"] for rl in read_lengths]
    colors = [
        "#2e85db" if by_rl[rl]["recommended"] else "#c9c9c9"
        for rl in read_lengths
    ]
    customdata = [
        [by_rl[rl].get("read_proportion", 0), by_rl[rl].get("offset", "n/a")]
        for rl in read_lengths
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=read_lengths,
            y=periodicity,
            marker_color=colors,
            customdata=customdata,
            hovertemplate=(
                "<b>Read length</b>: %{x}"
                "<br><b>Periodicity</b>: %{y:.3f}"
                "<br><b>Read proportion</b>: %{customdata[0]:.3f}"
                "<br><b>P-site offset</b>: %{customdata[1]}<extra></extra>"
            ),
        )
    )
    # 1/3 random-frame baseline: present the floor so the user reads the raw
    # periodicity value against it rather than having it hidden in scaling.
    fig.add_hline(
        y=1 / 3,
        line_dash="dash",
        line_color="#aaaaaa",
        annotation_text="Random-frame baseline (1/3 ≈ 0.33)",
        annotation_position="right",
        annotation_font_size=12,
    )
    fig.update_layout(
        title=(
            "Recommended read lengths "
            f"({recommended.get('n_recommended', 0)} selected, "
            f"{recommended.get('recommended_read_proportion', 0):.0%} of reads)"
        ),
        xaxis_title="Read Length",
        yaxis_title="Dominant-frame fraction (periodicity)",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    return {
        "name": "Recommended Read Lengths",
        "description": (
            "Per-read-length 3-nt periodicity. Blue bars are recommended for "
            "downstream P-site assignment / ORF calling (periodic enough and "
            "carrying enough reads); grey bars are not."
        ),
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(
            fig, config["plots"]["image_size"][0],
            config["plots"]["image_size"][1]
        ),
    }


def plot_gene_body_coverage(gene_body: dict, config: dict) -> dict:
    """Line plot of A-site density across relative CDS position (0-100%)."""
    profile = gene_body.get("profile", [])
    x = [round(100 * i / max(1, len(profile) - 1), 2) for i in range(len(profile))]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x, y=profile, mode="lines", name="",
            hovertemplate="<b>CDS position</b>: %{x:.0f}%"
                          "<br><b>Relative density</b>: %{y:.2f}<extra></extra>",
        )
    )
    fig.add_hline(y=1.0, line_dash="dash", line_color=config["plots"]["base_color"])
    fig.update_layout(
        title="Gene-body Coverage (5'-&gt;3')",
        xaxis_title="Relative position along CDS (%)",
        yaxis_title="Mean-normalised A-site density",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    return {
        "name": "Gene-body Coverage",
        "description": (
            "A-site density across relative CDS position. A 5' bump is the "
            "translation ramp; a 3' decline indicates drop-off. A flat line at "
            "1.0 is perfectly uniform elongation."
        ),
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(
            fig, config["plots"]["image_size"][0],
            config["plots"]["image_size"][1]
        ),
    }


def plot_library_complexity(complexity: dict, config: dict) -> dict:
    """Rarefaction (saturation) curve of distinct A-site positions vs depth."""
    fractions = complexity.get("fractions", [])
    distinct = complexity.get("distinct_positions", [])
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=fractions, y=distinct, mode="lines+markers", name="",
            hovertemplate="<b>Fraction of reads</b>: %{x:.0%}"
                          "<br><b>Distinct positions</b>: %{y}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Library Complexity (rarefaction)",
        xaxis_title="Fraction of reads sampled",
        yaxis_title="Expected distinct A-site positions",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    return {
        "name": "Library Complexity",
        "description": (
            "Expected distinct A-site positions recovered vs sequencing depth. "
            "A curve that has flattened by full depth is saturated (extra "
            "sequencing adds little); a still-rising curve is under-sequenced."
        ),
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(
            fig, config["plots"]["image_size"][0],
            config["plots"]["image_size"][1]
        ),
    }


_IMG_WARNED = False
_TRANSPARENT_PNG_1x1 = (
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMB/axPpVwAAAAASUVORK5CYII="
)


def plotly_to_image(fig: go.Figure, width: int, height: int) -> str:
    """Return a base64 PNG for PDF export, or a tiny placeholder if imaging is unavailable.

    - Respects RIBOMETRIC_SKIP_IMAGES=1 to force placeholder (CI/headless).
    - Catches Kaleido/Chrome RuntimeError and falls back silently to a 1x1 PNG.
    """
    global _IMG_WARNED
    if os.environ.get("RIBOMETRIC_SKIP_IMAGES") == "1":
        return _TRANSPARENT_PNG_1x1
    try:
        _prepare_report_figure(fig)
        img = pio.to_image(fig, format="png", width=width, height=height)
        return base64.b64encode(img).decode("ascii")
    except Exception as e:
        # Avoid noisy tracebacks when Chrome/Kaleido is missing; warn once.
        if not _IMG_WARNED:
            print("Note: static image export unavailable (", str(e).split("\n")[0], ") — continuing without images.")
            _IMG_WARNED = True
        return _TRANSPARENT_PNG_1x1


def plot_read_length_distribution(
    read_length_dict: dict, config: dict
) -> dict:
    """
    Generate a plot of the read length distribution for the full dataset

    Inputs:
        read_length_df: Dataframe containing the read length distribution
        config: Dictionary containing the configuration information

    Outputs:
        plot_read_length_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    hovertemplate = "<b>Read length</b>: %{x}" + "<br><b>Count</b>: %{y}"
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=list(read_length_dict.keys()),
            y=list(read_length_dict.values()),
            name="",
            hovertemplate=hovertemplate,
        )
    )
    fig.update_layout(
        title="Read Length Distribution",
        xaxis_title="Read Length",
        yaxis_title="Read Count",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    plot_read_length_dict = {
        "name": "Read Length Distribution",
        "description": "Distribution of read lengths for the full dataset",
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_read_length_dict


def plot_terminal_nucleotide_bias_distribution(
    terminal_nucleotide_bias_dict: dict, config: dict
) -> dict:
    """
    Generate a plot of ligation bias distribution for the full dataset

    Inputs:
        read_length_df: Dataframe containing the read length distribution
        config: Dictionary containing the configuration information

    Outputs:
        plot_terminal_nucleotide_bias_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    columns = 1
    if terminal_nucleotide_bias_dict["five_prime"] != {}:
        target_loop = ["five_prime"]
        if terminal_nucleotide_bias_dict["three_prime"] != {}:
            target_loop.append("three_prime")
            columns = 2
    else:
        target_loop = ["three_prime"]

    fig = make_subplots(
        rows=1,
        cols=int(columns),
        shared_yaxes=True,
        subplot_titles=["5' Ligation Bias", "3' Ligation Bias"]
        if len(target_loop) > 1
        else ["5' Ligation Bias"]
        if target_loop == ["five_prime"]
        else ["3' Ligation Bias"],
    )
    count = 0
    for current_target in target_loop:
        count += 1
        # Set colors according to the value
        color = ["#636efa" if value > 0 else "#da5325"
                 for value in terminal_nucleotide_bias_dict[current_target].values()]
        fig.add_trace(
            go.Bar(
                x=list(terminal_nucleotide_bias_dict[current_target].keys()),
                y=list(terminal_nucleotide_bias_dict[current_target].values()),
                name="",
                hovertemplate="<b>Nucleotides</b>:%{x}<br><b>Bias</b>:%{y}",
                marker=dict(color=color),
            ),
            row=1,
            col=count,
        )

        fig.add_hline(y=0)

        fig.update_xaxes(
            title_text="Read Start",
            row=1,
            col=count,
            title_font=dict(size=18)
        )
    fig.update_layout(
        title="Ligation Bias Distribution",
        yaxis_title="Proportion",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
        showlegend=False,
    )
    distance = 0.2
    for prime in terminal_nucleotide_bias_dict:
        for bias in terminal_nucleotide_bias_dict[prime].values():
            if abs(bias) > distance:
                distance = abs(bias)
    max_range = [-distance, distance]
    fig.update_yaxes(range=max_range)
    if columns > 1:
        fig.update_layout(
            xaxis=dict(
                domain=[0, 0.48], zeroline=False
            ),  # Adjust domain and remove x-axis zeroline for subplot 1
            xaxis2=dict(
                domain=[0.52, 1], zeroline=False
            ),  # Adjust domain and remove x-axis zeroline for subplot 2
        )
    else:
        fig.update_layout(
            xaxis=dict(domain=[0, 1], zeroline=False),
        )
    plot_terminal_nucleotide_bias_dict = {
        "name": "Ligation Bias Distribution",
        "description": "Distribution of end bases for the full dataset",
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_terminal_nucleotide_bias_dict


def plot_nucleotide_composition(
    nucleotide_composition_dict: dict, config: dict
) -> dict:
    """
    Generate a plot of the nucleotide composition for the full dataset

    Inputs:
        nucleotide_composition_dict: Dictionary containing the distribution of
        nucleotides per read position
        config: Dictionary containing the configuration information

    Outputs:
        plot_nucleotide_composition_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    colors = config["plots"]["nucleotide_colors"]
    fig = go.Figure()
    for nucleotide, distribution in nucleotide_composition_dict.items():
        fig.add_trace(
            go.Scatter(
                y=distribution, name=nucleotide, line_color=colors[nucleotide]
            )
        )
    fig.update_layout(
        title="Nucleotide Composition",
        xaxis_title="Position (nucleotides)",
        yaxis_title="Proportion",
        yaxis_range=[0, 1],
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    plot_nucleotide_composition_dict = {
        "name": "Nucleotide Composition",
        "description": "Nucleotide composition of the reads",
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_nucleotide_composition_dict


def plot_nucleotide_distribution(
    nucleotide_composition_dict: dict, config: dict
) -> dict:
    plot_data = []
    nt_start, nt_count = (
        config["plots"]["nucleotide_proportion"]["nucleotide_start"],
        config["plots"]["nucleotide_proportion"]["nucleotide_count"],
    )
    for nt in reversed(nucleotide_composition_dict):
        plot_data.append(
            go.Bar(
                name=nt,
                x=[*range(nt_start + 1, nt_start + nt_count + 1)],
                y=nucleotide_composition_dict[nt][
                    nt_start: nt_start + nt_count
                ],
                marker=dict(color=config["plots"]["nucleotide_colors"][nt]),
                # Set the text in the hovertemplate to proportion or count
                # depending on config
                hovertemplate="Proportion: %{y:.2%}"
                if not config["plots"]["mRNA_distribution"]["absolute_counts"]
                else "Count: %{x}",
            )
        )
    fig = go.Figure(plot_data)
    fig.update_layout(
        barmode="stack",
        title="Nucleotide Proportion",
        xaxis_title="",
        yaxis_title="Proportion",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
        legend={"traceorder": "reversed"},
    )
    plot_nucleotide_distribution_dict = {
        "name": "Nucleotide Distribution",
        "description": "Nucleotide distribution across specified reads \
(default: first 15 read)",
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_nucleotide_distribution_dict


def plot_read_frame_distribution(read_frame_dict: dict, config: dict) -> dict:
    """
    Generate a plot of the read frame distribution

    Inputs:
        read_frame_dict: Dataframe containing the read frame distribution
        config: Dictionary containing the configuration information

    Outputs:
        plot_read_frame_dict: Dictionary containing the plot name, description
        and plotly figure for html and pdf export
    """
    culled_read_frame_dict = read_frame_cull(read_frame_dict, config)
    # Calculates the read frame scores if 'show_scores' option
    # in config is not 'none'
    scored_read_frame_dict = (
        read_frame_score_trips_viz(culled_read_frame_dict)
        if config["plots"]["read_frame_distribution"]["show_scores"] != "none"
        else None
    )

    # Set minimum and maximum font sizes
    min_font_size, max_font_size = 5, 30

    # Calculate font size based on number of data points
    num_data_points = len(culled_read_frame_dict)
    font_size = max_font_size - (max_font_size - min_font_size) * (
        num_data_points / 50
    )
    # Generate plot
    plot_data = []
    for i in range(0, 3):
        plot_data.append(
            go.Bar(
                name="Frame " + str(i + 1),
                x=list(culled_read_frame_dict.keys()),
                #
                y=[
                    culled_read_frame_dict[x][y]
                    for x in culled_read_frame_dict
                    for y in culled_read_frame_dict[x]
                    if y == i
                ],
            )
        )
    fig = go.Figure(data=plot_data)
    fig.update_layout(barmode="group")
    fig.update_layout(
        title="Read Frame Distribution",
        xaxis_title="Read Length",
        yaxis_title="Read Count",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    # Place scores dynamically over bars
    if scored_read_frame_dict is not None:
        if config["plots"]["read_frame_distribution"]["show_scores"] == "all":
            for idx in enumerate(culled_read_frame_dict):
                if idx[1] != "global":
                    y_buffer = (
                        max(fig.data[0].y + fig.data[1].y + fig.data[2].y)
                        * 0.05
                    )

                    ymax = max(
                        fig.data[0].y[idx[0]],
                        fig.data[1].y[idx[0]],
                        fig.data[2].y[idx[0]],
                    )

                    if (
                        fig.data[0].y[idx[0]]
                        + fig.data[1].y[idx[0]]
                        + fig.data[2].y[idx[0]]
                        > y_buffer
                    ):
                        fig.add_annotation(
                            x=idx[1],
                            y=ymax + y_buffer,
                            text=round(scored_read_frame_dict[idx[1]], 2),
                            showarrow=False,
                            xanchor="center",
                            font={"size": font_size},
                        )
        fig.add_annotation(
            text=f'Score: {round(scored_read_frame_dict["global"], 2)}',
            showarrow=False,
            xref="paper",
            yref="paper",
            y=0.64,
            x=1.03,
            xanchor="left",
        )

    # Trim the x-axis to the populated read-length range. Guard against an
    # empty frame dict (e.g. no reads survived annotation/CDS filtering) and the
    # case where every bar is below the buffer, both of which would otherwise
    # leave lower_limit/upper_limit unbound.
    all_y = fig.data[0].y + fig.data[1].y + fig.data[2].y
    y_buffer = (max(all_y) * 0.05) if all_y else 0
    lower_limit = upper_limit = None
    for idx in enumerate(culled_read_frame_dict):
        count_sum = (
            fig.data[0].y[idx[0]]
            + fig.data[1].y[idx[0]]
            + fig.data[2].y[idx[0]])
        if count_sum > y_buffer:
            lower_limit = (idx[1])
            break
    for idx in list(enumerate((culled_read_frame_dict)))[::-1]:
        count_sum = (
            fig.data[0].y[idx[0]]
            + fig.data[1].y[idx[0]]
            + fig.data[2].y[idx[0]])
        if count_sum > y_buffer:
            upper_limit = (idx[1])
            break
    if lower_limit is not None and upper_limit is not None:
        fig.update_xaxes(range=[lower_limit-0.5, upper_limit+0.5])
    plot_read_frame_dict = {
        "name": "Read Frame Distribution",
        "description": "Frame distribution per read length",
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_read_frame_dict


def plot_mRNA_distribution(mRNA_distribution_dict: dict, config: dict) -> dict:
    """
    Generate a bar plot of the mRNA distribution

    Inputs:
        mRNA_distribution_dict: Dictionary containing the mRNA distribution
        over the read lengths
        config: Dictionary containing the configuration information

    Outputs:
        plot_mRNA_distribution_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    sum_mRNA_dict = sum_mRNA_distribution(mRNA_distribution_dict, config)
    plotting_order = config["plots"]["mRNA_read_breakdown"]["plotting_order"]
    plot_data = []

    for category in plotting_order:
        if category in sum_mRNA_dict:
            value = sum_mRNA_dict[category]
            plot_data.append(
                go.Bar(
                    name=category,
                    x=[value],
                    y=[""],
                    width=[0.3],
                    hovertemplate="Proportion: %{x:.2%}"
                    if not config[
                        "plots"][
                        "mRNA_distribution"][
                        "absolute_counts"]
                    else "Count: %{x}",
                    orientation="h",
                )
            )

    fig = go.Figure(plot_data)
    fig.update_layout(
        barmode="stack",
        title="mRNA Reads Breakdown",
        xaxis_title="Proportion"
        if not config["plots"]["mRNA_read_breakdown"]["absolute_counts"]
        else "Counts",
        yaxis_title="",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
        legend={"traceorder": "normal"},
    )
    fig.update_xaxes(range=[0, 1])
    plot_mRNA_distribution_dict = {
        "name": "mRNA Reads Breakdown",
        "description": "Shows the proportion of the different transcript \
regions represented in the reads",
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_mRNA_distribution_dict


from typing import Dict, List


def plot_mRNA_read_breakdown(
    mRNA_distribution_dict: Dict[int, Dict[str, int]], config: dict
) -> dict:
    """
    Generate a line plot of the mRNA distribution over the read lengths

    Inputs:
        mRNA_distribution_dict: Dictionary containing the mRNA distribution
        over the read lengths
        config: Dictionary containing the configuration information

    Outputs:
        plot_mRNA_distribution_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    plot_data: Dict[str, List[float]] = {}
    for read_length_dict in mRNA_distribution_dict.values():
        for category, count in read_length_dict.items():
            plot_data.setdefault(category, []).append(float(count))
    if not config["plots"]["mRNA_read_breakdown"]["absolute_counts"]:
        sum_data = {k: sum(v) for k, v in plot_data.items()}
        plot_data = {
            k: [x / sum(sum_data.values()) for x in v]
            for k, v in plot_data.items()
        }

    fig = go.Figure()
    for k, v in plot_data.items():
        fig.add_trace(
            go.Scatter(
                name=k,
                x=list(mRNA_distribution_dict.keys()),
                y=v,
                hovertemplate="Proportion: %{y:.2%}"
                if not config["plots"]["mRNA_read_breakdown"][
                    "absolute_counts"
                ]
                else "Count: %{x}",
            )
        )

    fig.update_layout(
        title="Nucleotide Distribution",
        xaxis_title="Position (nucleotides)",
        yaxis_title="Proportion"
        if not config["plots"]["mRNA_read_breakdown"]["absolute_counts"]
        else "Counts",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    fig.update_layout(
        title="Nucleotide Distribution",
        xaxis_title="Read length",
        yaxis_title="Proportion"
        if not config["plots"]["mRNA_read_breakdown"]["absolute_counts"]
        else "Counts",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
        legend={"traceorder": "normal"},
    )
    plot_mRNA_read_breakdown_dict = {
        "name": "mRNA Reads Breakdown over Read Length",
        "description": "Shows the proportion of the different transcript \
regions represented in the reads over the different read lengths.",
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_mRNA_read_breakdown_dict


def _metagene_stats_label(stats: dict | None) -> str:
    if not stats:
        return ""
    if stats.get("multimap_filter") == "unique_only":
        population = "MAPQ=255 unique reads"
    else:
        population = "all primary alignments"
    start_n = stats.get("start_window_reads", 0)
    stop_n = stats.get("stop_window_reads", 0)
    profile_n = stats.get("profile_reads", 0)
    input_n = stats.get("input_reads", 0)
    return (
        f"{population}; profile reads n={profile_n:,}/{input_n:,}; "
        f"window reads: start n={start_n:,}, stop n={stop_n:,}"
    )


def _metagene_description(base: str, stats: dict | None) -> str:
    stats_label = _metagene_stats_label(stats)
    if not stats_label:
        return base
    return f"{base} {stats_label}."


def plot_metagene_profile(
    metagene_profile_dict: dict,
    config: dict,
    stats: dict | None = None,
) -> dict:
    """
    Generate a plot of the distribution of reads depending on their distance
    to a target (default: start codon)

    Inputs:
        metagene_dict: Dictionary containing the counts as values and distance
        from target as keys
        config: Dictionary containing the configuration information

    Outputs:
        plot_metagene_profile_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    count = 0
    columns = 1
    frame_colors = {0: "#636efa", 1: "#ef553b", 2: "#00cc96"}
    if metagene_profile_dict["start"] != {}:
        target_loop = ["start"]
        if metagene_profile_dict["stop"] != {}:
            target_loop.append("stop")
            columns = 2
    else:
        target_loop = ["stop"]

    fig = make_subplots(
        rows=1,
        cols=columns,
        shared_yaxes=config["plots"]["metagene_profile"]["shared_yaxis"],
        subplot_titles=["Distance from 5'", "Distance from 3'"]
        if len(target_loop) > 1
        else ["Distance from 5'"]
        if target_loop == ["start"]
        else ["Distance from 3'"],
    )
    for current_target in target_loop:
        count += 1
        metagene_dict: Dict[int, int] = {}
        for inner_dict in metagene_profile_dict[current_target].values():
            for inner_key, inner_value in inner_dict.items():
                position = int(inner_key)
                count_value = inner_value if inner_value is not None else 0
                if (
                    position in metagene_dict
                    and metagene_dict[position] is not None
                ):
                    metagene_dict[position] += count_value
                else:
                    metagene_dict[position] = count_value
        sorted_positions = sorted(metagene_dict)
        color_frames: List[int] = [position % 3 for position in sorted_positions]
        color: List[str] = [frame_colors[i] for i in color_frames]

        fig.add_trace(
            go.Bar(
                x=sorted_positions,
                y=[metagene_dict[position] for position in sorted_positions],
                name="Distance from 5'"
                if current_target == "start"
                else "Distance from 3'",
                marker=dict(color=color),
            ),
            row=1,
            col=count,
        )
        fig.update_xaxes(
            title_text="Relative position (nt)",
            row=1,
            col=count,
            title_font=dict(size=18),
            range=config["plots"]["metagene_profile"]["distance_range"]
        )

    fig.update_layout(
        title="Metagene Profile",
        yaxis_title="Read Count",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
        bargap=0,
        showlegend=False,
    )
    if columns > 1:
        fig.update_layout(
            xaxis=dict(
                domain=[0, 0.48], zeroline=False
            ),  # Adjust domain and remove x-axis zeroline for subplot 1
            xaxis2=dict(
                domain=[0.52, 1], zeroline=False
            ),  # Adjust domain and remove x-axis zeroline for subplot 2
        )
    else:
        fig.update_layout(
            xaxis=dict(domain=[0, 1], zeroline=False),
        )

    fig.update_xaxes(
        range=config["plots"]["metagene_profile"]["distance_range"]
    )
    plot_metagene_profile_dict = {
        "name": "Metagene Profile",
        "description": _metagene_description(
            "Metagene profile showing read counts by distance from the selected "
            "target. The plotted read population follows the multimap filter.",
            stats,
        ),
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_metagene_profile_dict


def plot_metagene_heatmap(
    metagene_profile_dict: dict,
    config: dict,
    stats: dict | None = None,
) -> dict:
    """
    Generate a heatmap of the reads depending on their distance
    to a target, read length and count

    Inputs:
        metagene_heatmap_dict: Dictionary containing the counts as values
            and distance from target as keys
        config: Dictionary containing the configuration information

    Outputs:
        plot_metagene_heatmap: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    count = 0
    columns = 1
    if metagene_profile_dict["start"] != {}:
        target_loop = ["start"]
        if metagene_profile_dict["stop"] != {}:
            target_loop.append("stop")
            columns = 2
    else:
        target_loop = ["stop"]

    fig = make_subplots(
        rows=1,
        cols=columns,
        shared_yaxes=True,
        subplot_titles=["Distance from 5'", "Distance from 3'"]
        if len(target_loop) > 1
        else ["Distance from 5'"]
        if target_loop == ["start"]
        else ["Distance from 3'"],
    )

    for current_target in target_loop:
        count += 1
        x_data: List[int] = []
        y_data: List[int] = []
        z_data: List[int] = []
        for read_length, position_counts in metagene_profile_dict[
                                                current_target
                                                ].items():
            for position, counts in position_counts.items():
                x_data.append(int(position))
                y_data.append(int(read_length))
                z_data.append(int(counts))
        if x_data:
            ordered = sorted(zip(x_data, y_data, z_data), key=lambda row: (row[1], row[0]))
            x_data = [row[0] for row in ordered]
            y_data = [row[1] for row in ordered]
            z_data = [row[2] for row in ordered]

        if config["plots"]["metagene_profile"]["max_colorscale"] is None:
            z_max = max(max(z_data, default=0), 1)
            # Keep integer typing by scaling and rounding to nearest int for display
            z_data = [int(round((z / z_max) * z_max)) for z in z_data]

        fig.add_trace(
            go.Heatmap(
                x=x_data,
                y=y_data,
                z=z_data,
                colorscale=config["plots"]["metagene_profile"]["colorscale"],
                zmin=0,
                zmax=config["plots"]["metagene_profile"]["max_colorscale"],
            ),
            row=1,
            col=count,
        )

        fig.update_xaxes(
            title_text="Relative position (nt)",
            row=1,
            col=count,
            title_font=dict(size=18),
            range=[-50, 50],
            constrain="domain",
            zeroline=False,
            # range=config["plots"]["metagene_profile"]["distance_range"]
        )

    fig.update_layout(
        title="Metagene Heatmap",
        yaxis_title="Read length",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
        legend={"traceorder": "normal"},
        showlegend=False,
    )

    plot_metagene_heatmap = {
        "name": "Metagene Heatmap",
        "description": _metagene_description(
            "Metagene heatmap showing A-site distance from the selected target "
            "by read length. The plotted read population follows the multimap "
            "filter.",
            stats,
        ),
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_metagene_heatmap


def plot_metrics_summary(results_dict: dict, config: dict) -> dict:
    '''
    Generate the metrics summary payload used by the HTML/PDF report.

    Accepts either the full results_dict (preferred) or a plain metrics dict
    for backwards compatibility.

    Returns a dict with:
        plot    — bar-chart of anchored 0-1 scores
        metrics — scored records (key, score, raw, status, gate, tier, decision)
        context — context-strip data for the report header
        diagnostics — non-scored raw values shown as plain labels
    '''
    # Accept either full results_dict or bare metrics dict (legacy callers).
    if "metrics" in results_dict:
        metrics_dict = results_dict["metrics"]
    else:
        metrics_dict = results_dict
        results_dict = {}

    # Anchored scores come from the single scoring resolver (scoring.py).
    scored = build_scored_metrics({"metrics": metrics_dict}, config)

    # Drop any metrics explicitly excluded from the summary plot.
    excluded = set(config.get("plots", {}).get("exclude_metrics", []))
    scored = [m for m in scored if m["key"] not in excluded]

    # Only metrics that produced a score are plottable on the 0-1 bar chart.
    plottable = [m for m in scored if m["score"] is not None]
    df = pd.DataFrame(
        [{"Metric": m["key"], "Score": round(m["score"], 3)} for m in plottable]
    )

    width = 750
    height = 320
    fig = px.bar(df, x='Metric', y='Score', title="Summary Scores",
                 labels={'Score': 'Score'},
                 hover_data={'Metric': True, 'Score': ':.3f'},
                 height=400)
    fig.update_yaxes(range=[0, 1])
    fig.update_layout(showlegend=False)
    fig_html = plotly_to_html(fig)
    fig_image = plotly_to_image(fig, width, height)

    # --- Context strip ---------------------------------------------------
    rl_dist = results_dict.get("read_length_distribution", {})
    total_reads = int(sum(rl_dist.values())) if rl_dist else None
    dominant_rl = (
        max(rl_dist, key=rl_dist.__getitem__) if rl_dist else None
    )
    lib_type_data = results_dict.get("library_type")
    library_type_label = lib_type_data["label"] if isinstance(lib_type_data, dict) else None
    context = {
        "library_type": library_type_label,
        "dominant_read_length": dominant_rl,
        "total_reads": total_reads,
        "annotation": (
            config.get("argument", {}).get("annotation")
            or config.get("argument", {}).get("gff")
        ),
        "bam": config.get("argument", {}).get("bam"),
    }

    # --- Diagnostics (non-scored values shown as plain labels) -----------
    # These are metrics whose "good direction" depends on experiment type or
    # that don't belong on a universal pass/fail scale.
    _DIAG_KEYS = [
        "disome_proportion",
        "read_length_distribution_IQR_metric",
        "read_length_distribution_coefficient_of_variation_metric",
        "read_length_distribution_maxprop_metric",
        "read_length_distribution_bimodality_metric",
        "start_codon_enrichment_ratio",
        "stop_codon_readthrough_ratio",
        "five_prime_ramp_ratio",
        "three_prime_drop_ratio",
    ]
    _DIAG_LABELS = {
        "disome_proportion": "Di-some proportion",
        "read_length_distribution_IQR_metric": "Read-length IQR score",
        "read_length_distribution_coefficient_of_variation_metric": "Read-length CV",
        "read_length_distribution_maxprop_metric": "Read-length max-proportion",
        "read_length_distribution_bimodality_metric": "Read-length bimodality",
        "start_codon_enrichment_ratio": "Start-codon enrichment ratio",
        "stop_codon_readthrough_ratio": "Stop-codon read-through ratio",
        "five_prime_ramp_ratio": "5′ ramp ratio",
        "three_prime_drop_ratio": "3′ drop ratio",
    }

    def _disome_caption(raw, lib_type):
        if lib_type == "disome":
            return "Expected/desirable in disome-focused experiments."
        if lib_type in ("elongation", "initiation", "low_quality"):
            return "High values may indicate disome contamination or incomplete size selection."
        return "Interpretation depends on library type (monosome: high = potential contamination; disome: expected)."

    diagnostics = []
    for key in _DIAG_KEYS:
        raw = metrics_dict.get(key)
        if raw is None:
            continue
        if isinstance(raw, dict):
            raw = raw.get("global", raw.get("ratio"))
        if raw is None:
            continue
        entry = {
            "key": key,
            "label": _DIAG_LABELS.get(key, key.replace("_", " ").capitalize()),
            "raw": raw,
            "raw_label": f"{raw:.3f}" if isinstance(raw, float) else str(raw),
            "caption": None,
        }
        if key == "disome_proportion":
            entry["caption"] = _disome_caption(raw, library_type_label)
        diagnostics.append(entry)

    return {
        "plot": {
            "name": "Summary of Metrics",
            "description": "Bar chart showing anchored 0-1 scores (higher is better).",
            "fig_html": fig_html,
            "fig_image": fig_image,
        },
        "metrics": [
            {
                "key": m["key"],
                "name": m["key"].replace("_", " ").capitalize(),
                "score": round(m["score"], 3) if m["score"] is not None else None,
                "raw": m["raw"],
                "status": m["status"],
                "gate": m["gate"],
                "tier": m["tier"],
                "decision": m["decision"],
            }
            for m in scored
        ],
        "context": context,
        "diagnostics": diagnostics,
    }


def plot_read_frame_triangle(
        read_frame_triangle_dict: dict, config: dict) -> dict:
    '''
    Generate the triangle plot of the read frame distribution

    Inputs:
        read_frame_triangle_dict: Dictionary containing the read frame
        distribution
        config: Dictionary containing the configuration information

    '''
    data = [[k, v[0], v[1], v[2]] for k, v in read_frame_triangle_dict.items()]
    df = pd.DataFrame(data, columns=['Transcript', 'F1', 'F2', 'F3'])
    fig = px.scatter_ternary(
        df,
        a="F1",
        b="F2",
        c="F3",
        hover_name="Transcript"
        )

    fig.update_layout(
        title="Read Frame Triangle",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    plot_read_frame_triangle_dict = {
        "name": "Read Frame Triangle",
        "description": "Triangle plot showing the distribution of read frames \
for the full dataset",
        "fig_html": plotly_to_html(fig),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }

    return plot_read_frame_triangle_dict
