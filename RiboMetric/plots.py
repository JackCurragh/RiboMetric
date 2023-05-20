"""
This script contains the code for generating the plots for
RiboMetric reports
"""

from plotly import graph_objects as go
from plotly.subplots import make_subplots
from .modules import read_frame_cull, read_frame_score, sum_mRNA_distribution
import tempfile  # logoplot
import subprocess  # logoplot
import plotly.io as pio
import base64


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
    plots_list = [plot_metrics_summary(dict(), config)]
    plots_list.extend(
        [
            plot_read_length_distribution(
                results_dict["read_length_distribution"], config
            ),
            plot_ligation_bias_distribution(
                results_dict["ligation_bias_distribution"], config
            ),
            plot_nucleotide_composition(
                results_dict["nucleotide_composition"], config
            ),
            plot_read_frame_distribution(
                results_dict["read_frame_distribution"], config
            ),
        ]
    )
    if config["plots"]["logoplot"]["enable"]:
        plots_list.append(plot_logoplot(results_dict["sequence_slice"], config))

    if results_dict["mode"] == "annotation":
        print("Finished basic plots, generating annotation plots")
        plots_list.extend(
            [
                plot_mRNA_distribution(
                    results_dict["mRNA_distribution"],
                    config,
                ),
                plot_mRNA_read_breakdown(
                    results_dict["mRNA_distribution"],
                    config,
                ),
                plot_metagene_profile(
                    results_dict["metagene_profile"],
                    config,
                ),
                plot_metagene_heatmap(
                    results_dict["metagene_profile"], config
                ),
            ]
        )
    return plots_list


def plotly_to_image(fig: go.Figure, width: int, height: int) -> str:
    base_64_plot = base64.b64encode(
        pio.to_image(
            fig,
            format="jpg",
            width=width,
            height=height,
        )
    ).decode("ascii")
    return base_64_plot


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
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_read_length_dict


def plot_ligation_bias_distribution(
    ligation_bias_dict: dict, config: dict
) -> dict:
    """
    Generate a plot of ligation bias distribution for the full dataset

    Inputs:
        read_length_df: Dataframe containing the read length distribution
        config: Dictionary containing the configuration information

    Outputs:
        plot_ligation_bias_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    # Remove nucleotide sequences with 'N' if the option in config is False
    if config["plots"]["ligation_bias_distribution"]["include_N"] is False:
        ligation_bias_dict = {
            k: v for k, v in ligation_bias_dict.items() if "N" not in k
        }

    # Set colors according to the value
    color = ["#636efa" if value > 0 else "#da5325" for value in ligation_bias_dict.values()]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=list(ligation_bias_dict.keys()),
            y=list(ligation_bias_dict.values()),
            name="",
            hovertemplate="<b>Nucleotides</b>:%{x}<br><b>Proportion</b>:%{y}",
            marker=dict(color=color),
        )
    )
    fig.add_hline(y=0)
    fig.update_layout(
        title="Ligation Bias Distribution",
        xaxis_title="Read Start",
        yaxis_title="Proportion",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
    )
    plot_ligation_bias_dict = {
        "name": "Ligation Bias Distribution",
        "description": "Distribution of end bases for the full dataset",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_ligation_bias_dict


def plot_nucleotide_composition(
    nucleotide_composition_dict: dict, config: dict
) -> dict:
    """
    Generate a plot of the nucleotide composition for the full dataset

    Inputs:
        read_length_df: Dataframe containing the read length distribution
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
        "fig_html": pio.to_html(fig, full_html=False),
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
        "fig_html": pio.to_html(fig, full_html=False),
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
        read_frame_score(culled_read_frame_dict)
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
                        + fig.data[0].y[idx[0]]
                        + fig.data[0].y[idx[0]]
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
        fig.update_layout()
        fig.add_annotation(
            text=f'Score: {round(scored_read_frame_dict["global"], 2)}',
            showarrow=False,
            xref="paper",
            yref="paper",
            y=0.64,
            x=1.03,
            xanchor="left",
        )

    plot_read_frame_dict = {
        "name": "Read Frame Distribution",
        "description": "Frame distribution per read length",
        "fig_html": pio.to_html(fig, full_html=False),
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
    plot_data = []
    for k, v in sum_mRNA_dict.items():
        plot_data.append(
            go.Bar(
                name=k.replace("_", " ").title(),
                x=[v],
                y=[""],
                width=[0.3],
                hovertemplate="Proportion: %{x:.2%}"
                if not config["plots"]["mRNA_distribution"]["absolute_counts"]
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
    plot_mRNA_distribution_dict = {
        "name": "mRNA Reads Breakdown",
        "description": "Shows the proportion of the different transcript \
regions represented in the reads",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_mRNA_distribution_dict


def plot_mRNA_read_breakdown(
    mRNA_distribution_dict: dict, config: dict
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
    plot_data = {}
    for read_length in mRNA_distribution_dict.values():
        for category, count in read_length.items():
            if category not in plot_data:
                plot_data[category] = []
            plot_data[category].append(count)
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
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_mRNA_read_breakdown_dict


def plot_metagene_profile(metagene_profile_dict: dict, config: dict) -> dict:
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
    colnums = 1
    frame_colors = {0: "#636efa", 1: "#ef553b", 2: "#00cc96"}
    if metagene_profile_dict["start"] != {}:
        target_loop = ["start"]
        if metagene_profile_dict["stop"] != {}:
            target_loop.append("stop")
            colnums = 2
    else:
        target_loop = ["stop"]

    fig = make_subplots(
        rows=1,
        cols=colnums,
        shared_yaxes=config["plots"]["metagene_profile"]["shared_yaxis"],
        subplot_titles=["Distance from 5'", "Distance from 3'"]
        if len(target_loop) > 1
        else ["Distance from 5'"]
        if target_loop == ["start"]
        else ["Distance from 3'"],
    )
    for current_target in target_loop:
        count += 1
        metagene_dict = {}
        for inner_dict in metagene_profile_dict[current_target].values():
            for inner_key, inner_value in inner_dict.items():
                if (
                    inner_key in metagene_dict
                    and metagene_dict[inner_key] is not None
                ):
                    metagene_dict[inner_key] += (
                        inner_value if inner_value is not None else 0
                    )
                else:
                    metagene_dict[inner_key] = (
                        inner_value if inner_value is not None else 0
                    )
        n = 0
        color = [(x % 3) for x in metagene_dict.keys()]
        for i in color:
            color[n] = frame_colors[i]
            n += 1

        fig.add_trace(
            go.Bar(
                x=list(metagene_dict.keys()),
                y=list(metagene_dict.values()),
                name="Distance from 5'"
                if current_target == "start"
                else "Distance from 3'",
                marker=dict(color=color),
            ),
            row=1,
            col=count,
        )

    fig.update_layout(
        title="Metagene Profile",
        xaxis_title="Relative position",
        yaxis_title="Read Count",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
        bargap=0,
        showlegend=False,
    )
    if colnums > 1:
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
        "description": "Metagene profile showing the distance count of \
reads per distance away from a target (default: start codon).",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_metagene_profile_dict


def plot_metagene_heatmap(metagene_profile_dict: dict, config: dict) -> dict:
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
    colnums = 1
    if metagene_profile_dict["start"] != {}:
        target_loop = ["start"]
        if metagene_profile_dict["stop"] != {}:
            target_loop.append("stop")
            colnums = 2
    else:
        target_loop = ["stop"]

    fig = make_subplots(
        rows=1,
        cols=colnums,
        shared_yaxes=True,
        subplot_titles=["Distance from 5'", "Distance from 3'"]
        if len(target_loop) > 1
        else ["Distance from 5'"]
        if target_loop == ["start"]
        else ["Distance from 3'"],
    )

    for current_target in target_loop:
        count += 1
        x_data = []
        y_data = []
        z_data = []
        for k1, v1 in metagene_profile_dict[current_target].items():
            for k2, v2 in v1.items():
                x_data.append(k2)
                y_data.append(k1)
                z_data.append(v2)

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
        range=config["plots"]["metagene_profile"]["distance_range"]
    )
    fig.update_layout(
        title="Metagene Heatmap",
        xaxis_title="Relative position (nt)",
        yaxis_title="Read length",
        font=dict(
            family=config["plots"]["font_family"],
            size=18,
            color=config["plots"]["base_color"],
        ),
        legend={"traceorder": "normal"},
        showlegend=False,
    )
    if colnums > 1:
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
    plot_metagene_heatmap = {
        "name": "Metagene Heatmap",
        "description": "Metagene heatmap showing the distance between the \
            A-site and a target per read length and the counts in colorscale.",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_metagene_heatmap


def plot_logoplot(sequence_slice_dict: dict, config: dict) -> dict:
    """
    Generate a logoplot of the nucleotide distribution of the reads

    Inputs:
        sequence_slice_dict: Dictionary containing the sequences as values
            and read length as keys
        config: Dictionary containing the configuration information

    Outputs:
        plot_logoplot_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    # nt_start, nt_count = (
    #     config["plots"]["nucleotide_proportion"]["nucleotide_start"],
    #     config["plots"]["nucleotide_proportion"]["nucleotide_count"],
    # )
    with tempfile.TemporaryDirectory() as tempdir:
        with open(f"{tempdir}/temp_fasta.fasta", "w+") as fasta:
            count = 0
            for n in sequence_slice_dict.values():
                count += 1
                fasta.write(f">{count}\n{n}\n")
            with open(f"{tempdir}/temp_plot.png", "w+b") as plot:
                weblogo_prompt = [
                    "weblogo",
                    f"-f{fasta.name}",
                    "-Dfasta",
                    f"-o{plot.name}",
                    "-Fpng_print",
                ]
                subprocess.run(weblogo_prompt)
                plot.seek(0)
                fig_image = base64.b64encode(plot.read()).decode("utf-8")

    plot_logoplot_dict = {
        "name": "Logoplot",
        "description": "Logoplot created with WebLogo ver 3.7.12",
        "fig_html": f'<img src="data:image/png;base64,{fig_image}" \
            alt="Logoplot created with WebLogo" width="100%" height="auto"">',
        "fig_image": fig_image,
    }
    return plot_logoplot_dict

# WIP
def plot_metrics_summary(metrics_dict: dict, config: dict) -> dict:
    width = 450
    height = 400
    fig = go.Figure(data=go.Scatterpolar(
    r=[0.8, .5, .2, .6, .3, 0.8],
    theta=['read frame periodicity','ligation bias','nucleotide composition', 'mRNA breakdown',
            'metagene profile', 'read frame periodicity'],
    fill='toself'
    ))
        
    fig.update_layout(
    title="[WIP] Summary Scores",
    # autosize=False,
    # width=width,
    # height=height,
    polar=dict(
        radialaxis=dict(
        visible=False,
        ),
        radialaxis_range = [0, 1],
        angularaxis_thetaunit = "radians"
    ),
    showlegend=False,
    )
    plot_metrics_summary_dict = {
        "name": "Summary of Metrics",
        "description": "Radar plot showing the difference scores ranging \
from 0 to 1.",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     width,
                                     height),
    }
    return plot_metrics_summary_dict
