"""
This script contains the code for generating the plots for
RiboMetric reports
"""

from plotly import graph_objects as go
from plotly.subplots import make_subplots
from .modules import read_frame_cull, read_frame_score_trips_viz, sum_mRNA_distribution
from .results_output import normalise_score
import plotly.io as pio
import base64

import plotly.express as px
import pandas as pd


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
    plots_list = [plot_metrics_summary(results_dict["metrics"].copy(), config)]
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
                ),
                plot_metagene_heatmap(
                    results_dict["metagene_profile"],
                    config
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
    if terminal_nucleotide_bias_dict["five_prime"] != {}:
        target_loop = ["five_prime"]
        if terminal_nucleotide_bias_dict["three_prime"] != {}:
            target_loop.append("three_prime")
            columns = 2
    else:
        target_loop = ["three_prime"]

    fig = make_subplots(
        rows=1,
        cols=columns,
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
        "fig_html": pio.to_html(fig, full_html=False),
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
    fig.update_xaxes(range=[lower_limit-0.5, upper_limit+0.5])
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
        color = [(int(x) % 3) for x in metagene_dict.keys()]
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
        x_data = []
        y_data = []
        z_data = []
        for read_length, position_counts in metagene_profile_dict[
                                                current_target
                                                ].items():
            for position, counts in position_counts.items():
                x_data.append(int(position))
                y_data.append(int(read_length))
                z_data.append(int(counts))

        if config["plots"]["metagene_profile"]["max_colorscale"] is None:
            z_max = max(z_data)
            z_data = [z/z_max for z in z_data]

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
        "description": "Metagene heatmap showing the distance between the \
            A-site and a target per read length and the counts in colorscale.",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }
    return plot_metagene_heatmap


def plot_metrics_summary(metrics_dict: dict, config: dict) -> dict:
    '''
    generate the metrics summary plot that is a bar chart of the metrics
    found at the top of the report

    Inputs:
        metrics_dict: Dictionary containing the metrics and their scores
        config: Dictionary containing the configuration information

    Outputs:
        plot_metrics_summary_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    '''
    # for any entry in metrics_dict with a dict as value (these are where the
    # metrics are per read length) get average of the top 3 values of that dict
    for key, value in metrics_dict.items():
        if isinstance(value, dict):
            metrics_dict[key] = sum(sorted(value.values(), reverse=True)[:3])/3

    # Convert the metrics_dict to a DataFrame for easier plotting
    df = pd.DataFrame(list(metrics_dict.items()), columns=['Metric', 'Score'])

    # drop any metric that is not in max_mins
    df = df[df['Metric'].apply(lambda x: any(x.startswith(key) for key in config["max_mins"].keys()))]
    # normalise metrics between max_min values
    for metric in config["max_mins"]:
        df.loc[df['Metric'].str.startswith(metric), 'Score'] = normalise_score(
            df.loc[df['Metric'].str.startswith(metric), 'Score'].values[0],
            config["max_mins"][metric][0],
            config["max_mins"][metric][1]
        )
    # drop any metrics that are in exclude list in config
    df = df[~df['Metric'].isin(config["plots"]["exclude_metrics"])]

    width = 750
    height = 320
    # Create a bar chart using Plotly Express
    fig = px.bar(df, x='Metric', y='Score', title="Summary Scores",
                 labels={'Score': 'Score'},
                 hover_data={'Metric': True, 'Score': ':.3f'},
                 height=400)

    # Customize the layout if needed
    fig.update_layout(showlegend=False)

    # Convert the plot to HTML or an image as needed
    fig_html = pio.to_html(fig, full_html=False)
    fig_image = plotly_to_image(fig, width, height)

    plot_metrics_summary_dict = {
        "plot": {
            "name": "Summary of Metrics",
            "description": "Bar chart showing the difference scores ranging from 0 to 1.",
            "fig_html": fig_html,
            "fig_image": fig_image,
        },
        "metrics": [{"name": k.replace("_", " ").capitalize(), "score": round(v, 3)}
                    for k, v in metrics_dict.items() if isinstance(v, float)]
    }

    return plot_metrics_summary_dict


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
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": plotly_to_image(fig,
                                     config["plots"]["image_size"][0],
                                     config["plots"]["image_size"][1]),
    }

    return plot_read_frame_triangle_dict
