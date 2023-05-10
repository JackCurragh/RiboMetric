"""
This script contains the code for generating the plots for
RibosomeProfiler reports
"""

from plotly import graph_objects as go
from .modules import read_frame_cull, read_frame_score, sum_mRNA_distribution
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
    plots_list = []
    plots_list.extend([
        plot_read_length_distribution(
            results_dict["read_length_distribution"],
            config),
        plot_ligation_bias_distribution(
            results_dict["ligation_bias_distribution"],
        config),
        plot_nucleotide_composition(
            results_dict["nucleotide_composition"],
            config),
        plot_read_frame_distribution(
            results_dict["read_frame_distribution"],
            config)
    ])
    if results_dict["mode"] == "annotation_mode":
        plots_list.extend([
            plot_mRNA_distribution(
                results_dict["mRNA_distribution"],
                config),
            plot_mRNA_read_breakdown(
                results_dict["mRNA_distribution"],
                config),

        ])
    return plots_list


def plot_read_length_distribution(read_length_dict: dict, config: dict
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
        "fig_image": base64.b64encode(pio.to_image(fig, format="jpg")
                                      ).decode("ascii"),
    }
    return plot_read_length_dict


def plot_ligation_bias_distribution(ligation_bias_dict: dict, config: dict
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
    if config["plots"]["ligation_bias_distribution"]["include_N"] is False:
        ligation_bias_dict = {
            k: v for k, v in ligation_bias_dict.items() if "N" not in k
        }
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=list(ligation_bias_dict.keys()),
            y=list(ligation_bias_dict.values()),
            name="",
            hovertemplate="<b>Nucleotides</b>:%{x}<br><b>Proportion</b>:%{y}",
        )
    )
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
        "fig_image": base64.b64encode(pio.to_image(fig, format="jpg")
                                      ).decode("ascii"),
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
            go.Scatter(y=distribution, name=nucleotide,
                       line_color=colors[nucleotide])
        )
    fig.update_layout(
        title="Nucleotide Distribution",
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
        "fig_image": base64.b64encode(pio.to_image(fig, format="jpg")
                                      ).decode("ascii"),
    }
    return plot_nucleotide_composition_dict


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
    scored_read_frame_dict = read_frame_score(culled_read_frame_dict) \
        if config["plots"]["read_frame_distribution"]["show_scores"] != "none" \
        else None
    
    # Set minimum and maximum font sizes
    min_font_size, max_font_size = 5, 30

    # Calculate font size based on number of data points
    num_data_points = len(culled_read_frame_dict)
    font_size = max_font_size - (max_font_size - min_font_size) * (num_data_points / 50)

    plot_data = []
    for i in range(0, 3):
        plot_data.append(
            go.Bar(
                name="Frame " + str(i + 1),
                x=list(culled_read_frame_dict.keys()),
                y=[culled_read_frame_dict[x][y]
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
    if scored_read_frame_dict != None:
        if config["plots"]["read_frame_distribution"]["show_scores"] == "all":
            for idx in enumerate(culled_read_frame_dict):
                if idx[1] != "global":
                    y_buffer = max(fig.data[0].y+
                                fig.data[1].y+
                                fig.data[2].y)*0.05
                    ymax = max(fig.data[0].y[idx[0]],
                                fig.data[1].y[idx[0]],
                                fig.data[2].y[idx[0]])
                    if fig.data[0].y[idx[0]]\
                          + fig.data[0].y[idx[0]]\
                          + fig.data[0].y[idx[0]]\
                              > y_buffer:
                        fig.add_annotation(
                            x=idx[1],
                            y=ymax+y_buffer,
                            text=round(scored_read_frame_dict[idx[1]],2),
                            showarrow=False,
                            xanchor='center',
                            font={"size":font_size}
                        )
        fig.update_layout()
        fig.add_annotation(text=f'Score: {round(scored_read_frame_dict["global"],2)}', 
                        showarrow=False,
                        xref='paper',
                        yref='paper',
                        y=0.64,
                        x=1.03,
                        xanchor="left")
    plot_read_frame_dict = {
        "name": "Read Frame Distribution",
        "description": "Frame distribution per read length",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": base64.b64encode(pio.to_image(fig, format="jpg")
                                      ).decode("ascii"),
    }
    return plot_read_frame_dict

#WIP
def plot_frame_score_distribution(read_frame_dict: dict, config: dict) -> dict:
    """
    Generate a plot of the read frame score distribution

    Inputs:
        read_frame_dict: Dataframe containing the read frame distribution
        config: Dictionary containing the configuration information

    Outputs:
        plot_frame_score_dict: Dictionary containing the plot name, description
        and plotly figure for html and pdf export
    """
    culled_read_frame_dict = read_frame_cull(read_frame_dict, config)
    scored_read_frame_dict = read_frame_score(culled_read_frame_dict)


def plot_mRNA_distribution(mRNA_distribution_dict: dict, config: dict) -> dict:
    """
    Generate a bar plot of the mRNA distribution

    Inputs:
        mRNA_distribution_dict: Dataframe containing the mRNA distribution
        over the read lengths
        config: Dictionary containing the configuration information

    Outputs:
        plot_mRNA_distribution_dict: Dictionary containing the plot name,
        description and plotly figure for html and pdf export
    """
    sum_mRNA_dict = sum_mRNA_distribution(mRNA_distribution_dict, config)
    plot_data = []
    for k,v in sum_mRNA_dict.items():
        plot_data.append(
            go.Bar(
                name=k.replace("_"," ").title(),
                x=[v],
                y=[""],
                width=[0.3],
                hovertemplate = "Proportion: %{x:.2%}"
                if not config["plots"]["mRNA_distribution"]["absolute_counts"]
                else "Count: %{x}",
                orientation='h'
                )
            )

    fig = go.Figure(plot_data)
    fig.update_layout(
        barmode='stack',
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
            legend={'traceorder':'normal'},
        )
    plot_mRNA_distribution_dict = {
        "name": "mRNA Reads Breakdown",
        "description": "Shows the proportion of the different transcript \
regions represented in the reads",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": base64.b64encode(pio.to_image(fig, format="jpg")
                                      ).decode("ascii"),
    }
    return plot_mRNA_distribution_dict


def plot_mRNA_read_breakdown(mRNA_distribution_dict: dict, config: dict) -> dict:
    """
    Generate a line plot of the mRNA distribution over the read lenghts

    Inputs:
        mRNA_distribution_dict: Dataframe containing the mRNA distribution
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
            plot_data = {k: [x/sum(sum_data.values()) for x in v] for k, v in plot_data.items()}
    fig = go.Figure()
    for k,v in plot_data.items(): 
        fig.add_trace(
                go.Scatter(
                    name = k,
                    x=list(mRNA_distribution_dict.keys()),
                    y=v,
                    hovertemplate = "Proportion: %{y:.2%}"
                    if not config["plots"]["mRNA_read_breakdown"]["absolute_counts"]
                    else "Count: %{x}",
                ))
    fig.update_layout(
            title="Nucleotide Distribution",
            xaxis_title="Position (nucleotides)",
            yaxis_title="Proportion"
            if not config["plots"]["mRNA_read_breakdown"]["absolute_counts"]
            else "Counts",
            #yaxis_range=[0, 1],
            font=dict(
                family=config["plots"]["font_family"],
                size=18,
                color=config["plots"]["base_color"],
            ),
            legend={'traceorder':'normal'},
        )
    plot_mRNA_read_breakdown_dict = {
        "name": "mRNA Reads Breakdown over Read Length",
        "description": "Shows the proportion of the different transcript \
regions represented in the reads over the different read lengths.",
        "fig_html": pio.to_html(fig, full_html=False),
        "fig_image": base64.b64encode(pio.to_image(fig, format="jpg")
                                      ).decode("ascii"),
    }
    return plot_mRNA_read_breakdown_dict