"""
This script contains the code for generating the plots for
RibosomeProfiler reports
"""

from plotly import graph_objects as go
import plotly.io as pio

def plot_read_length_distribution(
    read_length_dict: dict, config: dict
) -> dict:
    """
    Generate a plot of the read length distribution for the full dataset

    Inputs:
        read_length_df: Dataframe containing the read length distribution
        config: Dictionary containing the configuration information

    Outputs:
        plot_read_length_dict: Dictionary containing the plot name, description and plotly figure in html for the read length distribution
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
            family="Helvetica Neue,Helvetica,Arial,sans-serif",
            size=18,
            color="#7f7f7f"
        ),
    )
    plot_read_length_dict = {
        "name": "Read Length Distribution",
        "description": "A plot showcasing the read length distribution of the reads",
        "fig": pio.to_html(fig, full_html=False)
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
        plot_ligation_bias_dict: Dictionary containing the plot name, description and plotly figure in html for the ligation bias distribution
    """
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=list(ligation_bias_dict.keys()),
            y=list(ligation_bias_dict.values()),
            name="",
            hovertemplate="<b>Nucleotides</b>: %{x}" + "<br><b>Proportion</b>: %{y}",
        )
    )
    fig.update_layout(
        title="Ligation Bias Distribution",
        xaxis_title="Read Start",
        yaxis_title="Proportion",
        font=dict(
            family="Helvetica Neue,Helvetica,Arial,sans-serif", size=18, color="#7f7f7f"
        ),
    )
    plot_ligation_bias_dict = {
        "name": "Ligation Bias Distribution",
        "description": "A plot showcasing the ligation bias distribution of the reads",
        "fig": pio.to_html(fig, full_html=False)
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
        plot_nucleotide_composition_dict: Dictionary containing the plot name, description and plotly figure in html for the nucleotide composition
    """
    colors = {"A": "#c93434", "C": "#2e85db", "G": "#f0de1f", "T": "#1fc24d"}
    fig = go.Figure()
    for nucleotide, distribution in nucleotide_composition_dict.items():
        fig.add_trace(
            go.Scatter(y=distribution, name=nucleotide, line_color=colors[nucleotide])
        )
    fig.update_layout(
        title="Nucleotide Distribution",
        xaxis_title="Position (nucleotides)",
        yaxis_title="Proportion",
        yaxis_range=[0, 1],
        font=dict(
            family="Helvetica Neue,Helvetica,Arial,sans-serif", size=18, color="#7f7f7f"
        ),
    )
    plot_nucleotide_composition_dict = {
        "name": "Nucleotide Composition",
        "description": "A plot showcasing the nucleotide composition of the reads",
        "fig": pio.to_html(fig, full_html=False)
    }
    return plot_nucleotide_composition_dict
