'''
This script contains the code for generating the plots for RibosomeProfiler reports 

'''

import plotly


def plot_read_length_distribution(read_length_dict: dict, config: dict) -> plotly.graph_objects.Figure:
    '''
    Generate a plot of the read length distribution for the full dataset

    Inputs:
        read_length_df: Dataframe containing the read length distribution
        config: Dictionary containing the configuration information

    Outputs:
        fig: Plotly figure containing the read length distribution
    '''
    fig = plotly.graph_objects.Figure()
    fig.add_trace(plotly.graph_objects.Bar(
        x=list(read_length_dict.keys()),
        y=list(read_length_dict.values()),
        name='',
        hovertemplate =
        '<b>Read length</b>: %{x}'+
        '<br><b>Count</b>: %{y}',
    ))
    fig.update_layout(
        title='Read Length Distribution',
        xaxis_title='Read Length',
        yaxis_title='Read Count',
        font=dict(
            family='Courier New, monospace',
            size=18,
            color='#7f7f7f'
        ))
    return fig

