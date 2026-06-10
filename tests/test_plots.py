#!/usr/bin/env python

"""Tests for `plots` package."""

from RiboMetric.plots import (
    plot_terminal_nucleotide_bias_distribution,
    plot_nucleotide_composition,
    plot_read_length_distribution,
    plot_read_frame_distribution,
)
from RiboMetric.modules import (
    read_length_distribution,
    terminal_nucleotide_bias_distribution,
    nucleotide_composition,
    a_site_calculation,
    read_frame_distribution,
)
import yaml
import pandas as pd


def is_plotly_html(fig_html):
    return fig_html.lstrip().startswith("<div") and "Plotly.newPlot" in fig_html


def assert_plotly_html(fig_html):
    assert is_plotly_html(fig_html)


def test_plot_read_length_distribution():
    with open("RiboMetric/config.yml", "r") as ymlfile:
        config = yaml.load(ymlfile, Loader=yaml.Loader)
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    plot_read_length = plot_read_length_distribution(
        read_length_distribution(read_df), config
    )
    assert_plotly_html(plot_read_length["fig_html"])


def test_plots():
    errors = []
    with open("RiboMetric/config.yml", "r") as ymlfile:
        config = yaml.load(ymlfile, Loader=yaml.Loader)
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    categories = ["first_dinucleotide", "last_dinucleotide"]
    read_df[categories] = read_df[categories].astype("category")
    sequence_data_single = {"A": [5, 2, 3, 1], "T": [0, 2, 1, 2]}
    if not is_plotly_html(plot_read_length_distribution(
            read_length_distribution(read_df),
            config)["fig_html"]):
        errors.append("Read length distribution plot html output error")
    if not is_plotly_html(plot_terminal_nucleotide_bias_distribution(
            terminal_nucleotide_bias_distribution(read_df), config)["fig_html"]):
        errors.append("Ligation bias distribution plot html output error")
    if not is_plotly_html(plot_nucleotide_composition(
            nucleotide_composition(sequence_data_single), config)["fig_html"]):
        errors.append("Nucleotide composition plot html output error")
    if not is_plotly_html(plot_read_frame_distribution(
            read_frame_distribution(
                a_site_calculation(read_df)
            ),
            config)["fig_html"]):
        errors.append("Read frame distribution plot html output error")
    assert not errors, "errors occurred:\n{}".format("\n".join(errors))


def test_plot_read_frame_distribution_empty():
    """An empty frame dict (e.g. no reads survived CDS filtering) must not
    crash the plot with an UnboundLocalError on the x-axis range."""
    with open("RiboMetric/config.yml", "r") as ymlfile:
        config = yaml.load(ymlfile, Loader=yaml.Loader)
    result = plot_read_frame_distribution({}, config)
    assert_plotly_html(result["fig_html"])


def test_plot_read_frame_distribution_all_below_buffer():
    """Frames present but all counts tiny/zero must also not crash."""
    with open("RiboMetric/config.yml", "r") as ymlfile:
        config = yaml.load(ymlfile, Loader=yaml.Loader)
    frame_dict = {29: {0: 0, 1: 0, 2: 0}, 30: {0: 0, 1: 0, 2: 0}}
    result = plot_read_frame_distribution(frame_dict, config)
    assert_plotly_html(result["fig_html"])
