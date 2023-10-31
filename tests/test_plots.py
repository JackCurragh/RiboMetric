#!/usr/bin/env python

"""Tests for `plots` package."""

from RiboMetric.plots import (
    plot_ligation_bias_distribution,
    plot_nucleotide_composition,
    plot_read_length_distribution,
    plot_read_frame_distribution,
)
from RiboMetric.modules import (
    read_length_distribution,
    ligation_bias_distribution,
    nucleotide_composition,
    a_site_calculation,
    read_frame_distribution,
)
import yaml
import pandas as pd


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
    assert "<div>" in plot_read_length["fig_html"]


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
    if (
        "<div>" not in plot_read_length_distribution(
            read_length_distribution(read_df),
            config)["fig_html"]
    ):
        errors.append("Read length distribution plot html output error")
    if ("<div>" not in plot_ligation_bias_distribution(
            ligation_bias_distribution(read_df), config)["fig_html"]):
        errors.append("Ligation bias distribution plot html output error")
    if ("<div>" not in plot_nucleotide_composition(
            nucleotide_composition(sequence_data_single), config)["fig_html"]):
        errors.append("Nucleotide composition plot html output error")
    if ("<div>" not in plot_read_frame_distribution(
            read_frame_distribution(
                a_site_calculation(read_df)
            ),
            config)["fig_html"]):
        errors.append("Read frame distribution plot html output error")
    assert not errors, "errors occurred:\n{}".format("\n".join(errors))
