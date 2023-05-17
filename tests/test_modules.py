"""
This script contains tests for the different functions found in modules.py
"""

from RiboMetric.modules import (
    read_length_distribution,
    ligation_bias_distribution,
    normalise_ligation_bias,
    nucleotide_composition,
    a_site_calculation,
    read_frame_distribution,
    metagene_profile,
)
import pandas as pd
import pytest


def test_read_length_distribution():
    """
    Test read length distribution calculation
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    read_length_dict = read_length_distribution(read_df)
    assert read_length_dict[29] == 4


@pytest.mark.parametrize("test_input,expected", [('ligation_bias_dict_2nt["AA"]', 0.25), ('ligation_bias_dict_3nt["AAA"]', 0.25)])
def test_ligation_bias_distribution(test_input, expected):
    """
    Test ligation bias distribution calculation
    """
    errors = []
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    ligation_bias_dict_2nt = ligation_bias_distribution(read_df)
    ligation_bias_dict_2nt = normalise_ligation_bias(read_df,
                                                    ligation_bias_dict_2nt)
    ligation_bias_dict_3nt = ligation_bias_distribution(read_df,
                                                    num_bases=3,
                                                    five_prime=False)
    ligation_bias_dict_3nt = normalise_ligation_bias(read_df,
                                                    ligation_bias_dict_3nt,
                                                    num_bases=3,
                                                    five_prime=False)
    assert eval(test_input) == expected


def test_nucleotide_composition():
    """
    Test nucleotide composition calculation
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    nucleotide_composition_dict = nucleotide_composition(read_df)
    assert nucleotide_composition_dict["A"] == [
        0.5,
        0.5,
        0.625,
        0.625,
        0.125,
        0,
        0,
        0.375,
        0.375,
        0.375,
        1
    ]


def test_a_site_calculation():
    """
    Test A-site calculation
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    a_site_df = a_site_calculation(read_df)
    assert a_site_df.a_site[0] == 366


def test_read_frame_distribution():
    """
    Test read frame labelling 
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    read_frame_dict = read_frame_distribution(a_site_calculation(read_df))
    assert read_frame_dict[33][0] == 1


def test_metagene_profile():
    """
    Test metagene distance calculations
    """
    # metagene_profile(annotated_read_df)