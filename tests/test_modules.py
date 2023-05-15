from RibosomeProfiler.modules import (
    read_length_distribution,
    ligation_bias_distribution,
    nucleotide_composition,
    a_site_calculation,
    read_frame_distribution,
)
import pandas as pd


def test_read_length_distribution():
    """Test read length distribution calculation"""
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    read_length_dict = read_length_distribution(read_df)
    assert read_length_dict[29] == 4


def test_ligation_bias_distribution():
    """Test ligation bias distribution calculation"""
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    ligation_bias_dict = ligation_bias_distribution(read_df)
    assert ligation_bias_dict["AA"] == 0.5


def test_nucleotide_composition():
    """Test nucleotide composition calculation"""
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
    ]


def test_a_site_calculation():
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    a_site_df = a_site_calculation(read_df)
    assert a_site_df.a_site[0] == 366


def test_read_frame_distribution():
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    read_frame_dict = read_frame_distribution(a_site_calculation(read_df))
    assert read_frame_dict[33][0] == 1
