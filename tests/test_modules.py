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
    annotate_reads,
    assign_mRNA_category,
    mRNA_distribution,
    metagene_profile,
)
import pandas as pd
import pytest


def test_a_site_calculation():
    """
    Test A-site calculation
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    a_site_df = a_site_calculation(read_df)
    assert a_site_df.a_site[0] == 353


def test_read_length_distribution():
    """
    Test read length distribution calculation
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    read_length_dict = read_length_distribution(read_df)
    assert read_length_dict[29] == 40


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ('AA_test', 0.5),
        ('ligation_bias_dict_norm["three_prime"]["TC"]', 0.275),
    ],
)
def test_ligation_bias_distribution(test_input, expected):
    """
    Test ligation bias distribution calculation
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    categories = ["first_dinucleotide", "last_dinucleotide"]
    read_df[categories] = read_df[categories].astype("category")
    sequence_background = {
        "5_prime_bg": {"AA": 0.25, "AG": 0.3, "TT": 0.45},
        "3_prime_bg": {"AA": 0.875, "TC": 0.125}}
    ligation_bias_dict = ligation_bias_distribution(read_df)
    AA_test = ligation_bias_dict["five_prime"]["AA"]
    ligation_bias_dict_norm = normalise_ligation_bias(
        ligation_bias_dict, sequence_background
    )

    # Just to satisfy the linter
    type(AA_test)
    type(ligation_bias_dict_norm)
    assert eval(test_input) == expected


def test_nucleotide_composition():
    """
    Test nucleotide composition calculation
    """
    sequence_data = {1: {
        "A": [4, 4, 5, 5, 1, 0, 0, 3, 3, 3, 2],
        "C": [1, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0],
        "G": [3, 3, 0, 0, 4, 5, 5, 0, 0, 0, 0],
        "T": [0, 1, 3, 0, 0, 0, 0, 5, 5, 5, 0]}}
    nucleotide_composition_dict = nucleotide_composition(sequence_data[1])
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
        1,
    ]


def test_read_frame_distribution():
    """
    Test read frame labelling
    """
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    read_frame_dict = read_frame_distribution(a_site_calculation(read_df))
    assert read_frame_dict[33][1] == 10


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ('mRNA_distribution_dict[29]["CDS"]', 10),
        ('mRNA_distribution_dict[21]["three_trailer"]', 40),
    ],
)
def test_mRNA_distribution(test_input, expected):
    """
    Test metagene distance calculations
    """
    annotation_df = pd.read_csv(
        "tests/test_data/test_annotation.tsv", sep="\t",
        dtype={"transcript_id": str, "cds_start": int,
               "cds_end": int, "transcript_length": int}
    )
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    a_site_df = a_site_calculation(read_df)
    annotated_read_df = annotate_reads(a_site_df, annotation_df)
    annotated_read_df = assign_mRNA_category(annotated_read_df)

    mRNA_distribution_dict = mRNA_distribution(annotated_read_df)

    # Just to satisfy the linter
    type(mRNA_distribution_dict)
    assert eval(test_input) == expected


def test_metagene_profile():
    """
    Test metagene distance calculations
    """
    annotation_df = pd.read_csv(
        "tests/test_data/test_annotation.tsv", sep="\t",
        dtype={"transcript_id": str, "cds_start": int,
               "cds_end": int, "transcript_length": int}
    )
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[
        read_df_pre.index.repeat(read_df_pre["count"])
    ].reset_index(drop=True)
    a_site_df = a_site_calculation(read_df)
    annotated_read_df = annotate_reads(a_site_df, annotation_df)
    metagene_profile_dict = metagene_profile(annotated_read_df)
    assert metagene_profile_dict["stop"][21][1] == 30
