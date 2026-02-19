"""
This script contains tests for the different functions found in modules.py
"""

from RiboMetric.modules import (
    read_length_distribution,
    terminal_nucleotide_bias_distribution,
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


def test_a_site_calculation(sample_read_df):
    """
    Test A-site calculation
    """
    a_site_df = a_site_calculation(sample_read_df)
    assert a_site_df.a_site[0] == 356


def test_read_length_distribution(sample_read_df):
    """
    Test read length distribution calculation
    """
    read_length_dict = read_length_distribution(sample_read_df)
    assert read_length_dict[29] == 40


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ('AA_test', 0.5),
        ('terminal_nucleotide_bias_dict_norm["three_prime"]["TC"]', 0.275),
    ],
)
def test_terminal_nucleotide_bias_distribution(test_input, expected, sample_read_df):
    """
    Test ligation bias distribution calculation
    """
    read_df = sample_read_df.copy()
    categories = ["first_dinucleotide", "last_dinucleotide"]
    read_df[categories] = read_df[categories].astype("category")
    sequence_background = {
        "5_prime_bg": {"AA": 0.25, "AG": 0.3, "TT": 0.45},
        "3_prime_bg": {"AA": 0.875, "TC": 0.125}}
    terminal_nucleotide_bias_dict = terminal_nucleotide_bias_distribution(read_df)
    AA_test = terminal_nucleotide_bias_dict["five_prime"]["AA"]
    terminal_nucleotide_bias_dict_norm = normalise_ligation_bias(
        terminal_nucleotide_bias_dict, sequence_background
    )

    # Just to satisfy the linter
    type(AA_test)
    type(terminal_nucleotide_bias_dict_norm)
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


def test_read_frame_distribution(sample_read_df):
    """
    Test read frame labelling
    """
    read_frame_dict = read_frame_distribution(a_site_calculation(sample_read_df))
    assert read_frame_dict[33][1] == 10


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ('mRNA_distribution_dict[29]["CDS"]', 30),
        ('mRNA_distribution_dict[21]["three_trailer"]', 40),
    ],
)
def test_mRNA_distribution(test_input, expected, sample_read_df, sample_annotation_df):
    """
    Test metagene distance calculations
    """
    a_site_df = a_site_calculation(sample_read_df)
    annotated_read_df = annotate_reads(a_site_df, sample_annotation_df)
    annotated_read_df = assign_mRNA_category(annotated_read_df)

    mRNA_distribution_dict = mRNA_distribution(annotated_read_df)

    # Just to satisfy the linter
    type(mRNA_distribution_dict)
    assert eval(test_input) == expected


def test_metagene_profile(sample_read_df, sample_annotation_df):
    """
    Test metagene distance calculations
    """
    a_site_df = a_site_calculation(sample_read_df)
    annotated_read_df = annotate_reads(a_site_df, sample_annotation_df)
    metagene_profile_dict = metagene_profile(annotated_read_df)

    assert metagene_profile_dict["stop"][21][4] == 30
