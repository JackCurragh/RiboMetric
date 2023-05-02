#!/usr/bin/env python

"""Tests for `RibosomeProfiler` package."""

from RibosomeProfiler.file_parser import parse_bam
from RibosomeProfiler.modules import *

def test_bam_parsing():
    """Test bam parsing"""
    bam = parse_bam(
        "tests/test_data/test.bam", 10000
    )
    assert len(bam) == 10001

def test_read_length_distribution():
    """Test read length distribution calculation"""
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[read_df_pre.index.repeat(read_df_pre['count'])].reset_index(drop=True)
    read_length_dict = read_length_distribution(read_df)
    assert read_length_dict[29] == 4

def test_ligation_bias_distribution():
    """Test ligation bias distribution calculation"""
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[read_df_pre.index.repeat(read_df_pre['count'])].reset_index(drop=True)
    ligation_bias_dict = ligation_bias_distribution(read_df)
    assert ligation_bias_dict["AA"] == 0.5

def test_nucleotide_composition():
    """Test nucleotide composition calculation"""
    read_df_pre = pd.read_csv("tests/test_data/test.csv")
    read_df = read_df_pre.loc[read_df_pre.index.repeat(read_df_pre['count'])].reset_index(drop=True)
    nucleotide_composition_dict = nucleotide_composition(read_df)
    assert nucleotide_composition_dict["A"] == [0.5,0.5,0.625,0.625,0.125,0,0,0.375,0.375,0.375]