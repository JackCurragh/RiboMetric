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
    read_df = pd.read_csv("tests/test_data/parsed_test.csv")
    read_length_dict = read_length_distribution(read_df)
    assert read_length_dict[16] == 41
