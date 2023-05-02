#!/usr/bin/env python

"""Tests for `RibosomeProfiler` package."""

from RibosomeProfiler.file_parser import parse_bam


def test_bam_parsing():
    """Test bam parsing"""
    bam = parse_bam(
        "tests/test_data/test.bam", 10000
    )
    assert len(bam) == 10001
