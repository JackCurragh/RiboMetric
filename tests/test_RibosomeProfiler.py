#!/usr/bin/env python

"""Tests for `RibosomeProfiler` package."""

from RibosomeProfiler.file_parser import parse_bam


def test_bam_parsing():
    """Test bam parsing"""
    bam = parse_bam(
        "tests/test_data/test.bam", 1000000
    )
    assert bam["total_reads"] == 1000000
    assert bam["mapped_reads"] == 1000000
    assert bam["unmapped_reads"] == 0
