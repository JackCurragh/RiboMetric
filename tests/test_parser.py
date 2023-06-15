
from RiboMetric.file_parser import parse_bam
import pandas as pd

def test_bam_parsing():
    """Test bam parsing"""
    bam = parse_bam(
        bam_file="tests/test_data/test.bam",
        num_reads=10000,
        num_processes=1,
    )[0]
    assert len(bam) == 9997
