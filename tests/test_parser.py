
from RiboMetric.file_parser import parse_bam
import pandas as pd

def test_bam_parsing():
    """Test bam parsing"""
    bam = pd.concat(parse_bam(
        "tests/test_data/test.bam", 10000
    ))
    assert len(bam) == 10000
