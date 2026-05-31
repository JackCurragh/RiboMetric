
from RiboMetric.file_parser import parse_bam
from RiboMetric.bam_processing import process_reads
import pandas as pd


def test_bam_parsing(test_data_dir):
    """Test bam parsing"""
    bam = parse_bam(
        bam_file=str(test_data_dir / "test.bam"),
        num_reads=10000,
        num_processes=1,
    )[0]
    assert len(bam) == 9997


def test_process_reads_accepts_query_name_alias():
    oxbow_df = pd.DataFrame({
        "query_name": ["read1_x3"],
        "seq": ["ACGT"],
        "cigar": ["4M"],
        "rname": ["tx1"],
        "pos": [10],
        "mapq": [255],
    })

    reads = process_reads(oxbow_df)

    assert reads["read_name"].iloc[0] == "read1_x3"
    assert reads["count"].iloc[0] == 3


def test_process_reads_empty_schema_returns_typed_empty_batch():
    reads = process_reads(pd.DataFrame())

    assert reads.empty
    assert "read_name" in reads.columns
    assert "reference_name" in reads.columns
