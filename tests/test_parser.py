
from RiboMetric.file_parser import parse_bam
from RiboMetric.bam_processing import process_reads
from RiboMetric.file_splitting import split_idxstats_df
import pandas as pd


def test_bam_parsing(test_data_dir):
    """Test bam parsing"""
    bam = parse_bam(
        bam_file=str(test_data_dir / "test.bam"),
        num_reads=10000,
        num_processes=1,
    )[0]
    assert 0 < len(bam) <= 10200


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


def test_split_idxstats_keeps_oversized_first_reference():
    idxstats_df = pd.DataFrame({
        "Reference": ["tx_empty", "tx_big", "tx_small"],
        "Length": ["50", "1000", "100"],
        "Mapped_Reads": ["0", "200", "1"],
        "Unmapped_Reads": ["0", "0", "0"],
    })

    batches = split_idxstats_df(idxstats_df, batch_size=100, num_reads=100)

    assert len(batches) == 1
    assert len(batches[0]) == 2
    assert batches[0]["Reference"].tolist() == ["tx_empty", "tx_big"]
