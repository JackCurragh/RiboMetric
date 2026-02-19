
from RiboMetric.file_parser import parse_bam


def test_bam_parsing(test_data_dir):
    """Test bam parsing"""
    bam = parse_bam(
        bam_file=str(test_data_dir / "test.bam"),
        num_reads=10000,
        num_processes=1,
    )[0]
    assert len(bam) == 9997
