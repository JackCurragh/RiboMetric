import pandas as pd
from RiboMetric.file_parser import check_annotation, parse_annotation


def test_check_annotation_minimal(tmp_path):
    p = tmp_path / "anno.tsv"
    p.write_text(
        "transcript_id\tcds_start\tcds_end\ttranscript_length\nTX1\t0\t90\t90\n"
    )
    assert check_annotation(str(p)) is True
    df = parse_annotation(str(p))
    assert list(df.columns[:4]) == [
        "transcript_id",
        "cds_start",
        "cds_end",
        "transcript_length",
    ]


def test_check_annotation_extended(tmp_path):
    p = tmp_path / "anno_ext.tsv"
    p.write_text(
        "\t".join(
            [
                "transcript_id",
                "cds_start",
                "cds_end",
                "transcript_length",
                "genomic_cds_starts",
                "genomic_cds_ends",
            ]
        )
        + "\nTX2\t5\t95\t120\t1,10\t20,30\n"
    )
    assert check_annotation(str(p)) is True
    df = parse_annotation(str(p))
    assert df.loc[0, "genomic_cds_starts"] == "1,10"
    assert df.loc[0, "genomic_cds_ends"] == "20,30"


def test_check_annotation_missing_required(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("transcript_id\tcds_start\nTX3\t0\n")
    assert check_annotation(str(p)) is False
    try:
        parse_annotation(str(p))
        assert False, "Expected ValueError for missing columns"
    except ValueError:
        pass

