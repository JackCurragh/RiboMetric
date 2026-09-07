import pytest

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


def test_parse_annotation_rejects_genomic_coordinates(tmp_path):
    """CDS coordinates in genomic/unspliced space (cds_end overruns the spliced
    transcript) must be rejected, since they silently scramble reading frame."""
    p = tmp_path / "genomic.tsv"
    rows = ["transcript_id\tcds_start\tcds_end\ttranscript_length"]
    # cds_end == cds_start + transcript_length -> overruns for any cds_start > 0
    for i in range(10):
        rows.append(f"TX{i}\t{100 + i}\t{3000 + i}\t{2900}")
    p.write_text("\n".join(rows) + "\n")
    with pytest.raises(ValueError, match="overruns|genomic|transcript-relative"):
        parse_annotation(str(p))


def test_parse_annotation_accepts_fully_coding_transcript(tmp_path):
    """A legitimate transcript that is entirely CDS (cds_start=0,
    cds_end=transcript_length) must NOT be flagged."""
    p = tmp_path / "coding.tsv"
    p.write_text(
        "transcript_id\tcds_start\tcds_end\ttranscript_length\n"
        "TX1\t0\t300\t300\n"
        "TX2\t30\t297\t400\n"
    )
    df = parse_annotation(str(p))
    assert len(df) == 2


def test_parse_annotation_tolerates_rare_overrun(tmp_path):
    """A single malformed row in an otherwise-valid file should not block the
    whole annotation (only systematic violations are rejected)."""
    rows = ["transcript_id\tcds_start\tcds_end\ttranscript_length"]
    for i in range(200):
        rows.append(f"TX{i}\t30\t297\t400")
    rows.append("BAD\t30\t9999\t400")  # one overrunning row (<1%)
    p = tmp_path / "mostly_ok.tsv"
    p.write_text("\n".join(rows) + "\n")
    df = parse_annotation(str(p))
    assert len(df) == 201

