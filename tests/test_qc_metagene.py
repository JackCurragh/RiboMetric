import pandas as pd

from RiboMetric.qc import (
    _frame_calibrated_offsets,
    _metagene_profile_stats,
    _metagene_report_reads,
    _offset_audit_record,
)


def test_metagene_report_reads_respects_unique_only_filter():
    df = pd.DataFrame(
        {
            "read_name": ["unique", "multi"],
            "mapq": [255, 0],
            "count": [1, 1],
        }
    )
    config = {"argument": {"multimap_filter": "unique_only"}}

    out = _metagene_report_reads(df, config)

    assert out["read_name"].tolist() == ["unique"]


def test_metagene_report_reads_keeps_all_when_multimap_filter_none():
    df = pd.DataFrame(
        {
            "read_name": ["unique", "multi"],
            "mapq": [255, 0],
            "count": [1, 1],
        }
    )
    config = {"argument": {"multimap_filter": "none"}}

    out = _metagene_report_reads(df, config)

    assert out["read_name"].tolist() == ["unique", "multi"]


def test_metagene_profile_stats_records_filter_and_window_support():
    annotated = pd.DataFrame({"count": [2, 3, 5]})
    profile_reads = pd.DataFrame({"count": [2, 3]})
    profile = {
        "start": {28: {-1: 2, 0: 3}},
        "stop": {28: {0: 1}},
    }
    config = {"argument": {"multimap_filter": "unique_only"}}

    stats = _metagene_profile_stats(annotated, profile_reads, profile, config)

    assert stats == {
        "multimap_filter": "unique_only",
        "input_reads": 10,
        "profile_reads": 5,
        "start_window_reads": 5,
        "stop_window_reads": 1,
    }


def test_offset_audit_uses_final_mapping_and_retains_raw_offset():
    # Simulates the pre-calibration dataframe still carrying the raw estimate.
    reads = pd.DataFrame({"read_length": [28] * 60, "offset": [14] * 60})
    config = {"argument": {"offset_calculation_method": "ribowaltz"}}

    audit = _offset_audit_record(
        reads,
        config,
        {28: 15},
        {"28": {"old_offset": 14, "new_offset": 15}},
        (8, 20),
        2 / 3,
        "a_site",
        15,
        raw_offsets={28: 14},
    )

    assert audit["raw_offsets"] == {"28": 14}
    assert audit["final_offsets"] == {"28": 15}
    assert audit["computed_offsets"] == {"28": 15}
    assert audit["applied_by_read_length"]["28"]["offsets"] == [15]


def test_external_offsets_are_explicitly_marked_final():
    reads = pd.DataFrame({"read_length": [28], "offset": [14]})
    audit = _offset_audit_record(
        reads,
        {"argument": {"offset_read_length": "offsets.tsv"}},
        {},
        {},
        (8, 20),
        2 / 3,
        "a_site",
        15,
    )

    assert audit["offsets_are_final"] is True
    assert audit["offset_calibration"] == "not_applied_external_offsets_treated_as_final"


def test_frame_calibrated_offsets_nudges_confident_nonzero_frame():
    df = pd.DataFrame(
        {
            "read_length": [28] * 60,
            "a_site": [102 + (i * 3) for i in range(60)],
            "cds_start": [100] * 60,
            "cds_end": [400] * 60,
            "mapq": [255] * 60,
            "count": [1] * 60,
        }
    )
    config = {
        "argument": {"multimap_filter": "unique_only", "global_offset": 15},
        "qc": {
            "read_frame_distribution": {
                "exclude_codons": 0,
                "offset_frame_correction_min_reads": 50,
                "offset_frame_correction_min_fraction": 0.6,
            }
        },
    }

    offsets, adjustments = _frame_calibrated_offsets(
        df,
        {28: 13},
        config,
        offset_bounds=(8, 20),
        max_read_length_fraction=2 / 3,
    )

    assert offsets[28] == 14
    assert adjustments["28"]["dominant_frame"] == 2
    assert adjustments["28"]["old_offset"] == 13
    assert adjustments["28"]["new_offset"] == 14


def test_frame_calibrated_offsets_ignores_low_count_lengths():
    df = pd.DataFrame(
        {
            "read_length": [28] * 10,
            "a_site": [102 + (i * 3) for i in range(10)],
            "cds_start": [100] * 10,
            "cds_end": [400] * 10,
            "mapq": [255] * 10,
            "count": [1] * 10,
        }
    )
    config = {
        "argument": {"multimap_filter": "unique_only", "global_offset": 15},
        "qc": {
            "read_frame_distribution": {
                "exclude_codons": 0,
                "offset_frame_correction_min_reads": 50,
                "offset_frame_correction_min_fraction": 0.6,
            }
        },
    }

    offsets, adjustments = _frame_calibrated_offsets(
        df,
        {28: 13},
        config,
        offset_bounds=(8, 20),
        max_read_length_fraction=2 / 3,
    )

    assert offsets[28] == 13
    assert adjustments == {}


def test_frame_calibrated_offsets_skips_when_nudge_is_clamped():
    """A nudge that falls out of bounds must not snap to the default offset.

    With ``dominant_frame == 1`` the shift is ``-1``. If the current offset sits
    on the lower bound, ``old - 1`` is invalid and ``sanitise_offset`` falls back
    to the (valid) default. That default is not a single-frame step, so applying
    it could land the dominant frame on frame 2 instead of frame 0. The offset
    must be left untouched.
    """
    df = pd.DataFrame(
        {
            # a_site - cds_start == 1 -> dominant frame 1 -> shift = -1
            "read_length": [28] * 60,
            "a_site": [101 + (i * 3) for i in range(60)],
            "cds_start": [100] * 60,
            "cds_end": [400] * 60,
            "mapq": [255] * 60,
            "count": [1] * 60,
        }
    )
    config = {
        "argument": {"multimap_filter": "unique_only", "global_offset": 15},
        "qc": {
            "read_frame_distribution": {
                "exclude_codons": 0,
                "offset_frame_correction_min_reads": 50,
                "offset_frame_correction_min_fraction": 0.6,
            }
        },
    }

    # old_offset == lower bound (8); old - 1 == 7 is out of bounds, so
    # sanitise_offset would otherwise return the default (15).
    offsets, adjustments = _frame_calibrated_offsets(
        df,
        {28: 8},
        config,
        offset_bounds=(8, 20),
        max_read_length_fraction=2 / 3,
    )

    assert offsets[28] == 8
    assert adjustments == {}


def test_frame_calibrated_offsets_nudges_frame_one_down():
    """In-bounds frame-1 dominance shifts the offset down by exactly one."""
    df = pd.DataFrame(
        {
            "read_length": [28] * 60,
            "a_site": [101 + (i * 3) for i in range(60)],
            "cds_start": [100] * 60,
            "cds_end": [400] * 60,
            "mapq": [255] * 60,
            "count": [1] * 60,
        }
    )
    config = {
        "argument": {"multimap_filter": "unique_only", "global_offset": 15},
        "qc": {
            "read_frame_distribution": {
                "exclude_codons": 0,
                "offset_frame_correction_min_reads": 50,
                "offset_frame_correction_min_fraction": 0.6,
            }
        },
    }

    offsets, adjustments = _frame_calibrated_offsets(
        df,
        {28: 13},
        config,
        offset_bounds=(8, 20),
        max_read_length_fraction=2 / 3,
    )

    assert offsets[28] == 12
    assert adjustments["28"]["dominant_frame"] == 1
    assert adjustments["28"]["new_offset"] == 12


def test_frame_calibrated_offsets_ignores_transcript_frame_ambiguous_fragments():
    df = pd.DataFrame(
        {
            "read_name": ["safe"] * 60 + ["ambig", "ambig"],
            "gene_id": ["g1"] * 60 + ["g2", "g2"],
            "transcript_id": ["tx1"] * 60 + ["tx2", "tx3"],
            "read_length": [28] * 62,
            "a_site": [100 + (i * 3) for i in range(60)] + [101, 102],
            "cds_start": [100] * 62,
            "cds_end": [400] * 62,
            "mapq": [255] * 62,
            "count": [1] * 62,
        }
    )
    config = {
        "argument": {"multimap_filter": "unique_only", "global_offset": 15},
        "qc": {
            "read_frame_distribution": {
                "exclude_codons": 0,
                "offset_frame_correction_min_reads": 50,
                "offset_frame_correction_min_fraction": 0.6,
            }
        },
    }

    offsets, adjustments = _frame_calibrated_offsets(
        df,
        {28: 13},
        config,
        offset_bounds=(8, 20),
        max_read_length_fraction=2 / 3,
    )

    assert offsets[28] == 13
    assert adjustments == {}
