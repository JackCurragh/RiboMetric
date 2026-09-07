"""
This script contains tests for the different functions found in modules.py
"""

import pandas as pd
import pytest

from RiboMetric.modules import (
    a_site_calculation,
    annotate_reads,
    asite_calculation_per_readlength,
    assign_mRNA_category,
    filter_unique_mappers,
    filter_unique_mappers_strict,
    metagene_profile,
    mRNA_distribution,
    normalise_ligation_bias,
    nucleotide_composition,
    read_frame_distribution,
    read_frame_distribution_annotated,
    read_length_distribution,
    representative_transcripts_for_offset_metagene,
    terminal_nucleotide_bias_distribution,
    unique_fragments_single_gene_for_offset_metagene,
)


def test_a_site_calculation(sample_read_df):
    """
    Test A-site calculation
    """
    a_site_df = a_site_calculation(sample_read_df)
    assert a_site_df.a_site[0] == 356


def test_read_length_distribution(sample_read_df):
    """
    Test read length distribution calculation
    """
    read_length_dict = read_length_distribution(sample_read_df)
    assert read_length_dict[29] == 40


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ("AA_test", 0.5),
        ('terminal_nucleotide_bias_dict_norm["three_prime"]["TC"]', 0.275),
    ],
)
def test_terminal_nucleotide_bias_distribution(test_input, expected, sample_read_df):
    """
    Test ligation bias distribution calculation
    """
    read_df = sample_read_df.copy()
    categories = ["first_dinucleotide", "last_dinucleotide"]
    read_df[categories] = read_df[categories].astype("category")
    sequence_background = {
        "5_prime_bg": {"AA": 0.25, "AG": 0.3, "TT": 0.45},
        "3_prime_bg": {"AA": 0.875, "TC": 0.125},
    }
    terminal_nucleotide_bias_dict = terminal_nucleotide_bias_distribution(read_df)
    AA_test = terminal_nucleotide_bias_dict["five_prime"]["AA"]
    terminal_nucleotide_bias_dict_norm = normalise_ligation_bias(
        terminal_nucleotide_bias_dict, sequence_background
    )

    # Just to satisfy the linter
    type(AA_test)
    type(terminal_nucleotide_bias_dict_norm)
    assert eval(test_input) == expected


def test_nucleotide_composition():
    """
    Test nucleotide composition calculation
    """
    sequence_data = {
        1: {
            "A": [4, 4, 5, 5, 1, 0, 0, 3, 3, 3, 2],
            "C": [1, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0],
            "G": [3, 3, 0, 0, 4, 5, 5, 0, 0, 0, 0],
            "T": [0, 1, 3, 0, 0, 0, 0, 5, 5, 5, 0],
        }
    }
    nucleotide_composition_dict = nucleotide_composition(sequence_data[1])
    assert nucleotide_composition_dict["A"] == [
        0.5,
        0.5,
        0.625,
        0.625,
        0.125,
        0,
        0,
        0.375,
        0.375,
        0.375,
        1,
    ]


def test_read_frame_distribution(sample_read_df):
    """
    Test read frame labelling
    """
    read_frame_dict = read_frame_distribution(a_site_calculation(sample_read_df))
    assert read_frame_dict[33][1] == 10


def test_read_frame_distribution_annotated_drops_transcript_ambiguous_fragments():
    df = pd.DataFrame(
        {
            "read_name": ["safe", "safe", "ambig", "ambig", "cross", "cross"],
            "read_length": [28, 28, 28, 28, 28, 28],
            "gene_id": ["g1", "g1", "g2", "g2", "g3", "g4"],
            "transcript_id": ["tx1", "tx2", "tx3", "tx4", "tx5", "tx6"],
            "a_site": [103, 103, 103, 104, 103, 103],
            "cds_start": [100, 100, 100, 100, 100, 100],
            "cds_end": [400, 400, 400, 400, 400, 400],
            "mapq": [255, 255, 255, 255, 255, 255],
            "mapq_available": [True, True, True, True, True, True],
            "count": [1, 1, 1, 1, 1, 1],
        }
    )

    out = read_frame_distribution_annotated(df, exclusion_length=0, unique_only=True)

    assert out == {}


def test_read_frame_distribution_annotated_keeps_single_transcript_weighted_fragment_once():
    df = pd.DataFrame(
        {
            "read_name": ["safe"],
            "read_length": [29],
            "gene_id": ["g1"],
            "transcript_id": ["tx1"],
            "a_site": [103],
            "cds_start": [100],
            "cds_end": [400],
            "mapq": [255],
            "mapq_available": [True],
            "count": [7],
        }
    )

    out = read_frame_distribution_annotated(df, exclusion_length=0, unique_only=True)

    assert out[29] == {0: 7, 1: 0, 2: 0}


def test_read_frame_distribution_annotated_handles_categorical_count():
    """Regression: ``count`` arrives as a non-ordered Categorical from BAM
    parsing. ``frame_safe_unique_fragments`` must not crash aggregating it
    (pandas refuses ``max`` over an unordered Categorical).
    """
    df = pd.DataFrame(
        {
            "read_name": ["safe", "safe"],
            "read_length": [29, 29],
            "gene_id": ["g1", "g1"],
            "transcript_id": ["tx1", "tx1"],
            "a_site": [103, 103],
            "cds_start": [100, 100],
            "cds_end": [400, 400],
            "mapq": [255, 255],
            "mapq_available": [True, True],
            "count": pd.Series([3, 7], dtype="category"),
        }
    )

    out = read_frame_distribution_annotated(df, exclusion_length=0, unique_only=True)

    # The duplicate read_name collapses to one fragment carrying max(count)=7.
    assert out[29] == {0: 7, 1: 0, 2: 0}


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ('mRNA_distribution_dict[29]["CDS"]', 30),
        ('mRNA_distribution_dict[21]["three_trailer"]', 40),
    ],
)
def test_mRNA_distribution(test_input, expected, sample_read_df, sample_annotation_df):
    """
    Test metagene distance calculations
    """
    a_site_df = a_site_calculation(sample_read_df)
    annotated_read_df = annotate_reads(a_site_df, sample_annotation_df)
    annotated_read_df = assign_mRNA_category(annotated_read_df)

    mRNA_distribution_dict = mRNA_distribution(annotated_read_df)

    # Just to satisfy the linter
    type(mRNA_distribution_dict)
    assert eval(test_input) == expected


def test_metagene_profile(sample_read_df, sample_annotation_df):
    """
    Test metagene distance calculations
    """
    a_site_df = a_site_calculation(sample_read_df)
    annotated_read_df = annotate_reads(a_site_df, sample_annotation_df)
    metagene_profile_dict = metagene_profile(annotated_read_df)

    assert metagene_profile_dict["stop"][21][4] == 30


def test_annotate_reads_extracts_gene_metadata_from_reference_name(sample_annotation_df):
    df = pd.DataFrame(
        {
            "read_name": ["r1"],
            "read_length": [29],
            "reference_name": ["tx1|gene1|havana|ott|txname1|genename1|100|protein_coding|"],
            "reference_start": [10],
            "a_site": [25],
            "count": [1],
        }
    )
    ann = pd.DataFrame(
        {
            "transcript_id": ["tx1"],
            "cds_start": [20],
            "cds_end": [80],
            "transcript_length": [100],
        }
    )

    out = annotate_reads(df, ann)

    assert out["gene_id"].iloc[0] == "gene1"
    assert out["gene_name"].iloc[0] == "genename1"


def test_representative_transcripts_for_offset_metagene_keeps_one_per_gene_and_length():
    df = pd.DataFrame(
        {
            "gene_id": ["g1", "g1", "g1", "g2"],
            "transcript_id": ["tx1", "tx1", "tx2", "tx3"],
            "read_length": [29, 29, 29, 29],
            "count": [1, 1, 5, 2],
        }
    )

    out = representative_transcripts_for_offset_metagene(df)

    assert sorted(out["transcript_id"].tolist()) == ["tx2", "tx3"]


def test_unique_fragments_single_gene_for_offset_metagene_drops_cross_gene_reads():
    df = pd.DataFrame(
        {
            "read_name": ["r1", "r1", "r2", "r2", "r3"],
            "gene_id": ["g1", "g1", "g2", "g3", "g4"],
            "transcript_id": ["tx1", "tx2", "tx3", "tx4", "tx5"],
        }
    )

    out = unique_fragments_single_gene_for_offset_metagene(df)

    assert sorted(out["read_name"].unique().tolist()) == ["r1", "r3"]


def test_filter_unique_mappers_requires_available_mapq():
    df = pd.DataFrame(
        {
            "mapq": [255, 0, 255],
            "mapq_available": [True, True, False],
            "count": [1, 1, 1],
        }
    )

    unique_df = filter_unique_mappers(df, enabled=True)

    assert unique_df.index.tolist() == [0]


def test_filter_unique_mappers_does_not_fallback_to_all_rows():
    df = pd.DataFrame(
        {
            "mapq": [0, 0],
            "mapq_available": [True, False],
            "count": [1, 1],
        }
    )

    unique_df = filter_unique_mappers(df, enabled=True)

    assert unique_df.empty


def test_filter_unique_mappers_prefers_mapq_for_frame_sensitive_filtering():
    df = pd.DataFrame(
        {
            "mapq": [255, 255, 255],
            "mapq_available": [True, True, True],
            "xa": [0, 1, 2],
            "count": [1, 1, 1],
        }
    )

    unique_df = filter_unique_mappers(df, enabled=True)

    assert unique_df.index.tolist() == [0, 1, 2]


def test_filter_unique_mappers_strict_removes_repeated_read_names():
    df = pd.DataFrame(
        {
            "read_name": ["read1", "read1", "read2"],
            "reference_name": ["tx1", "tx2", "tx1"],
            "reference_start": [10, 10, 20],
            "mapq": [255, 255, 255],
            "mapq_available": [True, True, True],
            "count": [1, 1, 1],
        }
    )

    unique_df = filter_unique_mappers_strict(df, enabled=True)

    assert unique_df["read_name"].tolist() == ["read2"]


def test_changepoint_offset_search_includes_upper_bound():
    rows = []
    for i in range(10):
        rows.append(
            {
                "read_name": f"peak18_{i}",
                "read_length": 30,
                "reference_start": 82,
                "cds_start": 100,
                "cds_end": 400,
                "mapq": 255,
                "mapq_available": True,
                "count": 1,
            }
        )
    for i in range(3):
        rows.append(
            {
                "read_name": f"peak17_{i}",
                "read_length": 30,
                "reference_start": 83,
                "cds_start": 100,
                "cds_end": 400,
                "mapq": 255,
                "mapq_available": True,
                "count": 1,
            }
        )
    df = pd.DataFrame(rows)

    offsets = asite_calculation_per_readlength(
        df,
        method="changepoint",
        offset_range=(10, 18),
        offset_target="p_site",
    )

    assert offsets[30] == 18


def test_offset_calculation_keeps_same_gene_repeated_read_names():
    df = pd.DataFrame(
        {
            "read_name": ["read1", "read1", "read2", "read2"],
            "read_length": [30, 30, 30, 30],
            "gene_id": ["g1", "g1", "g1", "g1"],
            "transcript_id": ["tx1", "tx2", "tx1", "tx2"],
            "reference_start": [82, 82, 82, 83],
            "cds_start": [100, 100, 100, 100],
            "cds_end": [400, 400, 400, 400],
            "mapq": [255, 255, 255, 255],
            "mapq_available": [True, True, True, True],
            "count": [1, 1, 1, 1],
        }
    )

    offsets = asite_calculation_per_readlength(
        df,
        method="changepoint",
        offset_range=(10, 18),
        offset_target="p_site",
    )

    assert offsets[30] == 18


def test_offset_calculation_drops_repeated_read_names_with_conflicting_peaks():
    rows = []
    for i in range(10):
        rows.append(
            {
                "read_name": f"unique18_{i}",
                "read_length": 30,
                "gene_id": "g1",
                "transcript_id": "tx1",
                "reference_start": 82,
                "cds_start": 100,
                "cds_end": 400,
                "mapq": 255,
                "mapq_available": True,
                "count": 1,
            }
        )
    for i in range(20):
        rows.extend(
            [
                {
                    "read_name": f"dup17_{i}",
                    "read_length": 30,
                    "gene_id": "g2",
                    "transcript_id": "tx2",
                    "reference_start": 83,
                    "cds_start": 100,
                    "cds_end": 400,
                    "mapq": 255,
                    "mapq_available": True,
                    "count": 1,
                },
                {
                    "read_name": f"dup17_{i}",
                    "read_length": 30,
                    "gene_id": "g3",
                    "transcript_id": "tx3",
                    "reference_start": 83,
                    "cds_start": 100,
                    "cds_end": 400,
                    "mapq": 255,
                    "mapq_available": True,
                    "count": 1,
                },
            ]
        )
    df = pd.DataFrame(rows)

    offsets = asite_calculation_per_readlength(
        df,
        method="changepoint",
        offset_range=(10, 18),
        offset_target="p_site",
    )

    assert offsets[30] == 18
