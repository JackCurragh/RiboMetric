"""Hand-computed answers for the coverage, region and library-type definitions
in docs/METRIC_CONTRACT.md (sections 3.2 and 5)."""

import pandas as pd
import pytest

from RiboMetric.metrics import (
    cds_coverage_metric,
    classify_library_type,
    proportion_of_reads_in_region,
    region_region_ratio_metric,
)
from RiboMetric.modules import gene_body_coverage_ramp

# CDS [10, 20): body positions 11..19 (9 nt); frame-0 body positions 13, 16, 19.
CDS_READS = pd.DataFrame(
    {
        "transcript_id": ["tx1"] * 4 + ["tx2"],
        "a_site": [11, 13, 16, 16, 55],
        "cds_start": [10, 10, 10, 10, 50],
        "cds_end": [20, 20, 20, 20, 60],
        "count": [1, 2, 1, 1, 1],
    }
)


@pytest.mark.parametrize(
    "min_reads,in_frame,expected",
    [
        (1, False, 3 / 9),  # positions 11, 13, 16 of 9
        (1, True, 2 / 3),  # frame-0 positions 13, 16 of 3
        (2, False, 2 / 9),  # positions with weight >= 2: 13 (2), 16 (1 + 1)
    ],
)
def test_cds_coverage_top_transcript(min_reads, in_frame, expected):
    value = cds_coverage_metric(
        CDS_READS, minimum_reads=min_reads, in_frame_coverage=in_frame, num_transcripts=1
    )
    assert value == pytest.approx(expected)


MRNA = {
    28: {"five_leader": 10, "start_codon": 0, "CDS": 80, "stop_codon": 0, "three_trailer": 10},
    45: {"five_leader": 0, "start_codon": 0, "CDS": 10, "stop_codon": 0, "three_trailer": 0},
}
MRNA["global"] = {
    region: MRNA[28][region] + MRNA[45][region] for region in MRNA[28]  # type: ignore[index]
}


def test_proportion_of_reads_in_region_per_length_and_pooled():
    prop = proportion_of_reads_in_region(MRNA, region="CDS")
    assert prop[28] == pytest.approx(0.8)
    assert prop[45] == pytest.approx(1.0)
    assert prop["global"] == pytest.approx(90 / 110)


def test_region_ratio_sums_all_observed_read_lengths():
    # Regional ratios use every observed length unless a caller explicitly
    # supplies a range.
    ratio = region_region_ratio_metric(MRNA, region1="CDS", region2="five_leader")
    assert ratio[28] == pytest.approx(8.0)
    assert ratio[45] is None
    assert ratio["global"] == pytest.approx(9.0)


def test_gene_body_ramp_and_drop():
    reads = pd.DataFrame(
        {"a_site": [5, 50], "cds_start": [0, 0], "cds_end": [100, 100], "count": [4, 1]}
    )
    ramp = gene_body_coverage_ramp(reads)
    # 5 reads over 100 bins, mean 0.05: bin 5 -> 80, bin 50 -> 20 after
    # normalisation. First tenth averages 8, the 40-59 middle averages 1.
    assert ramp["five_prime_ramp_ratio"] == pytest.approx(8.0)
    assert ramp["three_prime_drop_ratio"] == pytest.approx(0.0)


@pytest.mark.parametrize(
    "periodicity,prop_cds,start_ratio,label",
    [
        (0.30, 0.90, 1.0, "low_quality"),  # weak periodicity
        (0.80, 0.40, 1.0, "low_quality"),  # little CDS signal
        (0.80, 0.70, 4.0, "initiation"),  # start peak >= 3x the body
        (0.80, 0.70, 1.0, "elongation"),
        (0.80, 0.70, None, "elongation"),
    ],
)
def test_classify_library_type(periodicity, prop_cds, start_ratio, label):
    result = classify_library_type(periodicity, prop_cds, start_ratio)
    assert result["label"] == label
    assert result["evidence"]["periodicity"] == pytest.approx(periodicity)
