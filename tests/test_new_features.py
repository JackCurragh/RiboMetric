"""Tests for the new sample-level QC features.

* recommend_read_lengths      — which read lengths to keep downstream
* classify_library_type       — elongation / initiation / low quality
* gene_body_coverage_ramp     — 5'->3' relative-CDS metagene
* library_complexity_curve    — analytic rarefaction / saturation
"""

import numpy as np
import pandas as pd

from RiboMetric.metrics import classify_library_type, recommend_read_lengths
from RiboMetric.modules import gene_body_coverage_ramp, library_complexity_curve


# --------------------------------------------------------------------------- #
# recommend_read_lengths
# --------------------------------------------------------------------------- #
def test_recommend_read_lengths_selects_periodic_abundant_lengths():
    # 28 nt: strongly periodic and abundant -> recommended
    # 21 nt: periodic but rare -> not recommended (below read proportion)
    # 35 nt: abundant but flat -> not recommended (below periodicity)
    read_frame = {
        28: {0: 900, 1: 50, 2: 50},
        21: {0: 9, 1: 0, 2: 1},
        35: {0: 340, 1: 330, 2: 330},
    }
    rld = {28: 1000, 21: 10, 35: 1000}
    offsets = {28: 12, 21: 11, 35: 15}
    out = recommend_read_lengths(
        read_frame, rld, offsets, min_periodicity=0.5, min_read_proportion=0.05
    )
    assert out["recommended_lengths"] == [28]
    assert out["n_recommended"] == 1
    assert out["by_read_length"][28]["offset"] == 12
    # 1000 of 2010 reads are in the recommended length
    assert abs(out["recommended_read_proportion"] - 1000 / 2010) < 1e-3


def test_recommend_read_lengths_threshold_changes_selection():
    read_frame = {28: {0: 700, 1: 150, 2: 150}}  # periodicity 0.7
    rld = {28: 1000}
    strict = recommend_read_lengths(read_frame, rld, min_periodicity=0.8)
    loose = recommend_read_lengths(read_frame, rld, min_periodicity=0.6)
    assert strict["recommended_lengths"] == []
    assert loose["recommended_lengths"] == [28]


# --------------------------------------------------------------------------- #
# classify_library_type
# --------------------------------------------------------------------------- #
def test_classify_elongation():
    out = classify_library_type(
        periodicity=0.8, prop_reads_cds=0.8, start_codon_enrichment_ratio=1.0
    )
    assert out["label"] == "elongation"


def test_classify_initiation():
    out = classify_library_type(
        periodicity=0.8, prop_reads_cds=0.8, start_codon_enrichment_ratio=8.0
    )
    assert out["label"] == "initiation"


def test_classify_low_quality():
    out = classify_library_type(
        periodicity=0.2, prop_reads_cds=0.2, start_codon_enrichment_ratio=None
    )
    assert out["label"] == "low_quality"


# --------------------------------------------------------------------------- #
# gene_body_coverage_ramp
# --------------------------------------------------------------------------- #
def test_gene_body_ramp_detects_five_prime_enrichment():
    # Pile reads near the CDS 5' end of a 300 nt CDS
    a_sites = list(range(1, 31)) + list(range(140, 160))
    df = pd.DataFrame(
        {
            "transcript_id": ["tx1"] * len(a_sites),
            "a_site": a_sites,
            "cds_start": [0] * len(a_sites),
            "cds_end": [300] * len(a_sites),
            "count": [1] * len(a_sites),
        }
    )
    out = gene_body_coverage_ramp(df, n_bins=10)
    assert len(out["profile"]) == 10
    assert out["five_prime_ramp_ratio"] is not None
    assert out["five_prime_ramp_ratio"] > 1.0  # enriched at 5'


def test_gene_body_ramp_empty():
    df = pd.DataFrame(columns=["transcript_id", "a_site", "cds_start", "cds_end", "count"])
    out = gene_body_coverage_ramp(df)
    assert out["five_prime_ramp_ratio"] is None


# --------------------------------------------------------------------------- #
# library_complexity_curve
# --------------------------------------------------------------------------- #
def test_complexity_curve_monotonic_and_saturating():
    rng = np.random.default_rng(0)
    n = 500
    df = pd.DataFrame(
        {
            "transcript_id": rng.integers(0, 50, n).astype(str),
            "a_site": rng.integers(0, 300, n),
            "count": np.ones(n, dtype=int),
            "read_name": [f"r{i}" for i in range(n)],
        }
    )
    out = library_complexity_curve(df, n_points=10)
    d = out["distinct_positions"]
    assert len(d) == 10
    # rarefaction curve is non-decreasing
    assert all(d[i] <= d[i + 1] + 1e-6 for i in range(len(d) - 1))
    assert out["total_distinct_positions"] > 0
    assert 0.0 <= out["marginal_discovery_rate"] <= 1.0


def test_complexity_curve_low_complexity_saturates():
    # Every read hits one of 3 positions -> highly saturated, low marginal rate
    df = pd.DataFrame(
        {
            "transcript_id": ["tx1"] * 300,
            "a_site": [10, 20, 30] * 100,
            "count": [1] * 300,
            "read_name": [f"r{i}" for i in range(300)],
        }
    )
    out = library_complexity_curve(df, n_points=10)
    assert out["total_distinct_positions"] == 3
    assert out["marginal_discovery_rate"] < 0.01


# --------------------------------------------------------------------------- #
# compute_codon_dwell_times
# --------------------------------------------------------------------------- #
def test_codon_dwell_times_detects_paused_codon():
    from RiboMetric.rust import compute_codon_dwell_times

    # Build a synthetic transcript: 60-codon CDS. Make codon "CGA" appear at a
    # few elongation positions and pile most reads onto those positions, so CGA
    # should show an elevated dwell-time.
    codons = []
    cga_positions = set()
    for i in range(60):
        if i in (20, 30, 40):
            codons.append("CGA")
            cga_positions.add(i)
        else:
            codons.append("GCT")  # Ala, the "background" codon
    cds = "ATG" + "".join(codons) + "TAA"  # start + body + stop
    fasta = {"tx1": cds}
    ann = pd.DataFrame({"transcript_id": ["tx1"], "cds_start": [0], "cds_end": [len(cds)]})

    rows = []
    # Heavy reads on the CGA codons, light elsewhere. cds_start=0, so the A-site
    # position of codon i (0-based, including the ATG at codon 0) is 3*(i+1).
    for i in range(1, 61):
        a_site = 3 * i
        weight = 50 if (i - 1) in cga_positions else 1
        rows.append(
            {
                "transcript_id": "tx1",
                "reference_name": "tx1",
                "a_site": a_site,
                "cds_start": 0,
                "cds_end": len(cds),
                "count": weight,
                "read_name": f"r{i}",
            }
        )
    df = pd.DataFrame(rows)

    out = compute_codon_dwell_times(df, ann, fasta, elongation_5_nt=9, elongation_3_nt=9)
    assert out["transcripts_used"] == 1
    assert "CGA" in out["dwell_times"]
    # CGA carries far more reads than its exposure warrants -> dwell > 1
    assert out["dwell_times"]["CGA"] > 1.5
    assert out["cga_dwell"] is not None


def test_codon_dwell_times_empty_without_cds_reads():
    from RiboMetric.rust import compute_codon_dwell_times

    df = pd.DataFrame(
        columns=["transcript_id", "reference_name", "a_site", "cds_start", "cds_end", "count"]
    )
    out = compute_codon_dwell_times(df, pd.DataFrame(), {})
    assert out["transcripts_used"] == 0
    assert out["dwell_times"] == {}


# --------------------------------------------------------------------------- #
# floss_library_heterogeneity
# --------------------------------------------------------------------------- #
def test_floss_flags_aberrant_transcript():
    from RiboMetric.modules import floss_library_heterogeneity

    rows = []
    # 9 "normal" transcripts: length 29; 1 aberrant transcript: length 21.
    for t in range(9):
        for _ in range(30):
            rows.append({"transcript_id": f"tx{t}", "read_length": 29, "count": 1})
    for _ in range(30):
        rows.append({"transcript_id": "tx_bad", "read_length": 21, "count": 1})
    df = pd.DataFrame(rows)
    df["read_name"] = [f"r{i}" for i in range(len(df))]
    out = floss_library_heterogeneity(df, min_reads_per_transcript=20, floss_cutoff=0.3)
    assert out["n_transcripts_scored"] == 10
    # The single 21 nt transcript departs from the 29-nt-dominated reference.
    assert out["floss_aberrant_transcript_fraction"] > 0
    assert out["floss_median"] is not None


def test_floss_empty():
    from RiboMetric.modules import floss_library_heterogeneity

    out = floss_library_heterogeneity(
        pd.DataFrame(columns=["transcript_id", "read_length", "count"])
    )
    assert out["n_transcripts_scored"] == 0
