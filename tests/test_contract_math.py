"""Hand-computed answers for every metric definition in docs/METRIC_CONTRACT.md.

Each test states its input small enough to compute the expected value by hand
and asserts that value, not merely a range.
"""

import math

import pandas as pd
import pytest

from RiboMetric.metrics import (
    cds_enrichment_ratio,
    fourier_transform,
    information_metric_cutoff,
    periodicity_autocorrelation,
    periodicity_dominance,
    read_frame_information_content,
    read_frame_information_weighted_score,
    read_length_bimodality_coefficient,
    read_length_cv,
    read_length_iqr_fraction,
    read_length_max_proportion,
    read_length_normality_pvalue,
    recommend_read_lengths,
    terminal_nucleotide_bias_KL_divergence,
    terminal_nucleotide_bias_max_deviation,
    uniformity_autocorrelation,
    uniformity_entropy,
    uniformity_gini_index,
    uniformity_theil_index,
)
from RiboMetric.modules import (
    _get_weights,
    floss_library_heterogeneity,
    library_complexity_curve,
    metagene_profile,
    read_frame_score_trips_viz,
    reading_frame_triangle,
)


def _start_profile(series_by_length):
    return {"start": series_by_length, "stop": {}}


# --- 3.1 frame and periodicity ---------------------------------------------


def test_periodicity_dominance_per_length_global_and_min_reads():
    frames = {
        28: {0: 80, 1: 10, 2: 10},
        29: {0: 5, 1: 90, 2: 5},
        30: {0: 10, 1: 0, 2: 0},  # 10 reads: below the 100-read minimum
    }
    d = periodicity_dominance(frames, min_reads=100)
    assert d[28] == pytest.approx(0.8)
    assert d[29] == pytest.approx(0.9)
    assert 30 not in d
    # global shares one frame across lengths: frame 1 holds 100 of 210 reads
    assert d["global"] == pytest.approx(100 / 210)
    assert d["global_by_read_length_max"] == pytest.approx((80 + 90 + 10) / 210)


@pytest.mark.parametrize(
    "frames,expected",
    [
        ({0: 1, 1: 0, 2: 0}, 1.0),  # single frame
        ({0: 1, 1: 1, 2: 1}, 0.0),  # uniform
        # d = 0.8, 0.1, 0.1: H = 0.921928 bits, 1 - H/log2(3) = 0.418329
        ({0: 80, 1: 10, 2: 10}, 0.418329),
    ],
)
def test_periodicity_information_single_length(frames, expected):
    info = read_frame_information_content({28: frames})
    assert info[28][0] == pytest.approx(expected, abs=1e-5)


def test_periodicity_information_cutoff_and_weighted_score():
    info = read_frame_information_content(
        {28: {0: 80, 1: 10, 2: 10}, 29: {0: 1, 1: 1, 2: 1}}  # 29: 3 of 103 reads
    )
    cut = information_metric_cutoff(info, min_count_threshold=0.05)
    assert 29 not in cut  # 3 <= 0.05 * 103
    assert cut["global"] == pytest.approx(0.418329, abs=1e-5)
    weighted = read_frame_information_weighted_score(info)
    assert weighted == pytest.approx((0.418329 * 100 + 0.0 * 3) / 103, abs=1e-5)


def test_recommend_read_lengths_all_three_conditions():
    frames = {
        28: {0: 900, 1: 50, 2: 50},  # dominance 0.9, recommended
        29: {0: 400, 1: 300, 2: 300},  # dominance 0.4 < 0.5
        30: {0: 50, 1: 0, 2: 0},  # 50 frame reads < 100: skipped
    }
    rld = {28: 1000, 29: 1000, 30: 50, 31: 3950}
    rec = recommend_read_lengths(frames, rld, min_periodicity=0.5, min_read_proportion=0.05)
    assert rec["recommended_lengths"] == [28]
    assert rec["n_recommended"] == 1
    assert rec["recommended_read_proportion"] == pytest.approx(round(1000 / 6000, 4))
    assert 30 not in rec["by_read_length"]


def test_trips_viz_excludes_single_frame_lengths_from_global():
    scores = read_frame_score_trips_viz({28: {0: 80, 1: 10, 2: 10}, 29: {0: 50, 1: 0, 2: 0}})
    assert scores[28] == pytest.approx(0.875)
    assert scores[29] == 1
    assert scores["global"] == pytest.approx(0.875)


def test_fourier_periodic_flat_and_other_period():
    periodic = {p: v for p, v in enumerate([100, 5, 5] * 20)}
    flat = {p: 7 for p in range(60)}
    period5 = {p: v for p, v in enumerate([0, 0, 0, 0, 100] * 12)}
    f = fourier_transform(_start_profile({28: periodic, 29: flat, 30: period5}))
    assert f[28] == pytest.approx(1.0)
    assert f[29] == 0.0
    assert f[30] == pytest.approx(0.0, abs=1e-12)


def test_fourier_reads_positions_in_order():
    ordered = {p: v for p, v in enumerate([100, 5, 5] * 20)}
    shuffled = dict(sorted(ordered.items(), key=lambda kv: (kv[0] % 7, kv[0])))
    f = fourier_transform(_start_profile({28: shuffled}))
    assert f[28] == pytest.approx(1.0)


def test_periodicity_autocorrelation_global_matches_per_length_statistic():
    series = {p: v for p, v in enumerate([1, 0, 0] * 20)}
    a = periodicity_autocorrelation(_start_profile({28: series}), lag=3)
    # r0 = 20, r3 = 19 for twenty [1, 0, 0] repeats
    assert a[28] == pytest.approx(19 / 20)
    assert a["global"] == pytest.approx(a[28])


# --- 3.2 coverage -----------------------------------------------------------


def test_uniformity_entropy_uses_position_order():
    # 28 nt reads at 30, 31, 33, 36, 39, 40, 42 in a 30..42 window, zero-filled
    # out of order. Codon bins from 30: [2, 1, 1, 2, 1] of 7 reads over K = 5.
    series = {
        30: 1,
        31: 1,
        33: 1,
        36: 1,
        39: 1,
        40: 1,
        42: 1,
        32: 0,
        34: 0,
        35: 0,
        37: 0,
        38: 0,
        41: 0,
    }
    p = [2 / 7, 1 / 7, 1 / 7, 2 / 7, 1 / 7]
    expected = -sum(x * math.log2(x) for x in p) / math.log2(5)
    u = uniformity_entropy(_start_profile({28: series}))
    assert u[28] == pytest.approx(expected)  # 0.962961
    assert u["global"] == pytest.approx(expected)


def test_uniformity_entropy_even_and_concentrated():
    even = {p: 1 for p in range(30)}
    one_codon = {p: (5 if p < 3 else 0) for p in range(30)}
    u = uniformity_entropy(_start_profile({28: even, 29: one_codon}))
    assert u[28] == pytest.approx(1.0)
    assert u[29] == 0


def test_metagene_profile_returns_position_ordered_profiles():
    rows = [
        dict(read_length=28, a_site=100 + p, cds_start=100, cds_end=400, count=1)
        for p in (30, 31, 33, 36, 39, 40, 42)
    ] + [
        dict(read_length=29, a_site=100 + p, cds_start=100, cds_end=400, count=1)
        for p in range(30, 43)
    ]
    profile = metagene_profile(pd.DataFrame(rows), target="start", distance_range=[30, 42])
    for series in profile["start"].values():
        assert list(series) == sorted(series)
    assert list(profile["start"]) == [28, 29]


def test_gini_and_theil_even_and_concentrated():
    even = {p: 1 for p in range(12)}  # codon bins [3, 3, 3, 3]
    single = {p: (4 if p == 0 else 0) for p in range(12)}  # bins [4, 0, 0, 0]
    gini = uniformity_gini_index(_start_profile({28: even, 29: single}))
    theil = uniformity_theil_index(_start_profile({28: even, 29: single}))
    assert gini[28] == pytest.approx(0.0)
    assert gini[29] == pytest.approx(3 / 4)  # (K - 1) / K
    assert theil[28] == pytest.approx(0.0)
    assert theil[29] == pytest.approx(math.log(4))  # ln K


def test_uniformity_autocorrelation_excludes_lag_zero():
    series = {p: 1 for p in range(30)}  # ten equal codon bins
    u = uniformity_autocorrelation(_start_profile({28: series}))
    # For n equal values r_k = (n - k) / n: lags 1-4 give 0.9, 0.8, 0.7, 0.6
    assert u[28] == pytest.approx((0.9 + 0.8 + 0.7 + 0.6) / 4)


def _uniform_annotated(tx_specs, weight=1):
    rows = []
    for tx_id, tx_len, cds_s, cds_e in tx_specs:
        for pos in range(tx_len):
            if pos < cds_s:
                cat = "five_leader"
            elif pos == cds_s:
                cat = "start_codon"
            elif pos < cds_e:
                cat = "CDS"
            elif pos == cds_e:
                cat = "stop_codon"
            else:
                cat = "three_trailer"
            rows.append(
                dict(
                    transcript_id=tx_id,
                    a_site=pos,
                    cds_start=cds_s,
                    cds_end=cds_e,
                    transcript_length=tx_len,
                    mRNA_category=cat,
                    count=weight,
                )
            )
    return pd.DataFrame(rows)


def test_cds_enrichment_uniform_coverage_is_one():
    df = _uniform_annotated([("tx1", 100, 20, 80), ("tx2", 60, 10, 50)])
    assert cds_enrichment_ratio(df) == pytest.approx(1.0)


def test_cds_enrichment_weighted_equals_expanded():
    df = _uniform_annotated([("tx1", 100, 20, 80)])
    df.loc[df["mRNA_category"] == "CDS", "count"] = 3  # CDS reads collapsed x3
    expanded = df.loc[df.index.repeat(df["count"])].drop(columns=["count"])
    assert cds_enrichment_ratio(df) == pytest.approx(cds_enrichment_ratio(expanded))
    # observed = 3*59 / (3*59 + 41) against expected 59/100
    assert cds_enrichment_ratio(df) == pytest.approx((177 / 218) / (59 / 100))


def test_cds_enrichment_ignores_reads_on_ineligible_transcripts():
    base = _uniform_annotated([("tx1", 100, 20, 80)])
    noncoding = _uniform_annotated([("nc1", 100, 0, 0)])  # no CDS body
    assert cds_enrichment_ratio(pd.concat([base, noncoding])) == pytest.approx(
        cds_enrichment_ratio(base)
    )


def test_reading_frame_triangle_uses_cds_frame_and_body():
    df = pd.DataFrame(
        {
            "transcript_id": ["tx", "tx", "tx"],
            "a_site": [13, 14, 5],  # frame 0, frame 1, leader
            "cds_start": [10, 10, 10],
            "cds_end": [100, 100, 100],
            "count": [2, 1, 7],
        }
    )
    assert reading_frame_triangle(df) == {"tx": [2, 1, 0]}


# --- 3.3 read length ----------------------------------------------------------

RLD = {27: 10, 28: 40, 29: 40, 30: 10}


def test_read_length_iqr_fraction():
    # cumulative 0.1, 0.5, 0.9, 1.0: Q.25=28, Q.75=29, Q.1=27, Q.9=29
    assert read_length_iqr_fraction(RLD) == pytest.approx(0.5)


def test_read_length_cv():
    # mean 28.5, variance 0.65
    assert read_length_cv(RLD) == pytest.approx(math.sqrt(0.65) / 28.5)


def test_read_length_max_proportion():
    assert read_length_max_proportion(RLD) == pytest.approx(0.4)
    assert read_length_max_proportion({}) is None


def test_read_length_bimodality_coefficient():
    # symmetric: skew 0; biased excess kurtosis 1.0625/0.4225 - 3 = -0.485207;
    # n = 100 correction 3 * 99^2 / (98 * 97) = 3.093099
    expected = 1 / (-0.485207 + 3.093099)
    assert read_length_bimodality_coefficient(RLD) == pytest.approx(expected, abs=1e-5)


def test_read_length_degenerate_inputs_are_none():
    assert read_length_bimodality_coefficient({28: 3}) is None  # n = 3
    assert read_length_bimodality_coefficient({28: 50}) is None  # one length
    assert read_length_normality_pvalue({28: 3, 29: 4}) is None  # n = 7


# --- 3.5 terminal bias ---------------------------------------------------------

DINUCS = [a + b for a in "ACGT" for b in "ACGT"]
UNIFORM_BG = {d: 1 / 16 for d in DINUCS}


def test_terminal_bias_zero_when_observed_equals_background():
    observed = {"five_prime": dict(UNIFORM_BG)}
    assert terminal_nucleotide_bias_KL_divergence(observed, UNIFORM_BG) == pytest.approx(0.0)
    assert terminal_nucleotide_bias_max_deviation(observed, UNIFORM_BG) == pytest.approx(0.0)


def test_terminal_bias_single_dinucleotide():
    observed = {"five_prime": {d: (1.0 if d == "AA" else 0.0) for d in DINUCS}}
    assert terminal_nucleotide_bias_KL_divergence(observed, UNIFORM_BG) == pytest.approx(4.0)
    assert terminal_nucleotide_bias_max_deviation(observed, UNIFORM_BG) == pytest.approx(15 / 16)


# --- 3.6 library -----------------------------------------------------------------


def test_library_complexity_marginal_rate():
    df = pd.DataFrame({"transcript_id": ["tx"] * 3, "a_site": [10, 20, 30], "count": [1, 1, 2]})
    c = library_complexity_curve(df)
    # D(1) = 3; D(0.95) = 0.95 + 0.95 + (1 - 0.05^2) = 2.8975; margin = 0.05 * 4
    assert c["total_distinct_positions"] == 3
    assert c["marginal_discovery_rate"] == pytest.approx(round((3 - 2.8975) / 0.2, 4))


def test_floss_against_library_aggregate():
    df = pd.DataFrame(
        {
            "transcript_id": ["t1"] * 20 + ["t2"] * 20,
            "read_length": [28] * 20 + [30] * 20,
            "count": [1] * 40,
        }
    )
    f = floss_library_heterogeneity(df, min_reads_per_transcript=20, floss_cutoff=0.3)
    # reference 0.5 / 0.5: each transcript 0.5 * (0.5 + 0.5) = 0.5
    assert f["floss_median"] == pytest.approx(0.5)
    assert f["floss_aberrant_transcript_fraction"] == pytest.approx(1.0)


# --- 1 weights ----------------------------------------------------------------------


def test_weights_apply_even_when_a_read_name_repeats():
    df = pd.DataFrame({"read_name": ["r1", "r1", "r2"], "count": [5, 5, 2]})
    assert _get_weights(df).tolist() == [5, 5, 2]
