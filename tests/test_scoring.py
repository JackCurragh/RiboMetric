"""
Tests for the unified scoring resolver (RiboMetric/scoring.py).

These pin the anchor maths from docs/METRICS_DESIGN.md so the spec and the
implementation cannot silently drift apart.
"""

import math

import pytest

import pandas as pd

from RiboMetric.scoring import (
    score_value,
    resolve_status,
    build_scored_metrics,
    overall_gate_status,
    DEFAULT_SCORING,
    DEFAULT_STATUS,
)
from RiboMetric.metrics import cds_enrichment_ratio


# --- score method anchors --------------------------------------------------

@pytest.mark.parametrize(
    "raw,expected",
    [
        (1.0 / 3.0, 0.0),    # random frame -> 0
        (2.0 / 3.0, 0.5),    # interior reference -> 0.5
        (1.0, 1.0),          # all in one frame -> 1
        (0.25, 0.0),         # below random floor clips to 0
        (0.73, pytest.approx(0.595, abs=1e-3)),  # ~PASS boundary
    ],
)
def test_frame_dominance_rescaled(raw, expected):
    assert score_value("frame_dominance_rescaled", raw) == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize(
    "raw,expected",
    [(0.0, 1.0), (0.1, 0.9), (0.5, 0.5), (1.0, 0.0), (1.2, 0.0)],
)
def test_one_minus_rate(raw, expected):
    assert score_value("one_minus_rate", raw) == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize("raw,expected", [(0.0, 0.0), (0.5, 0.5), (1.0, 1.0), (1.5, 1.0), (-0.2, 0.0)])
def test_identity_clips(raw, expected):
    assert score_value("identity", raw) == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize(
    "raw,expected",
    [(0.5, 0.0), (1.0, 0.0), (2.0, 0.5), (4.0, 0.75)],
)
def test_enrichment_ratio(raw, expected):
    # Provided for 1C; anchored 1 - 1/E with E<=1 -> 0.
    assert score_value("enrichment_ratio", raw) == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize(
    "raw,expected",
    [(0.0, 1.0), (1.0, 0.5), (2.0, 0.0), (3.0, 0.0)],
)
def test_inverse_linear_default_kl_max(raw, expected):
    # Provided for 1D; default max_value 2.0.
    assert score_value("inverse_linear", raw) == pytest.approx(expected, abs=1e-6)


def test_unknown_method_returns_none():
    assert score_value("nonsense", 0.5) is None


# --- status resolution -----------------------------------------------------

@pytest.mark.parametrize(
    "score,expected",
    [(0.60, "PASS"), (0.75, "PASS"), (0.59, "WARNING"), (0.30, "WARNING"), (0.29, "FAIL"), (0.0, "FAIL")],
)
def test_resolve_status_default_thresholds(score, expected):
    assert resolve_status(score, DEFAULT_STATUS) == expected


def test_resolve_status_none_is_info():
    assert resolve_status(None, DEFAULT_STATUS) == "INFO"


# --- end-to-end over a results dict ----------------------------------------

def _results(metrics):
    return {"metrics": metrics}


def test_build_scored_metrics_handles_scalar_and_global_dict():
    metrics = {
        "periodicity_dominance": {"global": 2.0 / 3.0, "30": 0.9},
        "duplicate_rate": 0.1,
        "recommended_read_proportion": 0.8,
    }
    by_key = {m["key"]: m for m in build_scored_metrics(_results(metrics))}
    # identity method: score == raw == 2/3 ≈ 0.667 (below pass=0.70 → WARNING)
    assert by_key["periodicity_dominance"]["score"] == pytest.approx(2.0 / 3.0, abs=1e-6)
    assert by_key["periodicity_dominance"]["status"] == "WARNING"
    assert by_key["periodicity_dominance"]["gate"] is True
    assert by_key["duplicate_rate"]["score"] == pytest.approx(0.9, abs=1e-6)
    assert by_key["duplicate_rate"]["gate"] is False
    assert by_key["recommended_read_proportion"]["status"] == "PASS"


def test_none_raw_yields_unscored_info():
    metrics = {"marginal_position_discovery_rate": None}
    rec = build_scored_metrics(_results(metrics))[0]
    assert rec["score"] is None
    assert rec["status"] == "INFO"


def test_nan_raw_yields_unscored_info():
    metrics = {"duplicate_rate": float("nan")}
    rec = build_scored_metrics(_results(metrics))[0]
    assert rec["score"] is None and rec["status"] == "INFO"


def test_metrics_absent_from_spec_are_ignored():
    metrics = {"some_unknown_metric": 0.5, "duplicate_rate": 0.0}
    keys = {m["key"] for m in build_scored_metrics(_results(metrics))}
    assert "some_unknown_metric" not in keys
    assert "duplicate_rate" in keys


# --- gate membership drives the overall verdict ----------------------------

def test_overall_verdict_uses_gated_only():
    # Tier-1 (gated) both strong; a Tier-3 caveat fails but must not fail sample.
    metrics = {
        "periodicity_dominance": {"global": 1.0},      # gated, PASS
        "periodicity_information": {"global": 0.9},     # gated, PASS
        "duplicate_rate": 0.95,                          # not gated, FAIL-ish
    }
    scored = build_scored_metrics(_results(metrics))
    assert overall_gate_status(scored) == "PASS"


def test_overall_verdict_fails_on_gated_fail():
    metrics = {
        "periodicity_dominance": {"global": 0.4},       # gated -> low score, FAIL
        "periodicity_information": {"global": 0.9},
    }
    scored = build_scored_metrics(_results(metrics))
    assert overall_gate_status(scored) == "FAIL"


def test_overall_verdict_info_when_no_gated_scores():
    metrics = {"duplicate_rate": 0.1}  # only a non-gated metric
    scored = build_scored_metrics(_results(metrics))
    assert overall_gate_status(scored) == "INFO"


def test_config_override_changes_gate_and_threshold():
    config = {"scoring": {"duplicate_rate": {"gate": True, "status": {"pass": 0.95, "warn": 0.9}}}}
    metrics = {"duplicate_rate": 0.2}  # score 0.8 -> below new pass(0.95)/warn(0.9) -> FAIL
    scored = build_scored_metrics(_results(metrics), config)
    rec = [m for m in scored if m["key"] == "duplicate_rate"][0]
    assert rec["gate"] is True
    assert rec["status"] == "FAIL"
    assert overall_gate_status(scored) == "FAIL"


# --- S1: periodicity_dominance uses identity, not frame_dominance_rescaled ----

def test_periodicity_dominance_default_uses_identity():
    assert DEFAULT_SCORING["periodicity_dominance"]["method"] == "identity"
    assert "params" not in DEFAULT_SCORING["periodicity_dominance"]


def test_periodicity_dominance_pass_threshold_is_0_70():
    assert DEFAULT_SCORING["periodicity_dominance"]["status"]["pass"] == pytest.approx(0.70)


def test_periodicity_dominance_does_not_use_frame_dominance_rescaled():
    assert DEFAULT_SCORING["periodicity_dominance"]["method"] != "frame_dominance_rescaled"


# --- S3: R-O1 self-consistency (uniform reads ⇒ E ≈ 1 ⇒ score ≈ 0) ----------

def _make_uniform_annotated_df():
    """Minimal annotated_read_df with reads distributed uniformly per-nt."""
    # Two transcripts: each nt gets exactly one read.
    # tx1: length 100, cds_start=20, cds_end=80 → CDS body len = 59 nt
    # tx2: length 60,  cds_start=10, cds_end=50 → CDS body len = 39 nt
    # Uniform: one read per nt for all five region buckets combined.
    rows = []
    for tx_id, tx_len, cds_s, cds_e in [("tx1", 100, 20, 80), ("tx2", 60, 10, 50)]:
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
            rows.append({
                "transcript_id": tx_id,
                "a_site": pos,
                "cds_start": cds_s,
                "cds_end": cds_e,
                "transcript_length": tx_len,
                "mRNA_category": cat,
            })
    return pd.DataFrame(rows)


def test_cds_enrichment_self_consistency_uniform_reads():
    df = _make_uniform_annotated_df()
    E = cds_enrichment_ratio(df)
    assert E is not None
    assert E == pytest.approx(1.0, abs=0.01), f"Expected E≈1 for uniform reads, got {E}"
    assert score_value("enrichment_ratio", E) == pytest.approx(0.0, abs=0.02)


def test_cds_enrichment_no_annotation_returns_none():
    df = pd.DataFrame({"a_site": [10, 20], "read_length": [28, 30]})
    assert cds_enrichment_ratio(df) is None


def test_cds_enrichment_empty_df_returns_none():
    assert cds_enrichment_ratio(pd.DataFrame()) is None


# --- direction invariants ---------------------------------------------------
# Regression cover for the v1.4.3 scoring defects found in the 2026-09-07 audit.


def _score_for(key, raw, config=None):
    records = build_scored_metrics({"metrics": {key: raw}}, config)
    assert len(records) == 1, f"{key} is not in the scoring spec"
    return records[0]


@pytest.mark.parametrize("key", [
    "terminal_bias_maxabs_5prime",
    "terminal_bias_maxabs_3prime",
])
def test_terminal_bias_maxabs_is_not_inverted(key):
    """metrics.terminal_nucleotide_bias_max_absolute_metric returns
    ``1 - max|observed - expected|``, i.e. it is already higher-is-better.

    Scoring it with ``one_minus_rate`` inverted it: a perfectly unbiased
    library scored 0.00 FAIL and a severely biased one scored 0.80 PASS.
    """
    clean = _score_for(key, 1.0)      # zero deviation
    biased = _score_for(key, 0.20)    # 0.80 deviation
    assert clean["score"] > biased["score"]
    assert clean["score"] == pytest.approx(1.0)
    assert clean["status"] == "PASS"
    assert biased["status"] == "FAIL"


DIRECTION_CASES = [
    ("periodicity_dominance", 0.90, 0.34),
    ("periodicity_information", 0.90, 0.00),
    ("cds_enrichment_ratio", 4.0, 1.0),
    ("recommended_read_proportion", 0.90, 0.05),
    ("uniformity_entropy", 0.95, 0.10),
    ("marginal_position_discovery_rate", 0.02, 0.98),
    ("duplicate_rate", 0.02, 0.95),
    ("rpf_multimapper_rate", 0.02, 0.95),
    ("multimapper_rate", 0.02, 0.95),
    ("alignment_multimapper_rate", 0.02, 0.95),
    ("unique_rpf_rate", 0.98, 0.05),
    ("soft_clip_rate_5prime", 0.01, 0.90),
    ("terminal_bias_kl_5prime_raw", 0.02, 1.9),
    ("terminal_bias_kl_3prime_raw", 0.02, 1.9),
    ("terminal_bias_maxabs_5prime", 0.99, 0.10),
    ("terminal_bias_maxabs_3prime", 0.99, 0.10),
]


@pytest.mark.parametrize("key,raw_good,raw_bad", DIRECTION_CASES)
def test_every_scored_metric_is_higher_is_better(key, raw_good, raw_bad):
    """The single invariant the whole report rests on: for every scored
    metric, the better library must score higher. This is the check that was
    missing when the maxabs metrics shipped inverted."""
    good = _score_for(key, raw_good)
    bad = _score_for(key, raw_bad)
    assert 0.0 <= good["score"] <= 1.0
    assert 0.0 <= bad["score"] <= 1.0
    assert good["score"] > bad["score"], (
        f"{key}: raw={raw_good} scored {good['score']:.3f} but raw={raw_bad} "
        f"scored {bad['score']:.3f} -- direction is inverted"
    )


def test_every_default_scoring_entry_is_covered_by_the_direction_test():
    """A new scored metric must be added to the direction test above."""
    covered = {case[0] for case in DIRECTION_CASES}
    assert set(DEFAULT_SCORING) == covered


def test_periodicity_information_thresholds_are_anchored_to_dominance():
    """Entropy reduction is far more compressive than the dominant-frame
    fraction. The v1.4.0 thresholds (0.60/0.30) demanded dominance ~0.90 to
    pass, silently making this the strictest Tier 1 gate."""
    def information(dominance):
        p = [dominance, (1 - dominance) / 2, (1 - dominance) / 2]
        entropy = -sum(x * math.log2(x) for x in p if x > 0)
        return (math.log2(3) - entropy) / math.log2(3)

    spec = DEFAULT_SCORING["periodicity_information"]["status"]
    dom = DEFAULT_SCORING["periodicity_dominance"]["status"]

    # A library sitting exactly on the dominance PASS/WARN boundaries must not
    # be failed harder by the information cross-check.
    assert information(dom["pass"]) >= spec["pass"] - 1e-6
    assert information(dom["warn"]) >= spec["warn"] - 1e-6


def test_config_scoring_overrides_reach_the_qc_gate():
    """The gate and the HTML report must resolve the same status from the same
    results and the same config. generate_qc_status previously dropped config,
    so a user override moved the report but not the verdict."""
    from RiboMetric.results_output import evaluate_qc_status

    results = {"metrics": {"periodicity_dominance": {"global": 0.55}}}
    config = {
        "scoring": {
            "periodicity_dominance": {
                "method": "identity",
                "status": {"pass": 0.50, "warn": 0.40},
                "gate": True,
                "tier": 1,
            }
        }
    }
    assert evaluate_qc_status(results, "s")["overall_status"] == "WARNING"
    assert evaluate_qc_status(
        results, "s", None, config
    )["overall_status"] == "PASS"
    resolver = _score_for("periodicity_dominance", 0.55, config)
    assert resolver["status"] == "PASS"
