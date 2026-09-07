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
    # Records are keyed by score name and carry the source metric.
    dominance = by_key["periodicity_dominance_score"]
    assert dominance["metric"] == "periodicity_dominance"
    # identity method: score == raw == 2/3 ≈ 0.667 (below pass=0.70 → WARNING)
    assert dominance["score"] == pytest.approx(2.0 / 3.0, abs=1e-6)
    assert dominance["status"] == "WARNING"
    assert dominance["gate"] is True
    uniqueness = by_key["fragment_uniqueness_score"]
    assert uniqueness["metric"] == "duplicate_rate"
    assert uniqueness["raw"] == pytest.approx(0.1, abs=1e-6)
    assert uniqueness["score"] == pytest.approx(0.9, abs=1e-6)
    assert uniqueness["gate"] is False
    assert by_key["usable_read_fraction_score"]["status"] == "PASS"


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
    assert "fragment_uniqueness_score" in keys


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
    config = {"scoring": {"fragment_uniqueness_score": {"gate": True, "status": {"pass": 0.95, "warn": 0.9}}}}
    metrics = {"duplicate_rate": 0.2}  # score 0.8 -> below new pass(0.95)/warn(0.9) -> FAIL
    scored = build_scored_metrics(_results(metrics), config)
    rec = [m for m in scored if m["key"] == "fragment_uniqueness_score"][0]
    assert rec["gate"] is True
    assert rec["status"] == "FAIL"
    assert overall_gate_status(scored) == "FAIL"


# --- S1: periodicity_dominance uses identity, not frame_dominance_rescaled ----

def test_periodicity_dominance_default_uses_identity():
    spec = DEFAULT_SCORING["periodicity_dominance_score"]
    assert spec["method"] == "identity"
    assert spec["metric"] == "periodicity_dominance"
    assert "params" not in spec


def test_periodicity_dominance_pass_threshold_is_0_70():
    assert DEFAULT_SCORING["periodicity_dominance_score"]["status"]["pass"] == pytest.approx(0.70)


def test_periodicity_dominance_does_not_use_frame_dominance_rescaled():
    assert DEFAULT_SCORING["periodicity_dominance_score"]["method"] != "frame_dominance_rescaled"


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
# The single property the whole report rests on. Nothing asserted it before,
# which is why terminal_bias_maxabs shipped inverted for four releases.

from RiboMetric.registry import (          # noqa: E402
    METRIC_REGISTRY,
    SCORED_METRICS,
    HIGHER_BETTER,
    LOWER_BETTER,
    CONTEXT,
)


def _score_for(metric_key, raw, config=None):
    records = build_scored_metrics({"metrics": {metric_key: raw}}, config)
    assert len(records) == 1, f"{metric_key} produced no score"
    return records[0]


# A representative good and bad value for every scored raw metric, in the
# metric's own natural units and direction.
DIRECTION_CASES = [
    ("periodicity_dominance", 0.90, 0.34),
    ("periodicity_information", 0.90, 0.00),
    ("cds_enrichment_ratio", 4.0, 1.0),
    ("recommended_read_proportion", 0.90, 0.05),
    ("uniformity_entropy", 0.95, 0.10),
    ("marginal_position_discovery_rate", 0.02, 0.98),
    ("duplicate_rate", 0.02, 0.95),
    ("rpf_multimapper_rate", 0.02, 0.95),
    ("alignment_multimapper_rate", 0.02, 0.95),
    ("soft_clip_rate_5prime", 0.01, 0.90),
    ("floss_aberrant_transcript_fraction", 0.02, 0.90),
    ("terminal_bias_kl_5prime", 0.02, 1.9),
    ("terminal_bias_kl_3prime", 0.02, 1.9),
    ("terminal_bias_max_deviation_5prime", 0.01, 0.80),
    ("terminal_bias_max_deviation_3prime", 0.01, 0.80),
]


@pytest.mark.parametrize("metric_key,raw_good,raw_bad", DIRECTION_CASES)
def test_every_scored_metric_is_higher_is_better(metric_key, raw_good, raw_bad):
    good = _score_for(metric_key, raw_good)
    bad = _score_for(metric_key, raw_bad)
    assert 0.0 <= good["score"] <= 1.0
    assert 0.0 <= bad["score"] <= 1.0
    assert good["score"] > bad["score"], (
        f"{metric_key}: raw={raw_good} scored {good['score']:.3f} but "
        f"raw={raw_bad} scored {bad['score']:.3f} -- direction is inverted"
    )


@pytest.mark.parametrize("metric_key,raw_good,raw_bad", DIRECTION_CASES)
def test_registry_direction_matches_the_scoring_behaviour(
    metric_key, raw_good, raw_bad
):
    """The registry's declared direction must match what the scorer does.

    This is the check that would have caught the maxabs inversion: the
    registry says terminal_bias_max_deviation_5prime is lower_better, so the
    good sample must be the smaller raw value.
    """
    spec = METRIC_REGISTRY[metric_key]
    assert spec.direction in (HIGHER_BETTER, LOWER_BETTER), (
        f"{metric_key} is scored, so it cannot be direction={spec.direction}"
    )
    if spec.direction == HIGHER_BETTER:
        assert raw_good > raw_bad
    else:
        assert raw_good < raw_bad


def test_every_scored_metric_has_a_direction_case():
    """A new scored metric must declare a direction case above."""
    assert set(SCORED_METRICS.values()) == {c[0] for c in DIRECTION_CASES}


def test_registry_and_scoring_spec_agree():
    """registry.scored_as and DEFAULT_SCORING must name the same pairs."""
    from_registry = {
        spec.scored_as: spec.key
        for spec in METRIC_REGISTRY.values()
        if spec.scored_as
    }
    from_spec = {
        score_key: entry["metric"]
        for score_key, entry in DEFAULT_SCORING.items()
    }
    assert from_registry == from_spec


def test_every_score_key_ends_in_score():
    for score_key in DEFAULT_SCORING:
        assert score_key.endswith("_score"), score_key


def test_no_metric_key_is_also_a_score_key():
    """A key is either a measurement or a score, never both."""
    assert not set(METRIC_REGISTRY) & set(DEFAULT_SCORING)


def test_context_metrics_are_never_scored():
    """Metrics whose good direction depends on protocol must not carry a
    pass/fail badge (METRICS_DESIGN.md decision O3)."""
    for key, spec in METRIC_REGISTRY.items():
        if spec.direction == CONTEXT:
            assert spec.scored_as is None, (
                f"{key} is context-dependent and must not be scored"
            )


def test_periodicity_information_thresholds_are_anchored_to_dominance():
    """Entropy reduction is far more compressive than the dominant-frame
    fraction. The v1.4.0 thresholds (0.60/0.30) demanded dominance ~0.90 to
    pass, silently making this the strictest Tier 1 gate."""
    def information(dominance):
        p = [dominance, (1 - dominance) / 2, (1 - dominance) / 2]
        entropy = -sum(x * math.log2(x) for x in p if x > 0)
        return (math.log2(3) - entropy) / math.log2(3)

    spec = DEFAULT_SCORING["periodicity_information_score"]["status"]
    dom = DEFAULT_SCORING["periodicity_dominance_score"]["status"]

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
            "periodicity_dominance_score": {
                "metric": "periodicity_dominance",
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
