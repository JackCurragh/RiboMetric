"""The naming and direction contract from docs/METRIC_NAMING.md.

These tests are the reason the contract holds. Before them, nothing in the
codebase asserted which way was up for any metric, and
``terminal_bias_maxabs_5prime`` shipped inverted for four releases.
"""

import subprocess
import sys
from pathlib import Path

import pytest

from RiboMetric.registry import (
    CONTEXT,
    HIGHER_BETTER,
    LEGACY_METRIC_ALIASES,
    LOWER_BETTER,
    METRIC_REGISTRY,
    SCORED_METRICS,
    build_legacy_metrics,
)
from RiboMetric.scoring import DEFAULT_SCORING

REPO_ROOT = Path(__file__).resolve().parents[1]

VALID_UNITS = {
    "fraction",
    "rate",
    "bits",
    "ratio",
    "count",
    "pvalue",
    "index",
}


# --- registry well-formedness ----------------------------------------------


def test_every_spec_is_internally_consistent():
    for key, spec in METRIC_REGISTRY.items():
        assert spec.key == key
        assert spec.unit in VALID_UNITS, f"{key}: unknown unit {spec.unit!r}"
        assert spec.direction in (HIGHER_BETTER, LOWER_BETTER, CONTEXT)
        assert spec.summary, key
        # Sentence case, allowing summaries that legitimately open with a
        # symbol or digit ("5' leader reads relative to...").
        assert not spec.summary[0].islower(), key
        assert spec.summary.endswith("."), key


def test_no_metric_key_carries_a_goodness_transform_in_its_name():
    """A key that says "metric" or "score" is not naming a quantity."""
    for key in METRIC_REGISTRY:
        assert not key.endswith("_metric"), (
            f"{key}: '_metric' meant 'transformed' on some keys and nothing "
            "on others; name the quantity instead"
        )
    # One grandfathered exception, kept because it is a genuine score that
    # predates the split and is not part of the scored set.
    scores = [k for k in METRIC_REGISTRY if k.endswith("_score")]
    assert scores == ["periodicity_information_weighted_score"]


def test_rates_are_lower_better():
    """A "rate" that is higher-is-better would mean the name lies."""
    for key, spec in METRIC_REGISTRY.items():
        if spec.unit == "rate":
            assert spec.direction == LOWER_BETTER, key


def test_context_metrics_are_never_scored():
    for key, spec in METRIC_REGISTRY.items():
        if spec.direction == CONTEXT:
            assert spec.scored_as is None, (
                f"{key} is context-dependent; a pass/fail badge on it would "
                "be a protocol judgement (METRICS_DESIGN.md decision O3)"
            )


# --- registry <-> scoring spec ---------------------------------------------


def test_scored_as_matches_the_scoring_spec():
    from_registry = {
        spec.scored_as: spec.key for spec in METRIC_REGISTRY.values() if spec.scored_as
    }
    from_spec = {k: v["metric"] for k, v in DEFAULT_SCORING.items()}
    assert from_registry == from_spec


def test_every_scored_source_metric_is_registered():
    for score_key, spec in DEFAULT_SCORING.items():
        assert (
            spec["metric"] in METRIC_REGISTRY
        ), f"{score_key} scores {spec['metric']}, which is not registered"


def test_score_keys_and_metric_keys_are_disjoint():
    assert not set(METRIC_REGISTRY) & set(SCORED_METRICS)


def test_lower_better_metrics_use_a_flipping_method():
    """The direction flip must live in `method`, not in the metric function."""
    flipping = {"one_minus_rate", "inverse_linear"}
    for score_key, spec in DEFAULT_SCORING.items():
        direction = METRIC_REGISTRY[spec["metric"]].direction
        if direction == LOWER_BETTER:
            assert spec["method"] in flipping, (
                f"{score_key}: {spec['metric']} is lower-is-better but "
                f"method {spec['method']!r} does not flip it"
            )
        elif direction == HIGHER_BETTER:
            assert spec["method"] not in flipping, (
                f"{score_key}: {spec['metric']} is higher-is-better but "
                f"method {spec['method']!r} flips it"
            )


# --- legacy shim ------------------------------------------------------------


def test_every_legacy_alias_points_at_a_registered_metric():
    for old_key, (new_key, _transform) in LEGACY_METRIC_ALIASES.items():
        assert new_key in METRIC_REGISTRY, f"{old_key} -> {new_key}"


def test_legacy_aliases_do_not_collide_with_canonical_keys():
    """An old spelling must not shadow a canonical metric with a different
    value. `terminal_bias_kl_5prime` is the trap: it used to hold 1/(1+KL)
    and now holds KL, so it may only live in metrics_legacy."""
    for old_key, (new_key, transform) in LEGACY_METRIC_ALIASES.items():
        if old_key in METRIC_REGISTRY and old_key != new_key:
            pytest.fail(f"{old_key} is both a legacy alias and a canonical key")
        if old_key == new_key:
            assert transform != "identity", f"{old_key} maps to itself unchanged; drop the alias"


def test_legacy_metrics_reproduce_the_pre_2_0_values():
    metrics = {
        "terminal_bias_kl_5prime": 1.0,
        "terminal_bias_max_deviation_5prime": 0.25,
        "read_length_iqr_fraction": 0.4,
        "rpf_multimapper_rate": 0.2,
        "cds_coverage": 0.5,
    }
    legacy = build_legacy_metrics(metrics)

    # 1/(1 + KL)
    assert legacy["terminal_bias_kl_5prime_score"] == pytest.approx(0.5)
    assert legacy["terminal_bias_kl_5prime"] == pytest.approx(0.5)
    # raw bits passed through unchanged
    assert legacy["terminal_bias_kl_5prime_raw"] == pytest.approx(1.0)
    # 1 - deviation
    assert legacy["terminal_bias_maxabs_5prime"] == pytest.approx(0.75)
    # 1 - spread
    assert legacy["read_length_distribution_IQR_metric"] == pytest.approx(0.6)
    # plain aliases
    assert legacy["multimapper_rate"] == pytest.approx(0.2)
    assert legacy["unique_rpf_rate"] == pytest.approx(0.8)
    assert legacy["CDS_coverage_metric"] == pytest.approx(0.5)


def test_legacy_metrics_skip_absent_and_non_numeric_sources():
    assert build_legacy_metrics({}) == {}
    assert build_legacy_metrics({"cds_coverage": None}) == {}
    assert build_legacy_metrics({"cds_coverage": {"global": 0.5}}) == {}


# --- documentation cannot drift --------------------------------------------


def test_metrics_doc_is_generated_from_the_registry():
    """docs/METRICS.md documented four metrics that were never emitted and
    omitted 32 that were. It is generated now, so it cannot drift."""
    result = subprocess.run(
        [sys.executable, "scripts/generate_metrics_doc.py", "--check"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
