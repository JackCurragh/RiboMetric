"""
Single source of truth for turning raw RiboMetric quantities into anchored
0-1 scores and PASS/WARNING/FAIL statuses.

This module implements the scoring contract from ``docs/METRICS_DESIGN.md``:

* every score is in [0, 1] and higher-is-better;
* the raw value travels alongside the score;
* one status resolver is shared by the summary cards, the per-metric badges,
  and the QC gate;
* "scored" does not mean "gated" -- only metrics flagged ``gate: true`` count
  towards the overall QC verdict.

Phase 1A/1B scope: the resolver plus the scores that need no new data
(frame-dominance rescale, de-inflated periodicity information, recommended-read
proportion, mapping-hygiene rates, coverage entropy, saturation). CDS
enrichment (1C) and terminal-bias linear scores (1D) plug in here later by
adding a method and a config entry -- no consumer changes required.
"""

from __future__ import annotations

import math
from typing import Any, Callable, Dict, List, Optional


def _clip(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


# --------------------------------------------------------------------------
# Score methods: raw quantity -> [0, 1], higher-is-better.
# Each takes (raw, **params) and returns a float, or None when the raw value
# is missing/not applicable (so the metric is shown without a score).
# --------------------------------------------------------------------------


def _m_identity(raw: float) -> float:
    """Raw is already a 0-1, higher-is-better fraction."""
    return _clip(float(raw))


def _m_one_minus_rate(raw: float) -> float:
    """Lower-is-better rate in [0, 1] -> higher-is-better score."""
    return _clip(1.0 - float(raw))


def _m_frame_dominance_rescaled(raw: float, zero: float = 1.0 / 3.0, one: float = 1.0) -> float:
    """Rescale dominant-frame fraction against its random-frame floor.

    Dominance has a real lower bound at 1/3 (random); treating 0 as the floor
    overstates quality. ``(d - 1/3) / (2/3)`` so 1/3 -> 0, 2/3 -> 0.5, 1 -> 1.
    """
    span = one - zero
    if span <= 0:
        return 0.0
    return _clip((float(raw) - zero) / span)


def _m_inverse_linear(raw: float, max_value: float = 2.0) -> float:
    """Lower-is-better unbounded quantity (e.g. KL bits) -> linear score.

    ``1 - raw/max_value`` clipped. Replaces the uninterpretable ``1/(1+x)``.
    Used by terminal-bias KL in Phase 1D; defined here so the method exists.
    """
    if max_value <= 0:
        return 0.0
    return _clip(1.0 - float(raw) / max_value)


def _m_enrichment_ratio(raw: float) -> float:
    """Enrichment ratio E (>=0) -> ``1 - 1/E``. E<=1 -> 0, E=2 -> 0.5.

    Used by CDS enrichment in Phase 1C; defined here so the method exists.
    """
    e = float(raw)
    if e <= 1.0:
        return 0.0
    return _clip(1.0 - 1.0 / e)


SCORE_METHODS: Dict[str, Callable[..., float]] = {
    "identity": _m_identity,
    "one_minus_rate": _m_one_minus_rate,
    "frame_dominance_rescaled": _m_frame_dominance_rescaled,
    "inverse_linear": _m_inverse_linear,
    "enrichment_ratio": _m_enrichment_ratio,
}


# --------------------------------------------------------------------------
# Default scoring spec. Lives in code so the resolver works without config,
# but ``config["scoring"]`` (if present) is merged over the top per metric.
#   method     : key into SCORE_METHODS
#   params     : kwargs for the method
#   status     : {"pass": float, "warn": float} thresholds on the SCORE
#   gate       : whether this metric counts toward the overall QC verdict
#   tier       : 1 (identity) / 2 (usability) / 3 (caveats)
#   decision   : one-line consequence of a bad value
# --------------------------------------------------------------------------

DEFAULT_STATUS = {"pass": 0.60, "warn": 0.30}

DEFAULT_SCORING: Dict[str, Dict[str, Any]] = {
    # ---- Tier 1: is this Ribo-seq-like? (gated) ------------------------
    "periodicity_dominance": {
        "method": "identity",
        "status": {"pass": 0.70, "warn": 0.50},
        "gate": True,
        "tier": 1,
        "decision": "Low score: weak triplet structure; P-site assignment and "
        "ORF calling unreliable.",
    },
    "cds_enrichment_ratio": {
        "method": "enrichment_ratio",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": True,
        "tier": 1,
        "decision": "Low score: reads not enriched over coding sequence; library "
        "may reflect degradation, RNA contamination, or poor nuclease "
        "protection.",
    },
    "periodicity_information": {
        "method": "identity",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": True,
        "tier": 1,
        "decision": "Cross-check on frame dominance; large disagreement signals "
        "frame mixing or unstable offsets.",
    },
    # ---- Tier 2: usable for my analysis? (not gated) -------------------
    "recommended_read_proportion": {
        "method": "identity",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 2,
        "decision": "Low score: little of the library survives recommended-read "
        "filtering for frame-sensitive work.",
    },
    "uniformity_entropy": {
        "method": "identity",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 2,
        "decision": "Low score: coverage dominated by a few hotspots; broad "
        "quantification may be unreliable.",
    },
    "marginal_position_discovery_rate": {
        "method": "one_minus_rate",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 2,
        "decision": "Low score: library under-sequenced; more reads would "
        "discover substantially more positions.",
    },
    # ---- Tier 3: technical caveats (not gated) -------------------------
    "duplicate_rate": {
        "method": "one_minus_rate",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 3,
        "decision": "Low score: usable molecule diversity much lower than read "
        "depth suggests (protocol-dependent).",
    },
    "rpf_multimapper_rate": {
        "method": "one_minus_rate",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 3,
        "decision": "Low score: reduced confidence in locus/transcript-level " "quantification.",
    },
    "multimapper_rate": {
        "method": "one_minus_rate",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 3,
        "decision": "Low score: reduced confidence in locus/transcript-level " "quantification.",
    },
    "alignment_multimapper_rate": {
        "method": "one_minus_rate",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 3,
        "decision": "Low score: many alignment rows have evidence of another "
        "reported alignment.",
    },
    "unique_rpf_rate": {
        "method": "identity",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 3,
        "decision": "Low score: few fragments map uniquely.",
    },
    "soft_clip_rate_5prime": {
        "method": "one_minus_rate",
        "status": {"pass": 0.60, "warn": 0.30},
        "gate": False,
        "tier": 3,
        "decision": "Low score: 5' read ends frequently clipped; offset and "
        "terminal-bias interpretation may be unreliable.",
    },
    # Terminal bias (S4): raw KL bits via inverse_linear; 1/(1+KL) _score
    # variants are retired from the scored set and kept only as legacy keys.
    "terminal_bias_kl_5prime_raw": {
        "method": "inverse_linear",
        "params": {"max_value": 2.0},
        "status": {"pass": 0.70, "warn": 0.40},
        "gate": False,
        "tier": 3,
        "decision": "Low score: 5' terminal sequence bias may distort count "
        "quantification; consider correction.",
    },
    "terminal_bias_kl_3prime_raw": {
        "method": "inverse_linear",
        "params": {"max_value": 2.0},
        "status": {"pass": 0.70, "warn": 0.40},
        "gate": False,
        "tier": 3,
        "decision": "Low score: 3' terminal sequence bias may distort count "
        "quantification; consider correction.",
    },
    "terminal_bias_maxabs_5prime": {
        "method": "one_minus_rate",
        "status": {"pass": 0.70, "warn": 0.40},
        "gate": False,
        "tier": 3,
        "decision": "Low score: at least one 5' terminal dinucleotide is strongly "
        "over- or under-represented.",
    },
    "terminal_bias_maxabs_3prime": {
        "method": "one_minus_rate",
        "status": {"pass": 0.70, "warn": 0.40},
        "gate": False,
        "tier": 3,
        "decision": "Low score: at least one 3' terminal dinucleotide is strongly "
        "over- or under-represented.",
    },
}


def get_scoring_spec(config: Optional[Dict[str, Any]] = None) -> Dict[str, Dict[str, Any]]:
    """Return the effective scoring spec: defaults with ``config["scoring"]``
    merged over the top per metric."""
    spec = {k: dict(v) for k, v in DEFAULT_SCORING.items()}
    if config and isinstance(config.get("scoring"), dict):
        for key, override in config["scoring"].items():
            if key in spec and isinstance(override, dict):
                spec[key].update(override)
            elif isinstance(override, dict):
                spec[key] = dict(override)
    return spec


def _extract_raw(metric_value: Any) -> Optional[float]:
    """Pull a scalar raw value from a metric entry (scalar or dict-with-global)."""
    if isinstance(metric_value, dict):
        if "global" in metric_value:
            metric_value = metric_value["global"]
        else:
            return None
    if isinstance(metric_value, bool):
        return None
    if isinstance(metric_value, (int, float)) and not math.isnan(float(metric_value)):
        return float(metric_value)
    return None


def score_value(
    method: str, raw: float, params: Optional[Dict[str, Any]] = None
) -> Optional[float]:
    """Apply a named score method to a raw value."""
    fn = SCORE_METHODS.get(method)
    if fn is None:
        return None
    try:
        return float(fn(raw, **(params or {})))
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def resolve_status(score: Optional[float], status_thresholds: Dict[str, float]) -> str:
    """Map a 0-1 higher-is-better score to PASS/WARNING/FAIL (or INFO if None)."""
    if score is None:
        return "INFO"
    if score >= status_thresholds.get("pass", DEFAULT_STATUS["pass"]):
        return "PASS"
    if score >= status_thresholds.get("warn", DEFAULT_STATUS["warn"]):
        return "WARNING"
    return "FAIL"


def build_scored_metrics(
    results_dict: Dict[str, Any],
    config: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Produce the canonical scored-metric records consumed everywhere.

    Returns one record per metric present in both the results and the scoring
    spec::

        {key, raw, score, status, gate, tier, decision}

    Metrics whose raw value is missing or not applicable (e.g. a ``None``
    saturation rate) are returned with ``score=None`` and ``status="INFO"`` so
    they can be shown without a misleading 0%.
    """
    spec = get_scoring_spec(config)
    metrics = results_dict.get("metrics", {})
    records: List[Dict[str, Any]] = []
    for key, mspec in spec.items():
        if key not in metrics:
            continue
        raw = _extract_raw(metrics[key])
        score = score_value(mspec["method"], raw, mspec.get("params")) if raw is not None else None
        status_thresholds = mspec.get("status", DEFAULT_STATUS)
        records.append(
            {
                "key": key,
                "raw": raw,
                "score": score,
                "status": resolve_status(score, status_thresholds),
                "gate": bool(mspec.get("gate", False)),
                "tier": mspec.get("tier"),
                "decision": mspec.get("decision", ""),
            }
        )
    return records


def overall_gate_status(scored_metrics: List[Dict[str, Any]]) -> str:
    """Overall QC verdict from the gated subset only.

    FAIL if any gated metric fails, else WARNING if any warns, else PASS.
    Returns INFO if no gated metric produced a score.
    """
    gated = [m for m in scored_metrics if m["gate"] and m["score"] is not None]
    if not gated:
        return "INFO"
    statuses = {m["status"] for m in gated}
    if "FAIL" in statuses:
        return "FAIL"
    if "WARNING" in statuses:
        return "WARNING"
    return "PASS"
