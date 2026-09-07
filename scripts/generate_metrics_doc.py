#!/usr/bin/env python
"""Regenerate docs/METRICS.md from the metric registry and scoring spec.

docs/METRICS.md previously drifted badly: it documented four default metric
names that were never emitted and omitted 32 that were. Generating it from the
same objects the code uses makes that impossible.

Usage:  python scripts/generate_metrics_doc.py [--check]
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from RiboMetric.registry import (  # noqa: E402
    METRIC_REGISTRY,
    LEGACY_METRIC_ALIASES,
    HIGHER_BETTER,
    LOWER_BETTER,
    CONTEXT,
)
from RiboMetric.scoring import DEFAULT_SCORING  # noqa: E402

DIRECTION_LABEL = {
    HIGHER_BETTER: "higher is better",
    LOWER_BETTER: "lower is better",
    CONTEXT: "context-dependent",
}

TIER_HEADING = {
    1: ("Tier 1 — is this Ribo-seq?",
        "Gated. A failure here means frame-dependent analysis should not "
        "proceed."),
    2: ("Tier 2 — is it usable for my analysis?",
        "Not gated. These describe whether enough usable signal survives "
        "filtering."),
    3: ("Tier 3 — technical caveats",
        "Not gated. These colour interpretation but do not fail a sample."),
}


def _table(rows, headers):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def build() -> str:
    parts = [
        "# RiboMetric Metrics",
        "",
        "<!-- GENERATED FILE - do not edit by hand.",
        "     Regenerate with: python scripts/generate_metrics_doc.py -->",
        "",
        "Two kinds of key appear in RiboMetric output, and they follow "
        "different rules (see [METRIC_NAMING.md](METRIC_NAMING.md)):",
        "",
        "* **Metrics** (`results[\"metrics\"]`) name the quantity that was "
        "measured, in its natural units and its natural direction. A rate "
        "goes up as the library gets worse.",
        "* **Scores** (in the report, the QC status file and the summary "
        "plot) name the good property, always lie in [0, 1], and are always "
        "higher-is-better.",
        "",
        "## Scored metrics",
        "",
        "Every score below is higher-is-better. The direction flip, where "
        "one is needed, lives in the `method` column and nowhere else.",
        "",
    ]

    for tier in (1, 2, 3):
        heading, blurb = TIER_HEADING[tier]
        entries = [
            (score_key, spec) for score_key, spec in DEFAULT_SCORING.items()
            if spec.get("tier") == tier
        ]
        if not entries:
            continue
        parts += [f"### {heading}", "", blurb, ""]
        rows = []
        for score_key, spec in entries:
            metric = METRIC_REGISTRY[spec["metric"]]
            rows.append([
                f"`{score_key}`",
                f"`{metric.key}`",
                metric.unit,
                DIRECTION_LABEL[metric.direction],
                f"`{spec['method']}`",
                f"{spec['status']['pass']:.2f} / {spec['status']['warn']:.2f}",
            ])
        parts += [_table(rows, [
            "Score", "From metric", "Unit", "Metric direction", "Method",
            "Pass / warn",
        ]), ""]
        for score_key, spec in entries:
            parts.append(f"- **`{score_key}`** — "
                         f"{METRIC_REGISTRY[spec['metric']].summary} "
                         f"{spec.get('decision', '')}".strip())
        parts.append("")

    parts += [
        "## Diagnostics",
        "",
        "Reported as raw measurements with no pass/fail badge. Either their "
        "good direction depends on the protocol, or they describe shape "
        "rather than quality.",
        "",
    ]
    diag_rows = [
        [f"`{spec.key}`", spec.unit, DIRECTION_LABEL[spec.direction],
         spec.summary]
        for spec in sorted(METRIC_REGISTRY.values(), key=lambda s: s.key)
        if spec.scored_as is None
    ]
    parts += [_table(diag_rows, ["Metric", "Unit", "Direction", "Summary"]), ""]

    parts += [
        "## Renamed and removed keys",
        "",
        "Pre-2.0 spellings are reproduced under `results[\"metrics_legacy\"]` "
        "for one minor cycle and removed at v2.1. Where the transform is not "
        "`identity`, the old key held a *different number*: it had a "
        "goodness transform baked in.",
        "",
    ]
    legacy_rows = [
        [f"`{old}`", f"`{new}`",
         {"identity": "same value",
          "one_minus": "old = 1 − new",
          "inverse_1p": "old = 1 / (1 + new)"}[transform]]
        for old, (new, transform) in sorted(LEGACY_METRIC_ALIASES.items())
    ]
    parts += [_table(legacy_rows, ["Pre-2.0 key", "Canonical key",
                                   "Relationship"]), ""]
    return "\n".join(parts) + "\n"


def main() -> int:
    target = Path(__file__).resolve().parents[1] / "docs" / "METRICS.md"
    content = build()
    if "--check" in sys.argv:
        current = target.read_text() if target.exists() else ""
        if current != content:
            print(f"{target} is out of date; run "
                  "scripts/generate_metrics_doc.py")
            return 1
        print(f"{target} is up to date")
        return 0
    target.write_text(content)
    print(f"Wrote {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
