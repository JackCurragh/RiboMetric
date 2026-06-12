#!/usr/bin/env python3
"""
Before/after assessment for the unified scoring redesign (Phase 1A/1B).

Reads a RiboMetric JSON output (or any dict with a top-level "metrics" key) and
prints, per scored metric, the raw value alongside the OLD displayed score +
status and the NEW anchored score + status from the resolver.

This needs no pipeline run and no heavy deps -- it imports only RiboMetric.scoring
(which is Python 3.9-clean), so it works in any environment.

Usage:
    python3 scripts/score_before_after.py path/to/<sample>_RiboMetric.json
"""

import json
import sys
from pathlib import Path

# Allow running from the repo root without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from RiboMetric.scoring import build_scored_metrics  # noqa: E402


# --- OLD scoring behaviour, reproduced for honest comparison ----------------
# Old summary "score" = raw value (max_mins ranges were [0,1] no-ops, so
# normalise_score was identity). Old status = html_report._metric_status:
#   higher-is-better: >=0.7 pass, >=0.5 warn, else fail
#   lower-is-better : <=0.2 pass, <=0.5 warn, else fail
OLD_LOWER_IS_BETTER = {
    "duplicate_rate", "multimapper_rate", "rpf_multimapper_rate",
    "alignment_multimapper_rate", "soft_clip_rate_5prime", "disome_proportion",
    "marginal_position_discovery_rate",
}


def old_status(key, raw):
    if raw is None:
        return "info"
    if key in OLD_LOWER_IS_BETTER:
        if raw <= 0.2:
            return "pass"
        if raw <= 0.5:
            return "warn"
        return "fail"
    if raw >= 0.7:
        return "pass"
    if raw >= 0.5:
        return "warn"
    return "fail"


def _extract(metric_value):
    if isinstance(metric_value, dict):
        return metric_value.get("global")
    if isinstance(metric_value, (int, float)) and not isinstance(metric_value, bool):
        return float(metric_value)
    return None


def main(path):
    data = json.loads(Path(path).read_text())
    results = data.get("results", data)
    scored = build_scored_metrics(results)

    header = f"{'metric':<34}{'raw':>9}{'OLD score':>11}{'OLD':>9}   {'NEW score':>10}{'NEW':>9}  {'tier':>4} {'gate':>6}"
    print(header)
    print("-" * len(header))
    for m in scored:
        key = m["key"]
        raw = m["raw"]
        raw_s = f"{raw:.3f}" if isinstance(raw, (int, float)) else "n/a"
        old_sc = f"{raw:.3f}" if isinstance(raw, (int, float)) else "n/a"
        new_sc = f"{m['score']:.3f}" if m["score"] is not None else "n/a"
        chg = " *" if (
            isinstance(raw, (int, float)) and m["score"] is not None
            and abs(m["score"] - raw) > 0.005
        ) else ""
        os_ = old_status(key, raw)
        ns_ = m["status"].lower()
        sflag = "  <-- status changed" if os_ != ns_ else ""
        print(
            f"{key:<34}{raw_s:>9}{old_sc:>11}{os_:>9}   {new_sc:>10}{ns_:>9}  "
            f"{str(m['tier']):>4} {str(m['gate']):>6}{chg}{sflag}"
        )

    from RiboMetric.scoring import overall_gate_status
    print()
    print(f"NEW overall verdict (gated/Tier-1 only): {overall_gate_status(scored)}")
    print("'*' = score value changed; raw values are unchanged and shown for both.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
