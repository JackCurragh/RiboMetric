# RiboMetric Scoring — Phase 1 Sonnet Task List

**For:** an implementation agent (Sonnet) picking this up cold.
**Design source of truth:** `docs/METRICS_DESIGN.md`. Read it before starting.
**Status of Phase 1A/1B:** already implemented (resolver `RiboMetric/scoring.py`,
`scoring:` block in `config.yml`, rewired summary plot / report context / QC gate,
`sqrt` dropped in `metrics.calculate_3nt_periodicity_score`, tests in
`tests/test_scoring.py`, assessment script `scripts/score_before_after.py`).
These tasks finish Phase 1.

---

## Invariants — do not violate (every task)

1. **All scores are 0–1 and higher-is-better.**
2. **The raw value travels with the score** — the report shows both; never drop the
   raw quantity in natural units (in-frame %, KL bits, duplicate rate, ratio E).
3. **Gate ≠ universal.** Only metrics with `gate: true` (Tier-1 identity) decide the
   overall verdict. Tier-2/3 surface as caveats, never as a sample-level FAIL.
4. **Floors are shown, not scaled.** Where a metric has a meaningful non-zero floor,
   present the raw value and display the floor as context — do NOT rescale it into
   the number. **The periodicity correction in `METRICS_DESIGN.md` §5 is the worked
   example of this trap; re-read it so you don't repeat it.**

**Acceptance check for every task:** `pytest tests/test_scoring.py` and the other
collectible suites pass, and `python3 scripts/score_before_after.py <a RiboMetric.json>`
runs clean. After S0, also confirm one HTML report renders.

---

## Resolved design decisions (use these as given — already decided by Opus/owner)

### R-O1 — CDS enrichment denominator/numerator (for S3)
Confirmed against `modules.assign_mRNA_category` (the five A-site buckets:
`<cds_start`→five_leader, `==cds_start`→start_codon[1 pos],
`cds_start<a<cds_end`→CDS, `==cds_end`→stop_codon[1 pos], `>cds_end`→three_trailer).

* **Observed CDS-body fraction** = reads in the `"CDS"` bucket ÷ reads in *all five*
  buckets (this is exactly today's `prop_reads_CDS` global).
* **Per-transcript lengths:**
  * `CDS_body_len_t = max(0, cds_end_t − cds_start_t − 1)` — the interior positions
    only, excluding the single start_codon (`==cds_start`) and stop_codon
    (`==cds_end`) positions. (Identical to `cds_coverage`'s `interior_cds_length`.)
  * `eligible_mRNA_len_t = transcript_length_t` — the full transcript; every A-site
    position lands in one of the five buckets, so the eligible length that matches
    the observed denominator is the whole transcript.
* **Weight** `reads_t` = total reads on transcript t (sum over all five buckets).
* **Eligibility:** include transcript t only if `transcript_length_t > 0` **and**
  `CDS_body_len_t > 0` **and** t appears in the reads.
* **Formula:**
  `expected_CDS_fraction = Σ_t reads_t·(CDS_body_len_t / eligible_mRNA_len_t) / Σ_t reads_t`
  `E = observed_CDS_body_fraction / expected_CDS_fraction`
* **Self-consistency (must hold, add as a test):** if reads are distributed uniformly
  per nucleotide, observed ≈ expected ⇒ `E ≈ 1` ⇒ score ≈ 0.

### R-O2 — Terminal-bias `KL_max` (for S4)
`KL_max = 2.0` bits, with the `inverse_linear` method (`score = clip(1 − KL/KL_max)`).
Rationale: KL = 1 bit → score 0.5, preserving the midpoint of the old `1/(1+KL)`
transform so verdicts don't lurch; KL ≥ 2 bits → 0 (one/few dinucleotides dominate
the terminal composition = gross bias). **Provisional** — flagged for empirical
recalibration in Phase 2 (likely tightening toward ~1.0). Surface raw KL in bits.

### R-O3 — Context-dependent diagnostics (for S5a)
Any diagnostic whose "good direction" depends on experiment type gets a **neutral
caption keyed off `classify_library_type`, never a pass/fail badge**:
* `disome_proportion`: show raw value + caption — monosome library → "high values may
  indicate disome contamination or incomplete size selection"; disome experiment →
  "expected/desirable"; unknown → present both readings.
* The context strip declares the inferred library type up front so the caption has
  something to key on.

---

## Tasks

> Critical path: **S0 → S1, S3 → S2, S4 → S5 → S6**. S1/S3 are independent; S2/S4
> are independent. S6 is blocked on the Phase-2 corpus (separate, owner-led).

### S0 — Fix the dev environment *(blocks end-to-end verification of all below)*
The repo uses Python ≥3.10 `X | Y` union syntax (e.g. `metrics.py:178`, `plots.py:890`)
but the default local interpreter is 3.9, and the 3.11+ interpreters lack
scipy/plotly. Read `setup.py`/`pyproject.toml` to confirm the target Python, then set
up a venv with the full deps (scipy, plotly, pandas, numpy, pyyaml, jinja2, pytest).
**Done when:** `pytest` *collects all* test files (today `test_metrics`, `test_plots`,
`test_metrics_improved` fail to collect on 3.9) and one HTML report renders.

### S1 — Apply the periodicity correction
In `config.yml` `scoring:` and `RiboMetric/scoring.py`: switch `periodicity_dominance`
from `frame_dominance_rescaled` to `method: identity`; set status on the raw fraction
(`pass: 0.70, warn: 0.50`). Leave `periodicity_information` as `identity` (its floor is
already 0). Update `tests/test_scoring.py`: the `frame_dominance_rescaled` anchor tests
become an identity check; keep one test asserting dominance does **not** use the
rescale. Re-run `score_before_after.py` — periodicity should now show the raw in-frame
fraction (e.g. 0.79), not 0.69.

### S2 — Periodicity 1/3 baseline marker *(depends on S1)*
In the periodicity plot/report, add a `0.333` reference line + caption
("random-frame baseline"). Read-only; no score change. Goal: the floor is visible so
the user reads the raw value against it.

### S3 — 1C CDS enrichment *(uses R-O1)*
Add metric `cds_enrichment_ratio` computed per R-O1 (compute in `qc.py`, alongside
`prop_reads_CDS`). Score with the existing `enrichment_ratio` method; `gate: true`,
`tier: 1`. Add an annotation guard: with no annotation, skip gracefully (raw `None`,
status INFO) — never crash, never 0. Keep raw `prop_reads_CDS` visible as a diagnostic.
Add the R-O1 self-consistency test (uniform reads ⇒ E≈1 ⇒ score≈0).

### S4 — 1D technical caveats *(uses R-O2)*
(a) Terminal bias: route `terminal_bias_kl_5prime_raw` / `_3prime_raw` (KL in bits)
through `inverse_linear` with `max_value: 2.0`; show raw bits beside the score; retire
the `1/(1+KL)` `_score` variants from the scored set. (b) Max-abs terminal deviation:
score `1 − maxdiff`. (c) Coverage: `uniformity_entropy` stays `identity` (already wired).
Add all to the `scoring:` config with tier/gate.

### S5 — 1E report restructure
(a) **Reclassify diagnostics** *(uses R-O3)*: remove `read_length_distribution_*`
shape proxies (IQR/CV/maxprop/bimodality/normality), `disome_proportion`,
`classify_library_type`, and start/stop/codon metrics from the **scored** set; render
them as plots/labels/captions only.
(b) **Layout** *(flag for owner review before merge — user-facing)*: implement the
Tier 1/2/3 structure + the context strip (library type, mapping mode genomic vs
transcriptomic, annotation source, dominant read-length, total/mapped reads) +
raw-beside-score display, per `METRICS_DESIGN.md` §4. Use the three-way framing:
failed identity / passed-but-weak / passed-with-caveats.

### S6 — Phase-2 data collection (mechanical only) *(blocked on corpus, owner-led)*
Once a labelled corpus exists, run RiboMetric over it and emit per-metric raw-value
distributions + a p05/p50/p95 table. **No threshold decisions** — tabulate only.

---

## Deferred (do NOT start without owner sign-off)
- **D3** Align the `evaluate` CLI to the resolver (currently legacy threshold path;
  careful — external-YAML override semantics + exit codes).
- **D4** Reconcile the two periodicity-information "global" definitions
  (`information_metric_cutoff` vs `read_frame_information_weighted_score`).
- **D5** Drop/replace the normality metric (saturates to ~1.0; now a diagnostic).
