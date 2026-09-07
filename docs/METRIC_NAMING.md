# Metric Naming and Direction

**Status:** implemented. This document defines the contract;
[`RiboMetric/registry.py`](../RiboMetric/registry.py) is the machine-readable
form of it, [`tests/test_registry.py`](../tests/test_registry.py) enforces it,
and [`METRICS.md`](METRICS.md) is generated from it. For how raw values become
scores see [`METRICS_DESIGN.md`](METRICS_DESIGN.md).

Outcome on the reference sample: **54 metric keys → 39**, plus 15 scores split
into their own namespace. 25 pre-2.0 spellings are reproduced under
`results["metrics_legacy"]`, verified to match a pre-rename run of the same
BAM exactly.

---

## 1. The problem

`METRICS_DESIGN.md` fixed the *score* contract: every score is 0–1 and
higher-is-better. It did not fix the *name* contract. The result is that a
metric's name and the number stored under it frequently disagree about which
direction is good, and there is no way to tell which from the key alone.

Three distinct failure modes are live in v1.4.3:

**(a) The name says "badness", the stored value is goodness.**
`terminal_nucleotide_bias_max_absolute_metric()` returns `1 - max|obs − exp|`.
The key is `terminal_bias_maxabs_5prime`. A reader — or a scorer author — sees
"bias" and reasonably assumes high = more bias. S4 wired it to
`one_minus_rate` on exactly that assumption, and the score inverted: a
perfectly unbiased library scores 0.00 FAIL, a severely biased one scores 0.80
PASS. This is not a typo; it is what the naming makes likely.

The same shape appears throughout the read-length family:

| key | function actually returns |
|---|---|
| `read_length_distribution_IQR_metric` | `1 − IQR/(P90−P10)` |
| `read_length_distribution_coefficient_of_variation_metric` | `1/(1 + CV)` |
| `read_length_distribution_bimodality_metric` | `1/(1 + BC)` |
| `read_length_distribution_normality_metric` | `1 − p` |
| `terminal_bias_maxabs_5prime` / `_3prime` | `1 − max deviation` |
| `terminal_bias_kl_5prime` / `_score` | `1/(1 + KL)` |

In every case the key names a badness and holds a goodness. In every case the
inversion is baked into the metric function, where a downstream scorer cannot
see it.

**(b) Three names for one number.** `terminal_bias_kl_5prime`,
`terminal_bias_kl_5prime_score` and
`terminal_nucleotide_bias_distribution_5_prime_metric` are all
`1/(1 + KL)`. All ten `CDS_coverage*` / `cds_coverage*` keys are five numbers.
A consumer cannot tell which key is canonical.

**(c) `_metric` and `_raw` suffixes carry no consistent meaning.**
`terminal_bias_kl_5prime_raw` is the honest raw quantity. `*_metric` sometimes
means "transformed" (`read_length_distribution_IQR_metric`) and sometimes
means nothing (`CDS_coverage_metric` == `cds_coverage`).

---

## 2. The rule

> **A metric key names the physical quantity that was measured, in its natural
> units and natural direction. A score key names the good property, always
> lies in [0, 1], and is always higher-is-better.**

Concretely:

1. **`results["metrics"]` holds raw measurements only.** No `1 − x`, no
   `1/(1 + x)`, no clipping-to-goodness inside a metric function. A rate is a
   rate. Bias is bias. KL is bits. If the natural direction of the quantity is
   "lower is better", the key says so plainly (`duplicate_rate`,
   `terminal_bias_max_deviation_5prime`) and the number goes up when the
   library gets worse.

2. **`results["scores"]` holds derived scores only.** Every key ends in
   `_score`, every value is 0–1, higher-is-better, and the key names the
   *good* property — `fragment_uniqueness_score`, not `duplicate_rate_score`.
   The direction is then legible from the name without consulting the config.

3. **The direction flip lives in exactly one place**: the `method` field in
   `scoring.py` / `config.yml`. That is the only code allowed to know that a
   quantity is lower-is-better.

4. **One number, one key.** Aliases are removed, not accumulated.

5. **The invariant is tested, not assumed.** A registry test asserts that every
   scored metric declares a direction, that every score is in [0, 1], and that
   perturbing a raw value in the bad direction never raises its score.

Rule 3 is what actually prevents a repeat of the `maxabs` inversion. Once no
metric function inverts anything, "which way is up?" has a single answer per
metric, declared next to the threshold that uses it.

---

## 3. What changes

### 3.1 De-invert (metric function changes direction; key renamed to match)

These are the real behaviour changes. The stored number changes meaning, so the
key must change with it.

| current key | current value | new metric key | new value | new score key | method |
|---|---|---|---|---|---|
| `terminal_bias_maxabs_5prime` | `1 − maxdev` | `terminal_bias_max_deviation_5prime` | `maxdev` | `terminal_evenness_maxdev_5prime_score` | `one_minus_rate` |
| `terminal_bias_maxabs_3prime` | `1 − maxdev` | `terminal_bias_max_deviation_3prime` | `maxdev` | `terminal_evenness_maxdev_3prime_score` | `one_minus_rate` |
| `terminal_bias_kl_5prime_raw` | KL bits | `terminal_bias_kl_5prime` | KL bits | `terminal_evenness_kl_5prime_score` | `inverse_linear` |
| `terminal_bias_kl_3prime_raw` | KL bits | `terminal_bias_kl_3prime` | KL bits | `terminal_evenness_kl_3prime_score` | `inverse_linear` |
| `read_length_distribution_IQR_metric` | `1 − IQR/range` | `read_length_iqr_fraction` | `IQR/(P90−P10)` | — (diagnostic) | — |
| `read_length_distribution_coefficient_of_variation_metric` | `1/(1+CV)` | `read_length_cv` | `CV` | — (diagnostic) | — |
| `read_length_distribution_bimodality_metric` | `1/(1+BC)` | `read_length_bimodality_coefficient` | `BC` | — (diagnostic) | — |
| `read_length_distribution_normality_metric` | `1 − p` | `read_length_normality_pvalue` | `p` | — (diagnostic) | — |

Note that the whole read-length family loses its scores. Per
`METRICS_DESIGN.md` §Phase 1E these are already diagnostics; a `1/(1+x)`
"score" for them was never anchored to anything, and none of them is currently
gated. Removing the transform removes an uninterpretable number rather than
information.

### 3.2 Rename only (value unchanged, name now states the direction)

| current key | new metric key | new score key | method |
|---|---|---|---|
| `duplicate_rate` | `duplicate_rate` | `fragment_uniqueness_score` | `one_minus_rate` |
| `rpf_multimapper_rate` | `rpf_multimapper_rate` | `rpf_unique_mapping_score` | `one_minus_rate` |
| `alignment_multimapper_rate` | `alignment_multimapper_rate` | `alignment_unique_mapping_score` | `one_minus_rate` |
| `soft_clip_rate_5prime` | `soft_clip_rate_5prime` | `terminal_integrity_5prime_score` | `one_minus_rate` |
| `marginal_position_discovery_rate` | `marginal_position_discovery_rate` | `library_saturation_score` | `one_minus_rate` |
| `floss_aberrant_transcript_fraction` | `floss_aberrant_transcript_fraction` | `footprint_homogeneity_score` | `one_minus_rate` |
| `periodicity_dominance` | `periodicity_dominance` | `periodicity_dominance_score` | `identity` |
| `periodicity_information` | `periodicity_information` | `periodicity_information_score` | `identity` |
| `cds_enrichment_ratio` | `cds_enrichment_ratio` | `cds_enrichment_score` | `enrichment_ratio` |
| `uniformity_entropy` | `uniformity_entropy` | `coverage_uniformity_score` | `identity` |
| `recommended_read_proportion` | `recommended_read_proportion` | `usable_read_fraction_score` | `identity` |
| `read_length_distribution_maxprop_metric` | `read_length_max_proportion` | — (diagnostic) | — |

The raw keys in this group are already honest — a rate that goes up as the
library gets worse, under a name that says "rate". Only the score gains a name.

### 3.3 Delete (alias of another key)

| removed key | canonical key |
|---|---|
| `terminal_bias_kl_5prime`, `terminal_bias_kl_5prime_score`, `terminal_nucleotide_bias_distribution_5_prime_metric` | `terminal_evenness_kl_5prime_score` |
| `terminal_bias_kl_3prime`, `terminal_bias_kl_3prime_score`, `terminal_nucleotide_bias_distribution_3_prime_metric` | `terminal_evenness_kl_3prime_score` |
| `terminal_nucleotide_bias_max_absolute_metric_5_prime_metric` | `terminal_bias_max_deviation_5prime` |
| `terminal_nucleotide_bias_max_absolute_metric_3_prime_metric` | `terminal_bias_max_deviation_3prime` |
| `multimapper_rate` | `rpf_multimapper_rate` |
| `unique_rpf_rate` | `1 − rpf_multimapper_rate`; keep the rate, drop the complement |
| `CDS_coverage_metric` and its four `_inframe_*` / `_not_inframe_*` variants | `cds_coverage[_inframe]_<minreads>read_<ntx>tx` |

### 3.4 Unchanged

`disome_proportion`, `start_codon_enrichment_ratio`,
`stop_codon_readthrough_ratio`, `five_prime_ramp_ratio`,
`three_prime_drop_ratio`, `prop_reads_CDS` / `_leader` / `_trailer`,
`ratio_cds:leader` / `ratio_cds:trailer` / `ratio_leader:trailer`,
`floss_median`, `complexity_distinct_positions`, `n_recommended_read_lengths`,
`periodicity_information_weighted_score`.

These are already honest raw quantities with no baked-in inversion. Most are
context-dependent diagnostics that `METRICS_DESIGN.md` §O3 says must never
carry a pass/fail badge, so they need no score and therefore no direction.

---

## 4. What it costs

**Key count.** 54 emitted metric keys today → **34 raw metrics + 16 scores**.
Twenty alias keys disappear; the raw/score split makes the remaining 34
unambiguous.

**JSON schema is a public contract.** `results["metrics"]` is consumed by
riboseq.org ingestion, the `evaluate --expected` YAMLs, any cohort table built
from `--summary-tsv` or `--comparison-csv`, and by TranslonScorer. A rename is
a breaking change for all of them. Two things follow:

1. **This is a major version.** v2.0.0, not a 1.4.x patch.
2. **Ship a compatibility shim for one minor cycle.** `results_output` gains a
   `LEGACY_KEY_ALIASES` map and emits both spellings under
   `results["metrics_legacy"]`, with a deprecation note in the JSON
   `provenance` block. Drop it at v2.1.

**Values change, not just names.** §3.1 changes eight stored numbers. Any
cohort table or threshold YAML carrying those keys must be regenerated —
crucially including the 6k-sample cohort QC outputs, where
`read_length_distribution_*_metric` columns would silently flip meaning if the
key were reused. That is precisely why §3.1 renames rather than reuses.

**Migration order (all done).**

1. Correctness fixes landed first on 1.4.x so the shipped scores were right
   under the old names.
2. `registry.py` + `tests/test_registry.py`: a table of
   `{key, unit, direction, scored_as}` with property tests asserting that every
   lower-is-better metric is scored through a flipping `method`, that no
   context-dependent metric carries a badge, and that no key is both a metric
   and a score.
3. The eight metric functions in §3.1 de-inverted and renamed.
4. Scores split into their own key namespace; scored records carry both the
   score key and the `metric` they came from, so the report shows a raw value
   in natural units beside a score in the good direction.
5. Aliases in §3.3 deleted, with `build_legacy_metrics` reproducing them.
6. `METRICS.md` generated by `scripts/generate_metrics_doc.py`; a test runs it
   with `--check` so the docs cannot drift again.

Steps 2 and 6 are the ones that stop this recurring. The rename alone would
not: the reason the `maxabs` inversion shipped is that nothing in the codebase
asserted which way was up.
