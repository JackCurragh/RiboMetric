# Metric Naming and Direction

**Status:** proposal for review. Nothing in this document is implemented yet —
it defines the target naming scheme so the migration can be costed before it is
started. For current behaviour see [`METRICS.md`](METRICS.md); for how raw
values become scores see [`METRICS_DESIGN.md`](METRICS_DESIGN.md).

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
| `read_length_distribution_IQR_metric` | `1 − IQR/range` | `read_length_iqr_fraction` | `IQR/(P90−P10)` | `read_length_concentration_score` | `one_minus_rate` |
| `read_length_distribution_coefficient_of_variation_metric` | `1/(1+CV)` | `read_length_cv` | `CV` | — (diagnostic) | — |
| `read_length_distribution_bimodality_metric` | `1/(1+BC)` | `read_length_bimodality_coefficient` | `BC` | — (diagnostic) | — |
| `read_length_distribution_normality_metric` | `1 − p` | `read_length_normality_pvalue` | `p` | — (diagnostic) | — |

Note that `read_length_cv`, `read_length_bimodality_coefficient` and
`read_length_normality_pvalue` lose their scores entirely. Per
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

**Migration order.**

1. Land the correctness fixes on 1.4.x first (this branch) so the shipped
   scores are right under the current names.
2. Add the registry + invariant test on the current names — a table of
   `{key, direction, unit, scored_as}` with a property test asserting
   monotonicity. This makes the rename mechanical rather than judgement-based.
3. De-invert the eight metric functions in §3.1, renaming as you go, one family
   per commit (terminal bias, then read length).
4. Split `results["scores"]` out of `results["metrics"]`.
5. Delete the aliases in §3.3.
6. Rewrite `METRICS.md` from the registry so the docs cannot drift again
   (`docs/METRICS.md` currently documents four default metric names that are
   not emitted and omits 32 that are).

Steps 2 and 6 are the ones that stop this recurring. The rename alone does not:
the reason the `maxabs` inversion shipped is that nothing in the codebase
asserts which way is up.
