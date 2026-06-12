# RiboMetric Score Design Spec

**Status:** proposal for review; not yet implemented
**Scope:** defines what each metric should mean to a user and how that meaning is turned into a 0–1 score before code changes are made.

This is a design document, not a description of current behaviour. For current behaviour, see `METRICS.md` and `AUDIT_NOTES.md`.

> **Design note:** This document should remain a score-design specification, not an implementation PR. The main goal of this first iteration is to define interpretable directionality, anchors, provisional status thresholds, and metric roles. Empirical validation and recalibration come later once a reference corpus exists.

---

## 1. Why this exists

An audit of the current scoring found three problems:

1. **Different parts of the report make different decisions.**
   The summary cards, per-metric colour badges (`html_report._metric_status`), and QC gate (`results_output.evaluate_qc_status`) currently use separate threshold logic. As a result, the same metric on the same sample can appear acceptable in one place and unacceptable in another.
2. **Displayed percentages are not consistently calibrated.**
   Some metrics are raw fractions, some are transformed statistics, and some pass through a normalisation layer whose configured ranges are mostly `[0, 1]`. Therefore a displayed “70%” does not necessarily mean “70% of the way from bad to good”, and it does not mean the same thing across metrics.
3. **Some transforms are monotonic but not interpretable.**
   Several current scores use forms such as `1/(1+x)` for quantities such as CV, KL divergence, or Theil index. These preserve ordering, but they do not define a meaningful biological or analytical scale.

This spec fixes the meaning first. For each metric, it defines the user question, the raw quantity, the score mapping, and the decision consequence before any code is changed.

> **Design note:** The point is not that every current raw value is useless. Some raw values are meaningful. The problem is when raw values or arbitrary transforms are displayed as comparable percentages or used in shared pass/fail logic without a common interpretation.

---

## 2. The scoring contract

A RiboMetric score is not just a normalised number. It is a compact decision aid.

Every score obeys these rules:

1. **Directional and anchored.**
   The score defines what value of the underlying quantity maps to 0 and what maps to 1. Where possible, 0 is a theoretical or null-model value and 1 is the best meaningful value.
2. **Higher is always better.**
   Scores are mapped so that increasing score means increasing confidence, usability, or technical quality for the stated purpose.
3. **The raw value travels with the score.**
   The report shows both the score used for ranking/gating and the underlying quantity in natural units: frame dominance, CDS enrichment ratio, KL divergence in bits, duplicate rate, usable-read proportion, and so on.
4. **One score supports one primary decision.**
   Each scored metric names the consequence of a bad value: unreliable P-site offsets, RNA-like contamination, insufficient usable reads, hotspot-dominated coverage, distorted quantification, or reduced usable depth.
5. **Status thresholds are provisional until empirically calibrated.**
   In Phase 1, PASS/WARNING/FAIL cut-offs are set from theoretical expectations and presumed problematic values. In Phase 2, a reference corpus of good, borderline, and bad libraries will be used to recalibrate those thresholds.
6. **Make floors and baselines explicit; do not hide them in scaling.**
   Where a metric has a meaningful non-zero floor (e.g. periodicity's 1/3 random-frame baseline), present the raw, directly interpretable value and surface the floor as displayed context — a reference marker, an annotation, a one-line explanation. Do **not** rescale the floor to 0. Rescaling re-abstracts a number the user could otherwise read directly ("79% of A-sites in frame" becomes a meaningless "0.69"), and it buries the very baseline the user needs to judge the result. Reserve normalisation/rescaling for quantities that are genuinely not user-interpretable in their raw form. The job is to *aid understanding of the floor*, not to absorb it into the metric.

**Phase-1 anchoring policy:** use theoretical or null-model anchors wherever possible. Where no clean theoretical anchor exists, use a documented provisional anchor with clear directionality and mark it for empirical recalibration in Phase 2. Anchoring a score's interpretation to a floor (rule 6) means **showing** the floor, not subtracting it out.

A metric that does not change a general decision should not be forced into a percentage score. It should be shown as a diagnostic value, plot, or context label instead.

> **Design note:** For Phase 1, directionality is already valuable. We do not need to prove that every PASS/WARNING/FAIL boundary is empirically optimal before fixing the score semantics. The important distinction is:
>
> * score anchors should be interpretable now;
> * score direction should be interpretable now;
> * raw values should remain visible now;
> * status thresholds are reasonable provisional cut-offs until validated later.

---

## 3. Unified status model

All scored metrics are mapped onto the same 0–1, higher-is-better scale. A single config-driven status resolver then maps each score to a user-facing status.

Default Phase-1 thresholds:

|          Score | Status  |
| -------------: | ------- |
|      `>= 0.60` | PASS    |
| `0.30 – <0.60` | WARNING |
|       `< 0.30` | FAIL    |

These thresholds are provisional decision cut-offs. They are intended to reflect reasonable presumed problematic values under the Phase-1 theoretical scoring model. They are not claimed to be empirically optimal until calibrated against a reference corpus in Phase 2.

Per-metric overrides live in one config block. This replaces the current split between `DEFAULT_QC_THRESHOLDS`, `html_report._metric_status`, and `max_mins`.

The same status resolver is used by:

* summary cards;
* per-metric badges;
* report sections;
* QC gate logic.

However, **unified status does not mean universal gating**. Not every scored metric contributes to the overall QC gate. The gate uses an explicit subset of decision-critical metrics, especially Tier 1 metrics that answer “is this Ribo-seq-like?”. Other scored metrics may produce warnings or caveats without causing the whole sample to fail.

> **Design note:** The corrected principle is: one score scale, one status resolver, one config source of truth, explicit gate membership. The gate and the report should use the same status definitions for any metric they both display or evaluate, but the gate should not automatically include every scored metric.

A conceptual config structure might look like this:

```yaml
metrics:
  periodicity:
    score:
      method: frame_dominance_rescaled
      anchors:
        zero: 0.333333
        one: 1.0
    status:
      pass: 0.60
      warn: 0.30
    gate: true
  recommended_read_proportion:
    score:
      method: identity
    status:
      pass: 0.60
      warn: 0.30
    gate: false
  terminal_bias_kl_5prime:
    score:
      method: inverse_linear
      kl_max: 2.0
    status:
      pass: 0.70
      warn: 0.40
    gate: false
```

This is illustrative, not final syntax.

---

## 4. Report structure

The report is organised around the order in which a user should interpret a Ribo-seq library.

### Context strip — what am I looking at?

Shown before any pass/fail interpretation:

* inferred library type: elongation, initiation, disome-enriched, low-quality, or unknown;
* total reads and mapped reads;
* organism / annotation source;
* dominant read-length mode;
* analysis mode, where relevant: genomic BAM, transcriptomic BAM, monosome, disome, initiation-enriched, etc.

This context determines how later metrics should be interpreted. For example, high disome signal may be a warning in a monosome library but expected in a disome experiment.

### Tier 1 — Is this Ribo-seq-like?

These are the decision-critical metrics for basic Ribo-seq identity and frame-dependent interpretation:

* periodicity / frame dominance;
* periodicity information;
* CDS enrichment.

These metrics form the default global QC gate. If they fail, downstream frame-dependent analyses such as P-site assignment, ORF calling, or codon-level interpretation should not proceed without an explicit warning.

### Tier 2 — Is it usable for my analysis?

These metrics describe whether the library has enough usable signal for common analyses:

* recommended-read proportion;
* coverage breadth / concentration;
* saturation or marginal position discovery rate.

A library can pass Tier 1 but still be weak for a particular analysis because too little of the library survives filtering, coverage is hotspot-dominated, or additional sequencing would still discover many new positions.

### Tier 3 — What caveats should I carry forward?

These metrics describe technical distortions and loss of usable depth:

* duplicate rate;
* multimapping rate;
* soft-clipping;
* terminal nucleotide / ligation bias.

These usually produce warnings or interpretation caveats rather than causing the whole library to fail. Their thresholds should be protocol-aware where possible.

> **Design note:** The report should distinguish:
>
> * failed Ribo-seq identity;
> * passed but weak for a specific analysis;
> * passed with caveats.
>
> This avoids treating every warning as the same kind of failure.

---

## 5. Per-metric score specs

Template:

* **Q:** user question;
* **Raw:** underlying quantity and units;
* **Anchors:** values mapped to 0 and 1, with at least one interior reference where possible;
* **Score:** score mapping;
* **Decision:** what changes when the score is bad;
* **Design note:** reviewer / implementation caution.

---

# Tier 1 — Is this Ribo-seq-like?

Tier 1 metrics answer whether the library has the basic properties expected of Ribo-seq and whether frame-dependent interpretation is defensible. These metrics form the default global QC gate. Failure here means that P-site assignment, ORF calling, codon-level analysis, and other frame-dependent interpretations should not proceed without an explicit warning.

---

## `periodicity` / frame dominance

* **Q:** Can I trust reading-frame structure and P-site assignment?
* **Raw:** dominant-frame fraction of in-frame CDS reads. The meaningful range is 1/3–1.0.
* **Score:** the raw fraction itself (identity map). The value is directly interpretable — "0.79 means 79% of A-sites are in the dominant frame" — and that interpretability is the whole point. **Do not rescale against the 1/3 floor.**
* **Floor / baseline (shown, not hidden):** the random-frame expectation is 1/3 (≈0.33). Surface it in the report as a reference marker so the user reads the value against it: 0.33 = random, ~0.70+ = strong periodicity, 1.0 = perfect. Where this floor is not obvious to a given user, the report should *teach* it (a caption / baseline line on the plot), not bury it inside the number.
* **Status (thresholds set on the raw fraction, with the floor in mind):** PASS at dominance ≥ ~0.70, WARNING ≥ ~0.50, FAIL below. Provisional, not yet empirically calibrated.
* **Decision:** low value means weak or absent triplet structure. Frame-dependent interpretation, including P-site assignment and ORF calling, is unreliable unless rescued by length-specific evidence.

> **Design note (correction, supersedes earlier draft):** An earlier version of this spec scored frame dominance by rescaling against the 1/3 floor: `clip((d − 1/3)/(2/3), 0, 1)`. **That is rejected for interpretation.** Rescaling turns the most interpretable quantity in the report (a literal in-frame percentage) back into an abstract one (0.79 → 0.69) and hides the random baseline the user needs. Per contract rule 6, the floor is *communicated* (displayed as a baseline), and the score stays on the interpretable raw axis. This is the general pattern wherever a metric has a known floor — make it explicit, don't normalise it away.
>
> Good periodicity does not by itself prove that offsets are correct; offset estimation still needs to be length-aware.
>
> **Implementation status:** Phase 1A/1B code currently scores this with `frame_dominance_rescaled` (the rejected approach). That is now a known divergence from this plan — when implementation resumes, switch `periodicity_dominance` to `method: identity`, move the status thresholds onto the raw fraction (pass ~0.70 / warn ~0.50), and add the 1/3 baseline marker to the report.

---

## `periodicity_information`

* **Q:** Is the frame distribution non-random, regardless of which frame dominates?
* **Raw:** entropy reduction of the frame distribution.
* **Anchors:**
  * 0 = random frame distribution, where `H = log2(3)`;
  * 1 = perfect single-frame occupancy, where `H = 0`.
* **Score:**
  `(log2(3) - H) / log2(3)`
* **Decision:** use as a cross-check against frame dominance. Large disagreement between dominance and information can indicate frame mixing, unstable offsets, or read-length classes with different frame behaviour.
* **Implementation note:** drop the current `sqrt` transform because it inflates middling information values without improving interpretability.

> **Design note:** This should not automatically be treated as a fully independent second gate duplicating frame dominance. It is better framed as an information-content cross-check. If the gate includes both dominance and information, their relationship should be made explicit.

---

## CDS enrichment

* **Q:** Are reads concentrated in annotated coding sequence more than expected under an RNA-seq-like transcript coverage model?
* **Raw:** enrichment ratio:
  `E = observed_CDS_body_fraction / expected_CDS_fraction`
  Also show the raw `observed_CDS_body_fraction`.
* **Expected fraction (resolved):** expression-weighted, boundary-aligned to the
  observed `"CDS"` bucket:
  `expected_CDS_fraction = Σ_t reads_t·(CDS_body_len_t / eligible_mRNA_len_t) / Σ_t reads_t`
  where `CDS_body_len_t = max(0, cds_end_t − cds_start_t − 1)` (interior CDS,
  excluding the single start_codon/stop_codon positions, matching
  `assign_mRNA_category`), `eligible_mRNA_len_t = transcript_length_t` (the full
  transcript — every A-site falls in one of the five region buckets, so the eligible
  length matching the observed denominator is the whole transcript), and `reads_t` is
  the total reads on transcript t. Include t only if `transcript_length_t > 0` and
  `CDS_body_len_t > 0` and t appears in reads. Self-consistency: uniform-per-nucleotide
  reads ⇒ observed ≈ expected ⇒ E ≈ 1 ⇒ score ≈ 0.
* **Anchors:**
  * 0 = `E <= 1`, meaning no enrichment over the RNA-like expectation;
  * 1 = strong CDS enrichment;
  * interior reference: `E = 2` → score 0.5.
* **Score:**
  `clip(1 - 1/E, 0, 1)`
* **Decision:** low score suggests the library is not strongly enriched for ribosome-protected coding signal and may reflect degradation, RNA contamination, poor nuclease protection, or annotation/assignment mismatch.
* **Caveat:** this is a model-based score. Interpretation depends on annotation quality, transcript selection, expressed-transcript definition, and whether the experiment is expected to capture substantial translation outside annotated CDS.

> **Design note:** Replacing raw `prop_reads_CDS` with an enrichment ratio is the right move because raw CDS fraction is organism- and annotation-confounded. But this metric introduces dependency on annotation and transcript selection. It should be Tier 1, but with explicit caveats and an `annotation_required` or equivalent flag in implementation.

---

# Tier 2 — Is it usable for my analysis?

Tier 2 metrics answer whether a Ribo-seq-like library is useful for common downstream analyses. These metrics should generally not decide whether the sample is “Ribo-seq” at all. Instead, they describe whether enough usable, sufficiently distributed signal remains for the user’s intended analysis.

---

## `recommended_read_proportion`

* **Q:** Which read lengths do I keep, and how much of my library survives?
* **Raw:** fraction of reads in lengths whose dominant-frame fraction is at least `min_periodicity` and whose read-length share is at least `min_read_proportion`.
* **Anchors:**
  * 0 = none of the library is in recommended read lengths;
  * 1 = all of the library is in recommended read lengths;
  * interior reference: 0.5 = half the library survives the recommended-read filter.
* **Score:**
  `recommended_read_proportion`
* **Decision:** low score means that even if the library has some good periodic signal, little of the dataset is usable for frame-sensitive work after filtering.
* **Status note:** this should usually be an analysis-suitability warning rather than a global Ribo-seq identity failure.

> **Design note:** This is better than scoring IQR, CV, max-proportion, bimodality, or normality as separate proxies for read-length shape. It links the score to a real consequence: how much usable library remains. The report can still show read-length distribution diagnostics, but this should be the headline read-length score.

---

## Coverage breadth / concentration

Candidate metric: `uniformity_entropy`, optionally paired with Gini or another concentration statistic.

* **Q:** Can I quantify broadly, or is signal dominated by a few hotspots?
* **Raw:** normalised metagene-window entropy or another bounded concentration statistic.
* **Key subtlety:** perfect uniformity is not necessarily the biological ideal. Real translation has pauses, local rate variation, and non-uniform coverage. The useful warning is primarily at the low end, where coverage is dominated by a small number of positions.
* **Anchors, Phase 1 provisional:**
  * 0 = all signal concentrated in one bin;
  * 1 = maximally even distribution across bins.
* **Score:**
  `uniformity_entropy` as-is for Phase 1, interpreted as a concentration detector.
* **Decision:** low score means the library may not support broad position-level or gene-level quantification because a small number of hotspots carry too much of the signal.
* **Phase-2 note:** replace the theoretical top anchor with empirical good-library bands. The goal is not mathematical flatness; it is avoiding pathological concentration.

> **Design note:** This is one of the weakest places for pure theory. Directionality is useful now, but the “good” range really needs a reference corpus. Do not imply that perfectly uniform coverage is biologically ideal. The score should mainly flag pathological concentration.

---

## Saturation / complexity

Candidate metric: `marginal_position_discovery_rate`.

* **Q:** Did I sequence deep enough?
* **Raw:** marginal rate at which additional reads discover new covered positions. Lower marginal discovery means greater saturation.
* **Anchors:**
  * 0 = every new read discovers a new position, indicating severe under-sampling;
  * 1 = additional reads discover no new positions, indicating complete saturation.
* **Score:**
  `1 - marginal_position_discovery_rate`
* **Decision:** low score means additional sequencing is likely to discover substantially more signal; high score means the library is approaching saturation for covered positions.
* **Status note:** this is primarily an experimental-design and budget metric, not a Ribo-seq identity metric.

> **Design note:** This is useful, but be careful with interpretation. A saturated bad library can still be bad; saturation does not mean biological quality. It answers “would more sequencing add breadth?” not “is the library good?”.

---

# Tier 3 — Technical caveats

Tier 3 metrics describe technical distortions, loss of usable depth, and protocol-specific caveats. They should usually colour the report and inform interpretation, but not automatically fail the whole sample unless explicitly configured to do so.

---

## Mapping hygiene

Metrics:

* `duplicate_rate`
* `rpf_multimapper_rate`
* `soft_clip_rate_5prime`

### `duplicate_rate`

* **Q:** How much of my apparent depth may be dominated by duplicated molecules?
* **Raw:** fraction of reads marked or inferred as duplicates.
* **Anchors:**
  * 0 = duplicate rate 1.0, all reads duplicated;
  * 1 = duplicate rate 0.0, no duplicates detected.
* **Score:**
  `1 - duplicate_rate`
* **Decision:** low score means usable molecule diversity may be much lower than read depth suggests.
* **Protocol caveat:** high duplication may be expected in low-input Ribo-seq or highly enriched libraries.

> **Design note:** This should not be a naive universal gate. Duplication is protocol- and input-dependent. Defaults should be generous, and the raw duplicate rate must be shown beside the score.

### `rpf_multimapper_rate`

* **Q:** How many reads cannot be assigned uniquely under the chosen mapping model?
* **Raw:** fraction of reads that multimapped.
* **Anchors:**
  * 0 = multimapper rate 1.0;
  * 1 = multimapper rate 0.0.
* **Score:**
  `1 - rpf_multimapper_rate`
* **Decision:** low score means reduced confidence in locus-level or transcript-level quantification.
* **Mapping caveat:** on a transcriptome BAM, multimapping often reflects isoform multiplicity rather than ambiguity between genomic loci.

> **Design note:** Interpretation depends strongly on whether the BAM is genomic or transcriptomic. The context strip should say this up front. A high transcriptome multimapper rate is not the same biological problem as a high genomic multimapper rate.

### `soft_clip_rate_5prime`

* **Q:** Are read ends being clipped in a way that may affect P-site assignment or terminal-bias interpretation?
* **Raw:** fraction of reads with 5′ soft clipping, or average 5′ soft-clipped bases depending on current implementation.
* **Anchors:**
  * 0 = severe or universal 5′ soft clipping;
  * 1 = no 5′ soft clipping.
* **Score:**
  If the raw value is a rate: `1 - soft_clip_rate_5prime`.
  If the raw value is a mean clipped length, define a documented `soft_clip_max` and use:
  `clip(1 - mean_soft_clip / soft_clip_max, 0, 1)`
* **Decision:** low score means read-end positions may be unreliable for offset inference, terminal-bias assessment, or fine-grained positional analysis.

> **Design note:** Confirm the raw definition before implementation. A rate and a mean clipped length need different anchors. Do not silently treat them as equivalent.

---

## Terminal nucleotide / ligation bias

Metrics:

* `terminal_bias_kl_*`
* `terminal_nucleotide_bias_max_absolute_*`

### `terminal_bias_kl_*`

* **Q:** Is sequence bias from ligation, nuclease preference, or library preparation distorting counts?
* **Raw:** KL divergence in bits between observed terminal dinucleotide frequencies and a defined background composition.
* **Anchors, Phase 1 provisional:**
  * 1 = 0 bits KL, no detectable divergence from background;
  * 0 = `KL_max`, a documented concern threshold. **Resolved provisional value:
    `KL_max = 2.0` bits** — KL = 1 bit → score 0.5 (preserves the midpoint of the old
    `1/(1+KL)` transform so verdicts don't lurch); KL ≥ 2 bits → 0 (one/few
    dinucleotides dominate). To be tightened (likely toward ~1.0) in Phase 2.
* **Score:**
  `clip(1 - KL / KL_max, 0, 1)`
* **Decision:** low score means terminal sequence bias may distort count-level quantification and should be considered for correction or caveated interpretation.
* **Phase-2 note:** set `KL_max` and warning thresholds from empirical distributions of clean and biased libraries.

> **Design note:** Replacing `1/(1+KL)` is correct because that transform is monotonic but not interpretable. However, the Phase-1 `KL_max` will be provisional. The report should show raw KL in bits prominently because that is the expert-interpretable value.

### `terminal_nucleotide_bias_max_absolute_*`

* **Q:** What is the largest absolute terminal dinucleotide deviation from expectation?
* **Raw:** maximum absolute difference between observed and expected terminal dinucleotide probability.
* **Anchors:**
  * 0 = maximum possible or configured unacceptable deviation;
  * 1 = no deviation from expectation.
* **Score:**
  If raw is already bounded 0–1 as max absolute difference:
  `1 - maxdiff`
* **Decision:** low score means at least one terminal dinucleotide is strongly over- or under-represented and may distort positional or count-level interpretation.

> **Design note:** This is more immediately interpretable than KL for some users because it points to the worst offending terminal composition difference. It should complement, not necessarily replace, KL.

---

# 6. Diagnostics: no score, no percent

These values describe shape, context, or experiment-specific behaviour. They should be shown as plots, labels, or raw values, but they should not be forced into a general 0–1 score unless a specific decision rule is later defined.

---

## Read-length distribution diagnostics

Metrics:

* `read_length_distribution_IQR_metric`
* `read_length_distribution_CV_metric`
* `read_length_distribution_max_prop_metric`
* `read_length_distribution_bimodality`
* `read_length_distribution_normality`

These should not be primary scored metrics in the default report.

Read-length distributions in good Ribo-seq libraries are not necessarily normal. They may be multi-modal for legitimate biological, technical, or protocol-specific reasons. Normality in particular can become uninformative at realistic read counts.

The report should instead show:

* read-length distribution plot;
* dominant read-length mode;
* recommended read lengths;
* proportion of reads retained under recommended-read filtering;
* optional warning for clearly pathological distributions.

> **Design note:** The main scored read-length metric should be `recommended_read_proportion`, because it connects shape to consequence. IQR, CV, max-proportion, bimodality, and normality are competing proxies unless tied to a specific downstream decision.

---

## `disome_proportion`

* **Role:** diagnostic / context label.
* **Interpretation:** directional by experiment type.
  * In a monosome library, high disome proportion may indicate contamination or incomplete size selection.
  * In a disome-focused experiment, high disome proportion may be expected or desirable.
* **Default status:** no universal score, no universal pass/fail.
* **Presentation (resolved):** render the raw value with a neutral caption keyed off
  `classify_library_type`, never a pass/fail badge — monosome library → "high values
  may indicate disome contamination or incomplete size selection"; disome experiment →
  "expected/desirable"; unknown → present both readings. This is the general rule for
  any diagnostic whose "good direction" depends on experiment type; the context strip
  declares the library type up front so the caption has something to key on.

> **Design note:** This is a clear example of why not every value should become a score. The same value can mean contamination or signal depending on the intended library type.

---

## `classify_library_type`

* **Role:** context label.
* **Output:** elongation, initiation, disome-enriched, low-quality, unknown, or other supported classes.
* **Default status:** no score.

This should appear in the context strip and modify interpretation of later diagnostics.

> **Design note:** A classifier label is not a QC score. It informs how scores and diagnostics should be read.

---

## Start / stop / codon-level diagnostics

Metrics:

* `start_codon_enrichment_ratio`
* `stop_codon_readthrough_ratio`
* codon dwell-time statistics;
* RUST-style diagnostics;
* other experiment-specific metagene or codon effects.

These should be shown as raw values, plots, or optional analysis modules. They should not be included in the default global QC score unless the user explicitly requests an analysis-specific gate.

> **Design note:** These values are biologically interesting, but their interpretation is experiment-specific. For example, start enrichment may be expected in initiation profiling but not in elongation profiling. Stop/readthrough metrics can indicate biology, artefact, or annotation mismatch. They should not silently affect generic sample pass/fail.

---

# 7. Phase 2 — empirical calibration

Phase 1 makes the scores interpretable; it does not claim that every PASS/WARNING/FAIL boundary is empirically optimal. Phase 2 uses a reference corpus to test whether provisional decision thresholds match real library quality distributions.

Phase 2 replaces provisional anchors or thresholds with empirical bands derived from a reference corpus of known-good, borderline, and bad libraries.

Priority targets for empirical calibration:

1. coverage breadth / concentration;
2. terminal nucleotide bias KL thresholds;
3. duplication tolerances;
4. multimapping tolerances;
5. recommended-read proportion thresholds;
6. CDS enrichment expectations across annotations and organisms.

Deliverables:

1. A curated reference corpus labelled at least as good, borderline, and bad.
2. Per-metric empirical distributions of raw quantities.
3. Per-metric empirical anchor or threshold table, for example p05 / p50 / p95.
4. Recalibrated score and status thresholds behind the same scoring API.
5. A changelog documenting which Phase-1 assumptions were retained, adjusted, or rejected.

> **Design note:** The reference corpus should not only contain “good” libraries. It needs deliberately bad, borderline, contaminated, weakly periodic, low-depth, biased, and protocol-diverse examples. Otherwise Phase 2 will calibrate only the centre of the good-library distribution and fail to define useful warning/fail zones.

Possible corpus sources include the RiboSeq.Org / Trips-Viz ecosystem, public benchmarking datasets, manually reviewed internal examples, and protocol-specific known-good libraries.

---

# 8. What changes in code

This section is high-level only. Do not implement until the score design is accepted.

Required architectural changes:

1. **Create one scoring module.**
   It should map:
   `raw_quantity -> score -> status`
   using anchors, transforms, thresholds, and gate membership from one config-driven source of truth.
2. **Retire split threshold logic.**
   Replace the current split between summary card thresholds, badge thresholds, QC gate thresholds, and `max_mins`.
3. **Retire arbitrary transforms.**
   Replace `1/(1+x)` transforms and inert normalisation with explicit anchored maps.
4. **Score frame dominance as the raw in-frame fraction (identity), and surface the 1/3 floor as displayed context.**
   Do NOT rescale the floor into the number. Set status thresholds on the raw fraction (pass ~0.70 / warn ~0.50) and add a baseline marker to the report so the floor is taught, not hidden. (Supersedes the earlier "rescale against 1/3" instruction; see the periodicity design note.)
5. **Drop the square-root transform in periodicity information.**
   Entropy reduction is already naturally anchored.
6. **Add CDS enrichment computation.**
   Compute expected CDS fraction from annotation and expressed transcript models. Show both observed CDS fraction and enrichment ratio.
7. **Reclassify diagnostics.**
   Remove default percent scores and pass/fail labels for read-length normality, bimodality, disome proportion, library classification, start/stop enrichment, codon dwell-time, and RUST-style outputs unless the user explicitly enables an analysis-specific gate.
8. **Add explicit gate membership.**
   Every scored metric can have a status. Only selected metrics contribute to global QC.
9. **Make protocol/context visible.**
   Report context should include library type, mapping mode, annotation source, and analysis mode because they change interpretation.

> **Design note:** Codex or another implementation agent should not infer that all metrics need to be implemented at once. The first safe implementation target is the shared scoring/status resolver and a small number of high-confidence metrics: frame dominance, periodicity information, recommended-read proportion, and simple mapping hygiene rates. CDS enrichment requires additional annotation-dependent design.

---

# 9. Suggested initial implementation order

This is not part of the biological design, but it may help avoid a messy implementation.

## Phase 1A — scoring infrastructure

Implement:

* one scoring/status resolver;
* one metric config block;
* explicit `gate: true/false`;
* raw value + score + status data structure;
* report components reading from the same status resolver.

Do not change metric biology yet except where required to route through the new scoring layer.

## Phase 1B — clean theoretical scores

Implement or update:

* `periodicity` scored as the raw in-frame fraction (identity), with the 1/3 floor surfaced in the report — not rescaled;
* `periodicity_information` without square-root inflation;
* `recommended_read_proportion` as the headline read-length score;
* `duplicate_rate`, `rpf_multimapper_rate`, and `soft_clip_rate_5prime` as `1 - rate` where the raw quantity is confirmed to be a rate.

## Phase 1C — model-dependent scores

Implement:

* CDS enrichment;
* expected CDS fraction calculation;
* annotation/context warnings.

## Phase 1D — provisional technical caveats

Implement:

* terminal KL linear score with documented `KL_max`;
* max absolute terminal nucleotide deviation;
* coverage concentration score.

## Phase 1E — report cleanup

Reclassify as diagnostics:

* read-length normality;
* read-length bimodality;
* disome proportion;
* library type classifier;
* start/stop/codon-level metrics unless explicitly requested.

> **Design note:** This order avoids letting implementation complexity drive the biological design. The score model can be fixed before all metrics are perfect.

---

# 10. Codex handoff prompt

Use this prompt when passing the document to Codex or another coding agent:

```text
Update the RiboMetric scoring design/spec according to the markdown below.
Do not implement code yet unless explicitly instructed. Preserve the distinction
between raw quantities, anchored scores, provisional status thresholds,
diagnostics, and gate membership.
The key design principles are:
1. all scored metrics are 0-1 and higher-is-better;
2. raw values must remain visible beside scores;
3. scores need explicit anchors and decision consequences;
4. Phase-1 thresholds are provisional presumed cut-offs, not empirically
   validated boundaries;
5. unified status does not mean universal gating;
6. the global QC gate uses explicit gate membership;
7. diagnostics should not be forced into percent scores or pass/fail labels.
After updating the spec, propose an implementation plan, but do not modify code
until asked.
```

---

# 11. Summary of key design corrections

The first-iteration corrected design makes these changes:

1. It separates score anchors from status thresholds.
2. It treats Phase-1 status thresholds as provisional but useful.
3. It makes all scored metrics higher-is-better.
4. It preserves raw values beside scores.
5. It introduces explicit gate membership.
6. It distinguishes failed Ribo-seq identity from weak analysis suitability and technical caveats.
7. It keeps frame dominance as a raw, directly interpretable in-frame fraction with the 1/3 floor surfaced as context (earlier floor-rescaling rejected: floors are taught, not hidden in scaling).
8. It treats periodicity information as a cross-check unless explicitly gated.
9. It replaces raw CDS proportion with model-based CDS enrichment.
10. It demotes read-length shape proxies to diagnostics.
11. It treats disome proportion and library type as context-dependent diagnostics.
12. It defers empirical calibration to a reference corpus without blocking Phase-1 score cleanup.
