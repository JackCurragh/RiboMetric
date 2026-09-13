# RiboMetric Metrics

<!-- GENERATED FILE - do not edit by hand.
     Regenerate with: python scripts/generate_metrics_doc.py -->

Two kinds of key appear in RiboMetric output, and they follow different rules (see [METRIC_NAMING.md](METRIC_NAMING.md)):

* **Metrics** (`results["metrics"]`) name the quantity that was measured, in its natural units and its natural direction. A rate goes up as the library gets worse.
* **Scores** (in the report, the QC status file and the summary plot) name the good property, always lie in [0, 1], and are always higher-is-better.

## Scored metrics

Every score below is higher-is-better. The direction flip, where one is needed, lives in the `method` column and nowhere else.

### Tier 1 — is this Ribo-seq?

Gated. A failure here means frame-dependent analysis should not proceed.

| Score | From metric | Unit | Metric direction | Method | Pass / warn |
|---|---|---|---|---|---|
| `periodicity_dominance_score` | `periodicity_dominance` | fraction | higher is better | `identity` | 0.70 / 0.50 |
| `cds_enrichment_score` | `cds_enrichment_ratio` | ratio | higher is better | `enrichment_ratio` | 0.60 / 0.30 |
| `periodicity_information_score` | `periodicity_information` | fraction | higher is better | `identity` | 0.25 / 0.05 |

- **`periodicity_dominance_score`** — Fraction of coding A-sites in the dominant reading frame; the global value uses one shared dominant frame. Low score: weak triplet structure; P-site assignment and ORF calling unreliable.
- **`cds_enrichment_score`** — Observed CDS-body read fraction over the length-weighted expected fraction (E). Low score: reads not enriched over coding sequence; library may reflect degradation, RNA contamination, or poor nuclease protection.
- **`periodicity_information_score`** — Entropy reduction of the frame distribution against a uniform three-frame null. Cross-check on frame dominance; large disagreement signals frame mixing or unstable offsets.

### Tier 2 — is it usable for my analysis?

Not gated. These describe whether enough usable signal survives filtering.

| Score | From metric | Unit | Metric direction | Method | Pass / warn |
|---|---|---|---|---|---|
| `usable_read_fraction_score` | `recommended_read_proportion` | fraction | higher is better | `identity` | 0.60 / 0.30 |
| `coverage_uniformity_score` | `uniformity_entropy` | fraction | higher is better | `identity` | 0.60 / 0.30 |
| `library_saturation_score` | `marginal_position_discovery_rate` | rate | lower is better | `one_minus_rate` | 0.60 / 0.30 |

- **`usable_read_fraction_score`** — Fraction of the library carried by read lengths recommended for frame-sensitive work. Low score: little of the library survives recommended-read filtering for frame-sensitive work.
- **`coverage_uniformity_score`** — Normalised entropy of the codon-binned start-codon metagene. Low score: coverage dominated by a few hotspots; broad quantification may be unreliable.
- **`library_saturation_score`** — Fraction of reads at the margin of sequencing depth landing on a position not already seen; high means under-sequenced. Low score: library under-sequenced; more reads would discover substantially more positions.

### Tier 3 — technical caveats

Not gated. These colour interpretation but do not fail a sample.

| Score | From metric | Unit | Metric direction | Method | Pass / warn |
|---|---|---|---|---|---|
| `fragment_uniqueness_score` | `duplicate_rate` | rate | lower is better | `one_minus_rate` | 0.60 / 0.30 |
| `rpf_unique_mapping_score` | `rpf_multimapper_rate` | rate | lower is better | `one_minus_rate` | 0.60 / 0.30 |
| `alignment_unique_mapping_score` | `alignment_multimapper_rate` | rate | lower is better | `one_minus_rate` | 0.60 / 0.30 |
| `terminal_integrity_5prime_score` | `soft_clip_rate_5prime` | rate | lower is better | `one_minus_rate` | 0.60 / 0.30 |
| `footprint_homogeneity_score` | `floss_aberrant_transcript_fraction` | fraction | lower is better | `one_minus_rate` | 0.60 / 0.30 |
| `terminal_evenness_kl_5prime_score` | `terminal_bias_kl_5prime` | bits | lower is better | `inverse_linear` | 0.70 / 0.40 |
| `terminal_evenness_kl_3prime_score` | `terminal_bias_kl_3prime` | bits | lower is better | `inverse_linear` | 0.70 / 0.40 |
| `terminal_evenness_maxdev_5prime_score` | `terminal_bias_max_deviation_5prime` | fraction | lower is better | `one_minus_rate` | 0.70 / 0.40 |
| `terminal_evenness_maxdev_3prime_score` | `terminal_bias_max_deviation_3prime` | fraction | lower is better | `one_minus_rate` | 0.70 / 0.40 |

- **`fragment_uniqueness_score`** — Fraction of reads that are collapsed duplicates. Low score: usable molecule diversity much lower than read depth suggests (protocol-dependent).
- **`rpf_unique_mapping_score`** — Fraction of weighted fragments reported at more than one alignment location. Low score: reduced confidence in locus/transcript-level quantification.
- **`alignment_unique_mapping_score`** — Fraction of alignment rows whose fragment has evidence of another reported alignment. Low score: many alignment rows have evidence of another reported alignment.
- **`terminal_integrity_5prime_score`** — Fraction of reads with 5' soft-clipping. Low score: 5' read ends frequently clipped; offset and terminal-bias interpretation may be unreliable.
- **`footprint_homogeneity_score`** — Fraction of transcripts whose footprint-length profile departs from the library aggregate beyond the FLOSS cutoff. Low score: many transcripts have footprint-length profiles unlike the library aggregate; heterogeneous or contaminated library.
- **`terminal_evenness_kl_5prime_score`** — Kullback-Leibler divergence of observed 5' terminal dinucleotide frequencies from the background. Low score: 5' terminal sequence bias may distort count quantification; consider correction.
- **`terminal_evenness_kl_3prime_score`** — Kullback-Leibler divergence of observed 3' terminal dinucleotide frequencies from the background. Low score: 3' terminal sequence bias may distort count quantification; consider correction.
- **`terminal_evenness_maxdev_5prime_score`** — Largest absolute deviation of a 5' terminal dinucleotide frequency from its background frequency. Low score: at least one 5' terminal dinucleotide is strongly over- or under-represented.
- **`terminal_evenness_maxdev_3prime_score`** — Largest absolute deviation of a 3' terminal dinucleotide frequency from its background frequency. Low score: at least one 3' terminal dinucleotide is strongly over- or under-represented.

## Diagnostics

Reported as raw measurements with no pass/fail badge. Either their good direction depends on the protocol, or they describe shape rather than quality.

| Metric | Unit | Direction | Summary |
|---|---|---|---|
| `cds_coverage` | fraction | higher is better | Proportion of CDS positions covered, using the configured in-frame and minimum-read settings. |
| `cds_coverage_100read_100tx` | fraction | higher is better | CDS coverage, any frame, >=100 reads, 100 transcripts. |
| `cds_coverage_1read_1000tx` | fraction | higher is better | CDS coverage, any frame, >=1 read, 1000 transcripts. |
| `cds_coverage_inframe_100read_100tx` | fraction | higher is better | In-frame CDS coverage, >=100 reads, 100 transcripts. |
| `cds_coverage_inframe_1read_1000tx` | fraction | higher is better | In-frame CDS coverage, >=1 read, 1000 transcripts. |
| `cga_dwell` | ratio | context-dependent | Relative A-site dwell signal on the CGA codon. |
| `codon_dwell_cv` | index | context-dependent | Coefficient of variation of A-site codon dwell-times. |
| `codon_dwell_p90_p10` | ratio | context-dependent | Ratio of the 90th to 10th percentile codon dwell-time. |
| `complexity_distinct_positions` | count | context-dependent | Distinct A-site positions observed at full depth. |
| `disome_proportion` | fraction | context-dependent | Fraction of reads in the di-some read-length window; expected in a di-some experiment, contamination in a monosome one. |
| `five_prime_ramp_ratio` | ratio | context-dependent | A-site density in the 5' portion of the CDS over the body. |
| `floss_median` | index | lower is better | Median per-transcript FLOSS score. |
| `n_recommended_read_lengths` | count | context-dependent | Number of read lengths recommended for frame-sensitive work. |
| `periodicity_autocorrelation` | index | higher is better | Optional: triplet periodicity from signal autocorrelation. |
| `periodicity_fourier` | index | higher is better | Optional: Fourier power at the codon frequency. |
| `periodicity_information_weighted_score` | fraction | higher is better | Read-depth weighted periodicity information content. |
| `periodicity_trips-viz` | index | higher is better | Optional: Trips-Viz style triplet periodicity score. |
| `proline_dwell` | ratio | context-dependent | Relative A-site dwell signal on proline codons. |
| `prop_reads_CDS` | fraction | context-dependent | Proportion of A-sites falling in the CDS body. |
| `prop_reads_leader` | fraction | context-dependent | Proportion of A-sites falling in the 5' leader. |
| `prop_reads_trailer` | fraction | context-dependent | Proportion of A-sites falling in the 3' trailer. |
| `ratio_cds:leader` | ratio | higher is better | CDS reads relative to 5' leader reads. |
| `ratio_cds:trailer` | ratio | higher is better | CDS reads relative to 3' trailer reads. |
| `ratio_leader:trailer` | ratio | context-dependent | 5' leader reads relative to 3' trailer reads. |
| `read_length_bimodality_coefficient` | index | lower is better | Sarle's bimodality coefficient of the read length distribution. |
| `read_length_cv` | index | lower is better | Coefficient of variation of the read length distribution. |
| `read_length_iqr_fraction` | fraction | lower is better | Interquartile range of the read length distribution as a fraction of its 10th-90th percentile range. |
| `read_length_max_proportion` | fraction | context-dependent | Proportion of reads at the most frequent read length. |
| `read_length_normality_pvalue` | pvalue | context-dependent | Normaltest p-value for the read length distribution; a footprint distribution is not expected to be normal. |
| `rust_mean_kl_divergence` | bits | context-dependent | Mean RUST codon-metagene KL divergence. |
| `start_codon_enrichment_ratio` | ratio | context-dependent | Reads near the start codon relative to the CDS body; very high values indicate initiation-stalling treatments. |
| `stop_codon_readthrough_ratio` | ratio | lower is better | Reads downstream of the stop codon relative to upstream. |
| `three_prime_drop_ratio` | ratio | context-dependent | A-site density in the 3' portion of the CDS over the body. |
| `uniformity_autocorrelation` | index | higher is better | Optional: autocorrelation-based coverage smoothness. |
| `uniformity_gini_index` | index | lower is better | Optional: Gini coefficient of coding-region coverage. |
| `uniformity_theil_index` | index | lower is better | Optional: Theil inequality index of coding-region coverage. |

## Renamed and removed keys

Pre-2.0 spellings are reproduced under `results["metrics_legacy"]` for one minor cycle and removed at v2.1. Where the transform is not `identity`, the old key held a *different number*: it had a goodness transform baked in.

| Pre-2.0 key | Canonical key | Relationship |
|---|---|---|
| `CDS_coverage_metric` | `cds_coverage` | same value |
| `CDS_coverage_metric_inframe_100read_100tx` | `cds_coverage_inframe_100read_100tx` | same value |
| `CDS_coverage_metric_inframe_1read_1000tx` | `cds_coverage_inframe_1read_1000tx` | same value |
| `CDS_coverage_metric_not_inframe_100read_100tx` | `cds_coverage_100read_100tx` | same value |
| `CDS_coverage_metric_not_inframe_1read_1000tx` | `cds_coverage_1read_1000tx` | same value |
| `cds_coverage_not_inframe_100read_100tx` | `cds_coverage_100read_100tx` | same value |
| `cds_coverage_not_inframe_1read_1000tx` | `cds_coverage_1read_1000tx` | same value |
| `multimapper_rate` | `rpf_multimapper_rate` | same value |
| `read_length_distribution_IQR_metric` | `read_length_iqr_fraction` | old = 1 − new |
| `read_length_distribution_bimodality_metric` | `read_length_bimodality_coefficient` | old = 1 / (1 + new) |
| `read_length_distribution_coefficient_of_variation_metric` | `read_length_cv` | old = 1 / (1 + new) |
| `read_length_distribution_maxprop_metric` | `read_length_max_proportion` | same value |
| `read_length_distribution_normality_metric` | `read_length_normality_pvalue` | old = 1 − new |
| `terminal_bias_kl_3prime` | `terminal_bias_kl_3prime` | old = 1 / (1 + new) |
| `terminal_bias_kl_3prime_raw` | `terminal_bias_kl_3prime` | same value |
| `terminal_bias_kl_3prime_score` | `terminal_bias_kl_3prime` | old = 1 / (1 + new) |
| `terminal_bias_kl_5prime` | `terminal_bias_kl_5prime` | old = 1 / (1 + new) |
| `terminal_bias_kl_5prime_raw` | `terminal_bias_kl_5prime` | same value |
| `terminal_bias_kl_5prime_score` | `terminal_bias_kl_5prime` | old = 1 / (1 + new) |
| `terminal_bias_maxabs_3prime` | `terminal_bias_max_deviation_3prime` | old = 1 − new |
| `terminal_bias_maxabs_5prime` | `terminal_bias_max_deviation_5prime` | old = 1 − new |
| `terminal_nucleotide_bias_distribution_3_prime_metric` | `terminal_bias_kl_3prime` | old = 1 / (1 + new) |
| `terminal_nucleotide_bias_distribution_5_prime_metric` | `terminal_bias_kl_5prime` | old = 1 / (1 + new) |
| `terminal_nucleotide_bias_max_absolute_metric_3_prime_metric` | `terminal_bias_max_deviation_3prime` | old = 1 − new |
| `terminal_nucleotide_bias_max_absolute_metric_5_prime_metric` | `terminal_bias_max_deviation_5prime` | old = 1 − new |
| `unique_rpf_rate` | `rpf_multimapper_rate` | old = 1 − new |

