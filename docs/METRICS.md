# RiboMetric Metrics Guide

RiboMetric provides two categories of metrics: **default** (standard Ribo-Seq QC) and **optional** (theoretical/experimental approaches).

The HTML report organises scored metrics into three report tiers — Tier 1 (Ribo-seq identity, gated), Tier 2 (usability), Tier 3 (technical caveats) — plus a Diagnostics section for context-dependent values that are shown as raw measurements without pass/fail scoring.

## Default Metrics

These metrics are calculated by default and represent the standard quality checks expected for Ribo-Seq data:

### Read Length Distribution
- `read_length_distribution_IQR_metric` - Interquartile range normalized metric
- `read_length_distribution_coefficient_of_variation_metric` - CV of read length distribution
- `read_length_distribution_maxprop_metric` - Proportion of reads at most frequent length

### Terminal Nucleotide Bias
- `terminal_bias_kl_5prime_score` - Normalized 5' terminal nucleotide bias score derived from KL divergence; 1 means observed matches expected
- `terminal_bias_kl_3prime_score` - Normalized 3' terminal nucleotide bias score derived from KL divergence; 1 means observed matches expected
- `terminal_bias_kl_5prime_raw` - Raw 5' terminal Kullback-Leibler divergence in bits
- `terminal_bias_kl_3prime_raw` - Raw 3' terminal Kullback-Leibler divergence in bits
- `terminal_nucleotide_bias_max_absolute_5prime` - Maximum deviation from expected for 5' end
- `terminal_nucleotide_bias_max_absolute_3prime` - Maximum deviation from expected for 3' end

### 3-nt Periodicity
- **`periodicity_dominance`** - Proportion of reads in the dominant reading frame; the global value uses one shared dominant frame across all read lengths (recommended)
- `periodicity_information` - Shannon information content-based periodicity score

### Metagene Uniformity
- **`uniformity_entropy`** - Entropy-based uniformity across the codon-binned start-codon metagene window, not whole-CDS coverage uniformity (recommended)

### CDS Enrichment
- `cds_enrichment_ratio` - Ratio of observed CDS-body read fraction to the length-weighted expected fraction (E = observed / expected). E > 1 means reads are enriched in CDS relative to random sampling across transcript lengths. Used as a Tier 1 gate metric in place of the raw `prop_reads_CDS` proportion, which is confounded by transcript-length distribution.

### Coverage & Regional Distribution
- `CDS_coverage` - Proportion of CDS covered by reads
- `region_ratios` - Ratios between different mRNA regions (CDS:5'UTR, CDS:3'UTR, etc.)
- `region_proportions` - Proportion of reads in each region

### Alignment Quality
- `duplicate_rate` - Fraction of reads that are PCR/library duplicates (from collapsed count column)
- `multimapper_rate` - Backwards-compatible alias for `rpf_multimapper_rate`
- `rpf_multimapper_rate` - Fraction of weighted ribosome protected fragments reported at more than one alignment location; uses `NH:i:N` when present, then XA, and falls back to STAR-style `MAPQ < 255`. On a **transcriptome** BAM this counts alignment/transcript multiplicity (a read mapping to multiple isoforms of one gene), not distinct genomic loci — interpret accordingly.
- `unique_rpf_rate` - Fraction of weighted ribosome protected fragments that are uniquely mapped
- `alignment_multimapper_rate` - Fraction of reported alignment rows whose read/fragment has evidence of another reported alignment (same caveat: transcript multiplicity on a transcriptome BAM)
- `soft_clip_rate_5prime` - Fraction of reads with 5′ soft-clipping (adapter/trimming artefact indicator)

### Di-some Detection
- `read_length_distribution_bimodality_metric` - Bimodality coefficient of the read length distribution
- `disome_proportion` - Fraction of reads with length 50–70 nt (di-some / collision footprint window)

### Codon-level Translation Metrics
- `stop_codon_readthrough_ratio` - Reads downstream of stop codon (positions +1…+30) relative to upstream (-30…-1); elevated values indicate readthrough or frameshifting
- `start_codon_enrichment_ratio` - Reads near start codon (-5…+20) relative to CDS body (positions +30…+50); very high values suggest harringtonine/LTM treatment or initiation stalling

## Optional Metrics

These metrics are **only calculated when explicitly enabled** using `--enable-optional-metrics` or `--enable-metric <name>`. They represent alternative or more theoretical approaches:

### Alternative Periodicity Metrics
- `periodicity_autocorrelation` - Autocorrelation-based periodicity detection
- `periodicity_fourier` - Fourier transform power at 3-nt frequency
- `periodicity_trips_viz` - Trips-Viz algorithm score

### Alternative Uniformity Metrics
- `uniformity_autocorrelation` - Autocorrelation-based uniformity
- `uniformity_theil_index` - Theil index across metagene profile
- `uniformity_gini_index` - Gini coefficient across metagene profile

### Additional Read Length Metrics
- `read_length_distribution_normality` - Normality test p-value

### RUST Metric (requires `--fasta`)
- `rust_mean_kl_divergence` - Mean KL divergence across the 60-codon RUST metagene window; measures the extent of codon-specific A-site accumulation (O'Connor et al., Nat Commun 2016)

## Usage Examples

### Default behavior (standard metrics only):
```bash
RiboMetric run -b sample.bam -a annotation.tsv
```

### Enable all optional metrics:
```bash
RiboMetric run -b sample.bam -a annotation.tsv --enable-optional-metrics
```

### Enable specific optional metrics:
```bash
RiboMetric run -b sample.bam -a annotation.tsv \
  --enable-metric periodicity_fourier \
  --enable-metric uniformity_gini_index
```

## Recommendations

For standard Ribo-Seq QC, the **default metrics** provide comprehensive coverage:
- **Periodicity**: Use `periodicity_dominance` (simple, interpretable)
- **Uniformity**: Use `uniformity_entropy` for local metagene-window uniformity; use `CDS_coverage` for whole-CDS coverage

The optional metrics are useful for:
- Method comparison studies
- Developing new quality standards
- Research on Ribo-Seq quality assessment
- Exploring alternative analytical approaches

## Philosophy

RiboMetric scores each metric on a 0–1 scale anchored to interpretable reference points (0 = random-chance baseline, 1 = ideal). Scores drive pass/warn/fail labels in the HTML report and QC gate, but the raw value is always shown alongside — the score is a convenience, not a substitute for reading the number. Thresholds are provisional and will be recalibrated against a reference corpus in a future release; the `scoring:` block in `config.yml` exposes every threshold for local adjustment.
