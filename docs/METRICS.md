# RiboMetric Metrics Guide

RiboMetric provides two tiers of metrics: **default** (standard Ribo-Seq QC) and **optional** (theoretical/experimental approaches).

## Default Metrics

These metrics are calculated by default and represent the standard quality checks expected for Ribo-Seq data:

### Read Length Distribution
- `read_length_distribution_IQR` - Interquartile range normalized metric
- `read_length_distribution_coefficient_of_variation` - CV of read length distribution
- `read_length_distribution_maxprop` - Proportion of reads at most frequent length

### Terminal Nucleotide Bias
- `terminal_nucleotide_bias_KL_5prime` - Kullback-Leibler divergence for 5' end
- `terminal_nucleotide_bias_KL_3prime` - Kullback-Leibler divergence for 3' end
- `terminal_nucleotide_bias_max_absolute_5prime` - Maximum deviation from expected for 5' end
- `terminal_nucleotide_bias_max_absolute_3prime` - Maximum deviation from expected for 3' end

### 3-nt Periodicity
- **`periodicity_dominance`** - Proportion of reads in dominant reading frame (recommended)
- `periodicity_information` - Shannon information content-based periodicity score

### Metagene Uniformity
- **`uniformity_entropy`** - Entropy-based uniformity across metagene profile (recommended)

### Coverage & Regional Distribution
- `CDS_coverage` - Proportion of CDS covered by reads
- `region_ratios` - Ratios between different mRNA regions (CDS:5'UTR, CDS:3'UTR, etc.)
- `region_proportions` - Proportion of reads in each region

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
- `read_length_distribution_bimodality` - Bimodality coefficient
- `read_length_distribution_normality` - Normality test p-value

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
- **Uniformity**: Use `uniformity_entropy` (standard information-theoretic approach)

The optional metrics are useful for:
- Method comparison studies
- Developing new quality standards
- Research on Ribo-Seq quality assessment
- Exploring alternative analytical approaches

## Philosophy

RiboMetric provides objective measurements without labeling datasets as "good" or "bad". Users and communities can develop their own interpretation standards based on these metrics.
