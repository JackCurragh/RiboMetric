# RiboMetric Reporting Guide

This guide explains RiboMetric's output formats and how to use them for both pipeline integration and detailed sample review.

## Output Formats Overview

RiboMetric provides multiple output formats optimized for different use cases:

| Format | Use Case | Best For |
|--------|----------|----------|
| **Summary TSV** | Pipeline decisions | Quick pass/fail, multi-sample tracking |
| **QC Status JSON** | Automated QC | Machine-readable status with thresholds |
| **Comparison CSV** | Sample comparison | Side-by-side analysis of multiple samples |
| **Metrics Table CSV** | Detailed analysis | All metrics with read-length breakdowns |
| **HTML Report** | Sample review | Interactive, visual exploration |
| **PDF Report** | Documentation | Shareable, archivable reports |
| **Full JSON** | Reanalysis | Complete data for regenerating reports |

## For Pipeline Integration

### 1. Summary TSV - Quick Decisions

**Best for:** Automated pipelines, quick QC checks, multi-sample tracking

**Location:** `{sample}_summary.tsv`

**Format:** One line per sample, easily concatenated

```text
sample	timestamp	mode	total_reads	duplicate_rate	…	uniformity_entropy	…	cds_enrichment_ratio	…	periodicity_dominance	…
Sample1	2025-01-15T10:30:00	annotation	1500000	0.31	…	0.78	…	3.12	…	0.85	…
Sample2	2025-01-15T11:45:00	annotation	1200000	0.44	…	0.65	…	1.80	…	0.72	…
```

One column per metric in `results["metrics"]` (abridged above), each holding the
whole-library value. Which metrics appear depends on the inputs: a run without
an annotation or FASTA has fewer columns.

**Usage in pipelines:**

```bash
# Generate summary for each sample
for sample in *.bam; do
    RiboMetric run -b $sample -a annotation.tsv \
        --summary-tsv
done

# Concatenate all summaries (only safe when every sample used the same inputs,
# so the columns match)
cat *_summary.tsv | head -1 > all_samples_summary.tsv
tail -n +2 -q *_summary.tsv >> all_samples_summary.tsv

# Filter by column NAME: positions change with the inputs and between releases
awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; print; next}
            $col["periodicity_dominance"] >= 0.7' all_samples_summary.tsv > passing_samples.tsv
```

### 2. QC Status JSON - Automated Decisions

**Best for:** Automated QC gates, pipeline branching, status tracking

**Location:** `{sample}_qc_status.json`

**Format:**

```json
{
  "sample": "Sample1",
  "timestamp": "2025-01-15T10:30:00",
  "overall_status": "PASS",
  "checks": [
    {
      "metric": "periodicity_dominance_score",
      "source_metric": "periodicity_dominance",
      "value": 0.85,
      "score": 0.85,
      "status": "PASS",
      "gate": true,
      "tier": 1
    }
  ],
  "summary": {
    "total_checks": 5,
    "passed": 4,
    "warnings": 1,
    "failed": 0
  },
  "recommendation": "Sample passed all QC checks. Proceed with downstream analysis."
}
```

Each check pairs a score (`metric`, 0–1, higher-is-better) with the raw metric
it came from (`source_metric`, and its natural-units `value`). Only checks with
`"gate": true` — Tier 1 — decide `overall_status`; the rest are reported.

**Usage in pipelines:**

```python
import json

# Load QC status
with open('sample_qc_status.json') as f:
    qc = json.load(f)

# Make pipeline decision
if qc['overall_status'] == 'PASS':
    proceed_to_translation_analysis(sample)
elif qc['overall_status'] == 'WARNING':
    flag_for_review(sample)
    proceed_with_caution(sample)
else:  # FAIL
    exclude_from_analysis(sample)
    log_failure(sample, qc['recommendation'])
```

```bash
# Shell script example
STATUS=$(jq -r '.overall_status' sample_qc_status.json)

if [ "$STATUS" = "PASS" ]; then
    echo "Sample passed QC, continuing pipeline..."
    nextflow run translation_analysis.nf
elif [ "$STATUS" = "WARNING" ]; then
    echo "Sample has warnings, flagging for review..."
    echo "sample" >> samples_for_review.txt
else
    echo "Sample failed QC, excluding from analysis"
    echo "sample" >> failed_samples.txt
fi
```

### 3. Comparison CSV - Multi-Sample Analysis

**Best for:** Comparing metrics across many samples

**Location:** `{sample}_comparison.csv`, one file per sample

**Format:** Wide format, one column per value. Dict-valued metrics are flattened
with a suffix — `_global` for the whole library, `_rl<N>` for read length N —
so samples with different read lengths have different columns. Combine them by
column name, never by position.

```text
sample,timestamp,mode,duplicate_rate,…,uniformity_entropy_global,uniformity_entropy_rl28,…,cds_enrichment_ratio,…,periodicity_dominance_global,…
Sample1,2025-01-15T10:30:00,annotation,0.31,…,0.78,0.74,…,3.12,…,0.85,…
Sample2,2025-01-15T11:45:00,annotation,0.44,…,0.65,0.61,…,1.80,…,0.72,…
```

**Usage:**

```r
# R analysis
library(tidyverse)

# Load every per-sample file; bind_rows() aligns columns by name, which matters
# because per-read-length columns differ between samples
metrics <- list.files("comparison_results", pattern = "_comparison\\.csv$",
                      full.names = TRUE) |>
  map(read_csv, show_col_types = FALSE) |>
  bind_rows()

# Quick overview
metrics %>%
  select(sample, periodicity_dominance_global,
         uniformity_entropy_global, cds_enrichment_ratio) %>%
  summary()

# Identify outliers
outliers <- metrics %>%
  filter(periodicity_dominance_global < 0.5 |
         uniformity_entropy_global < 0.5)

# Plot distributions
metrics %>%
  ggplot(aes(x = periodicity_dominance_global)) +
  geom_histogram() +
  geom_vline(xintercept = 0.7, color = "red", linetype = "dashed")
```

## For Sample Review

### 1. Improved HTML Report - Interactive Review

**Best for:** Detailed sample inspection, identifying issues, generating figures

**Features:**
- **Executive Summary** - QC status at a glance
- **Key Metrics Dashboard** - Critical metrics with pass/warn/fail indicators
- **Searchable Metrics Table** - Find any metric quickly
- **Interactive Plots** - Zoom, pan, hover for details
- **Smooth Navigation** - Jump to any section
- **Print-Friendly** - Export clean PDFs via browser

**Access:** Open `{sample}_RiboMetric.html` in web browser

**Navigation:**
- Left sidebar links to all sections
- Click any metric in summary to jump to details
- Use search box to filter metrics table
- Scroll triggers automatic section highlighting

### 2. Metrics Table CSV - Detailed Analysis

**Best for:** Spreadsheet analysis, custom visualizations, read-length specific investigation

**Location:** `{sample}_metrics_table.csv`

**Format:**

```text
sample,metric,read_length_or_region,value,description
Sample1,periodicity_dominance,global,0.85,Fraction of coding A-sites in the dominant reading frame; the global value uses one shared dominant frame.
Sample1,periodicity_dominance,28,0.82,Fraction of coding A-sites in the dominant reading frame; the global value uses one shared dominant frame.
Sample1,periodicity_dominance,29,0.87,Fraction of coding A-sites in the dominant reading frame; the global value uses one shared dominant frame.
Sample1,uniformity_entropy,global,0.78,Normalised entropy of the codon-binned start-codon metagene.
...
```

Descriptions come from the metric registry, the same source that generates
[METRICS.md](METRICS.md). Besides `global` and per-read-length rows, some metrics
add summary rows such as `global_by_read_length_max`.

**Usage:**

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load detailed metrics
metrics = pd.read_csv('sample_metrics_table.csv')

# Analyze read-length specific periodicity
periodicity = metrics[metrics['metric'] == 'periodicity_dominance']
# per-read-length rows have a numeric label; 'global' and summary rows do not
per_length = periodicity['read_length_or_region'].astype(str).str.isdigit()
periodicity = periodicity[per_length].copy()
periodicity['read_length'] = periodicity['read_length_or_region'].astype(int)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(periodicity['read_length'], periodicity['value'], marker='o')
plt.axhline(y=0.7, color='r', linestyle='--', label='Pass threshold')
plt.xlabel('Read Length')
plt.ylabel('Periodicity Dominance')
plt.title('Periodicity by Read Length')
plt.legend()
plt.savefig('periodicity_by_readlength.png')
```

## Generating Improved Outputs

### Command-Line Usage

```bash
# HTML, JSON, CSV and PDF
RiboMetric run -b sample.bam -a annotation.tsv --all

# Add improved outputs
RiboMetric run -b sample.bam -a annotation.tsv \
    --summary-tsv \
    --qc-status \
    --comparison-csv \
    --metrics-table

# All four pipeline outputs at once
RiboMetric run -b sample.bam -a annotation.tsv --improved-outputs
```

### Python API Usage

The pipeline outputs can be regenerated from any saved result, for example
with different thresholds, without re-running the analysis:

```python
import json

from RiboMetric.results_output import generate_all_outputs

with open("sample_RiboMetric.json") as fh:
    data = json.load(fh)

generate_all_outputs(
    results_dict=data["results"],
    config=data["config"],
    sample_name="Sample1",
    output_directory="./results",
)
```

Passing `thresholds=` makes the QC status an explicit required-check policy:
every metric it names must be present, or that check fails.

## Customizing QC Thresholds

Report thresholds live in the `scoring:` section of `config.yml` (or a custom
config passed via `--config`). Each entry is keyed by the score it produces,
names the raw `metric` it reads, and sets `pass`/`warn` on the 0–1 score scale;
`gate: true` puts it in the Tier 1 verdict. These drive both the HTML report and
`qc_status.json`, so the two always agree. The entries below are copied from the
shipped `config.yml`:

```yaml
scoring:
  periodicity_dominance_score:
    metric: periodicity_dominance
    method: identity
    status: {pass: 0.70, warn: 0.50}
    gate: true
    tier: 1
  cds_enrichment_score:
    metric: cds_enrichment_ratio
    method: enrichment_ratio
    status: {pass: 0.60, warn: 0.30}
    gate: true
    tier: 1
  coverage_uniformity_score:
    metric: uniformity_entropy
    method: identity
    status: {pass: 0.60, warn: 0.30}
    gate: false
    tier: 2
  terminal_evenness_kl_5prime_score:
    metric: terminal_bias_kl_5prime
    method: inverse_linear
    params: {max_value: 2.0}
    status: {pass: 0.70, warn: 0.40}
    gate: false
    tier: 3
```

`RiboMetric evaluate -e thresholds.yml` is separate: its YAML sets pass/warn on
*raw* metrics, and each metric's direction comes from the metric registry.

Use with:

```bash
RiboMetric run -b sample.bam -a annotation.tsv \
    --config my_config.yml
```

## Best Practices

### For Pipeline Developers

1. **Use Summary TSV** for quick pass/fail decisions
2. **Use QC Status JSON** for structured pipeline logic
3. **Archive HTML reports** for later review of failed samples
4. **Log all outputs** to a central database for tracking
5. **Set appropriate thresholds** based on your protocol

### For Researchers

1. **Start with HTML Executive Summary** - Get overview instantly
2. **Review Key Metrics Dashboard** - Identify problem areas
3. **Use Comparison CSV** - Compare across samples/conditions
4. **Deep dive with Metrics Table** - Investigate specific read lengths
5. **Save HTML reports** - Document quality for publications

### For Core Facilities

1. **Provide HTML reports** to users for review
2. **Use QC Status** for automated sample acceptance
3. **Track metrics over time** with Comparison CSV
4. **Standardize thresholds** across projects
5. **Generate PDF reports** for archival records

## Interpreting Metrics

Every score, the metric it comes from, its tier and its pass/warn thresholds are
listed in [METRICS.md](METRICS.md). That page is generated from the metric
registry and the live scoring spec, so it cannot drift from the code — which the
hand-kept tables that used to be here had.

In short, **Tier 1** (CDS enrichment, periodicity dominance, periodicity information) is gated: a Tier 1 failure fails the sample.
**Tier 2** (coverage uniformity, library saturation, usable read fraction) says whether the library suits a given analysis. **Tier 3**
(alignment unique mapping, footprint homogeneity, fragment uniqueness, RPF unique mapping, terminal evenness KL 3prime, terminal evenness KL 5prime, terminal evenness maxdev 3prime, terminal evenness maxdev 5prime, terminal integrity 5prime) is informational.

### Warning Signs

Raw-value equivalents of the shipped thresholds:

- **Weak periodicity:** `periodicity_dominance` below 0.50 fails the Tier 1 gate (0.70 to pass); one third is the random baseline.
- **Low CDS enrichment:** E below ~1.4 fails, and E of about 2.5 is needed to pass.
- **Uneven coverage:** `uniformity_entropy` below 0.30 fails (0.60 to pass).
- **Strong terminal bias:** 5′ KL above 1.2 bits fails and above 0.6 bits warns; 2.0 bits or more scores zero.

## Example Workflows

### Workflow 1: Quick QC in Pipeline

```bash
#!/bin/bash
# Quick QC check before expensive analysis

for bam in data/*.bam; do
    sample=$(basename $bam .bam)

    # Run RiboMetric
    RiboMetric run -b $bam -a annotation.tsv \
        --qc-status \
        -o qc_results/

    # Check status
    status=$(jq -r '.overall_status' qc_results/${sample}_qc_status.json)

    if [ "$status" = "PASS" ]; then
        echo "$sample: PASS - submitting to translation analysis"
        sbatch run_translation.sh $bam
    else
        echo "$sample: $status - skipping"
    fi
done
```

### Workflow 2: Batch Sample Comparison

```bash
#!/bin/bash
# Each run writes its own comparison_results/{sample}_comparison.csv
for bam in data/*.bam; do
    RiboMetric run -b $bam -a annotation.tsv \
        --comparison-csv \
        -o comparison_results/
done

# Combine by column name and analyse (see the R example above)
Rscript compare_samples.R comparison_results/
```

### Workflow 3: Detailed Review

```bash
# Generate full reports for final samples
RiboMetric run -b sample.bam -a annotation.tsv \
    --improved-outputs \
    --html --pdf \
    -o final_reports/

# Open HTML for review
open final_reports/sample_RiboMetric.html

# Export key metrics
cat final_reports/sample_summary.tsv >> project_summary.tsv
```

## Troubleshooting

### Issue: Metrics seem off

**Solution:** Check the metrics table CSV for read-length specific values. Some metrics may be skewed by specific read lengths.

### Issue: QC status doesn't match expectations

**Solution:** Adjust the `scoring:` block of your config — it drives both the HTML report and `qc_status.json`. For `RiboMetric evaluate`, pass a thresholds YAML with `-e`. Different protocols may need different cutoffs.

### Issue: HTML report too large

**Solution:** Use `--subsample` to analyze fewer reads, or export to PDF for sharing.

### Issue: Can't compare samples

**Solution:** Combine the per-sample comparison CSVs by column name (`bind_rows()` in R, `pd.concat()` in pandas). Per-read-length columns (`_rl<N>`) differ between samples, so never concatenate them positionally.

## Further Reading

- [METRICS.md](METRICS.md) - Every metric and score, generated from the registry
- [METRIC_NAMING.md](METRIC_NAMING.md) - The 2.0 naming contract and the pre-2.0 key mapping
- `RiboMetric --help` - Command-line options
- [Documentation](https://ribometric.readthedocs.io) - Full online docs
