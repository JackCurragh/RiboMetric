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

```tsv
sample	timestamp	mode	total_reads	periodicity_dominance	uniformity_entropy	prop_reads_CDS
Sample1	2025-01-15T10:30:00	annotation	1500000	0.85	0.78	0.82
Sample2	2025-01-15T11:45:00	annotation	1200000	0.72	0.65	0.75
```

**Usage in pipelines:**

```bash
# Generate summary for each sample
for sample in *.bam; do
    RiboMetric run -b $sample -a annotation.tsv \
        --output-summary-tsv
done

# Concatenate all summaries
cat *_summary.tsv | head -1 > all_samples_summary.tsv
tail -n +2 -q *_summary.tsv >> all_samples_summary.tsv

# Filter passing samples
awk -F'\t' '$5 > 0.7 && $6 > 0.7' all_samples_summary.tsv > passing_samples.tsv
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
      "metric": "periodicity_dominance",
      "value": 0.85,
      "status": "PASS",
      "threshold_pass": 0.7,
      "threshold_warn": 0.5
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

**Location:** `{sample}_comparison.csv` (appends to same file)

**Format:** Wide format with all metrics as columns

```csv
sample,timestamp,periodicity_dominance_global,uniformity_entropy_global,prop_reads_CDS_global,...
Sample1,2025-01-15T10:30:00,0.85,0.78,0.82,...
Sample2,2025-01-15T11:45:00,0.72,0.65,0.75,...
Sample3,2025-01-15T13:00:00,0.91,0.82,0.88,...
```

**Usage:**

```r
# R analysis
library(tidyverse)

# Load comparison data
metrics <- read_csv("all_samples_comparison.csv")

# Quick overview
metrics %>%
  select(sample, periodicity_dominance_global,
         uniformity_entropy_global, prop_reads_CDS_global) %>%
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

```csv
sample,metric,read_length_or_region,value,description
Sample1,periodicity_dominance,global,0.85,Proportion of reads in dominant reading frame
Sample1,periodicity_dominance,28,0.82,Proportion of reads in dominant reading frame
Sample1,periodicity_dominance,29,0.87,Proportion of reads in dominant reading frame
Sample1,uniformity_entropy,global,0.78,Shannon entropy of metagene distribution
...
```

**Usage:**

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load detailed metrics
metrics = pd.read_csv('sample_metrics_table.csv')

# Analyze read-length specific periodicity
periodicity = metrics[metrics['metric'] == 'periodicity_dominance']
periodicity = periodicity[periodicity['read_length_or_region'] != 'global']
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
# Standard outputs (HTML, JSON, CSV)
RiboMetric run -b sample.bam -a annotation.tsv --all

# Add improved outputs
RiboMetric run -b sample.bam -a annotation.tsv \
    --output-summary-tsv \
    --output-qc-status \
    --output-comparison \
    --output-metrics-table

# Generate all formats
RiboMetric run -b sample.bam -a annotation.tsv --all-outputs
```

### Python API Usage

```python
from RiboMetric import run_analysis
from RiboMetric.results_output_improved import generate_all_outputs

# Run analysis
results = run_analysis('sample.bam', 'annotation.tsv')

# Generate all improved outputs
generate_all_outputs(
    results_dict=results,
    config=config,
    sample_name='Sample1',
    output_directory='./results',
    thresholds={
        'periodicity_dominance': {'pass': 0.7, 'warn': 0.5},
        'uniformity_entropy': {'pass': 0.7, 'warn': 0.5},
    }
)
```

## Customizing QC Thresholds

Create a `qc_thresholds.yaml` file:

```yaml
periodicity_dominance:
  pass: 0.7  # > 0.7 = PASS
  warn: 0.5  # 0.5-0.7 = WARNING, < 0.5 = FAIL

uniformity_entropy:
  pass: 0.7
  warn: 0.5

prop_reads_CDS:
  pass: 0.7
  warn: 0.5

read_length_distribution_IQR_metric:
  pass: 0.3
  warn: 0.2

terminal_nucleotide_bias_distribution_5_prime_metric:
  pass: 1.0  # Lower is better (less bias)
  warn: 2.0
```

Use with:

```bash
RiboMetric run -b sample.bam -a annotation.tsv \
    --qc-thresholds qc_thresholds.yaml
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

### Critical Metrics (Always Review)

| Metric | Good Range | Interpretation |
|--------|------------|----------------|
| `periodicity_dominance` | > 0.7 | Strong 3-nt periodicity, good ribosome footprints |
| `uniformity_entropy` | > 0.7 | Uniform CDS coverage, no major biases |
| `prop_reads_CDS` | > 0.7 | Most reads map to coding regions |
| `terminal_nucleotide_bias_5prime` | < 1.0 | Low adapter ligation bias |

### Warning Signs

- **Low periodicity** (< 0.5): RNA contamination, poor digestion, wrong read lengths
- **Low uniformity** (< 0.5): Biased coverage, PCR artifacts, degradation
- **Low CDS proportion** (< 0.5): rRNA contamination, poor mapping
- **High terminal bias** (> 2.0): Adapter ligation artifacts

## Example Workflows

### Workflow 1: Quick QC in Pipeline

```bash
#!/bin/bash
# Quick QC check before expensive analysis

for bam in data/*.bam; do
    sample=$(basename $bam .bam)

    # Run RiboMetric
    RiboMetric run -b $bam -a annotation.tsv \
        --output-qc-status \
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
# Generate comparison across all samples

# Initialize comparison file
> all_samples_comparison.csv

# Process each sample
for bam in data/*.bam; do
    sample=$(basename $bam .bam)
    RiboMetric run -b $bam -a annotation.tsv \
        --output-comparison \
        -o comparison_results/
done

# Analyze in R
Rscript compare_samples.R all_samples_comparison.csv
```

### Workflow 3: Detailed Review

```bash
# Generate full reports for final samples
RiboMetric run -b sample.bam -a annotation.tsv \
    --all-outputs \
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

**Solution:** Review and adjust thresholds in `qc_thresholds.yaml`. Different protocols may need different cutoffs.

### Issue: HTML report too large

**Solution:** Use `--subsample` to analyze fewer reads, or export to PDF for sharing.

### Issue: Can't compare samples

**Solution:** Use the comparison CSV format which has consistent columns across samples.

## Further Reading

- [METRICS.md](METRICS.md) - Detailed metric descriptions
- [QUICK_REFERENCE.md](QUICK_REFERENCE.md) - Command-line options
- [Documentation](https://ribometric.readthedocs.io) - Full online docs
