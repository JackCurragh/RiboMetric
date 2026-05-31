---
title: 'RiboMetric: A comprehensive quality assessment tool for ribosome profiling datasets'
tags:
  - Python
  - ribosome profiling
  - Ribo-Seq
  - translation
  - quality control
  - bioinformatics
authors:
  - name: Jack Tierney
    affiliation: 1
  - name: Lukas Wierdsma
    affiliation: 2
affiliations:
  - name: LAPTI Lab, School of Biochemistry and Cell Biology, University College Cork, Ireland
    index: 1
  - name: Howest University of Applied Sciences, Bruges, Belgium
    index: 2
date: 31 May 2026
bibliography: paper.bib
---

# Summary

Ribosome profiling (Ribo-Seq) is a sequencing technique that captures the
positions of actively translating ribosomes at nucleotide resolution
[@ingolia2009]. The resulting data encodes the translational landscape of a
cell, but like all high-throughput sequencing assays it is subject to a range
of technical artefacts—library preparation biases, ligation artefacts,
suboptimal ribosome protection, and alignment ambiguities—that must be
evaluated before downstream analysis. Poor-quality Ribo-Seq libraries yield
attenuated 3-nucleotide periodicity, atypical read-length distributions, or
elevated off-CDS signal, all of which can confound translation efficiency
estimates, ORF discovery, and codon-usage analyses.

`RiboMetric` is a Python command-line tool that automates this assessment.
Given a coordinate-sorted BAM file and an optional transcript annotation, it
produces a rich HTML report, machine-readable metric outputs (JSON, CSV, TSV),
and pipeline-ready QC gates, covering all major aspects of Ribo-Seq data
quality.

# Statement of Need

Several existing tools assess Ribo-Seq quality in isolation: RiboProfiling
[@popa2016] provides a basic R-based workflow, and riboWaltz [@lauria2018]
optimises P-site offsets and reports frame occupancy as a by-product.
Tools such as RiboTaper [@calviello2016] compute translated ORF evidence
but are not primarily QC-focused. However, no single tool delivers the full
spectrum of QC metrics—from raw read properties through periodicity and
regional distribution to automated pass/fail thresholds—in a
pipeline-native Python package.

`RiboMetric` fills this gap. It is the only tool that combines:

1. **Comprehensive metrics** in a single run: read-length distribution,
   terminal nucleotide ligation bias, 3-nt periodicity (multiple methods),
   metagene uniformity, CDS coverage, and regional read distribution.
2. **Multimapper-aware calculations** — for STAR transcriptome alignments,
   frame-sensitive metrics are restricted to uniquely-mapped reads
   (`MAPQ = 255`) by default, eliminating a known source of periodicity
   underestimation.
3. **Pipeline integration** — an `evaluate` subcommand accepts a result file
   and a YAML of expected thresholds and exits with a coded status (0/1/2 for
   PASS/WARNING/FAIL), enabling downstream gating in Nextflow or Snakemake
   workflows.
4. **Interactive exploration** — a `view` subcommand launches a terminal UI
   for navigating results without a browser.
5. **Accessibility** — pure Python, installable via `pip install RiboMetric`,
   with all system dependencies documented and a Docker image published on
   every tagged release.

# Implementation

`RiboMetric` is structured around three subcommands:

- `prepare` — converts a GFF3/GTF annotation into a compact TSV used by the
  `run` subcommand. The conversion is fully vectorised using `numpy` and
  `pandas`, and supports multi-threaded processing of large annotations.
- `run` — parses a BAM file with `oxbow` (Rust-backed, zero-copy) and
  `pysam`, calculates metrics, and writes requested outputs. Read lengths are
  derived from CIGAR consumed-bases rather than sequence length, ensuring
  correctness for soft-clipped alignments.
- `evaluate` — scores a previously produced result (JSON or metrics-table CSV)
  against user-specified thresholds and returns an exit code appropriate for
  pipeline automation.

Metric calculations are grouped into default (standard QC) and optional
(theoretical/experimental) sets. Optional metrics—including alternative
periodicity estimators (autocorrelation, Fourier, Trips-Viz) and uniformity
indices (Theil, Gini)—are enabled with `--enable-optional-metrics` or
`--enable-metric <name>`. All calculations respect read deduplication weights
stored in the `count` column, avoiding double-counting in collapsed libraries.

Outputs include interactive HTML reports (rendered with `Jinja2` and
`Plotly`), archivable PDF reports, and multiple tabular formats suited to
multi-sample comparison or integration with laboratory information management
systems.

# Availability

`RiboMetric` is available on PyPI (`pip install RiboMetric`) and as a
multi-architecture Docker image (`ghcr.io/jackcurragh/ribometric`). Source
code, documentation, and issue tracker are at
<https://github.com/JackCurragh/RiboMetric> under the MIT licence.
Documentation is hosted at <https://ribometric.readthedocs.io>.

# Acknowledgements

The authors thank the LAPTI lab at University College Cork for discussion and
test data, and all contributors who raised issues or submitted pull requests.

# References
