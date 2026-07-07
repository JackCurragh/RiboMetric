# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [1.4.1] — 2026-07-07

### Added

- **Production run provenance** — JSON outputs now include a `results["provenance"]`
  block with timestamp, package/Python/platform details, effective config SHA-256,
  config file fingerprint, and input file fingerprints. Small files are fully
  SHA-256 hashed; large files receive an explicit sampled fingerprint unless
  `RIBOMETRIC_FULL_INPUT_HASH=1` is set.
- **Offset audit record** — JSON outputs now include `results["offsets"]`, recording
  the offset source, target, bounds, method, computed offsets, global/external
  offset inputs, frame-adjustment details, and the offsets actually applied by
  read length.
- **Applied offsets TSV** — RiboMetric now writes `{sample}_offsets.tsv` by default
  (`offsets_tsv: True`). The table is suitable for cohort-level offset checks and
  includes applied offsets, read counts, computed/global offsets, and frame
  adjustment fields.

### Changed

- **JSON output is on by default** so production runs retain the complete audit
  object without requiring `--json`.
- **`--output-offsets` now writes the applied-offset audit TSV schema** instead of
  only writing internally computed offsets. This makes the explicit path useful
  for global, external read-length, and read-specific offset modes as well.

## [1.4.0] — 2026-06-11

### Added

- **Unified scoring module** (`RiboMetric/scoring.py`) — all metrics now resolved
  through a single config-driven scorer. Every score is anchored 0 (random/null
  baseline) → 1 (ideal). Raw value travels with the score in all outputs. Score
  methods: `identity`, `one_minus_rate`, `inverse_linear`, `enrichment_ratio`.
  Gate (`gate: true`) determines which metrics contribute to the overall pass/fail
  verdict; only Tier 1 metrics are gated by default.

- **CDS enrichment ratio** (`cds_enrichment_ratio`) — replaces raw `prop_reads_CDS`
  as the Tier 1 CDS gate. Formula: E = observed_CDS_body_fraction /
  length-weighted_expected_fraction. Normalises for transcript-length distribution,
  removing the organism-level confound present in the raw proportion.

- **New terminal-bias scored metrics** — `terminal_bias_kl_5prime_raw` and
  `terminal_bias_kl_3prime_raw` now scored via `inverse_linear` (KL_max = 2.0 bits;
  KL = 1 bit → score 0.5, ≥ 2 bits → score 0). Raw bits shown in report.
  `terminal_bias_maxabs_5prime` / `_3prime` scored via `one_minus_rate`.

- **Tier 1/2/3 report layout** — HTML report reorganised into three labelled
  sections: Tier 1 (Ribo-seq identity, gated), Tier 2 (usability), Tier 3
  (technical caveats). Raw value and anchored score shown side-by-side.

- **Context strip** — report header now shows library type, dominant read length,
  total reads, and annotation source at a glance.

- **Diagnostics section** — nine context-dependent metrics demoted from scored
  table to a plain-value diagnostics section: `disome_proportion`,
  `read_length_distribution_*_metric` (IQR, CV, max-prop, bimodality),
  `start_codon_enrichment_ratio`, `stop_codon_readthrough_ratio`,
  `five_prime_ramp_ratio`, `three_prime_drop_ratio`. No pass/fail badge; captions
  are library-type-aware (e.g. disome_proportion caption differs for disome vs
  elongation experiments).

### Changed

- **`periodicity_dominance` scoring** — switched from `frame_dominance_rescaled`
  (`(d − 1/3) / (2/3)`) to `identity` (raw in-frame fraction). Score now equals
  the interpretable number a user already reads (0.79 = 79 % of A-sites in frame).
  The 1/3 random-frame baseline is shown as a dashed reference line on the
  periodicity plot rather than subtracted into the score.

- **`scoring:` block in `config.yml`** — all default thresholds now live here and
  override code defaults. User-supplied configs are merged over this block, so
  individual thresholds can be changed without touching code.

- **Overall QC gate** — only Tier 1 metrics (`gate: true`) determine
  `overall_status`. Tier 2 / Tier 3 failures are visible in the report but do not
  block a PASS verdict.

## [1.3.0] — 2026-06-10

### Added

- **A-site codon dwell-times / pause sites** (requires ``--fasta``) — per-codon
  A-site occupancy relative to codon abundance, summarised as ``codon_dwell_cv``,
  ``codon_dwell_p90_p10`` (dynamic range), ``proline_dwell`` and ``cga_dwell``.
  Surfaces translation pausing distinct from RUST's information-divergence view,
  with a sorted dwell-time plot (proline / CGA highlighted).
- **FLOSS read-length heterogeneity** — per-transcript Fragment Length
  Organization Similarity Score vs the library aggregate, summarised at the
  sample level (``floss_median``, ``floss_aberrant_transcript_fraction``) with a
  histogram. A library-level homogeneity QC, not an ORF/translation classifier.
- **Recommended read lengths** — RiboMetric now reports which read lengths carry
  clean 3-nt periodicity (and their P-site offsets), the count, and the fraction
  of the library they represent (`recommended_read_proportion`). The periodicity
  cutoff is configurable (`qc.recommend.min_periodicity` / `--min-periodicity`).
  A "Recommended Read Lengths" plot highlights the selected lengths.
- **Gene-body coverage ramp** — a metagene of A-site density across *relative*
  CDS position (0–100%), exposing the 5′ translation ramp and 3′ drop-off, with
  `five_prime_ramp_ratio` / `three_prime_drop_ratio` sample metrics and a plot.
- **Library complexity / saturation** — an analytic rarefaction curve of distinct
  A-site positions vs depth (`marginal_position_discovery_rate`,
  `complexity_distinct_positions`) plus a plot, answering "was this sequenced
  deeply enough?" without random subsampling.
- **Library-type classification** — heuristic `elongation` / `initiation` /
  `low_quality` label from periodicity, CDS enrichment and start-codon
  enrichment, with the supporting evidence recorded.
- `evaluate` ships a default threshold for `recommended_read_proportion`; all new
  scalar metrics are thresholdable via the existing `-e` YAML.

### Fixed

- **GFF coordinate off-by-one** (`file_parser.gff_df_to_cds_df`) — exon lengths
  were computed as `end - start` despite GFF3 being 1-based inclusive. This
  undercounted `transcript_length` by one nt per exon and, because the same
  arithmetic fed the leader/trailer sums, shifted `cds_start` by the number of
  complete UTR exons — scrambling the reading frame (and therefore periodicity)
  for multi-exon transcripts. CDS spans are now whole codons as expected.
  **Annotation TSVs produced by `prepare` should be regenerated.**
- **Terminal-bias 3′ background** (`bam_processing.calculate_background`) — the
  `five_prime` flag was ignored, so the 5′ and 3′ ligation-bias backgrounds were
  identical and the 3′ metric used the wrong reference distribution.
- **`--offset-global` was silently ignored** — the CLI wrote `offset_global`
  while the pipeline only read `global_offset`. The flag now applies a fixed
  offset to all read lengths as documented; omit it to auto-calculate.
- **`evaluate` ignored metric directionality** — lower-is-better metrics
  (`duplicate_rate`, `multimapper_rate`, `disome_proportion`, …) were scored as
  higher-is-better, so a bad sample could PASS. Direction is now handled per
  metric and overridable with `direction: lower|higher` in the thresholds YAML.
- **Bare `RiboMetric` crashed** with an AttributeError instead of printing help.
- **RUST KL divergence** used a non-standard per-term `abs()` and voided most
  window positions; it now computes a standard KL divergence over the observed
  support.
- Metagene zero-fill no longer drops the maximum read length.

### Changed

- Sequence-composition / terminal-bias backgrounds now use **all** subsampled
  reads by default instead of a silent 1-in-10 sample
  (`RIBOMETRIC_SEQUENCE_CHUNK_STRIDE` to opt back into sampling).
- Spectral metrics (Fourier, Theil, Gini) default to **all observed read
  lengths** rather than a hardcoded 25–35 nt human window; the di-some window is
  configurable via `qc.disome`.
- CI now runs on pull requests, not only pushes to protected branches.

### Removed

- Dead code: `multitaper`, `wavelet_transform` (used SciPy APIs removed in
  ≥1.13), `change_point_analysis`, `calculate_expected_dinucleotide_freqs`,
  and the unused `sequence_mode`.

---

## [1.2.0] — 2026-05-31

### Added

- **Alignment quality metrics** — calculated for every run, no flags required:
  - `duplicate_rate` — fraction of reads that are PCR/library duplicates, derived
    from the collapsed `count` column (`1 − unique_seqs / total_weighted`).
  - `multimapper_rate` — backwards-compatible alias for the weighted RPF-level
    multimapper rate.
  - `rpf_multimapper_rate` — fraction of weighted ribosome protected fragments
    with evidence of multiple genomic loci; uses `NH:i:N` when present and
    falls back to STAR-style `MAPQ < 255`.
  - `unique_rpf_rate` — fraction of weighted ribosome protected fragments that
    are uniquely mapped.
  - `alignment_multimapper_rate` — fraction of reported alignment rows whose
    read/fragment has evidence of another genomic alignment.
  - `soft_clip_rate_5prime` — fraction of reads with 5′ soft-clipping; an indicator
    of incomplete adapter trimming or read-through artefacts.
  - These are also surfaced in `results["alignment_stats"]` alongside the
    flagstat totals (`total_reads`, `mapped_reads`, `unmapped_reads`).
- **Di-some detection metrics** (default):
  - `disome_proportion` — fraction of reads in the 50–70 nt window (di-some /
    collision footprints).
  - `read_length_distribution_bimodality` — bimodality coefficient; moved from
    optional → default so di-some presence is always flagged.
- **Codon-level translation quality metrics** (default, annotation required):
  - `stop_codon_readthrough_ratio` — reads at positions +1…+30 relative to stop
    codon divided by reads at −30…−1; elevated values indicate stop-codon
    readthrough or frameshifting.
  - `start_codon_enrichment_ratio` — reads at positions −5…+20 relative to start
    codon divided by reads at +30…+50; very high values are a proxy for
    harringtonine/LTM inhibitor treatment or initiation stalling.
- **RUST metric** (optional, requires `--fasta`):
  - `rust_mean_kl_divergence` — mean KL divergence across the 60-codon RUST
    metagene window (O'Connor et al., *Nat Commun* 2016). Measures codon-specific
    A-site accumulation along the ribosome elongation region. Enable with
    `--enable-metric rust_mean_kl_divergence --fasta transcriptome.fa.gz`.
  - Full RUST output (`metagene`, `kl_divergence`, `transcripts_used`) stored
    under `results["rust"]` in the JSON.
- **`--fasta` CLI flag** — path to a transcriptome FASTA (gzip supported); passed
  to annotation mode for sequence-dependent metrics (RUST, future codon usage).
- **Gzip FASTA support** — `parse_fasta()` now transparently handles `.gz`/`.bgz`
  files; previously only uncompressed FASTA was accepted.

### Changed

- `read_length_distribution_bimodality` moved from the optional metric list to
  default. Use `--disable-metric read_length_distribution_bimodality` to skip.

---

## [1.1.0] — 2026-05-31

### Added

- **`ribometric evaluate` subcommand** — gates a pipeline on QC thresholds.
  Accepts a RiboMetric JSON output or a metrics-table CSV, plus an optional
  YAML of per-metric `pass`/`warn` thresholds (built-in defaults used when
  omitted). Exits 0 / 1 / 2 for PASS / WARNING / FAIL so it can be branched
  on in Nextflow / Snakemake / shell pipelines. Closes #124.
- **`--skip-sequence-metrics` flag** — skips terminal nucleotide bias and
  nucleotide composition calculations. Useful for BAMs with no stored sequences.
  `parse_bam` is also now robust to sequence-less BAMs (no ZeroDivisionError).
  Closes #122.
- **`--multimap-filter` CLI flag** — controls whether frame-sensitive
  calculations (P-site offset detection, periodicity) operate on uniquely-mapped
  reads only (`unique_only`) or all primary reads (`none`). Default: `unique_only`.
  Controlled by `multimap_filter` in `config.yml`.
- **`ribometric view` TUI subcommand** — interactive terminal viewer for
  exploring a RiboMetric JSON result.
- **Improved output flags** — `--summary-tsv`, `--qc-status`, `--comparison-csv`,
  `--metrics-table`, `--improved-outputs`, `--output-offsets`.
- **Docker image** — `ghcr.io/jackcurragh/ribometric` built and pushed on every
  tagged release via GitHub Actions, for `linux/amd64` and `linux/arm64`.
- **`pixi run test`** — the pixi environment now includes all test/dev
  dependencies so `pixi run test` (and `test-cov`, `lint`, `typecheck`) work
  without a separate virtualenv.
- **`ribometric` (lowercase) entry point** — alias for `RiboMetric`, so the
  command matches what you'd type naturally.

### Fixed

- **MAPQ multimapper bug** — Oxbow returns MAPQ=255 as `NaN` for uniquely-mapped
  reads (BAM-spec sentinel). The old `fillna(0)` silently treated every
  uniquely-mapped read as a multimapper. Fixed to `fillna(255)`.
  Effect: `periodicity_dominance_global` increases substantially on STAR
  transcriptome BAMs (example: 0.53 → 0.70).
- **Read length from CIGAR** — `process_reads()` now derives read length from
  CIGAR consumed-bases (M+=/X+I) rather than SEQ length, making it robust to
  soft-clipped reads.
- **Soft-clip NaN misalignment in A-site calc** — removed a `dropna()` that
  caused read-position misalignment; now uses `astype("float")`.
- **Shell injection** — `os.system("samtools index …")` replaced with
  `subprocess.run(["samtools", "index", …])`.
- **`pattern_to_index` invalid-base handling** — previously returned `0` for
  invalid bases (silently assigned them to the first bucket); now returns `-1`.
- **Count-weighted metrics** — reads collapsed by deduplication carry a `count`
  column; all frame and periodicity metrics now use these as weights to avoid
  double-counting.
- **`prepare` subcommand race condition** — per-transcript `append`-mode file
  writes inside a multiprocessing pool could interleave. Fixed by collecting
  DataFrames from workers and writing once at the end.
- Various `apply(lambda)` → vectorised operations in `metrics.py`.

### Changed

- **Version source of truth** — `__version__` in `RiboMetric/__init__.py` is
  now the single authoritative version; `setup.py` reads it dynamically.
- **Build system modernised** — `python -m build` (PEP 517) replaces the
  deprecated `python setup.py sdist bdist_wheel`. `pyproject.toml` added.
- **PyPI publish via Trusted Publishing** — new `release.yml` workflow publishes
  on every `v*` tag using OIDC (no stored API token). CI `ci.yml` continues
  to build and test on PRs/pushes.
- **Dependency caps added** — `numpy<2` and `pandas<2.3` reflect the
  tested-compatible range with `oxbow`/`pyarrow`; caps will be relaxed once
  validated against newer releases.
- **`Development Status` classifier** — changed from `Pre-Alpha` to `Beta`.
- Server mode removed entirely from the codebase (was non-functional).
- Dead `sequence_mode` commented block removed from `qc.py`.

### Notes for upgraders

- The PyPI package was `0.1.15`; this release jumps to `1.1.0` to align the
  version line with the existing git tags (`v1.0.0`/`v1.0.1` from 2023 are
  now superseded).
- The `multimap_filter` config key defaults to `unique_only`. STAR
  transcriptome-mapped BAMs will show higher periodicity scores than with
  `<= 0.1.15` where all primary reads were used. Set `--multimap-filter none`
  to reproduce pre-1.1.0 behaviour.

---

## [0.1.15] — 2025-12-10

- Fix dependencies: remove unused `click`, bump `rich` to `>=13.3.3`.

## [0.1.14] and earlier

See the [git log](https://github.com/JackCurragh/RiboMetric/commits/main) for
the full history of the `0.1.x` release series.
