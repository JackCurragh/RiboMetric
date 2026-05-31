# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
