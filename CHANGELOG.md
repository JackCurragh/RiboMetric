# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Fixed

- **QC gating could pass with no evidence** — in the explicit-thresholds path
  (`RiboMetric evaluate --expected`), a metric that was absent from the results
  or non-numeric was silently skipped. A policy naming only metrics the run
  never produced returned `PASS` with zero checks performed and exit code 0, so
  a single typo in a thresholds YAML green-lit everything. This is reachable
  from an ordinary cohort: `qc.py` only populates the terminal nucleotide bias
  metrics when a sequence background is available, so a run without `--fasta`
  omits twelve metric keys.
- **Summary TSV could label values with the wrong metric** — `generate_summary_tsv`
  appended using the current row's key order without reading the existing
  header, so a sample whose metric set differed from the first sample's had its
  values written under the previous sample's column names.
- **Type and lint errors that shipped in v1.4.3** — an under-annotated record
  dict in `qc.py`, an `Any` return from `_config_path_used`, and an unused
  `cds_start` read in `file_parser.py`. These reached PyPI because `ci.yml`
  watches `main`/`dev` while releases were cut from `release/v1.3.0`, and
  `release.yml` publishes on tag push with no test dependency.

### Changed

- **An explicit threshold policy is now a required-check contract.** Every
  metric it names must be present and finite; otherwise the check fails with a
  diagnostic `reason`. Empty and malformed policies are rejected outright, and
  the `evaluate` command handler reports those as errors rather than as a
  verdict on the sample. A `FAIL` caused by an absent metric means incomplete QC
  evidence, not a poor biological sample. Callers that relied on missing metrics
  being skipped must trim their policy to the metrics the run actually produces.
- **QC status checks carry a `reason` field**, `None` when the metric was
  evaluated normally.
- **Appending a new metric column to an existing summary TSV now raises**
  instead of writing a row inconsistent with the header. Use a fresh file, or
  the comparison CSV writer, to widen the schema. Concurrent multi-process
  appending remains unsupported.

### Infrastructure

Brings the repo in line with `REPO_CONTRACT.md`, the shared maintenance
contract for the all-RiboSeq tools. None of this changes RiboMetric's
behaviour; it changes what can silently go wrong when releasing it.

- **The distribution name is lowercase (`ribometric`).** It was `RiboMetric`,
  and published only because `pyproject.toml` leaves `setuptools>=64` unpinned
  so CI happens to resolve a version that normalises the wheel filename itself.
  Pin setuptools, or build on an older one, and PyPI rejects the upload with a
  400 that `twine check` does not catch and that reads like an auth failure.
  The import package is still `RiboMetric` and both console scripts are
  unchanged.
- **Lint is enforced.** flake8 was commented out in `ci.yml` while the Makefile
  and pixi both invoked it, so nothing ran. The toolchain is now ruff (lint and
  import order, replacing flake8 *and* isort) plus black at line-length 100,
  both fatal in CI. The findings this exposed -- unused imports, mid-file
  imports, a duplicated import in the test suite -- are fixed, and the codebase
  is black-formatted.
- **`.bumpversion.cfg` replaces hand-editing.** bump2version now rewrites
  `RiboMetric/__init__.py`, `pixi.toml`, `CITATION.cff` and the conda recipe in
  one commit, and tags `vX.Y.Z`. `conda-ribometric/meta.yaml` had drifted to
  0.1.9 while PyPI shipped 1.4.3 because nothing kept them together.
- **`release.yml` calls `ci.yml` instead of carrying its own copy of the
  build.** A tagged release now re-runs the full gate -- lint, mypy, tests on
  3.10 and 3.12, build -- and publishes the artifact the gate built, so the
  shipped bits are the checked bits. It also creates a GitHub Release, which
  eleven tags so far have not.
- **CI asserts the artefact filenames**, and `make preflight` does the same
  locally.
- **One owner per container tag.** `build-container.yml` fired on `tags: '*'`
  with no test gate and wrote `:latest` on every push to `main`. `:latest` now
  means the most recent release and is written only by `release.yml`; `:main`
  and `:sha-<short>` come from `ci.yml`, gated on lint and tests.
- **The dead `lukasdev` CI trigger is removed** (the branch is
  `lukasdev-archive` now).
- **`make release` (a bare `twine upload`) is removed** -- it was a second
  publish path that ran no checks. Pushing the tag is the only one.

### Documentation

- **The Read the Docs build was broken.** `.readthedocs.yaml` used
  `python.version`, a key removed from config v2 years ago, and had no `build:`
  block — which is required, so the build failed before Sphinx ran. It now
  pins `ubuntu-24.04` and Python 3.12.
- **`docs/contributing.rst` and `docs/history.rst` included files that do not
  exist** (`CONTRIBUTING.rst`, `HISTORY.rst`), and both are in the toctree.
  `CONTRIBUTING.rst` is now written; `docs/history.rst` points at
  `CHANGELOG.md`, which is where the history actually lives.
- **The metrics documentation was unreachable.** `METRICS.md`,
  `METRICS_DESIGN.md` and `REPORTING_GUIDE.md` were not in any toctree, and
  `source_suffix` was `.rst` only, so Sphinx could not have rendered them
  regardless. Added myst-parser and put them in the toctree, along with
  `results`, `functions` and the new `RELEASE.md`. Internal working notes
  (`AUDIT_NOTES.md`, `SCORING_PHASE1_TASKS.md`, `TESTING.md`,
  `V1.2_REPORT_AND_PREPRINT_PLAN.md`) are explicitly excluded rather than
  published.
- **The API reference never built on Read the Docs.** `modules` is produced by
  `sphinx-apidoc`, which `make docs` ran locally but RTD did not, so the
  toctree entry dangled on the published site. `docs/conf.py` now generates it
  from a `builder-inited` hook, so local and RTD builds take the same path.
- **`make docs` could not run**, because it delegated to `docs/Makefile`, which
  calls a bare `python` — the same "command not found" failure the top-level
  `$(PYTHON)` convention exists to prevent. It now invokes Sphinx directly.
  Added `make docs-strict` for a warnings-fatal build.
- **The docs build is now warning-free and warnings are fatal**, in Read the
  Docs (`fail_on_warning: true`), in CI (a new `docs` job) and locally (`make
  docs-strict`). It went 43 -> 25 -> 0. The first pass fixed configuration and
  doc sources; the second fixed the docstrings behind the remaining 25:
  - Six functions wrote `Inputs:`/`Outputs:` with the body at the *same* indent
    as the label, which docutils reads as a block quote. 132 other docstrings
    were already indented correctly, so these were the outliers, not the norm.
  - The singular `Input:`/`Output:` spelling, used throughout
    `results_output.py` and `plots.py`, was not registered with napoleon, so 23
    more docstrings rendered as an undifferentiated paragraph rather than as
    parameter fields. All four spellings are now registered.
  - Three functions that return a multi-key dict documented the keys as a
    column-aligned block under `Returns:`. Google-style `Returns:` takes a
    single `type: description`, so napoleon read the first key as the return
    *type* and orphaned the rest. They now read `dict: ...` with the keys as a
    bullet list.
  - `autocorrelate_counts` mixed NumPy-style underlines with Google-style
    colons; it now matches the rest of its module. (`rust.py`'s genuine
    NumPy-style docstrings are untouched -- they were always valid.)
  - The `qc` and `rust` module docstrings had prose and an algorithm listing
    indented under a trailing colon, which is a block quote in RST. They are
    now a definition list and an enumerated list.
  - Six of the 25 were not RiboMetric's at all: `HelpScreen.compose` has no
    docstring, so autodoc inherited **textual's**, which is MkDocs-flavoured
    Markdown and cannot parse as RST. `autodoc_inherit_docstrings` is off.

  No executable code changed -- verified by checking every modified line in the
  package falls inside a docstring node.
- **`docs/RELEASE.md`** is new: the release runbook, the container-tag
  ownership split, and the traps that have actually cost this project releases
  (wheel-filename normalisation, 403-vs-400 diagnosis, the reusable-workflow
  permission startup failure, why there is no TestPyPI rehearsal).

### Packaging

- **The sdist shrank from 33 MB to 184 KB.** `MANIFEST.in` had
  `recursive-include tests *`, which shipped 47 MB of BAM/BAI fixtures for a
  package whose own code is 0.4 MB. Worse, it also swept in whatever test
  *outputs* were in the working tree — `tests/test_data/test_RiboMetric.html`
  is gitignored and was being published anyway — so the sdist varied depending
  on whether the suite had been run before building. Nothing consumed the
  shipped tests: the conda recipe checks imports and `--help` only, and CI runs
  the suite from a git checkout. Verified by installing the sdist into a clean
  environment: both entry points, `config.yml` and the report templates are all
  present.
- Test-run outputs written into `tests/test_data/` are now gitignored by name.
  They are ignored individually rather than by extension because that directory
  also holds tracked `.tsv`/`.json`/`.csv` *inputs*.

### Notes

- `main` now contains the v1.4.1-v1.4.3 releases, which had been cut from
  `release/v1.3.0` and never merged back. `dev` is the integration branch.

## [1.4.3] — 2026-08-19

### Fixed

- **Canonical final offsets** — frame-calibrated per-read-length offsets are now
  consistently used by annotation-dependent metrics and downstream profiles,
  and are reported in JSON and TSV audit outputs.
- **Offset provenance** — JSON now distinguishes raw estimates from final
  offsets; applied offsets no longer come from a stale pre-calibration dataframe.
- **External offset semantics** — externally supplied offsets are explicitly
  marked as final offsets because internal frame calibration is bypassed.

## [1.4.2] — 2026-07-10

### Fixed

- **`--output-offsets` CLI crash** — removed a function-local `Path` import that
  shadowed the module-level import and caused `UnboundLocalError` on normal
  `RiboMetric run --output-offsets ...` executions.

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
