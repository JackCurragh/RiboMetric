# RiboMetric Roadmap

Status as of the v1.1.0 release effort.

> **See also:** [`docs/AUDIT_NOTES.md`](docs/AUDIT_NOTES.md) — release-readiness
> audit (2026-06-09): what remains unscrutinised, optimisation targets, and
> code-cleanliness debt from the annotation-coordinate / offset-calibration work.

---

## Resolved (shipping in v1.1.0)

- **MAPQ multimapper fix** — Oxbow returns MAPQ=255 as `NaN`; the old `fillna(0)`
  zeroed uniquely-mapped reads. Now `fillna(255)`.
- **Multimapper filter** — frame-sensitive metrics operate on MAPQ=255 reads by
  default (`multimap_filter: unique_only`, `--multimap-filter`).
- **Read length from CIGAR** — `process_reads()` derives length from CIGAR
  (M+I+=/X), robust to soft-clipping.
- **Count-weighted metrics** — collapsed reads carry a `count` weight.
- **`prepare` performance** — vectorised `str.extract`, numpy-only CDS building,
  removed the multiprocessing file-write race.
- **Server mode removed** — `parse_bam` no longer carries a server path.
- **kaleido** — pinned `>=1.0.0` (modern wheels).
- **Security/correctness** — shell injection (`os.system` → `subprocess.run`),
  `pattern_to_index` invalid-base handling, soft-clip NaN A-site fix.

---

## In progress for v1.1.0

- Optional sequence-based metrics for sequence-less BAMs (#122).
- `RiboMetric evaluate` subcommand: gate a result against expected metrics (#124).
- Packaging: single-source version, PyPI Trusted Publishing, CHANGELOG.
- JOSS submission materials (`paper.md`, `CITATION.cff`, Zenodo archival).

---

## Future (v1.2+)

- **pyranges1 GFF parsing** — replace `gffpandas` + `str.extract` with
  `pr.read_gff3()` (Rust) for a further `prepare` speedup. Scope: `file_parser.py`,
  `setup.py`.
- **Metric-name registry** — `should_calculate_metric` (`qc.py`) still matches
  config keys via string variants. Replace with an explicit metric-ID → config-key
  registry.
- **Dependency caps** — decide on upper bounds for `pandas`/`scipy` once usage in
  production is understood.
