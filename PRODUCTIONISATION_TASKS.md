# RiboMetric productionisation task list

This checklist tracks the four work packages from the handoff supplied with this task. Items are marked complete only when the implementation and verification evidence support them.

## Task 1 — Annotation semantics

- [x] Support `cds_start == 0` and explicit `cds_start_complete` evidence.
- [x] Handle complete/incomplete CDSs, GTF/GFF3 semantics, and GFF3 stop-codon normalization.
- [x] Preserve all observed read lengths in frame tables; keep plot limits presentation-only.
- [x] Add and pass local fixtures for zero-start/incomplete CDS evidence, transcript coordinates, and malformed/ambiguous records.
- [x] Add dedicated GTF/GFF3 stop-codon fixture coverage for included/excluded stop-codon variants.
- [x] Update contract and productionisation notes where test-supported.

## Task 2 — Missing values and QC

- [x] Emit JSON `null` for uncomputable metrics, zero denominators, and specified empty inputs across measured producers.
- [x] Make null/absent gated metrics fail annotation-mode QC; align `evaluate` with `qc_status.json`.
- [x] Preserve `duplicate_rate` distinction between uncollapsed (`null`) and genuine zero (`0`).
- [x] Add regression coverage for empty inputs, missing evidence, null gates, and collapsed/uncollapsed BAMs.

## Task 3 — Reproducibility and auditability

- [x] Make `-S` seeded/reproducible and record seed, fraction, and realised count.
- [x] Calculate FLOSS from CDS-body reads for both reference and sample distributions.
- [x] Record effective offset execution, external/calibrated status, and terminal-bias weighting.
- [x] Verify collapsed-read weighting is applied exactly once and add deterministic provenance tests.

## Task 4 — Release candidate and Ensembl integration

- [x] Validate registry tests and strict documentation build.
- [x] Complete contract/implementation/golden-fixture validation.
- [x] Run Python 3.10/3.12 tests, strict docs, and Docker CLI smoke test.
- [x] Run package build/install smoke tests in a clean Python 3.12 environment.
- [x] Prepare a distinguishable local 1.5.0 release candidate without publishing/tagging.
- [x] Verify the existing Ensembl integration after validation; no cherry-pick or branch mutation is needed. Keep one digest-pinned image for `prepare` and `run`.

## Evidence log

- Initial state: existing uncommitted implementation in `RiboMetric/` was preserved.
- Task 2 evidence: focused annotation/module tests and full regression coverage pass, including missing evidence, null gates, duplicate-rate semantics, and unavailable mapping signals.
- Task 1 evidence: annotation and frame-distribution regression tests pass, including zero-start CDSs and preservation of observed read lengths.
- Task 2 evidence: default annotation-mode null-gate and duplicate-rate regressions pass; explicit `evaluate --expected` missing/non-finite checks pass.
- Full-suite evidence: 410 tests pass on `.test-venv` (Python 3.12) and `.test-venv310` (Python 3.10); Ruff and Black pass; strict Sphinx build passes.
- Ensembl evidence: `riboseq_unique_reads` is clean at `4f7fd45` and contains the same one-image integration as `5ac01ab` on `feat/ribometric-image-under-test`; no branch mutation was needed.
- Task 3 evidence: deterministic read-name subsampling, weighted terminal backgrounds, CDS-body FLOSS, offset audit provenance, and weighted-equivalence regressions pass in the focused and full Python 3.12 suites.
- Release-candidate evidence: clean Python 3.12 build produced `ribometric-1.5.0.tar.gz` and `ribometric-1.5.0-py3-none-any.whl`; `twine check`, wheel install, `RiboMetric --version`, `RiboMetric --help`, and template packaging all passed. No Docker image was built or pushed.
