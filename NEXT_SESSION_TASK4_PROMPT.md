# Next-session launch prompt — Task 4

Continue productionising RiboMetric in `/Users/jackt/projects/all-RiboSeq/RiboMetric`.

## Objective

Complete the release-candidate and Ensembl-integration checkpoint without discarding any existing uncommitted work. Validate the contract, package, and documentation artifacts; inspect the full diff and working tree; and prepare a distinguishable release candidate locally only. Do not commit, push, tag, publish, build/push a container, or delete branches unless explicitly asked.

## Current implementation state

Task 1 annotation/QC work and Task 3 reproducibility/audit work are present in the working tree. `-S/--subsample` now uses a stable SHA-256 ordering of read-name groups with default seed 42, records requested seed/fraction/realised count under `provenance.subsampling`, and recomputes sequence summaries from selected rows. FLOSS is called on CDS-body reads and applies count weights to both transcript and reference distributions. Offset provenance records requested/effective methods, external-supply status, raw/final mappings, and frame adjustments. Terminal-bias backgrounds are count-weighted, and weighted-equivalence tests cover collapsed-read weighting.

Direct legacy callers that construct a Namespace without the new `seed` field retain old fixture behavior for compatibility; normal CLI parsing always supplies seed 42. Revisit this only if the public API contract requires direct Namespace callers to receive the new sampling semantics.

## Files changed or inspected

Task 3 additionally touched `RiboMetric/RiboMetric.py`, `RiboMetric/arg_parser.py`, `RiboMetric/config.yml`, `RiboMetric/bam_processing.py`, `RiboMetric/file_parser.py`, `RiboMetric/modules.py`, `RiboMetric/qc.py`, `tests/test_audit_fixes.py`, `tests/golden/synthetic_tripsviz_a_site.json`, and `PRODUCTIONISATION_TASKS.md`. Other pre-existing modified files include `CHANGELOG.md`, `metrics.py`, `results_output.py`, `scoring.py`, `docs/METRIC_CONTRACT.md`, `docs/PRODUCTIONISATION.md`, and related tests/fixtures. Preserve all unrelated edits.

## Evidence completed

- `./.test-venv/bin/pytest -q`: 409 passed, Python 3.12.
- `./.test-venv310/bin/pytest -q`: 409 passed, Python 3.10.
- `./.test-venv/bin/ruff check RiboMetric tests`: passed.
- `./.test-venv/bin/black --check RiboMetric tests`: passed.
- `./.test-venv/bin/python -m sphinx -b html -W --keep-going docs /tmp/ribometric-sphinx-check`: passed.
- `git diff --check`: passed.

The synthetic golden output was intentionally regenerated because weighted terminal backgrounds and CDS-body FLOSS change contract values.

## Remaining tasks and acceptance criteria

1. Inspect `git diff` and `git status`; confirm no accidental generated files or unrelated reversions.
2. Reconcile the checklist, contract docs, productionisation notes, and changelog with only behavior supported by tests.
3. Investigate the package-build limitation recorded in the checklist (`setuptools.build_meta` unavailable in `.test-venv`) using a safe environment check. Do not mutate the release version to solve it.
4. Verify the Ensembl repository state read-only. `riboseq_unique_reads` already contains commit `5ac01ab` via `4f7fd45`; do not duplicate or cherry-pick it.
5. Update the checklist only when each item has direct command/test evidence.
6. Leave a final concise handoff listing exact commands, results, files, and unresolved limitations. No publication or git-history mutation is part of this task.

## Exact commands/tests

```bash
cd /Users/jackt/projects/all-RiboSeq/RiboMetric
git status --short
git diff --check
./.test-venv/bin/pytest -q
./.test-venv310/bin/pytest -q
./.test-venv/bin/ruff check RiboMetric tests
./.test-venv/bin/black --check RiboMetric tests
./.test-venv/bin/python -m sphinx -b html -W --keep-going docs /tmp/ribometric-sphinx-check
```

## Non-goals and traps

- Do not reset, commit, push, tag, publish, delete branches, or build/push a Docker image.
- Do not bump the release version.
- Do not modify Ensembl or duplicate the already-present integration.
- Do not treat externally supplied offsets as calibrated.
- Do not reintroduce BAM-order prefix subsampling, unweighted backgrounds, or double application of collapsed-read counts.
