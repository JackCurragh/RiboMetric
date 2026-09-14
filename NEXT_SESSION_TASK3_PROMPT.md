# Next-session launch prompt — Task 3

Continue productionising RiboMetric in `/Users/jackt/projects/all-RiboSeq/RiboMetric`.

The previous session completed and verified the annotation/QC checkpoint. Preserve all existing uncommitted work; do not reset, commit, push, tag, publish, or delete branches. The current working tree intentionally contains prior implementation changes and tests.

## Objective

Implement the remaining reproducibility and audit semantics:

- D-4: make `-S/--subsample` a seeded, reproducible sample and record seed, sampling fraction, and realised count;
- D-5: calculate FLOSS using CDS-body reads for both reference and sample distributions;
- D-8: record what offset calculation actually ran, including external offsets and frame calibration;
- D-9: weight terminal-bias backgrounds consistently with observed reads;
- verify collapsed-read weighting is applied exactly once;
- add tests for repeated names, count weights, external offsets, calibrated offsets, fixed offsets, and repeated runs with the same seed.

Extend provenance only as needed so a result can identify package version, command, configuration hashes, input hashes, subsampling seed/fraction/realised count, requested/effective offset method, and whether offsets were externally supplied or calibrated.

## Current implementation state

Annotation and QC changes are present in the working tree. `cds_start == 0` and explicit `cds_start_complete` are supported; observed read lengths are preserved; annotation-mode null/absent gated evidence fails; uncollapsed `duplicate_rate` is `null` while genuine collapsed zero remains `0`. Existing count weighting in `modules.py` was changed to avoid dropping weights when read names repeat. FLOSS currently uses annotated reads and count weights, but its exact CDS-body/reference semantics still need verification and tests. Offset provenance is partially present in `qc.py` and needs audit against actual execution paths. Subsampling is currently capped by BAM parsing order and has no complete seeded-sampling/provenance implementation.

## Invariants and decision rules

- Preserve legacy output keys and compatibility where the contract explicitly requires it.
- Do not apply collapsed-read weights twice; one row represents `count` reads unless the code explicitly proves otherwise.
- FLOSS reference and per-transcript distributions must use the same CDS-body eligibility rules.
- A plot/read-length display limit must never change numerical metrics.
- External offsets must not be described as calibrated; record requested and effective methods separately.
- Do not bump the release version, build/push a container, or modify Ensembl in this task.

## Evidence already completed

- `.test-venv/bin/pytest -q`: 404 passed (Python 3.12).
- `.test-venv310/bin/pytest -q`: 404 passed (Python 3.10).
- `.test-venv/bin/ruff check RiboMetric tests`: passed.
- `.test-venv/bin/black --check RiboMetric tests`: passed.
- strict Sphinx build: passed.
- Ensembl `riboseq_unique_reads` already contains commit `5ac01ab` via `4f7fd45`; do not duplicate that integration.

## Files touched or already modified

Existing modified files include `RiboMetric/bam_processing.py`, `RiboMetric/file_parser.py`, `RiboMetric/metrics.py`, `RiboMetric/modules.py`, `RiboMetric/qc.py`, `RiboMetric/results_output.py`, `RiboMetric/scoring.py`, `docs/METRIC_CONTRACT.md`, `docs/PRODUCTIONISATION.md`, `CHANGELOG.md`, and related tests/fixtures. The checklist is [PRODUCTIONISATION_TASKS.md](PRODUCTIONISATION_TASKS.md). Do not overwrite unrelated changes.

## Required commands

Run focused tests while developing, then:

```bash
./.test-venv/bin/pytest -q
./.test-venv310/bin/pytest -q
./.test-venv/bin/ruff check RiboMetric tests
./.test-venv/bin/black --check RiboMetric tests
./.test-venv/bin/python -m sphinx -b html -W --keep-going docs /tmp/ribometric-sphinx-check
```

Before finishing, inspect `git diff` and `git status`, update `PRODUCTIONISATION_TASKS.md` only for items genuinely supported by tests, and write the complete Task 4 continuation prompt including objective, implementation state, files changed, tests/results, unresolved decisions, exact commands, and non-goals. Do not merely say “continue”.
