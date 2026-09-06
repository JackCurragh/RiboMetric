## Summary

Fix three reproducible correctness failures in terminal-bias reporting and downstream QC/export. Base: `b2a0f4aedff8d15ec8081af0a1ef5701b1c64ac4`.

### 1. Zero-background terminal KL can appear perfect

For observed `{AA: 1, CC: 0}` and expected `{AA: 0, CC: 1}`, the current raw function returns 0 and the goodness score returns 1. The raw divergence should be positive infinity; the existing transform should yield 0.

The fix preserves zero-observation contributions, detects unsupported positive mass, validates finite probability values and avoids overflow in a probability ratio. It does not add smoothing or claim an empirical background zero is biologically impossible. Report raw values use JSON null plus the explicit status `zero_background_support` for this case. Score aliases remain numeric and finite.

### 2. Explicit evaluation policies can pass without required evidence

`evaluate_qc_status({"metrics": {}}, "sample", {"periodicity_dominance": {"pass": 0.7, "warn": 0.5}})` returns PASS with no checks. The evaluate command handler consequently returns 0. Null/string metrics, a misspelled required metric, partially missing policies and nonnumeric CSV values can produce the same false-PASS outcome.

An explicit policy now evaluates every required metric. Missing, nonnumeric and unexplained nonfinite values produce FAIL with a diagnostic reason. Empty/malformed policies are rejected; the command handler catches those errors and returns 2. Finite values retain pass/warn ordering and directionality. `thresholds=None` still delegates to the existing scored resolver. A known infinite raw KL represented by its JSON status is evaluated rather than silently skipped.

A FAIL caused by unavailable evidence is documented as incomplete QC evidence, not proof that the biological sample is poor.

### 3. Appending summary rows can swap metric labels

Write a first sample with metrics ordered `periodicity_dominance, uniformity_entropy`. Append a second sample with the order reversed and values `{uniformity_entropy: 0.21, periodicity_dominance: 0.74}`. The existing header labels the second row as periodicity 0.21 / uniformity 0.74.

The writer now honors the existing header order, retains null columns and leaves blanks for absent existing metrics. New columns or malformed headers are rejected before writing; use a fresh file or the comparison-CSV writer for schema expansion. Concurrent multi-process appending remains unsupported.

## Validation supplied with this contribution

104 focused regression cases passed in the authoring environment; 16 submission/patch safety tests passed. A further 500 seeded finite-value comparisons run within one regression case. These are not 500 additional pytest cases.

The authoring harness used the **complete, Git-blob-verified `evaluate.py`** and executable excerpts of the upstream metric, reporting and QC assignment code. It did not install the full upstream package or process a BAM. Tests cover actual JSON/CSV/YAML I/O, the command handler (including a subprocess exit test), KL-to-JSON-to-gate behavior, finite controls, missing inputs, malformed policies and column identity.

Before publication, the supplied submitter requires a fresh full checkout, exact commit/blob verification, executable excerpt/whole-source AST equivalence, failing unpatched controls and a passing complete focused suite against the real source files. The upstream CI workflow runs this focused suite on Python 3.11/3.12/3.13. CI has not been run by the authoring environment.

## Compatibility and review

- Explicit thresholds become a required-check contract, rather than best-effort checks. Users relying on silent omission must change their policy.
- JSON adds raw-KL status keys. Raw null is not zero divergence; downstream consumers should read the status or use the bounded score.
- Numeric-only CSV writers may omit null raw values/status strings; explicitly requiring an omitted raw value now fails safely rather than passing. JSON preserves the status.
- Existing summary columns are retained; appending a newly introduced metric now raises instead of writing an inconsistent row.
- The automatic scored resolver and metric gate memberships are unchanged. Adding a diagnostic does not automatically make it a Tier-1 gate.
- Full BAM/FASTA workflows, all other report consumers, every supported dependency version and original publication datasets have not been revalidated. No claim about published biological conclusions is made.

Reference semantics: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rel_entr.html
