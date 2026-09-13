# RiboMetric productionisation programme

**Status:** adopted 2026-09-13. Phase 0 not started.

The goal is to reach a point where RiboMetric is trusted enough that later
analyses with it are observations, not debugging exercises. This is a staged
programme with exit criteria, not an issue list: a phase is finished when its
exit criterion holds, not when its tasks are ticked off.

It is deliberately independent of any manuscript. Candidate manuscript
directions live in `../RiboMetric-Manuscript/PLAN.md` and are **not
requirements on RiboMetric**. Anything unexpected found during this programme
is logged in `../RiboMetric-Manuscript/OBSERVATIONS.md` before it is fixed or
interpreted. Decisions about metrics go in `../RiboMetric-Manuscript/decisions/`.

## The distinction this programme rests on

**Measurement correctness** is not **interpretation correctness**. We can prove
"frame-0 proportion = 0.81" was computed correctly before deciding whether 0.81
deserves PASS. Phases 0–4 establish measurement. Scoring thresholds are not
calibrated in this programme: empirical calibration (`METRICS_DESIGN.md`
Phase 2) comes after the measurement system is frozen.

## Where RiboMetric runs

Real-data runs go through the **ensembl-genes-nf Ribo-seq pipeline**
(`ensembl-genes-nf/pipelines/riboseq`), which is where RiboMetric is actually
used. It has two RiboMetric steps: `RIBOMETRIC_PREPARE` (organism setup: GTF →
annotation TSV) and `RIBOMETRIC` (per sample, on the STAR transcriptome BAM).

**riboseqorg-nf is legacy** and is not used by this programme. Its RiboMetric
module's defects (O-015) are out of scope.

ensembl-genes-nf is an Ensembl repository (REPO_CONTRACT tier C): changes to it
are proposed, not imposed.

## Phase 0 — Freeze the intended contract

Write down what RiboMetric is supposed to do before fixing anything else.

- **Every emitted number** (metrics, scores, `alignment_stats`, per-read-length
  tables, offsets, anything written to JSON/CSV/TSV) gets:
  `input → mathematical definition → expected range → interpretation →
  aggregation level → pathological cases`.
- **The pipeline**, stage by stage:
  `BAM + annotation → preparation → filtering → offset inference → assignment →
  raw measurements → scoring → evaluation/report`.
  Each stage states its inputs, outputs, defaults and what it discards.
- Output: `docs/METRIC_CONTRACT.md`, authoritative. `METRICS.md` stays generated
  from the registry; a later check should verify the registry against the
  contract. `METRICS_DESIGN.md` remains the interpretation/scoring design.
- Open questions to settle here: whether XA tags are supported (O-017); what
  "multimapper" means on a transcriptome BAM (O-005); one evaluation path,
  i.e. whether `evaluate`'s legacy threshold table (item D3) survives.

**Exit criterion:** every number emitted by RiboMetric has an explicit
definition.

## Phase 1 — Exhaustive correctness testing

This comes before any cohort run. Three test layers:

**Mathematical unit tests.** Every metric gets tiny inputs whose answer can be
computed by hand, and the test asserts that answer, not just "between 0 and 1".
For periodicity:

```text
perfect frame-0 periodicity      1/3 : 1/3 : 1/3
perfect frame-1 periodicity      flat coverage
perfect frame-2 periodicity      single spike
zero reads    one read    very short vector
```

Do the equivalent for every uniformity, bias, enrichment and read-length
statistic. Known target: the Fourier score selects the zero-frequency bin
(O-001). A perfectly periodic profile scores 0.40 and a flat one 1.00.

**Biological synthetic fixtures.** A tiny deterministic transcriptome and BAM,
used to test the whole pipeline rather than single functions:

```text
transcript A  perfect 28-nt periodic footprints
transcript B  perfect 29-nt footprints needing a different offset
transcript C  uniform RNA-like coverage
transcript D  no CDS reads
transcript E  multimappers
transcript F  secondary / supplementary alignments
```

This is where known per-length P-site offsets are asserted for every offset
method, including the ±1 frame calibration at the offset bounds (AUDIT_NOTES
§1.2–1.3). Transcript F covers alignment-tag recovery by row order (O-016).

**Golden integration tests.** A few very small fixtures with their expected
RiboMetric JSON committed. CI fails if a change alters those outputs. An
intentional change updates the golden file, and the commit says why.

**Exit criterion:** every metric and every major processing decision has a
focused unit test, plus an end-to-end fixture where appropriate.

## Phase 2 — Real-data validation passes

Not the whole Portal. A deliberately selected **validation panel** of 30–50
libraries:

- textbook high-quality monosome libraries;
- weak periodicity;
- extreme terminal bias;
- unusual length distributions;
- initiation profiling;
- disome / long-footprint profiling;
- very shallow and very deep libraries;
- several organisms and annotations;
- datasets already known well enough to spot nonsense;
- a few matched samples within studies.

Candidate members can be drawn from the classes in O-014.

**Running the panel.** Runs go through the ensembl-genes-nf Ribo-seq pipeline.
Before pass 1, that pipeline needs three changes (O-022):

- Both RiboMetric steps use the **image under test**: one image, referenced by
  digest, for `RIBOMETRIC_PREPARE` and `RIBOMETRIC`. Today `prepare` pulls the
  legacy, unpinned `riboseqorg-nf-ribometric:latest`, and `run` is pinned to
  1.4.3.
- Failed samples are reported, not dropped. `ANALYSIS:RIBOMETRIC` currently
  retries once and then ignores the failure.
- Each pass records the image digest, the pipeline commit and the
  `ribometric_*` parameters (offset method, offset target, sample size).

The pipeline's own `validation/VALIDATION_MATRIX.md` compares processing
choices (adapters, rRNA removal, alignment, multimappers). It is separate from
this programme: panel runs hold those settings fixed and record them.

**Inspecting the panel.** Inspect every report by hand. For each sample record:

```text
Expected characteristics
Observed characteristics
Metrics that make sense
Metrics that look surprising
Report/visualisation problems
Potential bugs
Interpretation uncertainty
```

Don't change a metric because one sample looks odd. Open an observation,
reproduce it, investigate, then decide.

- **Pass 1:** the panel. Fix what it reveals.
- **Pass 2:** the same panel, rerun exactly.
- **Pass 3:** expand to 100–200 diverse libraries.

Where panel records live is still open; proposed:
`../RiboMetric-Manuscript/validation/`, created when Phase 2 starts.

**Exit criterion:** pass 3 completes with every surprising result either
explained or recorded as an open observation. Only then is a Portal-scale run
trusted.

## Phase 3 — Stress and invariance testing

This tries to break the software. It is not the manuscript's perturbation
experiment.

- **Should not materially change results:** downsampling (within stated
  tolerance); BAM ordering; presence of secondary alignments; annotation order;
  transcript order; compression; thread count; repeated seeded runs; read names;
  irrelevant extra annotation; multimapper handling where the metric is
  defined not to depend on it.
- **Should change results, in the expected direction:** wrong P-site offsets;
  ±1 annotation shifts; removed CDSs; frame shifts; removal of particular read
  lengths.

Historical bugs become regression tests here. Annotation handling gets
disproportionate attention, because coordinate and frame errors have produced
plausible-looking but wrong periodicity before (O-002, O-003). Subsampling must
be seeded and the seed recorded.

**Exit criterion:** we know which outputs are invariant to irrelevant
implementation details and which are sensitive to biologically meaningful
changes, and both are encoded as tests.

## Phase 4 — Metric audit

Only now ask whether each metric deserves to survive:

- Is it correct?
- Does it measure something distinguishable from another metric?
- Can we explain what high and low values mean?
- Does it behave sensibly across real datasets?
- Does a user need it?

Each metric ends up in one of three statuses:

- **Primary:** exposed prominently; eventually eligible for scoring.
- **Diagnostic:** useful for understanding a sample, but not part of the
  headline assessment.
- **Deprecated / experimental:** kept temporarily for compatibility, or
  removed.

Eight defensible measurements beat twenty partially redundant ones. Thesis-era
estimators (the several periodicity, uniformity, terminal-bias and read-length
variants) have to earn their place; they don't keep it for having existed. The
corrected Fourier metric is compared properly here, not retained just because
it was fixed. Each outcome is recorded in `../RiboMetric-Manuscript/decisions/`.

**Exit criterion:** every emitted metric has a status and a decision record.

## Phase 5 — Report and UX audit

Once the numbers are trusted, use RiboMetric as if it were someone else's tool.
Given a BAM, can a user tell:

- what kind of library they appear to have;
- what looks good and what looks problematic;
- which read lengths to use;
- what evidence produced each conclusion;
- what they should *not* conclude from the report?

Fitness-for-purpose presentation is revisited here, not before.

Test bad input deliberately: wrong BAM type (genomic vs transcriptomic), empty
BAM, missing annotation, mismatched chromosome/transcript names, zero CDS reads,
malformed GFF/GTF, unsupported tags, no usable read lengths, absurdly deep BAM.

**Exit criterion:** every bad-input case produces an informative error rather
than a stack trace, and the report answers the questions above on the
validation panel.

## Phase 6 — Reproducibility and production engineering

- **Provenance block in every result JSON:** tool version and git SHA, command
  line, input names/sizes (optional checksum), annotation path and hash,
  subsample N and seed, offset method, filters, platform, timestamp. Currently
  `generate_json` writes only `{results, config}`.
- CI across supported Python versions, lint and type checking, installation
  from a clean environment, container build, CLI smoke tests, deterministic
  fixtures.
- Runtime and memory benchmarks, covering `prepare` memory (AUDIT_NOTES §2.2)
  and the extra pysam decode pass (§2.1).
- Migration notes for the renamed 2.0 metrics, including the one reused key
  (`terminal_bias_kl_5prime` now holds KL in bits).
- **ensembl-genes-nf integration:** propose that both RiboMetric steps take
  their image from one parameter, so `prepare` and `run` can't drift apart
  (O-022). Also check the module's `conda` directive, which doesn't list
  RiboMetric.

**Exit criterion:** a result JSON alone identifies exactly how it was
produced, and a clean install from a built package or container reproduces
the golden fixtures.

## Completion

The programme is complete when **two consecutive real-data validation passes
produce no correctness change to any metric definition**. The passes run on a
built package or container, not the development checkout. Visual and report
fixes don't reset the count; a changed definition (a periodicity formula, say)
does.

At that point the measurement system is frozen, and calibration and manuscript
work can follow.

## Known items, by phase

| item | source | phase |
|---|---|---|
| Fourier score selects DC bin | O-001 | 1 (test), 4 (survival) |
| No tests for offset values | AUDIT_NOTES §1.2–1.3 | 1 (fixtures A/B) |
| Alignment tags recovered by row order | O-016, AUDIT_NOTES §1.4 | 1 (fixture F) |
| XA tag no-op | O-017 | 0 |
| Transcriptome multimapper semantics | O-005 | 0, 5 |
| `evaluate` legacy threshold path (D3) | deferred scoring item | 0, 1 |
| Annotation coordinate handling | O-002, O-003 | 3 |
| CDS enrichment on strongly periodic libraries | O-021 | 2 |
| Aggregate vs per-length periodicity | O-013 | 2 |
| ensembl-genes-nf: legacy `prepare` image, 1.4.3 pin, ignored failures | O-022 | 2 (before pass 1), 6 |
| Provenance, seeded subsampling | — | 3, 6 |
| `prepare` memory, pysam double decode | AUDIT_NOTES §2.1–2.2 | 6 |
| Dead code, duplicated uniqueness logic | AUDIT_NOTES §3 | 4, 6 |
| riboseqorg-nf module | O-015 | out of scope (legacy) |
