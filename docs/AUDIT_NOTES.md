# Audit & Optimisation Notes

Working notes for a release-readiness pass. Captures what has **not** yet been
properly scrutinised, where the code can be **optimised**, and **cleanliness**
debt. Written 2026-06-09 against `main` (baseline commit `79c79b2`) with a large
in-progress, uncommitted feature set in the tree (offset-target / frame
calibration, NH/XA-aware alignment stats, MAPQ-sentinel recovery).

The goal is that a fresh session can pick any item up cold. Each item lists the
file(s)/symbol and concrete next step.

---

## Context: what was already done in this audit

- **Root cause of "poor periodicity"** found and fixed: the shipped
  `sample_data/gencode.v25.annotation_RiboMetric.tsv` had CDS coordinates in
  genomic/unspliced space, not transcript-relative. Regenerated correctly.
  Validated end-to-end on a real BAM (`SRR11005875`, v49) → dominance 0.79.
- **Guard added**: `validate_annotation_coordinates` in `file_parser.py` rejects
  annotations where `cds_end > transcript_length` for >1% of rows.
- **Crash fixed**: `plot_read_frame_distribution` (`plots.py`) on an empty/all-
  zero frame dict (`UnboundLocalError`).
- **Two hot loops vectorised** (see Optimisation §"Done").
- Full suite green: **195 passing**.

---

## 1. Not yet properly scrutinised

Ordered roughly by risk.

1. **HTML reporting layer — not read.**
   `html_report.py` (~250 changed lines) and `templates/base.html` (~370 changed
   lines) were *not* reviewed line-by-line; only confirmed that a full report
   renders (6.8 MB HTML) on the real BAM. Need: a read-through for missing-key
   access on metrics that may be absent (e.g. metrics gated behind `--fasta` or
   optional flags), and XSS/escaping of any user-derived strings injected into
   the template.

2. **Offset prediction methods — logic read, not unit-tested on real signal.**
   `modules.py`: `asite_calculation_per_readlength` and the three branches
   `changepoint` (~:1440), `ribowaltz` (~:1460), `tripsviz` (~:1474), plus
   `ribowaltz_psite_prediction` (~:1247). The frame-calibration that sits on top
   (`qc.py::_frame_calibrated_offsets`) was reasoned through and is correct for
   the ±1 nudge, but none of these have a focused test asserting the *offset
   values* produced on a known-periodic input. Need: a fixture BAM with a known
   P-site offset per read length and assertions per method.

3. **`offset_frame_adjustments` interaction with bounds.**
   `qc.py::_frame_calibrated_offsets`: when the ±1 nudge collides with
   `sanitise_offset` bounds/clamping, the new offset can differ from
   `old ± shift` yet still be applied (only the `new == old` case is skipped).
   Confirm a clamped nudge cannot land the dominant frame on the *wrong* non-zero
   frame. Single-pass only (no iteration) — confirm one pass always suffices.

4. **`_recover_alignment_tags_from_pysam` positional assumption.**
   `bam_processing.py`: recovers MAPQ/NH/XA by iterating pysam primary reads and
   aligning them **positionally** to oxbow rows, guarded only by total-length
   equality. Verified safe for STAR transcriptome BAMs (0 secondary alignments,
   per-split temp files). **Not** verified for: genome BAMs, BAMs containing
   secondary alignments, or any reader path where oxbow's row order differs from
   `fetch(until_eof=True)` minus secondaries. Need: an explicit per-read identity
   check (e.g. match on `qname`) instead of trusting order, or at least a test
   with secondary alignments present.

5. **XA tag handling is effectively a no-op on string XA.**
   `bam_processing.py` stores `XA` then `pd.to_numeric(..., errors="coerce")`;
   real BWA-style `XA:Z:` is a semicolon list of loci (a string) → coerced to
   `NaN` → the XA branch in `qc.py::calculate_alignment_stats` never fires. If XA
   support is intended, parse the locus count from the string. If not, drop the
   XA plumbing to avoid implying support that isn't there.

6. **New alignment-rate metrics: semantics & thresholds.**
   `rpf_multimapper_rate` / `alignment_multimapper_rate` descriptions say
   "genomic loci"; on a **transcriptome** BAM these are really alignment/
   transcript multiplicity. On `SRR11005875` the number (0.94) is genuine
   because STAR's NH and MAPQ agree there — but the wording can alarm users and
   could be wrong on a BAM where multi-isoform reads keep MAPQ=255. Also the
   `max_mins` ranges (`[0,1]`) and any QC pass/fail thresholds for these three
   metrics have **not** been calibrated against real distributions.

7. **Example reports are stale.**
   `example-reports/*` were generated against the *broken* annotation (this is
   the source of the 0.366 dominance in the repo). Regenerate against the fixed
   annotation before release or they misrepresent the tool.

8. **`prepare` output dropped the genomic CDS columns.**
   The current `gff_df_to_cds_df` emits only 4 columns; `parse_annotation` still
   special-cases optional `genomic_cds_starts`/`genomic_cds_ends`. Confirm
   nothing depends on them (a grep suggests they are unused in computation) and
   either remove the optional handling or have `prepare` emit them — pick one for
   consistency.

9. **Behaviour change in `filter_unique_mappers` not covered by a test.**
   It now returns an **empty** frame (not the full df) when there is no
   uniqueness signal (no NH, no MAPQ). This is the scenario that produced the
   `plots.py` crash now fixed. Add a regression test asserting graceful
   empty-result behaviour through the metric layer, not just the plot.

---

## 2. Optimisation opportunities

### Done this session
- `modules.py::frame_safe_unique_fragments` — Python groupby loop + per-group
  `concat` → vectorised `groupby.transform('nunique')` + `drop_duplicates`.
  **~411× faster** (54 s → 0.13 s for 150k fragments), identical output.
- `modules.py::reading_frame_triangle` — per-transcript Python groupby loop →
  single grouped `unstack`. **~135× faster** (35 s → 0.26 s for 600k rows),
  identical output. Removed dead `if len(proportion) < 3` branch.

### Open
1. **`_recover_alignment_tags_from_pysam` doubles BAM decode.**
   `bam_processing.py`: a full extra pysam pass per split, building three Python
   lists row-by-row. Meaningful cost on deep libraries. Options: skip entirely
   when oxbow already returns real (non-NaN) MAPQ; vectorise tag pull; or read
   tags in the same pass that oxbow uses.

2. **`prepare` holds the whole GFF/GTF in memory.**
   `file_parser.py::parse_gff` via `gffpandas.read_gff3` → ~2.5 GB resident for a
   3.1 GB v49 GTF. Will OOM on small machines/CI. Roadmap already notes
   `pyranges1`/`pr.read_gff3()` as the streaming replacement — worth doing, also
   gives a speedup. Document a minimum-RAM requirement in the meantime.

3. **`rust.py::run_rust` builds the CDS lookup with `iterrows()`** (~:150) over
   the full annotation (~200k rows). Replace with
   `dict(zip(annotation_df["transcript_id"], zip(cds_start, cds_end)))`. Only
   runs under `--fasta`, so lower priority.

4. **Minor `iterrows`/`apply`** at `file_splitting.py:68` (idxstats, small) and
   `plots.py:1172` (`apply` over a small metrics frame) — low impact, leave
   unless touched.

---

## 3. Code cleanliness

1. **Dead code — remove.**
   - `modules.py::proportion_of_kmer` is now unused (was only called by the old
     `reading_frame_triangle`). 
   - `modules.py::get_cart_point` is unused (the ternary→cartesian conversion is
     commented out in `reading_frame_triangle`).
   Removing both also drops the only consumers of the mid-file numpy-typing
   import.

2. **Mid-file imports.**
   `modules.py:1151` `from numpy.typing import NDArray` and `:1154`
   `from typing import cast as _cast` should move to the top import block (or
   disappear with `get_cart_point`).

3. **Duplicated uniqueness logic.**
   The NH/XA/MAPQ branching in `qc.py::calculate_alignment_stats` and
   `modules.py::filter_unique_mappers` independently re-implement "what counts as
   a unique fragment." Factor into one helper so the report's multimapper rate
   and the frame-sensitive filter can never drift apart.

4. **Magic numbers.**
   - `file_parser.py::validate_annotation_coordinates` uses a hard-coded 1%
     overrun tolerance — promote to a named constant with a comment.
   - `qc.py::_default_offset_for_target` hard-codes 12 (p_site) / 15 (a_site).

5. **`should_calculate_metric` string-variant matching** (`qc.py`) — already
   flagged in `TODO.md`; replace with an explicit metric-ID → config-key
   registry. Touches every new metric added this cycle.

6. **`reading_frame_triangle` docstring** still describes returning "cartesian
   coordinates of the triangle plot"; it now returns raw per-frame counts. Fix
   the docstring (the cartesian step happens in `plots.py`).

7. **Naming.** `calculate_alignment_stats` returns `unique_read_sequences` and
   `reported_alignment_rows` that are equal when `read_name` is unique — confirm
   the intended distinction and document it, or collapse.

---

## Suggested order for the next session
1. Remove dead code + fix mid-file imports (§3.1–3.2) — fast, low-risk.
2. Add offset/calibration unit tests (§1.2–1.3) — highest correctness value.
3. Harden `_recover_alignment_tags_from_pysam` (§1.4) — correctness on non-STAR.
4. Regenerate example reports (§1.7) — release blocker.
5. Optimise the pysam recovery pass and `prepare` memory (§2.1–2.2).
