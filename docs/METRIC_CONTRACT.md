# Metric contract

**Status:** reviewed 2026-09-14; authoritative. Phase 0 draft 2026-09-13,
with the Phase 1 fixes applied (§7) and the decisions of §8 taken.

This document says what every number RiboMetric emits is *supposed* to be:
its input, its mathematical definition, its range, what it means, the level
it is aggregated at, and how it behaves on degenerate input. It is the
reference that tests are written against.

It covers **measurement**, not judgement. Whether 0.81 deserves PASS belongs
to the scoring design (`METRICS_DESIGN.md`) and to calibration, which comes
later. This contract only fixes what 0.81 *is*.

Where the implementation does not yet match the contract, the entry says so
and names the observation in `RiboMetric-Manuscript/OBSERVATIONS.md` (IDs
`O-NNN`). Section 7 collects those divergences. Section 8 records the
decisions that settled open definitions, marked **[D-n]** where they apply;
each has a record in `RiboMetric-Manuscript/decisions/`.

`registry.py` records name, unit and direction for each metric, and
`METRICS.md` is generated from it. Both must agree with this contract; a test
that checks the registry against it is a Phase 1 item.

## 1. Conventions

These apply to every entry below unless the entry says otherwise.

**Coordinates.** Positions are 0-based offsets into the spliced transcript,
i.e. the reference sequence of a transcriptome BAM.

- `cds_start` is the index of the first nucleotide of the start codon.
- `cds_end` is the 0-based index of the first nucleotide of the stop codon,
  so the stop codon occupies `[cds_end, cds_end + 3)`, outside the CDS
  **[D-1]**. GTF `CDS` features already exclude the stop codon. For GFF3,
  `prepare` removes the final 3 nt when a `stop_codon` feature is covered by
  the CDS span; an already-exclusive CDS is left unchanged.
- `cds_start_complete` is true when the CDS has an annotated start codon and
  its 5′-most `CDS` feature has phase 0 **[D-3]**. `cds_start == 0` is a
  valid coordinate and is not itself evidence of an incomplete CDS.
- `transcript_length` is the sum of exon lengths.

**Read 5′ end.** `p5` is the 0-based transcript position of the first base
of the read *including* 5′ soft-clipped bases:
`p5 = POS − 1 − softclip5`, where POS is the 1-based SAM position.
Fixed in Phase 1: the implementation used the 1-based value,
`POS − softclip5`, so every reported offset was 1 nt low (O-027).

**Read length.** `L` is the number of query bases consumed by the alignment:
the sum of CIGAR `M`, `I`, `=` and `X` operations. Soft- and hard-clipped
bases are not part of the footprint. If the CIGAR is missing, the sequence
length is used.

**Offsets.** An offset `o(L)` is the 0-based distance from the read 5′ end to
the first nucleotide of the target codon: `site = p5 + o(L)`. The A-site
offset is the P-site offset plus 3. For a canonical 28-nt monosome footprint
these are P 12 and A 15. The analysis site is called `a_site` throughout,
even when the target is the P-site.

**Weights.** A BAM record whose read name ends `_xN` stands for `N` identical
reads (a collapsed library) and carries weight `w = N`; otherwise `w = 1`.
Every count in this document is a sum of weights unless it says rows.
Fixed in Phase 1: weights used to be dropped library-wide whenever any read
name occurred twice (O-037), and `cds_enrichment_ratio` never used them
(O-031).

**Frame.** For a read assigned to a transcript with a CDS,
`frame = (a_site − cds_start) mod 3`. Frame 0 places the site on the first
nucleotide of a codon.

**Read sets.** Four sets recur:

| set | definition |
|---|---|
| parsed reads | primary alignments (`-F 256`) of the subsample (§2.1) |
| annotated reads | parsed reads whose reference (the part of the name before the first `\|`) is a transcript in the annotation |
| unique reads | parsed reads with NH = 1 when an NH tag is present, else MAPQ = 255 where MAPQ is known; all parsed reads when `--multimap-filter none` |
| frame-safe reads | unique annotated reads whose fragment resolves to one gene, one transcript and one frame, collapsed to one row per fragment and read length |

**Aggregation.** Most metrics are reported per read length (integer key `L`)
and as `global`, which pools all read lengths by summing weights, not by
averaging per-length values. Entries that differ say so.

**Missing is not zero.** A quantity is `null` exactly when its input set is
empty or its denominator is zero, and absent when its input was not supplied
**[D-6]**. A computed 0 is reported as 0. Any remaining producer-specific
exception is listed with the metric entry.

**Units.** A *fraction* lies in [0, 1]. A *ratio* is ≥ 0 and unbounded. *Bits*
are base-2 information. A *count* is an integer. An *index* is a
dimensionless statistic with its own range, stated per entry.

## 2. Pipeline stages

A run is a fixed sequence of stages. Each stage lists what it consumes,
produces, defaults to and discards.

### 2.1 Parse

- **In:** a coordinate-sorted, indexed BAM aligned to the transcriptome.
- **Out:** one row per primary alignment: read name, `L`, reference,
  `p5`, 5′ soft-clip length, first and last dinucleotide, weight, MAPQ, NH, XA.
- **Discards:** secondary alignments (`samtools view -F 256`); unmapped reads.
  Supplementary alignments are kept.
- **Subsample (`-S N`):** a seeded random sample of reads **[D-4]**. A read is
  kept when a hash of its name and the seed falls below N / M, where M is the
  number of mapped primary records. Every record of a read is kept or
  dropped together, and the sample does not depend on record order. The seed
  (default 42), the fraction and the realised count go into `provenance`.
  Implemented and covered by deterministic pipeline tests.
- **Alignment tags:** MAPQ 255 is recovered from pysam, because oxbow reads
  STAR's 255 as null. Tags are joined positionally when both readers list the
  same records in the same order, and on read name otherwise. Fixed in
  Phase 1: a mismatch used to leave MAPQ null, which emptied the unique set
  when there was no NH tag (O-036).

### 2.2 Sequence background

- **In:** the read sequences of parsed reads.
- **Out:** for each read end, the frequency of every dinucleotide at all
  positions of the read except the terminal one being tested (the first
  position for the 5′ background, the last for the 3′), weighted by `w` like
  the observed frequencies **[D-9]**. Implemented and covered by weighted
  equivalence tests.
- **Skipped when:** the BAM has no stored sequences, or with
  `--skip-sequence-metrics`. Only the terminal-bias and
  nucleotide-composition outputs depend on this stage. Fixed in Phase 1: the
  whole annotation block used to be gated on it, and the run crashed without
  it (O-026).

### 2.3 Annotation join

- **In:** parsed reads; the annotation TSV (`transcript_id`, `cds_start`,
  `cds_end`, `transcript_length`, and `cds_start_complete` once `prepare`
  writes it) produced by `prepare`.
- **Out:** annotated reads.
- **Discards:** reads on transcripts absent from the annotation. They still
  count toward stage-2.1 statistics (read-length distribution, mapping
  hygiene, terminal bias).
- **Mode:** `annotation` when at least one read has its site in the CDS
  body; otherwise `annotation_free`. Without an annotation, every
  annotation-based metric is absent, including the frame table and every
  frame metric. Fixed in Phase 1: annotation-free runs with read sequences
  crashed (O-025).

### 2.4 Offset inference

Skipped when offsets are supplied (`--offset-read-length`,
`--offset-read-specific`, `--offset-global`), which are applied as given.

1. Take unique annotated reads that resolve to one gene, keeping one
   representative transcript per gene and read length (the one with the most
   weight).
2. Per read length, build the 5′-end metagene `p5 − cds_start` over
   [−50, 20].
3. Take the P-site offset as the position of the maximum count in the window
   −18…−10, negated. `changepoint`, `tripsviz` and `ribowaltz` all reduce to
   this argmax by default; `ribowaltz` adds an optional prominence check and
   a consensus fallback (O-030). The offset audit record names the
   computation that ran, including whether those options were active
   **[D-8]**. Implemented; the audit also records requested/effective method,
   external-supply status, raw/final mappings and frame adjustments.
4. Add 3 for the A-site target.
5. Keep the offset only if it lies within [8, 20], is at most ⌊2L/3⌋ and is
   below L. Otherwise use the default (P 12 / A 15), itself clipped to fit
   the read.
6. **Frame calibration:** for a read length with at least 50 frame-safe CDS
   reads whose dominant frame is not 0 and holds at least 0.6 of those reads,
   shift the offset by ±1 toward frame 0. Apply the shift only if it stays in
   bounds.

Read lengths with no reads in the offset set receive the clipped default.

### 2.5 Assignment

- `a_site = p5 + o(L)`.
- Region of each annotated read, by `a_site`:

| region | condition | proposed (D-2) |
|---|---|---|
| `five_leader` | `a_site < cds_start` | same |
| `start_codon` | `a_site == cds_start` | `cds_start ≤ a_site < cds_start + 3` |
| `CDS` (body) | `cds_start < a_site < cds_end` | `cds_start + 3 ≤ a_site < cds_end` |
| `stop_codon` | `a_site == cds_end` | `cds_end ≤ a_site < cds_end + 3` |
| `three_trailer` | `a_site > cds_end` | `a_site ≥ cds_end + 3` |

`start_codon` and `stop_codon` are single positions, not codons, so the
second and third bases of the start codon count as CDS body. Whole codons are
proposed, not accepted **[D-2]**; until they are, the middle column is the
contract. Under D-1, `stop_codon` is the first nucleotide of the stop codon.

### 2.6 The frame table

`F[L][f]`, the weighted count of frame-safe reads with read length `L` in
frame `f`, restricted to reads whose site lies in `(cds_start + 9,
cds_end − 9)` on transcripts with `cds_start_complete` true. That excludes
9 nt (three codons) at each CDS end, and CDSs whose frame is unknown. Every
observed read length is kept **[D-3]**. Emitted as
`read_frame_distribution`. Every frame metric is computed from it.

Plot limits in `plots.read_frame_distribution` do not change the numerical
frame table or periodicity metrics; they are applied only to presentation
metrics such as `periodicity_trips-viz`.

### 2.7 Measurement, scoring, evaluation

- Raw metrics (§3) are written to `results["metrics"]`.
- Scores (§4) are derived from them into `results["scores"]` by the rules in
  `scoring.py`, merged with the `scoring:` block of the config.
- The QC verdict (`qc_status.json`) is FAIL if any gated score fails, else
  WARNING if any warns, else PASS; INFO if no gated score applies. In
  `annotation` mode a gated score whose metric is `null` or absent counts as
  FAIL, for missing evidence **[D-6]**. Implemented and covered by QC and
  evaluation regression tests.
- `RiboMetric evaluate --expected` applies an explicit threshold policy.
  Every metric it names must be present and finite. Without `--expected`,
  `evaluate` applies the same scoring rules and gate as `qc_status.json`
  **[D-10]**. Explicit `--expected` policies remain an opt-in override.

## 3. Metrics

Entries give **input**, **definition**, **range**, **meaning**,
**aggregation**, **degenerate input** and **status**. "Conforms" means the
implementation matches the definition, as far as reading the code shows. It
is confirmed only when the Phase 1 unit test for that entry exists.

### 3.1 Frame and periodicity

#### `periodicity_dominance`

- **Input:** the frame table `F`, every observed read length **[D-3]**. Plot
  limits (15–40) affect presentation only.
- **Definition:** per read length with `n_L = Σ_f F[L][f] ≥ 100`:
  `max_f F[L][f] / n_L`. `global`: `max_f Σ_L F[L][f] / Σ_L n_L`, one frame
  shared by all lengths. `global_by_read_length_max`:
  `Σ_L max_f F[L][f] / Σ_L n_L`, each length at its own best frame.
- **Range:** [1/3, 1] for any length with reads. 1/3 is the random-frame
  floor, 1 is perfect single-frame occupancy.
- **Meaning:** the fraction of CDS footprints whose site falls in the
  dominant frame. It says whether frame-dependent analysis is defensible.
  `global` below `global_by_read_length_max` means lengths disagree on the
  dominant frame, i.e. offsets are wrong for some lengths.
- **Aggregation:** per read length, `global`, `global_by_read_length_max`.
- **Degenerate input:** a length with 0–99 reads is omitted. An empty table
  gives `global: null` (D-6); dominance cannot be below 1/3.
- **Status:** conforms.

#### `periodicity_information`

- **Input:** `F`, as for dominance.
- **Definition:** per read length,
  `I_L = 1 − H(p_L) / log₂3`, where `p_L` is the frame distribution and `H`
  its Shannon entropy in bits. A length is reported only if it carries more
  than `qc.read_frame_distribution.3nt_count_cutoff` (0.05) of all
  frame-table reads. `global` is the `n_L`-weighted mean of the reported
  `I_L`.
- **Range:** [0, 1]. 0 means a uniform frame distribution, 1 a single frame.
  For a dominant fraction d with the remainder split evenly, `I` is 0.054 at
  d = 0.5 and 0.255 at d = 0.7.
- **Meaning:** a cross-check on dominance that uses the whole frame
  distribution rather than only the top frame.
- **Aggregation:** per read length and `global`.
- **Degenerate input:** `global` is `null` when no length passes the cutoff
  (D-6).
- **Status:** conforms.

#### `periodicity_information_weighted_score`

- **Definition:** `Σ_L I_L n_L / Σ_L n_L` over *every* read length in `F`,
  with no share cutoff.
- **Range:** [0, 1]. **Aggregation:** scalar.
- **Degenerate input:** `null` when `F` is empty (D-6).
- **Status:** conforms. It duplicates the `global` of
  `periodicity_information`, except for the cutoff (a Phase 4 question).

#### `recommended_read_proportion`, `n_recommended_read_lengths`

- **Input:** `F` and the read-length distribution `R` (§5).
- **Definition:** read length `L` is recommended iff `n_L ≥ 100`,
  `max_f F[L][f] / n_L ≥ min_periodicity` (0.5, `--min-periodicity`), and
  `R[L] / Σ R ≥ min_read_proportion` (0.05).
  `recommended_read_proportion = Σ_{recommended L} R[L] / Σ R`;
  `n_recommended_read_lengths` is the count of recommended lengths.
- **Range:** fraction; count.
- **Meaning:** how much of the whole library, not only its CDS reads, lies in
  read lengths fit for frame-sensitive work.
- **Aggregation:** scalar. The per-length detail is in
  `recommended_read_lengths` (§5).
- **Degenerate input:** no recommended lengths gives 0 and 0. These are
  genuine values, not missing.
- **Status:** conforms. Changes under D-3, once `F` keeps lengths outside
  20–39.

#### `periodicity_trips-viz` *(optional)*

- **Definition:** per read length, with `t1 ≥ t2` the two largest frame
  counts: 0 if `t1 = 0`, 1 if `t2 = 0`, else `1 − t2/t1`.
  `global = 1 − Σ t2 / Σ t1`, summed over lengths with `t1 > 0` and `t2 > 0`.
- **Range:** [0, 1].
- **Degenerate input:** lengths with a single occupied frame are left out of
  `global`, so a library of perfect lengths has a `global` computed from its
  imperfect ones only.
- **Status:** conforms to this definition. Its fate is decided in Phase 4.

#### `periodicity_fourier` *(optional)*

- **Input:** the start-codon metagene. It is P-site aligned when
  `periodicity.use_psite_aligned` is true (the default), otherwise the
  A-site metagene over [30, 117].
- **Definition:** for the count series ordered by position, with its mean
  subtracted: the power at the frequency bin nearest 1/3 cycle per
  nucleotide, as a fraction of total power. Per read length, and `global` on
  the position-wise sum.
- **Range:** [0, 1]; 0 means no triplet component, 1 means all variation is
  triplet.
- **Status:** conforms. Fixed in Phase 1: the frequency grid was
  `np.fft.fftfreq(n, 1/n)` (integer frequencies), so the selected bin was
  always DC; a perfectly periodic profile scored 0.40 and a flat one 1.00
  (O-001).

#### `periodicity_autocorrelation` *(optional)*

- **Definition:** the autocorrelation at lag 3 of the start-codon metagene
  count series (position-ordered), normalised by lag 0. Per read length, and
  `global` on the position-wise sum.
- **Status:** conforms. Fixed in Phase 1: `global` was
  `(r₃ − mean r)/mean r`, a different statistic from the per-length `r₃/r₀`,
  and positions were out of order when there were gaps (O-024).

### 3.2 Coverage and regions

#### `uniformity_entropy`

- **Input:** the A-site metagene of annotated reads (all mapping classes,
  weighted) at distances `d = a_site − cds_start` in [30, 117] from the start
  codon.
- **Definition:** group the positions into consecutive 3-nt bins starting at
  30, giving K bins; the last may be partial. With `p_k` each bin's share of
  the counts, `U = −Σ p_k log₂ p_k / log₂ K`. Per read length, and `global`
  on the position-wise sum.
- **Range:** [0, 1]. 1 means evenly spread over codons, low means
  concentrated in few codons. Perfect evenness is *not* the biological ideal;
  the metric is meant for flagging concentration.
- **Aggregation:** per read length and `global`.
- **Degenerate input:** no counts gives `null` (D-6).
- **Status:** conforms. Fixed in Phase 1: positions were out of order when
  some had no reads, so bins grouped non-adjacent positions and the global
  sum misaligned (O-024).

#### `uniformity_autocorrelation`, `uniformity_gini_index`, `uniformity_theil_index` *(optional)*

- **Definitions:** over the same codon-binned start window, per read length
  and `global`.
  - Autocorrelation: the mean normalised autocorrelation at lags 1–4.
  - Gini: the Gini coefficient `G` of the bin counts. 0 means even; (K − 1)/K
    means all counts in one of the K bins. Lower is better.
  - Theil: the Theil T index `(1/K) Σ (x_k/μ) ln(x_k/μ)`, in [0, ln K].
    Lower is better.
- **Status:** conforms. Fixed in Phase 1:
  - autocorrelation included lag 0;
  - Gini was reported as `1 − G`, with a rank offset;
  - the "Theil" value was `1/(1 + Σ_codons H(within-codon frame proportions))`,
    which measures lack of periodicity, not inequality of coverage (O-032);
  - all three read positions out of order (O-024).

#### `cds_enrichment_ratio`

- **Input:** annotated reads, weighted, with their regions.
- **Definition:** restrict to *eligible* transcripts, those with
  `transcript_length > 0` and CDS body length `b_t = cds_end − cds_start − 1 > 0`.
  With `W_t` the weight of reads on transcript t and `C_t` the weight of its
  CDS-body reads:
  `observed = Σ C_t / Σ W_t`,
  `expected = Σ W_t (b_t / transcript_length_t) / Σ W_t`,
  `E = observed / expected`.
- **Range:** ratio ≥ 0. E = 1 is what uniform-per-nucleotide (RNA-seq-like)
  coverage would give; E > 1 means CDS enrichment.
- **Meaning:** CDS enrichment relative to what transcript geometry alone
  predicts. This removes the organism and annotation confound in the raw CDS
  fraction.
- **Aggregation:** scalar.
- **Degenerate input:** `null` when no eligible transcript has reads.
- **Status:** conforms. Fixed in Phase 1: `observed` counted rows rather than
  weights, over all annotated reads including ineligible transcripts (O-031).

#### `prop_reads_CDS`, `prop_reads_leader`, `prop_reads_trailer`

- **Input:** the region distribution `D[L][region]` (§5), i.e. annotated
  reads, all mapping classes, weighted.
- **Definition:** per read length, `D[L][region] / Σ_regions D[L]`.
  `global`: the same over all lengths pooled. "CDS" means the body region of
  §2.5.
- **Range:** fraction. **Aggregation:** per read length and `global`.
- **Degenerate input:** `null` for a length with no reads (D-6).
- **Status:** conforms.

#### `ratio_cds:leader`, `ratio_cds:trailer`, `ratio_leader:trailer`

- **Definition:** per read length, `D[L][A] / D[L][B]`. `global`: `Σ_L D[L][A] / Σ_L D[L][B]`.
- **Range:** ratio. **Aggregation:** per read length and `global`.
- **Degenerate input:** a zero denominator gives `null` (D-6).
- **Status:** conforms; all observed read lengths are included by default.

#### `cds_coverage` and its four variants

- **Input:** annotated reads with a CDS-body site, weighted.
- **Definition:** take the top N transcripts by CDS-body weight. The
  denominator is the total body length `Σ (cds_end − cds_start − 1)`, or
  `Σ ⌊(cds_end − cds_start − 1)/3⌋` for in-frame coverage. The numerator is
  the number of distinct `(transcript, a_site)` positions whose weight is at
  least m, keeping only frame-0 positions for in-frame coverage.

  | key | in-frame | m | N |
  |---|---|---|---|
  | `cds_coverage` | `qc.cds_coverage.in_frame_coverage` (true) | 1 | 100 |
  | `cds_coverage_1read_1000tx` | no | 1 | 1000 |
  | `cds_coverage_100read_100tx` | no | 100 | 100 |
  | `cds_coverage_inframe_1read_1000tx` | yes | 1 | 1000 |
  | `cds_coverage_inframe_100read_100tx` | yes | 100 | 100 |

- **Range:** fraction.
- **Meaning:** breadth of coverage on the best-covered coding sequences.
- **Degenerate input:** with no CDS reads no transcript is selected and the
  denominator is 0, so the value is `null` (D-6).
- **Status:** conforms.

#### `start_codon_enrichment_ratio`, `stop_codon_readthrough_ratio`

- **Input:** the metagene of unique annotated reads (all reads under
  `--multimap-filter none`), weighted, at distance from `cds_start` (start)
  or `cds_end` (stop).
- **Definition:**
  - Start: counts at distances −5…20 over counts at 30…50.
  - Stop: counts at 3…32, the 30 nt after the stop codon, over counts at
    −30…−1, relative to `cds_end` **[D-1]**.
  - Both are summed over all read lengths.
- **Range:** ratio.
- **Meaning:**
  - Start: a strong start peak points to initiation-inhibitor libraries.
  - Stop: reads past the stop codon, from readthrough, a wrong annotation or
    RNA contamination.
- **Degenerate input:** `null` when the denominator is 0.
- **Status:** conforms.

#### `five_prime_ramp_ratio`, `three_prime_drop_ratio`

- **Input:** annotated reads with a CDS-body site, weighted.
- **Definition:** bin the relative position `(a_site − cds_start)/(cds_end − cds_start)`
  into 100 bins and normalise to mean 1. Ramp = mean of bins 0–9 over the
  mean of bins 40–59. Drop = mean of bins 90–99 over the same middle.
- **Range:** ratio; 1 means flat.
- **Degenerate input:** `null` when the middle is 0 or there are no reads.
- **Status:** conforms.

### 3.3 Read length

All computed from the read-length distribution `R` of parsed reads (§5).

| key | definition | range | degenerate input |
|---|---|---|---|
| `read_length_iqr_fraction` | `(Q₀.₇₅ − Q₀.₂₅)/(Q₀.₉ − Q₀.₁)`, where `Q_x` is the smallest L with cumulative share ≥ x | fraction | `null` when `Q₀.₉ = Q₀.₁` (D-6) |
| `read_length_cv` | weighted standard deviation over weighted mean (population) | index ≥ 0 | |
| `read_length_max_proportion` | `max_L R[L] / Σ R` | fraction | `null` on empty `R` |
| `read_length_bimodality_coefficient` | Sarle's `(g₁² + 1)/(g₂ + 3(n−1)²/((n−2)(n−3)))`, with `g₁`, `g₂` the (biased) skewness and excess kurtosis of the expanded lengths, n = `Σ R`, floored at 0 | index ≥ 0; > 5/9 suggests bimodality | `null` for n ≤ 3 or a single read length |
| `read_length_normality_pvalue` *(optional)* | D'Agostino–Pearson test p-value | [0, 1] | `null` for n < 8 |
| `disome_proportion` | `Σ_{L ∈ [50, 70]} R[L] / Σ R`; window in `qc.disome` | fraction | |

These are diagnostics. No score is attached, and the direction is context
(`registry.py`).

### 3.4 Mapping hygiene

All computed from parsed reads.

- **`duplicate_rate`:** `1 − U / W`, where U is the number of distinct read
  names and W their total weight. It measures the share of reads that are
  collapse duplicates. For an uncollapsed BAM (every weight 1) it is 0 by
  construction, which is not a measurement, so it is `null` when no record
  carries a `_xN` suffix **[D-7]**.
- **`rpf_multimapper_rate`:** the weighted fraction of fragments flagged
  multi-mapping. The flag comes from the first available signal: NH > 1, else
  XA > 0, else MAPQ < 255 where MAPQ is known. The method used is recorded in
  `alignment_stats.multimapper_detection_method`. When no signal exists the
  method is `unavailable` and the rate is `null` (D-6).
- **`alignment_multimapper_rate`:** the same flag, as a fraction of alignment
  rows.
- **`soft_clip_rate_5prime`:** the weighted fraction of reads with a 5′ soft
  clip.

On a transcriptome BAM, "multi-mapping" includes one read reported on
several isoforms of the same gene, so these rates measure transcript
multiplicity, not distinct genomic loci (O-005).

### 3.5 Terminal bias

Requires the sequence background (§2.2).

- `obs5[x]`: the weighted fraction of parsed reads whose first two bases
  (5′→3′) are dinucleotide x. `obs3[x]`: the same for the last two bases,
  read 5′→3′. Dinucleotides containing N are left out of the numerators but
  stay in the denominator. Fixed in Phase 1: the 3′ dinucleotide was read in
  reverse (O-028).
- **`terminal_bias_kl_5prime`, `terminal_bias_kl_3prime`:**
  `Σ_x obs(x) log₂(obs(x)/bg(x))` over x with `obs > 0` and `bg > 0`, floored
  at 0. Bits; 0 means no departure from the background.
- **`terminal_bias_max_deviation_5prime`, `terminal_bias_max_deviation_3prime`:** `max_x |obs(x) − bg(x)|`.
  Fraction.
- Both are absent when no sequence background exists.

### 3.6 Library

- **`marginal_position_discovery_rate`:** with `c_i` the weight at each
  distinct `(transcript, a_site)` position of annotated reads,
  `D(f) = Σ_i (1 − (1 − f)^{c_i})` is the expected number of distinct
  positions at sampling fraction f. The rate is
  `(D(1) − D(0.95)) / (0.05 · Σ c_i)`, the share of the last 5% of reads that
  land on a new position. Fraction; 0 means saturated. `null` with no reads.
- **`complexity_distinct_positions`:** `D(1)`, the number of distinct
  positions. Count.
- **`floss_median`, `floss_aberrant_transcript_fraction`:** for each
  transcript with at least 20 weighted CDS-body reads,
  `FLOSS_t = ½ Σ_L |f_t(L) − f_ref(L)|`, where `f_t` is the read-length
  distribution of the transcript's CDS-body reads and `f_ref` the aggregate
  over all CDS-body reads **[D-5]**, as in Ingolia et al. 2014. The outputs
  are the median FLOSS and the fraction of transcripts with FLOSS > 0.3.
  Fraction; `null` when no transcript qualifies. Implemented on CDS-body
  reads for both distributions and covered by the contract pipeline tests.

### 3.7 Sequence-dependent *(require `--fasta`)*

- **`rust_mean_kl_divergence`:** the mean over window positions of the KL
  divergence of the RUST codon-occupancy profile (O'Connor et al. 2016) from
  its expectation. Window of 60 codons, A-site at codon 40, CDS trimmed by
  120 nt (5′) and 60 nt (3′). Bits.
- **`codon_dwell_cv`, `codon_dwell_p90_p10`, `proline_dwell`, `cga_dwell`:**
  for each codon c, `dwell(c) = (obs_c / Σ obs) / (exp_c / Σ exp)`. Here
  `obs_c` is the weighted A-site count on codon c and `exp_c` its frequency,
  both in CDS regions trimmed by 15 nt at each end. The outputs are the
  coefficient of variation across codons, the 90th/10th percentile ratio, the
  mean over the four proline codons, and CGA's dwell.

## 4. Scores

`results["scores"]` holds one record per score whose source metric is
present: `{metric, raw, score, status, gate, tier}`. `raw` is the metric's
`global` value (or its scalar value). `score = method(raw)`:

| method | definition |
|---|---|
| `identity` | `clip(raw, 0, 1)` |
| `one_minus_rate` | `clip(1 − raw, 0, 1)` |
| `inverse_linear` | `clip(1 − raw / max_value, 0, 1)` |
| `enrichment_ratio` | 0 if `raw ≤ 1`, else `1 − 1/raw` |

`status` is PASS if the score is at least `pass`, WARNING if at least `warn`,
else FAIL, and INFO when there is no score. The 15 scores, with their
methods, thresholds and gate membership, are defined by `DEFAULT_SCORING` in
`scoring.py` together with the `scoring:` block of `config.yml`. Thresholds
are provisional until calibration, which is outside this programme.

## 5. Other outputs

| key | content |
|---|---|
| `mode` | `annotation` or `annotation_free` (§2.3) |
| `read_length_distribution` | `R[L]`: weighted parsed reads per read length |
| `read_frame_distribution` | the frame table `F[L][f]` (§2.6) |
| `mRNA_distribution` | `D[L][region]` and `D["global"][region]`: weighted annotated reads per region (§2.5) |
| `metagene_profile` | `{start, stop}[L][d]`: weighted unique annotated reads at distance d over [−50, 50] from `cds_start` / `cds_end` |
| `metagene_profile_stats` | reads entering the metagene, and the multimap filter used |
| `recommended_read_lengths` | per length with ≥ 100 frame reads: `periodicity`, `read_proportion`, `n_frame_reads`, `recommended`, `offset`; plus the summaries in §3.1 |
| `computed_offsets`, `computed_offset_target` | final offset per read length (§2.4) and its target |
| `offsets` | offset audit record: source, target, default, bounds, method, raw and final offsets, frame adjustments, applied offsets per read length |
| `alignment_stats` | inputs to §3.4, `multimapper_detection_method`, `mapq_available_rate`, and `samtools flagstat` totals |
| `terminal_nucleotide_bias_distribution` | `obs − bg` per dinucleotide and end when `background_freq` is true, otherwise `obs` |
| `nucleotide_composition` | per read position, the fraction of A/C/G/T |
| `reading_frame_triangle` | per transcript, weighted counts of CDS-body reads by the §1 frame. Fixed in Phase 1: it used `a_site mod 3` over every read. |
| `gene_body_coverage` | the 100-bin profile behind §3.2's ramp and drop |
| `library_complexity` | `D(f)` at f = 0.1 … 1.0, plus §3.6 |
| `library_type` | a label with its evidence: `low_quality` if `global` dominance < 0.4 or `global` `prop_reads_CDS` < 0.5; else `initiation` if `start_codon_enrichment_ratio ≥ 3`; else `elongation` |
| `floss` | per-transcript FLOSS scores and §3.6 |
| `rust`, `codon_dwell_times` | full outputs behind §3.7 |
| `metrics_legacy` | pre-2.0 keys derived from the canonical ones (`registry.LEGACY_METRIC_ALIASES`); removed at 2.1 |
| `provenance` | timestamp, command, package version, Python, platform, config file and effective-config hashes, input paths, sizes and SHA-256; under D-4 also the subsample request, seed, realised fraction and count; offset audit separately records requested/effective method and external/calibrated status |

## 6. Invariants

These hold for every run and are what the Phase 3 invariance tests check.

1. **Order.** Results do not depend on BAM record order, annotation row
   order, transcript order, thread count or read names, except the `_xN`
   suffix.
2. **Reproducibility.** Two runs with the same inputs, config and seed give
   identical metrics.
3. **Weighting.** Collapsing identical reads into one `_xN` record changes no
   metric except `duplicate_rate`.
4. **Frame shift.** Shifting every offset by +1 moves `F` by one frame and
   leaves `periodicity_dominance` per length unchanged.
5. **Isolation.** Removing sequence data changes only §3.5 and
   `nucleotide_composition`.

## 7. Divergence register

| # | where | contract | implementation before Phase 1 | obs. | status |
|---|---|---|---|---|---|
| 1 | §1, §2.1 | `p5` is 0-based | 1-based; reported offsets 1 nt low | O-027 | fixed |
| 2 | §3.5 | 3′ dinucleotide read 5′→3′ | reversed | O-028 | fixed |
| 3 | §3.2 | metagene position-ordered | zero positions appended at the end | O-024 | fixed |
| 4 | §3.1 | Fourier at 1/3 cycle per nt | DC bin | O-001 | fixed |
| 5 | §2.2 | only §3.5 needs sequences | annotation block gated on sequences; crash | O-026 | fixed |
| 6 | §2.3 | annotation-free runs complete | crash | O-025 | fixed |
| 7 | §3.2 | `cds_enrichment_ratio` weighted, one read set | unweighted, mixed sets | O-031 | fixed |
| 8 | §3.2 | Gini, Theil as named | `1 − G`; "Theil" measured frame entropy | O-032 | fixed |
| 9 | §3.1 | autocorrelation one statistic | per-length ≠ global | this doc | fixed |
| 10 | §2.1 | random seeded subsample | header-order prefix | O-029 | fixed (D-4) |
| 11 | §2.1 | tag recovery robust | unique set could empty; XA:Z crash | O-036 | fixed |
| 12 | §1 | weights always applied | dropped if any name repeated | O-037 | fixed |
| 13 | §1 | missing is `null` | several metrics report 0 | this doc | fixed (D-6) |
| 14 | §2.7 | one threshold system | `evaluate` default uses legacy table | D3 | fixed (D-10) |
| 15 | §5 | triangle uses §1 frame | `a_site mod 3` | this doc | fixed |
| 16 | §3.3 | defined on tiny inputs | division by zero | this doc | fixed |
| 17 | §1 | stop codon outside the CDS for every format | GFF3 CDS keeps the stop | O-038 | fixed (D-1) |
| 18 | §1, §2.3 | `cds_start_complete` in the annotation | not written | O-039 | fixed (D-3) |
| 19 | §2.6 | frame table keeps every read length | 20–39 only | O-033 | fixed (D-3) |
| 20 | §2.6 | CDSs excluded by completeness | excluded by `cds_start == 0` | O-034, O-039 | fixed (D-3) |
| 21 | §3.1 | dominance over every read length | plot limits 15–40 | this doc | fixed (D-3) |
| 22 | §3.2 | stop window after the stop codon | 1…30 from `cds_end` | O-038 | fixed (D-1) |
| 23 | §2.2 | background weighted | one count per record | O-040 | fixed (D-9) |
| 24 | §2.7 | `null` gated metric fails the gate | skipped | O-041 | fixed (D-6) |
| 25 | §3.4 | `duplicate_rate` `null` when uncollapsed | 0 | O-042 | fixed (D-7) |
| 26 | §3.6 | FLOSS on CDS-body reads | all annotated reads | O-043 | fixed (D-5) |
| 27 | §2.4 | offset audit names the computation | method name only | O-030 | fixed (D-8) |

Phase 1 and Task 3 (2026-09-14) fixed the items marked fixed above, each
pinned by focused or end-to-end regression tests. Region granularity D-2
remains proposed and requires biological validation before changing the
current single-base assignment semantics.

## 8. Decisions

Taken 2026-09-14 on the owner's delegation. Each has a record in
`RiboMetric-Manuscript/decisions/` with its evidence and the alternatives
considered; a later change gets a new record.

| | question | decision | record | status |
|---|---|---|---|---|
| D-1 | the stop codon | always outside the CDS: `cds_end` is its first nucleotide; `prepare` normalises GFF3 | 0001 | accepted |
| D-2 | region granularity | whole-codon `start_codon` and `stop_codon` | 0002 | proposed |
| D-3 | read-length windows; `cds_start == 0` | every read length (per-length values still need 100 reads); exclude 5′-incomplete CDSs by `cds_start_complete` | 0003 | accepted |
| D-4 | subsampling | seeded, by read-name hash; seed, fraction and realised count recorded | 0004 | accepted |
| D-5 | FLOSS reference | CDS-body reads, for the reference and each transcript | 0005 | accepted |
| D-6 | missing values | `null` for an empty input or a zero denominator; a `null` gated metric fails the gate in `annotation` mode | 0006 | accepted |
| D-7 | duplicate rate on uncollapsed BAMs | `null` | 0007 | accepted |
| D-8 | offset methods | one computation under three names; the audit record names it; Phase 4 settles the names | 0008 | accepted (interim) |
| D-9 | background weighting | weight the background like the observed frequencies | 0009 | accepted |
| D-10 | one evaluation path | `evaluate` without `--expected` uses the scoring rules and gate | 0010 | accepted |

D-2 stays proposed until there is evidence, as the decisions log requires. The
Phase 2 initiation libraries will show how much signal sits on the second and
third nucleotides of the start codon.
