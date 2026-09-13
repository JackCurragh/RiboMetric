# Metric contract

**Status:** Phase 0 draft, 2026-09-13. Authoritative once reviewed.

This document says what every number RiboMetric emits is *supposed* to be:
its input, its mathematical definition, its range, what it means, the level
it is aggregated at, and how it behaves on degenerate input. It is the
reference that tests are written against.

It covers **measurement**, not judgement. Whether 0.81 deserves PASS belongs
to the scoring design (`METRICS_DESIGN.md`) and to calibration, which comes
later. This contract only fixes what 0.81 *is*.

Where the implementation does not yet match the contract, the entry says so
and names the observation in `RiboMetric-Manuscript/OBSERVATIONS.md` (IDs
`O-NNN`). Section 7 collects those divergences. Section 8 lists the decisions
that still need the owner, marked **[decision D-n]** where they arise.

`registry.py` records name, unit and direction for each metric, and
`METRICS.md` is generated from it. Both must agree with this contract; a test
that checks the registry against it is a Phase 1 item.

## 1. Conventions

These apply to every entry below unless the entry says otherwise.

**Coordinates.** Positions are 0-based offsets into the spliced transcript,
i.e. the reference sequence of a transcriptome BAM.

- `cds_start` is the index of the first nucleotide of the start codon.
- `cds_end` is the exclusive end of the annotated CDS: one past the last
  nucleotide covered by the annotation's `CDS` features. Whether the stop
  codon lies inside `[cds_start, cds_end)` therefore depends on the
  annotation source. GENCODE GTF `CDS` features exclude it, so there the stop
  codon occupies `[cds_end, cds_end + 3)`. **[decision D-1]**
- `transcript_length` is the sum of exon lengths.

**Read 5′ end.** `p5` is the 0-based transcript position of the first base
of the read *including* 5′ soft-clipped bases:
`p5 = POS − 1 − softclip5`, where POS is the 1-based SAM position.
*Divergence:* the implementation uses the 1-based value, `POS − softclip5`
(O-027).

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
*Divergence:* if any read name occurs twice, the implementation drops all
weights library-wide (O-037), and `cds_enrichment_ratio` never uses them
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

**Missing is not zero.** A quantity that cannot be computed (no reads, a zero
denominator, an input that was not supplied) is `null` or absent. It is never
reported as 0, because 0 is a legitimate value for most of these metrics.
*Divergence:* several metrics report 0 in that situation; listed per entry.
**[decision D-6]**

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
- **Subsample (`-S N`):** intended as a seeded random sample of N reads drawn
  across all references, with the seed recorded. *Divergence:* the
  implementation takes whole references in `samtools idxstats` order until N
  reads, i.e. a header-order prefix of the transcriptome (O-029).
  **[decision D-4]**
- **Alignment tags:** MAPQ 255 is recovered from pysam, because oxbow reads
  STAR's 255 as null. That recovery only proceeds when read names match row
  for row. *Divergence:* on a mismatch, MAPQ stays null and the unique set
  becomes empty when there is no NH tag (O-036).

### 2.2 Sequence background

- **In:** the read sequences of parsed reads.
- **Out:** for each read end, the frequency of every dinucleotide at all
  positions of the read except the terminal one being tested (the first
  position for the 5′ background, the last for the 3′). Computed over unique
  sequences, unweighted. **[decision D-9]**
- **Skipped when:** the BAM has no stored sequences, or with
  `--skip-sequence-metrics`. *Only* the terminal-bias and
  nucleotide-composition outputs should depend on this stage. *Divergence:* the
  whole annotation block is currently gated on it, and the run crashes without
  it (O-026).

### 2.3 Annotation join

- **In:** parsed reads; the annotation TSV (`transcript_id`, `cds_start`,
  `cds_end`, `transcript_length`) written by `prepare`.
- **Out:** annotated reads.
- **Discards:** reads on transcripts absent from the annotation. They still
  count toward stage-2.1 statistics (read-length distribution, mapping
  hygiene, terminal bias).
- **Mode:** `annotation` when at least one read has its site in the CDS
  body; otherwise `annotation_free`. Without an annotation, every
  annotation-based metric is absent. *Divergence:* annotation-free runs with
  read sequences crash (O-025).

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
   a consensus fallback (O-030).
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

| region | condition |
|---|---|
| `five_leader` | `a_site < cds_start` |
| `start_codon` | `a_site == cds_start` |
| `CDS` (body) | `cds_start < a_site < cds_end` |
| `stop_codon` | `a_site == cds_end` |
| `three_trailer` | `a_site > cds_end` |

`start_codon` and `stop_codon` are single positions, not codons, so the
second and third bases of the start codon count as CDS body. The meaning of
`stop_codon` also depends on D-1. **[decision D-2]**

### 2.6 The frame table

`F[L][f]`, the weighted count of frame-safe reads with read length `L` in
frame `f`, restricted to reads whose site lies in `(cds_start + 9,
cds_end − 9)`. That excludes 9 nt (three codons) at each CDS end. Emitted as
`read_frame_distribution`. Every frame metric is computed from it.

*Divergences:* only lengths 20–39 are kept (O-033), and transcripts with
`cds_start == 0` are excluded (O-034). Both are **[decision D-3]**.

### 2.7 Measurement, scoring, evaluation

- Raw metrics (§3) are written to `results["metrics"]`.
- Scores (§4) are derived from them into `results["scores"]` by the rules in
  `scoring.py`, merged with the `scoring:` block of the config.
- The QC verdict (`qc_status.json`) is FAIL if any gated score fails, else
  WARNING if any warns, else PASS; INFO if no gated score exists.
- `RiboMetric evaluate --expected` applies an explicit threshold policy.
  Every metric it names must be present and finite. *Divergence:* its default
  policy (no `--expected`) still uses a separate legacy threshold table
  (deferred item D3 in `SCORING_PHASE1_TASKS.md`).

## 3. Metrics

Entries give **input**, **definition**, **range**, **meaning**,
**aggregation**, **degenerate input** and **status**. "Conforms" means the
implementation matches the definition, as far as reading the code shows. It
is confirmed only when the Phase 1 unit test for that entry exists.

### 3.1 Frame and periodicity

#### `periodicity_dominance`

- **Input:** the frame table `F`, restricted to read lengths within
  `plots.read_frame_distribution` limits (15–40).
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
- **Degenerate input:** a length with 1–99 reads is omitted. *Divergence:* a
  length with 0 reads is reported as 0, and an empty table gives `global` 0;
  both should be absent or `null`.
- **Status:** conforms, apart from the zero cases.

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
- **Degenerate input:** *Divergence:* `global` is 0 when no length passes the
  cutoff; it should be `null`.
- **Status:** conforms.

#### `periodicity_information_weighted_score`

- **Definition:** `Σ_L I_L n_L / Σ_L n_L` over *every* read length in `F`,
  with no share cutoff.
- **Range:** [0, 1]. **Aggregation:** scalar.
- **Degenerate input:** *Divergence:* 0 when `F` is empty; should be `null`.
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
- **Status:** conforms. Depends on D-3, because `F` excludes lengths ≥ 40.

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
- **Status:** *Divergence:* the frequency grid is built with
  `np.fft.fftfreq(n, 1/n)` (integer frequencies), so the selected bin is
  always DC. A perfectly periodic profile scores 0.40 and a flat one 1.00
  (O-001).

#### `periodicity_autocorrelation` *(optional)*

- **Intended definition:** the autocorrelation at lag 3 of the start-codon
  metagene count series (position-ordered), normalised by lag 0. Per read
  length and `global`.
- **Status:** *Divergence:* per-length values are `r₃/r₀`, but `global` is
  `(r₃ − mean r)/mean r`, a different statistic. Positions are also not in
  order when there are gaps (O-024).

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
- **Degenerate input:** *Divergence:* no counts gives 0; should be `null`.
- **Status:** *Divergence:* positions are not position-ordered when some
  have no reads, so bins group non-adjacent positions and the per-length sum
  misaligns (O-024).

#### `uniformity_autocorrelation`, `uniformity_gini_index`, `uniformity_theil_index` *(optional)*

- **Intended definitions:** over the same codon-binned start window.
  - Autocorrelation: the mean normalised autocorrelation at lags 1–4. The
    implementation includes lag 0, so it is the mean of `[1, r₁ … r₄]`.
  - Gini: the Gini coefficient `G` of the bin counts. 0 means even, 1 means
    all counts in one bin; lower is better.
  - Theil: the Theil T index `(1/K) Σ (x_k/μ) ln(x_k/μ)`, in [0, ln K];
    lower is better.
- **Status:** *Divergence:* Gini is reported as `1 − G`. The "Theil" value is
  `1/(1 + Σ_codons H(within-codon frame proportions))`, which measures lack of
  periodicity, not inequality of coverage (O-032). All three are also subject
  to O-024.

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
- **Status:** *Divergence:* `observed` counts rows rather than weights, over
  all annotated reads including ineligible transcripts (O-031).

#### `prop_reads_CDS`, `prop_reads_leader`, `prop_reads_trailer`

- **Input:** the region distribution `D[L][region]` (§5), i.e. annotated
  reads, all mapping classes, weighted.
- **Definition:** per read length, `D[L][region] / Σ_regions D[L]`.
  `global`: the same over all lengths pooled. "CDS" means the body region of
  §2.5.
- **Range:** fraction. **Aggregation:** per read length and `global`.
- **Degenerate input:** a length with no reads gives 0.
- **Status:** conforms.

#### `ratio_cds:leader`, `ratio_cds:trailer`, `ratio_leader:trailer`

- **Definition:** per read length, `D[L][A] / D[L][B]`. `global`: `Σ_L D[L][A] / Σ_L D[L][B]`.
- **Range:** ratio. **Aggregation:** per read length and `global`.
- **Degenerate input:** *Divergence:* a zero denominator gives 0; should be
  `null`.
- **Status:** *Divergence:* only lengths 20–39 are summed (O-033, D-3).

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
- **Degenerate input:** no CDS reads gives 0.
- **Status:** conforms.

#### `start_codon_enrichment_ratio`, `stop_codon_readthrough_ratio`

- **Input:** the metagene of unique annotated reads (all reads under
  `--multimap-filter none`), weighted, at distance from `cds_start` (start)
  or `cds_end` (stop).
- **Definition:**
  - Start: counts at distances −5…20 over counts at 30…50.
  - Stop: counts at 1…30 over counts at −30…−1, relative to `cds_end`.
  - Both are summed over all read lengths.
- **Range:** ratio.
- **Meaning:**
  - Start: a strong start peak points to initiation-inhibitor libraries.
  - Stop: reads past the stop codon, from readthrough, a wrong annotation or
    RNA contamination.
- **Degenerate input:** `null` when the denominator is 0.
- **Status:** conforms. With D-1 unresolved, the stop window's first two
  nucleotides may be the stop codon itself.

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
| `read_length_iqr_fraction` | `(Q₀.₇₅ − Q₀.₂₅)/(Q₀.₉ − Q₀.₁)`, where `Q_x` is the smallest L with cumulative share ≥ x | fraction | 0 when `Q₀.₉ = Q₀.₁` |
| `read_length_cv` | weighted standard deviation over weighted mean (population) | index ≥ 0 | |
| `read_length_max_proportion` | `max_L R[L] / Σ R` | fraction | *Divergence:* division by zero on empty `R` |
| `read_length_bimodality_coefficient` | Sarle's `(g₁² + 1)/(g₂ + 3(n−1)²/((n−2)(n−3)))`, with `g₁`, `g₂` the (biased) skewness and excess kurtosis of the expanded lengths, n = `Σ R`, floored at 0 | index ≥ 0; > 5/9 suggests bimodality | *Divergence:* division by zero for n ≤ 3 |
| `read_length_normality_pvalue` *(optional)* | D'Agostino–Pearson test p-value | [0, 1] | undefined for n < 8 |
| `disome_proportion` | `Σ_{L ∈ [50, 70]} R[L] / Σ R`; window in `qc.disome` | fraction | |

These are diagnostics. No score is attached, and the direction is context
(`registry.py`).

### 3.4 Mapping hygiene

All computed from parsed reads.

- **`duplicate_rate`:** `1 − U / W`, where U is the number of distinct read
  names and W their total weight. It measures the share of reads that are
  collapse duplicates. For an uncollapsed BAM (every weight 1) it is 0 by
  construction, which is not a measurement. **[decision D-7]**
- **`rpf_multimapper_rate`:** the weighted fraction of fragments flagged
  multi-mapping. The flag comes from the first available signal: NH > 1, else
  XA > 0, else MAPQ < 255 where MAPQ is known. The method used is recorded in
  `alignment_stats.multimapper_detection_method`. When no signal exists the
  method is `unavailable` and the rate should be `null`. *Divergence:* it is
  reported as 0.
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
  stay in the denominator. *Divergence:* the 3′ dinucleotide is read in
  reverse (O-028).
- **`terminal_bias_kl_5prime` / `_3prime`:**
  `Σ_x obs(x) log₂(obs(x)/bg(x))` over x with `obs > 0` and `bg > 0`, floored
  at 0. Bits; 0 means no departure from the background.
- **`terminal_bias_max_deviation_5prime` / `_3prime`:** `max_x |obs(x) − bg(x)|`.
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
  transcript with at least 20 weighted annotated reads,
  `FLOSS_t = ½ Σ_L |f_t(L) − f_ref(L)|`, where `f_t` is the transcript's
  read-length distribution and `f_ref` the aggregate over all annotated
  reads. The outputs are the median FLOSS and the fraction of transcripts
  with FLOSS > 0.3. Fraction; `null` when no transcript qualifies. The
  docstring says the reference is the CDS aggregate, but the implementation
  uses all annotated reads. **[decision D-5]**

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
| `reading_frame_triangle` | per transcript, weighted counts by `a_site mod 3`. *Divergence:* this is not the §1 frame, because it ignores `cds_start` and includes UTR reads. |
| `gene_body_coverage` | the 100-bin profile behind §3.2's ramp and drop |
| `library_complexity` | `D(f)` at f = 0.1 … 1.0, plus §3.6 |
| `library_type` | a label with its evidence: `low_quality` if `global` dominance < 0.4 or `global` `prop_reads_CDS` < 0.5; else `initiation` if `start_codon_enrichment_ratio ≥ 3`; else `elongation` |
| `floss` | per-transcript FLOSS scores and §3.6 |
| `rust`, `codon_dwell_times` | full outputs behind §3.7 |
| `metrics_legacy` | pre-2.0 keys derived from the canonical ones (`registry.LEGACY_METRIC_ALIASES`); removed at 2.1 |
| `provenance` | timestamp, command, package version, Python, platform, config file and effective-config hashes, input paths, sizes and SHA-256 |

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

| # | where | contract | implementation | obs. |
|---|---|---|---|---|
| 1 | §1, §2.1 | `p5` is 0-based | 1-based; reported offsets are 1 nt low | O-027 |
| 2 | §3.5 | 3′ dinucleotide read 5′→3′ | reversed | O-028 |
| 3 | §3.2 | metagene position-ordered | zero positions appended at the end | O-024 |
| 4 | §3.1 | Fourier at 1/3 cycle per nt | DC bin | O-001 |
| 5 | §2.2 | only §3.5 needs sequences | annotation block gated on sequences; crash | O-026 |
| 6 | §2.3 | annotation-free runs complete | crash | O-025 |
| 7 | §3.2 | `cds_enrichment_ratio` weighted, one read set | unweighted, mixed sets | O-031 |
| 8 | §3.2 | Gini, Theil as named | `1 − G`; "Theil" measures frame entropy | O-032 |
| 9 | §3.1 | autocorrelation one statistic | per-length ≠ global | this doc |
| 10 | §2.1 | random seeded subsample | header-order prefix | O-029 |
| 11 | §2.1 | tag recovery robust | unique set can empty; XA:Z crash | O-036 |
| 12 | §1 | weights always applied | dropped if any name repeats | O-037 |
| 13 | §1 | missing is `null` | several metrics report 0 | this doc |
| 14 | §2.7 | one threshold system | `evaluate` default uses legacy table | D3 |
| 15 | §5 | triangle uses §1 frame | `a_site mod 3` | this doc |
| 16 | §3.3 | defined on tiny inputs | division by zero | this doc |

Items 1–8 and 11 are correctness bugs, fixed in Phase 1 against a failing
test first. Items 9, 12, 13, 15 and 16 are fixed in Phase 1 once D-6 is
settled. Item 10 waits on D-4 and item 14 on the `evaluate` decision.

## 8. Decisions needed

- **D-1 — the stop codon.** Should `prepare` normalise `cds_end` so the stop
  codon is always inside, or always outside, the CDS, regardless of
  annotation source?
- **D-2 — region granularity.** Should `start_codon` and `stop_codon` be whole
  codons rather than single positions?
- **D-3 — read-length windows.** Keep the frame table's 20–39 nt and the
  region ratios' 20–39 nt, or use all observed lengths with the
  `dominance_min_reads` guard? And keep excluding `cds_start == 0`
  transcripts?
- **D-4 — subsampling.** Random seeded sampling across references, with the
  seed recorded?
- **D-5 — FLOSS reference.** CDS reads or all annotated reads?
- **D-6 — missing values.** Adopt `null` for every "cannot compute" case
  (a breaking change for consumers that expect 0)?
- **D-7 — duplicate rate on uncollapsed BAMs.** Report `null`?
- **D-8 — offset methods.** Three names for one computation (O-030). Keep
  one, or implement genuinely different methods? This is a Phase 4 question.
- **D-9 — background weighting.** The 5′ and 3′ backgrounds count unique
  sequences, but the observed frequencies are weighted. Weight both?
