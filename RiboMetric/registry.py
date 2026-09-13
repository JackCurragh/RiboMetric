"""Canonical description of every metric RiboMetric emits.

This is the single source of truth for the naming and direction contract in
``docs/METRIC_NAMING.md``:

* a **metric** key names the physical quantity that was measured, in its
  natural units and its natural direction. Nothing here bakes in a ``1 - x``
  or ``1/(1 + x)`` transform, so a rate goes up as the library gets worse and
  the key says so;
* a **score** key names the good property, always lies in [0, 1] and is always
  higher-is-better;
* the direction flip lives in exactly one place, the ``method`` field of the
  scoring spec in :mod:`RiboMetric.scoring`.

The reason this module exists rather than the direction being implicit: before
it, ``terminal_bias_maxabs_5prime`` held ``1 - max_deviation`` under a key
named for bias, was scored as if it were the deviation, and shipped inverted
for four releases. Nothing in the codebase asserted which way was up. The
registry plus ``tests/test_registry.py`` now does.
"""

from __future__ import annotations

from typing import Dict, List, NamedTuple, Optional

# Direction of the raw quantity.
HIGHER_BETTER = "higher_better"
LOWER_BETTER = "lower_better"
CONTEXT = "context"  # good direction depends on protocol/library type


class MetricSpec(NamedTuple):
    """One raw metric: what it measures, in what units, which way is good."""

    key: str
    unit: str
    direction: str
    summary: str
    scored_as: Optional[str] = None  # score key, or None for a diagnostic


def _m(
    key: str,
    unit: str,
    direction: str,
    summary: str,
    scored_as: Optional[str] = None,
) -> MetricSpec:
    return MetricSpec(key, unit, direction, summary, scored_as)


METRIC_REGISTRY: Dict[str, MetricSpec] = {
    spec.key: spec
    for spec in [
        # -- Tier 1: is this Ribo-seq? ------------------------------------
        _m(
            "periodicity_dominance",
            "fraction",
            HIGHER_BETTER,
            "Fraction of coding A-sites in the dominant reading frame; the global "
            "value uses one shared dominant frame.",
            "periodicity_dominance_score",
        ),
        _m(
            "periodicity_information",
            "fraction",
            HIGHER_BETTER,
            "Entropy reduction of the frame distribution against a uniform " "three-frame null.",
            "periodicity_information_score",
        ),
        _m(
            "cds_enrichment_ratio",
            "ratio",
            HIGHER_BETTER,
            "Observed CDS-body read fraction over the length-weighted expected " "fraction (E).",
            "cds_enrichment_score",
        ),
        # -- Tier 2: usable for my analysis? ------------------------------
        _m(
            "recommended_read_proportion",
            "fraction",
            HIGHER_BETTER,
            "Fraction of the library carried by read lengths recommended for "
            "frame-sensitive work.",
            "usable_read_fraction_score",
        ),
        _m(
            "uniformity_entropy",
            "fraction",
            HIGHER_BETTER,
            "Normalised entropy of the codon-binned start-codon metagene.",
            "coverage_uniformity_score",
        ),
        _m(
            "marginal_position_discovery_rate",
            "rate",
            LOWER_BETTER,
            "Fraction of reads at the margin of sequencing depth landing on a "
            "position not already seen; high means under-sequenced.",
            "library_saturation_score",
        ),
        # -- Tier 3: technical caveats ------------------------------------
        _m(
            "duplicate_rate",
            "rate",
            LOWER_BETTER,
            "Fraction of reads that are collapsed duplicates.",
            "fragment_uniqueness_score",
        ),
        _m(
            "rpf_multimapper_rate",
            "rate",
            LOWER_BETTER,
            "Fraction of weighted fragments reported at more than one " "alignment location.",
            "rpf_unique_mapping_score",
        ),
        _m(
            "alignment_multimapper_rate",
            "rate",
            LOWER_BETTER,
            "Fraction of alignment rows whose fragment has evidence of another "
            "reported alignment.",
            "alignment_unique_mapping_score",
        ),
        _m(
            "soft_clip_rate_5prime",
            "rate",
            LOWER_BETTER,
            "Fraction of reads with 5' soft-clipping.",
            "terminal_integrity_5prime_score",
        ),
        _m(
            "terminal_bias_kl_5prime",
            "bits",
            LOWER_BETTER,
            "Kullback-Leibler divergence of observed 5' terminal dinucleotide "
            "frequencies from the background.",
            "terminal_evenness_kl_5prime_score",
        ),
        _m(
            "terminal_bias_kl_3prime",
            "bits",
            LOWER_BETTER,
            "Kullback-Leibler divergence of observed 3' terminal dinucleotide "
            "frequencies from the background.",
            "terminal_evenness_kl_3prime_score",
        ),
        _m(
            "terminal_bias_max_deviation_5prime",
            "fraction",
            LOWER_BETTER,
            "Largest absolute deviation of a 5' terminal dinucleotide "
            "frequency from its background frequency.",
            "terminal_evenness_maxdev_5prime_score",
        ),
        _m(
            "terminal_bias_max_deviation_3prime",
            "fraction",
            LOWER_BETTER,
            "Largest absolute deviation of a 3' terminal dinucleotide "
            "frequency from its background frequency.",
            "terminal_evenness_maxdev_3prime_score",
        ),
        _m(
            "floss_aberrant_transcript_fraction",
            "fraction",
            LOWER_BETTER,
            "Fraction of transcripts whose footprint-length profile departs "
            "from the library aggregate beyond the FLOSS cutoff.",
            "footprint_homogeneity_score",
        ),
        # -- Diagnostics: reported raw, never scored ----------------------
        _m(
            "read_length_iqr_fraction",
            "fraction",
            LOWER_BETTER,
            "Interquartile range of the read length distribution as a fraction "
            "of its 10th-90th percentile range.",
        ),
        _m(
            "read_length_cv",
            "index",
            LOWER_BETTER,
            "Coefficient of variation of the read length distribution.",
        ),
        _m(
            "read_length_bimodality_coefficient",
            "index",
            LOWER_BETTER,
            "Sarle's bimodality coefficient of the read length distribution.",
        ),
        _m(
            "read_length_normality_pvalue",
            "pvalue",
            CONTEXT,
            "Normaltest p-value for the read length distribution; a footprint "
            "distribution is not expected to be normal.",
        ),
        _m(
            "read_length_max_proportion",
            "fraction",
            CONTEXT,
            "Proportion of reads at the most frequent read length.",
        ),
        _m(
            "disome_proportion",
            "fraction",
            CONTEXT,
            "Fraction of reads in the di-some read-length window; expected in "
            "a di-some experiment, contamination in a monosome one.",
        ),
        _m("floss_median", "index", LOWER_BETTER, "Median per-transcript FLOSS score."),
        _m(
            "start_codon_enrichment_ratio",
            "ratio",
            CONTEXT,
            "Reads near the start codon relative to the CDS body; very high "
            "values indicate initiation-stalling treatments.",
        ),
        _m(
            "stop_codon_readthrough_ratio",
            "ratio",
            LOWER_BETTER,
            "Reads downstream of the stop codon relative to upstream.",
        ),
        _m(
            "five_prime_ramp_ratio",
            "ratio",
            CONTEXT,
            "A-site density in the 5' portion of the CDS over the body.",
        ),
        _m(
            "three_prime_drop_ratio",
            "ratio",
            CONTEXT,
            "A-site density in the 3' portion of the CDS over the body.",
        ),
        _m(
            "complexity_distinct_positions",
            "count",
            CONTEXT,
            "Distinct A-site positions observed at full depth.",
        ),
        _m(
            "n_recommended_read_lengths",
            "count",
            CONTEXT,
            "Number of read lengths recommended for frame-sensitive work.",
        ),
        _m(
            "periodicity_information_weighted_score",
            "fraction",
            HIGHER_BETTER,
            "Read-depth weighted periodicity information content.",
        ),
        _m("prop_reads_CDS", "fraction", CONTEXT, "Proportion of A-sites falling in the CDS body."),
        _m(
            "prop_reads_leader",
            "fraction",
            CONTEXT,
            "Proportion of A-sites falling in the 5' leader.",
        ),
        _m(
            "prop_reads_trailer",
            "fraction",
            CONTEXT,
            "Proportion of A-sites falling in the 3' trailer.",
        ),
        _m("ratio_cds:leader", "ratio", HIGHER_BETTER, "CDS reads relative to 5' leader reads."),
        _m("ratio_cds:trailer", "ratio", HIGHER_BETTER, "CDS reads relative to 3' trailer reads."),
        _m(
            "ratio_leader:trailer",
            "ratio",
            CONTEXT,
            "5' leader reads relative to 3' trailer reads.",
        ),
        _m(
            "cds_coverage",
            "fraction",
            HIGHER_BETTER,
            "Proportion of CDS positions covered, using the configured "
            "in-frame and minimum-read settings.",
        ),
        _m(
            "cds_coverage_1read_1000tx",
            "fraction",
            HIGHER_BETTER,
            "CDS coverage, any frame, >=1 read, 1000 transcripts.",
        ),
        _m(
            "cds_coverage_100read_100tx",
            "fraction",
            HIGHER_BETTER,
            "CDS coverage, any frame, >=100 reads, 100 transcripts.",
        ),
        _m(
            "cds_coverage_inframe_1read_1000tx",
            "fraction",
            HIGHER_BETTER,
            "In-frame CDS coverage, >=1 read, 1000 transcripts.",
        ),
        _m(
            "cds_coverage_inframe_100read_100tx",
            "fraction",
            HIGHER_BETTER,
            "In-frame CDS coverage, >=100 reads, 100 transcripts.",
        ),
        _m("rust_mean_kl_divergence", "bits", CONTEXT, "Mean RUST codon-metagene KL divergence."),
        _m(
            "codon_dwell_cv",
            "index",
            CONTEXT,
            "Coefficient of variation of A-site codon dwell-times.",
        ),
        _m(
            "codon_dwell_p90_p10",
            "ratio",
            CONTEXT,
            "Ratio of the 90th to 10th percentile codon dwell-time.",
        ),
        _m("proline_dwell", "ratio", CONTEXT, "Relative A-site dwell signal on proline codons."),
        _m("cga_dwell", "ratio", CONTEXT, "Relative A-site dwell signal on the CGA codon."),
        _m(
            "periodicity_autocorrelation",
            "index",
            HIGHER_BETTER,
            "Optional: triplet periodicity from signal autocorrelation.",
        ),
        _m(
            "periodicity_fourier",
            "index",
            HIGHER_BETTER,
            "Optional: Fourier power at the codon frequency.",
        ),
        _m(
            "periodicity_trips-viz",
            "index",
            HIGHER_BETTER,
            "Optional: Trips-Viz style triplet periodicity score.",
        ),
        _m(
            "uniformity_autocorrelation",
            "index",
            HIGHER_BETTER,
            "Optional: autocorrelation-based coverage smoothness.",
        ),
        _m(
            "uniformity_theil_index",
            "index",
            LOWER_BETTER,
            "Optional: Theil inequality index of coding-region coverage.",
        ),
        _m(
            "uniformity_gini_index",
            "index",
            LOWER_BETTER,
            "Optional: Gini coefficient of coding-region coverage.",
        ),
    ]
}


# ---------------------------------------------------------------------------
# Backwards compatibility.
#
# Every key RiboMetric <= 1.4.x emitted, mapped to the canonical metric it is
# now derived from plus the transform that reproduces its old value. Emitted
# under ``results["metrics_legacy"]`` for one minor cycle so existing consumers
# (riboseq.org ingestion, --expected policies, cohort tables) keep working
# while they migrate. Remove at v2.1.
# ---------------------------------------------------------------------------

IDENTITY = "identity"
ONE_MINUS = "one_minus"  # old value was 1 - new value
INVERSE_1P = "inverse_1p"  # old value was 1 / (1 + new value)

LEGACY_METRIC_ALIASES: Dict[str, tuple] = {
    # Terminal bias: three names for 1/(1+KL), two for 1-max_deviation.
    "terminal_bias_kl_5prime_raw": ("terminal_bias_kl_5prime", IDENTITY),
    "terminal_bias_kl_3prime_raw": ("terminal_bias_kl_3prime", IDENTITY),
    "terminal_bias_kl_5prime": ("terminal_bias_kl_5prime", INVERSE_1P),
    "terminal_bias_kl_3prime": ("terminal_bias_kl_3prime", INVERSE_1P),
    "terminal_bias_kl_5prime_score": ("terminal_bias_kl_5prime", INVERSE_1P),
    "terminal_bias_kl_3prime_score": ("terminal_bias_kl_3prime", INVERSE_1P),
    "terminal_nucleotide_bias_distribution_5_prime_metric": ("terminal_bias_kl_5prime", INVERSE_1P),
    "terminal_nucleotide_bias_distribution_3_prime_metric": ("terminal_bias_kl_3prime", INVERSE_1P),
    "terminal_bias_maxabs_5prime": ("terminal_bias_max_deviation_5prime", ONE_MINUS),
    "terminal_bias_maxabs_3prime": ("terminal_bias_max_deviation_3prime", ONE_MINUS),
    "terminal_nucleotide_bias_max_absolute_metric_5_prime_metric": (
        "terminal_bias_max_deviation_5prime",
        ONE_MINUS,
    ),
    "terminal_nucleotide_bias_max_absolute_metric_3_prime_metric": (
        "terminal_bias_max_deviation_3prime",
        ONE_MINUS,
    ),
    # Read length family.
    "read_length_distribution_IQR_metric": ("read_length_iqr_fraction", ONE_MINUS),
    "read_length_distribution_coefficient_of_variation_metric": ("read_length_cv", INVERSE_1P),
    "read_length_distribution_bimodality_metric": (
        "read_length_bimodality_coefficient",
        INVERSE_1P,
    ),
    "read_length_distribution_normality_metric": ("read_length_normality_pvalue", ONE_MINUS),
    "read_length_distribution_maxprop_metric": ("read_length_max_proportion", IDENTITY),
    # Mapping hygiene.
    "multimapper_rate": ("rpf_multimapper_rate", IDENTITY),
    "unique_rpf_rate": ("rpf_multimapper_rate", ONE_MINUS),
    # CDS coverage: ten keys for five numbers.
    "CDS_coverage_metric": ("cds_coverage", IDENTITY),
    "CDS_coverage_metric_not_inframe_1read_1000tx": ("cds_coverage_1read_1000tx", IDENTITY),
    "cds_coverage_not_inframe_1read_1000tx": ("cds_coverage_1read_1000tx", IDENTITY),
    "CDS_coverage_metric_not_inframe_100read_100tx": ("cds_coverage_100read_100tx", IDENTITY),
    "cds_coverage_not_inframe_100read_100tx": ("cds_coverage_100read_100tx", IDENTITY),
    "CDS_coverage_metric_inframe_1read_1000tx": ("cds_coverage_inframe_1read_1000tx", IDENTITY),
    "CDS_coverage_metric_inframe_100read_100tx": ("cds_coverage_inframe_100read_100tx", IDENTITY),
}


def _apply_legacy_transform(value: float, transform: str) -> Optional[float]:
    if transform == IDENTITY:
        return float(value)
    if transform == ONE_MINUS:
        return 1.0 - float(value)
    if transform == INVERSE_1P:
        return 1.0 / (1.0 + float(value))
    return None


def build_legacy_metrics(metrics: Dict[str, object]) -> Dict[str, float]:
    """Reproduce the pre-2.0 metric keys from the canonical ones.

    Only scalar metrics are translated; per-read-length dictionaries kept their
    keys across the rename and need no alias.
    """
    legacy: Dict[str, float] = {}
    for old_key, (new_key, transform) in LEGACY_METRIC_ALIASES.items():
        if new_key not in metrics:
            continue
        value = metrics[new_key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        converted = _apply_legacy_transform(float(value), transform)
        if converted is not None:
            legacy[old_key] = converted
    return legacy


def metrics_by_direction(direction: str) -> List[str]:
    """All registered metric keys with the given direction."""
    return sorted(key for key, spec in METRIC_REGISTRY.items() if spec.direction == direction)


def score_key_for(metric_key: str) -> Optional[str]:
    spec = METRIC_REGISTRY.get(metric_key)
    return spec.scored_as if spec else None


SCORED_METRICS: Dict[str, str] = {
    spec.scored_as: spec.key for spec in METRIC_REGISTRY.values() if spec.scored_as is not None
}
