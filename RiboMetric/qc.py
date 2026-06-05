"""
Main script for running qc analysis

Three main modes:
    annotation free: no gff file provided just use the bam file
    annotation based: gff file provided and use the bam file
    sequence based: gff file and transcriptome fasta file
                    provided and use the bam file

"""

import pandas as pd

from .modules import (
    chunked_annotate_reads,
    assign_mRNA_category,
    read_length_distribution,
    read_df_to_cds_read_df,
    terminal_nucleotide_bias_distribution,
    normalise_ligation_bias,
    nucleotide_composition,
    read_frame_distribution,
    read_frame_distribution_annotated,
    mRNA_distribution,
    metagene_profile,
    reading_frame_triangle,
    read_frame_score_trips_viz,
    read_frame_cull,
    asite_calculation_per_readlength,
    a_site_calculation_variable_offset,
    a_site_calculation,
    ribowaltz_psite_prediction,
    gene_body_coverage_ramp,
    library_complexity_curve,
    floss_library_heterogeneity,
    DEFAULT_OFFSET_BOUNDS,
    DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION,
)

from .metrics import (
    read_length_distribution_IQR_normalised_metric as rld_metric,
    read_length_distribution_coefficient_of_variation_metric as rldv_metric,
    read_length_distribution_normality_metric as rldn_metric,
    terminal_nucleotide_bias_KL_divergence as lbd_raw,
    terminal_nucleotide_bias_max_absolute_metric as lbmp_metric,
    read_frame_information_content as rf_info_metric,
    read_frame_information_weighted_score,
    information_metric_cutoff,
    cds_coverage_metric,
    region_region_ratio_metric,
    read_length_distribution_max_prop_metric as rldpp_metric,
    periodicity_autocorrelation,
    uniformity_autocorrelation,
    uniformity_entropy,
    uniformity_theil_index,
    uniformity_gini_index,
    periodicity_dominance,
    fourier_transform,
    read_length_distribution_bimodality,
    proportion_of_reads_in_region,
    recommend_read_lengths,
    classify_library_type,
)
from typing import Any, Dict, Optional


def should_calculate_metric(metric_name: str, config: dict) -> bool:
    """
    Check if a metric should be calculated based on the configuration

    Inputs:
        metric_name: Name of the metric to check
        config: Configuration dictionary containing enabled metrics

    Outputs:
        bool: True if the metric should be calculated, False otherwise
    """
    enabled_metrics = config.get("enabled_metrics", [])

    # If no enabled_metrics specified, calculate all (backwards compatibility)
    if not enabled_metrics:
        return True

    # Check various name formats the metric might be stored as
    metric_variants = [
        metric_name,
        metric_name.replace("_metric", ""),
        metric_name.replace("periodicity_", "").replace("uniformity_", ""),
    ]

    return any(variant in enabled_metrics for variant in metric_variants)


def annotation_mode(
    read_df: pd.DataFrame,
    sequence_data: dict,
    sequence_background: dict,
    annotation_df: pd.DataFrame = pd.DataFrame(),
    config: dict = {},
    fasta_dict: Optional[dict] = None,
) -> dict:
    """
    Run the annotation mode of the qc analysis

    Inputs:
        read_df: Dataframe containing the read information
                (keys are the read names)
        sequence_data: Dictionary containing the sequence data
        sequence_background: Dictionary containing the sequence background
        annotation_df: Dataframe containing the annotation information
        config: Dictionary containing the configuration information
        fasta_dict: Optional FASTA dict from parse_fasta; enables RUST metrics
                    and RNase context analysis when provided.

    Outputs:
        results_dict: Dictionary containing the results of the qc analysis
    """
    print("Calculating A site information...")
    offset_bounds = (
        int(config.get("argument", {}).get("offset_min", DEFAULT_OFFSET_BOUNDS[0])),
        int(config.get("argument", {}).get("offset_max", DEFAULT_OFFSET_BOUNDS[1])),
    )
    max_read_length_fraction = config.get("argument", {}).get(
        "offset_max_read_length_fraction",
        DEFAULT_OFFSET_MAX_READ_LENGTH_FRACTION,
    )
    max_read_length_fraction = (
        None
        if max_read_length_fraction is None
        else float(max_read_length_fraction)
    )
    if ("offset_read_length" in config["argument"]):
        print("Applying specified read length specific offsets")
        read_df = a_site_calculation(read_df,
                                     offset_file=config["argument"][
                                            "offset_read_length"],
                                     offset_type="read_length",
                                     offset_bounds=offset_bounds,
                                     max_read_length_fraction=max_read_length_fraction)

    elif ("offset_read_specific" in config['argument']):
        print("Applying read specific offsets")
        read_df = a_site_calculation(read_df,
                                     offset_file=config["argument"][
                                            "offset_read_specific"],
                                     offset_type="read_specific",
                                     offset_bounds=offset_bounds,
                                     max_read_length_fraction=max_read_length_fraction,
                                     )
    elif ("global_offset" in config["argument"]):
        print("Applying global offset")
        read_df = a_site_calculation(read_df,
                                     global_offset=config["argument"][
                                        "global_offset"],
                                     offset_type="global",
                                     offset_bounds=offset_bounds,
                                     max_read_length_fraction=max_read_length_fraction,
                                     )

    computed_offsets = {}
    if len(annotation_df) > 0:
        annotation = True
        print("Merging annotation and reads")
        if "a_site" not in read_df.columns:
            read_df = a_site_calculation(
                read_df,
                offset_type="global",
                offset_bounds=offset_bounds,
                max_read_length_fraction=max_read_length_fraction,
            )
            annotated_read_df = chunked_annotate_reads(read_df, annotation_df)

            print("assigning mRNA categories")
            annotated_read_df = assign_mRNA_category(annotated_read_df)
            # Optional safeguard for noisy metagenes: require a minimum
            # peak prominence when calling riboWaltz-style offsets.
            min_prom = (
                config.get("qc", {})
                .get("read_frame_distribution", {})
                .get("psite_min_prominence")
            )
            unique_only = (
                config.get("argument", {}).get("multimap_filter", "unique_only")
                == "unique_only"
            )
            offsets = asite_calculation_per_readlength(
                annotated_read_df,
                method=config["argument"]["offset_calculation_method"],
                default_offset=config["argument"].get("global_offset", 15),
                offset_bounds=offset_bounds,
                max_read_length_fraction=max_read_length_fraction,
                min_prominence=min_prom,
                unique_only=unique_only,
            )
            computed_offsets = offsets
            annotated_read_df = a_site_calculation_variable_offset(
                annotated_read_df,
                offsets,
                default_offset=config["argument"].get("global_offset", 15),
                offset_bounds=offset_bounds,
                max_read_length_fraction=max_read_length_fraction,
                validate_offsets=True,
                )
        else:
            annotated_read_df = chunked_annotate_reads(read_df, annotation_df)
            annotated_read_df = assign_mRNA_category(annotated_read_df)

        print("Subsetting to CDS reads")
        cds_read_df = read_df_to_cds_read_df(annotated_read_df)
        if cds_read_df.empty:
            print("No CDS reads found")
            annotation = False

    else:
        annotation = False
        read_df = a_site_calculation(
            read_df,
            offset_type="global",
            offset_bounds=offset_bounds,
            max_read_length_fraction=max_read_length_fraction,
        )

    print("Running modules")

    results_dict: Dict[str, Any] = {
        "mode": ("annotation" if annotation else "annotation_free"),
        "metrics": {}
    }

    if computed_offsets:
        results_dict["computed_offsets"] = {
            str(k): v for k, v in computed_offsets.items()
        }

    #######################################################################
    # ALIGNMENT STATS  (from read_df columns computed during BAM parsing)
    #######################################################################
    _count_arr = read_df["count"].astype(int)
    _total_weighted = int(_count_arr.sum())
    _unique_reads = len(read_df)
    _dup_rate = (
        (1.0 - _unique_reads / _total_weighted) if _total_weighted > 0 else 0.0
    )

    _mapq_arr = read_df["mapq"].astype(int)
    _multimap_mask = (_mapq_arr < 255).astype(int)
    _multimapper_rate = (
        float((_multimap_mask * _count_arr).sum() / _total_weighted)
        if _total_weighted > 0 else 0.0
    )

    results_dict["alignment_stats"] = {
        "total_reads_analysed": _total_weighted,
        "unique_read_sequences": _unique_reads,
        "duplicate_rate": round(_dup_rate, 4),
        "multimapper_rate": round(_multimapper_rate, 4),
    }
    results_dict["metrics"]["duplicate_rate"] = _dup_rate
    results_dict["metrics"]["multimapper_rate"] = _multimapper_rate

    # duplicate_rate is derived from per-read `count` weights, which are only
    # >1 when read names carry a collapse suffix (e.g. ``..._x12``). For
    # un-collapsed BAMs every count is 1, so the rate is necessarily 0 and not
    # informative — flag that explicitly rather than reporting a misleading 0.
    if _total_weighted == _unique_reads:
        print(
            "Note: reads are not collapsed (no '_xN' suffix); duplicate_rate "
            "is reported as 0 and should be treated as not applicable."
        )

    if "soft_clip_5" in read_df.columns:
        _sc5 = read_df["soft_clip_5"].astype(int)
        _sc5_mask = (_sc5 > 0).astype(int)
        _sc5_rate = (
            float((_sc5_mask * _count_arr).sum() / _total_weighted)
            if _total_weighted > 0 else 0.0
        )
        results_dict["alignment_stats"]["soft_clip_rate_5prime"] = round(_sc5_rate, 4)
        results_dict["metrics"]["soft_clip_rate_5prime"] = _sc5_rate

    #######################################################################
    # READ LENGTH DISTRIBUTION
    #######################################################################
    print("> read_length_distribution")
    results_dict["read_length_distribution"] = read_length_distribution(
        read_df
    )

    # Default read length metrics
    if should_calculate_metric("read_length_distribution_IQR", config):
        results_dict["metrics"][
            "read_length_distribution_IQR_metric"
            ] = rld_metric(
            results_dict["read_length_distribution"]
        )

    if should_calculate_metric("read_length_distribution_coefficient_of_variation", config):
        results_dict["metrics"][
            "read_length_distribution_coefficient_of_variation_metric"
            ] = rldv_metric(
            results_dict["read_length_distribution"]
        )

    if should_calculate_metric("read_length_distribution_maxprop", config):
        results_dict["metrics"][
            "read_length_distribution_maxprop_metric"] = rldpp_metric(
            results_dict["read_length_distribution"],
            num_top_readlens=1
        )

    # Optional read length metrics
    if should_calculate_metric("read_length_distribution_bimodality", config):
        results_dict["metrics"][
            "read_length_distribution_bimodality_metric"
            ] = read_length_distribution_bimodality(
                results_dict["read_length_distribution"]
            )

    if should_calculate_metric("read_length_distribution_normality", config):
        results_dict["metrics"][
            "read_length_distribution_normality_metric"
            ] = rldn_metric(
                results_dict["read_length_distribution"]
            )

    # Di-some proportion — fraction of reads in the di-some read-length window.
    # A secondary peak here indicates di-some contamination or a mixed library.
    # The window defaults to 50-70 nt (typical for mammalian di-somes) but is
    # configurable for other organisms / nucleases via qc.disome.
    _disome_cfg = config.get("qc", {}).get("disome", {}) if config else {}
    _disome_lo = int(_disome_cfg.get("lower_limit", 50))
    _disome_hi = int(_disome_cfg.get("upper_limit", 70))
    _rld = results_dict["read_length_distribution"]
    _total_rld = sum(_rld.values()) or 1
    _disome_count = sum(
        v for k, v in _rld.items() if _disome_lo <= int(k) <= _disome_hi
    )
    results_dict["metrics"]["disome_proportion"] = round(
        _disome_count / _total_rld, 4
    )

    #######################################################################
    # TERMINAL NUCLEOTIDE BIAS
    #######################################################################
    if sequence_background:
        print("> terminal_nucleotide_bias_distribution")
        results_dict[
            "terminal_nucleotide_bias_distribution"
            ] = terminal_nucleotide_bias_distribution(
            read_df,
            pattern_length=config[
                "plots"][
                "terminal_nucleotide_bias_distribution"][
                "nucleotide_count"
            ],
        )
        # Terminal nucleotide bias (standard + consolidated key aliases)
        _lbd5_raw = lbd_raw(
            results_dict["terminal_nucleotide_bias_distribution"],
            sequence_background["5_prime_bg"],
            prime="five_prime",
        )
        _lbd3_raw = lbd_raw(
            results_dict["terminal_nucleotide_bias_distribution"],
            sequence_background["3_prime_bg"],
            prime="three_prime",
        )
        _lbd5 = 1 / (1 + _lbd5_raw)
        _lbd3 = 1 / (1 + _lbd3_raw)
        _lbm5 = lbmp_metric(
            results_dict["terminal_nucleotide_bias_distribution"],
            sequence_background["5_prime_bg"],
            prime="five_prime",
        )
        _lbm3 = lbmp_metric(
            results_dict["terminal_nucleotide_bias_distribution"],
            sequence_background["3_prime_bg"],
            prime="three_prime",
        )
        # Legacy keys (preserve)
        results_dict["metrics"]["terminal_nucleotide_bias_distribution_5_prime_metric"] = _lbd5
        results_dict["metrics"]["terminal_nucleotide_bias_distribution_3_prime_metric"] = _lbd3
        results_dict["metrics"]["terminal_nucleotide_bias_max_absolute_metric_5_prime_metric"] = _lbm5
        results_dict["metrics"]["terminal_nucleotide_bias_max_absolute_metric_3_prime_metric"] = _lbm3
        # Consolidated keys (new)
        results_dict["metrics"]["terminal_bias_kl_5prime"] = _lbd5
        results_dict["metrics"]["terminal_bias_kl_3prime"] = _lbd3
        results_dict["metrics"]["terminal_bias_kl_5prime_score"] = _lbd5
        results_dict["metrics"]["terminal_bias_kl_3prime_score"] = _lbd3
        results_dict["metrics"]["terminal_bias_kl_5prime_raw"] = _lbd5_raw
        results_dict["metrics"]["terminal_bias_kl_3prime_raw"] = _lbd3_raw
        results_dict["metrics"]["terminal_bias_maxabs_5prime"] = _lbm5
        results_dict["metrics"]["terminal_bias_maxabs_3prime"] = _lbm3
        if config["plots"][
                "terminal_nucleotide_bias_distribution"][
                "background_freq"]:
            results_dict[
                "terminal_nucleotide_bias_distribution"
                ] = normalise_ligation_bias(
                results_dict["terminal_nucleotide_bias_distribution"],
                sequence_background=sequence_background,
                pattern_length=config[
                    "plots"][
                    "terminal_nucleotide_bias_distribution"][
                    "nucleotide_count"
                ],
            )

        print("> nucleotide_composition")
        results_dict["nucleotide_composition"] = nucleotide_composition(
            sequence_data)

        print("> read_frame_distribution")
        if annotation:
            coding_metagene = metagene_profile(
                annotated_read_df,
                target="start",
                distance_range=[30, 117],
                extend=True,
            )

            # Optionally build P-site aligned metagene for spectral metrics
            psite_aligned = None
            if config.get("periodicity", {}).get("use_psite_aligned", False):
                # Build 5'-end metagene around start codon for P-site detection
                metagene_5p = metagene_profile(
                    annotated_read_df,
                    target="start",
                    distance_range=[-50, 20],
                    position="reference_start",
                    extend=False,
                )
                # Predict P-site offsets (positive ints) per read length
                try:
                    psite_offsets = ribowaltz_psite_prediction(metagene_5p["start"])
                except Exception:
                    psite_offsets = {}

                # Shift start metagene counts by offset so P-site at start codon maps to 0
                aligned_start: Dict[int, Dict[int, int]] = {}
                for rl, counts in metagene_5p["start"].items():
                    off = psite_offsets.get(rl)
                    if off is None:
                        continue
                    shifted: Dict[int, int] = {}
                    for pos, val in counts.items():
                        new_pos = pos + off
                        shifted[new_pos] = shifted.get(new_pos, 0) + val
                    aligned_start[rl] = shifted
                if aligned_start:
                    psite_aligned = {"start": aligned_start, "stop": {}}

            # read_frame_dist must be computed before periodicity metrics
            exclude_nt = config["qc"]["read_frame_distribution"].get("exclude_codons", 9)
            unique_only = (
                config.get("argument", {}).get("multimap_filter", "unique_only")
                == "unique_only"
            )
            read_frame_dist = (
                read_frame_distribution_annotated(cds_read_df, exclusion_length=exclude_nt, unique_only=unique_only)
                if config["qc"]["use_cds_subset"]["read_frame_distribution"]
                and annotation
                else read_frame_distribution_annotated(annotated_read_df, exclusion_length=exclude_nt, unique_only=unique_only)
            )

            # Spectral metagene: prefer P-site aligned if available
            spectral_metagene = psite_aligned if psite_aligned else coding_metagene

            #######################################################################
            # Periodicity
            #######################################################################
            if should_calculate_metric("periodicity_autocorrelation", config):
                selected = select_read_lengths_for_global(
                    read_frame_dist,
                    results_dict["read_length_distribution"],
                    config,
                )
                metagene_for_autocorr = spectral_metagene
                if selected:
                    reduced: Dict[str, Dict[int, Dict[int, int]]] = {"start": {}, "stop": {}}
                    for rl in selected:
                        if rl in spectral_metagene["start"]:
                            reduced["start"][rl] = spectral_metagene["start"][rl]
                        if rl in spectral_metagene["stop"]:
                            reduced["stop"][rl] = spectral_metagene["stop"][rl]
                    metagene_for_autocorr = reduced
                results_dict["metrics"]["periodicity_autocorrelation"] = periodicity_autocorrelation(
                    metagene_for_autocorr.copy()
                )
            if should_calculate_metric("periodicity_fourier", config):
                selected = select_read_lengths_for_global(
                    read_frame_dist,
                    results_dict["read_length_distribution"],
                    config,
                )
                metagene_for_fourier = spectral_metagene
                if selected:
                    reduced = {"start": {}, "stop": {}}
                    for rl in selected:
                        if rl in spectral_metagene["start"]:
                            reduced["start"][rl] = spectral_metagene["start"][rl]
                        if rl in spectral_metagene["stop"]:
                            reduced["stop"][rl] = spectral_metagene["stop"][rl]
                    metagene_for_fourier = reduced
                results_dict["metrics"]["periodicity_fourier"] = fourier_transform(
                    metagene_for_fourier.copy()
                )

            results_dict["reading_frame_triangle"] = reading_frame_triangle(
                    annotated_read_df
                )
            # Compute entropy-based periodicity on culled read lengths
            culled_for_entropy = read_frame_cull(read_frame_dist, config)
            frame_info_content_dict = rf_info_metric(culled_for_entropy)
            results_dict["read_frame_distribution"] = read_frame_dist
            results_dict["metrics"]["periodicity_information"] =\
                information_metric_cutoff(
                    frame_info_content_dict,
                    config['qc']['read_frame_distribution']['3nt_count_cutoff']
                )

            results_dict["metrics"]["periodicity_information_weighted_score"] = \
                read_frame_information_weighted_score(
                    frame_info_content_dict,
                )

        if annotation:
            print("> mRNA_distribution")
            results_dict["mRNA_distribution"] = mRNA_distribution(
                annotated_read_df
                )

            print("> metagene_profile")
            results_dict["metagene_profile"] = metagene_profile(
                annotated_read_df,
                config["plots"]["metagene_profile"]["distance_target"],
                config["plots"]["metagene_profile"]["distance_range"],
            )

            ###############################################################
            # STOP CODON READTHROUGH RATIO
            ###############################################################
            _stop_meta = results_dict["metagene_profile"].get("stop", {})
            if _stop_meta:
                _before_stop = sum(
                    v for rl_d in _stop_meta.values()
                    for pos, v in rl_d.items()
                    if -30 <= pos < 0
                )
                _after_stop = sum(
                    v for rl_d in _stop_meta.values()
                    for pos, v in rl_d.items()
                    if 0 < pos <= 30
                )
                results_dict["metrics"]["stop_codon_readthrough_ratio"] = (
                    round(_after_stop / _before_stop, 4)
                    if _before_stop > 0 else None
                )
            else:
                results_dict["metrics"]["stop_codon_readthrough_ratio"] = None

            ###############################################################
            # START CODON ENRICHMENT RATIO
            # (proxy for translation-inhibitor treatment; high values
            #  indicate harringtonine/LTM use or stalled initiation)
            ###############################################################
            _start_meta = results_dict["metagene_profile"].get("start", {})
            if _start_meta:
                _near_start = sum(
                    v for rl_d in _start_meta.values()
                    for pos, v in rl_d.items()
                    if -5 <= pos <= 20
                )
                _body_start = sum(
                    v for rl_d in _start_meta.values()
                    for pos, v in rl_d.items()
                    if 30 <= pos <= 50
                )
                results_dict["metrics"]["start_codon_enrichment_ratio"] = (
                    round(_near_start / _body_start, 4)
                    if _body_start > 0 else None
                )
            else:
                results_dict["metrics"]["start_codon_enrichment_ratio"] = None

        else:
            results_dict["mRNA_distribution"] = {}
            results_dict["metagene_profile"] = {"start": {}, "stop": {}}

        #######################################################################
        # UNIFORMITY
        #######################################################################
        # Default + optional uniformity metrics only when annotation-derived metagene exists
        if annotation:
            if should_calculate_metric("uniformity_entropy", config):
                results_dict["metrics"]["uniformity_entropy"] = uniformity_entropy(
                    coding_metagene.copy()
                )

            if should_calculate_metric("uniformity_autocorrelation", config):
                results_dict["metrics"][
                    "uniformity_autocorrelation"
                    ] = uniformity_autocorrelation(
                    coding_metagene.copy()
                )
            if should_calculate_metric("uniformity_theil_index", config):
                results_dict["metrics"][
                    "uniformity_theil_index"
                    ] = uniformity_theil_index(
                    coding_metagene.copy()
                )
            if should_calculate_metric("uniformity_gini_index", config):
                results_dict["metrics"][
                    "uniformity_gini_index"
                    ] = uniformity_gini_index(
                    coding_metagene.copy()
                )

        #######################################################################
        # COVERAGE
        #######################################################################
        print("> cds_coverage_metric")
        _cds_cov = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=config["qc"]["cds_coverage"]["in_frame_coverage"]
            )
        results_dict["metrics"]["CDS_coverage_metric"] = _cds_cov
        results_dict["metrics"]["cds_coverage"] = _cds_cov
        results_dict["metrics"]["CDS_coverage_metric_not_inframe_1read_1000tx"] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=False,
            num_transcripts=1000
            )
        results_dict["metrics"]["cds_coverage_not_inframe_1read_1000tx"] = results_dict["metrics"]["CDS_coverage_metric_not_inframe_1read_1000tx"]
        results_dict["metrics"]["CDS_coverage_metric_not_inframe_100read_100tx"] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=100,
            in_frame_coverage=False,
            num_transcripts=100
            )
        results_dict["metrics"]["cds_coverage_not_inframe_100read_100tx"] = results_dict["metrics"]["CDS_coverage_metric_not_inframe_100read_100tx"]
        results_dict["metrics"]["CDS_coverage_metric_inframe_1read_1000tx"] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=1,
            in_frame_coverage=True,
            num_transcripts=1000
            )
        results_dict["metrics"]["cds_coverage_inframe_1read_1000tx"] = results_dict["metrics"]["CDS_coverage_metric_inframe_1read_1000tx"]
        results_dict["metrics"]["CDS_coverage_metric_inframe_100read_100tx"] = cds_coverage_metric(
            cds_read_df,
            minimum_reads=100,
            in_frame_coverage=True,
            num_transcripts=100
            )
        results_dict["metrics"]["cds_coverage_inframe_100read_100tx"] = results_dict["metrics"]["CDS_coverage_metric_inframe_100read_100tx"]

        #######################################################################
        # RNA REGIONAL SUPPORT
        #######################################################################
        results_dict["metrics"]["ratio_cds:leader"] =\
            region_region_ratio_metric(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region1="CDS",
                region2="five_leader",
            )
        results_dict["metrics"]["ratio_cds:trailer"] =\
            region_region_ratio_metric(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region1="CDS",
                region2="three_trailer",
            )
        results_dict["metrics"]["ratio_leader:trailer"] =\
            region_region_ratio_metric(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region1="five_leader",
                region2="three_trailer",
            )

        results_dict["metrics"]["prop_reads_CDS"] =\
            proportion_of_reads_in_region(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region="CDS",
            )

        results_dict["metrics"]["prop_reads_leader"] =\
            proportion_of_reads_in_region(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region="five_leader",
            )

        results_dict["metrics"]["prop_reads_trailer"] =\
            proportion_of_reads_in_region(
                mRNA_distribution=results_dict["mRNA_distribution"],
                region="three_trailer",
            )
        
        # Ensure read_frame_dist exists when annotation is False but sequence_background is present
        # so downstream culling and periodicity dominance always have input.
        if not annotation:
            read_frame_dist = read_frame_distribution(read_df)
            results_dict["read_frame_distribution"] = read_frame_dist
    else:
        read_frame_dist = read_frame_distribution(read_df)
        results_dict["read_frame_distribution"] = read_frame_dist

    culled_read_frame_dict = read_frame_cull(read_frame_dist, config)

    # Default: periodicity dominance (standard frame preference metric)
    if should_calculate_metric("periodicity_dominance", config):
        results_dict["metrics"]["periodicity_dominance"] = periodicity_dominance(
            culled_read_frame_dict
        )

    # Optional: trips-viz metric
    if should_calculate_metric("periodicity_trips_viz", config):
        results_dict["metrics"][
            "periodicity_trips-viz"
            ] = read_frame_score_trips_viz(
            culled_read_frame_dict)

    ###########################################################################
    # RUST METRIC (requires FASTA + annotation)
    ###########################################################################
    if (
        annotation
        and fasta_dict is not None
        and len(annotation_df) > 0
        and should_calculate_metric("rust_mean_kl_divergence", config)
    ):
        print("> rust (Ribo-seq Unit Step Transformation)")
        from .rust import compute_rust_metrics
        try:
            rust_result = compute_rust_metrics(
                annotated_read_df,
                annotation_df,
                fasta_dict,
            )
            results_dict["rust"] = rust_result
            results_dict["metrics"]["rust_mean_kl_divergence"] = (
                rust_result["mean_kl_divergence"]
            )
        except Exception as exc:
            print(f"Warning: RUST computation failed: {exc}")
            results_dict["rust"] = {}
            results_dict["metrics"]["rust_mean_kl_divergence"] = None

    ###########################################################################
    # RECOMMENDED READ LENGTHS (which lengths to keep for downstream analysis)
    ###########################################################################
    _rec_cfg = config.get("qc", {}).get("recommend", {}) if config else {}
    # --min-periodicity (config["argument"]) overrides the qc.recommend default.
    _cli_min_periodicity = config.get("argument", {}).get("min_periodicity")
    _min_periodicity = (
        float(_cli_min_periodicity)
        if _cli_min_periodicity is not None
        else float(_rec_cfg.get("min_periodicity", 0.5))
    )
    recommended = recommend_read_lengths(
        results_dict["read_frame_distribution"],
        results_dict["read_length_distribution"],
        offsets=computed_offsets or None,
        min_periodicity=_min_periodicity,
        min_read_proportion=float(_rec_cfg.get("min_read_proportion", 0.05)),
    )
    results_dict["recommended_read_lengths"] = recommended
    results_dict["metrics"]["n_recommended_read_lengths"] = (
        recommended["n_recommended"]
    )
    results_dict["metrics"]["recommended_read_proportion"] = (
        recommended["recommended_read_proportion"]
    )

    ###########################################################################
    # ANNOTATION-DERIVED SAMPLE DIAGNOSTICS
    ###########################################################################
    if annotation:
        # Gene-body 5'->3' coverage ramp (relative CDS position metagene)
        print("> gene_body_coverage_ramp")
        ramp = gene_body_coverage_ramp(cds_read_df)
        results_dict["gene_body_coverage"] = ramp
        results_dict["metrics"]["five_prime_ramp_ratio"] = (
            ramp["five_prime_ramp_ratio"]
        )
        results_dict["metrics"]["three_prime_drop_ratio"] = (
            ramp["three_prime_drop_ratio"]
        )

        # Library complexity / saturation (analytic rarefaction)
        print("> library_complexity_curve")
        complexity = library_complexity_curve(annotated_read_df)
        results_dict["library_complexity"] = complexity
        results_dict["metrics"]["marginal_position_discovery_rate"] = (
            complexity["marginal_discovery_rate"]
        )
        results_dict["metrics"]["complexity_distinct_positions"] = (
            complexity["total_distinct_positions"]
        )

        # Library-type classification (elongation / initiation / low quality)
        _dom = results_dict["metrics"].get("periodicity_dominance", {})
        _periodicity = (
            _dom.get("global", 0.0) if isinstance(_dom, dict) else float(_dom)
        )
        _prop_cds = results_dict["metrics"].get("prop_reads_CDS", {})
        _prop_cds = (
            _prop_cds.get("global", 0.0)
            if isinstance(_prop_cds, dict) else float(_prop_cds)
        )
        library_type = classify_library_type(
            periodicity=_periodicity,
            prop_reads_cds=_prop_cds,
            start_codon_enrichment_ratio=results_dict["metrics"].get(
                "start_codon_enrichment_ratio"
            ),
        )
        results_dict["library_type"] = library_type
        print(f"  Library type classified as: {library_type['label']}")

        # FLOSS-style read-length heterogeneity (no FASTA required)
        print("> floss_library_heterogeneity")
        _floss_cfg = config.get("qc", {}).get("floss", {}) if config else {}
        floss = floss_library_heterogeneity(
            annotated_read_df,
            min_reads_per_transcript=int(
                _floss_cfg.get("min_reads_per_transcript", 20)
            ),
            floss_cutoff=float(_floss_cfg.get("cutoff", 0.3)),
        )
        results_dict["floss"] = floss
        results_dict["metrics"]["floss_median"] = floss["floss_median"]
        results_dict["metrics"]["floss_aberrant_transcript_fraction"] = (
            floss["floss_aberrant_transcript_fraction"]
        )

        # A-site codon dwell-times / pause sites (requires FASTA)
        if fasta_dict is not None and len(annotation_df) > 0:
            print("> codon_dwell_times")
            from .rust import compute_codon_dwell_times
            _dwell_cfg = config.get("qc", {}).get("codon_dwell", {}) if config else {}
            try:
                dwell = compute_codon_dwell_times(
                    annotated_read_df,
                    annotation_df,
                    fasta_dict,
                    elongation_5_nt=int(_dwell_cfg.get("elongation_5_nt", 15)),
                    elongation_3_nt=int(_dwell_cfg.get("elongation_3_nt", 15)),
                )
            except Exception as exc:
                print(f"Warning: codon dwell-time computation failed: {exc}")
                dwell = {}
            results_dict["codon_dwell_times"] = dwell
            results_dict["metrics"]["codon_dwell_cv"] = dwell.get("codon_dwell_cv")
            results_dict["metrics"]["codon_dwell_p90_p10"] = dwell.get(
                "codon_dwell_p90_p10"
            )
            results_dict["metrics"]["proline_dwell"] = dwell.get("proline_dwell")
            results_dict["metrics"]["cga_dwell"] = dwell.get("cga_dwell")

    return results_dict


def select_read_lengths_for_global(read_frame_dist: dict, read_length_distribution: dict, config: dict) -> list:
    """Select read lengths to contribute to global spectral metrics.

    Supports either top-N by mass or a cumulative mass threshold.
    """
    opts = config.get("qc", {}).get("read_frame_distribution", {})
    top_n = opts.get("periodicity_top_n")
    mass_threshold = opts.get("periodicity_mass_threshold")
    if not top_n and not mass_threshold:
        return sorted(read_frame_dist.keys())

    # Build (len, mass) from read_length_distribution
    pairs = sorted(read_length_distribution.items(), key=lambda kv: kv[1], reverse=True)
    selected = []
    if top_n:
        selected = [int(k) for k, _ in pairs[: int(top_n)]]
    elif mass_threshold:
        total = sum(v for _, v in pairs)
        acc = 0
        for k, v in pairs:
            selected.append(int(k))
            acc += v
            if acc / total >= float(mass_threshold):
                break
    # Intersect with available keys in read_frame_dist
    return [l for l in selected if l in read_frame_dist]
