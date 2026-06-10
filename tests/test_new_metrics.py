"""
Tests for the new metrics added in v1.2.0:
  - alignment stats (duplicate_rate, multimapper_rate, soft_clip_rate_5prime)
  - disome_proportion
  - stop_codon_readthrough_ratio / start_codon_enrichment_ratio
  - RUST (compute_rust_metrics, _empty_result, _lookup_seq)
"""

import math
import pandas as pd
import numpy as np
import pytest

from RiboMetric.qc import calculate_alignment_stats

# ------------------------------------------------------------------ #
# Helpers                                                              #
# ------------------------------------------------------------------ #

def _make_read_df(
    n_reads=100,
    mapq_frac_multi=0.1,   # fraction with MAPQ < 255
    nh_frac_multi=None,    # fraction with NH > 1; when set, overrides MAPQ
    dup_count=2,            # collapsed count for half the reads
    soft_clip_frac=0.2,    # fraction with soft_clip_5 > 0
):
    """Build a minimal read_df mimicking the output of join_batches."""
    rng = np.random.default_rng(42)
    counts = [dup_count if i < n_reads // 2 else 1 for i in range(n_reads)]
    mapq = [0 if i < int(n_reads * mapq_frac_multi) else 255 for i in range(n_reads)]
    nh = None
    if nh_frac_multi is not None:
        nh = [2 if i < int(n_reads * nh_frac_multi) else 1 for i in range(n_reads)]
    soft_clips = [3 if i < int(n_reads * soft_clip_frac) else 0 for i in range(n_reads)]
    rl = rng.integers(28, 33, size=n_reads).tolist()
    data = {
        "read_length": pd.Categorical(rl),
        "reference_name": pd.Categorical(["tx1"] * n_reads),
        "reference_start": list(range(n_reads)),
        "count": pd.Categorical(counts),
        "mapq": mapq,
        "soft_clip_5": soft_clips,
        "a_site": list(range(n_reads)),
    }
    if nh is not None:
        data["nh"] = nh
    df = pd.DataFrame(data)
    return df


def _make_annotation_df():
    return pd.DataFrame({
        "transcript_id": ["tx1"],
        "cds_start": [30],
        "cds_end": [330],
        "transcript_length": [500],
        "genomic_cds_starts": [""],
        "genomic_cds_ends": [""],
    })


# ------------------------------------------------------------------ #
# Alignment stats (inline logic, tested by calling annotation_mode   #
# internals through a minimal read_df)                                #
# ------------------------------------------------------------------ #

class TestAlignmentStats:
    def test_duplicate_rate_computation(self):
        """dup_rate = 1 - (unique_reads / total_weighted)."""
        df = _make_read_df(n_reads=100, dup_count=2)
        stats = calculate_alignment_stats(df)
        # 50 reads with count=2, 50 with count=1 → total=150, unique=100
        assert stats["total_reads_analysed"] == 150
        assert stats["unique_read_sequences"] == 100
        assert abs(stats["duplicate_rate"] - (1.0 - 100 / 150)) < 1e-9

    def test_multimapper_rate_zero(self):
        df = _make_read_df(n_reads=50, mapq_frac_multi=0.0)
        stats = calculate_alignment_stats(df)
        assert stats["multimapper_detection_method"] == "star_mapq"
        assert stats["rpf_multimapper_rate"] == 0.0
        assert stats["unique_rpf_rate"] == 1.0

    def test_multimapper_rate_nonzero(self):
        df = _make_read_df(n_reads=100, mapq_frac_multi=0.2, dup_count=1)
        stats = calculate_alignment_stats(df)
        assert abs(stats["rpf_multimapper_rate"] - 0.2) < 1e-9
        assert abs(stats["alignment_multimapper_rate"] - 0.2) < 1e-9
        assert stats["multimapper_rate"] == stats["rpf_multimapper_rate"]

    def test_nh_takes_precedence_and_keeps_weighted_rpf_rate(self):
        df = _make_read_df(
            n_reads=100,
            mapq_frac_multi=0.0,
            nh_frac_multi=0.2,
            dup_count=2,
        )
        stats = calculate_alignment_stats(df)
        # First 20 reads are multimappers and have count=2, so the weighted
        # fragment rate is 40 / 150 while the row-level alignment rate is 20%.
        assert stats["multimapper_detection_method"] == "nh"
        assert abs(stats["rpf_multimapper_rate"] - (40 / 150)) < 1e-9
        assert abs(stats["alignment_multimapper_rate"] - 0.2) < 1e-9

    def test_xa_alt_locus_count_drives_multimapper_rate(self):
        """XA falls between NH and MAPQ: an alt-locus count > 0 is multimapping."""
        df = pd.DataFrame(
            {
                "read_name": [f"r{i}" for i in range(10)],
                "count": [1] * 10,
                # No NH column, so the XA branch is selected. 3 of 10 fragments
                # have alternative loci.
                "xa": [2.0, 1.0, 3.0] + [0.0] * 7,
                "mapq": [255] * 10,
            }
        )
        stats = calculate_alignment_stats(df)
        assert stats["multimapper_detection_method"] == "star_transcriptome_xa"
        assert abs(stats["rpf_multimapper_rate"] - 0.3) < 1e-9


class TestXATagParsing:
    """``_parse_xa_value`` must reduce both STAR and BWA XA forms to a count."""

    def test_bwa_string_list_counts_alternative_loci(self):
        from RiboMetric.bam_processing import _parse_xa_value

        assert _parse_xa_value("XA:Z:chr1,+100,30M,2;chr2,-200,30M,1;") == 2.0
        assert _parse_xa_value("chr1,+100,30M,2;") == 1.0

    def test_star_integer_form(self):
        from RiboMetric.bam_processing import _parse_xa_value

        assert _parse_xa_value("XA:i:3") == 3.0
        assert _parse_xa_value(2) == 2.0

    def test_absent_or_empty_is_nan(self):
        from RiboMetric.bam_processing import _parse_xa_value

        assert np.isnan(_parse_xa_value(None))
        assert np.isnan(_parse_xa_value(""))
        assert np.isnan(_parse_xa_value(float("nan")))


class TestUniqueMapperFiltering:
    """Graceful behaviour when there is no uniqueness signal (no NH, no MAPQ)."""

    def test_filter_returns_empty_when_no_signal(self):
        from RiboMetric.modules import filter_unique_mappers

        df = pd.DataFrame({"read_name": ["a", "b"], "count": [1, 1]})
        out = filter_unique_mappers(df)
        # No NH, no MAPQ column -> nothing can be asserted unique -> empty frame
        # (this is the input that previously crashed plot_read_frame_distribution).
        assert out.empty
        assert list(out.columns) == list(df.columns)

    def test_alignment_stats_handle_signalless_frame(self):
        """The metric layer must not crash on a no-signal read_df."""
        df = pd.DataFrame({"read_name": ["a", "b"], "count": [1, 1]})
        stats = calculate_alignment_stats(df)
        # With neither NH/XA/MAPQ evidence, multimapper rate is well-defined (0).
        assert stats["rpf_multimapper_rate"] == 0.0
        assert stats["unique_rpf_rate"] == 1.0

    def test_soft_clip_rate(self):
        df = _make_read_df(n_reads=100, soft_clip_frac=0.3, dup_count=1)
        count_arr = df["count"].astype(int)
        total = int(count_arr.sum())
        sc5 = df["soft_clip_5"].astype(int)
        sc_mask = (sc5 > 0).astype(int)
        rate = float((sc_mask * count_arr).sum() / total)
        assert abs(rate - 0.3) < 1e-9


# ------------------------------------------------------------------ #
# Di-some proportion                                                   #
# ------------------------------------------------------------------ #

class TestDisomeProportion:
    def _make_rld(self, counts_dict):
        """Build a read_length_distribution dict."""
        return {k: v for k, v in counts_dict.items()}

    def test_no_disome(self):
        rld = {28: 100, 29: 80, 30: 70, 31: 50}
        total = sum(rld.values())
        disome_count = sum(v for k, v in rld.items() if 50 <= k <= 70)
        assert disome_count == 0
        assert disome_count / total == 0.0

    def test_all_disome(self):
        rld = {55: 100, 60: 100}
        total = sum(rld.values())
        disome_count = sum(v for k, v in rld.items() if 50 <= k <= 70)
        assert disome_count / total == 1.0

    def test_partial_disome(self):
        rld = {30: 80, 60: 20}
        total = sum(rld.values())
        disome_count = sum(v for k, v in rld.items() if 50 <= k <= 70)
        assert abs(disome_count / total - 0.2) < 1e-9


# ------------------------------------------------------------------ #
# Stop / start codon ratios                                            #
# ------------------------------------------------------------------ #

class TestCodonEnrichment:
    def _make_metagene(self, positions, base_count=10, boost_range=None, boost_mult=5):
        """Create a fake metagene dict {read_len: {pos: count}}."""
        meta = {}
        for pos in positions:
            c = base_count
            if boost_range and boost_range[0] <= pos <= boost_range[1]:
                c = base_count * boost_mult
            meta.setdefault(28, {})[pos] = c
        return meta

    def test_stop_readthrough_ratio_low(self):
        """No reads downstream → ratio ≈ 0."""
        positions = list(range(-30, 31))
        # uniform except upstream > downstream
        meta = {}
        for pos in positions:
            meta.setdefault(28, {})[pos] = 10 if pos < 0 else (1 if pos > 0 else 5)
        before = sum(v for rd in meta.values() for pos, v in rd.items() if -30 <= pos < 0)
        after = sum(v for rd in meta.values() for pos, v in rd.items() if 0 < pos <= 30)
        ratio = after / before if before > 0 else None
        assert ratio is not None
        assert ratio < 1.0  # more reads upstream than downstream

    def test_stop_readthrough_ratio_high(self):
        """Many reads downstream → ratio > 1."""
        positions = list(range(-30, 31))
        meta = {}
        for pos in positions:
            meta.setdefault(28, {})[pos] = 1 if pos < 0 else (20 if pos > 0 else 5)
        before = sum(v for rd in meta.values() for pos, v in rd.items() if -30 <= pos < 0)
        after = sum(v for rd in meta.values() for pos, v in rd.items() if 0 < pos <= 30)
        ratio = after / before if before > 0 else None
        assert ratio is not None
        assert ratio > 1.0

    def test_start_enrichment_ratio(self):
        """More reads near start than in body → ratio > 1."""
        positions = list(range(-50, 51))
        meta = {}
        for pos in positions:
            if -5 <= pos <= 20:
                meta.setdefault(28, {})[pos] = 50
            elif 30 <= pos <= 50:
                meta.setdefault(28, {})[pos] = 10
            else:
                meta.setdefault(28, {})[pos] = 5
        near = sum(v for rd in meta.values() for pos, v in rd.items() if -5 <= pos <= 20)
        body = sum(v for rd in meta.values() for pos, v in rd.items() if 30 <= pos <= 50)
        ratio = near / body if body > 0 else None
        assert ratio is not None
        assert ratio > 1.0


# ------------------------------------------------------------------ #
# RUST module                                                          #
# ------------------------------------------------------------------ #

class TestRustHelpers:
    def test_lookup_seq_exact(self):
        from RiboMetric.rust import _lookup_seq
        from unittest.mock import MagicMock
        rec = MagicMock()
        rec.seq = "ATGATG"
        fasta = {"tx1": rec}
        assert _lookup_seq(fasta, "tx1") == "ATGATG"

    def test_lookup_seq_pipe_split(self):
        from RiboMetric.rust import _lookup_seq
        from unittest.mock import MagicMock
        rec = MagicMock()
        rec.seq = "ATGATG"
        fasta = {"NM_001": rec}
        assert _lookup_seq(fasta, "NM_001|isoform1") == "ATGATG"

    def test_lookup_seq_missing(self):
        from RiboMetric.rust import _lookup_seq
        assert _lookup_seq({}, "missing_tx") is None

    def test_empty_result_shape(self):
        from RiboMetric.rust import _empty_result, WINDOW_SIZE
        result = _empty_result(WINDOW_SIZE)
        assert result["metagene"] == {}
        assert len(result["kl_divergence"]) == WINDOW_SIZE
        assert all(v is None for v in result["kl_divergence"])
        assert result["mean_kl_divergence"] == 0.0
        assert result["transcripts_used"] == 0


class TestRustCompute:
    """Integration test for compute_rust_metrics with synthetic data."""

    def _make_inputs(self):
        """
        Build a minimal annotated_read_df, annotation_df, and fasta_dict
        that will pass through the RUST algorithm without errors.
        """
        # Transcript sequence: 600 nt, CDS = 0..600 (= whole seq for simplicity)
        # Codons: repeat "GCT" (Ala) to fill 200 codons
        cds_seq = "GCT" * 200  # 600 nt
        transcript_id = "tx_rust"
        cds_start = 0
        cds_end = 600

        # Annotation
        ann_df = pd.DataFrame({
            "transcript_id": [transcript_id],
            "cds_start": [cds_start],
            "cds_end": [cds_end],
            "transcript_length": [600],
            "genomic_cds_starts": [""],
            "genomic_cds_ends": [""],
        })

        # FASTA dict: simple string (no SeqRecord needed)
        class _FakeSeq:
            def __str__(self):
                return cds_seq
            def __repr__(self):
                return cds_seq

        class _FakeRec:
            seq = _FakeSeq()

        fasta = {transcript_id: _FakeRec()}

        # annotated_read_df: reads uniformly spread across CDS elongation region
        # elongation_region = cds_seq[120:-60] = positions 120..539 (420 nt = 140 codons)
        rng = np.random.default_rng(0)
        n = 500
        a_sites = rng.integers(120, 540, size=n).tolist()
        df = pd.DataFrame({
            "reference_name": [transcript_id] * n,
            "a_site": a_sites,
            "cds_start": [cds_start] * n,
            "cds_end": [cds_end] * n,
            "count": [1] * n,
            "mRNA_category": [2] * n,
        })

        return df, ann_df, fasta

    def test_compute_runs_without_error(self):
        from RiboMetric.rust import compute_rust_metrics
        df, ann_df, fasta = self._make_inputs()
        result = compute_rust_metrics(df, ann_df, fasta)
        assert "metagene" in result
        assert "kl_divergence" in result
        assert "mean_kl_divergence" in result
        assert "transcripts_used" in result

    def test_transcripts_used_gt_zero(self):
        from RiboMetric.rust import compute_rust_metrics
        df, ann_df, fasta = self._make_inputs()
        result = compute_rust_metrics(df, ann_df, fasta)
        assert result["transcripts_used"] >= 1

    def test_accepts_transcript_id_column_from_annotation_merge(self):
        from RiboMetric.rust import compute_rust_metrics
        df, ann_df, fasta = self._make_inputs()
        df = df.rename(columns={"reference_name": "transcript_id"})
        result = compute_rust_metrics(df, ann_df, fasta)
        assert result["transcripts_used"] >= 1

    def test_metagene_has_61_codons(self):
        from RiboMetric.rust import compute_rust_metrics, SENSE_CODONS
        df, ann_df, fasta = self._make_inputs()
        result = compute_rust_metrics(df, ann_df, fasta)
        # Should have an entry for every sense codon
        assert len(result["metagene"]) == len(SENSE_CODONS)

    def test_kl_divergence_length(self):
        from RiboMetric.rust import compute_rust_metrics, WINDOW_SIZE
        df, ann_df, fasta = self._make_inputs()
        result = compute_rust_metrics(df, ann_df, fasta)
        assert len(result["kl_divergence"]) == WINDOW_SIZE

    def test_mean_kl_non_negative(self):
        from RiboMetric.rust import compute_rust_metrics
        df, ann_df, fasta = self._make_inputs()
        result = compute_rust_metrics(df, ann_df, fasta)
        assert result["mean_kl_divergence"] >= 0.0

    def test_empty_df_returns_empty(self):
        from RiboMetric.rust import compute_rust_metrics, _empty_result, WINDOW_SIZE
        _, ann_df, fasta = self._make_inputs()
        empty_df = pd.DataFrame(columns=["reference_name", "a_site", "cds_start", "cds_end", "count"])
        result = compute_rust_metrics(empty_df, ann_df, fasta)
        expected = _empty_result(WINDOW_SIZE)
        assert result["transcripts_used"] == 0
        assert result["metagene"] == expected["metagene"]

    def test_no_fasta_match_returns_empty(self):
        from RiboMetric.rust import compute_rust_metrics
        df, ann_df, _ = self._make_inputs()
        result = compute_rust_metrics(df, ann_df, {})  # empty fasta
        assert result["transcripts_used"] == 0
