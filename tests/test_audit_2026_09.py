"""Regression tests for the defects found in the 2026-09-07 v1.4.3 audit.

Each test names the observable symptom it prevents, so a future change that
reintroduces one fails with an explanation rather than a bare assertion.
"""

import pandas as pd
import pytest

from RiboMetric.metrics import (
    DOMINANCE_MIN_READS,
    periodicity_dominance,
    recommend_read_lengths,
)
from RiboMetric.modules import (
    a_site_calculation_variable_offset,
    sum_mRNA_distribution,
)
from RiboMetric.plots import plot_mRNA_read_breakdown


# --- offsets ---------------------------------------------------------------


def test_default_offset_is_clipped_to_the_read_length():
    """A read length missing from offset_dict used to receive the bare
    default_offset, which is sanitised once against an effectively unbounded
    read length. A 15 nt read was given offset 15, placing the A-site one base
    past the read's 3' end -- and it appeared that way in the offsets TSV.
    """
    read_df = pd.DataFrame({
        "read_length": [15, 16, 30],
        "reference_start": [0, 0, 0],
    })
    out = a_site_calculation_variable_offset(
        read_df, offset_dict={30: 14}, default_offset=15
    )
    applied = dict(zip(out["read_length"], out["offset"]))

    assert applied[30] == 14, "explicit offsets must be preserved"
    for read_length in (15, 16):
        assert applied[read_length] < read_length, (
            f"offset {applied[read_length]} on a {read_length} nt read puts "
            "the A-site at or past the 3' end"
        )
    # The read-length ceiling (2/3 by default) applies to the fallback too.
    assert applied[15] <= 10
    assert applied[16] <= 10


def test_explicit_offsets_are_not_silently_revalidated():
    """validate_offsets=False must still mean 'trust the caller'."""
    read_df = pd.DataFrame({"read_length": [30], "reference_start": [0]})
    out = a_site_calculation_variable_offset(
        read_df, offset_dict={30: 25}, default_offset=15, validate_offsets=False
    )
    assert out["offset"].iloc[0] == 25


# --- per-read-length periodicity -------------------------------------------


def test_low_count_read_lengths_do_not_report_a_dominance():
    """Read lengths backed by one or two reads reported a dominant-frame
    fraction of exactly 1.000, which then reached cohort tables and was drawn
    as a full-height bar in the recommended-read-lengths plot."""
    read_frame_dict = {
        20: {0: 1, 1: 0, 2: 0},          # n=1  -> pure noise
        21: {0: 2, 1: 0, 2: 0},          # n=2  -> pure noise
        30: {0: 900, 1: 60, 2: 40},      # n=1000 -> real
    }
    dominance = periodicity_dominance(read_frame_dict)

    assert 20 not in dominance
    assert 21 not in dominance
    assert dominance[30] == pytest.approx(0.9)


def test_global_dominance_still_uses_every_read():
    """The min_reads floor must not move the gated Tier-1 value."""
    read_frame_dict = {
        20: {0: 1, 1: 0, 2: 0},
        30: {0: 900, 1: 60, 2: 40},
    }
    with_floor = periodicity_dominance(read_frame_dict, min_reads=100)
    without_floor = periodicity_dominance(read_frame_dict, min_reads=0)
    assert with_floor["global"] == without_floor["global"]
    assert (
        with_floor["global_by_read_length_max"]
        == without_floor["global_by_read_length_max"]
    )


def test_recommendations_skip_read_lengths_with_too_few_frame_reads():
    read_frame_dict = {
        20: {0: 1, 1: 0, 2: 0},
        30: {0: 900, 1: 60, 2: 40},
    }
    rld = {20: 1, 30: 1000}
    rec = recommend_read_lengths(read_frame_dict, rld)

    assert 20 not in rec["by_read_length"]
    assert rec["by_read_length"][30]["n_frame_reads"] == 1000
    assert rec["recommended_lengths"] == [30]


def test_dominance_floor_default_is_shared():
    assert DOMINANCE_MIN_READS == 100


# --- mRNA distribution ------------------------------------------------------


MRNA_DIST = {
    "global": {"CDS": 300, "five_leader": 100},
    28: {"CDS": 100, "five_leader": 40},
    29: {"CDS": 200, "five_leader": 60},
}


def test_sum_mrna_distribution_does_not_double_count_global():
    """mRNA_distribution carries a 'global' entry that is already the sum over
    read lengths. Including it doubled every absolute count (invisible in the
    default proportional view because the doubling cancels)."""
    config = {"plots": {"mRNA_distribution": {"absolute_counts": True}}}
    assert sum_mRNA_distribution(MRNA_DIST, config) == {
        "CDS": 300, "five_leader": 100,
    }


def test_sum_mrna_distribution_proportions_still_sum_to_one():
    config = {"plots": {"mRNA_distribution": {"absolute_counts": False}}}
    result = sum_mRNA_distribution(MRNA_DIST, config)
    assert sum(result.values()) == pytest.approx(1.0)
    assert result["CDS"] == pytest.approx(0.75)


def test_read_breakdown_plot_excludes_the_global_aggregate(monkeypatch):
    """'global' was passed to plotly as the first x category alongside numeric
    read lengths, making the aggregate the tallest point on the chart and
    doubling the normalisation denominator."""
    import RiboMetric.plots as plots

    captured = {}

    class _FakeFig:
        def add_trace(self, trace):
            captured.setdefault("x", trace.x)
            captured.setdefault("traces", []).append(trace)

        def update_layout(self, **kwargs):
            captured["layout"] = kwargs

    monkeypatch.setattr(plots.go, "Figure", lambda: _FakeFig())
    monkeypatch.setattr(plots, "plotly_to_html", lambda fig: "")
    monkeypatch.setattr(plots, "plotly_to_image", lambda fig, w, h: "")

    config = {
        "plots": {
            "mRNA_read_breakdown": {"absolute_counts": False},
            "font_family": "Arial",
            "base_color": "#000",
            "image_size": [720, 370],
        }
    }
    plot_mRNA_read_breakdown(MRNA_DIST, config)

    assert "global" not in captured["x"]
    assert list(captured["x"]) == [28, 29]
    # Proportions must be normalised over the read lengths only.
    total = sum(sum(t.y) for t in captured["traces"])
    assert total == pytest.approx(1.0)
    assert captured["layout"]["title"] == "mRNA Reads Breakdown over Read Length"


def test_global_offset_is_clipped_to_the_read_length():
    """The offset_type='global' path assigned the global offset to every read
    with no bounds check at all, so short reads got an A-site past their own
    3' end -- and that value was what the offsets audit TSV reported."""
    from RiboMetric.modules import a_site_calculation

    read_df = pd.DataFrame({
        "read_length": [15, 16, 30],
        "reference_start": [100, 100, 100],
    })
    out = a_site_calculation(read_df, offset_type="global", global_offset=15)
    applied = dict(zip(out["read_length"], out["offset"]))

    assert applied[30] == 15, "a valid global offset must be used as given"
    for read_length in (15, 16):
        assert applied[read_length] < read_length
        assert (
            out.loc[out["read_length"] == read_length, "a_site"].iloc[0]
            < 100 + read_length
        )
