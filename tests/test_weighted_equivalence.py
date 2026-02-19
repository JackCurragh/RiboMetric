"""
Weighted vs expanded equivalence tests on small fixtures.

These ensure that weighted implementations produce the same results
as explicit row expansion for key modules.
"""

import pandas as pd

from RiboMetric.modules import (
    a_site_calculation,
    read_frame_distribution,
    annotate_reads,
    assign_mRNA_category,
    mRNA_distribution,
    metagene_profile,
)


def expand(df: pd.DataFrame) -> pd.DataFrame:
    """Return a row-expanded copy; ensure multiples of ten to match fixtures.

    Test fixtures often use counts that are multiples of 10/100; the expanded
    dataframe should reflect exact multiplicities. This helper keeps the
    reference path explicit and leaves weighted path to use 'count' values.
    """
    return df.loc[df.index.repeat(df["count"].astype(int))].reset_index(drop=True)


def test_read_frame_distribution_weighted_equivalence(sample_read_df, sample_annotation_df):
    # Ensure weighted path truly uses weights, and expanded path has no weights
    wdf = sample_read_df.copy()
    edf = expand(sample_read_df).drop(columns=["count"]) if "count" in expand(sample_read_df).columns else expand(sample_read_df)
    w = a_site_calculation(wdf)
    e = a_site_calculation(edf)

    dw = read_frame_distribution(w)
    de = read_frame_distribution(e)
    assert dw == de


def test_mrna_distribution_weighted_equivalence(sample_read_df, sample_annotation_df):
    wdf = sample_read_df.copy()
    edf = expand(sample_read_df).drop(columns=["count"]) if "count" in expand(sample_read_df).columns else expand(sample_read_df)
    w = a_site_calculation(wdf)
    w = annotate_reads(w, sample_annotation_df)
    w = assign_mRNA_category(w)

    e = a_site_calculation(edf)
    e = annotate_reads(e, sample_annotation_df)
    e = assign_mRNA_category(e)

    dw = mRNA_distribution(w)
    de = mRNA_distribution(e)
    assert dw == de


def test_metagene_profile_weighted_equivalence(sample_read_df, sample_annotation_df):
    wdf = sample_read_df.copy()
    edf = expand(sample_read_df).drop(columns=["count"]) if "count" in expand(sample_read_df).columns else expand(sample_read_df)
    w = a_site_calculation(wdf)
    w = annotate_reads(w, sample_annotation_df)

    e = a_site_calculation(edf)
    e = annotate_reads(e, sample_annotation_df)

    pw = metagene_profile(w)
    pe = metagene_profile(e)
    assert pw == pe
