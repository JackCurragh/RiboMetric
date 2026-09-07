import pandas as pd

from RiboMetric.modules import (
    a_site_calculation,
    a_site_calculation_variable_offset,
    asite_calculation_per_readlength,
)


def make_df(refs, counts=None):
    return pd.DataFrame(
        {
            "read_name": [f"r{i}" for i in range(len(refs))],
            "reference_start": refs,
            "read_length": [28] * len(refs),
            "reference_name": ["tx1"] * len(refs),
            "first_dinucleotide": ["AA"] * len(refs),
            "last_dinucleotide": ["TT"] * len(refs),
            "count": counts if counts is not None else [1] * len(refs),
        }
    )


def test_asite_global_str_reference_start():
    df = make_df(["10", "20", "30"])
    out = a_site_calculation(df, offset_type="global", global_offset=15)
    assert (out["a_site"] == pd.Series([25, 35, 45])).all()


def test_asite_read_specific_str_offset(tmp_path):
    df = make_df([10, 20, 30])
    offsets = tmp_path / "offs.tsv"
    offsets.write_text("r0\t1\n" "r1\t2\n" "r2\t3\n")
    out = a_site_calculation(df, offset_type="read_specific", offset_file=str(offsets))
    assert (out["a_site"] == pd.Series([11, 22, 33])).all()


def test_asite_read_length_offsets(tmp_path):
    df = make_df([10, 20, 30])
    # read_length all 28, provide mapping file
    rl = tmp_path / "rl.tsv"
    rl.write_text("28\t5\n")
    out = a_site_calculation(df, offset_type="read_length", offset_file=str(rl))
    assert (out["a_site"] == pd.Series([15, 25, 35])).all()


def test_validated_variable_offsets_replace_implausible_offsets():
    df = pd.DataFrame(
        {
            "read_name": ["r0", "r1", "r2"],
            "reference_start": [10, 20, 30],
            "read_length": [27, 31, 32],
            "reference_name": ["tx1", "tx1", "tx1"],
            "first_dinucleotide": ["AA", "AA", "AA"],
            "last_dinucleotide": ["TT", "TT", "TT"],
            "count": [1, 1, 1],
        }
    )

    out = a_site_calculation_variable_offset(
        df,
        {27: 20, 31: 31, 32: 44},
        validate_offsets=True,
    )

    assert out["offset"].tolist() == [15, 15, 15]
    assert out["a_site"].tolist() == [25, 35, 45]


def test_tripsviz_offset_calculation_falls_back_for_extreme_offsets():
    df = pd.DataFrame(
        {
            "read_name": ["r0", "r1", "r2"],
            "reference_start": [56, 56, 88],
            "read_length": [32, 32, 32],
            "cds_start": [100, 100, 100],
            "cds_end": [400, 400, 400],
            "mapq": [255, 255, 255],
            "count": [10, 5, 1],
        }
    )

    offsets = asite_calculation_per_readlength(df, method="tripsviz")

    assert offsets[32] == 15


def test_automatic_offset_methods_share_a_site_target():
    df = pd.DataFrame(
        {
            "read_name": [f"r{i}" for i in range(12)],
            "reference_start": [-11] * 10 + [-13, -10],
            "read_length": [30] * 12,
            "cds_start": [0] * 12,
            "cds_end": [300] * 12,
            "mapq": [255] * 12,
            "count": [1] * 12,
        }
    )

    offsets = {
        method: asite_calculation_per_readlength(
            df,
            method=method,
            offset_target="a_site",
        )[30]
        for method in ["ribowaltz", "changepoint", "tripsviz"]
    }

    assert offsets == {
        "ribowaltz": 14,
        "changepoint": 14,
        "tripsviz": 14,
    }


def test_automatic_offset_methods_share_p_site_target():
    df = pd.DataFrame(
        {
            "read_name": [f"r{i}" for i in range(12)],
            "reference_start": [-11] * 10 + [-13, -10],
            "read_length": [30] * 12,
            "cds_start": [0] * 12,
            "cds_end": [300] * 12,
            "mapq": [255] * 12,
            "count": [1] * 12,
        }
    )

    offsets = {
        method: asite_calculation_per_readlength(
            df,
            method=method,
            offset_target="p_site",
            default_offset=12,
        )[30]
        for method in ["ribowaltz", "changepoint", "tripsviz"]
    }

    assert offsets == {
        "ribowaltz": 11,
        "changepoint": 11,
        "tripsviz": 11,
    }


def test_asite_calculate_uses_variable_offset():
    # Just ensure it runs and produces numeric a_site
    df = make_df([10, 20, 30])
    out = a_site_calculation(df, offset_type="calculate")
    assert pd.api.types.is_numeric_dtype(out["a_site"])
