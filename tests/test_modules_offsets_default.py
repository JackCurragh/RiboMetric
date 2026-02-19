import pandas as pd
from RiboMetric.modules import a_site_calculation_variable_offset


def test_asite_variable_offset_default_is_stable():
    df = pd.DataFrame(
        {
            "reference_start": [100, 200, 300],
            "read_length": pd.Series([28, 29, 30], dtype="category"),
        }
    )

    out1 = a_site_calculation_variable_offset(df.copy(), None)
    out2 = a_site_calculation_variable_offset(df.copy(), None)

    assert out1["a_site"].tolist() == out2["a_site"].tolist() == [115, 215, 315]

