import pandas as pd
from RiboMetric.modules import a_site_calculation


def make_df(refs, counts=None):
    return pd.DataFrame({
        'read_name': [f'r{i}' for i in range(len(refs))],
        'reference_start': refs,
        'read_length': [28]*len(refs),
        'reference_name': ['tx1']*len(refs),
        'first_dinucleotide': ['AA']*len(refs),
        'last_dinucleotide': ['TT']*len(refs),
        'count': counts if counts is not None else [1]*len(refs),
    })


def test_asite_global_str_reference_start():
    df = make_df(['10', '20', '30'])
    out = a_site_calculation(df, offset_type='global', global_offset=15)
    assert (out['a_site'] == pd.Series([25, 35, 45])).all()


def test_asite_read_specific_str_offset(tmp_path):
    df = make_df([10, 20, 30])
    offsets = tmp_path / 'offs.tsv'
    offsets.write_text('r0\t1\n' 'r1\t2\n' 'r2\t3\n')
    out = a_site_calculation(df, offset_type='read_specific', offset_file=str(offsets))
    assert (out['a_site'] == pd.Series([11, 22, 33])).all()


def test_asite_read_length_offsets(tmp_path):
    df = make_df([10, 20, 30])
    # read_length all 28, provide mapping file
    rl = tmp_path / 'rl.tsv'
    rl.write_text('28\t5\n')
    out = a_site_calculation(df, offset_type='read_length', offset_file=str(rl))
    assert (out['a_site'] == pd.Series([15, 25, 35])).all()


def test_asite_calculate_uses_variable_offset():
    # Just ensure it runs and produces numeric a_site
    df = make_df([10, 20, 30])
    out = a_site_calculation(df, offset_type='calculate')
    assert pd.api.types.is_numeric_dtype(out['a_site'])
