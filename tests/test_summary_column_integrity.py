"""Appending a new sample must not relabel its measurements."""
import csv
import pytest
from RiboMetric.results_output import generate_summary_tsv


def append(path,name,metrics):
    generate_summary_tsv({'metrics':metrics},{},name,name=str(path))

def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def test_reordered_metric_dictionary_preserves_column_identity(tmp_path):
    p=tmp_path/'summary.tsv'
    append(p,'first',{'periodicity_dominance':.91,'uniformity_entropy':.82})
    append(p,'second',{'uniformity_entropy':.21,'periodicity_dominance':.74})
    row=read(p)[1]
    assert float(row['periodicity_dominance'])==.74
    assert float(row['uniformity_entropy'])==.21

@pytest.mark.parametrize('second',[{'b':.2},{'a':None,'b':.2}])
def test_missing_or_null_metric_does_not_shift_following_values(tmp_path,second):
    p=tmp_path/'summary.tsv';append(p,'first',{'a':.8,'b':.9});append(p,'second',second)
    row=read(p)[1]
    assert row['a']=='' and float(row['b'])==.2


def test_new_metric_is_rejected_without_modifying_existing_file(tmp_path):
    p=tmp_path/'summary.tsv';append(p,'first',{'a':.8});before=p.read_bytes()
    with pytest.raises(ValueError,match='new metrics'):append(p,'second',{'a':.7,'b':.3})
    assert p.read_bytes()==before


def test_null_column_exists_from_first_row(tmp_path):
    p=tmp_path/'summary.tsv';append(p,'first',{'a':None,'b':.9});append(p,'second',{'a':.7,'b':.3})
    rows=read(p);assert rows[0]['a']=='' and float(rows[1]['a'])==.7

@pytest.mark.parametrize('header',['garbage\n','sample\ttimestamp\tmode\ttotal_reads\ta\ta\n'])
def test_invalid_saved_header_is_not_appended(tmp_path,header):
    p=tmp_path/'summary.tsv';p.write_text(header)
    with pytest.raises(ValueError):append(p,'second',{'a':.7})
    assert p.read_text()==header


def test_empty_file_receives_header(tmp_path):
    p=tmp_path/'summary.tsv';p.touch();append(p,'first',{'a':.8})
    assert read(p)[0]['sample']=='first'


def test_metadata_collision_is_rejected(tmp_path):
    p=tmp_path/'summary.tsv'
    with pytest.raises(ValueError):append(p,'first',{'sample':.8})
    assert not p.exists()
