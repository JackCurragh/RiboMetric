"""Required evidence must not disappear from explicit QC policy evaluation."""
import copy
import json
import math
import numpy as np
import pytest
import RiboMetric.results_output as output
from RiboMetric.results_output import evaluate_qc_status

P = {'periodicity_dominance': {'pass': .7, 'warn': .5}}

@pytest.mark.parametrize('metrics,reason', [
    ({}, 'missing_metric'),
    ({'periodicity_dominance': None}, 'unavailable_metric'),
    ({'periodicity_dominance': 'not-computed'}, 'nonnumeric_metric'),
    ({'periodicity_dominance': True}, 'nonnumeric_metric'),
    ({'periodicity_dominance': []}, 'nonnumeric_metric'),
    ({'periodicity_dominance': {'28': .9}}, 'missing_global_value'),
    ({'periodicity_dominance': {'global': None}}, 'unavailable_metric'),
    ({'periodicity_dominance': float('nan')}, 'nonfinite_metric'),
    ({'periodicity_dominance': float('inf')}, 'nonfinite_metric'),
    ({'periodicity_dominance': float('-inf')}, 'nonfinite_metric'),
    ({'periodicity_dominance': 10 ** 1000}, 'nonfinite_metric'),
])
def test_missing_or_invalid_measurement_fails(metrics,reason):
    result=evaluate_qc_status({'metrics':metrics},'sample',P)
    assert result['overall_status']=='FAIL'
    assert result['summary']['total_checks']==1
    assert result['summary']['unavailable']==1
    assert result['checks'][0]['reason']==reason
    assert result['checks'][0]['value'] is None
    json.dumps(result,allow_nan=False)


def test_partial_evidence_cannot_silently_pass():
    policy=dict(P,uniformity_entropy={'pass':.7,'warn':.5})
    result=evaluate_qc_status({'metrics':{'periodicity_dominance':.9}},'sample',policy)
    assert [c['status'] for c in result['checks']]==['PASS','FAIL']
    assert result['overall_status']=='FAIL'
    assert result['summary']['total_checks']==len(policy)

@pytest.mark.parametrize('direction,passval,warnval,values,statuses', [
    ('higher',.7,.5,[.8,.7,.6,.5,.4],['PASS','PASS','WARNING','WARNING','FAIL']),
    ('lower',.1,.2,[.05,.1,.15,.2,.3],['PASS','PASS','WARNING','WARNING','FAIL']),
])
def test_finite_contract_unchanged(direction,passval,warnval,values,statuses):
    for value,status in zip(values,statuses):
        result=evaluate_qc_status({'metrics':{'m':value}},'x',{'m':{'pass':passval,'warn':warnval,'direction':direction}})
        assert result['overall_status']==status

@pytest.mark.parametrize('value',[np.float64(.9),np.int64(1),{'global':.9}])
def test_numeric_and_global_measurements(value):
    result=evaluate_qc_status({'metrics':{'periodicity_dominance':value}},'x',P)
    assert result['overall_status']=='PASS'
    json.dumps(result,allow_nan=False)

@pytest.mark.parametrize('policy', [
    {}, [], {'m':{}}, {'m':{'pass':.5}}, {'m':{'pass':True,'warn':.1}},
    {'m':{'pass':'0.5','warn':.1}}, {'m':{'pass':math.nan,'warn':.1}},
    {'m':{'pass':math.inf,'warn':.1}}, {'m':{'pass':.5,'warn':math.inf}},
    {'m':{'pass':.5,'warn':.1,'direction':'Lower'}},
    {'m':{'pass':.5,'warn':.1,'direction':'lower'}},
    {'m':{'pass':.1,'warn':.5,'direction':'higher'}},
    {'m':{'pass':.5,'warn':.1,'directon':'lower'}}, {3:{'pass':.5,'warn':.1}},
    {'m':None}, {'m':{'pass':10**1000,'warn':.1}},
])
def test_invalid_policy_is_rejected(policy):
    with pytest.raises(ValueError):evaluate_qc_status({'metrics':{}},'x',policy)

@pytest.mark.parametrize('results',[{},[],{'metrics':None},{'metrics':[]}])
def test_invalid_results_mapping_rejected(results):
    with pytest.raises(ValueError):evaluate_qc_status(results,'x',P)


def test_no_explicit_policy_keeps_scored_resolver(monkeypatch):
    sentinel={'unchanged':'resolver'}
    monkeypatch.setattr(output,'_evaluate_qc_status_scored',lambda r,n:sentinel)
    assert evaluate_qc_status({'metrics':{}},'x',None) is sentinel

@pytest.mark.parametrize('prime',['5prime','3prime'])
def test_reported_infinite_raw_kl_is_not_skipped(prime):
    key=f'terminal_bias_kl_{prime}_raw'
    metrics={key:None,key+'_status':'zero_background_support'}
    result=evaluate_qc_status({'metrics':metrics},'x',{key:{'pass':.1,'warn':.2}})
    assert result['overall_status']=='FAIL'
    assert result['checks'][0]['value_status']=='zero_background_support'
    assert result['summary']['unavailable']==0
    json.dumps(result,allow_nan=False)


def test_unknown_null_is_not_equated_with_infinite_divergence():
    key='terminal_bias_kl_5prime_raw'
    result=evaluate_qc_status({'metrics':{key:None}},'x',{key:{'pass':.1,'warn':.2}})
    assert result['checks'][0]['reason']=='unavailable_metric'


def test_inputs_not_mutated():
    data={'metrics':{'periodicity_dominance':{'global':.9}}};saved=copy.deepcopy(data)
    policy=copy.deepcopy(P)
    evaluate_qc_status(data,'x',policy)
    assert data==saved and policy==P


def test_seeded_finite_oracle():
    rng=np.random.default_rng(6062026)
    for direction in ('lower','higher'):
        for _ in range(250):
            a,b=sorted(rng.uniform(-2,2,size=2));v=float(rng.uniform(-3,3))
            passed,warned=(a,b) if direction=='lower' else (b,a)
            expected=('PASS' if v<=passed else 'WARNING' if v<=warned else 'FAIL') if direction=='lower' else ('PASS' if v>=passed else 'WARNING' if v>=warned else 'FAIL')
            result=evaluate_qc_status({'metrics':{'m':v}},'x',{'m':{'pass':passed,'warn':warned,'direction':direction}})
            assert result['overall_status']==expected
