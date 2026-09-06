"""Modified KL -> actual QC assignments -> JSON writer -> command handler.

The assignment block is AST-selected from qc.py, avoiding BAM dependencies.
This is a cross-function regression, NOT a full BAM processing test.
"""
from argparse import Namespace
import ast
import json
import math
from pathlib import Path
import numpy as np
import pytest
import yaml
from RiboMetric.metrics import terminal_nucleotide_bias_KL_divergence
from RiboMetric.results_output import generate_json
from RiboMetric.evaluate import evaluate


def assign_raw_values(raw5,raw3):
    source=Path(__file__).resolve().parents[1]/'RiboMetric/qc.py'
    names={'_lbd5','_lbd3'}
    keys={f'terminal_bias_kl_{p}{suffix}' for p in ('5prime','3prime') for suffix in ('','_raw','_score','_raw_status')}
    keys.update(f'terminal_nucleotide_bias_distribution_{p}_metric' for p in ('5_prime','3_prime'))
    selected=[]
    for n in ast.walk(ast.parse(source.read_text())):
        if not isinstance(n,ast.Assign):continue
        t=n.targets[0]
        if (isinstance(t,ast.Name) and t.id in names) or (isinstance(t,ast.Subscript) and isinstance(t.slice,ast.Constant) and t.slice.value in keys):
            selected.append(n)
    assert len(selected)==12
    env={'np':np,'results_dict':{'metrics':{}},'_lbd5_raw':raw5,'_lbd3_raw':raw3}
    exec(compile(ast.Module(body=sorted(selected,key=lambda n:n.lineno),type_ignores=[]),str(source),'exec'),env)
    return env['results_dict']

@pytest.mark.parametrize('prime',['5prime','3prime'])
@pytest.mark.parametrize('score_gate',[False,True])
def test_disjoint_support_fails_after_written_json(tmp_path,prime,score_gate):
    raw=terminal_nucleotide_bias_KL_divergence({'five_prime':{'AA':1.,'CC':0.}},{'AA':0.,'CC':1.})
    assert math.isinf(raw)
    results=assign_raw_values(raw,raw)
    data=tmp_path/'results.json';generate_json(results,{},name=str(data))
    json.loads(data.read_text(),parse_constant=lambda token:pytest.fail(token))
    k=f'terminal_bias_kl_{prime}'+('_score' if score_gate else '_raw')
    thresholds={k:({'pass':.7,'warn':.5} if score_gate else {'pass':.1,'warn':.2})}
    cfg=tmp_path/'policy.yml';cfg.write_text(yaml.safe_dump(thresholds))
    out=tmp_path/'gate.json'
    assert evaluate(Namespace(input=str(data),expected=str(cfg),output=str(out),name='x'))==2
    report=json.loads(out.read_text());assert report['overall_status']=='FAIL'
    assert report['summary']['total_checks']==1
    assert results['metrics'][f'terminal_bias_kl_{prime}_score']==0.


def test_matching_distributions_keep_pass(tmp_path):
    raw=terminal_nucleotide_bias_KL_divergence({'five_prime':{'AA':1.,'CC':0.}},{'AA':1.,'CC':0.})
    results=assign_raw_values(raw,raw)
    assert results['metrics']['terminal_bias_kl_5prime_raw']==0.
    assert results['metrics']['terminal_bias_kl_5prime_score']==1.
    json.dumps(results,allow_nan=False)
