"""Exercise real JSON/CSV/YAML I/O and the evaluate command handler."""
from argparse import Namespace
import json
from pathlib import Path
import subprocess
import sys
import pytest
import yaml
from RiboMetric.evaluate import evaluate

@pytest.mark.parametrize('fmt,content,expected',[
    ('json',{'metrics':{}},2),
    ('json',{'metrics':{'periodicity_dominance':None}},2),
    ('json',{'metrics':{'periodicity_dominance':'missing'}},2),
    ('json',{'metrics':{'periodicity_dominance':.9}},0),
    ('json',{'metrics':{'periodicity_dominance':.6}},1),
    ('json',{'metrics':{'periodicity_dominance':.1}},2),
    ('csv','metric,value\nperiodicity_dominance,not-computed\n',2),
    ('csv','metric,value\nperiodicity_dominance,nan\n',2),
    ('csv','metric,value\nperiodicity_dominance,inf\n',2),
    ('csv','metric,value\nperiodicity_dominance,0.9\n',0),
])
def test_handler_file_roundtrip(tmp_path,fmt,content,expected,capsys):
    source=tmp_path/f'input.{fmt}';thresholds=tmp_path/'policy.yml';report=tmp_path/'qc.json'
    source.write_text(json.dumps({'results':content}) if fmt=='json' else content)
    thresholds.write_text(yaml.safe_dump({'periodicity_dominance':{'pass':.7,'warn':.5}}))
    code=evaluate(Namespace(input=str(source),expected=str(thresholds),output=str(report),name='sample'))
    assert code==expected
    status=json.loads(report.read_text(),parse_constant=lambda v:pytest.fail('Nonstandard JSON constant '+v))
    assert status['overall_status']=={0:'PASS',1:'WARNING',2:'FAIL'}[expected]
    assert status['summary']['total_checks']==1

@pytest.mark.parametrize('policy',['{}','thresholds: {}','[','m:\n  pass: 1\n','m:\n  pass: .nan\n  warn: 0\n'])
def test_bad_policy_returns_failure(tmp_path,policy):
    data=tmp_path/'data.json';data.write_text('{"metrics":{}}')
    cfg=tmp_path/'policy.yml';cfg.write_text(policy)
    assert evaluate(Namespace(input=str(data),expected=str(cfg),output=None,name=None))==2


def test_missing_results_returns_failure(tmp_path):
    assert evaluate(Namespace(input=str(tmp_path/'absent.json'),expected=None))==2


def test_empty_metrics_with_cli_defaults_fails(tmp_path):
    data=tmp_path/'data.json';data.write_text('{"metrics":{}}')
    assert evaluate(Namespace(input=str(data),expected=None,output=None,name=None))==2


def test_known_infinite_raw_value_survives_handler(tmp_path,capsys):
    k='terminal_bias_kl_5prime_raw';data=tmp_path/'data.json';cfg=tmp_path/'policy.yml';out=tmp_path/'report.json'
    data.write_text(json.dumps({'metrics':{k:None,k+'_status':'zero_background_support'}}))
    cfg.write_text(yaml.safe_dump({k:{'pass':.1,'warn':.2}}))
    assert evaluate(Namespace(input=str(data),expected=str(cfg),output=str(out),name='x'))==2
    assert '+infinity (zero background support)' in capsys.readouterr().out
    status=json.loads(out.read_text())
    assert status['checks'][0]['value'] is None


def test_subprocess_handler_exit_is_nonzero(tmp_path):
    data=tmp_path/'data.json';data.write_text('{"metrics":{}}')
    cfg=tmp_path/'policy.yml';cfg.write_text('periodicity_dominance:\n  pass: 0.7\n  warn: 0.5\n')
    command=('import sys; from argparse import Namespace; '
             'from RiboMetric.evaluate import evaluate; '
             'sys.exit(evaluate(Namespace(input=sys.argv[1],expected=sys.argv[2],output=None,name="fixture")))')
    run=subprocess.run([sys.executable,'-c',command,str(data),str(cfg)],text=True,capture_output=True)
    assert run.returncode==2,run.stdout+run.stderr
    assert 'missing_metric' in run.stdout
