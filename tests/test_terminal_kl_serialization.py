"""Execute the actual QC assignments and verify strict JSON round trips.

This isolates the serialization boundary without constructing a BAM or importing
unrelated optional QC dependencies. It is not a BAM/FASTA end-to-end test.
"""
import ast
import json
import math
from pathlib import Path

import numpy as np
import pytest

QC = Path(__file__).resolve().parents[1] / "RiboMetric" / "qc.py"


@pytest.mark.parametrize("raw_value,status,json_value", [
    (math.inf, "zero_background_support", None),
    (0.0, "finite", 0.0),
    (2.5, "finite", 2.5),
])
def test_raw_kl_status_and_json(raw_value, status, json_value):
    source = ast.parse(QC.read_text())
    statements = []
    keys = {"terminal_bias_kl_5prime_raw", "terminal_bias_kl_3prime_raw",
            "terminal_bias_kl_5prime_raw_status", "terminal_bias_kl_3prime_raw_status"}
    for node in ast.walk(source):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Subscript) and isinstance(target.slice, ast.Constant) and target.slice.value in keys:
                statements.append(node)
    assert len(statements) == 4, "Both raw values and both explicit status assignments must exist"
    block = ast.Module(body=sorted(statements, key=lambda n: n.lineno), type_ignores=[])
    env = {"math": math, "np": np, "_lbd5_raw": raw_value, "_lbd3_raw": raw_value,
           "results_dict": {"metrics": {}}}
    exec(compile(block, str(QC), "exec"), env)
    result = json.loads(json.dumps(env["results_dict"], allow_nan=False))
    for prime in ("5prime", "3prime"):
        assert result["metrics"][f"terminal_bias_kl_{prime}_raw"] == json_value
        assert result["metrics"][f"terminal_bias_kl_{prime}_raw_status"] == status
