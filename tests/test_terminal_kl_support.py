"""Zero-background regressions for the raw function and its goodness score."""
import math
import numpy as np
import pytest
from scipy.special import rel_entr
from RiboMetric.metrics import (
    terminal_nucleotide_bias_KL_divergence as raw,
    terminal_nucleotide_bias_KL_metric as score,
)


@pytest.mark.parametrize("prime", ["five_prime", "three_prime"])
@pytest.mark.parametrize("nested", [False, True])
def test_disjoint_support_is_not_perfect(prime, nested):
    observed = {prime: {"AA": 1.0, "CC": 0.0}}
    expected = {"AA": 0.0, "CC": 1.0}
    if nested:
        expected = {prime: expected}
    assert raw(observed, expected, prime) == math.inf
    assert score(observed, expected, prime) == 0.0


@pytest.mark.parametrize("nested", [False, True])
def test_zero_over_zero_contributes_zero(nested):
    observed = {"five_prime": {"AA": 1.0, "CC": 0.0}}
    expected = {"AA": 1.0, "CC": 0.0}
    if nested:
        expected = {"five_prime": expected}
    assert raw(observed, expected) == 0.0
    assert score(observed, expected) == 1.0


def test_partial_support_mismatch_dominates_other_terms():
    assert raw({"five_prime": {"AA": .1, "CC": .9}}, {"AA": 0, "CC": 1}) == math.inf


@pytest.mark.parametrize("p,q", [
    ([.2, .8], [.5, .5]), ([.1, .2, .7], [.4, .3, .3]),
    ([1.0, 0.0], [.1, .9]), ([.2, .8], [.2, .8]),
])
def test_finite_distributions_match_independent_reference(p, q):
    keys = [str(i) for i in range(len(p))]
    actual = raw({"five_prime": dict(zip(keys, p))}, dict(zip(keys, q)))
    expected = float(np.sum(rel_entr(p, q)) / np.log(2))
    assert actual == pytest.approx(expected, abs=2e-14)


def test_extreme_positive_background_does_not_overflow_ratio():
    expected = -math.log2(1e-320)
    assert raw({"five_prime": {"AA": 1}}, {"AA": 1e-320}) == pytest.approx(expected)


@pytest.mark.parametrize("value", [-.1, 1.1, math.nan, math.inf, -math.inf])
@pytest.mark.parametrize("which", ["observed", "expected"])
def test_invalid_values_raise(value, which):
    observed = {"five_prime": {"AA": 1}}
    expected = {"AA": 1}
    if which == "observed":
        observed["five_prime"]["AA"] = value
    else:
        expected["AA"] = value
    with pytest.raises(ValueError):
        raw(observed, expected)


def test_later_invalid_value_is_not_hidden_by_earlier_infinity():
    with pytest.raises(ValueError):
        raw({"five_prime": {"AA": 1, "CC": math.nan}}, {"AA": 0, "CC": 1})


def test_missing_expected_key_remains_an_error():
    with pytest.raises(KeyError):
        raw({"five_prime": {"AA": 1}}, {"CC": 1})


def test_inputs_not_mutated():
    observed = {"five_prime": {"AA": 0.4, "CC": 0.6}}
    expected = {"AA": .5, "CC": .5}
    raw(observed, expected)
    assert observed == {"five_prime": {"AA": .4, "CC": .6}}
    assert expected == {"AA": .5, "CC": .5}
