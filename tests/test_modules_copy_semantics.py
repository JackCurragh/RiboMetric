from RiboMetric.modules import normalise_ligation_bias


def test_normalise_ligation_bias_does_not_mutate_input():
    observed = {
        "five_prime": {"AA": 0.6, "AC": 0.4},
        "three_prime": {"AA": 0.5, "AC": 0.5},
    }
    expected = {
        "5_prime_bg": {"AA": 0.5, "AC": 0.5},
        "3_prime_bg": {"AA": 0.5, "AC": 0.5},
    }

    observed_copy = {k: v.copy() for k, v in observed.items()}
    out = normalise_ligation_bias(observed, expected)

    # Input unchanged
    assert observed == observed_copy

    # Output differs where expected
    assert out["five_prime"]["AA"] == observed["five_prime"]["AA"] - 0.5
    assert out["three_prime"]["AC"] == observed["three_prime"]["AC"] - 0.5

