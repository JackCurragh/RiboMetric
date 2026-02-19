"""
Pytest configuration and fixtures for RiboMetric tests
"""

import pytest
import pandas as pd
import os
from pathlib import Path


@pytest.fixture(scope="session")
def test_data_dir():
    """Return the path to the test data directory"""
    return Path(__file__).parent / "test_data"


@pytest.fixture(scope="session")
def config_path():
    """Return the path to the default config file"""
    return Path(__file__).parent.parent / "RiboMetric" / "config.yml"


@pytest.fixture
def sample_read_df(test_data_dir):
    """Load and return a sample read dataframe"""
    # Provide both expanded and weighted fixtures equivalently.
    read_df_pre = pd.read_csv(test_data_dir / "test.csv")
    # For pipeline code, we now use weights; keep original columns and do not expand.
    # Tests that expect expanded semantics either aggregate or use modules that now accept weights.
    return read_df_pre


@pytest.fixture
def sample_annotation_df(test_data_dir):
    """Load and return a sample annotation dataframe"""
    return pd.read_csv(
        test_data_dir / "test_annotation.tsv",
        sep="\t",
        dtype={
            "transcript_id": str,
            "cds_start": int,
            "cds_end": int,
            "transcript_length": int
        }
    )


@pytest.fixture
def sample_read_length_dict():
    """Return a sample read length distribution for testing"""
    return {
        20: 10,
        21: 15,
        25: 50,
        26: 80,
        27: 150,
        28: 200,
        29: 250,
        30: 200,
        31: 150,
        32: 80,
        33: 50,
        34: 30,
        35: 15,
        40: 10
    }


@pytest.fixture
def sample_read_frame_dict():
    """Return a sample read frame distribution for testing"""
    return {
        27: {0: 100, 1: 5, 2: 5},
        28: {0: 1000, 1: 50, 2: 50},
        29: {0: 800, 1: 20, 2: 10},
        30: {0: 10000, 1: 500, 2: 500},
        31: {0: 900, 1: 50, 2: 50},
        32: {0: 1000, 1: 100, 2: 100},
    }


@pytest.fixture
def sample_metagene_profile():
    """Return a sample metagene profile for testing"""
    profile = {"start": {}, "stop": {}}
    for read_len in [28, 29, 30, 31, 32]:
        profile["start"][read_len] = {
            i: (100 if i % 3 == 0 else 10) for i in range(-30, 90)
        }
        profile["stop"][read_len] = {
            i: (100 if i % 3 == 0 else 10) for i in range(-30, 90)
        }
    return profile


@pytest.fixture
def sample_sequence_background():
    """Return sample sequence background frequencies"""
    nucleotides = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
                   'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    freq = 1.0 / 16
    return {
        "5_prime_bg": {nt: freq for nt in nucleotides},
        "3_prime_bg": {nt: freq for nt in nucleotides}
    }


@pytest.fixture
def default_config():
    """Return a default configuration for testing"""
    return {
        "argument": {
            "bam": None,
            "annotation": None,
            "gff": None,
            "fasta": None,
            "subsample": None,
            "transcripts": None,
            "threads": 2,
            "offset_calculation_method": "changepoint",
        },
        "plots": {
            "read_frame_distribution": {
                "upper_limit": 40,
                "lower_limit": 20
            },
            "terminal_nucleotide_bias_distribution": {
                "nucleotide_count": 2,
                "background_freq": True
            },
            "metagene_profile": {
                "distance_target": "both",
                "distance_range": [-50, 50]
            },
            "mRNA_distribution": {
                "absolute_counts": False
            }
        },
        "qc": {
            "use_cds_subset": {
                "read_frame_distribution": True
            },
            "read_frame_distribution": {
                "3nt_count_cutoff": 0.05
            },
            "cds_coverage": {
                "in_frame_coverage": True
            }
        },
        "metrics": {
            "default": [
                "read_length_distribution_IQR",
                "read_length_distribution_coefficient_of_variation",
                "read_length_distribution_maxprop",
                "terminal_nucleotide_bias_KL_5prime",
                "terminal_nucleotide_bias_KL_3prime",
                "periodicity_dominance",
                "periodicity_information",
                "uniformity_entropy",
                "CDS_coverage",
                "region_ratios",
                "region_proportions"
            ],
            "optional": [
                "read_length_distribution_bimodality",
                "read_length_distribution_normality",
                "periodicity_autocorrelation",
                "periodicity_fourier",
                "periodicity_trips_viz",
                "uniformity_autocorrelation",
                "uniformity_theil_index",
                "uniformity_gini_index"
            ]
        },
        "enabled_metrics": []  # Will be populated by tests
    }
