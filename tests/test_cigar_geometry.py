import pandas as pd

from RiboMetric.cigar_geometry import (
    aligned_query_length,
    five_prime_soft_clip,
    is_profile_eligible,
    project_five_prime_offset,
)
from RiboMetric.modules import _project_sites


def test_aligned_query_length_excludes_soft_clips():
    assert aligned_query_length("3S28M4S") == 28
    assert is_profile_eligible(20, 0)
    assert is_profile_eligible(24, 0)
    assert not is_profile_eligible(20, 1)
    assert not is_profile_eligible(19, 0)


def test_forward_and_reverse_five_prime_origins():
    assert project_five_prime_offset("3S28M", 100, False, 0) == 100
    assert project_five_prime_offset("3S28M", 100, False, 27) == 127
    assert project_five_prime_offset("28M4S", 100, True, 0) == 127
    assert five_prime_soft_clip("3S28M", False) == 3
    assert five_prime_soft_clip("28M4S", True) == 4


def test_query_offset_projection_through_insertions_deletions_and_splices():
    assert project_five_prime_offset("10M2I10M", 100, False, 10) is None
    assert project_five_prime_offset("10M2I10M", 100, False, 12) == 110
    assert project_five_prime_offset("10M2D10M", 100, False, 10) == 112
    assert project_five_prime_offset("10M100N10M", 100, False, 10) == 210


def test_ribometric_site_projection_matches_the_shared_geometry():
    reads = pd.DataFrame(
        {
            "cigar": ["28M", "28M", "10M2I10M", "10M2D10M", "10M100N10M"],
            "reference_leftmost": [100, 100, 100, 100, 100],
            "reference_start": [100, 127, 100, 100, 100],
            "is_reverse": [False, True, False, False, False],
            "offset": [0, 0, 10, 10, 10],
        }
    )
    sites = _project_sites(reads)
    assert sites.iloc[[0, 1, 3, 4]].tolist() == [100, 127, 112, 210]
    assert pd.isna(sites.iloc[2])
