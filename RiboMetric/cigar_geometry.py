"""CIGAR rules shared by offset eligibility and genomic site projection."""

from __future__ import annotations

import re
from typing import Iterable, Optional, Tuple

_OPS = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar(cigar: str) -> Tuple[Tuple[int, str], ...]:
    pairs = tuple((int(n), op) for n, op in _OPS.findall(cigar or ""))
    if not pairs or "".join(f"{n}{op}" for n, op in pairs) != (cigar or ""):
        raise ValueError(f"Malformed CIGAR: {cigar!r}")
    return pairs


def aligned_query_length(cigar: str) -> int:
    return sum(n for n, op in parse_cigar(cigar) if op in "M=XI")


def five_prime_soft_clip(cigar: str, is_reverse: bool = False) -> int:
    pairs = parse_cigar(cigar)
    oriented = tuple(reversed(pairs)) if is_reverse else pairs
    for n, op in oriented:
        if op in "HP":
            continue
        return n if op == "S" else 0
    return 0


def is_profile_eligible(aligned_length: int, five_prime_clip: int) -> bool:
    return five_prime_clip == 0 and 20 <= aligned_length <= 40


def project_five_prime_offset(
    cigar: str, reference_start: int, is_reverse: bool, offset: int
) -> Optional[int]:
    """Map a 0-based aligned-query offset to a 0-based reference coordinate.

    Insertions consume query offset but have no reference base. Deletions and
    introns advance the reference without consuming query. Clips are outside
    the aligned footprint and consume neither offset nor reference here.
    """
    if offset < 0:
        return None
    pairs = parse_cigar(cigar)
    ref = reference_start
    coordinates = []
    for n, op in pairs:
        if op in "M=X":
            coordinates.extend(range(ref, ref + n))
            ref += n
        elif op == "I":
            coordinates.extend([None] * n)
        elif op in "DN":
            ref += n
    if is_reverse:
        coordinates.reverse()
    return coordinates[offset] if offset < len(coordinates) else None
