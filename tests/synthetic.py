"""Deterministic synthetic transcriptome and BAM with known geometry.

Used by the contract tests (docs/METRIC_CONTRACT.md). Every transcript has the
same layout: 900 nt, CDS ``[60, 660)`` in 0-based transcript coordinates. Each
carries a read population chosen so that one processing decision has one
observable consequence:

========  ==========================================================================
``TXA``   28 nt footprints, P-site 12 nt from the 5' end, every read in frame 0,
          plus a start-codon pile-up that the offset callers lock onto
``TXB``   29 nt footprints, P-site 13 nt from the 5' end, in frame 0
``TXC``   30 nt reads at every position: uniform, RNA-like coverage
``TXD``   28 nt reads only in the 5' leader and 3' trailer: no CDS reads
``TXE``   28 nt multimappers (NH 2, MAPQ 3) placed in frame 2: must not reach
          any frame-sensitive metric
``TXF``   28 nt unique reads in frame 0, each with a secondary copy (flag 256)
          shifted into frame 1: the secondaries must be dropped at parse
========  ==========================================================================

Read names carry a collapse suffix ``_xN`` so every record has a weight.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import pysam

TX_LEN = 900
CDS_START = 60
CDS_END = 660
TRUE_P_OFFSET = {28: 12, 29: 13}
A_SITE_SHIFT = 3
TRANSCRIPTS = ("TXA", "TXB", "TXC", "TXD", "TXE", "TXF")
N_CODONS = 190


@dataclass(frozen=True)
class Record:
    tid: int
    start: int  # 0-based leftmost aligned position
    length: int
    name: str
    flag: int = 0
    mapq: int = 255
    nh: int = 1

    @property
    def weight(self) -> int:
        return int(self.name.rsplit("_x", 1)[1])


def _codon_reads(
    tid: int,
    length: int,
    prefix: str,
    first_weight: int,
    weight: int,
    frame_shift: int = 0,
    flag: int = 0,
    mapq: int = 255,
    nh: int = 1,
) -> List[Record]:
    """One read per codon k whose P-site sits on codon k (shifted by frame_shift)."""
    p_offset = TRUE_P_OFFSET[length]
    reads = []
    for k in range(N_CODONS):
        start = CDS_START + 3 * k - p_offset + frame_shift
        w = first_weight if k == 0 else weight
        reads.append(Record(tid, start, length, f"{prefix}{k}_x{w}", flag, mapq, nh))
    return reads


def records() -> List[Record]:
    recs: List[Record] = []
    recs += _codon_reads(0, 28, "a", first_weight=200, weight=5)
    recs += _codon_reads(1, 29, "b", first_weight=200, weight=5)
    recs += [Record(2, p, 30, f"c{p}_x1") for p in range(0, TX_LEN - 30 + 1)]
    recs += [Record(3, p, 28, f"d{p}_x3") for p in list(range(0, 20)) + list(range(700, 860))]
    recs += _codon_reads(4, 28, "e", first_weight=5, weight=5, frame_shift=2, mapq=3, nh=2)
    primaries_f = _codon_reads(5, 28, "f", first_weight=5, weight=5)
    recs += primaries_f
    recs += [Record(5, r.start + 1, 28, r.name, flag=256) for r in primaries_f]
    return recs


def transcript_sequences(seed: int = 7) -> Dict[str, str]:
    rng = random.Random(seed)
    seqs = {}
    for tx in TRANSCRIPTS:
        s = "".join(rng.choice("ACGT") for _ in range(TX_LEN))
        seqs[tx] = s[:CDS_START] + "ATG" + s[CDS_START + 3 :]
    return seqs


def build_fixture(directory: Path, seed: int = 7) -> Tuple[Path, Path]:
    """Write ``synth.bam`` (+ index) and ``synth.tsv`` into directory."""
    directory.mkdir(parents=True, exist_ok=True)
    seqs = transcript_sequences(seed)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": tx, "LN": TX_LEN} for tx in TRANSCRIPTS],
    }
    bam = directory / "synth.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for r in sorted(records(), key=lambda r: (r.tid, r.start, r.flag, r.name)):
            seq = seqs[TRANSCRIPTS[r.tid]][r.start : r.start + r.length]
            a = pysam.AlignedSegment()
            a.query_name = r.name
            a.query_sequence = seq
            a.flag = r.flag
            a.reference_id = r.tid
            a.reference_start = r.start
            a.mapping_quality = r.mapq
            a.cigarstring = f"{r.length}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * r.length)
            a.set_tag("NH", r.nh)
            out.write(a)
    pysam.index(str(bam))
    annotation = directory / "synth.tsv"
    with open(annotation, "w") as fh:
        fh.write("transcript_id\tcds_start\tcds_end\ttranscript_length\n")
        for tx in TRANSCRIPTS:
            fh.write(f"{tx}\t{CDS_START}\t{CDS_END}\t{TX_LEN}\n")
    return bam, annotation


def primary_records() -> List[Record]:
    return [r for r in records() if not r.flag & 256]


def read_length_distribution() -> Dict[int, int]:
    dist: Dict[int, int] = {}
    for r in primary_records():
        dist[r.length] = dist.get(r.length, 0) + r.weight
    return dist
