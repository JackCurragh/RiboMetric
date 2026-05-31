"""
RUST (Ribo-seq Unit Step Transformation) metric for RiboMetric.

Implements the codon-level RUST metagene from:
    O'Connor et al., Nat Commun 2016. "Comparative survey of the expression of
    translation machinery components across diverse bacterial species."

Algorithm (O'Connor et al.):
    For each transcript with a valid CDS annotation:
    1.  Extract the elongation region (CDS[120:-60], skipping first 40 codons
        and last 20 to avoid initiation/termination artefacts).
    2.  Build a per-nucleotide A-site density profile.
    3.  Compute average_gene_density (mean density in the elongation region).
    4.  For each codon in the elongation region slide a 60-codon window
        (centred such that the A-site codon is at window position 40, i.e.
        40 codons upstream and 19 codons downstream are captured):
          - Increment the "total" counter for every codon in the window.
          - If the current codon is "enriched" (3-nt density > average),
            also increment the "enriched" counter.
    5.  RUST ratio per (codon, window_position) = enriched / total.
    6.  KL divergence per window position using observed RUST ratios as
        the observed distribution and per-codon expected enrichment rates
        as the expected distribution.

The scalar summary metric is ``rust_mean_kl_divergence`` — the mean KL
divergence across the 60 window positions (NA positions excluded).
"""

from __future__ import annotations

import math
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# Stop codons are excluded from the metagene output.
STOP_CODONS = {"TAA", "TAG", "TGA"}

# All valid DNA codons (order preserved for reproducibility).
ALL_CODONS: List[str] = [
    f"{a}{b}{c}"
    for a in "ACGT"
    for b in "ACGT"
    for c in "ACGT"
]
SENSE_CODONS: List[str] = [c for c in ALL_CODONS if c not in STOP_CODONS]

WINDOW_SIZE = 60       # window length in codons
ASITE_IN_WINDOW = 40   # 0-based position of the A-site codon within the window
ELONGATION_5_NT = 120  # nt to skip at CDS 5' end  (40 codons)
ELONGATION_3_NT = 60   # nt to skip at CDS 3' end  (20 codons)
MIN_PROFILE_LEN = 50   # minimum elongation-region nt to include a transcript


def _lookup_seq(
    fasta_dict: dict,
    name: str,
    base_index: Optional[Dict[str, str]] = None,
) -> Optional[str]:
    """Return the sequence string for *name*, trying several key variants.

    Lookup order:
    1. Exact match                        (e.g. ``ENST00000233.10``)
    2. Pipe-stripped exact match          (e.g. ``ENST00000233.10|iso1`` -> ``ENST00000233.10``)
    3. Version-stripped match via base_index (``ENST00000233.10`` -> ``ENST00000233``)
       - needed when annotation and FASTA come from different releases.
    """
    def _seq(rec):
        return str(rec.seq) if hasattr(rec, "seq") else str(rec)

    # 1. Exact
    if name in fasta_dict:
        return _seq(fasta_dict[name])
    # 2. Pipe-stripped
    base_pipe = name.split("|")[0]
    if base_pipe in fasta_dict:
        return _seq(fasta_dict[base_pipe])
    # 3. Version-stripped (cross-release matching)
    if base_index is not None:
        base_ver = base_pipe.split(".")[0]
        if base_ver in base_index:
            return base_index[base_ver]
    return None


def compute_rust_metrics(
    annotated_read_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    fasta_dict: dict,
    window: int = WINDOW_SIZE,
    asite_pos: int = ASITE_IN_WINDOW,
    elongation_5_nt: int = ELONGATION_5_NT,
    elongation_3_nt: int = ELONGATION_3_NT,
    min_profile_len: int = MIN_PROFILE_LEN,
) -> Dict:
    """
    Compute the RUST codon-level metagene from an annotated read DataFrame.

    Parameters
    ----------
    annotated_read_df : pd.DataFrame
        Output of ``chunked_annotate_reads`` + ``a_site_calculation_variable_offset``.
        Must have columns: reference_name, a_site, cds_start, cds_end, count.
    annotation_df : pd.DataFrame
        Annotation table with columns: transcript_id, cds_start, cds_end.
    fasta_dict : dict
        ``{record_id: SeqRecord | str}`` from ``parse_fasta``.
    window : int
        Number of codons in the sliding window (default 60).
    asite_pos : int
        0-based position of the A-site codon within the window (default 40).
    elongation_5_nt : int
        Nucleotides to trim from the CDS 5' end (default 120 = 40 codons).
    elongation_3_nt : int
        Nucleotides to trim from the CDS 3' end (default 60 = 20 codons).
    min_profile_len : int
        Minimum elongation-region length (nt) to include a transcript (default 50).

    Returns
    -------
    dict with keys:
        ``metagene``         : {codon: [rust_ratio_0 … rust_ratio_{window-1}]}
        ``kl_divergence``    : list of float or None, length == window
        ``mean_kl_divergence``: mean KL over non-NA positions (scalar)
        ``transcripts_used`` : int
    """
    # ------------------------------------------------------------------ #
    # 0.  Build version-stripped FASTA index for cross-release matching  #
    # ------------------------------------------------------------------ #
    # Annotation and FASTA may come from different Ensembl/GENCODE releases
    # (e.g. annotation from Ensembl 114, FASTA from GENCODE v49).  Build a
    # {base_id: seq_str} dict keyed on ENST without the ".N" version so
    # _lookup_seq can fall back to it.
    base_fasta_index: Dict[str, str] = {}
    for fid, rec in fasta_dict.items():
        base_id = fid.split("|")[0].split(".")[0]
        if base_id not in base_fasta_index:
            seq_str = str(rec.seq) if hasattr(rec, "seq") else str(rec)
            base_fasta_index[base_id] = seq_str

    # ------------------------------------------------------------------ #
    # 1.  Build per-transcript CDS lookup from annotation_df              #
    # ------------------------------------------------------------------ #
    ann_index: Dict[str, Tuple[int, int]] = {}
    for _, row in annotation_df.iterrows():
        tid = str(row["transcript_id"])
        ann_index[tid] = (int(row["cds_start"]), int(row["cds_end"]))

    # ------------------------------------------------------------------ #
    # 2.  Pre-compute per-transcript A-site density (pandas groupby)      #
    # ------------------------------------------------------------------ #
    # We need CDS-only reads with integer A-site positions.
    transcript_col = (
        "reference_name" if "reference_name" in annotated_read_df.columns
        else "transcript_id" if "transcript_id" in annotated_read_df.columns
        else None
    )
    needed = {"a_site", "cds_start", "cds_end", "count"}
    missing = needed - set(annotated_read_df.columns)
    if transcript_col is None:
        missing.add("reference_name/transcript_id")
    if missing:
        log.warning("RUST: annotated_read_df missing columns %s; skipping.", missing)
        return _empty_result(window)

    cds_df = annotated_read_df.dropna(subset=["a_site", "cds_start", "cds_end"])
    cds_df = cds_df[
        (cds_df["a_site"] >= cds_df["cds_start"])
        & (cds_df["a_site"] < cds_df["cds_end"])
    ].copy()

    if cds_df.empty:
        log.warning("RUST: no CDS reads found; skipping.")
        return _empty_result(window)

    cds_df["a_site"] = cds_df["a_site"].astype(int)
    cds_df["cds_start"] = cds_df["cds_start"].astype(int)
    cds_df["cds_end"] = cds_df["cds_end"].astype(int)
    cds_df["count"] = cds_df["count"].astype(float)

    # Group: (transcript, a_site) -> total weighted count
    density_series = (
        cds_df.groupby([transcript_col, "a_site"], observed=True)["count"]
        .sum()
    )

    # Per-transcript CDS boundaries (first occurrence suffices; they are constant)
    tx_info = (
        cds_df.groupby(transcript_col, observed=True)[["cds_start", "cds_end"]]
        .first()
    )

    # ------------------------------------------------------------------ #
    # 3.  Initialise enrichment accumulators                              #
    # ------------------------------------------------------------------ #
    # codon_enrichment[codon][window_pos] = [total, enriched]
    codon_enrichment: Dict[str, List[List[float]]] = {
        codon: [[0.0, 0.0] for _ in range(window)] for codon in ALL_CODONS
    }
    # per-codon expected_codon_density values (list per codon for averaging)
    codon_expected: Dict[str, List[float]] = {codon: [] for codon in ALL_CODONS}

    transcripts_used = 0

    # ------------------------------------------------------------------ #
    # 4.  Per-transcript RUST computation                                 #
    # ------------------------------------------------------------------ #
    for transcript in tx_info.index:
        cds_start = int(tx_info.loc[transcript, "cds_start"])
        cds_end = int(tx_info.loc[transcript, "cds_end"])

        # Sequence lookup (tries exact -> pipe-stripped -> version-stripped)
        seq = _lookup_seq(fasta_dict, str(transcript), base_index=base_fasta_index)
        if seq is None:
            continue

        # CDS sub-sequence (0-based, half-open on transcript coords)
        cds_seq = seq[cds_start:cds_end]
        if len(cds_seq) < elongation_5_nt + elongation_3_nt + min_profile_len:
            continue
        if len(cds_seq) % 3 != 0:
            continue

        # Elongation region
        elong_len = len(cds_seq) - elongation_5_nt - elongation_3_nt
        if elong_len < min_profile_len:
            continue
        if set(cds_seq) - set("ACGTacgt"):
            # Skip transcripts with ambiguous bases for simplicity
            cds_seq = cds_seq.upper()
            if set(cds_seq) - set("ACGT"):
                continue
        cds_seq = cds_seq.upper()

        # Build density profile (indexed within elongation region)
        profile = np.zeros(elong_len, dtype=float)
        tx_key = transcript
        if (tx_key, cds_start) in density_series.index or tx_key in density_series.index.get_level_values(0):
            try:
                tx_density = density_series.loc[tx_key]
            except KeyError:
                continue
            for a_abs, cnt in tx_density.items():
                # Convert absolute A-site to elongation-region index
                idx = int(a_abs) - cds_start - elongation_5_nt
                if 0 <= idx < elong_len:
                    profile[idx] += cnt
        else:
            continue

        if profile.sum() == 0:
            continue

        avg_density = float(profile.sum()) / elong_len

        # Expected codon density for this transcript (fraction of codons enriched)
        num_enriched = sum(
            1
            for n in range(0, elong_len - 2, 3)
            if (profile[n] + profile[n + 1] + profile[n + 2]) / 3.0 > avg_density
        )
        total_codons = elong_len // 3
        if total_codons == 0:
            continue
        expected_codon_density = float(num_enriched) / total_codons

        # ---------------------------------------------------------------- #
        # Sliding window over codons in the elongation region              #
        # Window starts at codon_start *in cds_seq* (not elongation_part). #
        # When we are at elongation position p=0, codon_start=0 in cds_seq #
        # => window = cds_seq[0:window*3] (first 60 codons of CDS).        #
        # A-site codon = window[asite_pos*3:(asite_pos+1)*3].               #
        # ---------------------------------------------------------------- #
        codon_start_in_cds = 0  # tracks where window begins in cds_seq
        for elong_pos in range(0, elong_len - 2, 3):
            codon_window = cds_seq[codon_start_in_cds: codon_start_in_cds + window * 3]
            if len(codon_window) < window * 3:
                codon_start_in_cds += 3
                continue
            # Skip windows with non-ACGT bases
            if set(codon_window) - set("ACGT"):
                codon_start_in_cds += 3
                continue

            # Is the current codon (A-site) enriched?
            codon_density_3nt = (
                profile[elong_pos]
                + (profile[elong_pos + 1] if elong_pos + 1 < elong_len else 0)
                + (profile[elong_pos + 2] if elong_pos + 2 < elong_len else 0)
            ) / 3.0
            enriched = codon_density_3nt > avg_density

            # Update accumulators for every codon in the window
            for win_pos in range(window):
                codon = codon_window[win_pos * 3: (win_pos + 1) * 3]
                codon_enrichment[codon][win_pos][0] += 1.0
                if enriched:
                    codon_enrichment[codon][win_pos][1] += 1.0

            # Record expected density for the A-site codon
            asite_codon = codon_window[asite_pos * 3: (asite_pos + 1) * 3]
            if asite_codon not in STOP_CODONS:
                codon_expected[asite_codon].append(expected_codon_density)

            codon_start_in_cds += 3

        transcripts_used += 1

    if transcripts_used == 0:
        log.warning("RUST: no transcripts contributed; returning empty result.")
        return _empty_result(window)

    # ------------------------------------------------------------------ #
    # 5.  Compute RUST ratios per (codon, window_position)                #
    # ------------------------------------------------------------------ #
    metagene: Dict[str, List[float]] = {}
    for codon in SENSE_CODONS:
        rust_ratios = []
        for win_pos in range(window):
            total, enriched = codon_enrichment[codon][win_pos]
            rust_ratios.append(enriched / total if total > 0 else 0.0)
        metagene[codon] = rust_ratios

    # ------------------------------------------------------------------ #
    # 6.  KL divergence per window position                               #
    # ------------------------------------------------------------------ #
    # Expected distribution: per-codon mean expected_codon_density
    exp_values = [
        float(np.mean(codon_expected[c])) if codon_expected[c] else 0.0
        for c in SENSE_CODONS
    ]
    exp_sum = sum(exp_values)
    if exp_sum == 0:
        q_values = [1.0 / len(SENSE_CODONS)] * len(SENSE_CODONS)
    else:
        q_values = [v / exp_sum for v in exp_values]

    kl_divergence: List[Optional[float]] = []
    for win_pos in range(window):
        obs = [metagene[c][win_pos] for c in SENSE_CODONS]
        obs_sum = sum(obs)
        if obs_sum == 0 or min(obs) == 0:
            kl_divergence.append(None)
            continue
        p_values = [v / obs_sum for v in obs]
        kl = sum(
            abs(p * math.log2(p / q))
            for p, q in zip(p_values, q_values)
            if p > 0 and q > 0
        )
        kl_divergence.append(kl)

    valid_kl = [v for v in kl_divergence if v is not None]
    mean_kl = float(np.mean(valid_kl)) if valid_kl else 0.0

    return {
        "metagene": metagene,
        "kl_divergence": kl_divergence,
        "mean_kl_divergence": mean_kl,
        "transcripts_used": transcripts_used,
        "window_size": window,
        "asite_position_in_window": asite_pos,
    }


def _empty_result(window: int) -> Dict:
    return {
        "metagene": {},
        "kl_divergence": [None] * window,
        "mean_kl_divergence": 0.0,
        "transcripts_used": 0,
        "window_size": window,
        "asite_position_in_window": ASITE_IN_WINDOW,
    }
