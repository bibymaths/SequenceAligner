"""FASTA parsing and alignment statistics utilities.

This module provides helper functions to parse aligned FASTA files and
compute residue‑level statistics. It also defines functions to build
coordinate mappings between alignment columns and ungapped residue
indices and to detect contiguous conserved blocks.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

logger = logging.getLogger(__name__)


def parse_alignment_fasta(path: Path) -> Dict[str, str]:
    """Parse a gapped alignment FASTA file.

    Parameters
    ----------
    path : Path
        Path to a FASTA file containing two aligned protein sequences
        including gap characters ('-').

    Returns
    -------
    Dict[str, str]
        A mapping from sequence identifier to sequence string. The
        sequences are returned as uppercase strings. If duplicate IDs
        exist, the last instance wins.
    """
    records = SeqIO.parse(str(path), "fasta")
    seqs: Dict[str, str] = {}
    for rec in records:
        seqs[str(rec.id)] = str(rec.seq).upper()
    if len(seqs) < 2:
        logger.warning("Expected at least two sequences in %s, found %d", path, len(seqs))
    return seqs


def parse_ungapped_fasta(path: Path) -> Dict[str, str]:
    """Parse an ungapped FASTA file.

    Parameters
    ----------
    path : Path
        Path to a FASTA file containing an ungapped sequence.

    Returns
    -------
    Dict[str, str]
        Mapping from sequence identifier to sequence string.
    """
    return parse_alignment_fasta(path)


def compute_alignment_stats(
        seq_a: str,
        seq_b: str,
        substitution_matrix: Optional[Dict[str, Dict[str, int]]] = None,
        similarity_threshold: int = 0,
) -> Dict[str, float]:
    """Compute summary statistics for a pairwise alignment.

    Parameters
    ----------
    seq_a : str
        First aligned sequence, including gaps ('-').
    seq_b : str
        Second aligned sequence, including gaps ('-').
    substitution_matrix : dict, optional
        A nested dictionary encoding substitution scores (e.g. BLOSUM62).
        When provided, positions with a positive substitution score are
        counted toward similarity. If ``None``, similarity metrics
        are not computed and returned values will be ``NaN``.
    similarity_threshold : int, default 0
        Score threshold above which a substitution is considered
        similar. Typical values are 0 or greater.

    Returns
    -------
    dict
        Dictionary of summary statistics, including alignment length,
        ungapped lengths, counts of matches, mismatches and gaps, and
        percent identity/similarity.

    Notes
    -----
    Percent identity is defined as the number of exact matches divided
    by the alignment length. Percent similarity counts positions with
    substitution score >= ``similarity_threshold`` (if a matrix is
    provided) and is divided by the alignment length. When no
    substitution matrix is given, ``percent_similarity`` will be NaN.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("Aligned sequences must have the same length")
    aln_len = len(seq_a)
    ungapped_a = seq_a.replace("-", "")
    ungapped_b = seq_b.replace("-", "")
    ungapped_a_len = len(ungapped_a)
    ungapped_b_len = len(ungapped_b)

    matches = 0
    similar = 0
    mismatches = 0
    gaps = 0

    for aa, bb in zip(seq_a, seq_b):
        if aa == '-' or bb == '-':
            gaps += 1
            continue
        if aa == bb:
            matches += 1
            if substitution_matrix is not None:
                similar += 1  # matches are also similar
        else:
            if substitution_matrix is not None:
                # safe lookup: unknown characters return 0
                score = substitution_matrix.get(aa, {}).get(bb, 0)
                if score >= similarity_threshold:
                    similar += 1
                else:
                    mismatches += 1
            else:
                mismatches += 1

    # compute percents
    percent_identity = matches / aln_len if aln_len > 0 else float('nan')
    percent_similarity = float('nan')
    if substitution_matrix is not None:
        percent_similarity = similar / aln_len if aln_len > 0 else float('nan')

    return {
        "alignment_length": aln_len,
        "ungapped_length_a": ungapped_a_len,
        "ungapped_length_b": ungapped_b_len,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "percent_identity": percent_identity,
        "percent_similarity": percent_similarity,
    }


def build_coordinate_maps(seq_a: str, seq_b: str) -> Tuple[List[Optional[int]], List[Optional[int]]]:
    """Generate coordinate maps from alignment positions to residue indices.

    Each element of the returned lists corresponds to an alignment
    column. For sequence A (first list), the value is the 0‑based
    residue index in the ungapped sequence or ``None`` if the aligned
    character at that position is a gap. The same holds for sequence B.

    Parameters
    ----------
    seq_a : str
        Aligned sequence for protein A.
    seq_b : str
        Aligned sequence for protein B.

    Returns
    -------
    Tuple[List[Optional[int]], List[Optional[int]]]
        Two lists of equal length mapping alignment columns to residue
        indices.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("Sequences must have the same length")
    a_map: List[Optional[int]] = []
    b_map: List[Optional[int]] = []
    a_idx = 0
    b_idx = 0
    for aa, bb in zip(seq_a, seq_b):
        if aa == '-':
            a_map.append(None)
        else:
            a_map.append(a_idx)
            a_idx += 1
        if bb == '-':
            b_map.append(None)
        else:
            b_map.append(b_idx)
            b_idx += 1
    return a_map, b_map


def extract_residue_pairs(
        seq_a: str, seq_b: str, a_map: List[Optional[int]], b_map: List[Optional[int]]
) -> List[Tuple[Optional[int], str, Optional[int], str]]:
    """Construct a list of residue pair information for each alignment column.

    Parameters
    ----------
    seq_a : str
        Aligned sequence for protein A.
    seq_b : str
        Aligned sequence for protein B.
    a_map : List[Optional[int]]
        Mapping from alignment column to residue index in sequence A.
    b_map : List[Optional[int]]
        Mapping from alignment column to residue index in sequence B.

    Returns
    -------
    List[Tuple[Optional[int], str, Optional[int], str]]
        Each tuple contains (index_a, residue_a, index_b, residue_b) for
        each alignment column. Gaps are represented by ``None`` for
        index and '-' for residue.
    """
    pairs: List[Tuple[Optional[int], str, Optional[int], str]] = []
    for idx, (aa, bb) in enumerate(zip(seq_a, seq_b)):
        pairs.append((a_map[idx], aa, b_map[idx], bb))
    return pairs


def contiguous_blocks(
        seq_a: str,
        seq_b: str,
        a_map: List[Optional[int]],
        b_map: List[Optional[int]],
        substitution_matrix: Dict[str, Dict[str, int]],
        min_block_length: int = 5,
        identity_threshold: float = 0.7,
        similarity_threshold: float = 0.8,
) -> List[Dict[str, object]]:
    """Detect contiguous conserved blocks in an alignment.

    A block is a contiguous stretch of alignment columns with no gaps in
    either sequence. Each block is classified based on the fraction of
    identical residues and substitution similarity.

    Parameters
    ----------
    seq_a, seq_b : str
        Aligned sequences.
    a_map, b_map : list
        Coordinate maps from alignment columns to residue indices for
        sequence A and B.
    substitution_matrix : dict
        Substitution scoring dictionary.
    min_block_length : int, optional
        Minimum length for a block to be reported.
    identity_threshold : float, optional
        Fraction of identical residues required to classify a block as
        ``"high_identity"``.
    similarity_threshold : float, optional
        Fraction of positions with substitution score >= 0 required to
        classify a block as ``"conservative"`` if not high identity.

    Returns
    -------
    List[dict]
        A list of block dictionaries, each containing block statistics
        and coordinates. Empty list if no block meets ``min_block_length``.
    """
    assert len(seq_a) == len(seq_b) == len(a_map) == len(b_map)
    blocks: List[Dict[str, object]] = []
    current_start: Optional[int] = None
    for i, (aa, bb) in enumerate(zip(seq_a, seq_b)):
        if aa != '-' and bb != '-':
            if current_start is None:
                current_start = i
        else:
            if current_start is not None:
                end = i  # non-inclusive
                if end - current_start >= min_block_length:
                    blocks.append(_summarize_block(
                        seq_a, seq_b, a_map, b_map, substitution_matrix,
                        current_start, end, identity_threshold, similarity_threshold,
                    ))
                current_start = None
    # handle tail block
    if current_start is not None:
        end = len(seq_a)
        if end - current_start >= min_block_length:
            blocks.append(_summarize_block(
                seq_a, seq_b, a_map, b_map, substitution_matrix,
                current_start, end, identity_threshold, similarity_threshold,
            ))
    return blocks


def _summarize_block(
        seq_a: str,
        seq_b: str,
        a_map: List[Optional[int]],
        b_map: List[Optional[int]],
        substitution_matrix: Dict[str, Dict[str, int]],
        start: int,
        end: int,
        identity_threshold: float,
        similarity_threshold: float,
) -> Dict[str, object]:
    """Compute statistics and classification for a contiguous block.

    Parameters
    ----------
    seq_a, seq_b, a_map, b_map : see ``contiguous_blocks``.
    start, end : int
        Alignment indices delimiting the block (end is exclusive).
    identity_threshold, similarity_threshold : float
        Classification thresholds.

    Returns
    -------
    dict
        Summary statistics for the block.
    """
    block_length = end - start
    identities = 0
    similarities = 0
    mismatches = 0
    # gather sequences for classification
    for i in range(start, end):
        aa = seq_a[i]
        bb = seq_b[i]
        if aa == bb:
            identities += 1
            similarities += 1  # identical is also similar
        else:
            score = substitution_matrix.get(aa, {}).get(bb, 0)
            if score >= 0:
                similarities += 1
            else:
                mismatches += 1
    frac_identity = identities / block_length if block_length > 0 else 0.0
    frac_similarity = similarities / block_length if block_length > 0 else 0.0
    # determine classification
    if frac_identity >= identity_threshold:
        classification = "high_identity"
    elif frac_similarity >= similarity_threshold:
        classification = "conservative"
    else:
        classification = "mismatch_rich"
    # derive coordinate ranges in residue indices
    a_start = a_map[start]
    a_end = a_map[end - 1]
    if a_start is not None and a_end is not None:
        # inclusive range; convert to residue range
        a_range = (a_start, a_end)
    else:
        a_range = None
    b_start = b_map[start]
    b_end = b_map[end - 1]
    if b_start is not None and b_end is not None:
        b_range = (b_start, b_end)
    else:
        b_range = None
    return {
        "start": start,
        "end": end,
        "length": block_length,
        "identities": identities,
        "similarities": similarities,
        "mismatches": mismatches,
        "frac_identity": frac_identity,
        "frac_similarity": frac_similarity,
        "classification": classification,
        "seqA_range": a_range,
        "seqB_range": b_range,
    }
