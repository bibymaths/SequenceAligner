"""Residue‑level support profile calculations.

This module exposes functions to compute per‑residue metrics across
different alignment strategies. Each residue in a protein is annotated
with participation flags in global, local and LCS alignments, partner
residue indices, dynamic programming support scores, strong conserved
block membership and gap proximity. The results are returned as a
pandas DataFrame that can be saved to TSV or plotted.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_residue_support(
        seq_len: int,
        seq_str: str,
        method_data: Dict[str, dict],
        window: int = 2,
) -> pd.DataFrame:
    """Compute residue support profiles for a protein across alignment methods.

    Parameters
    ----------
    seq_len : int
        Length of the ungapped protein sequence.
    seq_str : str
        Ungapped protein sequence (1‑letter codes) for the protein.
    method_data : dict
        Dictionary keyed by method name (``"global"``, ``"local"``, ``"lcs"``)
        containing the following items for each method:

        ``a_map``: list[Optional[int]]
            Mapping from alignment columns to residue indices in the protein.
        ``b_map``: list[Optional[int]]
            Mapping from alignment columns to residue indices in the partner
            protein (for partner index extraction).
        ``aligned_a``: str
            The gapped alignment string for this protein.
        ``aligned_b``: str
            The gapped alignment string for the partner protein.
        ``dp_matrix``: np.ndarray, optional
            DP scoring matrix for this method. If provided, support
            scores are extracted from the cell corresponding to each
            residue and partner index. If absent, support scores are
            NaN.
        ``blocks``: pd.DataFrame, optional
            DataFrame of conserved blocks, as returned by
            ``block_detection.detect_blocks_to_df``. Residues inside
            blocks classified as ``high_identity`` or ``conservative``
            are considered part of a strong block.

    window : int, optional
        Number of alignment columns around the residue to consider when
        computing gap proximity and local DP support. A window of 2
        examines ±2 alignment columns.

    Returns
    -------
    pandas.DataFrame
        DataFrame with one row per residue and columns describing
        participation and support in each method. Columns include:
        ``residue_index``, ``residue``, and per‑method features:
        ``<method>_participates``, ``<method>_partner_index``,
        ``<method>_dp_score``, ``<method>_local_support``,
        ``<method>_strong_block``, ``<method>_gap_proximity``.
    """
    # Initialise DataFrame
    df = pd.DataFrame({
        'residue_index': list(range(seq_len)),
        'residue': list(seq_str),
    })

    for method, data in method_data.items():
        a_map: List[Optional[int]] = data.get('a_map')
        b_map: List[Optional[int]] = data.get('b_map')
        aligned_a: str = data.get('aligned_a', '')
        aligned_b: str = data.get('aligned_b', '')
        dp_matrix: Optional[np.ndarray] = data.get('dp_matrix')
        blocks_df: Optional[pd.DataFrame] = data.get('blocks')

        participates = [False] * seq_len
        partner_indices = [None] * seq_len
        dp_scores = [np.nan] * seq_len
        local_support_scores = [np.nan] * seq_len
        strong_block_flags = [False] * seq_len
        gap_proximities = [0] * seq_len

        # Precompute mapping of residue index to alignment positions
        residue_to_aln_positions: Dict[int, List[int]] = {}
        if a_map:
            for aln_idx, res_idx in enumerate(a_map):
                if res_idx is not None:
                    residue_to_aln_positions.setdefault(res_idx, []).append(aln_idx)

        for res_idx in range(seq_len):
            aln_positions = residue_to_aln_positions.get(res_idx, [])
            if aln_positions:
                participates[res_idx] = True
                # choose first alignment column for partner / dp mapping
                aln_col = aln_positions[0]
                # partner index
                if b_map and aln_col < len(b_map):
                    partner_indices[res_idx] = b_map[aln_col]
                # dp score: mapping to DP cell (i+1, j+1) since DP includes zero row/col
                if dp_matrix is not None:
                    i_dp = res_idx + 1
                    j_dp = (partner_indices[res_idx] + 1) if partner_indices[res_idx] is not None else None
                    if j_dp is not None and i_dp < dp_matrix.shape[0] and j_dp < dp_matrix.shape[1]:
                        dp_scores[res_idx] = dp_matrix[i_dp, j_dp]
                # local support: search neighbourhood around dp cell
                if dp_matrix is not None:
                    # compute max score in window around cell; if partner index unknown, use row
                    i_dp = res_idx + 1
                    if partner_indices[res_idx] is not None:
                        j_dp = partner_indices[res_idx] + 1
                        i_start = max(0, i_dp - window)
                        i_end = min(dp_matrix.shape[0], i_dp + window + 1)
                        j_start = max(0, j_dp - window)
                        j_end = min(dp_matrix.shape[1], j_dp + window + 1)
                        local_region = dp_matrix[i_start:i_end, j_start:j_end]
                        if local_region.size > 0:
                            local_support_scores[res_idx] = float(np.max(local_region))
                    else:
                        # no partner: take row slice
                        i_start = max(0, i_dp - window)
                        i_end = min(dp_matrix.shape[0], i_dp + window + 1)
                        local_region = dp_matrix[i_start:i_end, :]
                        local_support_scores[res_idx] = float(np.max(local_region))
                # strong block membership
                if blocks_df is not None and not blocks_df.empty:
                    # check if residue index falls into any block range of seqA_range
                    for _, block in blocks_df.iterrows():
                        if block.get('seqA_range') is not None:
                            start, end = block['seqA_range']
                            # inclusive range in residue coordinates
                            if start <= res_idx <= end and block['classification'] in ('high_identity', 'conservative'):
                                strong_block_flags[res_idx] = True
                                break
                # gap proximity
                if aligned_a and aligned_b and a_map:
                    # examine alignment columns within +/- window of each aln_position
                    gap_count = 0
                    for aln_col in aln_positions:
                        for offset in range(-window, window + 1):
                            col = aln_col + offset
                            if 0 <= col < len(aligned_a):
                                if aligned_a[col] == '-' or aligned_b[col] == '-':
                                    gap_count += 1
                    gap_proximities[res_idx] = gap_count
        # assign columns to dataframe
        df[f"{method}_participates"] = participates
        df[f"{method}_partner_index"] = partner_indices
        df[f"{method}_dp_score"] = dp_scores
        df[f"{method}_local_support"] = local_support_scores
        df[f"{method}_strong_block"] = strong_block_flags
        df[f"{method}_gap_proximity"] = gap_proximities
    return df
