"""High‑level block detection and classification utilities.

This module provides a thin abstraction over the contiguous block
functions in ``fasta_utils`` to generate DataFrames and summarise
conserved blocks for output. It allows direct detection of blocks
using simple parameters, writing results to TSV if required.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from . import fasta_utils

logger = logging.getLogger(__name__)


def detect_blocks_to_df(
        seq_a: str,
        seq_b: str,
        a_map: List[Optional[int]],
        b_map: List[Optional[int]],
        substitution_matrix: Dict[str, Dict[str, int]],
        min_block_length: int,
        identity_threshold: float,
        similarity_threshold: float,
) -> pd.DataFrame:
    """Detect conserved blocks and return a DataFrame of results.

    Parameters
    ----------
    seq_a, seq_b, a_map, b_map : see ``fasta_utils.contiguous_blocks``.
    substitution_matrix : dict
        Substitution scores used for similarity classification.
    min_block_length, identity_threshold, similarity_threshold : see
        ``fasta_utils.contiguous_blocks``.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with one row per conserved block and the
        statistics computed for each block.
    """
    blocks = fasta_utils.contiguous_blocks(
        seq_a,
        seq_b,
        a_map,
        b_map,
        substitution_matrix,
        min_block_length=min_block_length,
        identity_threshold=identity_threshold,
        similarity_threshold=similarity_threshold,
    )
    df = pd.DataFrame(blocks)
    return df
