"""Utilities specific to longest common subsequence (LCS) analysis.

This module contains functions to parse the DP length matrix and
traceback pointers used for computing the LCS. It also provides
helpers for interpreting the LCS within the context of the aligned
proteins.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from .dp_matrix import load_dp_matrix

logger = logging.getLogger(__name__)


def load_lcs_dp_lengths(
    bin_path: Optional[Path],
    txt_path: Optional[Path],
    shape: Tuple[int, int],
    dtype: str = "int32",
) -> np.ndarray:
    """Load the DP matrix of LCS lengths.

    Parameters
    ----------
    bin_path, txt_path : Path, optional
        Locations of the binary and text files storing the LCS length
        matrix. At least one must exist.
    shape : tuple[int, int]
        Expected shape of the matrix.
    dtype : str, optional
        Data type of the matrix. LCS lengths are integers and thus
        default to ``int32``.

    Returns
    -------
    np.ndarray
        The loaded DP matrix.
    """
    return load_dp_matrix(bin_path, txt_path, shape, dtype=dtype)


def load_traceback_pointers(path: Path, shape: Tuple[int, int]) -> np.ndarray:
    """Load an LCS traceback pointer matrix from text.

    The pointer symbols are expected to be one of 'D' (diagonal), 'U'
    (up) or 'L' (left) and correspond to the direction of the next
    pointer when performing the traceback.

    Parameters
    ----------
    path : Path
        Path to the pointer matrix in text form.
    shape : tuple[int, int]
        The expected shape of the pointer matrix.

    Returns
    -------
    np.ndarray
        A 2‑D array of strings with shape ``shape`` containing
        direction symbols. Missing entries are filled with empty
        strings.
    """
    rows, cols = shape
    pointers = np.full(shape, "", dtype=object)
    with path.open("r") as fh:
        for r, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            for c, symbol in enumerate(parts):
                if r < rows and c < cols:
                    pointers[r, c] = symbol
    return pointers


def traceback_lcs(
    pointers: np.ndarray,
    seq_a: str,
    seq_b: str,
) -> List[Tuple[int, int]]:
    """Reconstruct the LCS path from a pointer matrix.

    Parameters
    ----------
    pointers : np.ndarray
        LCS traceback pointer matrix with elements 'D', 'U', 'L'.
    seq_a, seq_b : str
        Ungapped sequences. Used to determine matrix dimensions.

    Returns
    -------
    list[tuple[int, int]]
        List of matrix coordinates visited during the traceback,
        starting from (len(seq_a), len(seq_b)) and moving towards (0,0).

    Notes
    -----
    This routine follows the pointers backwards and stops when both
    indices reach 0. If multiple pointers exist (ties), this simple
    implementation prefers the order D, U, L. If pointers are missing,
    it falls back to decrementing indices in favour of diagonal moves.
    """
    i, j = len(seq_a), len(seq_b)
    path: List[Tuple[int, int]] = [(i, j)]
    while i > 0 or j > 0:
        symbol = (
            pointers[i, j] if (i < pointers.shape[0] and j < pointers.shape[1]) else ""
        )
        if symbol == "D":
            i -= 1
            j -= 1
        elif symbol == "U":
            i -= 1
        elif symbol == "L":
            j -= 1
        else:
            # fallback: prefer diagonal if both indices > 0
            if i > 0 and j > 0:
                i -= 1
                j -= 1
            elif i > 0:
                i -= 1
            elif j > 0:
                j -= 1
        path.append((i, j))
    path.reverse()
    return path
