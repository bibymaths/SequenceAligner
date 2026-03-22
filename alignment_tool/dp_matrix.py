"""Dynamic programming matrix loading and validation.

Functions in this module facilitate loading scoring matrices generated
by global, local and LCS alignment algorithms. Where both binary and
text representations exist, the binary file is preferred for
performance using NumPy memory mapping.

References
==========
Dynamic programming for sequence alignment constructs a table where
rows correspond to positions in the first sequence and columns
correspond to positions in the second sequence【489419667959862†L180-L205】.
The table has dimensions (m+1)×(n+1) for sequences of length m and n,
respectively, with the extra row and column representing the empty
prefix.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


def infer_shape(seq_a_len: int, seq_b_len: int) -> Tuple[int, int]:
    """Infer the dynamic programming matrix shape for two sequences.

    According to standard alignment algorithms, the DP matrix has
    dimensions (len(seq_a) + 1, len(seq_b) + 1)【489419667959862†L180-L205】.

    Parameters
    ----------
    seq_a_len : int
        Length of the first (ungapped) sequence.
    seq_b_len : int
        Length of the second (ungapped) sequence.

    Returns
    -------
    tuple[int, int]
        The inferred matrix shape.
    """
    return (seq_a_len + 1, seq_b_len + 1)


def load_dp_matrix(
        bin_path: Optional[Path],
        txt_path: Optional[Path],
        shape: Tuple[int, int],
        dtype: str = "float64",
) -> np.ndarray:
    """Load a DP matrix from binary or text representation.

    Parameters
    ----------
    bin_path : Path, optional
        Path to binary file containing the matrix in row‑major order.
    txt_path : Path, optional
        Path to text file containing whitespace separated rows of
        scores. Used as a fallback if the binary file is missing.
    shape : tuple[int, int]
        The expected shape of the matrix. This must match the number
        of elements in the file.
    dtype : str, optional
        The NumPy data type to use when reading the binary file. The
        default of ``float64`` should accommodate most scoring
        matrices. If the binary file was created with 32‑bit ints, you
        should pass ``dtype='int32'``.

    Returns
    -------
    np.ndarray
        The loaded DP matrix with shape ``shape``. If a binary file
        exists, the returned array may be a ``memmap`` backed by the
        file; otherwise, a regular NumPy array loaded from the text file.

    Raises
    ------
    FileNotFoundError
        If neither binary nor text file is provided.
    ValueError
        If the loaded data does not match the expected shape.
    """
    if bin_path and bin_path.exists():
        logger.debug("Loading DP matrix from binary file %s", bin_path)
        # Use memmap to avoid loading entire matrix into memory when not needed
        total_elements = shape[0] * shape[1]
        file_size = bin_path.stat().st_size
        dtype_np = np.dtype(dtype)
        expected_bytes = dtype_np.itemsize * total_elements
        if file_size < expected_bytes:
            logger.warning(
                "Binary DP matrix %s is smaller than expected (%d bytes < %d bytes)",
                bin_path, file_size, expected_bytes,
            )
        try:
            memmap = np.memmap(bin_path, dtype=dtype_np, mode='r', shape=shape)
        except Exception as e:
            logger.error("Failed to memory map %s: %s", bin_path, e)
            raise
        return memmap
    elif txt_path and txt_path.exists():
        logger.debug("Loading DP matrix from text file %s", txt_path)
        data = np.loadtxt(txt_path, dtype=float)
        if data.shape != shape:
            raise ValueError(
                f"Text DP matrix at {txt_path} has shape {data.shape}, expected {shape}"
            )
        return data
    else:
        raise FileNotFoundError(
            f"No DP matrix file found. Checked binary: {bin_path}, text: {txt_path}"
        )


def validate_matrix_shape(matrix: np.ndarray, shape: Tuple[int, int]) -> None:
    """Validate that a DP matrix has the expected dimensions.

    Parameters
    ----------
    matrix : np.ndarray
        The loaded matrix to validate.
    shape : tuple[int, int]
        The expected shape.

    Raises
    ------
    ValueError
        If the matrix shape does not match the expected shape.
    """
    if matrix.shape != shape:
        raise ValueError(
            f"Matrix shape mismatch: expected {shape}, got {matrix.shape}"
        )
