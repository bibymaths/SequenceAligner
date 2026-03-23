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
    """
    Load a DP matrix from binary or text representation with support for int32 headers.

    This function attempts, in order:
    1. int32 header (rows, cols) + int32 matrix
    2. raw int32
    3. raw float32
    4. raw float64
    5. fallback to txt

    If header is detected, it OVERRIDES the provided shape.

    Notes
    -----
    - Header format assumed: first 8 bytes = two int32 values (rows, cols)
    - Matrix stored in row-major order
    """

    if bin_path and bin_path.exists():
        logger.debug("Loading DP matrix from binary file %s", bin_path)

        file_size = bin_path.stat().st_size

        # --------------------------------------------------
        # 1. Try int32 header (rows, cols) + int32 matrix
        # --------------------------------------------------
        try:
            with open(bin_path, "rb") as fh:
                header = np.fromfile(fh, dtype=np.int32, count=2)

            if len(header) == 2:
                rows, cols = int(header[0]), int(header[1])

                if rows > 0 and cols > 0:
                    expected_bytes = 8 + rows * cols * 4

                    if expected_bytes == file_size:
                        logger.info(
                            "Detected int32 header: shape=(%d, %d)", rows, cols
                        )

                        return np.memmap(
                            bin_path,
                            dtype=np.int32,
                            mode="r",
                            offset=8,
                            shape=(rows, cols),
                        )
        except Exception as e:
            logger.debug("Header detection failed: %s", e)

        # --------------------------------------------------
        # 2. Try raw int32
        # --------------------------------------------------
        if file_size % 4 == 0:
            total_elements = file_size // 4
            if total_elements == shape[0] * shape[1]:
                logger.info("Loading raw int32 matrix with expected shape")
                return np.memmap(
                    bin_path,
                    dtype=np.int32,
                    mode="r",
                    shape=shape,
                )

        # --------------------------------------------------
        # 3. Try raw float32
        # --------------------------------------------------
        if file_size % 4 == 0:
            total_elements = file_size // 4
            if total_elements == shape[0] * shape[1]:
                logger.info("Loading raw float32 matrix with expected shape")
                return np.memmap(
                    bin_path,
                    dtype=np.float32,
                    mode="r",
                    shape=shape,
                )

        # --------------------------------------------------
        # 4. Try raw float64 (your original assumption)
        # --------------------------------------------------
        if file_size % 8 == 0:
            total_elements = file_size // 8
            if total_elements == shape[0] * shape[1]:
                logger.info("Loading raw float64 matrix with expected shape")
                return np.memmap(
                    bin_path,
                    dtype=np.float64,
                    mode="r",
                    shape=shape,
                )

        # --------------------------------------------------
        # If nothing worked → error
        # --------------------------------------------------
        logger.error(
            "Failed to infer binary DP matrix format for %s (file size: %d bytes)",
            bin_path, file_size
        )
        raise ValueError(f"Unsupported or inconsistent DP matrix format: {bin_path}")

    # --------------------------------------------------
    # 5. Fallback to text
    # --------------------------------------------------
    elif txt_path and txt_path.exists():
        logger.debug("Loading DP matrix from text file %s", txt_path)

        data = np.loadtxt(txt_path, dtype=float)

        if data.shape != shape:
            logger.warning(
                "Text matrix shape %s does not match expected %s — using text shape",
                data.shape, shape
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
