"""
Binary matrix parser and downsampling utilities.

The C++ alignment tools write dynamic programming (DP) matrices to
binary files (``.bin``). These matrices are typically stored in
row‑major order with 32‑bit integers representing scores. To support
visualization in the web UI without overwhelming the client, we
downsample matrices larger than 1000×1000 by selecting a grid of
rows/columns. The full matrix is preserved on disk for offline
analysis or high‑resolution plotting.

This module uses NumPy if available for efficient parsing; otherwise
falls back to Python's built‑in array module.
"""

import math
from pathlib import Path
from typing import List, Tuple
import numpy as np


def parse_bin_matrix(file_path: Path, shape: Tuple[int, int]) -> List[List[int]]:
    """Read a binary DP matrix file into a 2D list using the provided shape."""
    data = np.fromfile(file_path, dtype=np.int32)
    if data.size != shape[0] * shape[1]:
        raise ValueError("Size mismatch for DP matrix")
    return data.reshape(shape).tolist()


def downsample_matrix(matrix: List[List[int]], max_dim: int = 1000) -> List[List[int]]:
    """Downsample a 2D list by selecting rows and columns at regular intervals."""
    n_rows = len(matrix)
    n_cols = len(matrix[0]) if n_rows > 0 else 0
    if n_rows <= max_dim and n_cols <= max_dim:
        return matrix
    row_step = math.ceil(n_rows / max_dim)
    col_step = math.ceil(n_cols / max_dim)
    return [row[::col_step] for row in matrix[::row_step]]
