"""Path and traceback utilities for sequence alignment analyses.

Traceback or path files typically contain a series of coordinate
positions in the dynamic programming matrix that were visited during
the reconstruction of the optimal alignment. This module provides
parsing utilities and basic statistics computed on these paths.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


def load_path(path_file: Path) -> List[Tuple[int, int]]:
    """Load a list of DP matrix coordinates from a path file.

    Parameters
    ----------
    path_file : Path
        Path to a text file where each line contains two integers
        separated by whitespace, representing the row and column index
        (0‑based) of a cell in the dynamic programming matrix.

    Returns
    -------
    list[tuple[int, int]]
        Sequence of (row, column) coordinates along the traceback
        path. The order corresponds to the order they appear in the
        file.
    """
    coords: List[Tuple[int, int]] = []
    with path_file.open('r') as fh:
        for line_no, line in enumerate(fh, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            parts = stripped.split()
            if len(parts) != 2:
                logger.warning(
                    "Skipping malformed line %d in %s: %s", line_no, path_file, stripped
                )
                continue
            try:
                i, j = int(parts[0]), int(parts[1])
            except ValueError:
                logger.warning(
                    "Non‑integer coordinate on line %d in %s: %s", line_no, path_file, stripped
                )
                continue
            coords.append((i, j))
    return coords


def validate_path_dimensions(path: List[Tuple[int, int]], shape: Tuple[int, int]) -> None:
    """Ensure all coordinates in a path lie within a DP matrix.

    Parameters
    ----------
    path : list[tuple[int, int]]
        The list of (row, column) coordinates.
    shape : tuple[int, int]
        Shape of the dynamic programming matrix.

    Raises
    ------
    ValueError
        If any coordinate lies outside the matrix bounds.
    """
    rows, cols = shape
    for (i, j) in path:
        if not (0 <= i < rows and 0 <= j < cols):
            raise ValueError(f"Path coordinate ({i}, {j}) is outside matrix shape {shape}")


def overlay_path_on_matrix(matrix: np.ndarray, path: List[Tuple[int, int]]) -> np.ndarray:
    """Create a copy of the DP matrix with the traceback path overlaid.

    The original matrix is not modified. Path positions are set to
    ``np.nan`` to distinguish them from the score landscape. This
    representation can be used for plotting.

    Parameters
    ----------
    matrix : np.ndarray
        DP scoring matrix.
    path : list[tuple[int, int]]
        Traceback path coordinates.

    Returns
    -------
    np.ndarray
        A copy of ``matrix`` where path cells are set to NaN.
    """
    overlay = np.array(matrix, copy=True, dtype=float)
    for (i, j) in path:
        if 0 <= i < overlay.shape[0] and 0 <= j < overlay.shape[1]:
            overlay[i, j] = np.nan
    return overlay


def compute_path_metrics(path: List[Tuple[int, int]]) -> dict:
    """Compute summary metrics describing a traceback path.

    This function attempts to capture simple geometric properties of
    the path, such as the number of diagonal, horizontal and vertical
    moves and the lengths of contiguous runs. These metrics help
    characterise gap runs and diagonal continuity.

    Parameters
    ----------
    path : list[tuple[int, int]]
        The ordered traceback coordinates.

    Returns
    -------
    dict
        Dictionary with metrics including counts of diagonal,
        horizontal and vertical steps, the number and lengths of gap
        runs and the number of direction changes.
    """
    if not path:
        return {
            "num_steps": 0,
            "diagonal_steps": 0,
            "horizontal_steps": 0,
            "vertical_steps": 0,
            "gap_runs": 0,
            "avg_gap_run_length": 0.0,
            "direction_changes": 0,
        }
    diagonal_steps = 0
    horizontal_steps = 0
    vertical_steps = 0
    direction_changes = 0
    gap_run_lengths: List[int] = []
    current_gap_run = 0
    prev_direction: Optional[str] = None
    for (prev_i, prev_j), (i, j) in zip(path, path[1:]):
        di = i - prev_i
        dj = j - prev_j
        if di == 1 and dj == 1:
            step = 'diag'
            diagonal_steps += 1
        elif di == 1 and dj == 0:
            step = 'vert'
            vertical_steps += 1
        elif di == 0 and dj == 1:
            step = 'horiz'
            horizontal_steps += 1
        else:
            # Non unit step: record as change
            step = 'other'
        # Track gap runs: vertical or horizontal contiguous segments
        if step in {'vert', 'horiz'}:
            if prev_direction == step:
                current_gap_run += 1
            else:
                if current_gap_run > 0:
                    gap_run_lengths.append(current_gap_run)
                current_gap_run = 1
        else:
            if current_gap_run > 0:
                gap_run_lengths.append(current_gap_run)
                current_gap_run = 0
        # Direction change detection
        if prev_direction and step != prev_direction:
            direction_changes += 1
        prev_direction = step
    # flush last gap run
    if current_gap_run > 0:
        gap_run_lengths.append(current_gap_run)

    avg_gap_run = float(np.mean(gap_run_lengths)) if gap_run_lengths else 0.0
    return {
        "num_steps": len(path) - 1,
        "diagonal_steps": diagonal_steps,
        "horizontal_steps": horizontal_steps,
        "vertical_steps": vertical_steps,
        "gap_runs": len(gap_run_lengths),
        "avg_gap_run_length": avg_gap_run,
        "direction_changes": direction_changes,
    }
