"""Plotting utilities for alignment analysis.

This module contains functions for generating visual representations of
dynamic programming matrices, residue support profiles and conserved
blocks. The plots are saved to files using matplotlib. Each
function accepts an output path and optional DPI setting to allow
customisation of image resolution.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def plot_dp_heatmap(
    matrix: np.ndarray,
    out_path: Path,
    path_coords: Optional[List[Tuple[int, int]]] = None,
    title: Optional[str] = None,
    dpi: int = 150,
) -> None:
    """Plot a heatmap of a dynamic programming matrix.

    Parameters
    ----------
    matrix : np.ndarray
        The DP scoring matrix to visualise. Values are displayed
        using a diverging colour map. NaN values are masked.
    out_path : Path
        Destination file for the PNG image.
    path_coords : list of tuple[int, int], optional
        A list of (row, column) coordinates representing a traceback
        path. If provided, the path is overlayed on the heatmap as
        a white line. Coordinates should be 0‑based and refer to
        matrix indices.
    title : str, optional
        Title for the plot.
    dpi : int, optional
        Output resolution in dots per inch.

    Returns
    -------
    None
        The plot is saved to ``out_path``. Any existing file will be
        overwritten.
    """
    logger.debug("Plotting DP heatmap to %s", out_path)
    # Mask NaNs for colormap
    data = np.array(matrix)
    masked = np.ma.masked_invalid(data)
    # Choose a colour map that highlights both high and low values
    cmap = plt.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(8, 6), dpi=dpi)
    im = ax.imshow(masked, aspect="auto", origin="lower", cmap=cmap)
    ax.set_xlabel("Sequence B index")
    ax.set_ylabel("Sequence A index")
    if title:
        ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    # Overlay path if provided
    if path_coords:
        # path_coords is list of (i,j); convert to arrays for plotting
        rows, cols = zip(*path_coords)
        ax.plot(cols, rows, color="white", linewidth=1)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def plot_residue_support(
    df: pd.DataFrame,
    methods: List[str],
    out_path: Path,
    title: Optional[str] = None,
    dpi: int = 150,
) -> None:
    """Plot residue support profiles for multiple alignment methods.

    The resulting figure contains one row per method, with panels
    showing DP scores, local support scores, strong block membership
    and gap proximity along the sequence. The x‑axis corresponds to
    residue index. For methods lacking certain metrics, missing
    values are plotted as NaNs.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame returned from ``residue_profiles.compute_residue_support``.
    methods : list[str]
        List of method names to include (e.g. ["global", "local", "lcs"]).
    out_path : Path
        Destination file for the figure.
    title : str, optional
        Overall title for the figure.
    dpi : int, optional
        Figure resolution.

    Returns
    -------
    None
        The plot is saved to ``out_path``.
    """
    num_methods = len(methods)
    fig, axes = plt.subplots(
        nrows=num_methods, ncols=4, figsize=(12, 3 * num_methods), dpi=dpi
    )
    if num_methods == 1:
        axes = np.expand_dims(axes, axis=0)
    x = df["residue_index"]
    for i, method in enumerate(methods):
        # DP score plot
        ax = axes[i, 0]
        y = df[f"{method}_dp_score"].astype(float)
        ax.plot(x, y, color="tab:blue")
        ax.set_ylabel(f"{method} DP score")
        ax.set_xlabel("Residue index")
        # Local support plot
        ax = axes[i, 1]
        y = df[f"{method}_local_support"].astype(float)
        ax.plot(x, y, color="tab:orange")
        ax.set_ylabel(f"{method} local support")
        ax.set_xlabel("Residue index")
        # Strong block membership
        ax = axes[i, 2]
        y = df[f"{method}_strong_block"].astype(int)
        ax.bar(x, y, color="tab:green", width=1.0)
        ax.set_ylabel(f"{method} strong block")
        ax.set_xlabel("Residue index")
        ax.set_ylim(-0.05, 1.05)
        # Gap proximity
        ax = axes[i, 3]
        y = df[f"{method}_gap_proximity"].astype(int)
        ax.plot(x, y, color="tab:red")
        ax.set_ylabel(f"{method} gap proximity")
        ax.set_xlabel("Residue index")
    # Set column titles
    col_titles = ["DP score", "Local support", "Strong block", "Gap proximity"]
    for j, col_title in enumerate(col_titles):
        axes[0, j].set_title(col_title)
    if title:
        fig.suptitle(title, y=1.02, fontsize=14)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def plot_conserved_blocks_comparison(
    blocks_dict: Mapping[str, pd.DataFrame],
    seq_length: int,
    out_path: Path,
    title: Optional[str] = None,
    dpi: int = 150,
) -> None:
    """Compare conserved blocks across alignment methods.

    This function creates a plot with one horizontal track per method.
    Each track spans the length of the ungapped sequence and highlights
    contiguous blocks classified as high identity or conservative. High
    identity blocks are coloured dark green and conservative blocks
    light green. Regions without blocks remain white.

    Parameters
    ----------
    blocks_dict : mapping[str, pandas.DataFrame]
        Mapping from method name to DataFrame of conserved blocks as
        returned by ``block_detection.detect_blocks_to_df``. Only
        blocks with a defined ``seqA_range`` are plotted.
    seq_length : int
        Length of the ungapped protein sequence used for the x‑axis.
    out_path : Path
        Path to save the image.
    title : str, optional
        Title for the figure.
    dpi : int, optional
        Resolution of the figure.

    Returns
    -------
    None
        The plot is saved to ``out_path``.
    """
    methods = list(blocks_dict.keys())
    num_methods = len(methods)
    fig, ax = plt.subplots(figsize=(12, 1.5 * num_methods), dpi=dpi)
    # Draw base rectangles
    for idx, method in enumerate(methods):
        y_base = num_methods - idx - 1
        ax.hlines(y_base, 0, seq_length, color="lightgray", linewidth=8)
        df = blocks_dict[method]
        if df is not None and not df.empty:
            for _, row in df.iterrows():
                a_range = row["seqA_range"]
                if not a_range:
                    continue
                start, end = a_range
                classification = row["classification"]
                if classification == "high_identity":
                    colour = "#006400"  # dark green
                elif classification == "conservative":
                    colour = "#66c2a5"  # light green
                else:
                    colour = "#cccccc"  # grey for mismatch rich
                ax.hlines(y_base, start, end + 1, color=colour, linewidth=8)
        ax.text(seq_length + 1, num_methods - idx - 1, method, va="center")
    ax.set_ylim(-1, num_methods)
    ax.set_xlim(0, seq_length + 5)
    ax.set_yticks([])
    ax.set_xlabel("Residue index (sequence A)")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def plot_alignment_method_comparison(
    category_series: pd.Series,
    out_path: Path,
    title: Optional[str] = None,
    dpi: int = 150,
) -> None:
    """Visualise alignment participation categories along the sequence.

    The input series should assign a category label (e.g. 'global_only',
    'local_only', 'global_local_shared', etc.) to each residue index. This
    function plots a coloured bar for each category along the x‑axis.

    Parameters
    ----------
    category_series : pandas.Series
        Series indexed by residue_index with categorical values.
    out_path : Path
        Location to save the figure.
    title : str, optional
        Title for the plot.
    dpi : int, optional
        Image resolution.

    Returns
    -------
    None
        The plot is saved to ``out_path``.
    """
    # Assign a colour to each category
    categories = category_series.unique()
    palette = {
        "global_only": "#1f77b4",
        "local_only": "#ff7f0e",
        "lcs_only": "#2ca02c",
        "global_local_shared": "#9467bd",
        "global_lcs_shared": "#17becf",
        "local_lcs_shared": "#e377c2",
        "all_shared": "#8c564b",
        "none": "#7f7f7f",
    }
    colours = [palette.get(cat, "#cccccc") for cat in category_series]
    fig, ax = plt.subplots(figsize=(12, 2), dpi=dpi)
    x = category_series.index
    ax.bar(x, [1] * len(x), color=colours, width=1.0)
    ax.set_yticks([])
    ax.set_xlim(min(x), max(x))
    ax.set_xlabel("Residue index")
    # Build legend
    handles = []
    labels = []
    for cat in categories:
        if cat in palette:
            handles.append(plt.Rectangle((0, 0), 1, 1, color=palette[cat]))
            labels.append(cat)
    ax.legend(handles, labels, bbox_to_anchor=(1.01, 1), loc="upper left")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
