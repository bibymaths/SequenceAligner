"""Alignment method comparison utilities.

This module defines functions to compare residue participation
across global, local and LCS alignment strategies. It supports
constructing category assignments for each residue, summarising
segments of contiguous categories and producing DataFrames for
reporting and plotting.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def assign_participation_categories(df: pd.DataFrame) -> pd.Series:
    """Assign participation categories to each residue.

    Given a residue support DataFrame containing boolean columns
    ``global_participates``, ``local_participates`` and
    ``lcs_participates``, this function creates a categorical series
    representing the combination of methods in which each residue is
    aligned. The categories are named as follows:

    - ``global_only``: residue participates only in global alignment.
    - ``local_only``: residue participates only in local alignment.
    - ``lcs_only``: residue participates only in LCS.
    - ``global_local_shared``: participates in global and local but not LCS.
    - ``global_lcs_shared``: participates in global and LCS but not local.
    - ``local_lcs_shared``: participates in local and LCS but not global.
    - ``all_shared``: participates in all three methods.
    - ``none``: participates in none of the methods.

    Parameters
    ----------
    df : pandas.DataFrame
        Residue support DataFrame with boolean participation columns.

    Returns
    -------
    pandas.Series
        Categorical series with the category for each residue index.
    """
    categories = []
    for _, row in df.iterrows():
        g = bool(row.get("global_participates", False))
        l = bool(row.get("local_participates", False))
        s = bool(row.get("lcs_participates", False))
        if g and not l and not s:
            categories.append("global_only")
        elif not g and l and not s:
            categories.append("local_only")
        elif not g and not l and s:
            categories.append("lcs_only")
        elif g and l and not s:
            categories.append("global_local_shared")
        elif g and not l and s:
            categories.append("global_lcs_shared")
        elif not g and l and s:
            categories.append("local_lcs_shared")
        elif g and l and s:
            categories.append("all_shared")
        else:
            categories.append("none")
    return pd.Series(categories, index=df.index, name="category")


def summarise_category_segments(category_series: pd.Series) -> pd.DataFrame:
    """Summarise contiguous segments of identical categories.

    This helper takes a series of category labels indexed by residue
    index and merges consecutive positions with the same category into
    segments. It returns a DataFrame with start and end indices
    (inclusive), the category label and the length of each segment.

    Parameters
    ----------
    category_series : pandas.Series
        Series mapping residue index to category label.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ``start``, ``end``, ``category`` and
        ``length`` summarising consecutive segments.
    """
    segments: List[Dict[str, int | str]] = []
    current_category: str = None  # type: ignore[assignment]
    start_idx: int = None  # type: ignore[assignment]
    for idx, cat in category_series.items():
        if current_category is None:
            current_category = cat
            start_idx = idx
        elif cat != current_category:
            # close previous segment
            segments.append({
                "start": start_idx,
                "end": idx - 1,
                "category": current_category,
                "length": (idx - 1) - start_idx + 1,
            })
            current_category = cat
            start_idx = idx
    # close final segment
    if current_category is not None and start_idx is not None:
        end_idx = category_series.index[-1]
        segments.append({
            "start": start_idx,
            "end": end_idx,
            "category": current_category,
            "length": end_idx - start_idx + 1,
        })
    return pd.DataFrame(segments)
