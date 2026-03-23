"""Alignment method comparison utilities.

This module defines functions to compare residue participation
across global, local and LCS alignment strategies. It supports
constructing category assignments for each residue, summarising
segments of contiguous categories and producing DataFrames for
reporting and plotting.
"""

from __future__ import annotations

import logging
from typing import Dict, List

import numpy as np
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
    is_global = df.get("global_participates", pd.Series(False, index=df.index)).astype(
        bool
    )
    is_local = df.get("local_participates", pd.Series(False, index=df.index)).astype(
        bool
    )
    is_lcs = df.get("lcs_participates", pd.Series(False, index=df.index)).astype(bool)

    conditions = [
        (is_global & ~is_local & ~is_lcs),  # global_only
        (~is_global & is_local & ~is_lcs),  # local_only
        (~is_global & ~is_local & is_lcs),  # lcs_only
        (is_global & is_local & ~is_lcs),  # global_local_shared
        (is_global & ~is_local & is_lcs),  # global_lcs_shared
        (~is_global & is_local & is_lcs),  # local_lcs_shared
        (is_global & is_local & is_lcs),  # all_shared
    ]

    choices = [
        "global_only",
        "local_only",
        "lcs_only",
        "global_local_shared",
        "global_lcs_shared",
        "local_lcs_shared",
        "all_shared",
    ]

    # 3. Use np.select to assign values based on conditions.
    # The 'default' handles the 'none' (False, False, False) case.
    res = np.select(conditions, choices, default="none")

    # 4. Return as a Series with a Categorical dtype for memory efficiency
    return pd.Series(
        pd.Categorical(res, categories=choices + ["none"]),
        index=df.index,
        name="category",
    )


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
            segments.append(
                {
                    "start": start_idx,
                    "end": idx - 1,
                    "category": current_category,
                    "length": (idx - 1) - start_idx + 1,
                }
            )
            current_category = cat
            start_idx = idx
    # close final segment
    if current_category is not None and start_idx is not None:
        end_idx = category_series.index[-1]
        segments.append(
            {
                "start": start_idx,
                "end": end_idx,
                "category": current_category,
                "length": end_idx - start_idx + 1,
            }
        )
    return pd.DataFrame(segments)
