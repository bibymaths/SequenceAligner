"""Summary report generation.

This module aggregates analysis results into a single JSON summary file.
The summary captures metadata about the input files, sequences, dynamic
programming matrices, conserved blocks, alignment statistics and
comparative observations. It also records any warnings or limitations
encountered during analysis.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

logger = logging.getLogger(__name__)


def generate_summary_json(summary_data: Mapping[str, Any], out_path: Path) -> None:
    """Write a summary dictionary to a JSON file.

    Parameters
    ----------
    summary_data : dict
        Arbitrary nested dictionary containing summary information.
    out_path : Path
        Destination path for the JSON file. Parent directories are
        created if necessary.

    Returns
    -------
    None
        The summary is written to disk.
    """
    logger.info("Writing summary JSON to %s", out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with out_path.open("w") as fh:
            json.dump(summary_data, fh, indent=2)
    except Exception as exc:
        logger.error("Failed to write summary JSON: %s", exc)


def build_summary_data(
    input_files: Mapping[str, str],
    sequence_ids: Tuple[str, str],
    sequence_lengths: Tuple[int, int],
    dp_shapes: Mapping[str, Tuple[int, int]],
    stats_metadata: Mapping[str, Any],
    blocks_top: Mapping[str, List[Mapping[str, Any]]],
    alignment_stats: Mapping[str, Mapping[str, float]],
    category_counts: Optional[Mapping[str, Any]] = None,
    warnings: Optional[List[str]] = None,
    notes: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Assemble summary information into a dictionary.

    Parameters
    ----------
    input_files : dict
        Mapping of input file descriptors to file paths used in analysis.
    sequence_ids : tuple[str, str]
        Identifiers of the two aligned proteins.
    sequence_lengths : tuple[int, int]
        Ungapped lengths of the two proteins.
    dp_shapes : mapping[str, tuple[int, int]]
        Dimensions of DP matrices for each alignment method.
    stats_metadata : mapping[str, any]
        Parsed contents of stats JSON files keyed by method ('global',
        'local'). LCS may not have such metadata.
    blocks_top : mapping[str, list]
        Top conserved blocks per method. Each list contains block
        dictionaries (as returned from block detection) sorted by
        decreasing identity fraction.
    alignment_stats : mapping[str, mapping[str, float]]
        Alignment summary statistics for each method (keys like
        'alignment_length', 'percent_identity', etc.).
    category_counts : mapping, optional
        Counts of residues falling into each participation category.
    warnings : list[str], optional
        Any warnings encountered during analysis.
    notes : list[str], optional
        Additional notes or limitations.

    Returns
    -------
    dict
        Structured summary dictionary.
    """
    summary: Dict[str, Any] = {
        "input_files": dict(input_files),
        "sequence_ids": list(sequence_ids),
        "sequence_lengths": list(sequence_lengths),
        "dp_shapes": {k: list(v) for k, v in dp_shapes.items()},
        "stats_metadata": stats_metadata,
        "top_blocks": {},
    }
    # Top blocks
    for method, blocks in blocks_top.items():
        summary["top_blocks"][method] = blocks
    summary["alignment_stats"] = alignment_stats
    if category_counts is not None:
        summary["participation_counts"] = category_counts
    if warnings:
        summary["warnings"] = list(warnings)
    if notes:
        summary["notes"] = list(notes)
    return summary
