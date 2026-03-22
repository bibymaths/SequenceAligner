"""Substitution‑aware analyses for protein alignments.

This module provides utilities to interpret residue substitutions in
aligned protein sequences. It supports loading the BLOSUM62
substitution matrix via Biopython, classifying amino acids by
physicochemical categories, distinguishing identical, conservative
and radical substitutions, and summarising these properties across
an alignment. These functions can be used to generate TSV summaries
or feed into higher level analyses.
"""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

import pandas as pd
from Bio.Align import substitution_matrices

logger = logging.getLogger(__name__)


def load_substitution_matrix(name: str | None) -> Optional[Mapping[str, Mapping[str, float]]]:
    """Load a substitution matrix by name.

    Parameters
    ----------
    name : str or None
        Name of the substitution matrix to load. Currently only
        ``"blosum62"`` (case insensitive) is supported. If ``None``
        or ``"none"``, no matrix is loaded and the function returns
        ``None``.

    Returns
    -------
    dict or None
        A nested dictionary mapping residue letters to their scores
        against all other residues. If the matrix cannot be loaded,
        ``None`` is returned.

    Notes
    -----
    The Biopython ``substitution_matrices`` module is used to load
    the BLOSUM62 matrix. Biopython represents matrices as objects
    that behave like dictionaries. This function converts the matrix
    to a plain dictionary of dictionaries for easier use.
    """
    if not name or name.lower() in {"none", "null", ""}:
        logger.info("No substitution matrix requested")
        return None
    name_lower = name.lower()
    if name_lower == "blosum62":
        try:
            matrix = substitution_matrices.load("BLOSUM62")
        except Exception as exc:
            logger.error("Failed to load BLOSUM62 matrix: %s", exc)
            return None
        # Convert to nested dict
        sub_dict: Dict[str, Dict[str, float]] = {}
        for aa in matrix.alphabet:
            sub_dict[aa] = {}
            for bb in matrix.alphabet:
                sub_dict[aa][bb] = float(matrix[aa, bb])
        return sub_dict
    else:
        logger.warning("Unsupported substitution matrix name: %s", name)
        return None


def classify_residue(residue: str) -> List[str]:
    """Assign physicochemical categories to an amino acid.

    Parameters
    ----------
    residue : str
        A single character representing an amino acid.

    Returns
    -------
    list[str]
        A list of categories. Possible categories include:
        ``"glycine"``, ``"proline"``, ``"cysteine"``, ``"aromatic"``,
        ``"positive"``, ``"negative"``. An empty list is returned
        for unknown or uncategorised residues.

    References
    ----------
    Amino acids can be grouped into classes based on their
    physicochemical properties. For example, glycine (G) and proline
    (P) have unique structural roles, cysteine (C) contains sulfur,
    aromatic residues include phenylalanine (F), tryptophan (W) and
    tyrosine (Y), and charged residues include arginine (R), histidine
    (H), lysine (K), aspartic acid (D) and glutamic acid (E)【915331229988753†L39-L50】【915331229988753†L111-L113】.
    """
    residue = residue.upper()
    categories: List[str] = []
    if residue == "G":
        categories.append("glycine")
    if residue == "P":
        categories.append("proline")
    if residue == "C":
        categories.append("cysteine")
    if residue in {"F", "W", "Y"}:
        categories.append("aromatic")
    if residue in {"R", "H", "K"}:
        categories.append("positive")
    if residue in {"D", "E"}:
        categories.append("negative")
    return categories


def summarise_substitutions(
        seq_a: str,
        seq_b: str,
        substitution_matrix: Optional[Mapping[str, Mapping[str, float]]],
        similarity_threshold: float = 0.0,
) -> pd.DataFrame:
    """Summarise substitution types and residue categories across an alignment.

    Parameters
    ----------
    seq_a, seq_b : str
        Aligned protein sequences (including gaps). Must be of equal
        length.
    substitution_matrix : dict or None
        Nested dictionary of substitution scores. If provided, scores
        >= ``similarity_threshold`` are considered conservative
        substitutions. If ``None``, only exact matches are counted and
        all non‑identical pairs are classified as ``"radical"``.
    similarity_threshold : float, optional
        Score threshold above which a substitution is considered
        conservative. Default is 0.0.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with columns summarising counts of identical,
        conservative and radical substitutions, number of gap pairs,
        and counts of notable categories (glycine, proline, cysteine,
        aromatic, positive, negative) among identical and conservative
        substitutions.

    Notes
    -----
    Gap positions (where either ``seq_a`` or ``seq_b`` has a '-' )
    are excluded from substitution classification and counted under
    ``gap_pairs``. Positions involving unknown residues that are not
    present in the substitution matrix are treated as radical
    substitutions.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("Aligned sequences must be of equal length")
    # Counters
    counts = {
        "identical": 0,
        "conservative": 0,
        "radical": 0,
        "gap_pairs": 0,
    }
    category_counts = {
        "identical": {
            "glycine": 0,
            "proline": 0,
            "cysteine": 0,
            "aromatic": 0,
            "positive": 0,
            "negative": 0,
        },
        "conservative": {
            "glycine": 0,
            "proline": 0,
            "cysteine": 0,
            "aromatic": 0,
            "positive": 0,
            "negative": 0,
        },
    }
    for aa, bb in zip(seq_a.upper(), seq_b.upper()):
        if aa == '-' or bb == '-':
            counts["gap_pairs"] += 1
            continue
        if aa == bb:
            counts["identical"] += 1
            for cat in classify_residue(aa):
                category_counts["identical"][cat] += 1
        else:
            if substitution_matrix is not None:
                score = substitution_matrix.get(aa, {}).get(bb, None)
                if score is not None and score >= similarity_threshold:
                    counts["conservative"] += 1
                    for cat in set(classify_residue(aa) + classify_residue(bb)):
                        category_counts["conservative"][cat] += 1
                else:
                    counts["radical"] += 1
            else:
                counts["radical"] += 1
    # Build DataFrame
    data = {
        "metric": ["identical", "conservative", "radical", "gap_pairs"],
        "count": [
            counts["identical"],
            counts["conservative"],
            counts["radical"],
            counts["gap_pairs"],
        ],
    }
    df = pd.DataFrame(data)
    # Append category counts as separate columns
    for group in ["identical", "conservative"]:
        for cat, val in category_counts[group].items():
            df.loc[df["metric"] == group, f"{group}_{cat}"] = val
    # Fill NaN with 0
    df = df.fillna(0)
    return df
