"""Inventory management for alignment output files.

This module defines data structures and helper functions for
discovering and validating the set of files produced by a pairwise
alignment pipeline. It assumes a fixed naming scheme corresponding
to the files described in the problem specification. Because this tool
is designed around a specific set of inputs, missing files are
considered errors unless optional according to the analysis being
performed.

Functions
=========

scan_results_dir(results_dir: Path) -> AlignmentFiles
    Traverse a directory and record the presence of expected files.

validate_files(alignment_type: str, files: AlignmentFiles) -> None
    Validate that all required files exist for the chosen alignment
    type. Raises a ``FileNotFoundError`` if mandatory files are
    missing.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class AlignmentFiles:
    """Container for all alignment related file paths.

    Attributes correspond to the various alignment files available.
    Not all fields need to be populated for every analysis. Binary
    versions take precedence over text versions where both exist.
    """

    # Alignment FASTA files
    global_alignment: Optional[Path] = None
    local_alignment: Optional[Path] = None
    lcs_alignment: Optional[Path] = None
    lcs: Optional[Path] = None

    # DP matrix files
    global_dp_bin: Optional[Path] = None
    global_dp_txt: Optional[Path] = None
    local_dp_bin: Optional[Path] = None
    local_dp_txt: Optional[Path] = None
    lcs_dp_bin: Optional[Path] = None
    lcs_dp_txt: Optional[Path] = None

    # Path coordinate files
    global_path: Optional[Path] = None
    local_path: Optional[Path] = None
    lcs_path: Optional[Path] = None

    # Traceback pointers for LCS
    lcs_traceback_bin: Optional[Path] = None
    lcs_traceback_txt: Optional[Path] = None

    # Statistics JSON
    global_stats: Optional[Path] = None
    local_stats: Optional[Path] = None

    def available_alignment_types(self) -> list[str]:
        """Return a list of alignment types detected in the inventory.

        Returns
        -------
        list[str]
            A list containing any of ``"global"``, ``"local"`` and
            ``"lcs"`` for which the requisite alignment FASTA exists.
        """
        types = []
        if self.global_alignment:
            types.append("global")
        if self.local_alignment:
            types.append("local")
        if self.lcs_alignment:
            types.append("lcs")
        return types


def scan_results_dir(results_dir: Path) -> AlignmentFiles:
    """Scan a directory for alignment output files.

    Parameters
    ----------
    results_dir : Path
        Path to the directory that contains the alignment outputs.

    Returns
    -------
    AlignmentFiles
        An instance describing the discovered files. Missing files are
        represented by ``None``.

    Notes
    -----
    Only the filenames specified in the problem statement are
    considered. Unexpected files are ignored.
    """
    results_dir = results_dir.expanduser().resolve()
    files = AlignmentFiles()
    if not results_dir.is_dir():
        raise FileNotFoundError(f"Results directory {results_dir} does not exist")

    # Map of expected filenames to object attributes
    filename_map = {
        "global_alignment.fasta": "global_alignment",
        "local_alignment.fasta": "local_alignment",
        "lcs_alignment.fasta": "lcs_alignment",
        "lcs.fasta": "lcs",
        "global_dp_matrix.bin": "global_dp_bin",
        "global_dp_matrix.txt": "global_dp_txt",
        "local_dp_matrix.bin": "local_dp_bin",
        "local_dp_matrix.txt": "local_dp_txt",
        "lcs_dp_lengths.bin": "lcs_dp_bin",
        "lcs_dp_lengths.txt": "lcs_dp_txt",
        "global_path.txt": "global_path",
        "local_path.txt": "local_path",
        "lcs_path.txt": "lcs_path",
        "lcs_traceback_pointers.bin": "lcs_traceback_bin",
        "lcs_traceback_pointers.txt": "lcs_traceback_txt",
        "global_stats.json": "global_stats",
        "local_stats.json": "local_stats",
    }

    for filename, attr in filename_map.items():
        path = results_dir / filename
        if path.exists():
            setattr(files, attr, path)
            logger.debug("Found %s at %s", attr, path)
        else:
            logger.debug("Missing expected file %s", path)

    return files


def validate_files(alignment_type: str, files: AlignmentFiles) -> None:
    """Validate that required files exist for a specific alignment type.

    Parameters
    ----------
    alignment_type : str
        One of ``"global"``, ``"local"``, or ``"lcs"``.
    files : AlignmentFiles
        The file inventory to validate.

    Raises
    ------
    FileNotFoundError
        If any mandatory file is missing.

    Notes
    -----
    The minimal requirements are:
    - For global: ``global_alignment``, ``global_dp_bin`` or ``global_dp_txt``,
      ``global_path``, and ``global_stats``.
    - For local: similar set with ``local_`` prefix.
    - For LCS: ``lcs_alignment``, ``lcs`` (ungapped), ``lcs_dp_bin`` or
      ``lcs_dp_txt``, ``lcs_path``, and ``lcs_traceback_*``.
    """
    missing: list[str] = []
    if alignment_type == "global":
        if not files.global_alignment:
            missing.append("global_alignment.fasta")
        if not (files.global_dp_bin or files.global_dp_txt):
            missing.append("global_dp_matrix.bin or global_dp_matrix.txt")
        if not files.global_path:
            missing.append("global_path.txt")
        if not files.global_stats:
            missing.append("global_stats.json")
    elif alignment_type == "local":
        if not files.local_alignment:
            missing.append("local_alignment.fasta")
        if not (files.local_dp_bin or files.local_dp_txt):
            missing.append("local_dp_matrix.bin or local_dp_matrix.txt")
        if not files.local_path:
            missing.append("local_path.txt")
        if not files.local_stats:
            missing.append("local_stats.json")
    elif alignment_type == "lcs":
        if not files.lcs_alignment:
            missing.append("lcs_alignment.fasta")
        if not files.lcs:
            missing.append("lcs.fasta")
        if not (files.lcs_dp_bin or files.lcs_dp_txt):
            missing.append("lcs_dp_lengths.bin or lcs_dp_lengths.txt")
        if not files.lcs_path:
            missing.append("lcs_path.txt")
        if not (files.lcs_traceback_bin or files.lcs_traceback_txt):
            missing.append("lcs_traceback_pointers.bin or lcs_traceback_pointers.txt")
    else:
        raise ValueError(f"Unknown alignment type: {alignment_type}")

    if missing:
        raise FileNotFoundError(
            f"Missing required files for {alignment_type} alignment: {', '.join(missing)}"
        )
