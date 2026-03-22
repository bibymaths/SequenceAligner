"""Top‑level package for the protein alignment analysis tool.

This package provides utilities to parse and analyse pairwise protein
alignment outputs derived from dynamic programming alignment algorithms.
The main entrypoint is via the CLI defined in ``alignment_tool.cli``.

The version is stored in ``__version__`` for programmatic access.
"""

__all__ = [
    "__version__",
    "cli",
    "file_inventory",
    "fasta_utils",
    "dp_matrix",
    "path_utils",
    "lcs_utils",
    "block_detection",
    "residue_profiles",
    "substitution_analysis",
    "plotting",
    "comparison",
    "summary",
]

#: Current version of the software
__version__ = "0.1.0"
