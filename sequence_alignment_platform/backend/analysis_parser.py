"""
Utilities for discovering and parsing downstream analysis outputs.

The analysis tools produce a variety of file types: TSV/CSV tables,
PNG plots, and summary JSON. This module scans an analysis
directory, groups outputs by categories (e.g. global, local, lcs,
method‑comparison), and provides functions to parse tabular data into
Python structures suitable for JSON responses.
"""

import csv
from pathlib import Path
from typing import Dict, List


def discover_analysis_outputs(analysis_dir: Path) -> Dict[str, Dict[str, List[str]]]:
    """Group analysis files by category and return a nested dict."""
    grouped: Dict[str, Dict[str, List[str]]] = {}
    for path in analysis_dir.iterdir():
        if not path.is_file():
            continue
        parts = path.name.split("_")
        if len(parts) < 2:
            key = "misc"
        else:
            key = parts[1]  # e.g. global, local, lcs, etc.
        ext_group = grouped.setdefault(key, {})
        ext = path.suffix.lstrip(".")
        ext_group.setdefault(ext, []).append(path.name)
    return grouped


def parse_tsv(file_path: Path) -> List[Dict[str, str]]:
    """Parse a TSV or CSV file and return a list of dict records."""
    delimiter = "\t" if file_path.suffix == ".tsv" else ","
    with open(file_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        return [dict(row) for row in reader]