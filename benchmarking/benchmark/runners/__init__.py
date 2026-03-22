"""
Runners for each alignment tool.

This package contains one module per supported alignment tool. Each module
exposes a `run` function that executes the tool on given input files, measures
runtime and peak memory consumption, captures stdout/stderr, and returns
alignment statistics using parsers defined in ``benchmark/parsers``. If a
tool is not installed or does not support the requested sequence type (DNA or
protein), the runner returns ``None`` and logs an appropriate message.
"""

from .blast_runner import run as run_blast  # noqa: F401
from .bowtie2_runner import run as run_bowtie2  # noqa: F401
from .bwa_runner import run as run_bwa  # noqa: F401
from .mafft_runner import run as run_mafft  # noqa: F401
from .clustal_runner import run as run_clustal  # noqa: F401

__all__ = [
    "run_blast",
    "run_bowtie2",
    "run_bwa",
    "run_mafft",
    "run_clustal",
]
