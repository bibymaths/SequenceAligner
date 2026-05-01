"""
Runner for BWA (Burrows–Wheeler Aligner).

BWA is a DNA aligner capable of mapping reads to a reference genome. Before
alignment the reference must be indexed using ``bwa index``. This runner
creates an index for the provided target FASTA (unless an existing index is
found), then runs ``bwa mem`` to align the query FASTA. The resulting SAM
output is parsed via ``benchmark.parsers.sam_parser``. BWA does not support
protein sequences, so this runner exits early for protein inputs.

The BWA manual outlines basic usage: create an index with ``bwa index
ref.fa`` and align sequences using ``bwa mem ref.fa reads.fq > aln.sam``【28747307351112†L32-L37】.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, Optional

from .. import utils
from ..parsers import sam_parser


def run(
    query_path: str,
    target_path: str,
    sequence_type: str,
    threads: int,
    timeout: Optional[int],
    work_dir: str,
    log_path: str,
) -> Optional[Dict[str, object]]:
    """Run BWA on the given DNA sequences.

    Parameters
    ----------
    query_path: str
        Path to FASTA file with query sequences.
    target_path: str
        Path to FASTA file with target sequences.
    sequence_type: str
        Should be 'dna'; BWA is not designed for proteins.
    threads: int
        Number of threads for BWA operations.
    timeout: int or None
        Maximum allowed runtime for each command.
    work_dir: str
        Directory for temporary index and SAM output.
    log_path: str
        Path to write combined logs of indexing and alignment.

    Returns
    -------
    Optional[Dict[str, object]]
        Dictionary with runtime, memory usage and alignment metrics or
        ``None`` if BWA cannot run.
    """
    logger = logging.getLogger("bwa_runner")
    if sequence_type != "dna":
        logger.warning("BWA only supports DNA alignment; skipping run")
        return None
    if not utils.check_executable("bwa"):
        logger.error("bwa executable not found in PATH; skipping BWA run")
        return None
    # Prepare index prefix
    base_name = os.path.basename(target_path)
    index_prefix = os.path.join(work_dir, f"bwa_index_{base_name}")
    # BWA index produces files with extensions .amb, .ann, .bwt, .pac, .sa
    index_files = [
        f"{index_prefix}.{ext}" for ext in ["amb", "ann", "bwt", "pac", "sa"]
    ]
    index_exists = all(os.path.exists(f) for f in index_files)
    build_cmd = None
    build_runtime = 0.0
    build_memory = None
    build_exit = 0
    build_stdout = ""
    build_stderr = ""
    if not index_exists:
        # Build index
        build_cmd = ["bwa", "index", "-p", index_prefix, target_path]
        build_runtime, build_memory, build_exit, build_stdout, build_stderr = (
            utils.run_subprocess_with_resource_tracking(
                build_cmd, timeout=timeout, capture_output=True
            )
        )
        # Write build log
        with open(log_path, "w", encoding="utf-8") as logf:
            logf.write("Command: " + " ".join(build_cmd) + "\n")
            logf.write("Exit code: " + str(build_exit) + "\n")
            logf.write("=== STDOUT (index build) ===\n")
            logf.write(build_stdout + "\n")
            logf.write("=== STDERR (index build) ===\n")
            logf.write(build_stderr + "\n")
    # Align
    sam_file = os.path.join(work_dir, f"bwa_{base_name}.sam")
    align_cmd = ["bwa", "mem"]
    if threads and threads > 1:
        align_cmd += ["-t", str(threads)]
    align_cmd += [index_prefix, query_path]
    # bwa mem writes SAM to stdout; we will capture via capture_output True
    align_runtime, align_memory, align_exit, align_stdout, align_stderr = (
        utils.run_subprocess_with_resource_tracking(
            align_cmd, timeout=timeout, capture_output=True
        )
    )
    # Write alignment logs
    with open(log_path, "a", encoding="utf-8") as logf:
        logf.write("Command: " + " ".join(align_cmd) + "\n")
        logf.write("Exit code: " + str(align_exit) + "\n")
        logf.write("=== STDOUT (alignment) ===\n")
        logf.write(align_stdout + "\n")
        logf.write("=== STDERR (alignment) ===\n")
        logf.write(align_stderr + "\n")
    # Write SAM file for parsing
    with open(sam_file, "w", encoding="utf-8") as sf:
        sf.write(align_stdout)
    # Parse metrics
    query_lengths = utils.read_fasta_lengths(query_path)
    metrics = sam_parser.parse_sam(align_stdout, query_lengths)
    # Summarise runtime and memory
    total_runtime = build_runtime + align_runtime
    peak_memory = None
    if build_memory is not None and align_memory is not None:
        peak_memory = max(build_memory, align_memory)
    elif build_memory is not None:
        peak_memory = build_memory
    else:
        peak_memory = align_memory
    return {
        "tool": "bwa",
        "command": " | ".join(
            [" ".join(build_cmd) if build_cmd else "", " ".join(align_cmd)]
        ).strip(" |"),
        "exit_code": align_exit,
        "runtime": total_runtime,
        "memory": peak_memory,
        "metrics": metrics,
    }
