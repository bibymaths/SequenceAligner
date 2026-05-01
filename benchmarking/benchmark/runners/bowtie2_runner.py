"""
Runner for Bowtie2.

Bowtie2 is a fast DNA read aligner that requires building an index from the
reference genome. This runner builds an index on-the-fly for the provided
target FASTA if no cached index is found, then aligns the query FASTA using
the Bowtie2 ``-f -U`` options for unpaired FASTA reads. Bowtie2 outputs
alignments in SAM format, which are parsed with ``benchmark.parsers.sam_parser``.

Bowtie2 does not support protein sequences, so this runner returns ``None``
when invoked for protein benchmarking.

Reference: The Bowtie2 manual describes how to build an index using
``bowtie2-build <reference_in> <bt2_base>`` and then align reads with
``bowtie2 -x <bt2_base> -U <reads.fq> -S <output.sam>``【354970811977329†L2124-L2160】.
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
    """Run Bowtie2 on the given DNA sequences.

    Parameters
    ----------
    query_path: str
        Path to FASTA file containing the query sequence(s).
    target_path: str
        Path to FASTA file containing the target (reference) sequence(s).
    sequence_type: str
        Expected to be 'dna'. Bowtie2 does not support proteins.
    threads: int
        Number of threads for Bowtie2 to use for both index building and
        alignment.
    timeout: int or None
        Maximum runtime in seconds for each of the Bowtie2 commands.
    work_dir: str
        Directory where temporary files (index and SAM outputs) will be stored.
    log_path: str
        Path to write the combined Bowtie2 build and alignment logs.

    Returns
    -------
    Optional[Dict[str, object]]
        Dictionary containing runtime, memory usage and alignment metrics,
        or ``None`` if Bowtie2 is unavailable or the sequence type is not DNA.
    """
    logger = logging.getLogger("bowtie2_runner")
    if sequence_type != "dna":
        logger.warning("Bowtie2 only supports DNA alignment; skipping run")
        return None
    if not (
        utils.check_executable("bowtie2") and utils.check_executable("bowtie2-build")
    ):
        logger.error("bowtie2 or bowtie2-build not found in PATH; skipping Bowtie2 run")
        return None
    # Prepare index directory and prefix
    base_name = os.path.basename(target_path)
    index_prefix = os.path.join(work_dir, f"bowtie2_index_{base_name}")
    # Determine if index already exists
    index_files = [
        f"{index_prefix}.{ext}.bt2" for ext in ["1", "2", "3", "4", "rev.1", "rev.2"]
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
        build_cmd = ["bowtie2-build", target_path, index_prefix]
        if threads and threads > 1:
            build_cmd += ["--threads", str(threads)]
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
    # Align reads
    # Bowtie2 expects reads; we treat query FASTA as unpaired reads with -f -U
    sam_file = os.path.join(work_dir, f"bowtie2_{base_name}.sam")
    align_cmd = [
        "bowtie2",
        "-x",
        index_prefix,
        "-f",
        "-U",
        query_path,
        "-S",
        sam_file,
    ]
    if threads and threads > 1:
        align_cmd += ["--threads", str(threads)]
    align_runtime, align_memory, align_exit, align_stdout, align_stderr = (
        utils.run_subprocess_with_resource_tracking(
            align_cmd, timeout=timeout, capture_output=True
        )
    )
    # Append alignment logs
    with open(log_path, "a", encoding="utf-8") as logf:
        logf.write("Command: " + " ".join(align_cmd) + "\n")
        logf.write("Exit code: " + str(align_exit) + "\n")
        logf.write("=== STDOUT (alignment) ===\n")
        logf.write(align_stdout + "\n")
        logf.write("=== STDERR (alignment) ===\n")
        logf.write(align_stderr + "\n")
    # Read SAM content
    try:
        with open(sam_file, "r", encoding="utf-8") as sf:
            sam_content = sf.read()
    except FileNotFoundError:
        sam_content = ""
    # Parse metrics
    query_lengths = utils.read_fasta_lengths(query_path)
    metrics = sam_parser.parse_sam(sam_content, query_lengths)
    # Combine runtime and memory: sum build and align runtime, take max memory
    total_runtime = build_runtime + align_runtime
    peak_memory = None
    if build_memory is not None and align_memory is not None:
        peak_memory = max(build_memory, align_memory)
    elif build_memory is not None:
        peak_memory = build_memory
    else:
        peak_memory = align_memory
    return {
        "tool": "bowtie2",
        "command": " | ".join(
            [" ".join(build_cmd) if build_cmd else "", " ".join(align_cmd)]
        ).strip(" |"),
        "exit_code": align_exit,
        "runtime": total_runtime,
        "memory": peak_memory,
        "metrics": metrics,
    }
