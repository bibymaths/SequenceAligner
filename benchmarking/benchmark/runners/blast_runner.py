"""
Runner for BLAST (blastn and blastp).

This module implements a function to execute NCBI BLAST+ for pairwise
alignments. For nucleotide sequences the ``blastn`` program is used, and
for proteins ``blastp``. The runner invokes BLAST in subject mode to
compare two FASTA files directly without building a database. Output is
generated in tabular format (outfmt 6) and parsed using
``benchmark.parsers.blast_parser``.

If BLAST is not installed or unavailable in the execution environment, the
runner returns ``None`` and logs the condition. It also handles unsupported
sequence types by returning ``None``.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, Optional

from .. import utils
from ..parsers import blast_parser


def run(
    query_path: str,
    target_path: str,
    sequence_type: str,
    threads: int,
    timeout: Optional[int],
    work_dir: str,
    log_path: str,
) -> Optional[Dict[str, object]]:
    """Run BLAST on the given query and target sequences.

    Parameters
    ----------
    query_path: str
        Path to the FASTA file containing the query sequence(s).
    target_path: str
        Path to the FASTA file containing the subject (target) sequence(s).
    sequence_type: str
        Either 'dna' or 'protein'. Determines which BLAST program to use.
    threads: int
        Number of threads to request from BLAST via ``-num_threads``.
    timeout: int or None
        Maximum number of seconds to allow the BLAST search to run.
    work_dir: str
        Directory where temporary files may be written.
    log_path: str
        File path where the raw stdout/stderr will be logged.

    Returns
    -------
    Optional[Dict[str, object]]
        A dictionary containing runtime, memory usage and alignment metrics. If
        BLAST is not installed or the sequence type is unsupported, returns
        ``None``.
    """
    logger = logging.getLogger('blast_runner')
    # Determine which program to use
    if sequence_type == 'dna':
        prog = 'blastn'
    elif sequence_type == 'protein':
        prog = 'blastp'
    else:
        logger.warning(f"BLAST does not support sequence type: {sequence_type}")
        return None
    if not utils.check_executable(prog):
        logger.error(f"{prog} not found in PATH; skipping BLAST run")
        return None
    # Compose command. Use subject mode to avoid database creation.
    cmd = [prog, '-query', query_path, '-subject', target_path, '-outfmt', '6']
    # Use threads if supported
    if threads and threads > 1:
        cmd += ['-num_threads', str(threads)]
    logger.debug(f"Running BLAST command: {' '.join(cmd)}")
    # Execute command
    elapsed, peak_mem, exit_code, stdout, stderr = utils.run_subprocess_with_resource_tracking(
        cmd, timeout=timeout, capture_output=True
    )
    # Write raw logs
    with open(log_path, 'w', encoding='utf-8') as logf:
        logf.write('Command: ' + ' '.join(cmd) + '\n')
        logf.write('Exit code: ' + str(exit_code) + '\n')
        logf.write('=== STDOUT ===\n')
        logf.write(stdout)
        logf.write('\n=== STDERR ===\n')
        logf.write(stderr)
    # Parse metrics
    # Precompute sequence lengths for coverage
    query_lengths = utils.read_fasta_lengths(query_path)
    subject_lengths = utils.read_fasta_lengths(target_path)
    metrics = blast_parser.parse_blast_outfmt6(stdout, query_lengths, subject_lengths)
    return {
        'tool': prog,
        'command': ' '.join(cmd),
        'exit_code': exit_code,
        'runtime': elapsed,
        'memory': peak_mem,
        'metrics': metrics,
    }
