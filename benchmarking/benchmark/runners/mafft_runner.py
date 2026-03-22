"""
Runner for MAFFT.

MAFFT is a general-purpose multiple sequence aligner for amino acid or nucleotide
sequences. For benchmarking pairwise alignments we use MAFFT on a two-sequence
FASTA file, letting MAFFT choose an appropriate algorithm via ``--auto``. The
aligned FASTA output is parsed to compute identity, mismatches, gaps and
coverage. MAFFT can align both proteins and nucleotides; however, in this
benchmarking framework we restrict its use to protein sequences to ensure
methodological consistency.

Usage is simple: ``mafft [arguments] input > output``【395817269289108†L36-L59】. The type of sequences
(amino acid or nucleotide) is automatically recognised【395817269289108†L29-L36】.
"""

from __future__ import annotations

import logging
import os
import tempfile
from typing import Dict, Optional

from .. import utils
from ..parsers import msa_parser


def _get_first_sequence_id(fasta_path: str) -> str:
    """Return the ID of the first sequence in a FASTA file."""
    with open(fasta_path, 'r', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                return line[1:].split()[0]
    raise ValueError(f"No sequence found in {fasta_path}")


def run(
    query_path: str,
    target_path: str,
    sequence_type: str,
    threads: int,
    timeout: Optional[int],
    work_dir: str,
    log_path: str,
) -> Optional[Dict[str, object]]:
    """Run MAFFT on a pair of sequences.

    Parameters
    ----------
    query_path: str
        Path to FASTA file containing the query sequence.
    target_path: str
        Path to FASTA file containing the target sequence.
    sequence_type: str
        Either 'protein' or 'dna'. For benchmarking we only execute MAFFT
        on proteins to avoid unfair comparisons.
    threads: int
        Number of threads for MAFFT. MAFFT uses ``--thread N``.
    timeout: int or None
        Maximum runtime in seconds.
    work_dir: str
        Directory to place temporary input and output files.
    log_path: str
        File path to write logs.

    Returns
    -------
    Optional[Dict[str, object]]
        Dictionary containing runtime, memory and alignment metrics, or
        ``None`` if MAFFT should not run or is unavailable.
    """
    logger = logging.getLogger('mafft_runner')
    if sequence_type != 'protein':
        # MAFFT can align DNA, but we limit its use to proteins for fair comparison
        logger.warning("Skipping MAFFT for non-protein sequences")
        return None
    if not utils.check_executable('mafft'):
        logger.error("mafft not found in PATH; skipping MAFFT run")
        return None
    # Read IDs and lengths
    query_id = _get_first_sequence_id(query_path)
    target_id = _get_first_sequence_id(target_path)
    query_lengths = utils.read_fasta_lengths(query_path)
    target_lengths = utils.read_fasta_lengths(target_path)
    query_len = query_lengths.get(query_id, 0)
    target_len = target_lengths.get(target_id, 0)
    # Create combined FASTA
    combined_path = os.path.join(work_dir, f"mafft_input_{query_id}_{target_id}.fasta")
    with open(query_path, 'r', encoding='utf-8') as qf, open(target_path, 'r', encoding='utf-8') as tf, open(combined_path, 'w', encoding='utf-8') as cf:
        cf.write(qf.read())
        if not qf.seekable():
            # Should not happen; just append newline
            cf.write('\n')
        cf.write(tf.read())
    # Prepare MAFFT command
    cmd = ['mafft', '--auto']
    if threads and threads > 1:
        cmd += ['--thread', str(threads)]
    cmd.append(combined_path)
    # Run MAFFT
    runtime, memory, exit_code, stdout, stderr = utils.run_subprocess_with_resource_tracking(
        cmd, timeout=timeout, capture_output=True
    )
    # Write logs
    with open(log_path, 'w', encoding='utf-8') as logf:
        logf.write('Command: ' + ' '.join(cmd) + '\n')
        logf.write('Exit code: ' + str(exit_code) + '\n')
        logf.write('=== STDOUT ===\n')
        logf.write(stdout + '\n')
        logf.write('=== STDERR ===\n')
        logf.write(stderr + '\n')
    # Parse alignment
    metrics = msa_parser.parse_aligned_fasta(stdout, query_id, target_id, query_len, target_len)
    return {
        'tool': 'mafft',
        'command': ' '.join(cmd),
        'exit_code': exit_code,
        'runtime': runtime,
        'memory': memory,
        'metrics': metrics,
    }
