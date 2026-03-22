"""
Runner for Clustal Omega.

Clustal Omega is a general purpose multiple sequence alignment tool for
proteins. According to its documentation, in its current form Clustal Omega
can only align protein sequences and not DNA or RNA【738991639517845†L64-L66】. This runner
therefore skips execution for nucleotide inputs. It combines the query and
target FASTA into a single file, runs Clustal Omega to produce an aligned
FASTA, and parses the result using ``benchmark.parsers.msa_parser``.

The command-line usage involves specifying an input file with ``-i`` and
an output file with ``-o``. We set ``--outfmt=fasta`` to obtain a FASTA
alignment and ``--force`` to allow overwriting output files【357385058027719†L314-L316】.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, Optional

from .. import utils
from ..parsers import msa_parser


def _get_first_sequence_id(fasta_path: str) -> str:
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
    """Run Clustal Omega on a pair of protein sequences.

    Parameters
    ----------
    query_path, target_path, sequence_type, threads, timeout, work_dir, log_path:
        Same as other runners. Only protein sequences are supported.

    Returns
    -------
    Optional[Dict[str, object]]
        Dictionary with runtime, memory usage and alignment metrics or
        ``None`` if Clustal Omega should not be run.
    """
    logger = logging.getLogger('clustal_runner')
    if sequence_type != 'protein':
        logger.warning("Clustal Omega only supports proteins; skipping run")
        return None
    if not utils.check_executable('clustalo'):
        logger.error("clustalo executable not found; skipping Clustal Omega run")
        return None
    # Read IDs and lengths
    query_id = _get_first_sequence_id(query_path)
    target_id = _get_first_sequence_id(target_path)
    query_lengths = utils.read_fasta_lengths(query_path)
    target_lengths = utils.read_fasta_lengths(target_path)
    query_len = query_lengths.get(query_id, 0)
    target_len = target_lengths.get(target_id, 0)
    # Combine sequences into single FASTA
    combined_path = os.path.join(work_dir, f"clustal_input_{query_id}_{target_id}.fasta")
    with open(query_path, 'r', encoding='utf-8') as qf, open(target_path, 'r', encoding='utf-8') as tf, open(combined_path, 'w', encoding='utf-8') as cf:
        cf.write(qf.read())
        cf.write(tf.read())
    # Output file
    output_path = os.path.join(work_dir, f"clustal_output_{query_id}_{target_id}.fasta")
    # Compose command
    cmd = [
        'clustalo',
        '-i', combined_path,
        '-o', output_path,
        '--outfmt=fasta',
        '--force',
    ]
    if threads and threads > 1:
        cmd += ['--threads', str(threads)]
    # Run command
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
    # Read aligned FASTA
    try:
        with open(output_path, 'r', encoding='utf-8') as outf:
            aligned_content = outf.read()
    except FileNotFoundError:
        aligned_content = ''
    metrics = msa_parser.parse_aligned_fasta(aligned_content, query_id, target_id, query_len, target_len)
    return {
        'tool': 'clustalo',
        'command': ' '.join(cmd),
        'exit_code': exit_code,
        'runtime': runtime,
        'memory': memory,
        'metrics': metrics,
    }
