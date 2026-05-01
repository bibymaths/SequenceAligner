"""
Parser for pairwise alignments produced by multiple sequence aligners.

MAFFT and Clustal Omega produce multiple sequence alignments in FASTA format
where all sequences are padded with gap characters ('-') to the same length.
For pairwise benchmarking we consider just two sequences: the query and the
target. This parser computes basic statistics such as the number of matches,
mismatches, gaps, percentage identity and coverage.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple


def _read_fasta_sequences(content: str) -> List[Tuple[str, str]]:
    """Parse FASTA content into a list of (identifier, sequence) tuples.

    Parameters
    ----------
    content: str
        Raw FASTA content.

    Returns
    -------
    List[Tuple[str, str]]
        List of (id, sequence) tuples in the order they appear in the file.
    """
    seqs: List[Tuple[str, str]] = []
    current_id: Optional[str] = None
    current_seq: List[str] = []
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                seqs.append((current_id, "".join(current_seq)))
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id is not None:
        seqs.append((current_id, "".join(current_seq)))
    return seqs


def parse_aligned_fasta(
    content: str,
    query_id: str,
    target_id: str,
    query_length: int,
    target_length: int,
) -> Dict[str, Optional[float]]:
    """Compute alignment statistics from aligned FASTA content.

    Parameters
    ----------
    content: str
        FASTA formatted alignment containing two sequences.
    query_id: str
        Identifier of the query sequence in the original input. The parser uses
        this to select the appropriate aligned sequence.
    target_id: str
        Identifier of the target sequence in the original input.
    query_length: int
        Length of the unaligned query sequence.
    target_length: int
        Length of the unaligned target sequence.

    Returns
    -------
    Dict[str, Optional[float]]
        Dictionary with keys ``identity``, ``alignment_length``, ``mismatches``,
        ``gap_count``, ``query_coverage``, ``target_coverage``. Values are floats
        or ``None`` when they cannot be determined.
    """
    sequences = _read_fasta_sequences(content)
    if len(sequences) < 2:
        return {
            "identity": None,
            "alignment_length": None,
            "mismatches": None,
            "gap_count": None,
            "query_coverage": None,
            "target_coverage": None,
        }
    # Find aligned query and target by ID
    aligned_query = None
    aligned_target = None
    for sid, seq in sequences:
        if sid == query_id:
            aligned_query = seq
        elif sid == target_id:
            aligned_target = seq
    # Fallback: if exact IDs not found, use first two sequences
    if aligned_query is None or aligned_target is None:
        aligned_query = sequences[0][1]
        aligned_target = sequences[1][1]
    # Ensure sequences are same length
    if len(aligned_query) != len(aligned_target):
        min_len = min(len(aligned_query), len(aligned_target))
        aligned_query = aligned_query[:min_len]
        aligned_target = aligned_target[:min_len]
    matches = 0
    mismatches = 0
    gaps = 0
    # Count non-gap positions for coverage
    query_aligned_residues = 0
    target_aligned_residues = 0
    for a, b in zip(aligned_query, aligned_target):
        if a != "-":
            query_aligned_residues += 1
        if b != "-":
            target_aligned_residues += 1
        if a == "-" or b == "-":
            gaps += 1
        elif a.upper() == b.upper():
            matches += 1
        else:
            mismatches += 1
    alignment_length = len(aligned_query)
    identity = (matches / alignment_length * 100.0) if alignment_length > 0 else None
    query_cov = (
        (query_aligned_residues / query_length * 100.0) if query_length > 0 else None
    )
    target_cov = (
        (target_aligned_residues / target_length * 100.0) if target_length > 0 else None
    )
    return {
        "identity": identity,
        "alignment_length": float(alignment_length),
        "mismatches": float(mismatches),
        "gap_count": float(gaps),
        "query_coverage": query_cov,
        "target_coverage": target_cov,
    }
