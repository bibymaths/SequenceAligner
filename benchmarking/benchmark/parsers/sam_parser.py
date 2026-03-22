"""
Parser for SAM alignments produced by tools like Bowtie2 and BWA.

SAM (Sequence Alignment/Map) format records detailed information about how
sequencing reads align to a reference. For benchmarking purposes we are
interested in summary statistics such as the percentage of identical bases,
alignment length, number of mismatches and gap events (insertions or
deletions). This module provides a simple parser that operates on the first
mapped alignment in a SAM file and derives these metrics.

The NM optional field in SAM records reports the edit distance between the
read and the reference (i.e., the sum of mismatches and gap openings). We
extract the NM tag and compute mismatches by subtracting the total number of
gap operations in the CIGAR string. Identity is calculated as the ratio of
matches to the total aligned length. Coverage is the aligned length divided
by the query sequence length.
"""

from __future__ import annotations

import re
from typing import Dict, Optional


def _parse_cigar(cigar: str) -> Dict[str, int]:
    """Parse a CIGAR string and count operation lengths.

    Returns a dictionary with operation codes as keys and the sum of lengths
    for each code as integer values.
    """
    op_counts: Dict[str, int] = {}
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        op_counts[op] = op_counts.get(op, 0) + int(length)
    return op_counts


def parse_sam(
    content: str,
    query_lengths: Dict[str, int],
) -> Dict[str, Optional[float]]:
    """Parse SAM content and compute alignment statistics.

    The parser considers the first mapped alignment (i.e., a record whose
    FLAG does not have bit 0x4 set). Header lines (starting with '@') are
    ignored. If no alignments are found, all metrics are set to ``None``.

    Parameters
    ----------
    content: str
        The full SAM file as a single string.
    query_lengths: Dict[str, int]
        Mapping of query sequence identifiers to their lengths. Required for
        coverage calculation.

    Returns
    -------
    Dict[str, Optional[float]]
        Dictionary with keys ``identity``, ``alignment_length``, ``mismatches``,
        ``gap_count``, ``query_coverage``. Values are floats or ``None`` when
        metrics could not be determined.
    """
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith('@'):
            continue
        fields = line.split('\t')
        if len(fields) < 11:
            # Malformed record
            continue
        flag = int(fields[1])
        # bit 0x4 indicates unmapped read
        if flag & 0x4:
            continue
        qname = fields[0]
        cigar = fields[5]
        # Compute query aligned length (excluding soft clipping and hard clipping)
        cigar_ops = _parse_cigar(cigar)
        aligned_query_length = 0
        for op_code, count in cigar_ops.items():
            if op_code in ('M', '=', 'X', 'I'):
                aligned_query_length += count
        # Total gap events (insertions and deletions)
        gap_events = cigar_ops.get('I', 0) + cigar_ops.get('D', 0)
        # Extract NM tag (edit distance)
        nm_value: Optional[int] = None
        for field in fields[11:]:
            if field.startswith('NM:i:'):
                try:
                    nm_value = int(field.split(':', 2)[2])
                except ValueError:
                    nm_value = None
                break
        mismatches: Optional[int] = None
        if nm_value is not None:
            # mismatches = edit distance - gap events
            mismatches = max(nm_value - gap_events, 0)
        # Compute identity and coverage
        identity: Optional[float] = None
        if aligned_query_length > 0 and mismatches is not None:
            matches = aligned_query_length - mismatches - gap_events
            if matches < 0:
                matches = 0
            identity = matches / aligned_query_length * 100.0
        # Coverage
        qlen = query_lengths.get(qname)
        query_cov: Optional[float] = None
        if qlen:
            query_cov = aligned_query_length / qlen * 100.0
        # Return metrics for the first alignment
        return {
            "identity": identity,
            "alignment_length": float(aligned_query_length) if aligned_query_length > 0 else None,
            "mismatches": float(mismatches) if mismatches is not None else None,
            "gap_count": float(gap_events),
            "query_coverage": query_cov,
        }
    # No mapped alignments found
    return {
        "identity": None,
        "alignment_length": None,
        "mismatches": None,
        "gap_count": None,
        "query_coverage": None,
    }
