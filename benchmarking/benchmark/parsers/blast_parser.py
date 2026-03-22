"""
Parser for BLAST tabular output (outfmt 6).

The BLAST+ suite (blastn/blastp) can output tab-separated values using
``-outfmt 6``, which contains one line per high-scoring pair (HSP). The default
set of columns includes important information such as percentage identity,
alignment length, mismatches, gap openings, e-value and bit score. This
module provides a function to parse such output and compute derived metrics
like coverage on the query and subject sequences.

References
----------
The metagenomics wiki describes the columns of BLAST tabular output format 6.
It lists the twelve default fields: ``qseqid, sseqid, pident, length, mismatch,
gapopen, qstart, qend, sstart, send, evalue, bitscore`` and explains the
meaning of each【633730938334422†L361-L399】.
"""

from __future__ import annotations

from typing import Dict, Optional


def parse_blast_outfmt6(
    content: str,
    query_lengths: Dict[str, int],
    subject_lengths: Optional[Dict[str, int]] = None,
) -> Dict[str, Optional[float]]:
    """Parse BLAST outfmt 6 content and compute summary metrics.

    Parameters
    ----------
    content: str
        Raw output from a BLAST command run with ``-outfmt 6``.
    query_lengths: Dict[str, int]
        Mapping of query sequence identifiers to their lengths.
    subject_lengths: Dict[str, int], optional
        Mapping of subject (target) sequence identifiers to their lengths.

    Returns
    -------
    Dict[str, Optional[float]]
        Dictionary with keys ``identity``, ``alignment_length``, ``mismatches``,
        ``gap_count``, ``bitscore``, ``evalue``, ``query_coverage`` and
        ``subject_coverage``. Values are floats or ``None`` if not available
        (for example subject coverage if subject_lengths is not provided or
        the alignment length cannot be determined).
    """
    lines = [line.strip() for line in content.splitlines() if line.strip()]
    if not lines:
        # No alignments found
        return {
            "identity": None,
            "alignment_length": None,
            "mismatches": None,
            "gap_count": None,
            "bitscore": None,
            "evalue": None,
            "query_coverage": None,
            "subject_coverage": None,
        }
    # Pick the first (best) alignment line
    parts = lines[0].split('\t')
    # Default BLAST outfmt 6 has 12 fields
    # Ensure at least 12 fields present
    if len(parts) < 12:
        raise ValueError("Unexpected BLAST outfmt6 format: fewer than 12 columns")
    qseqid = parts[0]
    sseqid = parts[1]
    pident = float(parts[2])
    alength = int(parts[3])
    mismatch = int(parts[4])
    gapopen = int(parts[5])
    qstart = int(parts[6])
    qend = int(parts[7])
    sstart = int(parts[8])
    send = int(parts[9])
    evalue = float(parts[10])
    bitscore = float(parts[11])
    # Compute coverage as aligned portion of query or subject divided by total length
    qlen = query_lengths.get(qseqid)
    if qlen:
        query_cov = (abs(qend - qstart) + 1) / qlen * 100.0
    else:
        query_cov = None
    subj_cov = None
    if subject_lengths is not None:
        slen = subject_lengths.get(sseqid)
        if slen:
            subj_cov = (abs(send - sstart) + 1) / slen * 100.0
    # Gap count from gap openings; total number of gaps is unknown but openings gives
    # a lower bound. We use the gap openings as gap count here.
    gap_count = gapopen
    return {
        "identity": pident,
        "alignment_length": float(alength),
        "mismatches": float(mismatch),
        "gap_count": float(gap_count),
        "bitscore": bitscore,
        "evalue": evalue,
        "query_coverage": query_cov,
        "subject_coverage": subj_cov,
    }
