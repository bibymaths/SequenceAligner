"""
Parsers for alignment outputs.

Submodules in this package define functions to extract useful metrics from
tool-specific output formats. For example, ``blast_parser`` handles tabular
BLAST output (outfmt 6), ``sam_parser`` parses SAM alignments from aligners
such as Bowtie2 and BWA, and ``msa_parser`` computes statistics from aligned
FASTA files produced by multiple sequence aligners like MAFFT or Clustal
Omega.
"""
