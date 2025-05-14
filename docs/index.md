# SequenceAligner Documentation

## Overview

SequenceAligner is a lightweight C++ tool designed for the comparative alignment of biological sequences. It supports three classic algorithms:

* **Longest Common Subsequence (LCS)**: Identifies the longest sequence of characters that appear left-to-right (not necessarily contiguously) in both sequences.
* **Global Alignment (Needleman-Wunsch)**: Finds the best full-length alignment between two sequences using dynamic programming.
* **Local Alignment (Smith-Waterman)**: Identifies the highest scoring local subsequence alignment useful for identifying conserved motifs.

This tool is optimized for command-line use with FASTA inputs and outputs colorful visual feedback.

## Command-Line Usage

```bash
./aligner_v1 files/seq1.fasta files/seq2.fasta 2
```

Where:

* First argument is the path to the first FASTA file.
* Second argument is the path to the second FASTA file.
* Third argument is the alignment type:

  * `1`: LCS
  * `2`: Global
  * `3`: Local

---

For complete API documentation, see [API Reference](https://bibymaths.github.io/SequenceAligner/api/index.html).
