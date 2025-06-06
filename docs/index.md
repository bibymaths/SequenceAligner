# SequenceAligner Documentation

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg) ![CMake](https://img.shields.io/badge/CMake-≥3.10-blue.svg)
[![Doxygen](https://img.shields.io/badge/docs-Doxygen-blue)](https://bibymaths.github.io/SequenceAligner/api/index.html)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15414690.svg)](https://doi.org/10.5281/zenodo.15414690)

## Overview

SequenceAligner is a lightweight C++ tool designed for the comparative alignment of biological sequences. It supports three classic algorithms:

* **Longest Common Subsequence (LCS)**: Identifies the longest sequence of characters that appear left-to-right (not necessarily contiguously) in both sequences.
* **Global Alignment (Needleman-Wunsch)**: Finds the best full-length alignment between two sequences using dynamic programming.
* **Local Alignment (Smith-Waterman)**: Identifies the highest scoring local subsequence alignment useful for identifying conserved motifs.

This tool is optimized for command-line use with FASTA inputs and outputs colorful visual feedback.

## Command-Line Usage

### Usage

```bash
./aligner --query <query.fasta> --target <target.fasta> --choice <1|2|3|4> [--mode dna|protein] [--outdir <directory>] [--verbose]
```

### Where:

* `--query` specifies the first input FASTA file.
* `--target` specifies the second input FASTA file.
* `--choice` selects the alignment method:

  * `1` = global
  * `2` = local
  * `3` = LCS
  * `4` = all methods
* `--mode` (optional) sets scoring mode: `dna` (default) or `protein`.
* `--outdir` (optional) sets the output directory (default is current directory).
* `--verbose` (optional) enables progress and detailed output.

---

Please refer to [API Documentation](https://bibymaths.github.io/SequenceAligner/api/index.html)
