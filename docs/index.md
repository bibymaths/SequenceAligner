<img src="assets/logo.png" alt="SequenceAligner" width="300">

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg) ![CMake](https://img.shields.io/badge/CMake-≥3.10-blue.svg)
[![Doxygen](https://img.shields.io/badge/docs-Doxygen-blue)](https://bibymaths.github.io/SequenceAligner/api/index.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15414690.svg)](https://doi.org/10.5281/zenodo.15414690)

SequenceAligner is a high-performance C++ tool for biological sequence alignment. It combines **classical dynamic
programming** with **modern indexing techniques (FM-index)** and **parallel computation (MPI + OpenMP)**.

It supports:

* Longest Common Subsequence (LCS)
* Global alignment (Needleman–Wunsch with affine gaps)
* Local alignment (Smith–Waterman with affine gaps)
* Seed-and-extend alignment using FM-index

The tool is optimized for **large-scale sequence comparison**, not just textbook examples.

---

## What makes this different

This is not a naive DP implementation.

From the code:

* FM-index enables fast substring search and seed generation
* Suffix arrays are constructed in (O(n \log n))
* Alignment uses affine gap penalties (GAP_OPEN, GAP_EXTEND)
* SIMD (`immintrin.h`) and MPI support large-scale execution
* Optional binary output for DP matrices

---

## Workflow

1. Parse FASTA input
2. Build FM-index on target
3. Generate k-mer seeds
4. Chain seeds into candidate regions
5. Run DP alignment (global/local)
6. Output formatted alignment + optional matrices

---

Please refer to [API Documentation](https://bibymaths.github.io/SequenceAligner/api/index.html)
