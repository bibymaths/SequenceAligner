## Output Format

The output consists of two aligned sequences printed line-by-line:

* **Green**: matched characters
* **Red**: gaps
* **Cyan**: mismatches

Each output block includes the base range for easier visual indexing.

!!! note "FM-Index Anchoring with Proteins"
If you see the message `"FM-index anchoring unavailable/failed. Falling back to MPI full DP"` during protein alignment, this is **expected behavior, not a bug**.

    * **Exact Matches vs. Similarity:** The FM-Index requires exact substring matches (k-mers, typically 5-8 characters) to build anchors. Distant protein sequences (<30% identity) often preserve *chemical similarity* rather than exact identity, meaning they may only share very short exact matches (2-3 amino acids).
    * **Smart Fallback:** Lacking exact matches long enough to safely anchor, the program intentionally skips the FM-Index phase. It gracefully falls back to the full Smith-Waterman or Needleman-Wunsch DP matrix (using BLOSUM62) to ensure a biologically accurate alignment. 
    * **Why not lower the k-mer size?** Forcing a tiny k-mer threshold (like `k=3`) on proteins would result in massive amounts of random, noisy seeds, completely destroying both accuracy and performance. 

---

## Additional Outputs (From Code)

Beyond alignment:

- DP matrices (text or binary)
- Traceback matrices
- LCS sequence output
- Indexed position ranges

Binary formats:

- Efficient storage for large matrices
- Row-major layout with metadata header

--- 

## Future Improvements

* Support multi-sequence alignment.
* Allow multiple FASTA entries.
* Export alignment results in standard formats (CLUSTAL, Stockholm).
* Web-based or GUI interface.

