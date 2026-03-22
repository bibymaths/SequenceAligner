Here is your content converted cleanly into Markdown, with the table replaced by numbering and structure improved for
readability.

---

# Alignment Analysis Tool Architecture

## Overview

The Python package `alignment_tool` provides a production-quality command-line tool for analysing pairwise protein
alignment outputs. It takes the concrete files produced by global, local and longest-common-subsequence (LCS) alignments
as its only data sources and produces biologically interpretable reports, TSV tables, plots and a JSON summary. The
design is modular, with type-hinted functions, clear doc-strings and logging throughout.

---

## Module Structure

1. **file_inventory.py**
   Scans a results directory and records the presence of expected files. It exposes an `AlignmentFiles` dataclass and a
   `validate_files` function to ensure that the minimum required files for each alignment type are present.

2. **fasta_utils.py**
   Parses aligned and unaligned FASTA files and computes alignment statistics (alignment length, matches, mismatches,
   gaps, percent identity and similarity). It also builds coordinate maps from alignment columns to residue indices and
   detects contiguous conserved blocks.

3. **dp_matrix.py**
   Loads dynamic programming matrices (score matrices for global/local and length matrices for LCS) from binary or text
   files and validates their shape. It uses memory mapping for large binary matrices. The DP matrix is expected to have
   dimensions `(len(seqA)+1, len(seqB)+1)` as per standard alignment algorithms [1].

4. **path_utils.py**
   Parses traceback path files, validates that path coordinates lie within the DP matrix, overlays a path onto a matrix
   and computes simple metrics such as numbers of diagonal/horizontal/vertical steps, gap runs and direction changes.

5. **lcs_utils.py**
   Loads DP length matrices and traceback pointer matrices specific to LCS and reconstructs the LCS path by following
   D/U/L pointers.

6. **block_detection.py**
   A thin wrapper that invokes conserved block detection from `fasta_utils` and returns a `pandas.DataFrame` for
   convenient output.

7. **residue_profiles.py**
   Computes per-residue support profiles for each protein across alignment methods. Each residue is annotated with
   participation flags, partner indices, DP scores, local support, strong block membership and gap proximity.

8. **substitution_analysis.py**
   Provides facilities to load the BLOSUM62 substitution matrix via Biopython, classify amino acids into biologically
   meaningful categories (glycine, proline, cysteine, aromatic, positive and negative residues [2][3]) and summarise
   identical, conservative and radical substitutions across an alignment.

9. **plotting.py**
   Generates plots using matplotlib, including DP score heatmaps (with or without path overlays), residue-support
   profiles, conserved-block comparison tracks and alignment-method participation maps.

10. **comparison.py**
    Assigns a category to each residue based on which methods align it (e.g. `global_only`, `local_only`,
    `global_local_shared`, etc.) and summarises contiguous segments of identical categories.

11. **summary.py**
    Assembles analysis outputs into a single JSON summary capturing input files, sequence IDs and lengths, matrix
    dimensions, method metadata, top conserved blocks and overlap statistics.

12. **cli.py**
    Implements the command-line interface. Subcommands (`global`, `local`, `lcs`, `compare` and `full`) perform the
    corresponding analyses by orchestrating the above modules, writing TSV/JSON/PNG outputs and logging progress.

13. **main.py**
    Provides a thin wrapper so that the package can be executed with `python -m alignment_tool`.

---

## CLI Usage

After installing or unzipping the package, the tool can be run from the command line.

The `--results-dir` option must point to the directory containing the files listed in the problem description.

Optional arguments control:

* output directory
* filename prefix
* block thresholds
* window size
* plot resolution
* substitution matrix

### Examples

```bash
# analyse only the global alignment
data/analyze_alignment.py global --results-dir results/ --outdir analysis_out --prefix run1

# analyse all three methods and perform comparative analysis
python -m alignment_tool.cli full --results-dir results/ --outdir analysis_out --prefix sample --blosum blosum62 --plot-dpi 200

# run comparison only (requires at least two alignments)
python -m alignment_tool.cli compare --results-dir results/ --outdir analysis_out --prefix sample
```

The `full` subcommand automatically runs global, local and LCS analyses (if files exist) and generates comparison
outputs and a summary report.

---

## Outputs

### Per-alignment outputs

For each alignment method, the tool writes:

1. `<prefix>_<method>_alignment_summary.tsv`
   Alignment length, ungapped lengths, matches, mismatches, gaps, percent identity and similarity.

2. `<prefix>_<method>_conserved_blocks.tsv`
   Statistics for contiguous gap-free blocks (length, identity, similarity, classification, residue ranges).

3. `<prefix>_<method>_path_metrics.tsv`
   Path structure metrics: diagonal/horizontal/vertical steps, gap runs, direction changes.

4. `<prefix>_<method>_residue_support_<sequenceID>.tsv`
   Per-residue support metrics (participation, partner index, DP score, local support, block membership, gap proximity).

5. `<prefix>_<method>_substitution_summary.tsv`
   Counts of identical, conservative, radical substitutions and residue category counts.

6. PNG plots

    * DP heatmaps (with/without path)
    * Residue support profiles

---

### Comparison outputs

1. `<prefix>_alignment_method_comparison_categories_<sequenceID>.tsv`
   Per-residue category assignments (e.g. `global_only`, `global_local_shared`).

2. `<prefix>_alignment_method_comparison_<sequenceID>.tsv`
   Contiguous category segments with start/end indices and lengths.

3. `<prefix>_alignment_method_comparison_<sequenceID>.png`
   Visual comparison of alignment coverage across methods.

---

### Full pipeline output

* `<prefix>_summary.json`
  Consolidated metadata:

    * input files
    * sequence info
    * DP matrix shapes
    * scoring parameters
    * top conserved blocks
    * identity/similarity stats
    * overlap across methods

---

## Biological Interpretation and Limitations

The tool focuses on biologically interpretable signals:

* Gap-free high-identity regions → potential conserved domains
* Conservative substitutions → functional similarity signals
* Residue support profiles → robustness of alignment at residue level
* Special residues tracked:

    * glycine
    * proline
    * cysteine
    * aromatic
    * charged residues [2][3]

DP heatmaps with traceback paths reveal:

* diagonal corridors → conserved regions
* long gaps → insertions/deletions

### Limitations

* No structural inference
* No functional claims
* LCS captures **only exact matches** → misses conservative substitutions
* Global alignment may force weak regions
* Local alignment isolates strongest core

This tool reports these differences—it does not interpret beyond evidence.

---

## Expected Output Directory Structure

Example (`--outdir analysis_out`, `--prefix sample`):

```
analysis_out/
├── sample_global_alignment_summary.tsv
├── sample_global_conserved_blocks.tsv
├── sample_global_path_metrics.tsv
├── sample_global_residue_support_P14416.tsv
├── sample_global_residue_support_P35462.tsv
├── sample_global_substitution_summary.tsv
├── sample_global_dp_heatmap.png
├── sample_global_dp_heatmap_with_path.png
├── sample_global_residue_support_P14416.png
├── sample_global_residue_support_P35462.png
├── sample_local_alignment_summary.tsv
├── ...
├── sample_lcs_alignment_summary.tsv
├── ...
├── sample_alignment_method_comparison_categories_P14416.tsv
├── sample_alignment_method_comparison_P14416.tsv
├── sample_alignment_method_comparison_P14416.png
├── sample_alignment_method_comparison_categories_P35462.tsv
├── sample_alignment_method_comparison_P35462.tsv
├── sample_alignment_method_comparison_P35462.png
└── sample_summary.json
```

The exact outputs depend on which alignment files are present in `results/`.

---

## References

[1] Chapter 3: Sequence Alignments – Applied Bioinformatics
[https://open.oregonstate.education/appliedbioinformatics/chapter/chapter-3/](https://open.oregonstate.education/appliedbioinformatics/chapter/chapter-3/)

[2][3] IMGT Education
[https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html](https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html)

---