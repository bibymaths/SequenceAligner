# Benchmarks

This section compares **SequenceAligner** against widely used sequence alignment tools.

Benchmarks are designed to evaluate:

- Runtime performance
- Memory usage
- Scalability
- Alignment accuracy (where applicable)

---

## Tools Compared

The following tools will be included:

- BLAST
- Bowtie2
- BWA
- MAFFT (for alignment comparison)
- Clustal Omega

These tools represent different classes:

| Tool          | Category                    |
|---------------|-----------------------------|
| BLAST         | Heuristic local alignment   |
| Bowtie2       | FM-index based aligner      |
| BWA           | Burrows-Wheeler aligner     |
| MAFFT         | Multiple sequence alignment |
| Clustal Omega | Progressive MSA             |

---

## Benchmark Design

### Datasets

Benchmarks will be performed on:

- Synthetic DNA sequences (controlled length and mutation rate)
- Real biological datasets (e.g. UniProt, GenBank)
- Long sequences (>100k bp)
- Protein sequences (BLOSUM62 mode)

---

### Metrics

#### 1. Runtime

Measured as:

- Total execution time
- Time per base pair

#### 2. Memory Usage

- Peak RAM usage
- Memory scaling with sequence length

#### 3. Accuracy

- Alignment score
- Identity percentage
- Gap distribution

---

## Methodology

Each tool is run under identical conditions:

- Same input sequences
- Same hardware
- Same number of threads (where applicable)

SequenceAligner configurations:

- FM-index enabled
- Affine gap penalties
- SIMD optimizations active
- MPI disabled (single-node baseline unless stated)

---

## Planned Experiments

### 1. Short Sequence Benchmark

- Length: 100–1,000 bp
- Focus: overhead vs classic tools

### 2. Medium Sequence Benchmark

- Length: 1k–50k bp
- Focus: DP vs indexed methods

### 3. Long Sequence Benchmark

- Length: 50k–1M bp
- Focus: scalability and memory

### 4. Protein Alignment Benchmark

- BLOSUM62 scoring
- Real protein datasets

---

## Expected Observations (Hypothesis)

Based on implementation:

- Faster than pure DP methods due to seed filtering
- Competitive with Bowtie2/BWA for exact matches
- Higher accuracy than heuristic tools in local alignment
- Better performance on long sequences due to:
    - FM-index seeding
    - Reduced DP search space

---

## Results

*(To be populated)*

### Runtime Comparison

| Tool            | Dataset | Time (s) |
|-----------------|---------|----------|
| SequenceAligner | TBD     | TBD      |
| BLAST           | TBD     | TBD      |
| Bowtie2         | TBD     | TBD      |

---

### Memory Usage

| Tool            | Peak RAM |
|-----------------|----------|
| SequenceAligner | TBD      |
| BLAST           | TBD      |

---

### Accuracy

| Tool            | Identity (%) | Score |
|-----------------|--------------|-------|
| SequenceAligner | TBD          | TBD   |

---

## Reproducibility

All benchmark scripts and datasets will be provided.

---

## Notes

SequenceAligner is not a drop-in replacement for all tools:

- It emphasizes **exact and controlled alignment**
- Not optimized for database-scale search (like BLAST)
- Not a full MSA tool (like MAFFT)

Instead, it targets:

- High-precision pairwise alignment
- Algorithmic transparency
- Research and development workflows