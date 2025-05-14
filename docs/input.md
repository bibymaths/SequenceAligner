
## Input Format

* Input files must be in standard **FASTA** format.
* Only the first sequence per file is read.
* FASTA headers (lines starting with '>') are skipped automatically.

Example:

```fasta
>sequence1
AGCTAGCTAGCTA
```

---

## Test Files

Located in the `files/` directory:

* `seq1.fasta`: Genomic reference (e.g. RefSeqGene)
* `seq3.fasta`: mRNA variant
* `seq2.fasta`: Additional test case
