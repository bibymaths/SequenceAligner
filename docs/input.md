
## Input Format

* Input files must be in standard **FASTA** format.
* Only the first sequence per file is read.
* FASTA headers (lines starting with '>') are skipped automatically.

Example:

For DNA sequences, use the following format: 

```fasta
>sequence1
AGCTAGCTAGCTA
>sequence2 
GCTAGCTAGCTAG
``` 

For protein sequences, use the following format:

```fasta 
>protein1
MKTAYIAKQRQISFVKSHFSRQDILDL
>protein2
MKTAAYIAKQRQISFVKSHFSRQDILDL
```

---

## Test Files

Located in the `files/` directory: 

* `dna1.fasta`: Contains DNA sequences for testing. 
* `protein1.fasta`: Contains protein sequences for testing.
* `dna2.fasta`: Contains DNA sequences for testing.
* `protein2.fasta`: Contains protein sequences for testing.
