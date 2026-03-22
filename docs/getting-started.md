## Installation

```bash
mkdir build && cd build
cmake ..
make
```

## Basic Usage

```bash
./aligner --query q.fasta --target t.fasta --choice 2 --mode dna
```

## Modes

* `dna`: EDNAFULL scoring
* `protein`: BLOSUM62 scoring

!!! Note  
The scoring matrices are embedded directly in the code.