# Sequence Alignment Program

## Description
This C++ program performs sequence alignment using three different methods:
1. **Longest Common Subsequence (LCS)**: Finds the length of the longest common subsequence between two sequences.
2. **Global Alignment**: Uses a scoring system with match, mismatch, and gap penalties to align two sequences globally.
3. **Local Alignment**: Identifies the best local alignment between two sequences based on a scoring system.

## Requirements
- C++ compiler (e.g., `g++`)
- Input sequence files (`seq1.txt` and `seq2.txt`)

## Compilation
To compile the program, run:

```sh
g++ -o sequence_alignment sequence_alignment.cpp
```

## Usage
Run the program from the terminal:

```sh
./sequence_alignment
```

Then, enter your choice of alignment method:
- `1` for LCS
- `2` for Global Alignment
- `3` for Local Alignment

## Input Format
The program reads two DNA or protein sequences from text files:
- `seq1.txt` (First sequence)
- `seq2.txt` (Second sequence)

Each file should contain a single line with the sequence.

## Output
- LCS method prints the length of the longest common subsequence.
- Global and Local alignment methods print the alignment score.

## Example
**Input files:**
```
seq1.txt:
ATGCTGA

seq2.txt:
ATGTTGA
```

**Running the program:**
```sh
./sequence_alignment
Select Alignment Method:
1. LCS
2. Global
3. Local
Choice: 2
```

**Output:**
```
Global Alignment Score: 10
```

## License
This program is open-source and free to use for educational purposes.
