
# SequenceAligner

A simple and efficient C++ implementation of three classic sequence alignment algorithms:  
- **Longest Common Subsequence (LCS)**
- **Global Alignment (Needleman-Wunsch)**
- **Local Alignment (Smith-Waterman)**

The alignments are visualized with colored base-level matches and mismatches for clarity.

---

## Project Structure

```

./
├── CMakeLists.txt         
├── files/              
│   ├── seq1.fasta
│   ├── seq2.fasta
│   └── seq3.fasta
├── LICENSE      
├── README.md        
└── src/            
├── main\_v1.cpp   
├── main\_v2.cpp            # (Optional) experimental version
└── main\_v3.cpp            # (Optional) experimental version

````

---

## Getting Started

### Prerequisites

- A C++17-compatible compiler (e.g. `g++`, `clang++`)
- CMake ≥ 3.10

---

Each version implements the same interface. Choose the appropriate binary based on the features or optimizations you want to test.

---

### Build Instructions

```bash
git clone https://github.com/yourusername/SequenceAligner.git
cd SequenceAligner
mkdir build && cd build
cmake ..
make
```

This will compile **three executables**:

* `aligner_v1` → builds from `src/main_v1.cpp`
* `aligner_v2` → builds from `src/main_v2.cpp`
* `aligner_v3` → builds from `src/main_v3.cpp`

---

### Usage

```bash
./aligner_v1 <seq1.fasta> <seq2.fasta> <choice>
```

---

## Input Format

The program accepts standard **FASTA** files. 

---

## History

This project was originally developed in 2014 as part of the *Biological Computation* course during my Bachelor's in Bioinformatics at JUIT, Solan. It was later revisited and optimized for better performance, readability, and maintainability.

---

## License

This project is licensed under the [BSD 3-Clause License](./LICENSE) © 2025 Abhinav Mishra.