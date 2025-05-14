
# SequenceAligner

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg) ![CMake](https://img.shields.io/badge/CMake-≥3.10-blue.svg)
[![Doxygen](https://img.shields.io/badge/docs-Doxygen-blue)](https://bibymaths.github.io/SequenceAligner/api/index.html)  
 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15414690.svg)](https://doi.org/10.5281/zenodo.15414690)

A simple and efficient C++ implementation of three classic sequence alignment algorithms:  
- **Longest Common Subsequence (LCS)**
- **Global Alignment (Needleman-Wunsch)**
- **Local Alignment (Smith-Waterman)**

The alignments are visualized with colored base-level matches and mismatches for clarity.

Please refer to [API Documentation](https://bibymaths.github.io/SequenceAligner/api/index.html)

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

### Installing OpenMPI for v2 and v3

#### Fedora

To use OpenMPI compilers (`mpicc`, `mpic++`, etc.) and `mpirun` on Fedora:

1. **Install OpenMPI and development headers**:

   ```bash
   sudo dnf install openmpi openmpi-devel
   ```

2. **Enable environment modules**:

   ```bash
   source /etc/profile.d/modules.sh
   ```

3. **Load the OpenMPI module**:

   ```bash
   module load mpi/openmpi-x86_64
   ```

4. **Persistent setup** (optional): Add the above two lines to your `~/.bashrc` to avoid repeating them in each session.

Reference: [OpenMPI on Fedora](https://brandonrozek.com/blog/openmpi-fedora/)

---

#### Debian / Ubuntu

Install OpenMPI and development files with:

```bash
sudo apt install libopenmpi-dev
```

---

#### macOS (Homebrew)

Install using Homebrew:

```bash
brew install openmpi
```

---

For official documentation, see:
[OpenMPI Quickstart Guide](https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html)

---

### Build Instructions
 
Each version implements the same interface. Choose the appropriate binary based on the features or optimizations you want to test. 

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

This project is licensed under the [BSD 3-Clause License](./LICENSE)  
© 2025 Abhinav Mishra.