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

## Getting Started

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

*NOTE*: **Sometimes, you need to run the above command before cmake and make commands.**

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
brew install llvm
```

---

For official documentation, see:
[OpenMPI Quickstart Guide](https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html)

---

### Build Instructions

Each version implements the same interface. Choose the appropriate binary based on the features or optimizations you
want to test.

For Debian/Ubuntu, Fedora, and the following commands will build the project:

```bash
git clone https://github.com/yourusername/SequenceAligner.git
cd SequenceAligner
mkdir build && cd build
cmake ..
make
cd ..
```

For macOS (Intel Chip), you can use the following commands:

```bash 
git clone https://github.com/yourusername/SequenceAligner.git
cd SequenceAligner
mkdir build && cd build
cmake .. \
  -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang \
  -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ \
  -DCMAKE_C_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/llvm/include" \
  -DCMAKE_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/llvm/include" \
  -DCMAKE_EXE_LINKER_FLAGS="-L/usr/local/opt/llvm/lib -lomp"
make
cd ..
``` 

For macOS (Apple Silicon), you can use the following commands:

```bash 
git clone https://github.com/yourusername/SequenceAligner.git
cd SequenceAligner
mkdir build && cd build
cmake .. \
  -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ \
  -DCMAKE_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/llvm/include" \
  -DCMAKE_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/llvm/include" \
  -DCMAKE_EXE_LINKER_FLAGS="-L/opt/homebrew/opt/llvm/lib -lomp"
make
cd ..
```

This will compile an executable named `aligner` in the project's directory:

* `aligner` → builds from `src/main.cpp`

---

### Usage

To run on cluster or server with OpenMPI, use the following command:

```bash
mpirun -np <num_processes> ./aligner \
  --query <query.fasta> \
  --target <target.fasta> \
  --choice <1|2|3|4> \
  --mode <dna|protein> \
  [--outdir <output_directory>] \ 
  [--binary <binary file DP>] \ 
  [--txt <text file DP>] \
  [--gap_open <float>] \
  [--gap_extend <float>] \
  [--verbose]
``` 

**Note**: Use `mpirun` when you want performance via parallelism—especially for long sequences or many alignments.  
It can run on multiple CPU cores or even multiple nodes if configured

To run on local machine, use the following command:

```bash 
./aligner \
  --query <query.fasta> \
  --target <target.fasta> \
  --choice <1|2|3|4> \
  --mode <dna|protein> \
  [--outdir <output_directory>] \ 
  [--binary <binary file DP>] \ 
  [--txt <text file DP>] \
  [--gap_open <float>] \
  [--gap_extend <float>] \
  [--verbose]
```

---

## 🚀 Running with CLion GUI

To run the aligner within CLion, you must configure the **Run/Debug Configurations** to handle the MPI environment and
the library paths.

### 1. Setup Environment & Pathing

Open **Run > Edit Configurations** and select the `aligner` target. Click the folder icon next to **Environment
variables** and add the following:

| Name               | Value                                | Purpose                            |
|:-------------------|:-------------------------------------|:-----------------------------------|
| `OMPI_MCA_patcher` | `^overwrite`                         | Prevents MPI/ASan memory conflicts |
| `ASAN_OPTIONS`     | `protect_shadow_gap=0`               | Allows MPI to use memory gaps      |
| `LIBRARY_PATH`     | `/home/abhinavmishra/micromamba/lib` | Links the SDSL-lite static library |

> **Note:** Ensure your **Working Directory** is set to the project root so the program can find the `/files` folder.

---

### 2. Test Case 1: DNA Global Alignment

Use this test case to verify the FM-Index and basic DNA sequence matching.

* **Program Arguments:**
    ```text
    --query files/dna1.fasta --target files/dna2.fasta --choice 1 --mode dna --verbose
    ```
* **Expected Output:** The console should show the MPI rank initialization, the building of the FM-Index for
  `dna2.fasta`, and the resulting alignment score.

---

### 3. Test Case 2: Protein Alignment

Use this test case to verify protein scoring matrices and sequence processing.

* **Program Arguments:**
    ```text
    --query files/prot1.fasta --target files/prot3.fasta --choice 1 --mode protein --verbose
    ```
* **Expected Output:** High-level logs showing the protein alphabet processing and the alignment of the amino acid
  sequences.

---

### 🛠 Troubleshooting "Illegal Instruction"

If the program crashes with `Caught signal 4 (Illegal instruction)`, ensure the `OMPI_MCA_patcher` variable is set
correctly. This error occurs because the **AddressSanitizer** and **OpenMPI** are attempting to manage the same memory
regions simultaneously.

### **Explanation of Options**

| Option                | Description                                                                        |
|-----------------------|------------------------------------------------------------------------------------|
| `--query`             | Path to the query FASTA file                                                       |
| `--target`            | Path to the target FASTA file                                                      |
| `--choice`            | Alignment method: <br> `1 = global` <br> `2 = local` <br> `3 = LCS` <br> `4 = all` |
| `--mode`              | Scoring mode: `dna` (uses EDNAFULL) or `protein` (uses BLOSUM62)                   |
| `--outdir` *(opt)*    | Output directory (default is current directory `.`)                                | 
| `--binary` *(opt)*    | Output binary file for dynamic programming matrix (default: `dp.bin`)              | 
| `--txt` *(opt)*       | Output text file for dynamic programming matrix (default: `dp.txt`)                |
| `--gap_open` *(opt)*  | Gap opening penalty (default: `-5.0`)                                              |
| `--gap_extend`\*(opt) | Gap extension penalty (default: `-1.0`)                                            |
| `--verbose` *(opt)*   | Show colored alignment and progress bars                                           |
| `--help`              | Show help and usage instructions                                                   |

---

## Input Format

The program accepts standard **FASTA** files.

---

## History

This project was originally developed in 2014 as part of the *Biological Computation* course during my Bachelor's in
Bioinformatics at JUIT, Solan. It was later revisited and optimized for better performance, readability, and
maintainability.

---

## License

This project is licensed under the [BSD 3-Clause License](./LICENSE)  
© 2025 Abhinav Mishra.