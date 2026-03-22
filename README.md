<img src="docs/assets/logo.png" alt="SequenceAligner" width="300">

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

### Installing [OpenMPI](https://www.open-mpi.org/software/ompi/v5.0/) for v2 and v3

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

# 🛠️ One-Time Library & Environment Setup

Follow these steps to ensure all dependencies are installed and the system is configured to handle the
AddressSanitizer (ASan) and `libdivsufsort` requirements.

---

## 1. Install Libraries via Micromamba

We use Micromamba to manage the `sdsl-lite` and `libdivsufsort` dependencies. This avoids manual compilation errors with
GCC 15.

```bash
# 1. Install the SDSL and DivSufSort libraries
mamba install -c conda-forge sdsl-lite libdivsufsort

# 2. Identify your environment path for CLion
echo $CONDA_PREFIX
````

**Note:**
The path returned (usually `/home/USER/micromamba`) is what you will use for the `CMAKE_PREFIX_PATH`.

---

## 2. Configure CLion Settings

To make the libraries visible to your project, update the CMake settings in the CLion GUI:

* Open **Settings (Ctrl+Alt+S)**
* Navigate to:
  `Build, Execution, Deployment > CMake`

In **CMake options**, paste:

```plaintext
-DCMAKE_PREFIX_PATH=/home/YOUR_USER/micromamba
```

Click **Apply**, then click the **Reload CMake Project** icon in the CMake tab.

---

## 3. Fix the libasan Version Bridge

Fedora 41+ uses a newer version of AddressSanitizer. Since the project links against version 8, you must create a
symbolic link (bridge):

```bash
# Find your actual system ASan version and link it to version 8
ACTUAL_ASAN=$(ls /usr/lib64/libasan.so.[1-9]* | head -n 1)
sudo ln -sf $ACTUAL_ASAN /usr/lib64/libasan.so.8.0.0
```

---

## 4. Mandatory Run Settings (Signal 4 Fix)

To prevent the **Illegal Instruction** crash caused by conflicts between OpenMPI and AddressSanitizer, add the following
environment variables to every Run Configuration (`aligner` and `fmindex`):

* Go to:
  `Run > Edit Configurations > Environment Variables`

Paste:

```plaintext
OMPI_MCA_patcher=^overwrite;ASAN_OPTIONS=protect_shadow_gap=0;LIBRARY_PATH=/home/YOUR_USER/micromamba/lib
```

---

# 🏁 Setup Verification Checklist

* **[SDSL](https://anaconda.org/channels/conda-forge/packages/sdsl-lite/overview) Check:**

  ```bash
  ls /home/YOUR_USER/micromamba/lib/libsdsl.a
  ```

* **[DivSufSort](https://github.com/y-256/libdivsufsort) Check:**

  ```bash
  ls /home/YOUR_USER/micromamba/lib/libdivsufsort.a
  ```

* **[ASan](https://gnu.googlesource.com/gcc/+/362cbc2d5f18c9f00dc3b945fb01c43bb0d36aae/libasan) Check:**

  ```bash
  ls -l /usr/lib64/libasan.so.8.0.0
  ```

  Should point to your actual system library.

* **MPI Check:**
  Ensure `OMPI_MCA_patcher=^overwrite` is present in your CLion environment variables.

---

## 🛠️ CLion Configuration Guide

This project consists of two tools: `fmindex` (for pre-processing) and `aligner` (for sequence matching). Both require
specific Environment Variables to handle the **OpenMPI + AddressSanitizer** memory conflict.

### 1. Global Environment Setup

For **ALL** run configurations below, copy and paste this single string into the **Environment variables** field in
CLion:

**Copy this:**
`OMPI_MCA_patcher=^overwrite;ASAN_OPTIONS=protect_shadow_gap=0;LIBRARY_PATH=/home/abhinavmishra/micromamba/lib`

---

### 2. Phase 1: Generating FM-Indexes

Before aligning, you must generate `.fmidx` files from your FASTA data. Create two **CMake Application** configurations
in CLion:

#### **Config A: Index DNA 1**

* **Target:** `fmindex`
* **Program arguments:** `files/dna1.fasta -s $`
* **Purpose:** Generates an index for the query sequence.

#### **Config B: Index DNA 2**

* **Target:** `fmindex`
* **Program arguments:** `files/dna2.fasta -s $`
* **Purpose:** Generates an index for the target sequence.

---

### 3. Phase 2: Running the Aligner

Once the `.fmidx` files are generated in your project root, switch to the `aligner` configuration:

* **Target:** `aligner`
* **Program arguments (DNA Test):**
  `--query files/dna1.fasta --target files/dna2.fasta --choice 1 --mode dna --verbose`
* **Working Directory:** Must be set to the project root (e.g., `/home/abhinavmishra/git/SequenceAligner`) so the app
  can find the `/files` folder.

---

### 📝 Summary for Quick Setup

| Configuration Name   | Target    | Program Arguments                                                                    |
|:---------------------|:----------|:-------------------------------------------------------------------------------------|
| **Generate Index 1** | `fmindex` | `files/dna1.fasta -s $`                                                              |
| **Generate Index 2** | `fmindex` | `files/dna2.fasta -s $`                                                              |
| **Run DNA Aligner**  | `aligner` | `--query files/dna1.fasta --target files/dna2.fasta --choice 1 --mode dna --verbose` |

---

### ⚠️ Common Error: Signal 4 (Illegal Instruction)

If the program crashes immediately upon launch, verify that `OMPI_MCA_patcher=^overwrite` is present in your *
*Environment variables**. This prevents OpenMPI from corrupting memory regions monitored by the AddressSanitizer.

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

### 📝 Note on FM-Index Anchoring with Proteins

If you see the message `"FM-index anchoring unavailable/failed. Falling back to MPI full DP"` during protein alignment,
this is **expected behavior, not a bug.**

* **Exact Matches vs. Similarity:** The FM-Index relies on finding exact substring matches (k-mers, typically 5-8
  characters long) to build anchors. Distant protein sequences (e.g., <30% identity) often preserve *chemical
  similarity* rather than exact character identity, meaning they may only share very short exact matches (2-3 amino
  acids).
* **Smart Fallback:** Because the sequences lack exact matches long enough to safely anchor the alignment, the program
  intentionally skips the FM-Index phase. It gracefully falls back to the full Smith-Waterman or Needleman-Wunsch DP
  matrix (using BLOSUM62) to ensure a biologically accurate alignment.
* **Why not lower the k-mer size?** Forcing a tiny k-mer threshold (like `k=3`) on proteins would result in massive
  amounts of random, noisy seeds, completely destroying both accuracy and performance.

---

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