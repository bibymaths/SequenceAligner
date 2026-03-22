# Sequence Alignment Benchmarking Framework

This repository contains a modular, **reproducible benchmarking system** for
comparing the performance and accuracy of widely used sequence alignment
tools. The framework focuses on **pairwise alignments** of DNA and
protein sequences and is designed to make fair comparisons between
functionally similar tools while clearly documenting limitations and
unsupported use cases.

## Supported tools

The benchmark evaluates the following external alignment programs:

| Tool              | Sequence type(s)                           | Notes                                                                                                                                                                                                                                                                                                                                                |
|-------------------|--------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **BLAST**         | DNA (via `blastn`), protein (via `blastp`) | The NCBI BLAST+ suite can align nucleotide or amino‑acid sequences directly between two FASTA files using the `-query` and `-subject` arguments. The benchmark captures tabular output (`-outfmt 6`) and parses the first high‑scoring pair to compute identity, alignment length, mismatches, gap openings and coverage【633730938334422†L361-L399】. |
| **Bowtie2**       | DNA only                                   | A fast DNA read aligner that requires building an FM‑index with `bowtie2‑build` and then aligning reads with `bowtie2`【354970811977329†L2124-L2160】. Bowtie2 does **not** support protein alignment and will be skipped for protein inputs.                                                                                                          |
| **BWA**           | DNA only                                   | The Burrows–Wheeler Aligner indexes a reference with `bwa index` and aligns reads with `bwa mem`【28747307351112†L32-L37】. It also lacks support for proteins and will be skipped for those runs.                                                                                                                                                     |
| **MAFFT**         | Protein                                    | A multiple sequence aligner that automatically detects sequence type and chooses an appropriate algorithm【395817269289108†L29-L36】. While MAFFT can align nucleotides, this benchmark restricts MAFFT to **protein** inputs to maintain methodological parity with Clustal Omega.                                                                    |
| **Clustal Omega** | Protein only                               | According to its documentation, Clustal Omega only aligns protein sequences and does not handle DNA/RNA inputs【738991639517845†L64-L66】. It accepts an input FASTA via `-i` and writes an aligned FASTA via `-o`, with `--outfmt=fasta` to enforce FASTA output【357385058027719†L312-L316】.                                                          |

Only standard, widely used tools are included. **Custom or in‑house aligners are deliberately excluded.**

## Directory structure

The benchmark code resides under the `benchmark/` package:

```
benchmark/
  benchmark.py        # CLI entrypoint
  utils.py            # Shared helpers for timing, memory tracking and FASTA I/O
  parsers/            # Output parsers for BLAST, SAM and aligned FASTA
  runners/            # Tool‑specific runners (BLAST, Bowtie2, BWA, MAFFT, Clustal)

configs/
  default.yaml        # Example configuration file (edit to specify your FASTA files)
results/              # Generated CSV/JSON results will be written here
logs/                 # Raw logs of each tool run
outputs/              # Temporary work directories for intermediate files
```

The `node_modules` and other assets are unrelated to the benchmarking framework and can be ignored.

## Installation and prerequisites

This framework is written in Python and depends on a few additional packages:

* **Python 3.8+**
* `PyYAML` for reading configuration files
* `psutil` for capturing peak memory usage

Install the Python dependencies (preferably in a virtual environment) with:

```bash
pip install pyyaml psutil
```

The benchmark **requires the alignment tools themselves** to be installed and
available on your system `PATH`. Specifically:

* `blastn` and `blastp` (BLAST+ package)
* `bowtie2` and `bowtie2-build`
* `bwa`
* `mafft`
* `clustalo` (Clustal Omega)

You need to install **5 external tools**. The cleanest way depends on your OS.

---

# 🔧 Recommended (Linux / WSL / HPC): use Conda

This is the most reproducible setup.

### 1. Install Miniconda (if not already)

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### 2. Create environment

```bash
conda create -n alignbench python=3.10 -y
conda activate alignbench
```

### 3. Install all tools

```bash
conda install -c bioconda -c conda-forge \
  blast bowtie2 bwa mafft clustalo -y
```

This installs:

* `blastn`, `blastp`
* `bowtie2`, `bowtie2-build`
* `bwa`
* `mafft`
* `clustalo`

---

### 4. Verify installation

```bash
blastn -version
bowtie2 --version
bwa
mafft --version
clustalo --version
```

If all respond → you're good.

---

# 🐧 Alternative: Ubuntu (APT + manual)

Some tools available directly:

```bash
sudo apt update
sudo apt install -y \
  ncbi-blast+ \
  bowtie2 \
  bwa \
  mafft \
  clustalo
```

⚠️ Note:

* Versions may be outdated
* Conda is preferred for consistency

---

# 🍎 macOS (Homebrew)

```bash
brew install blast bowtie2 bwa mafft clustal-omega
```

---

# 🧪 HPC cluster (common scenario for you)

Use module system or conda:

### Option A: modules

```bash
module load blast
module load bowtie2
module load bwa
module load mafft
module load clustal-omega
```

### Option B: user conda (preferred)

Same as Conda instructions above.

---

# 📦 Python dependencies (required)

```bash
pip install pyyaml psutil
```

---

# 🚀 Quick sanity test

Run:

```bash
which blastn
which bowtie2
which bwa
which mafft
which clustalo
```

If paths show → CLI will work.

---

# ⚠️ Important pitfalls (don’t ignore)

* `clustalo` binary name ≠ `clustal` → must be **clustalo**
* `blast` package gives multiple binaries (`blastn`, `blastp`)
* Bowtie2 needs BOTH:

    * `bowtie2`
    * `bowtie2-build`
* BWA requires indexing → handled automatically in your framework
* MAFFT writes to stdout → already handled in runner

---

# ✔ Final check with your framework

After installing:

```bash
python benchmark/benchmark.py --config configs/default.yaml
```

If tools are correctly installed:

* No “not found in PATH” errors
* Results will populate in `/results`

---

## Configuring inputs

Create a configuration file (YAML) specifying the paths to your input FASTA
files and benchmarking parameters. A sample file is provided at
`configs/default.yaml`:

```yaml
dna:
  query: path/to/query_dna.fasta
  target: path/to/target_dna.fasta
protein:
  query: path/to/query_protein.fasta
  target: path/to/target_protein.fasta
runs: 3
timeout: 300
threads: 4
```

* `dna.query`, `dna.target`: FASTA files containing one or more DNA sequences.
  Only BLAST (blastn), Bowtie2 and BWA will run on these.
* `protein.query`, `protein.target`: FASTA files with protein sequences. BLAST
  (blastp), MAFFT and Clustal Omega will be executed on these.
* `runs`: Number of times each tool is executed; results are aggregated
  (mean, median, standard deviation, minimum and maximum).
* `timeout`: Maximum seconds allowed for each tool invocation. Set to `null` to
  disable timeouts.
* `threads`: Number of threads passed to tools that support multithreading.

## Running the benchmark

Execute the benchmark from the repository root using the provided CLI script:

```bash
python benchmark/benchmark.py --config configs/default.yaml
```

The script will:

1. Create `results/`, `logs/` and `outputs/` directories as needed.
2. Loop over the specified number of runs.
3. For each sequence type (DNA and protein), call the appropriate runner
   functions according to the tool map:
    * DNA: BLAST (blastn), Bowtie2, BWA
    * Protein: BLAST (blastp), MAFFT, Clustal Omega
4. Measure **wall‑clock runtime** and **peak memory usage** using
   `psutil`, capturing the command, exit code, stdout and stderr.
5. Parse the output to compute accuracy metrics such as percentage identity,
   alignment length, mismatches, gap count and coverage. BLAST outputs are
   parsed from tabular format【633730938334422†L361-L399】; SAM outputs (Bowtie2/BWA) are parsed
   via the `NM` tag and CIGAR string; MAFFT/Clustal outputs are parsed from
   aligned FASTA.
6. Write aggregated runtime and memory statistics to CSV files
   (`results/runtime.csv`, `results/memory.csv`). Accuracy metrics are
   summarised per tool and metric in `results/accuracy.csv`.
7. Store the full per‑run results, including commands and raw metrics, in
   `results/full_results.json`, and record environment details (OS, CPU,
   RAM, tool versions) in `results/environment.json`.
8. Log raw stdout/stderr for each tool invocation under `logs/`.

### Interpreting the output

The generated CSV files have the following structure:

* **runtime.csv** – each row reports the mean, median, standard deviation,
  minimum and maximum wall‑clock time (in seconds) for a tool on a given
  sequence type.
* **memory.csv** – similar statistics for peak memory usage (in megabytes).
* **accuracy.csv** – for each metric (`identity`, `alignment_length`,
  `mismatches`, `gap_count`, `query_coverage`, `subject_coverage` and
  `target_coverage`), statistics are summarised across runs per tool and
  sequence type. Metrics that are not applicable (e.g. BWA on proteins) are
  recorded as `NA`.

The JSON files contain richer information:

* **full_results.json** – nested dictionary of per‑run results including
  commands, exit codes, raw runtime/memory numbers and parsed metrics.
* **environment.json** – details about the machine used for benchmarking,
  such as operating system, CPU cores, total memory, Python version and
  discovered tool versions.

## Methodological considerations

* **Fair comparisons:** DNA and protein aligners solve fundamentally
  different problems. The framework separates runs into DNA and protein
  categories and does not compare tools across these categories. For
  example, Bowtie2 and BWA results appear only under DNA.
* **Unsupported modes:** Bowtie2 and BWA only align DNA and will not be run
  on protein sequences. Clustal Omega cannot align DNA【738991639517845†L64-L66】. MAFFT is
  capable of aligning nucleotides but is limited to proteins in this
  benchmark to keep comparisons with Clustal Omega fair. When a tool
  does not support a given sequence type, its metrics are reported as
  `NA`.
* **First alignment only:** For BLAST, the framework considers only the
  first high‑scoring pair from the tabular output. MAFFT and Clustal
  produce pairwise MSAs; the parser computes identity and coverage by
  counting matches, mismatches and gaps across the aligned sequences.
* **Index caching:** Bowtie2 and BWA build indexes into `outputs/` and
  reuse them on subsequent runs to avoid redundant work. If you change
  the target FASTA, remove the corresponding files in `outputs/` to
  force index rebuilding.
* **Error handling:** If any command fails (non‑zero exit code or timeout),
  the exit code is recorded and all metrics for that run are set to `NA`.

## Extending the framework

The code is intentionally modular. To add a new tool:

1. Create a new module under `benchmark/runners/` that exposes a `run()`
   function following the existing signature. The function should build
   any required indices, invoke the tool via `utils.run_subprocess_with_resource_tracking`,
   write logs and return a dictionary with `runtime`, `memory` and
   `metrics` keys.
2. Add a parser under `benchmark/parsers/` if the tool’s output format is
   not already handled.
3. Register the runner in `benchmark/runners/__init__.py` and update the
   tool map in `benchmark/benchmark.py` to indicate which sequence types it
   supports.

## Licensing

This benchmark wraps third‑party alignment tools that each have their own
licenses. Consult the respective tool documentation for license terms.
The Python code in this repository is released under an open license (see
LICENSE file if provided).

---

Feel free to contribute improvements or report issues. Happy benchmarking!