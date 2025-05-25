import platform

# Detect OS and architecture
OS = platform.system()
ARCH = platform.machine()

# Build directory
BUILD_DIR = "build"

# Compiler configuration based on platform
if OS == "Linux":
    CMAKE_ARGS = []
elif OS == "Darwin" and ARCH == "x86_64":
    CMAKE_ARGS = [
        "-DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang",
        "-DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++",
        '-DCMAKE_C_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/llvm/include"',
        '-DCMAKE_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/llvm/include"',
        '-DCMAKE_EXE_LINKER_FLAGS="-L/usr/local/opt/llvm/lib -lomp"'
    ]
elif OS == "Darwin" and ARCH == "arm64":
    CMAKE_ARGS = [
        "-DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang",
        "-DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++",
        '-DCMAKE_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/llvm/include"',
        '-DCMAKE_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/llvm/include"',
        '-DCMAKE_EXE_LINKER_FLAGS="-L/opt/homebrew/opt/llvm/lib -lomp"'
    ]
else:
    raise RuntimeError("Unsupported platform or architecture")

rule all:
    input:
        "results/summary.png"

rule build:
    output:
        "aligner"
    params:
        cmake_args=" ".join(CMAKE_ARGS)
    shell:
        """
        mkdir -p build
        cd build && cmake .. {params.cmake_args} && make
        """

rule run_alignment:
    input:
        aligner="aligner",
        query="files/dna1.fasta",
        target="files/dna2.fasta"
    output:
        "results/global_dp_matrix.txt",
        "results/local_dp_matrix.txt",
        "results/lcs_traceback.txt"
    shell:
        """
        ./aligner \
          --query {input.query} \
          --target {input.target} \
          --choice 4 \
          --mode dna \
          --outdir results \
          --verbose \
          --gap_open -5 \
          --gap_extend -1 \
          --txt
        """

rule visualize:
    input:
        script="plotter/plotDP.sh",
        lcs="results/lcs_traceback.txt",
        globalalign="results/global_dp_matrix.txt",
        localalign="results/local_dp_matrix.txt"
    output:
        "results/summary.png"
    shell:
        """
        chmod +x {input.script}
        {input.script} {input.lcs} {input.globalalign} {input.localalign} results
        """

rule clean:
    shell:
        """
        rm -rf build aligner results
        """