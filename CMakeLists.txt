cmake_minimum_required(VERSION 3.15)
project(SequenceAligner LANGUAGES C CXX)

# Enforce C++17 and disable GNU extensions
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find MPI
find_package(MPI REQUIRED)

# macOS: OpenMP via Homebrew LLVM
if(APPLE)
  message(STATUS "Detected macOS: setting OpenMP flags for Homebrew LLVM")
  set(OpenMP_C_FLAGS   "-fopenmp")
  set(OpenMP_CXX_FLAGS "-fopenmp")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  # Ensure this path is correct for your Homebrew LLVM installation
  set(OpenMP_omp_LIBRARY "/usr/local/opt/llvm/lib/libomp.dylib") # Or use find_library
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)

# --- SDSL-lite Configuration ---
# Define the paths to your sdsl-lite installation (mambaforge in your case)
set(SDSL_INCLUDE_DIR "/home/abhinavmishra/mambaforge/include")
set(SDSL_LIBRARY_DIR "/home/abhinavmishra/mambaforge/lib")

# Create executable (assuming your source file is src/fm.cpp for the aligner)
# If your main source file is still src/test.cpp, adjust accordingly.
# Let's assume you want to build an executable named 'aligner' from 'src/fm.cpp'
add_executable(aligner src/test.cpp) # Changed 'aligner' to 'aligner' and 'src/test.cpp' to 'src/fm.cpp'

# Include directories
target_include_directories(aligner PRIVATE
  ${MPI_CXX_INCLUDE_DIRS}
  ${SDSL_INCLUDE_DIR}      # Add SDSL include directory
)

# Link directories (so CMake knows where to find -lsdsl etc.)
target_link_directories(aligner PRIVATE
  ${SDSL_LIBRARY_DIR}      # Add SDSL library directory
)

# Link libraries
target_link_libraries(aligner PRIVATE
  ${MPI_CXX_LIBRARIES}
  OpenMP::OpenMP_CXX
  sdsl                 # Link against libsdsl.a or libsdsl.so
  divsufsort           # Link against libdivsufsort.a or libdivsufsort.so
  divsufsort64         # Link against libdivsufsort64.a or libdivsufsort64.so
  pthread              # Link against pthread
)

# Common compiler warnings and flags
target_compile_options(aligner PRIVATE
  -Wall
  -Wextra
  -Wno-cast-function-type
  # -O3                  # Moved -O3 to be applied with AVX2 or generally
  # -DNDEBUG             # Usually set by CMAKE_BUILD_TYPE=Release
)

# Enable AddressSanitizer in Debug builds
# Note: This target_compile_options and target_link_options should apply to aligner
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Enabling AddressSanitizer for aligner")
  target_compile_options(aligner PRIVATE -fsanitize=address -fno-omit-frame-pointer -g)
  target_link_options(aligner PRIVATE -fsanitize=address)
endif()


# Optional SIMD acceleration if supported
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
if(COMPILER_SUPPORTS_AVX2)
  target_compile_options(aligner PRIVATE -mavx2 -march=native -funroll-loops -O3)
else()
  target_compile_options(aligner PRIVATE -O3) # Apply -O3 even if AVX2 is not supported
  message(WARNING "Compiler does not support AVX2; skipping AVX2-specific flags but keeping -O3.")
endif()

# Place the built executable in the project root
set_target_properties(aligner PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)