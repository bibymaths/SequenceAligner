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
  set(OpenMP_omp_LIBRARY "/usr/local/opt/llvm/lib/libomp.dylib")
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find LZMA (XZ) compression library
find_library(LZMA_LIBRARIES NAMES lzma HINTS /usr/local/lib /usr/lib /usr/bin)
find_path(LZMA_INCLUDE_DIRS lzma.h HINTS /usr/local/include /usr/include /usr/bin)

if(NOT LZMA_LIBRARIES OR NOT LZMA_INCLUDE_DIRS)
  message(FATAL_ERROR "liblzma not found. Make sure XZ is installed (e.g., brew install xz)")
endif()

# Enable AddressSanitizer in Debug builds
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Enabling AddressSanitizer")
  target_compile_options(aligner_private_sanitizers PRIVATE -fsanitize=address -fno-omit-frame-pointer -g)
  target_link_options(aligner_private_sanitizers PRIVATE -fsanitize=address)
endif()

# Create executable
add_executable(aligner src/test.cpp)

target_include_directories(aligner PRIVATE
  ${MPI_CXX_INCLUDE_DIRS}
  ${LZMA_INCLUDE_DIRS}
)

# Include directories
target_link_libraries(aligner PRIVATE
  ${MPI_CXX_LIBRARIES}
  OpenMP::OpenMP_CXX
  ${LZMA_LIBRARIES}
)

# Common compiler warnings and flags
target_compile_options(aligner PRIVATE
  -Wall
  -Wextra
  -Wno-cast-function-type
)

# Optional SIMD acceleration if supported
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
if(COMPILER_SUPPORTS_AVX2)
  target_compile_options(aligner PRIVATE -mavx2 -O3 -march=native -funroll-loops)
else()
  message(WARNING "Compiler does not support AVX2; skipping related flags.")
endif()

# Place the built executable in the project root
set_target_properties(aligner PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)
