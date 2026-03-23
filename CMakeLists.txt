cmake_minimum_required(VERSION 3.15)
project(SequenceAligner LANGUAGES C CXX)
option(DOCKER_BUILD "Build configuration for Docker image" OFF)

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

# Optional SIMD acceleration if supported
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)

# ------------------------------------------------------------------------------
# Shared settings helper
# ------------------------------------------------------------------------------
function(configure_sequence_target target_name)
  target_include_directories(${target_name} PRIVATE
          ${MPI_CXX_INCLUDE_DIRS}
  )

  target_link_directories(${target_name} PRIVATE
          ${MPI_CXX_LIBRARIES}
  )

  target_link_libraries(${target_name} PRIVATE
          ${MPI_CXX_LIBRARIES}
          OpenMP::OpenMP_CXX
          divsufsort
          divsufsort64
          pthread
  )

  target_compile_options(${target_name} PRIVATE
          -Wall
          -Wextra
          -Wno-cast-function-type
  )

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug build for ${target_name}: keeping symbols for Valgrind")
    target_compile_options(${target_name} PRIVATE -g -fno-omit-frame-pointer)
  endif()

  if(COMPILER_SUPPORTS_AVX2)
    target_compile_options(${target_name} PRIVATE -mavx2 -march=native -funroll-loops -O3)
  else()
    target_compile_options(${target_name} PRIVATE -O3)
    message(WARNING "Compiler does not support AVX2 for ${target_name}; using -O3 only.")
  endif()

  set_target_properties(${target_name} PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
endfunction()

# ------------------------------------------------------------------------------
# fmindex.cpp -> object file target
# ------------------------------------------------------------------------------
add_library(fmindex_obj OBJECT src/fmindex.cpp)
configure_sequence_target(fmindex_obj)

set_target_properties(fmindex_obj PROPERTIES
        OBJECT_OUTPUTS "${CMAKE_BINARY_DIR}/fmindex.o"
)

# ------------------------------------------------------------------------------
# aligner executable from main.cpp
# ------------------------------------------------------------------------------
add_executable(aligner src/main.cpp)
configure_sequence_target(aligner)

# ------------------------------------------------------------------------------
# seed_aligner executable
# In Docker, avoid linking fmindex_obj if fmindex.cpp also defines main().
# Outside Docker, preserve the current behavior.
# ------------------------------------------------------------------------------
if(DOCKER_BUILD)
  add_executable(seed_aligner src/main_fmindex.cpp)
else()
  add_executable(seed_aligner
          src/main_fmindex.cpp
          $<TARGET_OBJECTS:fmindex_obj>
  )
endif()
configure_sequence_target(seed_aligner)