cmake_minimum_required(VERSION 3.15)
project(SequenceAligner LANGUAGES C CXX)
option(DOCKER_BUILD "Build configuration for Docker image" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)

if(APPLE)
  message(STATUS "Detected macOS: setting OpenMP flags for Homebrew libomp")
  set(OpenMP_C_FLAGS   "-Xpreprocessor -fopenmp")
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
  set(OpenMP_C_LIB_NAMES   "omp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  # Resolve libomp path dynamically via brew — works on both Intel and Apple Silicon
  execute_process(
          COMMAND brew --prefix libomp
          OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
          OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set(OpenMP_omp_LIBRARY "${HOMEBREW_LIBOMP_PREFIX}/lib/libomp.dylib")
  message(STATUS "libomp path: ${OpenMP_omp_LIBRARY}")
endif()

find_package(OpenMP REQUIRED)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)

function(configure_sequence_target target_name)
  target_include_directories(${target_name} PRIVATE
          ${MPI_CXX_INCLUDE_DIRS}
  )

  if(APPLE)
    execute_process(
            COMMAND brew --prefix libomp
            OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    target_include_directories(${target_name} PRIVATE
            "${HOMEBREW_LIBOMP_PREFIX}/include"
    )
    target_link_directories(${target_name} PRIVATE
            "${HOMEBREW_LIBOMP_PREFIX}/lib"
    )
  endif()

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

  # AVX2 + immintrin.h only on x86_64 — not available on ARM64
  if(COMPILER_SUPPORTS_AVX2 AND NOT CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64|ARM64")
    target_compile_options(${target_name} PRIVATE -mavx2 -march=native -funroll-loops -O3)
  else()
    target_compile_options(${target_name} PRIVATE -O3)
    message(WARNING "Skipping AVX2 for ${target_name} (ARM or unsupported); using -O3 only.")
  endif()

  set_target_properties(${target_name} PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
endfunction()

add_executable(fmindex src/fmindex.cpp)
configure_sequence_target(fmindex)

add_executable(aligner src/main.cpp)
configure_sequence_target(aligner)

add_executable(seed_aligner src/main_fmindex.cpp)
configure_sequence_target(seed_aligner)