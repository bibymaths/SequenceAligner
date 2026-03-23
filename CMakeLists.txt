cmake_minimum_required(VERSION 3.15)
project(SequenceAligner LANGUAGES C CXX)
option(DOCKER_BUILD "Build configuration for Docker image" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)

if(APPLE)
  message(STATUS "Detected macOS: setting OpenMP flags for Homebrew LLVM")
  set(OpenMP_C_FLAGS "-fopenmp")
  set(OpenMP_CXX_FLAGS "-fopenmp")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY "/usr/local/opt/llvm/lib/libomp.dylib")
endif()

find_package(OpenMP REQUIRED)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)

function(configure_sequence_target target_name)
  target_include_directories(${target_name} PRIVATE
          ${MPI_CXX_INCLUDE_DIRS}
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

add_executable(fmindex src/fmindex.cpp)
configure_sequence_target(fmindex)

add_executable(aligner src/main.cpp)
configure_sequence_target(aligner)

add_executable(seed_aligner src/main_fmindex.cpp)
configure_sequence_target(seed_aligner)