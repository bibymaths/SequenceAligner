cmake_minimum_required(VERSION 3.10)
project(SequenceAligner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(MPI REQUIRED)

# macOS: manually configure OpenMP if using LLVM from Homebrew
if(APPLE)
  message(STATUS "Detected macOS: manually setting OpenMP flags for LLVM")
  set(OpenMP_C_FLAGS "-fopenmp")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY "/usr/local/opt/llvm/lib/libomp.dylib")
  set(OpenMP_CXX_FLAGS "-fopenmp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_CXX_LIBRARIES "/usr/local/opt/llvm/lib/libomp.dylib")
endif()

find_package(OpenMP REQUIRED)

# Optional: Debug build sanitizer
if (CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "Enabling AddressSanitizer")
  set(ASAN_FLAGS "-fsanitize=address -fno-omit-frame-pointer -g")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ASAN_FLAGS}")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${ASAN_FLAGS}")
endif()

add_executable(aligner src/main.cpp)

target_include_directories(aligner PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(aligner PRIVATE ${MPI_CXX_LIBRARIES} OpenMP::OpenMP_CXX)

# Optional SIMD acceleration
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
if(COMPILER_SUPPORTS_AVX2)
  target_compile_options(aligner PRIVATE -mavx2 -O3 -march=native -funroll-loops)
else()
  message(WARNING "Compiler does not support AVX2. Skipping -mavx2.")
endif()

target_compile_options(aligner PRIVATE -Wall -Wextra -Wno-cast-function-type)

set_target_properties(aligner PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)
