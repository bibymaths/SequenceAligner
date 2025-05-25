cmake_minimum_required(VERSION 3.10)
project(SequenceAligner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
#find_package(ZLIB REQUIRED)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "Enabling AddressSanitizer")
  set(ASAN_FLAGS "-fsanitize=address -fno-omit-frame-pointer -g")
  set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} ${ASAN_FLAGS}")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${ASAN_FLAGS}")
endif()

add_executable(aligner src/main.cpp)

target_include_directories(aligner PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(aligner PRIVATE ${MPI_CXX_LIBRARIES} OpenMP::OpenMP_CXX)
# target_link_libraries(aligner PRIVATE ${MPI_CXX_LIBRARIES} OpenMP::OpenMP_CXX ${ZLIB_LIBRARIES})

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
if(COMPILER_SUPPORTS_AVX2)
  target_compile_options(aligner PRIVATE -mavx2 -O3 -march=native -funroll-loops)
else()
  message(WARNING "Compiler does not support AVX2. Skipping -mavx2.")
endif()

target_compile_options(aligner PRIVATE -Wall -Wextra)
target_compile_options(aligner PRIVATE
    -Wno-cast-function-type
)

set_target_properties(aligner PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)