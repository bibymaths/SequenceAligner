# CMakeLists.txt for the SequenceAligner project
# This file is used to configure the build process using CMake
# CMake version 3.10 or higher is required
cmake_minimum_required(VERSION 3.10)
project(SequenceAligner)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find MPI
find_package(MPI REQUIRED)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Add executables
add_executable(aligner src/main.cpp)

# MPI + OpenMP + AVX2
target_include_directories(aligner PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(aligner PRIVATE ${MPI_CXX_LIBRARIES} OpenMP::OpenMP_CXX)
target_compile_options(aligner PRIVATE -mavx2)