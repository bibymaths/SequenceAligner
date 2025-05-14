cmake_minimum_required(VERSION 3.10)
project(SequenceAligner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(aligner_v1 src/main_v1.cpp)
add_executable(aligner_v2 src/main_v2.cpp)
add_executable(aligner_v3 src/main_v3.cpp)

