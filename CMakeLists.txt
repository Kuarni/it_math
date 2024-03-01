cmake_minimum_required(VERSION 3.22)
project(poisson C)

set(CMAKE_C_STANDARD 11)

add_subdirectory(argparse)
include_directories(argparse)

add_link_options(-fopenmp)
add_executable(poisson main.c)
target_link_libraries(poisson PRIVATE argparse)
