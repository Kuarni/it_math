cmake_minimum_required(VERSION 3.22)
project(poisson C)

set(CMAKE_C_STANDARD 11)

add_subdirectory(argparse)
include_directories(argparse)

add_subdirectory(tasks)
include_directories(tasks)

add_link_options(-fopenmp)
add_compile_options(-fopenmp)

add_executable(poisson it_math_1.c
        it_math_1.h)

set(CMAKE_EXE_LINKER_FLAGS "-static")
target_link_libraries(poisson PRIVATE argparse)
target_link_libraries(poisson PUBLIC tasks)


