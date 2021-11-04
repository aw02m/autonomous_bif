cmake_minimum_required(VERSION 3.21)

project(pd CXX)

add_executable(main main.cpp newton.cpp dynamical_system.cpp sys_func.cpp)
target_compile_options(main PUBLIC -O3 -mtune=native -march=native -mfpmath=both)