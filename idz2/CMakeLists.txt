cmake_minimum_required(VERSION 3.22)
project(Fluid)

set(CMAKE_CXX_STANDARD 20)



add_executable(Fluid
        idz2.cpp)

add_compile_definitions(TYPES=$ENV{TYPES})
add_compile_definitions(SIZES=$ENV{SIZES})
