cmake_minimum_required(VERSION 3.25)
project(src)

set(CMAKE_CXX_STANDARD 20)

add_executable(src
        main.cpp "fourier_transforms.hpp" fourier_transforms.cpp)
