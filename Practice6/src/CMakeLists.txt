cmake_minimum_required(VERSION 3.25)
project(NumericalIntegrationExample)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(NumericalIntegrationExample
        Integration_Scheme.h
        Integration_Scheme_Interval.cpp
        Integration_Scheme_Interval.h
        main.cpp
        Point.cpp
        Point.h)
