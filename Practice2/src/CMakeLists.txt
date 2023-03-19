cmake_minimum_required(VERSION 3.20)

project(CubicInterpolationSpline)
set(${CMAKE_CXX_STANDARD} 23)

set(headers
        Cubic_Interpolation_Spline_1D.h
        Point.h
        Spline.h
        regular_partions.h
        )

set(source
        Cubic_Interpolation_Spline_1D.cpp
        Point.cpp
        regular_partions.cpp
        main.cpp)

add_executable(${CubicInterpolationSpline}
                ${headers}
                ${source})