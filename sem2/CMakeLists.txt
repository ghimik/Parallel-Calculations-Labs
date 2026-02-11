cmake_minimum_required(VERSION 3.10)
project(MatrixInverse)

find_package(OpenMP REQUIRED)


set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(main
    src/main.cpp
    src/matrix.cpp
    src/inverse_sequential.cpp
    src/inverse_openmp.cpp
)

target_link_libraries(main OpenMP::OpenMP_CXX)

