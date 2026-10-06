cmake_minimum_required(VERSION 3.20)

project(dsmc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

add_executable(dsmc
    src/main.cpp
)

target_include_directories(dsmc PRIVATE include)

target_link_libraries(dsmc PRIVATE MPI::MPI_CXX)
