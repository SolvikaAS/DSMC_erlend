cmake_minimum_required(VERSION 3.20)

project(dsmc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)

add_executable(dsmc
    src/main.cpp
    src/Simulation.cpp
    src/Initialization.cpp
    src/Collision.cpp
    src/Boundary.cpp
    src/Geometry.cpp
    src/Sampling.cpp
    src/MPIManager.cpp
    src/Gamma.cpp
    rng/BRNG.cxx
    rng/DiscreteUniformRNG.cxx
    rng/PoissonRNG.cxx
    rng/RayleighRNG.cxx
    rng/IsotropicVectorRNG.cxx
)

target_include_directories(dsmc PRIVATE include rng)

target_compile_definitions(dsmc PRIVATE
    DSMC_RESULTS_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}/results"
)

target_link_libraries(dsmc PRIVATE MPI::MPI_CXX)
