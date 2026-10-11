cmake_minimum_required(VERSION 3.20)

project(dsmc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)

add_library(dsmc_core STATIC
    src/Simulation.cpp
    src/Initialization.cpp
    src/Collision.cpp
    src/Boundary.cpp
    src/Geometry.cpp
    src/Sampling.cpp
    src/MPIManager.cpp
    rng/BRNG.cxx
    rng/DiscreteUniformRNG.cxx
    rng/PoissonRNG.cxx
    rng/RayleighRNG.cxx
    rng/IsotropicVectorRNG.cxx
)

target_include_directories(dsmc_core PUBLIC include rng)

target_compile_definitions(dsmc_core PRIVATE
    DSMC_RESULTS_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}/results"
)

target_link_libraries(dsmc_core PUBLIC MPI::MPI_CXX)

add_executable(dsmc src/main.cpp)
target_link_libraries(dsmc PRIVATE dsmc_core)

include(CTest)
if(BUILD_TESTING)
    add_executable(geometry_tests tests/GeometryTests.cpp)
    target_link_libraries(geometry_tests PRIVATE dsmc_core)
    add_test(NAME box_geometry_and_reflection COMMAND geometry_tests)
endif()
