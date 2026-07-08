cmake_minimum_required(VERSION 3.16)
project(MeshThings VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_HEAT_DEMO "Build MPI heat conduction demo" ON)

include(FetchContent)
set(OVM_ENABLE_APPLICATIONS OFF CACHE BOOL "" FORCE)
set(OVM_ENABLE_UNITTESTS OFF CACHE BOOL "" FORCE)
set(OVM_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(OVM_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    OpenVolumeMesh
    GIT_REPOSITORY https://github.com/OpenVolumeMesh/OpenVolumeMesh.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(OpenVolumeMesh)

find_package(MPI REQUIRED)
find_library(METIS_LIBRARY metis REQUIRED)
find_library(HYPRE_LIBRARY HYPRE REQUIRED)
find_path(HYPRE_INCLUDE_DIR NAMES HYPRE.h PATH_SUFFIXES hypre)

add_library(meshthings_core STATIC
    src/plcParser.cpp
    src/concreteDelaunayTriangulation.cpp
    src/delaunayTriangulationFactory.cpp
    src/constrainedDelaunay.cpp
    src/wedgeExtruder.cpp
    src/openVolumeMeshAdapter.cpp
    src/meshPartitioner.cpp
    src/heatFEM.cpp
    src/hypreSolver.cpp
)

target_include_directories(meshthings_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${HYPRE_INCLUDE_DIR}
)

target_link_libraries(meshthings_core
    PUBLIC
        OpenVolumeMesh::OpenVolumeMesh
        MPI::MPI_CXX
        ${METIS_LIBRARY}
        ${HYPRE_LIBRARY}
)

if(BUILD_HEAT_DEMO)
    add_executable(heat_conduction_mpi examples/heat_conduction_mpi.cpp)
    target_link_libraries(heat_conduction_mpi PRIVATE meshthings_core MPI::MPI_CXX)
endif()
