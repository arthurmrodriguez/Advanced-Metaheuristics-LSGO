project (GAEDALib)

cmake_minimum_required (VERSION 2.8)

set (CMAKE_C_COMPILER   mpicc)
set (CMAKE_CXX_COMPILER mpicxx)

cmake_policy(SET CMP0042 NEW)
SET(CMAKE_MACOSX_RPATH OFF)

set (CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".dylib" ".a")

find_package (MPI REQUIRED)

find_library(LIBCONFIG config++)
find_path   (LIBCONFIG_INCLUDE libconfig.h++)

add_subdirectory (gaeda)      # generate libgaeda
add_subdirectory (gaedaexec)  # generate gaedaexec
add_subdirectory (problems)   # generate problems
add_subdirectory (utils)      # generate several utilities
