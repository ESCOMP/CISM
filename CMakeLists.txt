cmake_minimum_required(VERSION 3.10)

# Brings in CIME CMake modules for tests
list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})

# Need to set up the build flags and pfunit stuff?
include(CIME_initial_setup)

project(CISM_sta_units LANGUAGES Fortran VERSION 0.1)

include(CIME_utils)

add_subdirectory(libglad)
add_subdirectory(libglimmer)
add_subdirectory(tests/pfunit)




