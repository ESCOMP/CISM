# cism_driver and CISM front-end build

# Need mod files from glimmer
include_directories (${CISM_BINARY_DIR}/include  ${PYTHON_INC_DIR})

link_directories (${CISM_HDF5_LIB_DIR} ${PYTHON_LIB_DIR} )

# DMR 6/2/14 -- Moved simple_forcing.F90 and testsfg.F90 to the glimmercismfortran lib.
# These are local source files needed to make the cism_driver executable and CISM front-end
add_executable(cism_driver cism_driver.F90 cism_front_end.F90 gcm_cism_interface.F90
               gcm_to_cism_glint.F90)

#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
IF (CISM_STATIC_LINKING)
  SET_TARGET_PROPERTIES(cism_driver PROPERTIES LINK_SEARCH_START_STATIC 1)
  SET_TARGET_PROPERTIES(cism_driver PROPERTIES LINK_SEARCH_END_STATIC 1)
ENDIF()

IF (CISM_FORCE_FORTRAN_LINKER)
  SET_PROPERTY(TARGET cism_driver PROPERTY LINKER_LANGUAGE Fortran)
ENDIF()

MESSAGE("CISM_HDF5_LIBS: " ${CISM_HDF5_LIBS})

# Executable depends on several glimmer libraries,
# and potentially an f90main.o file (on jaguar with PGI at least)

target_link_libraries(cism_driver
    ${CISM_FMAIN}
    glimmercismfortran
    ${PYTHON_LIBS}
    ${CISM_NETCDF_LIBS}
    ${CISM_HDF5_LIBS}
    ${CISM_MPI_LIBS}
    ${CISM_IMPLICIT_LINK_LIBS}
    ${CISM_EXTRA_LIBS}
    ${CISM_GPTL_LIB}
)

# Helpful(?) message near end of configuration step
MESSAGE("")
MESSAGE("   Executable cism_driver should appear in dir: build_dir/cism_driver")
MESSAGE("")
