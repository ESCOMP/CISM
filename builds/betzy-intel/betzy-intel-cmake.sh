# Run this script by typing: source betzy-intel-cmake
# After this script completes, type: make -j 8
# If rebuilding, type 'make clean' before running 'make -j 8'

# This cmake configuration script is set up to perform a parallel build without Trilinos

# Set path to top cism directory
# Note, this is an easy way to build out of source.
# In directory you want to build in, run:
#  source $CISM/builds/fram-intel/fram-intel-cmake $CISM
# where $CISM is the path to the top level cism directory.
if [ $# -eq 0 ]
then
    cism_top="../.." 
else
    cism_top=${1%/}
fi

echo CISM: "${cism_top}"

#module purge
#module load StdEnv
#module load intel/2020a
#module load netCDF/4.7.1-iimpi-2019b
#module load netCDF-Fortran/4.5.2-iimpi-2019b

module purge
module load StdEnv
module load intel/2020a
module load netCDF/4.7.4-iimpi-2020a
module load netCDF-Fortran/4.5.2-iimpi-2020a
module load CMake/3.16.4-GCCcore-9.3.0

## Resulting module list 

# remove old build data:
rm -f ./CMakeCache.txt
rm -fr ./CMakeFiles

echo
echo "Doing CMake Configuration step"

# A few non-intuitive things:
#
# - CISM_FORCE_FORTRAN_LINKER: without this, cmake tries to use a C++ linker, which doesn't work
#
# - CISM_INCLUDE_IMPLICIT_LINK_LIBRARIES: if this is on (the default), some
#   libraries are included on the link line which can't be found (e.g.,
#   hdf5). This may be related to the fact that trilinos on yellowstone is old,
#   and/or the fact that cmake wants to use a C++ linker but we're telling it to
#   use a fortran linker.

cmake \
  -D CISM_BUILD_CISM_DRIVER:BOOL=ON \
  -D CISM_ENABLE_BISICLES=OFF \
  -D CISM_ENABLE_FELIX=OFF \
\
  -D CISM_USE_TRILINOS:BOOL=OFF \
  -D CISM_MPI_MODE:BOOL=ON \
  -D CISM_SERIAL_MODE:BOOL=OFF \
\
  -D CISM_USE_GPTL_INSTRUMENTATION:BOOL=OFF \
  -D CISM_COUPLED:BOOL=OFF \
  -D CISM_USE_CISM_FRONT_END:BOOL=OFF \
\
  -D CISM_TRILINOS_DIR=$NOTSPECIFIED_TRILINOS_ROOT \
  -D CISM_NETCDF_DIR=$EBROOTNETCDF \
  -D CISM_FORCE_FORTRAN_LINKER:BOOL=OFF \
  -D CISM_INCLUDE_IMPLICIT_LINK_LIBRARIES:BOOL=ON \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
\
  -D CMAKE_CXX_COMPILER=mpiicpc \
  -D CMAKE_C_COMPILER=mpiicc \
  -D CMAKE_Fortran_COMPILER=mpiifort \
\
  -D CMAKE_Fortran_FLAGS:STRING="-fp-model source -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -xHost -O2 -mkl" \
  -D CMAKE_C_FLAGS:STRING="-O2 -fp-model precise -xHost" \
  -D CMAKE_CXX_FLAGS:STRING="-O2 -fp-model precise -xHost" \
\
  -D CISM_EXTRA_LIBS:STRING="-lnetcdff" \
  "${cism_top}"

# Note: last argument above  "../.." or ${cism_top} is path to top cism directory
