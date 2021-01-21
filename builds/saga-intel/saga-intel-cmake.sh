# Run this script by typing: source saga-intel-cmake
# After this script completes, type: make -j 8
# If rebuilding, type 'make clean' before running 'make -j 8'

# This cmake configuration script is set up to perform a parallel build without Trilinos

# Set path to top cism directory
# Note, this is an easy way to build out of source.
# In directory you want to build in, run:
#  source $CISM/builds/saga-intel/saga-intel-cmake $CISM
# where $CISM is the path to the top level cism directory.
if [ $# -eq 0 ]
then
    cism_top="../.." 
else
    cism_top=${1%/}
fi

echo CISM: "${cism_top}"

module purge
module load netCDF/4.7.1-iimpi-2019b
module load netCDF-Fortran/4.5.2-iimpi-2019b
module load CMake/3.15.3-GCCcore-8.3.0
module load OpenBLAS/0.3.7-GCC-8.3.0

## Resulting module list 
#  1) StdEnv                              (S)
#  2) GCCcore/8.3.0
#  3) zlib/1.2.11-GCCcore-8.3.0           (H)
#  4) binutils/2.32-GCCcore-8.3.0         (H)
#  5) iccifort/2019.5.281
#  6) impi/2018.5.288-iccifort-2019.5.281
#  7) iimpi/2019b
#  8) Szip/2.1.1-GCCcore-8.3.0            (H)
#  9) HDF5/1.10.5-iimpi-2019b
# 10) cURL/7.66.0-GCCcore-8.3.0           (H)
# 11) netCDF/4.7.1-iimpi-2019b
# 12) netCDF-Fortran/4.5.2-iimpi-2019b
# 13) ncurses/6.1-GCCcore-8.3.0           (H)
# 14) bzip2/1.0.8-GCCcore-8.3.0           (H)
# 15) CMake/3.15.3-GCCcore-8.3.0
# 16) GCC/8.3.0
# 17) OpenBLAS/0.3.7-GCC-8.3.0

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
