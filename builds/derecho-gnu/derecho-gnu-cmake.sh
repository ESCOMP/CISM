#!/bin/sh

# After this script completes, type: make -j 8
# If rebuilding, type 'make clean' before running 'make -j 8'

# Set path to top cism directory
# Note, this is an easy way to build out of source.
# In directory you want to build in, run:
#   $ source $CISM/builds/derecho-gnu/derecho-gnu-cmake.sh $CISM
# where $CISM is the path to the top level cism directory.
if [ $# -eq 0 ]
then
    cism_top="../.." 
else
    cism_top=${1%/}
fi

source /etc/profile.d/z00_modules.sh

source "${cism_top}/builds/derecho-gnu/derecho-gnu-modules"

echo CISM: "${cism_top}"

# remove old build data:
rm -f ./CMakeCache.txt
rm -rf ./CMakeFiles

echo
echo "Doing CMake Configuration step"

# Note: the compilation flags were taken from the defaults for a CESM build on
# cheyenne-intel (using cime at 84aafd5). Some of these options are probably
# unnecessary for a standalone cism build, but I am keeping things consistent
# with the CESM build for simplicity.

# GPTL timing is off by default. To use it, first build libgptl.a in
# utils/libgptl (see utils/libgptl/README.CISM), then set
# CISM_USE_GPTL_INSTRUMENTATION=ON and CISM_GPTL_DIR=<path to libgptl> below.

cmake \
  -D CISM_BUILD_CISM_DRIVER:BOOL=ON \
\
  -D CISM_USE_GPTL_INSTRUMENTATION:BOOL=OFF \
  -D CISM_COUPLED:BOOL=OFF \
  -D CISM_USE_CISM_FRONT_END:BOOL=OFF \
\
  -D CISM_NETCDF_DIR="$NETCDF" \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
\
  -D CMAKE_C_COMPILER=cc \
  -D CMAKE_Fortran_COMPILER=ftn \
\
  -D CMAKE_EXE_LINKER_FLAGS="" \
\
  -D CMAKE_Fortran_FLAGS:STRING="-fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -fallow-argument-mismatch  -fallow-invalid-boz  -g -Wall -O0 -g -fbacktrace -ffpe-trap=zero,overflow -fcheck=bounds " \
  -D CMAKE_C_FLAGS:STRING="-std=gnu99 -g -Wall -O0 -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds " \
  "${cism_top}"

