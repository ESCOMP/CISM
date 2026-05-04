# Run this script by typing: source olivia-gnu-cmake.sh
# After this script completes, type: make
# If rebuilding, type 'make clean' before running 'make'
# or call ./clean.sh

## December 2025
#module purge
#source /opt/cray/pe/lmod/lmod/init/profile
#export MODULEPATH=/cluster/software/modules/Core/
#module load NRIS/CPU
#module load CMake/3.26.3-GCCcore-12.3.0
#module load netCDF-Fortran/4.6.1-gompi-2023a
#module load OpenBLAS/0.3.23-GCC-12.3.0

# source modules
. arch-modules

# remove old build data:
rm -f ./CMakeCache.txt
rm -rf ./CMakeFiles

echo
echo "Doing CMake Configuration step"

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
  -D CISM_NETCDF_DIR=$EBROOTNETCDFMINFORTRAN \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
\
  -D CMAKE_CXX_COMPILER=mpiicpc \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_Fortran_COMPILER=mpif90 \
\
  -D CMAKE_EXE_LINKER_FLAGS="-Wl,-rpath=${EBROOTNETCDFMINFORTRAN}/lib,-rpath=${EBROOTNETCDF}/lib" \
\
  -D CMAKE_Fortran_FLAGS:STRING="-fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -ffree-form -O" \
  -D CMAKE_C_FLAGS:STRING="-std=gnu99 -O" \
  -D CISM_EXTRA_LIBS:STRING="-lopenblas" \
 ../..

# Note: last argument above  "../.."  is path to top-level cism directory

