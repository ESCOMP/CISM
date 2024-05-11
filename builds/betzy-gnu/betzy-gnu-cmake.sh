# Run this script by typing: source cheyenne-intel-cmake
# After this script completes, type: make -j 8
# If rebuilding, type 'make clean' before running 'make -j 8'

# This cmake configuration script is set up to perform a parallel build with Trilinos

#module reset
#module load ncarenv/1.2
#module load gnu/8.3.0
#module load openblas/0.3.6
#module load openmpi/3.1.4
#module load netcdf/4.7.1
#module load ncarcompilers/0.5.0
#module load cmake/3.14.4
#module load python/2.7.13
#module load numpy/1.12.0
#module load netcdf4-python/1.2.7

module purge
module load StdEnv
module load GCCcore/11.2.0
module load CMake/3.22.1-GCCcore-11.2.0
module load netCDF/4.8.1-gompi-2021b
module load netCDF-Fortran/4.5.3-gompi-2021b 
module load OpenBLAS/0.3.18-GCC-11.2.0

# remove old build data:
rm -f ./CMakeCache.txt
rm -rf ./CMakeFiles

echo
echo "Doing CMake Configuration step"

# Note: the compilation flags were taken from the defaults for a CESM build on
# cheyenne-intel (using cime at 84aafd5). Some of these options are probably
# unnecessary for a standalone cism build, but I am keeping things consistent
# with the CESM build for simplicity.
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
  -D CMAKE_EXE_LINKER_FLAGS="-ldl" \
\
  -D CMAKE_Fortran_FLAGS:STRING="-fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -ffree-form -O" \
  -D CMAKE_C_FLAGS:STRING="-std=gnu99 -O" \
  -D CISM_EXTRA_LIBS:STRING="-lopenblas" \
 ../..

# Note: last argument above  "../.."  is path to top-level cism directory

#  -D CMAKE_C_FLAGS:STRING="-qno-opt-dynamic-align -fp-model precise -std=gnu99 -qopt-report -O2 -debug minimal " \
#  -D CMAKE_CXX_FLAGS:STRING="-qno-opt-dynamic-align -fp-model precise -std=gnu99 -qopt-report -O2 -debug minimal " \#  -D CMAKE_Fortran_FLAGS:STRING="-qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source -qopt-report -O2 -debug minimal " \


#  -D CMAKE_Fortran_FLAGS:STRING="-fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -ffree-form -O -L/cluster/software/OpenBLAS/0.3.18-GCC-11.2.0/lib/libopenblas.a" \
