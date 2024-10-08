Macport environment from 10.2024 (Apple M3)  

For CISM:
port install gcc12
port select --set gcc mp-gcc12
port install netcdf
port install netcdf-fortran
port install openmpi 
port select --set mpi openmpi-mp-fortran

./clean.sh
. mac-gnu-cmake
make

