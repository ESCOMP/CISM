!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_parallel.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2018
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module cism_parallel

  ! Parallel communication (MPI) for CISM.
  ! The code is in the modules below, grouped by theme. Other modules use these names
  !  through cism_parallel (all names in these modules are public):
  !   cism_parallel_core            parallel_type, shared variables (communicator, ranks, halo widths),
  !                                 MPI setup and shutdown, global and local index utilities
  !   cism_parallel_broadcast       broadcasts from the main task
  !   cism_parallel_reduce          global sums and reductions
  !   cism_parallel_halo            halo updates
  !   cism_parallel_gather_scatter  gathers, scatters, distributed reads and writes, printing
  !   cism_parallel_netcdf          netCDF wrappers (I/O on the main task)
  !   cism_parallel_decomp          domain decomposition, row and column communicators

  use cism_parallel_core
  use cism_parallel_broadcast
  use cism_parallel_reduce
  use cism_parallel_halo
  use cism_parallel_gather_scatter
  use cism_parallel_netcdf
  use cism_parallel_decomp

  implicit none

end module cism_parallel

!=======================================================================
