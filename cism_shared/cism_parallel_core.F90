!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_parallel_core.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_parallel_core

  ! Part of the cism_parallel module:
  !  the parallel_type and the shared variables (communicator, ranks, halo widths);
  !  MPI setup and shutdown; global and local index utilities
  ! Other modules use these through cism_parallel.

  use glimmer_global, only : dp, sp
  use glimmer_paramets, only: iulog

  implicit none

  ! integers associated with the main global communicator
  ! (for communication among all tasks)
  ! Note: When running CISM with multiple ice sheet instances, these values are assumed
  !       to be the same for all instances.
  integer :: comm           ! integer ID for the main global communicator
  integer :: tasks          ! total number of tasks
  integer :: this_rank      ! integer ID for the local task
  integer :: main_rank      ! integer ID for the master task
  logical :: main_task      ! true if this_rank = main_rank

  ! Halo information
  ! Note: The glissade dycore currently requires nhalo = 2, and the glide dycore requires nhalo = 0.
  !       For glide simulations, nhalo is set to 0 by calling distributed_grid with optional argument nhalo = 0.
  !       The halo parameters are outside the parallel derived type, on the assumption
  !        that all instances will use the same dycore.

  integer, save :: nhalo = 2

  !TODO - Define lhalo and uhalo in terms of nhalo?

  integer, save :: lhalo = 2
  integer, save :: uhalo = 2

  ! halo widths for staggered grid
!  integer,parameter :: staggered_lhalo = lhalo
!  integer,parameter :: staggered_uhalo = uhalo-1
  integer, save :: staggered_lhalo = 2
  integer, save :: staggered_uhalo = 1

  ! Information needed to carry out parallel operations (halo, broadcast, gather, scatter, etc.)
  !  for a particular ice sheet instance
  ! Note: At one time, these were independent module-level variables in this module.
  !       This works only for runs with a single ice sheet instance.
  !       With multiple instances, each instance requires its own grid and decomposition.

  type :: parallel_type

     ! global boundary conditions
     logical :: periodic_bc = .true.  ! doubly periodic
     logical :: outflow_bc = .false.  ! if true, set scalars in global halo to zero, and set
                                      ! staggered data to zero beyond the global boundary
     logical :: no_ice_bc = .false.   ! if true, scalars adjacent to the global boundary are set to zero;
                                      ! this includes halo cells and one row of locally owned cells

     ! grid parameters on the distributed grid
     integer :: global_ewn, global_nsn  ! no. of cells in x and y dimensions on global grid
     integer :: local_ewn,  local_nsn   ! no. of cells in x and y dimensions on local grid (including halos)
     integer :: own_ewn,    own_nsn     ! no. of locally owned cells in x and y dimensions (excluding halos) 
     integer :: global_row_offset, global_col_offset  ! no. of row and columns preceding the local rows and columns

     integer :: ewlb, ewub, nslb, nsub    ! indices for lower and upper cells in x and y dimensions
     integer :: east, west, north, south  ! integer ID for neighboring processes in each direction
     integer :: ewtasks, nstasks          ! no. of tasks in x and y directions
     integer :: ewrank, nsrank            ! this processor's rank in x and y directions

     !WHL - added to handle gathers and scatters correctly when computing on active blocks only
     integer :: global_minval_ewlb, global_maxval_ewub
     integer :: global_minval_nslb, global_maxval_nsub

     ! logical variables to identify corner tasks.
     ! A southeast corner task is a task with active south and east neighbors but an
     !  inactive southeast neighbor; and similarly for other directions.
     ! Subroutine distributed_grid_active_blocks typically distributes tasks
     !  in a way that there are corner tasks.
     logical :: southwest_corner   ! true if a task has inactive SW neighbor, active S and W 
     logical :: southeast_corner   ! true if a task has inactive SE neighbor, active S and E
     logical :: northeast_corner   ! true if a task has inactive NE neighbor, active N and E
     logical :: northwest_corner   ! true if a task has inactive NW neighbor, active N and W

     ! bounds of locally owned vertices on staggered grid
     ! For periodic BC, staggered_ilo = staggered_jlo = staggered_lhalo+1.
     ! For outflow BC the locally owned vertices include the southern and western rows
     !  of the global domain, so staggered_ilo = staggered_jlo = staggered_lhalo on
     !  processors that include these rows.
     ! TODO: Carrying around these variables adds many lines of code.
     !       Would it be possible to support outflow BCs without these variables?
     !       Might want to revisit whether the staggered grid is needed.
     integer :: staggered_ilo
     integer :: staggered_jlo
     integer :: staggered_ihi
     integer :: staggered_jhi

     ! optional row-based and column-based communicators
     ! These can be used for global tridiagonal solves

     ! integers associated with the row-based communicator
     ! (for communication among tasks with the same value of nsrank)
     integer :: comm_row       ! integer ID for the row-based communicator
     integer :: tasks_row      ! total number of tasks on the local row
     integer :: this_rank_row  ! integer ID for the local task in the row
     integer :: main_rank_row  ! integer ID for the master task on the row
     logical :: main_task_row  ! true if this_rank_row = main_rank_row

     ! integers associated with the column-based communicator
     ! (for communication among tasks with the same value of ewrank)
     integer :: comm_col       ! integer ID for the column-based communicator
     integer :: tasks_col      ! total number of tasks on the local column
     integer :: this_rank_col  ! integer ID for the local task in the column
     integer :: main_rank_col  ! integer ID for the master task on the column
     logical :: main_task_col  ! true if this_rank_col = main_rank_col

     ! option to compute reproducible sums
     logical :: reprosum       ! if true, compute reproducible global sums
                               ! (interface parallel_reduce_sum)

  end type parallel_type

  interface parallel_boundary_value
     module procedure parallel_boundary_value_real8_2d
     module procedure parallel_boundary_value_real8_3d
  end interface parallel_boundary_value

  interface parallel_convert_haloed_to_nonhaloed
     module procedure parallel_convert_haloed_to_nonhaloed_real4_2d
     module procedure parallel_convert_haloed_to_nonhaloed_real8_2d
  end interface parallel_convert_haloed_to_nonhaloed

contains

!=======================================================================

  function distributed_execution()
     ! Returns if running distributed or not.
     logical distributed_execution

     distributed_execution = .true.
  end function distributed_execution

!=======================================================================

  function distributed_isparallel()

     implicit none
     logical :: distributed_isparallel

     distributed_isparallel = .true.

  end function distributed_isparallel

!=======================================================================

  function distributed_owner(ew, ewn, ns, nsn, parallel)

    implicit none
    logical :: distributed_owner
    integer :: ew,ewn,ns,nsn
    type(parallel_type) :: parallel

    associate(  &
         local_ewn  => parallel%local_ewn,    &
         local_nsn  => parallel%local_nsn     &
         )

    ! begin
    distributed_owner = (ew > lhalo .and. ew <= local_ewn-uhalo .and.&
         ns > lhalo .and. ns <= local_nsn-uhalo)

    end associate

  end function distributed_owner

!=======================================================================

  subroutine not_parallel(file, line)

    implicit none
    integer :: line
    character(len=*) :: file

    ! begin
    call parallel_stop(file,line)

  end subroutine not_parallel

!=======================================================================

  subroutine parallel_barrier

    use mpi_mod
    implicit none

    integer :: ierror

    ! begin
    call mpi_barrier(comm,ierror)

  end subroutine parallel_barrier

!=======================================================================

  function parallel_boundary(ew, ns, parallel)

    implicit none
    logical :: parallel_boundary
    integer :: ew,ewn,ns,nsn
    type(parallel_type) :: parallel

    ! begin
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         global_ewn  => parallel%global_ewn,   &
         global_nsn  => parallel%global_nsn,   &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub          &
         )

    parallel_boundary = (ewlb<1 .and. ew==1+lhalo) .or. &
               (ewub>global_ewn .and. ew==local_ewn-uhalo) .or. &
                        (nslb<1 .and. ns==1+lhalo) .or. &
               (nsub>global_nsn .and. ns==local_nsn-uhalo)

    end associate

  end function parallel_boundary

!=======================================================================

  ! subroutines belonging to the parallel_boundary_value interface

  subroutine parallel_boundary_value_real8_2d(&
       field,          &
       boundary_value, &
       parallel)

    ! Insert a specified value into cells on the global boundary.
    ! Typically called with value = 0.0 or an special value.

    real(dp), dimension(:,:), intent(inout) :: field
    real(dp), intent(in) :: boundary_value
    type(parallel_type) :: parallel

    integer :: ew, ns

    ! begin
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn     &
         )

    do ns = 1, local_nsn
       do ew = 1, local_ewn
          if (parallel_boundary(ew, ns, parallel)) then
             field(ew,ns) = boundary_value
          endif
       enddo
    enddo

    end associate

  end subroutine parallel_boundary_value_real8_2d


  subroutine parallel_boundary_value_real8_3d(&
       field,          &
       boundary_value, &
       parallel)

    ! Insert a specified value into cells on the global boundary.
    ! Typically called with value = 0.0 or an special value.

    real(dp), dimension(:,:,:), intent(inout) :: field
    real(dp), intent(in) :: boundary_value
    type(parallel_type) :: parallel

    integer :: ew, ns

    ! begin
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn     &
         )

    do ns = 1, local_nsn
       do ew = 1, local_ewn
          if (parallel_boundary(ew, ns, parallel)) then
             field(:,ew,ns) = boundary_value
          endif
       enddo
    enddo

    end associate

  end subroutine parallel_boundary_value_real8_3d

!=======================================================================

  ! subroutines belonging to the parallel_convert_haloed_to_nonhaloed interface

  subroutine parallel_convert_haloed_to_nonhaloed_real4_2d(input_with_halo, output_no_halo, parallel)

    ! Given an input array that has halo cells, return an output array without halo cells
    real(sp),dimension(:,:), intent(in)  :: input_with_halo
    real(sp),dimension(:,:), intent(out) :: output_no_halo
    type(parallel_type) :: parallel

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn       &
         )

    if (size(input_with_halo,1) /= local_ewn .or. size(input_with_halo,2) /= local_nsn) then
       write(iulog,*) "Unexpected size for input_with_halo: ", &
            size(input_with_halo,1), size(input_with_halo,2)
       write(iulog,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_no_halo,1) /= own_ewn .or. size(output_no_halo,2) /= own_nsn) then
       write(iulog,*) "Unexpected size for output_no_halo: ", &
            size(output_no_halo,1), size(output_no_halo,2)
       write(iulog,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_no_halo(1:own_ewn, 1:own_nsn) = &
         input_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

    end associate

  end subroutine parallel_convert_haloed_to_nonhaloed_real4_2d


  subroutine parallel_convert_haloed_to_nonhaloed_real8_2d(input_with_halo, output_no_halo, parallel)

    ! Given an input array that has halo cells, return an output array without halo cells
    real(dp),dimension(:,:), intent(in)  :: input_with_halo
    real(dp),dimension(:,:), intent(out) :: output_no_halo
    type(parallel_type) :: parallel

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn       &
         )

    if (size(input_with_halo,1) /= local_ewn .or. size(input_with_halo,2) /= local_nsn) then
       write(iulog,*) "Unexpected size for input_with_halo: ", &
            size(input_with_halo,1), size(input_with_halo,2)
       write(iulog,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_no_halo,1) /= own_ewn .or. size(output_no_halo,2) /= own_nsn) then
       write(iulog,*) "Unexpected size for output_no_halo: ", &
            size(output_no_halo,1), size(output_no_halo,2)
       write(iulog,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_no_halo(1:own_ewn, 1:own_nsn) = &
         input_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

    end associate

  end subroutine parallel_convert_haloed_to_nonhaloed_real8_2d

!=======================================================================

  subroutine parallel_finalise

    use mpi_mod
    implicit none
    integer :: ierror

    ! begin
    call mpi_finalize(ierror)

  end subroutine parallel_finalise

!=======================================================================

  !TODO - Is function parallel_globalID still needed?  No longer called except from glissade_test_halo.

  function parallel_globalID(locns, locew, upstride, parallel)

    ! Returns a unique ID for a given row and column reference that is identical across all processors.
    ! For instance if Proc 0: (17,16) is the same global cell as Proc 3: (17,1), then the globalID will be the same for both.
    ! These IDs are spaced upstride apart.  upstride = number of vertical layers.
    ! Typically (upn) + number of ghost layers (2 = top and bottom)

    integer,intent(in) :: locns, locew, upstride
    integer :: parallel_globalID
    type(parallel_type) :: parallel

    ! locns is local NS (row) grid index
    ! locew is local EW (col) grid index
    integer :: global_row, global_col, global_ID
    character(len=40) :: local_coord

    associate(  &
         periodic_bc       => parallel%periodic_bc,        &
         global_ewn        => parallel%global_ewn,         &
         global_nsn        => parallel%global_nsn,         &
         global_row_offset => parallel%global_row_offset,  &
         global_col_offset => parallel%global_col_offset   &
         )

    ! including global domain halo adds lhalo to offsets
    global_row = (locns - lhalo) + (global_row_offset + lhalo)
    global_col = (locew - lhalo) + (global_col_offset + lhalo)

    ! if halo cell and if using periodic boundary conditions,
    ! define global ID to be associated non-halo cell
    !WHL - Commented out the deprecated options horiz_bcs_type_south/north/west/east.
    !      Replaced with periodic_bc logic, but not tested.

    if (global_row <= lhalo) then
!       if (horiz_bcs_type_south == HORIZ_BCS_CYCLIC) then
       if (periodic_bc) then
          global_row = global_row + global_nsn
       endif
    endif

    if (global_row > (global_nsn+lhalo)) then
!       if  (horiz_bcs_type_north == HORIZ_BCS_CYCLIC) then
       if (periodic_bc) then
          global_row = global_row - global_nsn
       endif
    endif

    if (global_col .le. lhalo) then
!       if (horiz_bcs_type_west == HORIZ_BCS_CYCLIC) then
       if (periodic_bc) then
          global_col = global_col + global_ewn
       endif
    endif

    if (global_col > (global_ewn+lhalo)) then
!       if (horiz_bcs_type_east == HORIZ_BCS_CYCLIC) then
       if (periodic_bc) then
          global_col = global_col - global_ewn
       endif
    endif

    ! including global domain halo adds (lhalo + uhalo) to global_ewn
    global_ID = ((global_row - 1) * (global_ewn + lhalo + uhalo) + (global_col - 1)) * upstride + 1

    ! Testing Code
    ! write(local_coord, "A13,I10.1,A2,I10.1,A1") " (NS, EW) = (", locns, ", ", locew, ")"
    ! write(iulog,*) "Processor reference ", this_rank, local_coord, " globalID = ", global_ID

    !return value
    parallel_globalID = global_ID

    end associate

  end function parallel_globalID

!=======================================================================

  function parallel_globalID_scalar(locew, locns, upstride, parallel)

    !WHL - This function is similar to parallel_globalID, but assigns 0's to cells outside the global domain

    ! Returns a unique ID for a given row and column reference that is identical across all processors.
    ! For instance if Proc 0: (17,16) is the same global cell as Proc 3: (17,1), then the globalID will be the same for both.
    ! These IDs are spaced upstride apart.  upstride = number of vertical layers.
    integer,intent(in) :: locns, locew, upstride
    integer :: parallel_globalID_scalar
    type(parallel_type) :: parallel

    ! locns is local NS (row) grid index
    ! locew is local EW (col) grid index
    integer :: global_row, global_col, global_ID
    character(len=40) :: local_coord

    associate(  &
         global_ewn        => parallel%global_ewn,         &
         global_nsn        => parallel%global_nsn,         &
         global_row_offset => parallel%global_row_offset,  &
         global_col_offset => parallel%global_col_offset   &
         )

    ! including global domain halo adds lhalo to offsets
    global_row = (locns - lhalo) + global_row_offset
    global_col = (locew - lhalo) + global_col_offset

    ! including global domain halo adds (lhalo + uhalo) to global_ewn
    global_ID = ((global_row - 1)*(global_ewn) + (global_col - 1)) * upstride + 1

    ! Testing Code
    ! write(local_coord, "A13,I10.1,A2,I10.1,A1") " (NS, EW) = (", locns, ", ", locew, ")"
    ! write(iulog,*) "Processor reference ", this_rank, local_coord, " globalID = ", global_ID

    !return value
    parallel_globalID_scalar = global_ID

    end associate

  end function parallel_globalID_scalar

!=======================================================================

  subroutine parallel_globalindex(ilocal, jlocal, iglobal, jglobal, parallel)

    ! Calculates the global i,j indices from the local i,j indices
    integer,intent(in)  :: ilocal,  jlocal  ! These include the halos
    integer,intent(out) :: iglobal, jglobal ! These do NOT include halos
    type(parallel_type) :: parallel

    associate(  &
         global_row_offset => parallel%global_row_offset,  &
         global_col_offset => parallel%global_col_offset   &
         )

    ! Note: if the local index is in a halo, still convert that to its location
    ! on the global grid (even though that location on the global grid is owned
    ! by a different processor!)
    ! No check is currently made for being located in the global (periodic) halo
    iglobal = (ilocal - lhalo) + global_col_offset
    jglobal = (jlocal - lhalo) + global_row_offset

    end associate

  end subroutine parallel_globalindex

!=======================================================================

  ! Note: parallel_initialise should generally be called only by standalone cism drivers
  ! When cism is nested inside a climate model (so mpi_init has already been called) use parallel_set_info instead

  subroutine parallel_initialise

    use mpi_mod 
    implicit none

    integer :: ierror 
    integer, parameter :: my_main_rank = 0

    ! begin 
    call mpi_init(ierror)
    call parallel_set_info(mpi_comm_world, my_main_rank)

  end subroutine parallel_initialise

!=======================================================================

  ! Note: parallel_set_info should be called directly when cism is nested inside a climate model
  !       (Since mpi_init has already been called, do NOT use parallel_initialise)

  subroutine parallel_set_info(my_comm, my_main_rank)

    use mpi_mod
    implicit none
    integer, intent(in) :: my_comm       ! CISM's global communicator
    integer, intent(in) :: my_main_rank  ! rank of the master task

    integer :: ierror 

    ! begin
    comm = my_comm
    main_rank = my_main_rank
    call mpi_comm_size(comm, tasks, ierror)
    call mpi_comm_rank(comm, this_rank, ierror)
    main_task = (this_rank==main_rank)

  end subroutine parallel_set_info

!=======================================================================

  subroutine parallel_stop(file, line)

    use mpi_mod
    implicit none
    integer :: line
    character(len=*) :: file

    integer :: ierror

    ! begin
    if (main_task) write(0,*) "PARALLEL STOP in ",file," at line ",line
    call mpi_abort(MPI_COMM_WORLD, 1001, ierror)
    stop "PARALLEL STOP"

  end subroutine parallel_stop

!=======================================================================

end module cism_parallel_core

!=======================================================================
