!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   parallel_mpi.F90 - part of the Community Ice Sheet Model (CISM)  
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

  use netcdf
  use glimmer_global, only : dp, sp

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

  ! Debug and Verification Level
  integer, parameter :: DEBUG_LEVEL = 1
	! If > 0, then debug code executed.  Added for parallel_halo_verify()

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

  ! common work space
  integer, dimension(4), save :: d_gs_mybounds
  integer, dimension(:,:), allocatable,save :: d_gs_bounds

  ! distributed gather flow control parameter
  integer,parameter :: max_gather_block_size = 64 ! max and default

  ! Information needed to carry out parallel operations (halo, broadcast, gather, scatter, etc.)
  !  for a particular ice sheet instance
  ! Note: At one time, these were independent module-level variables in parallel_mpi.F90.
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

  end type parallel_type

  ! Information on the local & global bounds of an array
  ! This is used to distinguish between arrays on the staggered vs. unstaggered grids
  type, private :: bounds_info_type
     ! Global number of points in each dimension
     integer :: global_ewn
     integer :: global_nsn

     ! Range of indices that this proc is responsible for (excludes halo cells)
     ! These are the indices in global index space
     integer :: mybounds_ew_lb
     integer :: mybounds_ew_ub
     integer :: mybounds_ns_lb
     integer :: mybounds_ns_ub

     ! Local indices that this proc is responsible for (excludes halo cells)
     ! These are the indices in local index space
     integer :: ilo
     integer :: ihi
     integer :: jlo
     integer :: jhi
  end type bounds_info_type

  ! Parallel interfaces (listed alphabetically for the most part)

  interface broadcast
     module procedure broadcast_character
     module procedure broadcast_integer
     module procedure broadcast_integer_1d
     module procedure broadcast_logical
     module procedure broadcast_logical_1d
     module procedure broadcast_real4
     module procedure broadcast_real4_1d
     module procedure broadcast_real8     
     module procedure broadcast_real8_1d
  end interface

  interface distributed_gather_var
     module procedure distributed_gather_var_integer_2d
     module procedure distributed_gather_var_logical_2d
     module procedure distributed_gather_var_real4_2d
     module procedure distributed_gather_var_real4_3d
     module procedure distributed_gather_var_real8_2d
     module procedure distributed_gather_var_real8_3d
  end interface

  interface distributed_gather_var_row
     module procedure distributed_gather_var_row_real8_2d
  end interface

  interface distributed_gather_all_var_row
     module procedure distributed_gather_all_var_row_real8_2d
  end interface

  interface distributed_gather_var_col
     module procedure distributed_gather_var_col_real8_2d
  end interface

  interface distributed_gather_all_var_col
     module procedure distributed_gather_all_var_col_real8_2d
  end interface

  interface distributed_get_var
     module procedure distributed_get_var_integer_2d
     module procedure distributed_get_var_real4_1d
     module procedure distributed_get_var_real4_2d
     module procedure distributed_get_var_real8_1d
     module procedure distributed_get_var_real8_2d
     module procedure distributed_get_var_real8_3d
  end interface

  interface distributed_print
     ! Gathers a distributed variable and writes to file
     module procedure distributed_print_integer_2d
     module procedure distributed_print_real8_2d
     module procedure distributed_print_real8_3d
  end interface

  interface distributed_put_var
     module procedure distributed_put_var_integer_2d
     module procedure distributed_put_var_real4_1d
     module procedure distributed_put_var_real4_2d
     module procedure distributed_put_var_real8_1d
     module procedure distributed_put_var_real8_2d
     module procedure distributed_put_var_real8_3d
  end interface

  interface distributed_scatter_var
     module procedure distributed_scatter_var_integer_2d
     module procedure distributed_scatter_var_logical_2d
     module procedure distributed_scatter_var_real4_2d
     module procedure distributed_scatter_var_real4_3d
     module procedure distributed_scatter_var_real8_2d
     module procedure distributed_scatter_var_real8_3d
  end interface

  interface distributed_scatter_var_row
     module procedure distributed_scatter_var_row_real8_2d
  end interface

  interface distributed_scatter_var_col
     module procedure distributed_scatter_var_col_real8_2d
  end interface

  interface parallel_boundary_value
     module procedure parallel_boundary_value_real8_2d
     module procedure parallel_boundary_value_real8_3d
  end interface parallel_boundary_value

  interface parallel_convert_haloed_to_nonhaloed
     module procedure parallel_convert_haloed_to_nonhaloed_real4_2d
     module procedure parallel_convert_haloed_to_nonhaloed_real8_2d
  end interface parallel_convert_haloed_to_nonhaloed

  interface parallel_convert_nonhaloed_to_haloed
     module procedure parallel_convert_nonhaloed_to_haloed_real4_2d
     module procedure parallel_convert_nonhaloed_to_haloed_real8_2d
  end interface parallel_convert_nonhaloed_to_haloed
  
  interface parallel_def_var
     module procedure parallel_def_var_dimids
     module procedure parallel_def_var_nodimids
  end interface

  interface parallel_get_att
     module procedure parallel_get_att_character
     module procedure parallel_get_att_real4
     module procedure parallel_get_att_real4_1d
     module procedure parallel_get_att_real8
     module procedure parallel_get_att_real8_1d
  end interface

  interface parallel_get_var
     module procedure parallel_get_var_integer
     module procedure parallel_get_var_real4
     module procedure parallel_get_var_real8
     module procedure parallel_get_var_integer_1d
     module procedure parallel_get_var_real4_1d
     module procedure parallel_get_var_real8_1d
     module procedure parallel_get_var_integer_2d
     module procedure parallel_get_var_real8_2d
  end interface

  interface parallel_global_sum
     module procedure parallel_global_sum_integer_2d
     module procedure parallel_global_sum_real4_2d
     module procedure parallel_global_sum_real8_2d
  end interface

  interface parallel_halo
     module procedure parallel_halo_integer_2d
     module procedure parallel_halo_logical_2d
     module procedure parallel_halo_real4_2d
     module procedure parallel_halo_real8_2d
     module procedure parallel_halo_real8_3d
     module procedure parallel_halo_real8_4d
  end interface

  interface parallel_halo_extrapolate
     module procedure parallel_halo_extrapolate_integer_2d
     module procedure parallel_halo_extrapolate_real8_2d
  end interface

  interface parallel_halo_tracers
     module procedure parallel_halo_tracers_real8_3d
     module procedure parallel_halo_tracers_real8_4d
  end interface

  interface parallel_halo_verify
     module procedure parallel_halo_verify_integer_2d
     module procedure parallel_halo_verify_real8_2d
     module procedure parallel_halo_verify_real8_3d
  end interface

  interface parallel_print
     module procedure parallel_print_integer_2d
     module procedure parallel_print_real8_2d
     module procedure parallel_print_real8_3d
  end interface

  interface parallel_put_att
     module procedure parallel_put_att_character
     module procedure parallel_put_att_integer
     module procedure parallel_put_att_real4
     module procedure parallel_put_att_real4_1d
     module procedure parallel_put_att_real8
     module procedure parallel_put_att_real8_1d
  end interface

  interface parallel_put_var
     module procedure parallel_put_var_integer
     module procedure parallel_put_var_integer_1d
     module procedure parallel_put_var_real4
     module procedure parallel_put_var_real8
     module procedure parallel_put_var_real8_1d
  end interface

  interface parallel_reduce_max
     module procedure parallel_reduce_max_integer
     module procedure parallel_reduce_max_real4
     module procedure parallel_reduce_max_real8
     module procedure parallel_reduce_max_real8_1d
  end interface

  ! This reduce interface determines the global max value and the processor on which it occurs
  interface parallel_reduce_maxloc
     module procedure parallel_reduce_maxloc_integer
     module procedure parallel_reduce_maxloc_real4
     module procedure parallel_reduce_maxloc_real8
  end interface

  interface parallel_reduce_min
     module procedure parallel_reduce_min_integer
     module procedure parallel_reduce_min_real4
     module procedure parallel_reduce_min_real8
     module procedure parallel_reduce_min_real8_1d
  end interface

  ! This reduce interface determines the global min value and the processor on which it occurs
  interface parallel_reduce_minloc
     module procedure parallel_reduce_minloc_integer
     module procedure parallel_reduce_minloc_real4
     module procedure parallel_reduce_minloc_real8
  end interface

  interface parallel_reduce_sum
     module procedure parallel_reduce_sum_integer
     module procedure parallel_reduce_sum_real4
     module procedure parallel_reduce_sum_real8
     module procedure parallel_reduce_sum_integer_nvar
     module procedure parallel_reduce_sum_real8_nvar
  end interface

  interface staggered_parallel_halo
     module procedure staggered_parallel_halo_integer_2d
     module procedure staggered_parallel_halo_integer_3d
     module procedure staggered_parallel_halo_real8_2d
     module procedure staggered_parallel_halo_real8_3d
     module procedure staggered_parallel_halo_real8_4d
  end interface

  interface staggered_parallel_halo_extrapolate
     module procedure staggered_parallel_halo_extrapolate_integer_2d
     module procedure staggered_parallel_halo_extrapolate_real8_2d
  end interface

contains

!=======================================================================

  ! subroutines belonging to the broadcast interface

  subroutine broadcast_character(c, proc)

    use mpi_mod
    implicit none
    character(len=*) :: c
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror, n
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    n = len(c)
    call mpi_bcast(c,n,mpi_character,source,comm,ierror)

  end subroutine broadcast_character


  subroutine broadcast_integer(i, proc)

    use mpi_mod
    implicit none
    integer :: i
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(i,1,mpi_integer,source,comm,ierror)

  end subroutine broadcast_integer


  subroutine broadcast_integer_1d(a, proc)

    use mpi_mod
    implicit none
    integer,dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(a,size(a),mpi_integer,source,comm,ierror)

  end subroutine broadcast_integer_1d


  subroutine broadcast_logical(l, proc)

    use mpi_mod
    implicit none
    logical :: l
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(l,1,mpi_logical,source,comm,ierror)

  end subroutine broadcast_logical


  subroutine broadcast_logical_1d(l, proc)

    use mpi_mod
    implicit none
    logical,dimension(:) :: l
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(l,size(l),mpi_logical,source,comm,ierror)

  end subroutine broadcast_logical_1d


  subroutine broadcast_real4(r, proc)

    use mpi_mod
    implicit none
    real(sp) :: r
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(r,1,mpi_real4,source,comm,ierror)

  end subroutine broadcast_real4


  subroutine broadcast_real4_1d(a, proc)

    use mpi_mod
    implicit none
    real(sp),dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(a,size(a),mpi_real4,source,comm,ierror)

  end subroutine broadcast_real4_1d


  subroutine broadcast_real8(r, proc)

    use mpi_mod
    implicit none
    real(dp) :: r
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(r,1,mpi_real8,source,comm,ierror)

  end subroutine broadcast_real8


  subroutine broadcast_real8_1d(a, proc)

    use mpi_mod
    implicit none
    real(dp),dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(a,size(a),mpi_real8,source,comm,ierror)

  end subroutine broadcast_real8_1d

!=======================================================================

  function distributed_execution()
     ! Returns if running distributed or not.
     logical distributed_execution

     distributed_execution = .true.
  end function distributed_execution

!=======================================================================

  ! subroutines belonging to the distributed_gather_var interface

  ! WHL, July 2019:
  ! There is an issue with allocating the global_values array in the distributed_gather_var_*,
  !  distributed_get_var_*, distributed_print_*, and distributed_put_var_* functions and subroutines
  !  when computing only on active blocks (compute_blocks = 1).
  ! This array is allocated based on the max and min of ewlb, ewub, nslb, and nsub over the global domain.
  ! Previously, this was done based on a bounds array computed in fc_gather_int, which gathers
  !  the bounds for all tasks (i.e., all active blocks).
  ! However, these global bounds will be incorrect if either the west, east, south, or north row
  !  of the global domain contains only inactive blocks.
  ! The fix is to allocate global values based on global_minval_ewlb, global_maxval_ewub,
  !  global_minval_nslb, and global_maxval_nsub, which are now computed at initialization
  !  based on the bounds in all blocks (including inactive blocks), not just active blocks.

  subroutine distributed_gather_var_integer_2d(values, global_values, parallel)

    ! Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    integer,dimension(:,:),intent(in) :: values
    integer,dimension(:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:),allocatable :: recvbuf
    integer,dimension(:,:),allocatable :: sendbuf

    associate(  &
         ewlb      => parallel%ewlb,       &
         ewub      => parallel%ewub,       &
         nslb      => parallel%nslb,       &
         nsub      => parallel%nsub,       &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         global_minval_ewlb => parallel%global_minval_ewlb,  &
         global_maxval_ewub => parallel%global_maxval_ewub,  &
         global_minval_nslb => parallel%global_minval_nslb,  &
         global_maxval_nsub => parallel%global_maxval_nsub   &
    )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_gather does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(&
!!                 minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
!!                 minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       allocate(global_values(&
                 global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
                 global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1) &
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if

    allocate(sendbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_int(sendbuf,size(sendbuf),mpi_integer,&
       recvbuf,recvcounts,displs,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_var_integer_2d


  subroutine distributed_gather_var_logical_2d(values, global_values, parallel)

    ! Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    logical,dimension(:,:),intent(in) :: values
    logical,dimension(:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    logical,dimension(:),allocatable :: recvbuf
    logical,dimension(:,:),allocatable :: sendbuf

    associate(  &
         ewlb      => parallel%ewlb,       &
         ewub      => parallel%ewub,       &
         nslb      => parallel%nslb,       &
         nsub      => parallel%nsub,       &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         global_minval_ewlb => parallel%global_minval_ewlb,  &
         global_maxval_ewub => parallel%global_maxval_ewub,  &
         global_minval_nslb => parallel%global_minval_nslb,  &
         global_maxval_nsub => parallel%global_maxval_nsub   &
    )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_gather does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(&
!!                 minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
!!                 minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       allocate(global_values(&
                 global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
                 global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:) = .false.
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_log(sendbuf,size(sendbuf),mpi_logical,&
         recvbuf,recvcounts,displs,mpi_logical,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_var_logical_2d


  subroutine distributed_gather_var_real4_2d(values, global_values, parallel)

    ! Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(sp),dimension(:,:),intent(in) :: values
    real(sp),dimension(:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    real(sp),dimension(:),allocatable :: recvbuf
    real(sp),dimension(:,:),allocatable :: sendbuf

    associate(  &
         ewlb      => parallel%ewlb,       &
         ewub      => parallel%ewub,       &
         nslb      => parallel%nslb,       &
         nsub      => parallel%nsub,       &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         global_minval_ewlb => parallel%global_minval_ewlb,  &
         global_maxval_ewub => parallel%global_maxval_ewub,  &
         global_minval_nslb => parallel%global_minval_nslb,  &
         global_maxval_nsub => parallel%global_maxval_nsub   &
    )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_gather does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(&
!!                 minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
!!                 minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       allocate(global_values(&
                 global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
                 global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:) = 0.0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1) &
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real4(sendbuf,size(sendbuf),mpi_real4,&
       recvbuf,recvcounts,displs,mpi_real4,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_var_real4_2d


  subroutine distributed_gather_var_real4_3d(values, global_values, parallel, ld1, ud1)

    ! Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(sp),dimension(:,:,:),intent(in) :: values
    real(sp),dimension(:,:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel
    integer,optional,intent(in) :: ld1, ud1

    integer :: i,ierror,j,k,d1l,d1u
    integer,dimension(:),allocatable :: displs,recvcounts
    real(sp),dimension(:),allocatable :: recvbuf
    real(sp),dimension(:,:,:),allocatable :: sendbuf

    associate(  &
         ewlb      => parallel%ewlb,       &
         ewub      => parallel%ewub,       &
         nslb      => parallel%nslb,       &
         nsub      => parallel%nsub,       &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         global_minval_ewlb => parallel%global_minval_ewlb,  &
         global_maxval_ewub => parallel%global_maxval_ewub,  &
         global_minval_nslb => parallel%global_minval_nslb,  &
         global_maxval_nsub => parallel%global_maxval_nsub   &
    )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_gather does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       if (present(ld1)) then
         d1l = ld1
       else
         d1l = 1
       endif
       if (present(ud1)) then
         d1u = ud1
       else
         d1u = size(values,1)-(d1l-1)
       endif
       if (size(values,1) /= d1u-d1l+1) then
          write(*,*) "size(values,1) .ne. d1u-d1l+1 in gather call"
          call parallel_stop(__FILE__, __LINE__)
       endif
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(d1l:d1u,&
!!                              minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
!!                              minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       allocate(global_values(d1l:d1u,&
                global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
                global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:,:) = 0.0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)&
                      *size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(size(values,1),&
                     d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:,:) = values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real4(sendbuf,size(sendbuf),mpi_real4,&
       recvbuf,recvcounts,displs,mpi_real4,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(:,&
                        d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/size(values,1),&
                       d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_var_real4_3d


  subroutine distributed_gather_var_real8_2d(values, global_values, parallel)

    ! Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:),intent(in) :: values
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:),allocatable :: sendbuf

    associate(  &
         ewlb      => parallel%ewlb,       &
         ewub      => parallel%ewub,       &
         nslb      => parallel%nslb,       &
         nsub      => parallel%nsub,       &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         global_minval_ewlb => parallel%global_minval_ewlb,  &
         global_maxval_ewub => parallel%global_maxval_ewub,  &
         global_minval_nslb => parallel%global_minval_nslb,  &
         global_maxval_nsub => parallel%global_maxval_nsub   &
    )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_gather does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(&
!!                 minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
!!                 minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       allocate(global_values(&
                 global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
                 global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:) = 0.0d0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if

    allocate(sendbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_var_real8_2d


  subroutine distributed_gather_var_real8_3d(values, global_values, parallel, ld1, ud1)

    ! Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:),intent(in) :: values
    real(dp),dimension(:,:,:),allocatable,intent(inout) :: global_values
    integer,optional,intent(in) :: ld1, ud1
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k,d1l,d1u
    integer,dimension(:),allocatable :: displs,recvcounts
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:,:),allocatable :: sendbuf

    associate(  &
         ewlb      => parallel%ewlb,       &
         ewub      => parallel%ewub,       &
         nslb      => parallel%nslb,       &
         nsub      => parallel%nsub,       &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         global_minval_ewlb => parallel%global_minval_ewlb,  &
         global_maxval_ewub => parallel%global_maxval_ewub,  &
         global_minval_nslb => parallel%global_minval_nslb,  &
         global_maxval_nsub => parallel%global_maxval_nsub   &
    )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_gather does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       if (present(ld1)) then
         d1l = ld1
       else
         d1l = 1
       endif
       if (present(ud1)) then
         d1u = ud1
       else
         d1u = size(values,1)-(d1l-1)
       endif
       if (size(values,1) /= d1u-d1l+1) then
          write(*,*) "size(values,1) .ne. d1u-d1l+1 in gather call"
          call parallel_stop(__FILE__, __LINE__)
       endif
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(d1l:d1u,&
!!                              minval(d_gs_bounds(1,:)):maxval(d_gs_bounds(2,:)),&
!!                              minval(d_gs_bounds(3,:)):maxval(d_gs_bounds(4,:))))
       allocate(global_values(d1l:d1u,&
                global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
                global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:,:) = 0.0d0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)&
                      *size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       if (allocated(global_values)) then
          deallocate(global_values)
       endif
       allocate(global_values(1,1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(size(values,1),&
                     d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    sendbuf(:,:,:) = values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(:,&
                        d_gs_bounds(1,i):d_gs_bounds(2,i),&
                        d_gs_bounds(3,i):d_gs_bounds(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/size(values,1),&
                       d_gs_bounds(2,i)-d_gs_bounds(1,i)+1,&
                       d_gs_bounds(4,i)-d_gs_bounds(3,i)+1/))
       end do
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_var_real8_3d

!=======================================================================

  ! subroutines belonging to the distributed_gather_var_row interface

  subroutine distributed_gather_var_row_real8_2d(values, global_values, parallel)

    ! Gather data along a row of tasks onto the main task for that row.
    ! Based on distributed_gather_var_real8_2d.
    ! Note: The first index represents a data dimension that is the same on each task,
    !        whose size generally is less than own_ewn.
    !       The second index represents the north-south dimension, and is assumed
    !        to have size own_nsn (i.e., the data extend over locally owned cells only).
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which main_task_row will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:),intent(in) :: values
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:),allocatable :: sendbuf

    integer,dimension(4) :: d_gs_mybounds_row
    integer,dimension(:,:),allocatable,save :: d_gs_bounds_row
    logical :: new_bounds_row

    associate(  &
         comm_row      => parallel%comm_row,       &
         tasks_row     => parallel%tasks_row,      &
         this_rank_row => parallel%this_rank_row,  &
         main_rank_row => parallel%main_rank_row,  &
         main_task_row => parallel%main_task_row,  &
         own_nsn       => parallel%own_nsn         &
         )

    if (size(values,2) /= own_nsn) then
       ! Note: Removing this restriction would require some recoding below.
       write(*,*) "ERROR: distributed_gather_var_row requires N-S array size of own_nsn"
       write(*,*) 'rank, own_nsn, size(values,2) =', this_rank, own_nsn, size(values,2)
       call parallel_stop(__FILE__, __LINE__)
    end if

    d_gs_mybounds_row(1) = this_rank_row*size(values,1) + 1
    d_gs_mybounds_row(2) = (this_rank_row+1)*size(values,1)
    d_gs_mybounds_row(3) = 1
    d_gs_mybounds_row(4) = own_nsn

    if (allocated(d_gs_bounds_row)) then
       ! d_gs_bounds_row already computed
       ! Recompute only if there is a mismatch between d_gs_bounds_row and size(values)
       if (d_gs_bounds_row(2,1) - d_gs_bounds_row(1,1) + 1 == size(values,1)) then
          new_bounds_row = .false.   ! use the saved value
       else
          new_bounds_row = .true.    ! recompute
          !WHL - debug
!          if (main_task_row) then
!             print*, this_rank, 'Recompute d_gs_bounds_row'
!             print*, '  current size =', (d_gs_bounds_row(2,:) - d_gs_bounds_row(1,:) + 1)
!             print*, '  size(values,1) =', size(values,1)
!          endif
       endif
    else
       new_bounds_row = .true.
!       if (main_task_row) print*, this_rank, 'Allocate d_gs_bounds_row'
    endif

    !Note: The d_gs_bounds_array is needed only on main_task_col, so memory goes unused on other tasks.
    !      But if we saved memory by allocating an array of size(1,1) on other tasks,
    !       then the new_bounds logic above would not work on all tasks.
    !      The alternative would be to always compute d_gs_bounds_col using fc_gather_int,
    !       but it is more efficient to call mpi_allgather only as needed.
    if (new_bounds_row) then
       if (allocated(d_gs_bounds_row)) deallocate(d_gs_bounds_row)

!       if (main_task_row) then
!          allocate(d_gs_bounds_row(4,tasks_row))
!       else
!          allocate(d_gs_bounds_row(1,1))
!       endif
!       call fc_gather_int(d_gs_mybounds_row,4,mpi_integer,d_gs_bounds_row,4,&
!            mpi_integer,main_rank_row,comm_row)
       allocate(d_gs_bounds_row(4,tasks_row))
       call mpi_allgather(d_gs_mybounds_row,4, mpi_integer,  &
                          d_gs_bounds_row, 4, mpi_integer,  &
                          comm_row, ierror)
    endif

    if (main_task_row) then
       if (allocated(global_values)) deallocate(global_values)
       allocate(global_values(size(values,1)*tasks_row, own_nsn))
       global_values(:,:) = 0.0d0
       allocate(displs(tasks_row+1))
       allocate(recvcounts(tasks_row))
       recvcounts(:) = (d_gs_bounds_row(2,:)-d_gs_bounds_row(1,:)+1)&
                      *(d_gs_bounds_row(4,:)-d_gs_bounds_row(3,:)+1)
       displs(1) = 0
       do i = 1,tasks_row
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks_row+1)))
    else
       if (allocated(global_values)) deallocate(global_values)
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if

    !Note: Would need to uncomment the following if sendbuf were not identical
    !      to the input values array
!!    allocate(sendbuf(d_gs_mybounds_row(1):d_gs_mybounds_row(2),&
!!                     d_gs_mybounds_row(3):d_gs_mybounds_row(4)))
!!    sendbuf(:,:) = values(1:size(values,1), 1:own_nsn)
!!    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
!!       recvbuf,recvcounts,displs,mpi_real8,main_rank_row,comm_row)

    call fc_gatherv_real8(values,size(values),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank_row,comm_row)

    if (main_task_row) then
       do i = 1, tasks_row
          global_values(d_gs_bounds_row(1,i):d_gs_bounds_row(2,i),&
                        d_gs_bounds_row(3,i):d_gs_bounds_row(4,i)) = &
                reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds_row(2,i)-d_gs_bounds_row(1,i)+1,&
                       d_gs_bounds_row(4,i)-d_gs_bounds_row(3,i)+1/))
       end do
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_var_row_real8_2d

!=======================================================================

  ! subroutines belonging to the distributed_gather_all_var_row interface

  subroutine distributed_gather_all_var_row_real8_2d(values, global_values, parallel)

    ! Gather global data along a row of tasks onto each task for that row.
    ! Based on distributed_gather_var_real8_2d.
    ! Note: The first index represents a data dimension that is the same on each task,
    !        whose size generally is less than own_ewn.
    !       The second index represents the north-south dimension, and is assumed
    !        to have size own_nsn (i.e., the data extend over locally owned cells only).
    ! values = local portion of distributed variable
    ! global_values = allocatable array in which each task will store values for that row.
    ! If global_values is allocated, then it will be deallocated and reallocated.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:),intent(in) :: values
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel

    integer :: i,ierror,k
    integer,dimension(:),allocatable :: displs,recvcounts
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:),allocatable :: sendbuf

    integer, dimension(4) :: d_gs_mybounds_row
    integer, dimension(:,:), allocatable, save :: d_gs_bounds_row
    logical :: new_bounds_row

    associate(  &
         comm_row      => parallel%comm_row,       &
         tasks_row     => parallel%tasks_row,      &
         this_rank_row => parallel%this_rank_row,  &
         main_rank_row => parallel%main_rank_row,  &
         main_task_row => parallel%main_task_row,  &
         own_nsn       => parallel%own_nsn         &
         )

    if (size(values,2) /= own_nsn) then
       ! Note: Removing this restriction would require some recoding below.
       ! TODO: Do this recoding.  This subroutine currently fails with outflow BC, because
       !       the southern and western rows of tasks have an extra locally owned vertex,
       !       giving size(values,2) = own_nsn + 1
       write(*,*) "ERROR: distributed_gather_var_row requires N-S array size of own_nsn"
       write(*,*) 'rank, own_nsn, size(values,2) =', this_rank, own_nsn, size(values,2)
       call parallel_stop(__FILE__, __LINE__)
    end if

    d_gs_mybounds_row(1) = this_rank_row*size(values,1) + 1
    d_gs_mybounds_row(2) = (this_rank_row+1)*size(values,1)
    d_gs_mybounds_row(3) = 1
    d_gs_mybounds_row(4) = own_nsn

    if (allocated(d_gs_bounds_row)) then
       ! d_gs_bounds_row already computed
       ! Recompute only if there is a mismatch between d_gs_bounds_row and size(values)
       if (d_gs_bounds_row(2,1) - d_gs_bounds_row(1,1) + 1 == size(values,1)) then
          new_bounds_row = .false.   ! use the saved value
       else
          new_bounds_row = .true.    ! recompute
          !WHL - debug
!          if (main_task_row) then
!             print*, this_rank, 'Recompute d_gs_bounds_row'
!             print*, '  current size =', (d_gs_bounds_row(2,:) - d_gs_bounds_row(1,:) + 1)
!             print*, '  size(values,1) =', size(values,1)
!          endif
       endif
    else
       new_bounds_row = .true.
!       if (main_task_row) print*, this_rank, 'Allocate d_gs_bounds_row'
    endif

    if (new_bounds_row) then
       if (allocated(d_gs_bounds_row)) deallocate(d_gs_bounds_row)
       allocate(d_gs_bounds_row(4,tasks_row))
       call mpi_allgather(d_gs_mybounds_row,4, mpi_integer,  &
                          d_gs_bounds_row, 4, mpi_integer,  &
                          comm_row, ierror)
    endif

    if (allocated(global_values)) deallocate(global_values)
    allocate(global_values(size(values,1)*tasks_row, own_nsn))
    global_values(:,:) = 0.0d0
    allocate(displs(tasks_row+1))
    allocate(recvcounts(tasks_row))
    recvcounts(:) = (d_gs_bounds_row(2,:)-d_gs_bounds_row(1,:)+1)&
                   *(d_gs_bounds_row(4,:)-d_gs_bounds_row(3,:)+1)
    displs(1) = 0
    do i = 1,tasks_row
       displs(i+1) = displs(i)+recvcounts(i)
    end do
    allocate(recvbuf(displs(tasks_row+1)))

    !Note: Would need to uncomment the following and call mpi_allgatherv
    !      with sendbuf arguments if sendbuf were not identical
    !      to the input values array
!!    allocate(sendbuf(d_gs_mybounds_row(1):d_gs_mybounds_row(2),&
!!                     d_gs_mybounds_row(3):d_gs_mybounds_row(4)))
!!    sendbuf(:,:) = values(1:size(values,1), 1:own_nsn)
!!    call mpi_allgatherv(sendbuf, size(sendbuf), mpi_real8, &
!!                        recvbuf, recvcounts, displs, mpi_real8, &
!!                        comm_row, ierror)

    call mpi_allgatherv(values, size(values), mpi_real8, &
                        recvbuf, recvcounts, displs, mpi_real8, &
                        comm_row, ierror)

    do i = 1, tasks_row
       global_values(d_gs_bounds_row(1,i):d_gs_bounds_row(2,i),&
                     d_gs_bounds_row(3,i):d_gs_bounds_row(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                   (/d_gs_bounds_row(2,i)-d_gs_bounds_row(1,i)+1,&
                     d_gs_bounds_row(4,i)-d_gs_bounds_row(3,i)+1/))
    end do

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_all_var_row_real8_2d

!=======================================================================

  ! subroutines belonging to the distributed_gather_var_col interface

  subroutine distributed_gather_var_col_real8_2d(values, global_values, parallel)

    ! Gather data along a column of tasks onto the main task for that column.
    ! Based on distributed_gather_var_real8_2d.
    ! Note: The first index represents a data dimension that is the same on each task,
    !        whose size generally is less than own_nsn.
    !       The second index represents the east-west dimension, and is assumed
    !        to have size own_ewn (i.e., the data extend over locally owned cells only).
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which main_task_col will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:),intent(in) :: values
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,recvcounts
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:),allocatable :: sendbuf

    integer,dimension(4) :: d_gs_mybounds_col
    integer,dimension(:,:),allocatable,save :: d_gs_bounds_col
    logical :: new_bounds_col

    associate(  &
         comm_col      => parallel%comm_col,       &
         tasks_col     => parallel%tasks_col,      &
         this_rank_col => parallel%this_rank_col,  &
         main_rank_col => parallel%main_rank_col,  &
         main_task_col => parallel%main_task_col,  &
         own_ewn       => parallel%own_ewn         &
         )

    if (size(values,2) /= own_ewn) then
       ! Note: Removing this restriction would require some recoding below.
       write(*,*) "ERROR: distributed_gather_var_row requires E-W array size of own_ewn"
       write(*,*) 'rank, own_ewn, size(values,2) =', this_rank, own_ewn, size(values,2)
       call parallel_stop(__FILE__, __LINE__)
    end if

    d_gs_mybounds_col(1) = this_rank_col*size(values,1) + 1
    d_gs_mybounds_col(2) = (this_rank_col+1)*size(values,1)
    d_gs_mybounds_col(3) = 1
    d_gs_mybounds_col(4) = own_ewn

    if (allocated(d_gs_bounds_col)) then
       ! d_gs_bounds_col already computed
       ! Recompute only if there is a mismatch between d_gs_bounds_col and size(values)
       if (d_gs_bounds_col(2,1) - d_gs_bounds_col(1,1) + 1 == size(values,1)) then
          new_bounds_col = .false.   ! use the saved value
       else
          new_bounds_col = .true.    ! recompute
          !WHL - debug
!          if (main_task_col) then
!             print*, this_rank, 'Recompute d_gs_bounds_col'
!             print*, '  current size =', (d_gs_bounds_col(2,:) - d_gs_bounds_col(1,:) + 1)
!             print*, '  size(values,1) =', size(values,1)
!          endif
       endif
    else
       new_bounds_col = .true.
!       if (main_task_col) print*, this_rank, 'Allocate d_gs_bounds_col'
    endif

    !Note: The d_gs_bounds_array is needed only on main_task_col, so memory goes unused on other tasks.
    !      But if we saved memory by allocating an array of size(1,1) on other tasks,
    !       then the new_bounds logic above would not work on all tasks.
    !      The alternative would be to always compute d_gs_bounds_col using fc_gather_int,
    !       but it is more efficient to call mpi_allgather only as needed.

    if (new_bounds_col) then
       if (allocated(d_gs_bounds_col)) deallocate(d_gs_bounds_col)

!       if (main_task_col) then
!          allocate(d_gs_bounds_col(4,tasks_col))
!       else
!          allocate(d_gs_bounds_col(1,1))
!       endif
!       call fc_gather_int(d_gs_mybounds_col,4,mpi_integer,d_gs_bounds_col,4,&
!            mpi_integer,main_rank_col,comm_col)
       allocate(d_gs_bounds_col(4,tasks_col))
       call mpi_allgather(d_gs_mybounds_col,4, mpi_integer,  &
                          d_gs_bounds_col, 4, mpi_integer,  &
                          comm_col, ierror)
    endif

    if (main_task_col) then
       if (allocated(global_values)) deallocate(global_values)
       allocate(global_values(size(values,1)*tasks_col, own_ewn))
       global_values(:,:) = 0.0d0
       allocate(displs(tasks_col+1))
       allocate(recvcounts(tasks_col))
       recvcounts(:) = (d_gs_bounds_col(2,:)-d_gs_bounds_col(1,:)+1)&
                      *(d_gs_bounds_col(4,:)-d_gs_bounds_col(3,:)+1)
       displs(1) = 0
       do i = 1,tasks_col
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks_col+1)))
    else
       if (allocated(global_values)) deallocate(global_values)
       allocate(global_values(1,1))  ! This prevents a problem with NULL pointers later.
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if

    !Note: Would need to uncomment the following if sendbuf were not identical
    !      to the input values array
!!    allocate(sendbuf(d_gs_mybounds_col(1):d_gs_mybounds_col(2),&
!!                     d_gs_mybounds_col(3):d_gs_mybounds_col(4)))
!!    sendbuf(:,:) = values(1:size(values,1), 1:own_ewn)
!!    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
!!       recvbuf,recvcounts,displs,mpi_real8,main_rank_col,comm_col)

    call fc_gatherv_real8(values,size(values),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank_col,comm_col)

    if (main_task_col) then
       do i = 1, tasks_col
          global_values(d_gs_bounds_col(1,i):d_gs_bounds_col(2,i),&
                        d_gs_bounds_col(3,i):d_gs_bounds_col(4,i)) = &
                reshape(recvbuf(displs(i)+1:displs(i+1)), &
                     (/d_gs_bounds_col(2,i)-d_gs_bounds_col(1,i)+1,&
                       d_gs_bounds_col(4,i)-d_gs_bounds_col(3,i)+1/))
       end do
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_var_col_real8_2d

!=======================================================================

  ! subroutines belonging to the distributed_gather_all_var_col interface

  subroutine distributed_gather_all_var_col_real8_2d(values, global_values, parallel)

    ! Gather global data along a column of tasks onto each task for that column.
    ! Based on distributed_gather_var_real8_2d.
    ! Note: The first index represents a data dimension that is the same on each task,
    !        whose size generally is less than own_nsn.
    !       The second index represents the east-west dimension, and is assumed
    !        to have size own_ewn (i.e., the data extend over locally owned cells only).
    ! values = local portion of distributed variable
    ! global_values = allocatable array in which each task will store values for that column.
    ! If global_values is allocated, then it will be deallocated and reallocated.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:),intent(in) :: values
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values
    type(parallel_type) :: parallel

    integer :: i,ierror,k
    integer,dimension(:),allocatable :: displs,recvcounts
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:),allocatable :: sendbuf

    integer,dimension(4) :: d_gs_mybounds_col
    integer,dimension(:,:),allocatable,save :: d_gs_bounds_col
    logical :: new_bounds_col

    associate(  &
         comm_col      => parallel%comm_col,       &
         tasks_col     => parallel%tasks_col,      &
         this_rank_col => parallel%this_rank_col,  &
         main_rank_col => parallel%main_rank_col,  &
         main_task_col => parallel%main_task_col,  &
         own_ewn       => parallel%own_ewn         &
         )

    if (size(values,2) /= own_ewn) then
       ! Note: Removing this restriction would require some recoding below.
       write(*,*) "ERROR: distributed_gather_var_row requires E-W array size of own_ewn"
       write(*,*) 'rank, own_ewn, size(values,2) =', this_rank, own_ewn, size(values,2)
       call parallel_stop(__FILE__, __LINE__)
    end if

    d_gs_mybounds_col(1) = this_rank_col*size(values,1) + 1
    d_gs_mybounds_col(2) = (this_rank_col+1)*size(values,1)
    d_gs_mybounds_col(3) = 1
    d_gs_mybounds_col(4) = own_ewn

    if (allocated(d_gs_bounds_col)) then
       ! d_gs_bounds_col already computed
       ! Recompute only if there is a mismatch between d_gs_bounds_col and size(values)
       if (d_gs_bounds_col(2,1) - d_gs_bounds_col(1,1) + 1 == size(values,1)) then
          new_bounds_col = .false.   ! use the saved value
       else
          new_bounds_col = .true.    ! recompute
          !WHL - debug
!          if (main_task_col) then
!             print*, this_rank, 'Recompute d_gs_bounds_col'
!             print*, '  current size =', (d_gs_bounds_col(2,:) - d_gs_bounds_col(1,:) + 1)
!             print*, '  size(values,1) =', size(values,1)
!          endif
       endif
    else
       new_bounds_col = .true.
!       if (main_task_col) print*, this_rank, 'Allocate d_gs_bounds_col'
    endif

    if (new_bounds_col) then
       if (allocated(d_gs_bounds_col)) deallocate(d_gs_bounds_col)
       allocate(d_gs_bounds_col(4,tasks_col))
       call mpi_allgather(d_gs_mybounds_col,4, mpi_integer,  &
                          d_gs_bounds_col, 4, mpi_integer,  &
                          comm_col, ierror)
    endif

    if (allocated(global_values)) deallocate(global_values)
    allocate(global_values(size(values,1)*tasks_col, own_ewn))
    global_values(:,:) = 0.0d0
    allocate(displs(tasks_col+1))
    allocate(recvcounts(tasks_col))
    recvcounts(:) = (d_gs_bounds_col(2,:)-d_gs_bounds_col(1,:)+1)&
                   *(d_gs_bounds_col(4,:)-d_gs_bounds_col(3,:)+1)
    displs(1) = 0
    do i = 1,tasks_col
       displs(i+1) = displs(i)+recvcounts(i)
    end do
    allocate(recvbuf(displs(tasks_col+1)))

    !Note: Would need to uncomment the following and call mpi_allgatherv
    !      with sendbuf arguments if sendbuf were not identical
    !      to the input values array
!!    allocate(sendbuf(d_gs_mybounds_col(1):d_gs_mybounds_col(2),&
!!                     d_gs_mybounds_col(3):d_gs_mybounds_col(4)))
!!    sendbuf(:,:) = values(1:size(values,1), 1:own_ewn)
!!    call mpi_allgatherv(sendbuf,size(sendbuf),mpi_real8,&
!!                        recvbuf, recvcounts, displs, mpi_real8, &
!!                        comm_col, ierror)

    call mpi_allgatherv(values, size(values), mpi_real8, &
                        recvbuf, recvcounts, displs, mpi_real8, &
                        comm_col, ierror)

    do i = 1, tasks_col
       global_values(d_gs_bounds_col(1,i):d_gs_bounds_col(2,i),&
                     d_gs_bounds_col(3,i):d_gs_bounds_col(4,i)) = &
             reshape(recvbuf(displs(i)+1:displs(i+1)), &
                   (/d_gs_bounds_col(2,i)-d_gs_bounds_col(1,i)+1,&
                     d_gs_bounds_col(4,i)-d_gs_bounds_col(3,i)+1/))
    end do

    end associate
    ! automatic deallocation

  end subroutine distributed_gather_all_var_col_real8_2d

!=======================================================================

  ! functions belonging to the distributed get_var interface
  ! Note: 'parallel' is before 'start' to be consistent with the distributed_put_var interface.
  !       For distributed_put_var, 'parallel' must go first since 'start' can be optional.

  function distributed_get_var_integer_2d(ncid, varid, values, parallel, start)

    use mpi_mod
    implicit none
    integer :: distributed_get_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values
    type(parallel_type) :: parallel

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: sendbuf
    integer,dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    associate(  &
         ewlb       => parallel%ewlb,       &
         ewub       => parallel%ewub,       &
         nslb       => parallel%nslb,       &
         nsub       => parallel%nsub,       &
         local_ewn  => parallel%local_ewn,  &
         local_nsn  => parallel%local_nsn,  &
         global_ewn => parallel%global_ewn, &
         global_nsn => parallel%global_nsn, &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(global_minval_ewlb:global_maxval_ewub, &
                              global_minval_nslb:global_maxval_nsub))
       global_values(:,:) = 0
       distributed_get_var_integer_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_integer_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_integer,&
         recvbuf,size(recvbuf),mpi_integer,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))

    end associate
    !automatic deallocation

  end function distributed_get_var_integer_2d


  function distributed_get_var_real4_1d(ncid, varid, values, parallel, start)

    use mpi_mod
    use netcdf
    implicit none
    integer :: distributed_get_var_real4_1d,ncid,varid
    integer,dimension(:) :: start
    real(sp),dimension(:) :: values
    type(parallel_type) :: parallel

    integer :: i,ierror,myn,status,x1id,y1id
    integer,dimension(2) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(sp),dimension(:),allocatable :: global_values,sendbuf

    ! begin

    associate(  &
         ewlb       => parallel%ewlb,       &
         ewub       => parallel%ewub,       &
         nslb       => parallel%nslb,       &
         nsub       => parallel%nsub,       &
         local_ewn  => parallel%local_ewn,  &
         local_nsn  => parallel%local_nsn,  &
         global_ewn => parallel%global_ewn, &
         global_nsn => parallel%global_nsn, &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x1id)
    call broadcast(y1id)
    if (varid==x1id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub
       myn = global_ewn
    else if (varid==y1id) then
       mybounds(1) = nslb
       mybounds(2) = nsub
       myn = global_nsn
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:))))
       if (varid==x1id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub))
       elseif (varid==y1id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub))
       endif
       global_values(:) = 0.0
       distributed_get_var_real4_1d = &
            nf90_get_var(ncid,varid,global_values(1:myn),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = bounds(2,:)-bounds(1,:)+1
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
               global_values(bounds(1,i):bounds(2,i))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real4_1d)
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         values,size(values),mpi_real4,main_rank,comm,ierror)

    end associate
    !automatic deallocation

  end function distributed_get_var_real4_1d


  function distributed_get_var_real4_2d(ncid, varid, values, parallel, start)

    use mpi_mod
    implicit none
    integer :: distributed_get_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(sp),dimension(:,:) :: values
    type(parallel_type) :: parallel

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(sp),dimension(:),allocatable :: sendbuf
    real(sp),dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    associate(  &
         ewlb       => parallel%ewlb,       &
         ewub       => parallel%ewub,       &
         nslb       => parallel%nslb,       &
         nsub       => parallel%nsub,       &
         local_ewn  => parallel%local_ewn,  &
         local_nsn  => parallel%local_nsn,  &
         global_ewn => parallel%global_ewn, &
         global_nsn => parallel%global_nsn, &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(global_minval_ewlb:global_maxval_ewub, &
                              global_minval_nslb:global_maxval_nsub))
       global_values(:,:) = 0.0
       distributed_get_var_real4_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real4_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         recvbuf,size(recvbuf),mpi_real4,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))

    end associate
    !automatic deallocation

  end function distributed_get_var_real4_2d


  function distributed_get_var_real8_1d(ncid, varid, values, parallel, start)

    use mpi_mod
    use netcdf
    implicit none
    integer :: distributed_get_var_real8_1d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:) :: values
    type(parallel_type) :: parallel

    integer :: i,ierror,myn,status,x1id,y1id
    integer,dimension(2) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(dp),dimension(:),allocatable :: global_values,sendbuf

    ! begin

    associate(  &
         ewlb       => parallel%ewlb,       &
         ewub       => parallel%ewub,       &
         nslb       => parallel%nslb,       &
         nsub       => parallel%nsub,       &
         local_ewn  => parallel%local_ewn,  &
         local_nsn  => parallel%local_nsn,  &
         global_ewn => parallel%global_ewn, &
         global_nsn => parallel%global_nsn, &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x1id)
    call broadcast(y1id)
    if (varid==x1id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub
       myn = global_ewn
    else if (varid==y1id) then
       mybounds(1) = nslb
       mybounds(2) = nsub
       myn = global_nsn
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:))))
       if (varid==x1id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub))
       elseif (varid==y1id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub))
       endif
       global_values(:) = 0.0d0
       distributed_get_var_real8_1d = &
            nf90_get_var(ncid,varid,global_values(1:myn),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = bounds(2,:)-bounds(1,:)+1
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
               global_values(bounds(1,i):bounds(2,i))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real8_1d)
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         values,size(values),mpi_real8,main_rank,comm,ierror)

    end associate
    !automatic deallocation

  end function distributed_get_var_real8_1d


  function distributed_get_var_real8_2d(ncid, varid, values, parallel, start)

    use mpi_mod
    implicit none
    integer :: distributed_get_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:,:) :: values
    type(parallel_type) :: parallel

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(dp),dimension(:),allocatable :: sendbuf
    real(dp),dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    associate(  &
         ewlb       => parallel%ewlb,       &
         ewub       => parallel%ewub,       &
         nslb       => parallel%nslb,       &
         nsub       => parallel%nsub,       &
         local_ewn  => parallel%local_ewn,  &
         local_nsn  => parallel%local_nsn,  &
         global_ewn => parallel%global_ewn, &
         global_nsn => parallel%global_nsn, &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)

    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(global_minval_ewlb:global_maxval_ewub, &
                              global_minval_nslb:global_maxval_nsub))
       global_values(:,:) = 0.0d0
       distributed_get_var_real8_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real8_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))

    end associate
    !automatic deallocation

  end function distributed_get_var_real8_2d


  function distributed_get_var_real8_3d(ncid, varid, values, parallel, start)

    use mpi_mod
    implicit none
    integer :: distributed_get_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:,:,:) :: values
    type(parallel_type) :: parallel

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(dp),dimension(:),allocatable :: sendbuf
    real(dp),dimension(:,:,:),allocatable :: global_values,recvbuf

    ! begin

    associate(  &
         ewlb       => parallel%ewlb,       &
         ewub       => parallel%ewub,       &
         nslb       => parallel%nslb,       &
         nsub       => parallel%nsub,       &
         local_ewn  => parallel%local_ewn,  &
         local_nsn  => parallel%local_nsn,  &
         global_ewn => parallel%global_ewn, &
         global_nsn => parallel%global_nsn, &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:)),size(values,3)))
       allocate(global_values(global_minval_ewlb:global_maxval_ewub, &
                              global_minval_nslb:global_maxval_nsub, &
                              size(values,3)))
       global_values(:,:,:) = 0.0d0
       distributed_get_var_real8_3d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns,:),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*&
            (bounds(4,:)-bounds(3,:)+1)*size(values,3)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(global_values(&
               bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i),:),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real8_3d)
    allocate(recvbuf(local_ewn,local_nsn,size(values,3)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(:,:,:) = recvbuf(:size(values,1),:size(values,2),:)

    end associate
    !automatic deallocation

  end function distributed_get_var_real8_3d

!=======================================================================

  subroutine distributed_grid(ewn,      nsn,        &
                              parallel,             &
                              nhalo_in, global_bc_in)

    ! Divide the global domain into blocks, with one task per block.
    ! Set various grid and domain variables for the local task.

    implicit none
    integer, intent(inout) :: ewn, nsn                  ! global grid dimensions
    type(parallel_type), intent(inout) :: parallel      ! info for parallel communication, computed here
    integer, intent(in), optional :: nhalo_in           ! number of rows of halo cells
    character(*), intent(in), optional :: global_bc_in  ! string indicating the global BC option

    integer :: best,i,j,metric
    real(dp) :: rewtasks,rnstasks
    integer :: ProcsEW

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         global_ewn  => parallel%global_ewn,   &
         global_nsn  => parallel%global_nsn,   &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn,      &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         ewtasks     => parallel%ewtasks,      &
         nstasks     => parallel%nstasks,      &
         ewrank      => parallel%ewrank,       &
         nsrank      => parallel%nsrank,       &
         global_col_offset  => parallel%global_col_offset,  &
         global_row_offset  => parallel%global_row_offset,  &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub, &
         southwest_corner   => parallel%southwest_corner,   &
         southeast_corner   => parallel%southeast_corner,   &
         northeast_corner   => parallel%northeast_corner,   &
         northwest_corner   => parallel%northwest_corner,   &
         staggered_ilo      => parallel%staggered_ilo,      &
         staggered_ihi      => parallel%staggered_ihi,      &
         staggered_jlo      => parallel%staggered_jlo,      &
         staggered_jhi      => parallel%staggered_jhi       &
         )

    ! set the boundary conditions (periodic by default)
    ! Note: The no-penetration BC is treated as periodic. This BC may need some more work.

    if (present(global_bc_in)) then
       if (trim(global_bc_in) == 'periodic') then
          periodic_bc = .true.
          outflow_bc = .false.
          no_ice_bc = .false.
          if (main_task) write(*,*) 'Setting periodic boundary conditions'
       elseif (trim(global_bc_in) == 'outflow') then
          periodic_bc = .false.
          outflow_bc = .true.
          no_ice_bc = .false.
          if (main_task) write(*,*) 'Setting outflow boundary conditions'
       elseif (trim(global_bc_in) == 'no_penetration') then
          periodic_bc = .true.   ! Currently use the same halo logic for no_penetration and periodic  
          outflow_bc = .false.
          no_ice_bc = .false.
          if (main_task) write(*,*) 'Setting no_penetration boundary conditions'
       elseif (trim(global_bc_in) == 'no_ice') then
          periodic_bc = .false.
          outflow_bc = .false.
          no_ice_bc = .true.
          if (main_task) write(*,*) 'Setting no_ice boundary conditions'
       else
          if (main_task) write(*,*) 'Error: Invalid global_bc option for distributed_grid subroutine'
          call parallel_stop(__FILE__, __LINE__)
       endif
    else   ! default to periodic
       periodic_bc = .true.
       outflow_bc = .false.
       no_ice_bc = .false.
    endif

    ! Optionally, change the halo values
    ! Note: The Glissade higher-order dycore requires nhalo = 2.
    !       The Glide SIA dycore requires nhalo = 0.
    ! The default halo values at the top of the module are appropriate for
    !  the higher-order dycores.  Here they can be reset to zero for Glide.

    if (present(nhalo_in)) then
       if (main_task) then
          write(*,*) 'Setting halo values: nhalo =', nhalo_in
          if (nhalo_in < 0) then
             write(*,*) 'ERROR: nhalo must be >= 0'
             call parallel_stop(__FILE__, __LINE__)
          elseif (nhalo_in /= 2) then
             write(*,*) 'WARNING: parallel dycores tested only with nhalo = 2'
          endif
       endif 
       nhalo = nhalo_in
       lhalo = nhalo
       uhalo = nhalo
       staggered_lhalo = lhalo
       staggered_uhalo = max(uhalo-1, 0)
    endif

    global_ewn = ewn
    global_nsn = nsn

    ewtasks = 0
    nstasks = 0
    best = huge(best)
    do i = 1,min(tasks,global_ewn)
       j = tasks/i
       if (j<=global_nsn.and.i*j==tasks) then ! try to use all tasks
          metric = abs(i*global_nsn-j*global_ewn) ! zero if ewn/nsn == i/j
          if (metric<best) then
             best = metric
             ewtasks = i
             nstasks = j
          end if
       end if
    end do
    if (ewtasks*nstasks/=tasks) call parallel_stop(__FILE__,__LINE__)

    ! Store critical value for creating global IDs.  Defines grid distribution.
    ProcsEW = ewtasks

    ! For globalID calculations, determine processor's global grid index offsets.
    ! Sum block sizes for row blocks preceding this_rank.
    ! Do not include halo offsets in global calculations.
    ! (There are ProcsEW processors per row.)
    global_col_offset = 0
    do ewrank=0,mod(this_rank, ProcsEW)-1
      rewtasks = 1/real(ewtasks,dp)
      ewlb = nint(ewrank*global_ewn*rewtasks)+1
      ewub = nint((ewrank+1)*global_ewn*rewtasks)
      own_ewn = ewub-ewlb+1
      global_col_offset = global_col_offset + own_ewn
    enddo

    ! Sum block sizes for column blocks preceding this_rank
    ! (Integer division required for this_rank/ProcsEW)
    global_row_offset = 0
    do nsrank=0,(this_rank/ProcsEW)-1
      rnstasks = 1/real(nstasks,dp)
      nslb = nint(nsrank*global_nsn*rnstasks)+1
      nsub = nint((nsrank+1)*global_nsn*rnstasks)
      own_nsn = nsub-nslb+1
      global_row_offset = global_row_offset + own_nsn
    enddo

    ! Set local processor's grid indices, including halo offsets
    ewrank = mod(this_rank,ewtasks)
    rewtasks = 1/real(ewtasks,dp)
    ewlb = nint(ewrank*global_ewn*rewtasks)+1-lhalo
    ewub = nint((ewrank+1)*global_ewn*rewtasks)+uhalo
    local_ewn = ewub-ewlb+1
    own_ewn = local_ewn-lhalo-uhalo
    ewn = local_ewn

    nsrank = this_rank/ewtasks
    rnstasks = 1/real(nstasks,dp)
    nslb = nint(nsrank*global_nsn*rnstasks)+1-lhalo
    nsub = nint((nsrank+1)*global_nsn*rnstasks)+uhalo
    local_nsn = nsub-nslb+1
    own_nsn = local_nsn-lhalo-uhalo
    nsn = local_nsn

    ! Determine the global bounds of ewlb, ewub, nslb, and nsub.
    ! These are used to allocate global arrays used in gathers and scatters.
    global_minval_ewlb = parallel_reduce_min(ewlb)
    global_maxval_ewub = parallel_reduce_max(ewub)
    global_minval_nslb = parallel_reduce_min(nslb)
    global_maxval_nsub = parallel_reduce_max(nsub)

    west = this_rank-1
    if ((west/ewtasks<this_rank/ewtasks).or.(west<0)) west = west+ewtasks
    east = this_rank+1
    if (east/ewtasks>this_rank/ewtasks) east = east-ewtasks
    south = this_rank-ewtasks
    if (south<0) south = south+tasks
    north = this_rank+ewtasks
    if (north>=tasks) north = north-tasks

    ! Specify that this is not a corner task.
    ! A southeast corner task is defined as a task with active south and east neighbors but an
    !  inactive southeast neighbor; and similarly for other directions.
    ! This subroutine assumes a standard rectangular layout of tasks, in which case
    !  there are no such corner tasks.
    ! Subroutine distributed_grid_active_blocks (below) typically distributes tasks
    !  such that there are corner tasks.
    southwest_corner = .false.
    southeast_corner = .false.
    northeast_corner = .false.
    northwest_corner = .false.

    ! Set the limits of locally owned vertices on the staggered grid.
    ! For periodic BC, staggered_ilo = staggered_jlo = staggered_lhalo+1.
    ! For outflow BC the locally owned vertices include the southern and western rows
    !  of the global domain, so staggered_ilo = staggered_jlo = staggered_lhalo on
    !  processors that include these rows.
    ! Note: For no_ice BC, we assume (uvel,vvel) = 0 along the global boundary.
    !       In this case, vertices along the southern and western edges of the global boundary
    !        are not considered to be locally owned by any task.

    if (outflow_bc .and. this_rank <= west) then  ! on west edge of global domain
       staggered_ilo = staggered_lhalo
    else
       staggered_ilo = staggered_lhalo+1
    endif
    staggered_ihi = ewn - 1 - staggered_uhalo

    if (outflow_bc .and. this_rank <= south) then  ! on south edge of global domain
       staggered_jlo = staggered_lhalo
    else
       staggered_jlo = staggered_lhalo+1
    endif
    staggered_jhi = nsn - 1 - staggered_uhalo
    
    ! Check that we have not split up the problem too much.  We do not want halos overlapping in either dimension.
    ! local_* - lhalo - uhalo is the actual number of non-halo cells on a processor.
    if ((local_nsn - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(*,*) "NS halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    if ((local_ewn  - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(*,*) "EW halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

!    call parallel_barrier
!    print*, 'task, west, east, south, north:', this_rank, west, east, south, north

    ! Uncomment to print grid geometry
!    write(*,*) "Process ", this_rank, " Total = ", tasks, " ewtasks = ", ewtasks, " nstasks = ", nstasks
!    write(*,*) "Process ", this_rank, " ewrank = ", ewrank, " nsrank = ", nsrank
!    write(*,*) "Process ", this_rank, " l_ewn = ", local_ewn, " o_ewn = ", own_ewn
!    write(*,*) "Process ", this_rank, " l_nsn = ", local_nsn, " o_nsn = ", own_nsn
!    write(*,*) "Process ", this_rank, " ewlb = ", ewlb, " ewub = ", ewub
!    write(*,*) "Process ", this_rank, " nslb = ", nslb, " nsub = ", nsub
!    write(*,*) "Process ", this_rank, " east = ", east, " west = ", west
!    write(*,*) "Process ", this_rank, " north = ", north, " south = ", south
!    write(*,*) "Process ", this_rank, " ew_vars = ", own_ewn, " ns_vars = ", own_nsn
!    write(*,*) "Process ", this_rank, " global_col_offset = ", global_col_offset, &
!                                      " global_row_offset = ", global_row_offset

    call distributed_print_grid(own_ewn, own_nsn)

    end associate

  end subroutine distributed_grid

!=======================================================================

  subroutine distributed_grid_active_blocks(ewn,      nsn,        &
                                            nx_block, ny_block,   &
                                            ice_domain_mask,      &
                                            parallel,             &
                                            inquire_only)

    ! Divide the global domain into blocks, setting various grid and domain variables
    !  for each block as in subroutine distributed_grid above.
    ! Then read a mask to identify which blocks are potentially active
    !  (i.e, ice can be present in some or all of the block).
    !  Assign a task to each active block.  Activate additional blocks as needed
    !  so that there is exactly one task per active block.
    ! Assume no_ice boundary conditions.  For these BCs, scalars (including ice thickness)
    !  are set to zero not only in the global halo (as for outflow BCs), but also
    !  along one row just inside the global boundary.  This ensures that velocity = 0
    !  for vertices along the global boundary, and allows halo routines to work
    !  correctly when some tasks are inactive.
    ! Set the neighbor task indices (west, east, south, north) such that if a task
    !  does not have an active west neighbor, it has west = this_rank (and similarly
    !  for other directions).
    ! Set logical variables for corner tasks.  A southeast corner task is defined as
    !  a task with active south and east neighbors but an inactive southeast neighbor
    !  (and similarly for other directions).
    !
    ! Note: The code does not yet compute the total ice mass zeroed out near the global boundary.
    !       This does not result in a mass conservation error, because the total dmass/dt term
    !        is computed by summing over dH/dt in locally owned cells only.
    !       There is a similar issue with outflow BC; the code does not compute the ice mass
    !        zeroed out in the global halo.
    ! TODO: Compute dmass/dt by summing over the global ice mass each time step?
    !       Then compute an outflow flux to bring mass conservation back into balance.
    !       This would require computing the total ice mass before and after halo updates.x

    ! Created by WHL, July 2019, based on subroutine distributed_grid above.

    implicit none

    integer, intent(inout) :: ewn, nsn              ! global grid dimensions
    integer, intent(in) :: nx_block, ny_block       ! block sizes in each direction
    integer, intent(in), dimension(:,:) :: &
         ice_domain_mask                            ! = 1 where ice is potentially present and active, else = 0
    type(parallel_type), intent(inout) :: parallel  ! info for parallel communication, computed here
    logical, intent(in), optional :: inquire_only   ! if true, then report the number of active blocks and abort

    integer :: i, j, nb, nt
    integer :: nblocks               ! number of blocks = ewtasks * nstasks
    real(dp) :: rewtasks, rnstasks
    integer :: ProcsEW

    integer :: nblocks_active        ! number of active blocks
    logical :: only_inquire          ! local version of inquire_only
    integer :: corner_block          ! block number for a corner neighbor

    ! arrays with dimension 'nblocks'

    integer, dimension(:), allocatable ::  &
         ewrank_block,             &  ! grid column in which a block lies, increasing from W to E
         nsrank_block,             &  ! grid row in which a block lies, increasing from S to N
         ewlb_block, ewub_block,   &  ! lower and upper bounds in E-W direction for each block, including halos
         nslb_block, nsub_block,   &  ! lower and upper bounds in N-S direction for each block, including halos
         own_ewn_block,            &  ! number of cells in E-W direction on each block   !TODO - rename?  Not "owned"
         own_nsn_block,            &  ! number of cells in N-S direction on each block   !TODO - rename?  Not "owned"
         global_col_offset_block,  &  ! column offset for each block
         global_row_offset_block,  &  ! row offset for each block
         local_ewn_block,          &  ! number of cells in E-W direction on each block, including halos
         local_nsn_block,          &  ! number of cells in N-S direction on each block, including halos
         ewn_block,                &  ! set to local_ewn_block
         nsn_block,                &  ! set to local_nsn_block
         east_block, west_block,   &  ! east and west neighbors of each block
         north_block, south_block     ! north and south neighbors of each block

    integer, dimension(:), allocatable ::  &
         block_to_task                ! task (if any) assigned to each block

    logical, dimension(:), allocatable :: &
         block_is_active,          &  ! true for active blocks
         block_is_active_new          ! augmented version of block_is_active

    ! arrays with dimension 'tasks'

    integer, dimension(:), allocatable ::  &
         task_to_block                ! block associated with each task

    logical :: verbose_active_blocks = .false.

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         global_ewn  => parallel%global_ewn,   &
         global_nsn  => parallel%global_nsn,   &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn,      &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         ewtasks     => parallel%ewtasks,      &
         nstasks     => parallel%nstasks,      &
         ewrank      => parallel%ewrank,       &
         nsrank      => parallel%nsrank,       &
         global_col_offset  => parallel%global_col_offset,  &
         global_row_offset  => parallel%global_row_offset,  &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub, &
         southwest_corner   => parallel%southwest_corner,   &
         southeast_corner   => parallel%southeast_corner,   &
         northeast_corner   => parallel%northeast_corner,   &
         northwest_corner   => parallel%northwest_corner,   &
         staggered_ilo      => parallel%staggered_ilo,      &
         staggered_ihi      => parallel%staggered_ihi,      &
         staggered_jlo      => parallel%staggered_jlo,      &
         staggered_jhi      => parallel%staggered_jhi       &
         )

    if (present(inquire_only)) then
       only_inquire = inquire_only
       if (only_inquire) verbose_active_blocks = .true.
    else
       only_inquire = .false.
    endif

    ! Set the boundary conditions.
    ! For the active_blocks option, only the no_ice BC is supported.

    no_ice_bc = .true.
    outflow_bc = .false.
    periodic_bc = .false.

    if (main_task) write(*,*) 'Setting no_ice boundary conditions, assigning tasks to active blocks only'

    ! Note: There is no option to change nhalo.
    !       The Glissade higher-order dycore requires nhalo = 2, which is set with other halo values at the top of this module.
    !       The Glide SIA dycore (with nhalo = 0) is serial-only and is not compatible with the active-block option.

    global_ewn = ewn
    global_nsn = nsn

    ! Given the global dimensions global_ewn and global_nsn, along with the specified block size nx_block, ny_block,
    ! determine the number of blocks needed to span the global domain in each direction.
    ! Note: I considered renaming ewtasks and nstasks to ewblocks and nsblocks, but kept the old names
    !       because they're used elsewhere in the code.  With the active_block logic, the domain of active tasks
    !       might not fully span the region defined by (ewtasks,nstasks).

    if (nx_block == 0 .or. ny_block == 0) then
       write(*,*) 'Error: Must have nx_block and ny_block > 0'
       call parallel_stop(__FILE__, __LINE__)
    endif

    if (mod(global_ewn, nx_block) > 0) then
       ewtasks = global_ewn/nx_block + 1
    else
       ewtasks = global_ewn/nx_block
    endif

    if (mod(global_nsn, ny_block) > 0) then
       nstasks = global_nsn/ny_block + 1
    else
       nstasks = global_nsn/ny_block
    endif

    nblocks = ewtasks * nstasks

    if (main_task .and. verbose_active_blocks) then
       print*, 'global_ewn, global_nsn =', global_ewn, global_nsn
       print*, 'nx_block, ny_block =', nx_block, ny_block
       print*, 'ewtasks, nstasks, nblocks =', ewtasks, nstasks, nblocks
    endif

    ! Allocate block variables
    allocate(ewrank_block(0:nblocks-1))
    allocate(global_col_offset_block(0:nblocks-1))
    allocate(ewlb_block(0:nblocks-1))
    allocate(ewub_block(0:nblocks-1))
    allocate(local_ewn_block(0:nblocks-1))
    allocate(own_ewn_block(0:nblocks-1))
    allocate(ewn_block(0:nblocks-1))

    allocate(nsrank_block(0:nblocks-1))
    allocate(global_row_offset_block(0:nblocks-1))
    allocate(nslb_block(0:nblocks-1))
    allocate(nsub_block(0:nblocks-1))
    allocate(local_nsn_block(0:nblocks-1))
    allocate(own_nsn_block(0:nblocks-1))
    allocate(nsn_block(0:nblocks-1))

    allocate(block_is_active(0:nblocks-1))
    allocate(block_is_active_new(0:nblocks-1))

    allocate(west_block(0:nblocks-1))
    allocate(east_block(0:nblocks-1))
    allocate(south_block(0:nblocks-1))
    allocate(north_block(0:nblocks-1))

    ! Determine properties of each block.
    ! Since not all blocks are active, we may not have a one-to-one correspondence between
    !  blocks and tasks, so we cannot yet compute task-specific properties on the local task.
    ! Note: These calculations could potentially be done on main_task only, then broadcast to local tasks.
    !       Doing them on all tasks avoids the need to broadcast arrays of size (nblocks) to all local tasks.

    ! For globalID calculations, determine each block's global grid index offsets.
    ! Do not include halo offsets in global calculations.

    ! Store critical value for creating global IDs.  Defines grid distribution.
    !Note: Currently, ProcsEW is used only in profile.F90; could simply replace with ewtasks?
    ProcsEW = ewtasks

    ! Loop over blocks.  Note zero-based indexing.
    do nb = 0, nblocks-1

       ! Sum block sizes for row blocks preceding each block.
       ! The number of blocks per row is ewtasks.

       global_col_offset_block(nb) = 0
       do ewrank = 0, mod(nb,ewtasks) - 1
          rewtasks = 1.0d0/real(ewtasks,dp)
          ewlb_block(nb) = nint( ewrank    * global_ewn * rewtasks) + 1
          ewub_block(nb) = nint((ewrank+1) * global_ewn * rewtasks)
          own_ewn_block(nb) = ewub_block(nb) - ewlb_block(nb) + 1
          global_col_offset_block(nb) = global_col_offset_block(nb) + own_ewn_block(nb)
       enddo

       ! Sum block sizes for column blocks preceding each block.
       ! The number of blocks per column is nstasks.
       ! (Integer division required for nb/ewtasks)

       global_row_offset_block(nb) = 0
       do nsrank = 0, (nb/ewtasks) - 1
          rnstasks = 1.0d0/real(nstasks,dp)
          nslb_block(nb) = nint( nsrank * global_nsn * rnstasks) + 1
          nsub_block(nb) = nint((nsrank+1) * global_nsn * rnstasks)
          own_nsn_block(nb) = nsub_block(nb) - nslb_block(nb) + 1
          global_row_offset_block(nb) = global_row_offset_block(nb) + own_nsn_block(nb)
       enddo

       ! Set each block's grid indices, including halo offsets

       ewrank_block(nb) = mod(nb, ewtasks)
       rewtasks = 1.0d0/real(ewtasks,dp)
       ewlb_block(nb) = nint( ewrank * global_ewn * rewtasks) + 1 - lhalo
       ewub_block(nb) = nint((ewrank+1) * global_ewn * rewtasks) + uhalo
       local_ewn_block(nb) = ewub_block(nb) - ewlb_block(nb) + 1
       own_ewn_block(nb) = local_ewn_block(nb) - lhalo - uhalo
       ewn_block(nb) = local_ewn_block(nb)

       nsrank_block(nb) = nb / ewtasks
       rnstasks = 1.0d0/real(nstasks,dp)
       nslb_block(nb) = nint( nsrank * global_nsn * rnstasks) + 1 - lhalo
       nsub_block(nb) = nint((nsrank+1) * global_nsn * rnstasks) + uhalo
       local_nsn_block(nb) = nsub_block(nb) - nslb_block(nb) + 1
       own_nsn_block(nb) = local_nsn_block(nb) - lhalo - uhalo
       nsn_block(nb) = local_nsn_block(nb)

       ! Identify the west, east, south, and north neighbors of each block
       !
       ! Periodic BCs satisfy the convention that for blocks on the western edge of the domain,
       !  'west' is the easternmost block of that row, and for blocks on the easter edge of the domain,
       ! 'east' is the westernmost block of that row.  And similiarly for south and north.
       !
       ! For no_ice BC, we have a different convention: For blocks along the western edge of the domain
       !  (where west_block/ewtasks < nb/ewtasks by integer division), west_block is set to nb.
       ! According to the halo logic, when this_rank <= west (and specifically when this_rank = west),
       !  no message is sent to the west, and none will be received from the west.
       ! And similarly for other directions.

       west_block(nb) = nb - 1
       if ((west_block(nb)/ewtasks < nb/ewtasks) .or. (west_block(nb) < 0)) west_block(nb) = nb

       east_block(nb) = nb + 1
       if (east_block(nb)/ewtasks > nb/ewtasks) east_block(nb) = nb

       south_block(nb) = nb - ewtasks
       if (south_block(nb) < 0) south_block(nb) = nb

       north_block(nb) = nb + ewtasks
       if (north_block(nb) >= nblocks) north_block(nb) = nb

    enddo   ! nblocks

    ! Determine the global bounds of ewlb, ewub, nslb, and nsub.
    ! These are used to allocate global arrays used in gathers and scatters.
    global_minval_ewlb = minval(ewlb_block)
    global_maxval_ewub = maxval(ewub_block)
    global_minval_nslb = minval(nslb_block)
    global_maxval_nsub = maxval(nsub_block)

    !WHL - debug
    if (main_task) then
       print*, 'global_minval_ewlb =', global_minval_ewlb
       print*, 'global_maxval_ewub =', global_maxval_ewub
       print*, 'global_minval_nslb =', global_minval_nslb
       print*, 'global_maxval_nsub =', global_maxval_nsub
    endif

    ! Determine which blocks are active.
    ! A block is active if one or more locally owned cells has ice_domain_mask = 1.
    ! Note: This calculation is done on main_task only, since this is where
    !        the global array ice_domain_mask is nonzero.
    !       The following loop could be expensive on fine global grids.

    block_is_active(:) = .false.
    nblocks_active = 0

    if (main_task) then
       do nb = 0, nblocks-1
          ij_outer: do j = nslb_block(nb) + lhalo, nsub_block(nb) - uhalo
             do i = ewlb_block(nb) + lhalo, ewub_block(nb) - uhalo
                if (ice_domain_mask(i,j) == 1) then
                   block_is_active(nb) = .true.
                   nblocks_active = nblocks_active + 1
                   exit ij_outer
                endif
             enddo
          enddo ij_outer
       enddo

       if (verbose_active_blocks) then
          print*, 'nblocks, nblocks_active:', nblocks, nblocks_active
          print*, ' '
          print*, 'Block layout:'
          do j = nstasks-1, 0, -1
             do i = 0, ewtasks-1
                nb = ewtasks*j + i
                write(6, '(i5)', advance='no') nb
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'block_is_active:'
          do j = nstasks-1, 0, -1
             do i = 0, ewtasks-1
                nb = ewtasks*j + i
                write(6, '(l5)', advance='no') block_is_active(nb)
             enddo
             print*, ' '
          enddo
       endif  ! verbose_active_blocks

       if (only_inquire) then  ! report the number of active blocks, then abort cleanly
          write(*,*)
          write(*,*) 'The number of active blocks with this domain and block layout is ', nblocks_active
          write(*,*) 'The total number of blocks is ', nblocks
          write(*,*) 'Please resubmit with nblocks_active <= tasks <= n_blocks'
          call parallel_stop(__FILE__, __LINE__)
       else   ! abort if tasks < nblocks_active or tasks > nblocks; otherwise proceed
          if (tasks < nblocks_active) then
             write(*,*)
             write(*,*) 'Fatal error: tasks < nblocks_active'
             write(*,*) 'Number of tasks =', tasks
             write(*,*) 'Minimum number of tasks to compute on all active blocks is ', nblocks_active
             call parallel_stop(__FILE__, __LINE__)
          elseif (tasks > nblocks) then
             write(*,*)
             write(*,*) 'Fatal error: tasks > nblocks'
             write(*,*) 'Number of tasks =', tasks
             write(*,*) 'Maximum number of tasks to compute on all blocks is ', nblocks
             call parallel_stop(__FILE__, __LINE__)
          endif
       endif   ! only_inquire

       ! If nblocks_active < tasks, then activate additional blocks along the ice sheet boundary until
       !  nblocks_active = tasks. Do not activate "orphan" blocks disconnected from other active blocks.

       if (nblocks_active < tasks) then

          if (verbose_active_blocks) then
             print*, 'nblocks_active, tasks =', nblocks_active, tasks
             print*, 'Adding more active blocks:'
          endif

          do while (nblocks_active < tasks)

             block_is_active_new(:) = block_is_active(:)

             do nb = 0, nblocks-1
                if (.not. block_is_active(nb)) then
                   if (block_is_active(west_block(nb))  .or. block_is_active(east_block(nb)) .or.  &
                       block_is_active(south_block(nb)) .or. block_is_active(north_block(nb))) then
                      if (verbose_active_blocks) print*, 'Activate block', nb
                      block_is_active_new(nb) = .true.
                      nblocks_active = nblocks_active + 1
                      if (nblocks_active == tasks) exit
                   endif
                endif
             enddo   ! nb

             block_is_active(:) = block_is_active_new(:)

          enddo   ! nblocks_active < tasks

          if (nblocks_active < tasks) then  ! should not happen, but check just in case
             write(*,*) 'Error, still have nblocks_active < tasks:', nblocks_active, tasks
             call parallel_stop(__FILE__, __LINE__)
          endif

          if (verbose_active_blocks) then
             print*, 'After activating more blocks: nblocks, nblocks_active:', nblocks, nblocks_active
             print*, ' '
             print*, 'Block layout:'
             do j = nstasks-1, 0, -1
                do i = 0, ewtasks-1
                   nb = ewtasks*j + i
                   write(6, '(i5)', advance='no') nb
                enddo
                print*, ' '
             enddo
             print*, ' '
             print*, 'block_is_active:'
             do j = nstasks-1, 0, -1
                do i = 0, ewtasks-1
                   nb = ewtasks*j + i
                   write(6, '(l5)', advance='no') block_is_active(nb)
                enddo
                print*, ' '
             enddo
          endif  ! verbose_active_blocks

       endif   ! nblocks_active < tasks

    endif   ! main_task

    ! Now that we have one task per active block, set up grid info on the local task.

    ! Broadcast active blocks to all tasks
    call broadcast(block_is_active)
    nblocks_active = count(block_is_active)

    ! Set up a correspondence between tasks and active blocks
    allocate(block_to_task(0:nblocks-1))
    allocate(task_to_block(0:tasks-1))

    block_to_task(:) = -1   ! do not set to 0, since 0 is a valid block and task number
    task_to_block(:) = -1

    nt = -1
    do nb = 0, nblocks-1
       if (block_is_active(nb)) then
          nt = nt + 1
          task_to_block(nt) = nb
          block_to_task(nb) = nt
       endif
    enddo

    if (main_task .and. verbose_active_blocks) then
       print*, ' '
       print*, 'Task layout:'
       do j = nstasks-1, 0, -1
          do i = 0, ewtasks-1
             nb = ewtasks*j + i
             write(6, '(i5)', advance='no') block_to_task(nb)
          enddo
          print*, ' '
       enddo
    endif

    ! Assign grid info for the local task (nt = this_rank)
    ! This info can be copied directly from the block arrays computed above.
    nt = this_rank
    nb = task_to_block(nt)

    ewrank = ewrank_block(nb)
    global_col_offset = global_col_offset_block(nb)
    ewlb = ewlb_block(nb)
    ewub = ewub_block(nb)
    own_ewn = own_ewn_block(nb)
    local_ewn = local_ewn_block(nb)
    ewn = ewn_block(nb)

    nsrank = nsrank_block(nb)
    global_row_offset = global_row_offset_block(nb)
    nslb = nslb_block(nb)
    nsub = nsub_block(nb)
    own_nsn = own_nsn_block(nb)
    local_nsn = local_nsn_block(nb)
    nsn = nsn_block(nb)

    ! Assign W, E, S and N neighbor tasks, based on the block layout.
    ! Note: The halo logic for no_ice BCs is designed so that if (say) west = this_rank,
    !  then no halo values are sent to the west neighbor.

    if (block_is_active(west_block(nb))) then   ! set 'west' to the task that owns west_block(nb)
       west = block_to_task(west_block(nb))
    else  ! set 'west' to the local task
       west = nt
    endif

    if (block_is_active(east_block(nb))) then   ! set 'east' to the task that owns east_block(nb)
       east = block_to_task(east_block(nb))
    else  ! set 'east' to the local task
       east = nt
    endif

    if (block_is_active(south_block(nb))) then   ! set 'south' to the task that owns south_block(nb)
       south = block_to_task(south_block(nb))
    else  ! set 'south' to the local task
       south = nt
    endif

    if (block_is_active(north_block(nb))) then   ! set 'north' to the task that owns north_block(nb)
       north = block_to_task(north_block(nb))
    else  ! set 'north' to the local task
       north = nt
    endif

    if (main_task .and. verbose_active_blocks) then
       print*, ' '
       print*, 'this_rank, nb, west(nb), east(nb), south(nb), north(nb):',  &
                this_rank, nb, west_block(nb), east_block(nb), south_block(nb), north_block(nb)
       print*, 'this_rank, nt, west(nt), east(nt), south(nt), north(nt):',  &
                this_rank, nt, west, east, south, north
    endif

    ! Identify blocks at the corner of the global domain.  Each such block has a locally owned
    ! corner cell where scalars should be set to zero in halo updates.

    southwest_corner = .false.
    southeast_corner = .false.
    northeast_corner = .false.
    northwest_corner = .false.

    if (this_rank > west .and. this_rank > south) then  ! west and south blocks are active
       corner_block = west_block(south_block(nb))
       if (.not.block_is_active(corner_block)) then     ! southwest block is not active
          southwest_corner = .true.
       endif
    endif

    if (this_rank < east .and. this_rank > south) then  ! east and south blocks are active
       corner_block = east_block(south_block(nb))
       if (.not.block_is_active(corner_block)) then     ! southeast block is not active
          southeast_corner = .true.
       endif
    endif

    if (this_rank < east .and. this_rank < north) then  ! east and north blocks are active
       corner_block = east_block(north_block(nb))
       if (.not.block_is_active(corner_block)) then     ! northeat block is not active
          northeast_corner = .true.
       endif
    endif

    if (this_rank > west .and. this_rank < north) then  ! west and north blocks are active
       corner_block = west_block(north_block(nb))
       if (.not.block_is_active(corner_block)) then     ! northwest block is not active
          northwest_corner = .true.
       endif
    endif

    if (verbose_active_blocks) then
!       if (southwest_corner) print*, 'Southwest corner, task =', this_rank
!       if (southeast_corner) print*, 'Southeast corner, task =', this_rank
!       if (northeast_corner) print*, 'Northeast corner, task =', this_rank
!       if (northwest_corner) print*, 'Northwest corner, task =', this_rank
    endif

    ! Set the limits of locally owned vertices on the staggered grid.
    ! For periodic BC, staggered_ilo = staggered_jlo = staggered_lhalo+1.
    ! For no_ice BC, we use the same values as for periodic.
    ! This means that vertices along the southern and western global boundaries
    !  are not considered to be locally owned by any task.  This is acceptable since
    !  there is no ice adjacent to the boundary, and thus (u,v) = 0 for boundary vertices.

    staggered_ilo = staggered_lhalo+1
    staggered_ihi = ewn - 1 - staggered_uhalo
    staggered_jlo = staggered_lhalo+1
    staggered_jhi = nsn - 1 - staggered_uhalo

    ! Check that we have not split up the problem too much.  We do not want halos overlapping in either dimension.
    ! local_* - lhalo - uhalo is the actual number of non-halo cells on a processor.
    if ((local_nsn - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(*,*) "NS halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    if ((local_ewn  - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(*,*) "EW halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    if (verbose_active_blocks) then
    call parallel_barrier
!       print*, 'task, west, east, south, north:', this_rank, west, east, south, north
!       print*, 'task, SW, SE, NE, NW:', this_rank, &
!            southwest_corner, southeast_corner, northwest_corner, northeast_corner
    endif

    ! Uncomment to print grid geometry
!    write(*,*) " "
!    write(*,*) "Process ", this_rank, " Total = ", tasks, " ewtasks = ", ewtasks, " nstasks = ", nstasks
!    write(*,*) "Process ", this_rank, " ewrank = ", ewrank, " nsrank = ", nsrank
!    write(*,*) "Process ", this_rank, " l_ewn = ", local_ewn, " o_ewn = ", own_ewn
!    write(*,*) "Process ", this_rank, " l_nsn = ", local_nsn, " o_nsn = ", own_nsn
!    write(*,*) "Process ", this_rank, " ewlb = ", ewlb, " ewub = ", ewub
!    write(*,*) "Process ", this_rank, " nslb = ", nslb, " nsub = ", nsub
!    write(*,*) "Process ", this_rank, " east = ", east, " west = ", west
!    write(*,*) "Process ", this_rank, " north = ", north, " south = ", south
!    write(*,*) "Process ", this_rank, " ew_vars = ", own_ewn, " ns_vars = ", own_nsn
!    write(*,*) "Process ", this_rank, " global_col_offset = ", global_col_offset, &
!                                      " global_row_offset = ", global_row_offset

    call distributed_print_grid(own_ewn, own_nsn)

    end associate

  end subroutine distributed_grid_active_blocks

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

  subroutine distributed_print_grid(l_ewn,l_nsn)

    ! Gathers and prints the overall grid layout by processor counts.
    use mpi_mod
    implicit none

    integer :: l_ewn, l_nsn

    integer :: i,j,curr_count
    integer,dimension(2) :: mybounds
    integer,dimension(:,:),allocatable :: bounds

    ! begin
    mybounds(1) = l_ewn
    mybounds(2) = l_nsn

    if (main_task) then
       allocate(bounds(2,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          if (bounds(1,i) .ne. -1) then
             ! total up number of processors with matching distribution
             curr_count = 1
             do j = i+1,tasks
                if ((bounds(1,i) .eq. bounds(1,j)) .and. (bounds(2,i) .eq. bounds(2,j))) then
                   ! if matching current distribution, increment counter
                   curr_count = curr_count + 1
                   bounds(1,j) = -1  ! mark so not counted later
                   bounds(2,j) = -1
                endif
             enddo
             write(*,*) "Layout(EW,NS) = ", bounds(1,i), bounds(2,i), " total procs = ", curr_count
          endif
       end do
    end if
    ! automatic deallocation

  end subroutine distributed_print_grid

!=======================================================================

  ! subroutines belonging to the distributed_print interface

  subroutine distributed_print_integer_2d(name, values, parallel)

    use mpi_mod
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: values
    type(parallel_type) :: parallel

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: recvbuf
    integer,dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_print does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(&
            global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
            global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_int(sendbuf,size(sendbuf),mpi_integer,&
       recvbuf,recvcounts,displs,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,1)<local_ewn) then
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_print_integer_2d


  subroutine distributed_print_real8_2d(name, values, parallel)

    use mpi_mod
    implicit none
    character(*) :: name
    real(dp),dimension(:,:) :: values
    type(parallel_type) :: parallel

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_print does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(&
            global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
            global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:) = 0.0d0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,1)<local_ewn) then
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_print_real8_2d


  subroutine distributed_print_real8_3d(name, values, parallel)

    use mpi_mod
    implicit none
    character(*) :: name
    real(dp),dimension(:,:,:) :: values
    type(parallel_type) :: parallel

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:,:),allocatable :: global_values,sendbuf

    ! begin

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_print does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(size(values,1),minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(size(values,1), &
            global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
            global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:,:) = 0.0d0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)*size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(size(values,1),mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:,:) = values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    sendbuf(:,mybounds(1):mybounds(2),mybounds(3):mybounds(4)) = sendbuf(:,mybounds(1):mybounds(2),mybounds(3):mybounds(4))
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(:,bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/size(values,1),bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,2)<local_ewn) then
          do j = lbound(global_values,3),ubound(global_values,3)
             do i = lbound(global_values,2),ubound(global_values,2)
                write(u,'(2i6,100g15.5e3)') j,i,global_values(:,i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,3),ubound(global_values,3)
             do i = lbound(global_values,2),ubound(global_values,2)
                write(u,'(2i6,100g15.5e3)') j,i,global_values(:,i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if

    end associate
    ! automatic deallocation

  end subroutine distributed_print_real8_3d

!=======================================================================

  ! functions belonging to the distributed_put_var interface

  function distributed_put_var_integer_2d(ncid, varid, values, parallel, start)

    use mpi_mod
    implicit none
    integer :: distributed_put_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values
    type(parallel_type) :: parallel

    type(bounds_info_type) :: bounds_info
    integer :: i,ierror
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: recvbuf
    integer,dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    associate(  &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    bounds_info = get_bounds_info(size(values,1), parallel)

    mybounds(1) = bounds_info%mybounds_ew_lb
    mybounds(2) = bounds_info%mybounds_ew_ub
    mybounds(3) = bounds_info%mybounds_ns_lb
    mybounds(4) = bounds_info%mybounds_ns_ub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(&
            global_minval_ewlb+lhalo:global_maxval_ewub-uhalo, &
            global_minval_nslb+lhalo:global_maxval_nsub-uhalo))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(bounds_info%ilo:bounds_info%ihi, &
                          bounds_info%jlo:bounds_info%jhi)
    call fc_gatherv_int(sendbuf,size(sendbuf),mpi_integer,&
       recvbuf,recvcounts,displs,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_integer_2d = nf90_put_var(ncid,varid,&
            global_values(1:bounds_info%global_ewn, 1:bounds_info%global_nsn),start)
    end if
    call broadcast(distributed_put_var_integer_2d)

    end associate
    !automatic deallocation

  end function distributed_put_var_integer_2d


  function distributed_put_var_real4_1d(ncid, varid, values, parallel, start)

    use mpi_mod
    use netcdf
    implicit none
    integer :: distributed_put_var_real4_1d,ncid,varid
    real(sp),dimension(:) :: values
    type(parallel_type) :: parallel
    integer,dimension(:),optional :: start

    integer :: i,ierror,myn,status,x0id,x1id,y0id,y1id
    integer,dimension(2) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(sp),dimension(:),allocatable :: global_values,recvbuf

    ! begin

    associate(  &
         global_ewn  => parallel%global_ewn,   &
         global_nsn  => parallel%global_nsn,   &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x0id)
    call broadcast(x1id)
    call broadcast(y0id)
    call broadcast(y1id)
    if (varid==x0id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub-1
       myn = global_ewn-1
    else if (varid==x1id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub
       myn = global_ewn
    else if (varid==y0id) then
       mybounds(1) = nslb
       mybounds(2) = nsub-1
       myn = global_nsn-1
    else if (varid==y1id) then
       mybounds(1) = nslb
       mybounds(2) = nsub
       myn = global_nsn
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:))))
       if (varid==x0id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub-1))
       elseif (varid==x1id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub))
       elseif (varid==y0id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub-1))
       elseif (varid==y1id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub))
       endif
       global_values(:) = 0.0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = bounds(2,:)-bounds(1,:)+1
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    call fc_gatherv_real4(values,size(values),mpi_real4,&
       recvbuf,recvcounts,displs,mpi_real4,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i)) = &
               recvbuf(displs(i)+1:displs(i+1))
       end do
       if (present(start)) then
         distributed_put_var_real4_1d = &
            nf90_put_var(ncid,varid,global_values(1:myn),start)
       else
         distributed_put_var_real4_1d = &
            nf90_put_var(ncid,varid,global_values(1:myn))
       endif
    end if
    call broadcast(distributed_put_var_real4_1d)

    end associate
    !automatic deallocation

  end function distributed_put_var_real4_1d


  function distributed_put_var_real4_2d(ncid, varid, values, parallel, start)

    use mpi_mod
    implicit none
    integer :: distributed_put_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(sp),dimension(:,:) :: values
    type(parallel_type) :: parallel

    type(bounds_info_type) :: bounds_info
    integer :: i,ierror
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(sp),dimension(:),allocatable :: recvbuf
    real(sp),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    associate(  &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    bounds_info = get_bounds_info(size(values,1), parallel)

    mybounds(1) = bounds_info%mybounds_ew_lb
    mybounds(2) = bounds_info%mybounds_ew_ub
    mybounds(3) = bounds_info%mybounds_ns_lb
    mybounds(4) = bounds_info%mybounds_ns_ub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(global_minval_ewlb:global_maxval_ewub, &
                              global_minval_nslb:global_maxval_nsub))
       global_values(:,:) = 0.0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(bounds_info%ilo:bounds_info%ihi, &
                          bounds_info%jlo:bounds_info%jhi)
    call fc_gatherv_real4(sendbuf,size(sendbuf),mpi_real4,&
       recvbuf,recvcounts,displs,mpi_real4,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_real4_2d = nf90_put_var(ncid,varid,&
            global_values(1:bounds_info%global_ewn, 1:bounds_info%global_nsn),start)
    end if
    call broadcast(distributed_put_var_real4_2d)

    end associate
    !automatic deallocation

  end function distributed_put_var_real4_2d


  function distributed_put_var_real8_1d(ncid, varid, values, parallel, start)

    use mpi_mod
    use netcdf
    implicit none
    integer :: distributed_put_var_real8_1d,ncid,varid
    real(dp),dimension(:) :: values
    type(parallel_type) :: parallel
    integer,dimension(:),optional :: start

    integer :: i,ierror,myn,status,x0id,x1id,y0id,y1id

    integer,dimension(2) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(dp),dimension(:),allocatable :: global_values,recvbuf

    ! begin

    associate(  &
         global_ewn  => parallel%global_ewn,   &
         global_nsn  => parallel%global_nsn,   &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x0id)
    call broadcast(x1id)
    call broadcast(y0id)
    call broadcast(y1id)
    if (varid==x0id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub-1
       myn = global_ewn-1
    else if (varid==x1id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub
       myn = global_ewn
    else if (varid==y0id) then
       mybounds(1) = nslb
       mybounds(2) = nsub-1
       myn = global_nsn-1
    else if (varid==y1id) then
       mybounds(1) = nslb
       mybounds(2) = nsub
       myn = global_nsn
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:))))
       if (varid==x0id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub-1))
       elseif (varid==x1id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub))
       elseif (varid==y0id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub-1))
       elseif (varid==y1id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub))
       endif
       global_values(:) = 0.0d0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = bounds(2,:)-bounds(1,:)+1
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    call fc_gatherv_real8(values,size(values),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i)) = &
               recvbuf(displs(i)+1:displs(i+1))
       end do
       if (present(start)) then
         distributed_put_var_real8_1d = &
            nf90_put_var(ncid,varid,global_values(1:myn),start)
       else
         distributed_put_var_real8_1d = &
            nf90_put_var(ncid,varid,global_values(1:myn))
       endif
    end if
    call broadcast(distributed_put_var_real8_1d)

    end associate
    !automatic deallocation

  end function distributed_put_var_real8_1d


  function distributed_put_var_real8_2d(ncid, varid, values, parallel, start)

    use mpi_mod
    implicit none
    integer :: distributed_put_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:,:) :: values
    type(parallel_type) :: parallel

    type(bounds_info_type) :: bounds_info
    integer :: i,ierror
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    associate(  &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    bounds_info = get_bounds_info(size(values,1), parallel)

    mybounds(1) = bounds_info%mybounds_ew_lb
    mybounds(2) = bounds_info%mybounds_ew_ub
    mybounds(3) = bounds_info%mybounds_ns_lb
    mybounds(4) = bounds_info%mybounds_ns_ub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:))))
       allocate(global_values(global_minval_ewlb:global_maxval_ewub, &
                              global_minval_nslb:global_maxval_nsub))
       global_values(:,:) = 0.0d0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(bounds_info%ilo:bounds_info%ihi, &
                          bounds_info%jlo:bounds_info%jhi)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_real8_2d = nf90_put_var(ncid,varid,&
            global_values(1:bounds_info%global_ewn, 1:bounds_info%global_nsn),start)
    end if
    call broadcast(distributed_put_var_real8_2d)

    end associate
    !automatic deallocation

  end function distributed_put_var_real8_2d


  function distributed_put_var_real8_3d(ncid, varid, values, parallel, start)

    use mpi_mod
    implicit none
    integer :: distributed_put_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(dp),dimension(:,:,:) :: values
    type(parallel_type) :: parallel

    type(bounds_info_type) :: bounds_info
    integer :: i,ierror,nz
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(dp),dimension(:),allocatable :: recvbuf
    real(dp),dimension(:,:,:),allocatable :: global_values,sendbuf

    ! begin

    !TODO: Add these global_minval/maxval terms to the bounds_info type?
    associate(  &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub  &
         )

    nz = size(values,3)
    bounds_info = get_bounds_info(size(values,1), parallel)

    mybounds(1) = bounds_info%mybounds_ew_lb
    mybounds(2) = bounds_info%mybounds_ew_ub
    mybounds(3) = bounds_info%mybounds_ns_lb
    mybounds(4) = bounds_info%mybounds_ns_ub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,4,mpi_integer,bounds,4,&
       mpi_integer,main_rank,comm)
    if (main_task) then
       !WHL - See comments above on allocating the global_values array
!!       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
!!            minval(bounds(3,:)):maxval(bounds(4,:)),nz))
       allocate(global_values(global_minval_ewlb:global_maxval_ewub, &
                              global_minval_nslb:global_maxval_nsub, nz))
       global_values(:,:,:) = 0.0d0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)&
            *nz
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4),nz))
    sendbuf(:,:,:) = values(bounds_info%ilo:bounds_info%ihi, &
                            bounds_info%jlo:bounds_info%jhi, &
                            :)
    call fc_gatherv_real8(sendbuf,size(sendbuf),mpi_real8,&
       recvbuf,recvcounts,displs,mpi_real8,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i),:) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1,nz/))
       end do
       distributed_put_var_real8_3d = nf90_put_var(ncid,varid,&
            global_values(1:bounds_info%global_ewn, 1:bounds_info%global_nsn ,:),start)
    end if
    call broadcast(distributed_put_var_real8_3d)

    end associate
    !automatic deallocation

  end function distributed_put_var_real8_3d

!=======================================================================

  ! subroutines belonging to the distributed_scatter_var interface

  subroutine distributed_scatter_var_integer_2d(values, global_values, parallel)

    ! Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! This subroutine expects global_values to be allocated on all tasks.
    ! It can be allocated with zero size on tasks other than main_task.
    use mpi_mod
    implicit none
    integer,dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    integer,dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:),allocatable :: sendbuf
    integer,dimension(:,:),allocatable :: recvbuf

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub          &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_scatter does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if

    allocate(recvbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_integer,&
         recvbuf,size(recvbuf),mpi_integer,main_rank,comm,ierror)
    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)

    end associate
    deallocate(global_values)   ! TODO - Is this deallocation necessary, here and below?
    ! automatic deallocation

  end subroutine distributed_scatter_var_integer_2d


  subroutine distributed_scatter_var_logical_2d(values, global_values, parallel)

    ! Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! This subroutine expects global_values to be allocated on all tasks.
    ! It can be allocated with zero size on tasks other than main_task.
    use mpi_mod
    implicit none
    logical,dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    logical,dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    logical,dimension(:),allocatable :: sendbuf
    logical,dimension(:,:),allocatable :: recvbuf

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub          &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_scatter does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_logical,&
         recvbuf,size(recvbuf),mpi_logical,main_rank,comm,ierror)
    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)

    end associate
    deallocate(global_values)
    ! automatic deallocation

  end subroutine distributed_scatter_var_logical_2d


  subroutine distributed_scatter_var_real4_2d(values, global_values, parallel)

    ! Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! This subroutine expects global_values to be allocated on all tasks.
    ! It can be allocated with zero size on tasks other than main_task.
    use mpi_mod
    implicit none
    real(sp),dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    real(sp),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(sp),dimension(:),allocatable :: sendbuf
    real(sp),dimension(:,:),allocatable :: recvbuf

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub          &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_scatter does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                     d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                     (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         recvbuf,size(recvbuf),mpi_real4,main_rank,comm,ierror)
    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)

    end associate
    deallocate(global_values)
    ! automatic deallocation

  end subroutine distributed_scatter_var_real4_2d


  subroutine distributed_scatter_var_real4_3d(values, global_values, parallel)

    ! Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! This subroutine expects global_values to be allocated on all tasks.
    ! It can be allocated with zero size on tasks other than main_task.
    use mpi_mod
    implicit none
    real(sp),dimension(:,:,:),intent(inout) :: values  ! populated from values on main_task
    real(sp),dimension(:,:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(sp),dimension(:),allocatable :: sendbuf
    real(sp),dimension(:,:,:),allocatable :: recvbuf

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub          &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_scatter does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)*size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(:,&
                                   d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(size(values,1),&
                     d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         recvbuf,size(recvbuf),mpi_real4,main_rank,comm,ierror)
    values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:,:)

    end associate
    deallocate(global_values)
    ! automatic deallocation

  end subroutine distributed_scatter_var_real4_3d


  subroutine distributed_scatter_var_real8_2d(values, global_values, parallel)

    ! Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! This subroutine expects global_values to be allocated on all tasks.
    ! It can be allocated with zero size on tasks other than main_task.
    use mpi_mod
    implicit none
    real(dp),dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(dp),dimension(:),allocatable :: sendbuf
    real(dp),dimension(:,:),allocatable :: recvbuf

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub          &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_scatter does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)

    end associate
    deallocate(global_values)
    ! automatic deallocation

  end subroutine distributed_scatter_var_real8_2d


  subroutine distributed_scatter_var_real8_3d(values, global_values, parallel, deallocflag)

    ! Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    ! This subroutine expects global_values to be allocated on all tasks.
    ! It can be allocated with zero size on tasks other than main_task.
    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:),intent(inout) :: values  ! populated from values on main_task
    real(dp),dimension(:,:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    type(parallel_type) :: parallel
    logical,optional :: deallocflag   ! TODO - Is this flag needed?  Not in other scatter subroutines.

    logical :: deallocmem
    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(dp),dimension(:),allocatable :: sendbuf
    real(dp),dimension(:,:,:),allocatable :: recvbuf

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub          &
         )

    if (uhalo==0 .and. size(values,1)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "distributed_scatter does not currently work for"
       write(*,*) "variables on the staggered grid when uhalo=0"
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (present(deallocflag)) then
       deallocmem = deallocflag
    else
       deallocmem = .true.
    endif

    ! first time
    if (.not. allocated(d_gs_bounds)) then
       if (main_task) then
          allocate(d_gs_bounds(4,tasks))
       else
          allocate(d_gs_bounds(1,1))
       endif
    endif

    ! Note: Originally, d_gs_mybounds and d_gs_bounds were computed only once.
    !       Now, they are recomputed each time to allow for multiple instances with different ewlb, etc.
    d_gs_mybounds(1) = ewlb+lhalo
    d_gs_mybounds(2) = ewub-uhalo
    d_gs_mybounds(3) = nslb+lhalo
    d_gs_mybounds(4) = nsub-uhalo
    call fc_gather_int(d_gs_mybounds,4,mpi_integer,d_gs_bounds,4,&
         mpi_integer,main_rank,comm)

    if (main_task) then
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (d_gs_bounds(2,:)-d_gs_bounds(1,:)+1)&
                      *(d_gs_bounds(4,:)-d_gs_bounds(3,:)+1)&
                      *size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))

       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(:,&
                                   d_gs_bounds(1,i):d_gs_bounds(2,i),&
                                   d_gs_bounds(3,i):d_gs_bounds(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    allocate(recvbuf(size(values,1),&
                     d_gs_mybounds(1):d_gs_mybounds(2),&
                     d_gs_mybounds(3):d_gs_mybounds(4)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:,:)

    end associate
    if (deallocmem) deallocate(global_values)
    ! automatic deallocation

  end subroutine distributed_scatter_var_real8_3d

!=======================================================================

  ! subroutines belonging to the distributed_scatter_var_row interface

  subroutine distributed_scatter_var_row_real8_2d(values, global_values, parallel)

    ! Scatter data to a row of tasks from the main task for that row.
    ! Based on distributed_scatter_var_real8_2d.
    ! Note: The first index represents a data dimension that is the same on each task,
    !        whose size generally is less than own_ewn.
    !       The second index represents the north-south dimension, and is assumed
    !        to have size own_nsn (i.e., the data extend over locally owned cells only).
    ! values = local portion of distributed variable
    ! global_values = reference to allocatable array in which main_task_row holds the variable.
    ! global_values is deallocated at the end.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:),intent(inout) :: values  ! populated from values on main_task_row
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task_row
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(dp),dimension(:),allocatable :: sendbuf
    real(dp),dimension(:,:),allocatable :: recvbuf

    integer,dimension(4) :: d_gs_mybounds_row
    integer,dimension(:,:),allocatable,save :: d_gs_bounds_row

    associate(  &
         comm_row      => parallel%comm_row,       &
         tasks_row     => parallel%tasks_row,      &
         this_rank_row => parallel%this_rank_row,  &
         main_rank_row => parallel%main_rank_row,  &
         main_task_row => parallel%main_task_row,  &
         own_nsn       => parallel%own_nsn         &
         )

    if (size(values,2) /= own_nsn) then
       ! Note: Removing this restriction would require some recoding below.
       write(*,*) "ERROR: distributed_scatter_var_row requires N-S array size of own_nsn"
       call parallel_stop(__FILE__, __LINE__)
    end if

    d_gs_mybounds_row(1) = this_rank_row*size(values,1) + 1
    d_gs_mybounds_row(2) = (this_rank_row+1)*size(values,1)
    d_gs_mybounds_row(3) = 1
    d_gs_mybounds_row(4) = own_nsn

    if (allocated(d_gs_bounds_row)) then
       ! do nothing; d_gs_bounds_row already computed
    else   ! first time
       if (main_task_row) then
          allocate(d_gs_bounds_row(4,tasks_row))
       else
          allocate(d_gs_bounds_row(1,1))
       endif

       call fc_gather_int(d_gs_mybounds_row,4,mpi_integer,d_gs_bounds_row,4,&
            mpi_integer,main_rank_row,comm_row)
    endif

    if (main_task_row) then
       allocate(displs(tasks_row+1))
       allocate(sendcounts(tasks_row))
       sendcounts(:) = (d_gs_bounds_row(2,:)-d_gs_bounds_row(1,:)+1)&
                      *(d_gs_bounds_row(4,:)-d_gs_bounds_row(3,:)+1)
       displs(1) = 0
       do i = 1,tasks_row
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks_row+1)))
       do i = 1,tasks_row
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds_row(1,i):d_gs_bounds_row(2,i),&
                                   d_gs_bounds_row(3,i):d_gs_bounds_row(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if

    !Note: Would need to uncomment the following if recvbuf were not identical
    !      to the output values array
!!    allocate(recvbuf(d_gs_mybounds_row(1):d_gs_mybounds_row(2),&
!!                     d_gs_mybounds_row(3):d_gs_mybounds_row(4)))
!!    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
!!         recvbuf,size(recvbuf),mpi_real8,main_rank_row,comm_row,ierror)
!!    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         values,size(values),mpi_real8,main_rank_row,comm_row,ierror)

    end associate
    ! automatic deallocation

  end subroutine distributed_scatter_var_row_real8_2d

!=======================================================================

  ! subroutines belonging to the distributed_scatter_var_col interface

  subroutine distributed_scatter_var_col_real8_2d(values, global_values, parallel)

    ! Scatter data to a column of tasks from the main task for that column
    ! Based on distributed_scatter_var_real8_2d.
    ! Note: The first index represents a data dimension that is the same on each task,
    !        whose size generally is less than own_nsn.
    !       The second index represents the east-west dimension, and is assumed
    !        to have size own_ewn (i.e., the data extend over locally owned cells only).
    ! values = local portion of distributed variable
    ! global_values = reference to allocatable array in which main_task_col holds the variable.
    ! global_values is deallocated at the end.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:),intent(inout) :: values  ! populated from values on main_task_col
    real(dp),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task_col
    type(parallel_type) :: parallel

    integer :: i,ierror,j,k
    integer,dimension(:),allocatable :: displs,sendcounts
    real(dp),dimension(:),allocatable :: sendbuf
    real(dp),dimension(:,:),allocatable :: recvbuf
    integer,dimension(4) :: d_gs_mybounds_col
    integer,dimension(:,:),allocatable,save :: d_gs_bounds_col

    associate(  &
         comm_col      => parallel%comm_col,       &
         tasks_col     => parallel%tasks_col,      &
         this_rank_col => parallel%this_rank_col,  &
         main_rank_col => parallel%main_rank_col,  &
         main_task_col => parallel%main_task_col,  &
         own_ewn       => parallel%own_ewn         &
         )

    if (size(values,2) /= own_ewn) then
       ! Note: Removing this restriction would require some recoding below.
       write(*,*) "ERROR: distributed_scatter_var_col requires E-W array size of own_nsn"
       call parallel_stop(__FILE__, __LINE__)
    end if

    d_gs_mybounds_col(1) = this_rank_col*size(values,1) + 1
    d_gs_mybounds_col(2) = (this_rank_col+1)*size(values,1)
    d_gs_mybounds_col(3) = 1
    d_gs_mybounds_col(4) = own_ewn

    if (allocated(d_gs_bounds_col)) then
       ! do nothing; d_gs_bounds_col already computed
    else   ! first time
       if (main_task_col) then
          allocate(d_gs_bounds_col(4,tasks_col))
       else
          allocate(d_gs_bounds_col(1,1))
       endif

       call fc_gather_int(d_gs_mybounds_col,4,mpi_integer,d_gs_bounds_col,4,&
            mpi_integer,main_rank_col,comm_col)
    endif

    if (main_task_col) then
       allocate(displs(tasks_col+1))
       allocate(sendcounts(tasks_col))
       sendcounts(:) = (d_gs_bounds_col(2,:)-d_gs_bounds_col(1,:)+1)&
                      *(d_gs_bounds_col(4,:)-d_gs_bounds_col(3,:)+1)
       displs(1) = 0
       do i = 1,tasks_col
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks_col+1)))
       do i = 1,tasks_col
          sendbuf(displs(i)+1:displs(i+1)) = &
             reshape(global_values(d_gs_bounds_col(1,i):d_gs_bounds_col(2,i),&
                                   d_gs_bounds_col(3,i):d_gs_bounds_col(4,i)),&
                                   (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if

    !Note: Would need to uncomment the following if recvbuf were not identical
    !      to the output values array
!!    allocate(recvbuf(d_gs_mybounds_col(1):d_gs_mybounds_col(2),&
!!                     d_gs_mybounds_col(3):d_gs_mybounds_col(4)))
!!    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
!!         recvbuf,size(recvbuf),mpi_real8,main_rank_col,comm_col,ierror)
!!    values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo) = recvbuf(:,:)
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         values,size(values),mpi_real8,main_rank_col,comm_col,ierror)

    end associate
    ! automatic deallocation

  end subroutine distributed_scatter_var_col_real8_2d

!=======================================================================

  function get_bounds_info(array_ew_size, parallel) result(bounds_info)

    ! Determine information on the local & global bounds of an array
    ! This is used to distinguish between arrays on the staggered vs. unstaggered grids
    
    integer, intent(in) :: array_ew_size  ! size of the array of interest in the EW direction
    type(parallel_type) :: parallel
    type(bounds_info_type) :: bounds_info

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

    if (array_ew_size == local_ewn) then
       bounds_info%global_ewn = global_ewn
       bounds_info%global_nsn = global_nsn

       bounds_info%mybounds_ew_lb = ewlb + lhalo
       bounds_info%mybounds_ew_ub = ewub - uhalo
       bounds_info%mybounds_ns_lb = nslb + lhalo
       bounds_info%mybounds_ns_ub = nsub - uhalo

       bounds_info%ilo = 1 + lhalo
       bounds_info%ihi = local_ewn - uhalo
       bounds_info%jlo = 1 + lhalo
       bounds_info%jhi = local_nsn - uhalo

    else if (array_ew_size == (local_ewn - 1)) then
       bounds_info%global_ewn = global_ewn - 1
       bounds_info%global_nsn = global_nsn - 1

       bounds_info%mybounds_ew_lb = ewlb + staggered_lhalo
       bounds_info%mybounds_ew_ub = (ewub - 1) - staggered_uhalo
       bounds_info%mybounds_ns_lb = nslb + staggered_lhalo
       bounds_info%mybounds_ns_ub = (nsub - 1) - staggered_uhalo

       bounds_info%ilo = 1 + staggered_lhalo
       bounds_info%ihi = (local_ewn - 1) - staggered_uhalo
       bounds_info%jlo = 1 + staggered_lhalo
       bounds_info%jhi = (local_nsn - 1) - staggered_uhalo

    else
       call parallel_stop(__FILE__, __LINE__)
    end if

    end associate

  end function get_bounds_info

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

  function parallel_close(ncid)

    implicit none
    integer :: ncid

    integer :: parallel_close

    ! begin
    if (main_task) parallel_close = nf90_close(ncid)
    call broadcast(parallel_close)

  end function parallel_close

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
       write(*,*) "Unexpected size for input_with_halo: ", &
            size(input_with_halo,1), size(input_with_halo,2)
       write(*,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_no_halo,1) /= own_ewn .or. size(output_no_halo,2) /= own_nsn) then
       write(*,*) "Unexpected size for output_no_halo: ", &
            size(output_no_halo,1), size(output_no_halo,2)
       write(*,*) "Expected size is: ", own_ewn, own_nsn
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
       write(*,*) "Unexpected size for input_with_halo: ", &
            size(input_with_halo,1), size(input_with_halo,2)
       write(*,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_no_halo,1) /= own_ewn .or. size(output_no_halo,2) /= own_nsn) then
       write(*,*) "Unexpected size for output_no_halo: ", &
            size(output_no_halo,1), size(output_no_halo,2)
       write(*,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_no_halo(1:own_ewn, 1:own_nsn) = &
         input_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo)

    end associate

  end subroutine parallel_convert_haloed_to_nonhaloed_real8_2d

!=======================================================================

  ! subroutines belonging to the parallel_convert_nonhaloed_to_haloed interface

  subroutine parallel_convert_nonhaloed_to_haloed_real4_2d(input_no_halo, output_with_halo, parallel)

    ! Given an input array without halo cells, return an output array with halo cells
    real(sp),dimension(:,:), intent(in)  :: input_no_halo
    real(sp),dimension(:,:), intent(out) :: output_with_halo
    type(parallel_type) :: parallel
    
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn       &
         )

    if (size(input_no_halo,1) /= own_ewn .or. size(input_no_halo,2) /= own_nsn) then
       write(*,*) "Unexpected size for input_no_halo: ", &
            size(input_no_halo,1), size(input_no_halo,2)
       write(*,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_with_halo,1) /= local_ewn .or. size(output_with_halo,2) /= local_nsn) then
       write(*,*) "Unexpected size for output_with_halo: ", &
            size(output_with_halo,1), size(output_with_halo,2)
       write(*,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = &
         input_no_halo(1:own_ewn, 1:own_nsn)

    call parallel_halo(output_with_halo, parallel)
    
    end associate

  end subroutine parallel_convert_nonhaloed_to_haloed_real4_2d


  subroutine parallel_convert_nonhaloed_to_haloed_real8_2d(input_no_halo, output_with_halo, parallel)

    ! Given an input array without halo cells, return an output array with halo cells
    real(dp),dimension(:,:), intent(in)  :: input_no_halo
    real(dp),dimension(:,:), intent(out) :: output_with_halo
    type(parallel_type) :: parallel
    
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn       &
         )

    if (size(input_no_halo,1) /= own_ewn .or. size(input_no_halo,2) /= own_nsn) then
       write(*,*) "Unexpected size for input_no_halo: ", &
            size(input_no_halo,1), size(input_no_halo,2)
       write(*,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_with_halo,1) /= local_ewn .or. size(output_with_halo,2) /= local_nsn) then
       write(*,*) "Unexpected size for output_with_halo: ", &
            size(output_with_halo,1), size(output_with_halo,2)
       write(*,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = &
         input_no_halo(1:own_ewn, 1:own_nsn)

    call parallel_halo(output_with_halo, parallel)
    
    end associate

  end subroutine parallel_convert_nonhaloed_to_haloed_real8_2d

!=======================================================================

  function parallel_create(path, cmode, ncid)

    implicit none
    integer :: cmode, ncid, parallel_create

    character(len=*) :: path

    ! begin
    if (main_task) parallel_create = nf90_create(path,cmode,ncid)
    call broadcast(parallel_create)
    call broadcast(ncid)

  end function parallel_create

!=======================================================================

  !WHL - Added the next two subroutines to create communicators for rows and columns of tasks.
  !      These are useful for a tridiagonal solve in 1D along a global row or column.
  !      Cf. parallel_initialise and parallel_set_info

  subroutine parallel_create_comm_row(comm, parallel)

    use mpi_mod
    implicit none
    integer, intent(in) :: comm          ! global communicator
    type(parallel_type) :: parallel

    integer, parameter :: my_main_rank_row = 0
    integer :: ierror

    associate(  &
         nsrank        => parallel%nsrank,        &
         comm_row      => parallel%comm_row,      &
         tasks_row     => parallel%tasks_row,     &
         this_rank_row => parallel%this_rank_row, &
         main_rank_row => parallel%main_rank_row, &
         main_task_row => parallel%main_task_row  &
         )

    ! Create a communicator (comm_row) that groups rows according to their value of nsrank,
    !  and assign a main_task to each row.

    call mpi_comm_split(comm, nsrank, this_rank, comm_row, ierror)
    call mpi_comm_size(comm_row, tasks_row, ierror)
    call mpi_comm_rank(comm_row, this_rank_row, ierror)
    main_rank_row = my_main_rank_row
    main_task_row = (this_rank_row==main_rank_row)

    if (main_task_row) print*, 'Create comm_row: this_rank, tasks_row, this_rank_row, main_task_row =', &
         this_rank, tasks_row, this_rank_row, main_task_row

    end associate

  end subroutine parallel_create_comm_row

!=======================================================================

  subroutine parallel_create_comm_col(comm, parallel)

    use mpi_mod
    implicit none
    integer, intent(in) :: comm          ! global communicator
    type(parallel_type) :: parallel

    integer, parameter :: my_main_rank_col = 0
    integer :: ierror

    associate(  &
         ewrank        => parallel%ewrank,        &
         comm_col      => parallel%comm_col,      &
         tasks_col     => parallel%tasks_col,     &
         this_rank_col => parallel%this_rank_col, &
         main_rank_col => parallel%main_rank_col, &
         main_task_col => parallel%main_task_col  &
         )

    ! Create a communicator (comm_col) that groups columns according to their value of ewrank,
    !  and assign a main_task to each column.

    call mpi_comm_split(comm, ewrank, this_rank, comm_col, ierror)
    call mpi_comm_size(comm_col, tasks_col, ierror)
    call mpi_comm_rank(comm_col, this_rank_col, ierror)
    main_rank_col = my_main_rank_col
    main_task_col = (this_rank_col==main_rank_col)

    if (main_task_col) print*, 'Create comm_col: this_rank, tasks_col, this_rank_col, main_task_col =', &
         this_rank, tasks_col, this_rank_col, main_task_col

    end associate

  end subroutine parallel_create_comm_col

!=======================================================================

  function parallel_def_dim(ncid, name, len, dimid)

    use netcdf
    implicit none
    integer :: dimid, len, ncid, parallel_def_dim
    type(parallel_type) :: parallel

    character(len=*) :: name

    ! begin
    if (main_task) parallel_def_dim = nf90_def_dim(ncid,name,len,dimid)
    call broadcast(parallel_def_dim)
    call broadcast(dimid)

  end function parallel_def_dim

!=======================================================================

  ! functions belonging to the parallel_def_var interface

  function parallel_def_var_dimids(ncid, name, xtype, dimids, varid)

    implicit none
    integer :: ncid, parallel_def_var_dimids, varid, xtype
    integer,dimension(:) :: dimids
    character(len=*) :: name

    ! begin
    if (main_task) parallel_def_var_dimids = &
         nf90_def_var(ncid,name,xtype,dimids,varid)
    call broadcast(parallel_def_var_dimids)
    call broadcast(varid)

  end function parallel_def_var_dimids


  function parallel_def_var_nodimids(ncid, name, xtype, varid)

    implicit none
    integer :: ncid, parallel_def_var_nodimids, varid, xtype
    character(len=*) :: name

    ! begin
    if (main_task) parallel_def_var_nodimids = &
         nf90_def_var(ncid,name,xtype,varid)
    call broadcast(parallel_def_var_nodimids)
    call broadcast(varid)

  end function parallel_def_var_nodimids

!=======================================================================

  function parallel_enddef(ncid)

    implicit none
    integer :: ncid,parallel_enddef

    ! begin
    if (main_task) parallel_enddef = nf90_enddef(ncid)
    call broadcast(parallel_enddef)

  end function parallel_enddef

!=======================================================================

  subroutine parallel_finalise

    use mpi_mod
    implicit none
    integer :: ierror

    ! begin
    call mpi_finalize(ierror)

  end subroutine parallel_finalise

!=======================================================================

  ! functions belonging to the parallel_get_att interface

  function parallel_get_att_character(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_get_att_character,varid
    character(len=*) :: name,values

    ! begin
    if (main_task) parallel_get_att_character = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_character)
    call broadcast(values)

  end function parallel_get_att_character


  function parallel_get_att_real4(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_get_att_real4,varid
    character(len=*) :: name
    real(sp) :: values

    ! begin
    if (main_task) parallel_get_att_real4 = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real4)
    call broadcast(values)

  end function parallel_get_att_real4


  function parallel_get_att_real4_1d(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_get_att_real4_1d, varid
    character(len=*) :: name
    real(sp),dimension(:) :: values

    ! begin
    if (main_task) parallel_get_att_real4_1d = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real4_1d)
    call broadcast(values)

  end function parallel_get_att_real4_1d


  function parallel_get_att_real8(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_get_att_real8, varid
    character(len=*) :: name
    real(dp) :: values

    ! begin
    if (main_task) parallel_get_att_real8 = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real8)
    call broadcast(values)

  end function parallel_get_att_real8


  function parallel_get_att_real8_1d(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_get_att_real8_1d, varid
    character(len=*) :: name
    real(dp),dimension(:) :: values

    ! begin
    if (main_task) parallel_get_att_real8_1d = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real8_1d)
    call broadcast(values)

  end function parallel_get_att_real8_1d

!=======================================================================

  ! functions belonging to the parallel_get_var interface

  !WHL - Added parallel_get_var functions in analogy to parallel_put_var functions.
  !      Similar to distributed_get_var, but they lack a 'start' argument.
  !      The scalar and 1D functions broadcast values from main_task to other tasks.
  !      The 2D functions do not broadcast; they only bring a global array to main_task.

  function parallel_get_var_integer(ncid, varid, values)

    implicit none
    integer :: ncid, parallel_get_var_integer, varid
    integer :: values

    ! begin
    if (main_task) parallel_get_var_integer = nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_integer)
    call broadcast(values)

  end function parallel_get_var_integer


  function parallel_get_var_real4(ncid, varid, values)

    implicit none
    integer :: ncid, parallel_get_var_real4, varid
    real(sp) :: values

    ! begin
    if (main_task) parallel_get_var_real4 = nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_real4)
    call broadcast(values)

  end function parallel_get_var_real4


  function parallel_get_var_real8(ncid, varid, values)

    implicit none
    integer :: ncid,parallel_get_var_real8,varid
    real(dp) :: values

    ! begin
    if (main_task) parallel_get_var_real8 = nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_real8)
    call broadcast(values)

  end function parallel_get_var_real8


  function parallel_get_var_integer_1d(ncid, varid, values)

    implicit none
    integer :: ncid,parallel_get_var_integer_1d,varid
    integer,dimension(:) :: values

    ! begin
    if (main_task) parallel_get_var_integer_1d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_integer_1d)
    call broadcast(values)

  end function parallel_get_var_integer_1d

  function parallel_get_var_real4_1d(ncid, varid, values)

    implicit none
    integer :: ncid,parallel_get_var_real4_1d,varid
    real(sp),dimension(:) :: values

    ! begin
    if (main_task) parallel_get_var_real4_1d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_real4_1d)
    call broadcast(values)

  end function parallel_get_var_real4_1d

  function parallel_get_var_real8_1d(ncid, varid, values)

    implicit none
    integer :: ncid,parallel_get_var_real8_1d,varid
    real(dp),dimension(:) :: values

    ! begin
    if (main_task) parallel_get_var_real8_1d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_real8_1d)
    call broadcast(values)

  end function parallel_get_var_real8_1d


  function parallel_get_var_integer_2d(ncid, varid, values)

    implicit none
    integer :: ncid,parallel_get_var_integer_2d,varid
    integer,dimension(:,:) :: values

    ! begin
    if (main_task) parallel_get_var_integer_2d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_integer_2d)
!!    call broadcast(values)  ! no broadcast subroutine for 2D arrays

  end function parallel_get_var_integer_2d


  function parallel_get_var_real8_2d(ncid, varid, values)

    implicit none
    integer :: ncid,parallel_get_var_real8_2d,varid
    real(dp),dimension(:,:) :: values

    ! begin
    if (main_task) parallel_get_var_real8_2d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_real8_2d)
!!    call broadcast(values)  ! no broadcast subroutine for 2D arrays

  end function parallel_get_var_real8_2d

!=======================================================================

  subroutine parallel_global_edge_mask(global_edge_mask, parallel)

    ! Create a mask = 1 in locally owned cells at the edge of the global domain,
    ! = 0 elsewhere

    integer, dimension(:,:), intent(out) :: global_edge_mask
    type(parallel_type) :: parallel

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south)

    ! Check array dimensions

    ! unknown grid
    if (size(global_edge_mask,1)/=local_ewn .or. size(global_edge_mask,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(global_edge_mask,1), ",", size(global_edge_mask,2), &
            ") and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! Identify cells at the edge of the global domain

    global_edge_mask = 0

    if (this_rank >= east) then  ! at east edge of global domain
       global_edge_mask(local_ewn-uhalo,:) = 1
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       global_edge_mask(lhalo+1,:) = 1
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       global_edge_mask(:,local_nsn-uhalo) = 1
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       global_edge_mask(:,lhalo+1) = 1
    endif

    call parallel_halo(global_edge_mask, parallel)

    end associate

  end subroutine parallel_global_edge_mask

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
    ! write(*,*) "Processor reference ", this_rank, local_coord, " globalID = ", global_ID

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
    ! write(*,*) "Processor reference ", this_rank, local_coord, " globalID = ", global_ID

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

  function parallel_global_sum_integer_2d(a, parallel)

    ! Calculates the global sum of a 2D integer field

    integer,dimension(:,:),intent(in) :: a
    type(parallel_type) :: parallel

    integer :: i, j
    integer :: local_sum
    integer :: parallel_global_sum_integer_2d

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn)

    local_sum = 0
    do j = nhalo+1, local_nsn-nhalo
       do i = nhalo+1, local_ewn-nhalo
          local_sum = local_sum + a(i,j)
       enddo
    enddo
    parallel_global_sum_integer_2d = parallel_reduce_sum(local_sum)

    end associate

  end function parallel_global_sum_integer_2d

!=======================================================================

  function parallel_global_sum_real4_2d(a, parallel)

    ! Calculates the global sum of a 2D single-precision field

    real(sp),dimension(:,:),intent(in) :: a
    type(parallel_type) :: parallel

    integer :: i, j
    real(sp) :: local_sum
    real(sp) :: parallel_global_sum_real4_2d

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn)

    local_sum = 0.0
    do j = nhalo+1, local_nsn-nhalo
       do i = nhalo+1, local_ewn-nhalo
          local_sum = local_sum + a(i,j)
       enddo
    enddo
    parallel_global_sum_real4_2d = parallel_reduce_sum(local_sum)

    end associate

  end function parallel_global_sum_real4_2d

!=======================================================================

  function parallel_global_sum_real8_2d(a, parallel)

    ! Calculates the global sum of a 2D double-precision field

    real(dp),dimension(:,:),intent(in) :: a
    type(parallel_type) :: parallel

    integer :: i, j
    real(dp) :: local_sum
    real(dp) :: parallel_global_sum_real8_2d

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn)

    local_sum = 0.0d0
    do j = nhalo+1, local_nsn-nhalo
       do i = nhalo+1, local_ewn-nhalo
          local_sum = local_sum + a(i,j)
       enddo
    enddo
    parallel_global_sum_real8_2d = parallel_reduce_sum(local_sum)

    end associate

  end function parallel_global_sum_real8_2d

!=======================================================================

  subroutine parallel_localindex(iglobal, jglobal, ilocal, jlocal, rlocal, parallel)

    ! Calculates the local i,j indices and rank from the global i,j indices
    integer,intent(in) :: iglobal, jglobal 
    integer,intent(out)  :: ilocal, jlocal, rlocal
    type(parallel_type) :: parallel

    integer :: flag, flag_out

    associate(  &
         own_ewn           => parallel%own_ewn,            &
         own_nsn           => parallel%own_nsn,            &
         global_row_offset => parallel%global_row_offset,  &
         global_col_offset => parallel%global_col_offset   &
         )

    flag = 0   ! This flag will be flipped on exactly one processor if the global point is valid
    ilocal = iglobal + lhalo - global_col_offset
    jlocal = jglobal + lhalo - global_row_offset

    ! Check whether these are valid values of ilocal and jlocal
    ! If so, then flip the flag and broadcast these values
    if ( (ilocal > lhalo .and. ilocal <= lhalo + own_ewn)  &
                         .and.                             &
         (jlocal > lhalo .and. jlocal <= lhalo + own_nsn) ) then
       flag = 1
    endif

    call parallel_reduce_maxloc(flag, flag_out, rlocal)

    if (flag_out==1) then
       call broadcast(ilocal, rlocal)
       call broadcast(jlocal, rlocal)
    else ! global indices are invalid
       if (main_task) then
          write(*,*) 'Invalid global indices: iglobal, jglobal =', iglobal, jglobal
          call parallel_stop(__FILE__,__LINE__)
       endif
    endif

    end associate

  end subroutine parallel_localindex

!=======================================================================

  ! subroutines belonging to the parallel_halo interface

  subroutine parallel_halo_integer_2d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel
    
    integer :: erequest, nrequest, srequest, wrequest, ierror
    integer,dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    integer,dimension(parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
      a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:) = 0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:) = 0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:) = 0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo) = 0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:) = 0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:) = 0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:) = 0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1) = 0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1) = 0
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1) = 0
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:) = 0
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:) = 0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_integer_2d


  subroutine parallel_halo_logical_2d(a, parallel)

    use mpi_mod
    implicit none
    logical,dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,nrequest,srequest,wrequest
    logical,dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    logical,dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    logical,dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    logical,dimension(parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_logical,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_logical,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_logical,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_logical,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_logical,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_logical,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_logical,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_logical,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                        ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:) = .false.
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:) = .false.
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:) = .false.
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo) = .false.
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:) = .false.
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:) = .false.
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:) = .false.
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1) = .false.
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1) = .false.
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1) = .false.
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:) = .false.
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:) = .false.

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_logical_2d


  subroutine parallel_halo_real4_2d(a, parallel)

    use mpi_mod
    implicit none
    real(sp),dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,nrequest,srequest,wrequest
    real(sp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(sp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(sp),dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    real(sp),dimension(parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real4,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real4,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real4,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real4,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real4,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real4,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real4,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real4,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                        ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:) = 0.
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:) = 0.
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:) = 0.
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo) = 0.
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:) = 0.
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:) = 0.
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:) = 0.
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1) = 0.
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1) = 0.
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1) = 0.
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:) = 0.
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:) = 0.

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_real4_2d


  subroutine parallel_halo_real8_2d(a, parallel, periodic_offset_ew, periodic_offset_ns)

    !WHL - added optional arguments for periodic offsets, to support ismip-hom test cases

    use mpi_mod
    implicit none
    real(dp),dimension(:,:) :: a
    type(parallel_type) :: parallel
    real(dp), intent(in), optional :: &
       periodic_offset_ew,  &! offset halo values by this amount
                             ! if positive, the offset is positive for W halo, negative for E halo
       periodic_offset_ns    ! offset halo values by this amount
                             ! if positive, the offset is positive for S halo, negative for N halo
    
    integer :: erequest,ierror,nrequest,srequest,wrequest
    real(dp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
       write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    if (present(periodic_offset_ew)) then
       if (periodic_offset_ew /= 0.d0) then
          if (this_rank <= west) then   ! this proc lies at the west edge of the global domain
!             print*, 'Offset at west edge: this_rank, west =', this_rank, west
             a(:lhalo,1+lhalo:local_nsn-uhalo) =   &
                a(:lhalo,1+lhalo:local_nsn-uhalo) + periodic_offset_ew
          endif
          if (this_rank >= east) then   ! this proc lies at the east edge of the global domain
!             print*, 'Offset at east edge: this_rank, east =', this_rank, east
             a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) =    &
                a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) - periodic_offset_ew
          endif
       endif
    endif

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (present(periodic_offset_ns)) then
       if (periodic_offset_ns /= 0.d0) then
          if (this_rank <= south) then  ! this proc lies at the south edge of the global domain
!             print*, 'Offset at south edge: this_rank, south =', this_rank, south
             a(:,:lhalo) = a(:,:lhalo) + periodic_offset_ns
          endif
          if (this_rank >= north) then  ! this proc lies at the north edge of the global domain
!             print*, 'Offset at north edge: this_rank, north =', this_rank, north
             a(:,local_nsn-uhalo+1:) = a(:,local_nsn-uhalo+1:) - periodic_offset_ns
          endif
       endif
    endif

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1) = 0.d0
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1) = 0.d0
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:) = 0.d0
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_real8_2d


  subroutine parallel_halo_real8_3d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(dp),dimension(size(a,1), lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(size(a,1), uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(size(a,1), parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1), parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,2)==local_ewn-1.and.size(a,3)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn.or.size(a,3)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:,:)

    nsend(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,local_nsn-uhalo+1:) = nrecv(:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,local_ewn-uhalo+1:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:lhalo,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,local_nsn-uhalo+1:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:lhalo) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,local_ewn-uhalo:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:lhalo+1,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,local_nsn-uhalo:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:lhalo+1) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:,:lhalo+1,:lhalo+1) = 0.d0
       if (southeast_corner) a(:,local_ewn-lhalo:,:lhalo+1) = 0.d0
       if (northeast_corner) a(:,local_ewn-lhalo:,local_nsn-lhalo:) = 0.d0
       if (northwest_corner) a(:,:lhalo+1,local_nsn-lhalo:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_real8_3d


  subroutine parallel_halo_real8_4d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(dp),dimension(size(a,1), size(a,2), lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(size(a,1), size(a,2), uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(size(a,1), size(a,2), parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1), size(a,2), parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,3)==local_ewn-1.and.size(a,4)==local_nsn-1) return

    ! unknown grid
    if (size(a,3)/=local_ewn.or.size(a,4)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ",", size(a,4), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:,:) = &
         a(:,:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:,:,:)

    nsend(:,:,:,:) = a(:,:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:,:) = a(:,:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,:,:lhalo) = srecv(:,:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,:,local_nsn-uhalo+1:) = nrecv(:,:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,:,local_ewn-uhalo+1:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:,:lhalo,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,:,local_nsn-uhalo+1:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:,:lhalo) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,:,local_ewn-uhalo:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:,:lhalo+1,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,:,local_nsn-uhalo:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:,:lhalo+1) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:,:,:lhalo+1,:lhalo+1) = 0.d0
       if (southeast_corner) a(:,:,local_ewn-lhalo:,:lhalo+1) = 0.d0
       if (northeast_corner) a(:,:,local_ewn-lhalo:,local_nsn-lhalo:) = 0.d0
       if (northwest_corner) a(:,:,:lhalo+1,local_nsn-lhalo:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_real8_4d

!=======================================================================

  ! subroutines belonging to the parallel_halo_extrapolate interface

  subroutine parallel_halo_extrapolate_integer_2d(a, parallel)

    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: i, j

    ! begin
    associate(  &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    ! Extrapolate the field into halo cells along the global boundary.

    ! First update the halos so that we are sure the interior halos are correct
    call parallel_halo(a, parallel)

! Useful for debugging small domains (the YYYY is just a tag for grepping the output,
!  particularly if you prepend the processor number, e.g. "0YYYY")
!  do j = 1, size(a,2)
!     write(6, "(i3, 'YYYY BEFORE row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    if (this_rank >= east) then  ! at east edge of global domain
       ! extrapolate eastward
       do i = size(a,1)-uhalo+1, size(a,1)
          a(i, :) = a(size(a,1)-uhalo, :)
       enddo
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! extrapolate westward
       do i = 1, lhalo
          a(i, :) = a(lhalo+1, :)
       enddo
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       ! extrapolate northward
       do j = size(a,2)-uhalo+1, size(a,2)
          a(:, j) = a(:, size(a,2)-uhalo)
       enddo
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! extrapolate southward
       do j = 1, lhalo
          a(:, j) = a(:, lhalo+1)
       enddo
    endif

! Useful for debugging small domains
!  do j = 1, size(a,2)
!     write(6, "(i3, 'YYYY AFTER  row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    end associate

  end subroutine parallel_halo_extrapolate_integer_2d


  subroutine parallel_halo_extrapolate_real8_2d(a, parallel)

    implicit none
    real(dp),dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: i, j

    ! begin
    associate(  &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    ! Extrapolate the field into halo cells along the global boundary.

    ! First update the halos so that we are sure the interior halos are correct
    call parallel_halo(a, parallel)

! Useful for debugging small domains (the YYYY is just a tag for grepping the output,
!  particularly if you prepend the processor number, e.g. "0YYYY")
!  do j = 1, size(a,2)
!     write(6, "(i3, 'YYYY BEFORE row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    if (this_rank >= east) then  ! at east edge of global domain
       ! extrapolate eastward
       do i = size(a,1)-uhalo+1, size(a,1)
          a(i, :) = a(size(a,1)-uhalo, :)
       enddo
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! extrapolate westward
       do i = 1, lhalo
          a(i, :) = a(lhalo+1, :)
       enddo
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       ! extrapolate northward
       do j = size(a,2)-uhalo+1, size(a,2)
          a(:, j) = a(:, size(a,2)-uhalo)
       enddo
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! extrapolate southward
       do j = 1, lhalo
          a(:, j) = a(:, lhalo+1)
       enddo
    endif

! Useful for debugging small domains
!  do j = 1, size(a,2)
!     write(6, "(i3, 'YYYY AFTER  row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    end associate

  end subroutine parallel_halo_extrapolate_real8_2d

!=======================================================================

  ! subroutines belonging to the parallel_halo_tracers interface

  subroutine parallel_halo_tracers_real8_3d(a, parallel)

    ! Custom halo routine for tracer arrays with dimension(nx,ny,ntracers)
    ! Will work for any 3D array with (nx,ny) in the first two slots

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: ierror, one, erequest, nrequest, srequest, wrequest
    real(dp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo, size(a,3)) :: esend,wrecv
    real(dp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo, size(a,3)) :: erecv,wsend
    real(dp),dimension(parallel%local_ewn, lhalo, size(a,3)) :: nsend,srecv
    real(dp),dimension(parallel%local_ewn, uhalo, size(a,3)) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn .or. size(a,2)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo,:)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo,:)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo,:) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo,:) = erecv(:,:,:)

    nsend(:,:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo,:)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,:)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo,:) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:,:) = nrecv(:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:,:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo,:) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:,:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1,:) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1,:) = 0.d0
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1,:) = 0.d0
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:,:) = 0.d0
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:,:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_tracers_real8_3d

  subroutine parallel_halo_tracers_real8_4d(a, parallel)

    ! Custom halo routine for tracer arrays with dimension(nx,ny,ntracers,nz)
    ! Will work for any 4D array with (nx,ny) in the first two slots

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: ierror, one, erequest, nrequest, srequest, wrequest
    real(dp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo, size(a,3), size(a,4)) :: esend,wrecv
    real(dp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo, size(a,3), size(a,4)) :: erecv,wsend
    real(dp),dimension(parallel%local_ewn, lhalo, size(a,3), size(a,4)) :: nsend,srecv
    real(dp),dimension(parallel%local_ewn, uhalo, size(a,3), size(a,4)) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn .or. size(a,2)/=local_nsn) then
         write(*,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ",", size(a,4), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo,:,:)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo,:,:)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo,:,:) = wrecv(:,:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo,:,:) = erecv(:,:,:,:)

    nsend(:,:,:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo,:,:)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,:,:)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo,:,:) = srecv(:,:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:,:,:) = nrecv(:,:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:,:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:,:,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:,:,:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo,:,:) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:,:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:,:,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:,:,:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1,:,:) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1,:,:) = 0.d0
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1,:,:) = 0.d0
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:,:,:) = 0.d0
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:,:,:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_tracers_real8_4d

!=======================================================================

  ! subroutines belonging to the parallel_halo_verify interface

  function parallel_halo_verify_integer_2d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: ierror, erequest, nrequest, srequest, wrequest
    integer,dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    integer,dimension(parallel%local_ewn, uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_integer_2d

    ! begin

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    ! ANY True if any value is true (LOGICAL)
    notverify_flag = ANY(a(:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:))

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:lhalo) /= srecv(:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,local_nsn-uhalo+1:) /= nrecv(:,:))

    ! if notverify_flag is TRUE, then there was some difference detected
    if (notverify_flag) then
         write(*,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_integer_2d = .NOT. notverify_flag  ! return if verified (True) or not verified (False)

    end associate

  end function parallel_halo_verify_integer_2d


  function parallel_halo_verify_real8_2d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:) :: a
    type(parallel_type) :: parallel
    
    integer :: ierror, erequest, nrequest, srequest, wrequest
    real(dp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(parallel%local_ewn, uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_real8_2d

    ! begin
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    notverify_flag = ANY(a(:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:))

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:lhalo) /= srecv(:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,local_nsn-uhalo+1:) /= nrecv(:,:))

    if (notverify_flag) then
         write(*,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_real8_2d = .NOT. notverify_flag

    end associate

  end function parallel_halo_verify_real8_2d


  function parallel_halo_verify_real8_3d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:) :: a
    type(parallel_type) :: parallel
    
    integer :: ierror, one, erequest, nrequest, srequest, wrequest
    real(dp),dimension(size(a,1), lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(size(a,1), uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(size(a,1), parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1), parallel%local_ewn, uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_real8_3d

    ! begin
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,2)==local_ewn-1.and.size(a,3)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn.or.size(a,3)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    notverify_flag = ANY(a(:,:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:,:))

    nsend(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:,:lhalo) /= srecv(:,:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:,local_nsn-uhalo+1:) /= nrecv(:,:,:))

    if (notverify_flag) then
         write(*,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_real8_3d = .NOT. notverify_flag

    end associate

  end function parallel_halo_verify_real8_3d

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

  function parallel_inq_attname(ncid, varid, attnum, name)

    implicit none
    integer :: attnum, ncid, parallel_inq_attname, varid
    character(len=*) :: name

    ! begin
    if (main_task) parallel_inq_attname = &
         nf90_inq_attname(ncid,varid,attnum,name)
    call broadcast(parallel_inq_attname)
    call broadcast(name)

  end function parallel_inq_attname

!=======================================================================

  function parallel_inq_dimid(ncid, name, dimid)

    implicit none
    integer :: dimid, ncid, parallel_inq_dimid
    character(len=*) :: name

    ! begin
    if (main_task) parallel_inq_dimid = nf90_inq_dimid(ncid,name,dimid)
    call broadcast(parallel_inq_dimid)
    call broadcast(dimid)

  end function parallel_inq_dimid

!=======================================================================

  function parallel_inq_varid(ncid, name, varid)

    implicit none
    integer :: ncid,parallel_inq_varid,varid
    character(len=*) :: name

    ! begin
    if (main_task) parallel_inq_varid = nf90_inq_varid(ncid,name,varid)
    call broadcast(parallel_inq_varid)
    call broadcast(varid)

  end function parallel_inq_varid

!=======================================================================

  function parallel_inquire(ncid, nvariables)

    implicit none
    integer :: ncid,parallel_inquire,nvariables

    ! begin
    if (main_task) parallel_inquire = nf90_inquire(ncid,nvariables=nvariables)
    call broadcast(parallel_inquire)
    call broadcast(nvariables)

  end function parallel_inquire

!=======================================================================

  function parallel_inquire_dimension(ncid, dimid, name, len)

    implicit none
    integer :: dimid, ncid, parallel_inquire_dimension
    integer,optional :: len
    character(len=*),optional :: name
    
    integer :: l
    
    ! begin
    if (present(name)) then
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,name,len=l)
       call broadcast(name)
    else
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,len=l)
    end if
    call broadcast(parallel_inquire_dimension)
    if (present(len)) then
       call broadcast(l)
       len = l
    end if

  end function parallel_inquire_dimension

!=======================================================================

  function parallel_inquire_variable(ncid, varid, &
                                     name, ndims, dimids, natts)

    implicit none
    integer :: ncid, parallel_inquire_variable, varid
    integer,optional :: ndims, natts
    character(len=*),optional :: name
    integer,dimension(:),optional :: dimids

    integer :: nd,na

    ! begin
    if (present(name)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,name=name)
       call broadcast(parallel_inquire_variable)
       call broadcast(name)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (present(dimids)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,dimids=dimids)
       call broadcast(parallel_inquire_variable)
       call broadcast(dimids)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (main_task) parallel_inquire_variable = &
         nf90_inquire_variable(ncid,varid,ndims=nd,natts=na)
    call broadcast(parallel_inquire_variable)
    if (present(ndims)) then
       call broadcast(nd)
       ndims = nd
    end if
    if (present(natts)) then
       call broadcast(na)
       natts = na
    end if

  end function parallel_inquire_variable

!=======================================================================

  function parallel_open(path, mode, ncid)

    implicit none
    integer :: mode, ncid, parallel_open
    character(len=*) :: path

    ! begin
    if (main_task) parallel_open = nf90_open(path,mode,ncid)
    call broadcast(parallel_open)

  end function parallel_open

!=======================================================================

  subroutine parallel_print_all(name, values, parallel)

    implicit none
    character(*) :: name
    real(dp),dimension(:,:,:) :: values
    type(parallel_type) :: parallel

    integer,parameter :: u = 33
    integer :: i,j,t

    ! begin
    associate(  &
         ewlb      => parallel%ewlb,        &
         nslb      => parallel%nslb         &
         )

    if (main_task) then
       open(unit=u,file=name,form="formatted",status="replace")
       close(u)
    end if
    do t = 0,tasks-1
       call parallel_barrier
       if (t==this_rank) then
          open(unit=u,file=name,form="formatted",position="append")
          do j = 1,size(values,3)
             do i = 1,size(values,2)
                write(u,'(2i5,100g15.5e3)') nslb+j-1,ewlb+i-1,values(:,i,j)
             end do
             write(u,'()')
          end do
          write(u,'(//)')
          close(u)
       end if
    end do

    end associate

  end subroutine parallel_print_all

!=======================================================================

  ! subroutines belonging to the parallel_print interface

  subroutine parallel_print_integer_2d(name, values)

    implicit none
    character(*) :: name
    integer,dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j

    ! begin

    if (main_task) then
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
       close(u)
    end if

    call parallel_barrier  ! Only the main_task writes the variable.  Rest wait here.

    ! automatic deallocation

  end subroutine parallel_print_integer_2d


  subroutine parallel_print_real8_2d(name, values)

    implicit none
    character(*) :: name
    real(dp),dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j

    ! begin
    if (main_task) then
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
       close(u)
    end if

    call parallel_barrier  ! Only the main_task writes the variable.  Rest wait here.

  end subroutine parallel_print_real8_2d


  subroutine parallel_print_real8_3d(name, values)

    implicit none
    character(*) :: name
    real(dp),dimension(:,:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,j

    ! begin
    if (main_task) then
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       do j = lbound(values,3),ubound(values,3)
          do i = lbound(values,2),ubound(values,2)
             write(u,'(2i6,100g15.5e3)') j,i,values(:,i,j)
          end do
          write(u,'()')
       end do
       close(u)
    end if

    call parallel_barrier  ! Only the main_task writes the variable.  Rest wait here.

  end subroutine parallel_print_real8_3d

!=======================================================================

  ! functions belonging to the parallel_put_att interface

  function parallel_put_att_character(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_put_att_character, varid
    character(len=*) :: name, values

    ! begin
    if (main_task) parallel_put_att_character = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_character)

  end function parallel_put_att_character


  function parallel_put_att_integer(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_put_att_integer, varid
    character(len=*) :: name
    integer :: values

    ! begin
    if (main_task) parallel_put_att_integer = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_integer)

  end function parallel_put_att_integer


  function parallel_put_att_real4(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_put_att_real4,varid
    character(len=*) :: name
    real(sp) :: values

    ! begin
    if (main_task) parallel_put_att_real4 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4)

  end function parallel_put_att_real4


  function parallel_put_att_real4_1d(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_put_att_real4_1d,varid
    character(len=*) :: name
    real(sp),dimension(:) :: values

    ! begin
    if (main_task) parallel_put_att_real4_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4_1d)

  end function parallel_put_att_real4_1d


  function parallel_put_att_real8(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_put_att_real8, varid
    character(len=*) :: name
    real(dp) :: values

    ! begin
    if (main_task) parallel_put_att_real8 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8)

  end function parallel_put_att_real8


  function parallel_put_att_real8_1d(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_put_att_real8_1d,varid
    character(len=*) :: name
    real(dp),dimension(:) :: values

    ! begin
    if (main_task) parallel_put_att_real8_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8_1d)

  end function parallel_put_att_real8_1d

!=======================================================================

  ! functions belonging to the parallel_put_var interface

  function parallel_put_var_integer(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_integer,varid
    integer :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_integer = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_integer = nf90_put_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_put_var_integer)

  end function parallel_put_var_integer


  function parallel_put_var_integer_1d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_integer_1d,varid
    integer,dimension(:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_integer_1d = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_integer_1d = nf90_put_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_put_var_integer_1d)

  end function parallel_put_var_integer_1d


  function parallel_put_var_real4(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_real4,varid
    real(sp) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real4 = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real4 = nf90_put_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_put_var_real4)

  end function parallel_put_var_real4


  function parallel_put_var_real8(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_real8,varid
    real(dp) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real8 = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real8 = nf90_put_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_put_var_real8)

  end function parallel_put_var_real8


  function parallel_put_var_real8_1d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_real8_1d,varid
    real(dp),dimension(:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values)
       end if
    end if
    call broadcast(parallel_put_var_real8_1d)

  end function parallel_put_var_real8_1d

!=======================================================================

  function parallel_redef(ncid)

    implicit none
    integer :: ncid,parallel_redef

    ! begin
    if (main_task) parallel_redef = nf90_redef(ncid)
    call broadcast(parallel_redef)

  end function parallel_redef

!=======================================================================

  ! functions for parallel reduction of logical variables
  ! * parallel_reduce_log_or returns 'true' iff x = 'true' on at least one processor
  ! * parallel_reduce_log_and returns 'true' iff x = 'true' on all processors

  function parallel_reduce_log_or(x)

    use mpi_mod
    implicit none
    logical :: x

    integer :: ierror
    logical :: recvbuf,sendbuf, parallel_reduce_log_or

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_logical,mpi_lor,comm,ierror)
    parallel_reduce_log_or = recvbuf

  end function parallel_reduce_log_or

  function parallel_reduce_log_and(x)

    use mpi_mod
    implicit none
    logical :: x

    integer :: ierror
    logical :: recvbuf,sendbuf, parallel_reduce_log_and

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_logical,mpi_land,comm,ierror)
    parallel_reduce_log_and = recvbuf

  end function parallel_reduce_log_and

!=======================================================================

  ! functions belonging to the parallel_reduce_sum interface

  function parallel_reduce_sum_integer(x)

    use mpi_mod
    implicit none
    integer :: x

    integer :: ierror
    integer :: recvbuf,sendbuf, parallel_reduce_sum_integer

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_integer,mpi_sum,comm,ierror)
    parallel_reduce_sum_integer = recvbuf

  end function parallel_reduce_sum_integer


  function parallel_reduce_sum_real4(x)

    use mpi_mod
    implicit none
    real(sp) :: x

    integer :: ierror
    real(sp) :: recvbuf,sendbuf, parallel_reduce_sum_real4

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real4,mpi_sum,comm,ierror)
    parallel_reduce_sum_real4 = recvbuf

  end function parallel_reduce_sum_real4


  function parallel_reduce_sum_real8(x)

    use mpi_mod
    implicit none
    real(dp) :: x

    integer :: ierror
    real(dp) :: recvbuf,sendbuf, parallel_reduce_sum_real8

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real8,mpi_sum,comm,ierror)
    parallel_reduce_sum_real8 = recvbuf

  end function parallel_reduce_sum_real8


  function parallel_reduce_sum_integer_nvar(x)

    use mpi_mod
    implicit none
    integer :: x(:)

    integer :: ierror, nvar
    integer, dimension(size(x)) :: recvbuf,sendbuf, parallel_reduce_sum_integer_nvar

    ! begin
    nvar = size(x)
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,nvar,mpi_integer,mpi_sum,comm,ierror)
    parallel_reduce_sum_integer_nvar = recvbuf

  end function parallel_reduce_sum_integer_nvar


  function parallel_reduce_sum_real8_nvar(x)

    use mpi_mod
    implicit none
    real(dp) :: x(:)

    integer :: ierror, nvar
    real(dp), dimension(size(x)) :: recvbuf,sendbuf, parallel_reduce_sum_real8_nvar

    ! begin
    nvar = size(x)
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,nvar,mpi_real8,mpi_sum,comm,ierror)
    parallel_reduce_sum_real8_nvar = recvbuf

  end function parallel_reduce_sum_real8_nvar

!=======================================================================

  ! functions belonging to the parallel_reduce_max interface

  function parallel_reduce_max_integer(x)

    use mpi_mod
    implicit none
    integer :: x

    integer :: ierror
    integer :: recvbuf,sendbuf, parallel_reduce_max_integer

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_integer,mpi_max,comm,ierror)
    parallel_reduce_max_integer = recvbuf

  end function parallel_reduce_max_integer


  function parallel_reduce_max_real4(x)

    use mpi_mod
    implicit none
    real(sp) :: x

    integer :: ierror
    real(sp) :: recvbuf,sendbuf, parallel_reduce_max_real4

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real4,mpi_max,comm,ierror)
    parallel_reduce_max_real4 = recvbuf

  end function parallel_reduce_max_real4


  function parallel_reduce_max_real8(x)

    use mpi_mod
    implicit none
    real(dp) :: x

    integer :: ierror
    real(dp) :: recvbuf,sendbuf, parallel_reduce_max_real8

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real8,mpi_max,comm,ierror)
    parallel_reduce_max_real8 = recvbuf

  end function parallel_reduce_max_real8

  function parallel_reduce_max_real8_1d(x)

    use mpi_mod
    implicit none
    real(dp), dimension(:) :: x

    integer :: ierror
    real(dp), dimension(size(x)) :: recvbuf,sendbuf, parallel_reduce_max_real8_1d

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,size(x),mpi_real8,mpi_max,comm,ierror)
    parallel_reduce_max_real8_1d = recvbuf

  end function parallel_reduce_max_real8_1d

!=======================================================================

  ! functions belonging to the parallel_reduce_maxloc interface

  subroutine parallel_reduce_maxloc_integer(xin, xout, xprocout)

    use mpi_mod
    implicit none
    integer, intent(in) :: xin         ! variable to reduce
    integer, intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    integer :: ierror
    integer, dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2INTEGER,mpi_maxloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1)

  end subroutine parallel_reduce_maxloc_integer


  subroutine parallel_reduce_maxloc_real4(xin, xout, xprocout)

    use mpi_mod
    implicit none
    real(sp), intent(in) :: xin         ! variable to reduce
    real(sp), intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout    ! processor on which reduced value occurs

    integer :: ierror
    real(sp), dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x (coerced to a real)
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2REAL,mpi_maxloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1) ! coerced back to integer

  end subroutine parallel_reduce_maxloc_real4


  subroutine parallel_reduce_maxloc_real8(xin, xout, xprocout)

    use mpi_mod
    implicit none
    real(dp), intent(in) :: xin         ! variable to reduce
    real(dp), intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout    ! processor on which reduced value occurs

    integer :: ierror
    real(dp), dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x (coerced to a real)
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION,mpi_maxloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1) ! coerced back to integer

  end subroutine parallel_reduce_maxloc_real8


!=======================================================================

  ! subroutines belonging to the parallel_reduce_min interface

  function parallel_reduce_min_integer(x)

    use mpi_mod
    implicit none
    integer :: x

    integer :: ierror
    integer :: recvbuf,sendbuf, parallel_reduce_min_integer

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_integer,mpi_min,comm,ierror)
    parallel_reduce_min_integer = recvbuf

  end function parallel_reduce_min_integer


  function parallel_reduce_min_real4(x)

    use mpi_mod
    implicit none
    real(sp) :: x

    integer :: ierror
    real(sp) :: recvbuf,sendbuf, parallel_reduce_min_real4

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real4,mpi_min,comm,ierror)
    parallel_reduce_min_real4 = recvbuf

  end function parallel_reduce_min_real4


  function parallel_reduce_min_real8(x)

    use mpi_mod
    implicit none
    real(dp) :: x

    integer :: ierror
    real(dp) :: recvbuf,sendbuf, parallel_reduce_min_real8

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real8,mpi_min,comm,ierror)
    parallel_reduce_min_real8 = recvbuf

  end function parallel_reduce_min_real8


  function parallel_reduce_min_real8_1d(x)

    use mpi_mod
    implicit none
    real(dp), dimension(:) :: x

    integer :: ierror
    real(dp), dimension(size(x)) :: recvbuf,sendbuf, parallel_reduce_min_real8_1d

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,size(x),mpi_real8,mpi_min,comm,ierror)
    parallel_reduce_min_real8_1d = recvbuf

  end function parallel_reduce_min_real8_1d

!=======================================================================

  ! subroutines belonging to the parallel_reduce_minloc interface

  subroutine parallel_reduce_minloc_integer(xin, xout, xprocout)

    use mpi_mod
    implicit none
    integer, intent(in) :: xin         ! variable to reduce
    integer, intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    integer :: ierror
    integer, dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2INTEGER,mpi_minloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1)

  end subroutine parallel_reduce_minloc_integer


  subroutine parallel_reduce_minloc_real4(xin, xout, xprocout)

    use mpi_mod
    implicit none
    real(sp), intent(in) :: xin        ! variable to reduce
    real(sp), intent(out) :: xout      ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    integer :: ierror
    real(sp), dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x (coerced to a real)
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2REAL,mpi_minloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1) ! coerced back to integer

  end subroutine parallel_reduce_minloc_real4


  subroutine parallel_reduce_minloc_real8(xin, xout, xprocout)

    use mpi_mod
    implicit none
    real(dp), intent(in) :: xin        ! variable to reduce
    real(dp), intent(out) :: xout      ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    integer :: ierror
    real(dp), dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x (coerced to a real)
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION,mpi_minloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1) ! coerced back to integer

  end subroutine parallel_reduce_minloc_real8

!=======================================================================

  subroutine parallel_show_minmax(label,values)

    use mpi_mod
    implicit none
    character(*) :: label
    real(dp),dimension(:,:,:) :: values
    
    integer :: ierror
    real(dp) :: allmin,allmax,mymin,mymax

    ! begin
    mymin = minval(values(:,1+lhalo:size(values,2)-uhalo,&
         1+lhalo:size(values,3)-uhalo))
    mymax = maxval(values(:,1+lhalo:size(values,2)-uhalo,&
         1+lhalo:size(values,3)-uhalo))
    call mpi_reduce(mymin,allmin,1,mpi_real8,mpi_min,main_rank,comm,ierror)
    call mpi_reduce(mymax,allmax,1,mpi_real8,mpi_max,main_rank,comm,ierror)
    if (main_task) print *,label,allmin,allmax

  end subroutine parallel_show_minmax

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

  function parallel_sync(ncid)

    implicit none
    integer :: ncid,parallel_sync

    ! begin
    if (main_task) parallel_sync = nf90_sync(ncid)
    call broadcast(parallel_sync)

  end function parallel_sync

!=======================================================================

  subroutine staggered_no_penetration_mask(umask, vmask, parallel)

    implicit none
    integer,dimension(:,:) :: umask, vmask  ! mask set to 1 wherever the outflow velocity should be zero
    type(parallel_type) :: parallel

    associate(  &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         east      => parallel%east,       &
         west      => parallel%west,       &
         north     => parallel%north,      &
         south     => parallel%south       &
         )

    ! initialize the no-penetration masks to 0
    umask(:,:) = 0
    vmask(:,:) = 0

    if (this_rank >= east) then  ! at east edge of global domain
       ! set u velocity mask = 1 at the east global boundary and vertices eastward
       umask(local_ewn-uhalo:,:) = 1
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! set u velocity mask = 1 at the west global boundary and vertices westward
       umask(:lhalo,:) = 1
    endif
    
    if (this_rank >= north) then  ! at north edge of global domain
       ! set v velocity mask = 1 at the north global boundary and vertices northward
       vmask(:,local_nsn-uhalo:) = 1
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! set v velocity mask = 1 at the south global boundary and vertices southward
       vmask(:,:lhalo) = 1
    endif

    call staggered_parallel_halo(umask, parallel)
    call staggered_parallel_halo(vmask, parallel)

    end associate

  end subroutine staggered_no_penetration_mask

!=======================================================================

  ! subroutines belonging to the staggered_parallel_halo interface

  !-----------------------------------------------------------------
  ! Comments on the staggered_parallel_halo subroutines:
  !
  ! The following subroutines implement a staggered grid halo update for integer and real(dp) arrays
  !  of various sizes.
  ! As the grid is staggered, the array 'a' is one smaller in both dimensions than an unstaggered array.
  !
  ! The grid is laid out from the SW, and the lower left corner is assigned to this_rank = 0.
  ! Its eastern neighbor is task_id = 1, proceeding rowwise and starting from the western edge.
  ! The South-most processes own one additional row of stagggered variables on the southern edge
  ! and have one less 'southern' halo row than other processes. Likewise, the West-most processes own one 
  ! additional column of staggered variables on the western edge and have one less 'western' halo column. 
  !
  ! The default BCs are periodic, in which case the southernmost and westernmost rows of the 
  ! global domain are filled with data from the northernmost and easternmost rows.
  !
  ! WHL: Comments on the implementation of outflow BC:
  !
  ! With outflow BCs, the southernmost and westernmost rows contain correct data (e.g., velocity
  ! at global boundaries) that should be passed to adjacent processors and should not be overwritten.
  ! This is implemented by defining variables called staggered_ilo, staggered_ihi, staggered_jlo
  ! and staggered_jhi, which denote the limits of locally owned staggered data on each processor.
  !
  ! For processors that include the western global boundary, staggered_ilo = staggered_lhalo; 
  !  for other processors (or for all processors when periodic_bc = T), staggered_ilo = staggered_lhalo+1. 
  ! Similarly, for processors that include the southern global boundary, staggered_jlo = staggered_lhalo; 
  !  for other processors (or for all processors when periodic_bc = T), staggered_jlo = staggered_lhalo+1.
  ! For all processors, staggered_ihi = (nx-1)-staggered_uhalo and staggered_jhi = (ny-1)-staggered_uhalo.
  ! Using these limits, we can correctly pass all and only the locally owned staggered data
  !  from each processor to its neighbors.
  ! Data are passed between the northern and southern global rows, and between the eastern and western rows,
  !  only when periodic_bc = T.
  !
  ! In addition to correctly passing staggered data along global boundaries, outflow BC imply
  !  that all staggered data are zero at points beyond the global boundary. The global boundary
  !  is delimited by:
  !    i = staggered_lhalo on the west boundary
  !    i = nx - 1 - staggered_uhalo on the east boundary
  !    j = staggered_lhalo on the south boundary
  !    j = ny - 1 - staggered_uhalo on the north boundary
  ! 
  ! TODO: Think about whether code simplifications are possible below using staggered_ilo/ihi/jlo/jhi.
  !
  !-----------------------------------------------------------------

  subroutine staggered_parallel_halo_integer_2d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a 2D integer field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    integer,dimension(staggered_lhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    integer,dimension(staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    integer,dimension(size(a,1),staggered_lhalo) :: nsend,srecv
    integer,dimension(size(a,1),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &       
         north     => parallel%north,          &    
         south     => parallel%south,          &      
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      wsend(:, 1:staggered_jhi-staggered_jlo+1) = &
           a(1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      esend(:, 1:staggered_jhi-staggered_jlo+1) = &
           a(size(a,1)-staggered_uhalo-staggered_lhalo+1:size(a,1)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(size(a,1)-staggered_uhalo+1:size(a,1), staggered_jlo:staggered_jhi) = &
          erecv(:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:) = &
        a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:) = &
        a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,size(a,2)-staggered_uhalo+1:size(a,2)) = nrecv(:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,1:staggered_lhalo) = srecv(:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary
    ! For no_ice BC, zero the field along and beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(1:staggered_lhalo-1, :) = 0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(size(a,1)-staggered_uhalo+1:size(a,1), :) = 0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, 1:staggered_lhalo-1) = 0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, size(a,2)-staggered_uhalo+1:size(a,2)) = 0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(1:staggered_lhalo, :) = 0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(size(a,1)-staggered_uhalo:size(a,1), :) = 0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, 1:staggered_lhalo) = 0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, size(a,2)-staggered_uhalo:size(a,2)) = 0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:staggered_lhalo, :staggered_lhalo) = 0
       if (southeast_corner) a(size(a,1)-staggered_uhalo:size(a,1), :staggered_lhalo) = 0
       if (northeast_corner) a(size(a,1)-staggered_uhalo:size(a,1), &
                               size(a,2)-staggered_uhalo:size(a,2)) = 0
       if (northwest_corner) a(:staggered_lhalo, size(a,2)-staggered_uhalo:size(a,2)) = 0

    endif   ! outflow or no_ice_bc

    end associate

  end subroutine staggered_parallel_halo_integer_2d


  subroutine staggered_parallel_halo_integer_3d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a 3D integer field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    integer,dimension(size(a,1), staggered_lhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    integer,dimension(size(a,1), staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    integer,dimension(size(a,1),size(a,2),staggered_lhalo) :: nsend,srecv
    integer,dimension(size(a,1),size(a,2),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &
         north     => parallel%north,          &
         south     => parallel%south,          &
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,2)/=local_ewn-1.or.size(a,3)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      wsend(:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
       esend(:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:, size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(:, size(a,2)-staggered_uhalo+1:size(a,2), staggered_jlo:staggered_jhi) = &
          erecv(:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(:, 1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:,:) = &
        a(:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:,:) = &
        a(:,:,size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,:,size(a,3)-staggered_uhalo+1:size(a,3)) = nrecv(:,:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,:,1:staggered_lhalo) = srecv(:,:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then
          a(:, 1:staggered_lhalo-1, :) = 0
       endif

       if (this_rank >= east) then
          a(:, size(a,2)-staggered_uhalo+1:size(a,2), :) = 0
       endif

       if (this_rank <= south) then
          a(:, :, 1:staggered_lhalo-1) = 0
       endif

       if (this_rank >= north) then
          a(:, :, size(a,3)-staggered_uhalo+1:size(a,3)) = 0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(:, 1:staggered_lhalo, :) = 0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(:, size(a,2)-staggered_uhalo:size(a,2), :) = 0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, :, 1:staggered_lhalo) = 0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, :, size(a,3)-staggered_uhalo:size(a,3)) = 0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:, :staggered_lhalo, :staggered_lhalo) = 0
       if (southeast_corner) a(:, size(a,2)-staggered_uhalo:size(a,2), :staggered_lhalo) = 0
       if (northeast_corner) a(:, size(a,2)-staggered_uhalo:size(a,2), &
                               size(a,3)-staggered_uhalo:size(a,3)) = 0
       if (northwest_corner) a(:, :staggered_lhalo, size(a,3)-staggered_uhalo:size(a,3)) = 0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine staggered_parallel_halo_integer_3d


  subroutine staggered_parallel_halo_real8_2d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a real(dp) 2D field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    real(dp),dimension(staggered_lhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    real(dp),dimension(staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    real(dp),dimension(size(a,1),staggered_lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &
         north     => parallel%north,          &
         south     => parallel%south,          &
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      wsend(:, 1:staggered_jhi-staggered_jlo+1) = &
           a(1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      esend(:, 1:staggered_jhi-staggered_jlo+1) = &
           a(size(a,1)-staggered_uhalo-staggered_lhalo+1:size(a,1)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(size(a,1)-staggered_uhalo+1:size(a,1), staggered_jlo:staggered_jhi) = &
          erecv(:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:) = &
        a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:) = &
        a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,size(a,2)-staggered_uhalo+1:size(a,2)) = nrecv(:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,1:staggered_lhalo) = srecv(:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(1:staggered_lhalo-1, :) = 0.0d0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(size(a,1)-staggered_uhalo+1:size(a,1), :) = 0.0d0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, 1:staggered_lhalo-1) = 0.0d0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, size(a,2)-staggered_uhalo+1:size(a,2)) = 0.0d0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(1:staggered_lhalo, :) = 0.0d0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(size(a,1)-staggered_uhalo:size(a,1), :) = 0.0d0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, 1:staggered_lhalo) = 0.0d0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, size(a,2)-staggered_uhalo:size(a,2)) = 0.0d0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:staggered_lhalo, :staggered_lhalo) = 0.0d0
       if (southeast_corner) a(size(a,1)-staggered_uhalo:size(a,1), :staggered_lhalo) = 0.0d0
       if (northeast_corner) a(size(a,1)-staggered_uhalo:size(a,1), &
                               size(a,2)-staggered_uhalo:size(a,2)) = 0.0d0
       if (northwest_corner) a(:staggered_lhalo, size(a,2)-staggered_uhalo:size(a,2)) = 0.0d0

    endif   ! outflow or no_ice_bc

    end associate

  end subroutine staggered_parallel_halo_real8_2d


  subroutine staggered_parallel_halo_real8_3d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a real(dp) 3D field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    real(dp),dimension(size(a,1), staggered_lhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    real(dp),dimension(size(a,1), staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    real(dp),dimension(size(a,1),size(a,2),staggered_lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1),size(a,2),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &
         north     => parallel%north,          &
         south     => parallel%south,          &
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,2)/=local_ewn-1 .or. size(a,3)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    ! Send and receive east-west messages

    if (this_rank > west .or. periodic_bc) then
      wsend(:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
       esend(:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:, size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(:, size(a,2)-staggered_uhalo+1:size(a,2), staggered_jlo:staggered_jhi) = &
          erecv(:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(:, 1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    ! Send and receive north-south messages

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:,:) = &
        a(:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:,:) = &
        a(:,:,size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,:,size(a,3)-staggered_uhalo+1:size(a,3)) = nrecv(:,:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,:,1:staggered_lhalo) = srecv(:,:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then
          a(:, 1:staggered_lhalo-1, :) = 0.0d0
       endif

       if (this_rank >= east) then
          a(:, size(a,2)-staggered_uhalo+1:size(a,2), :) = 0.0d0
       endif

       if (this_rank <= south) then
          a(:, :, 1:staggered_lhalo-1) = 0.0d0
       endif

       if (this_rank >= north) then
          a(:, :, size(a,3)-staggered_uhalo+1:size(a,3)) = 0.0d0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(:, 1:staggered_lhalo, :) = 0.0d0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(:, size(a,2)-staggered_uhalo:size(a,2), :) = 0.0d0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, :, 1:staggered_lhalo) = 0.0d0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, :, size(a,3)-staggered_uhalo:size(a,3)) = 0.0d0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:, :staggered_lhalo, :staggered_lhalo) = 0.0d0
       if (southeast_corner) a(:, size(a,2)-staggered_uhalo:size(a,2), :staggered_lhalo) = 0.0d0
       if (northeast_corner) a(:, size(a,2)-staggered_uhalo:size(a,2), &
                               size(a,3)-staggered_uhalo:size(a,3)) = 0.0d0
       if (northwest_corner) a(:, :staggered_lhalo, size(a,3)-staggered_uhalo:size(a,3)) = 0.0d0

    endif   ! outflow or no_ice_bc

    end associate

  end subroutine staggered_parallel_halo_real8_3d


  subroutine staggered_parallel_halo_real8_4d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a real(dp) 4D field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    real(dp),dimension(size(a,1), size(a,2), staggered_lhalo,&
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    real(dp),dimension(size(a,1), size(a,2), staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    real(dp),dimension(size(a,1),size(a,2),size(a,3),staggered_lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1),size(a,2),size(a,3),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &
         north     => parallel%north,          &
         south     => parallel%south,          &
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,3)/=local_ewn-1 .or. size(a,4)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      wsend(:,:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:,:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
       esend(:,:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:,:, size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(:,:, size(a,3)-staggered_uhalo+1:size(a,3), staggered_jlo:staggered_jhi) = &
          erecv(:,:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(:,:, 1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:,:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:,:,:) = &
        a(:,:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:,:,:) = &
        a(:,:,:,size(a,4)-staggered_uhalo-staggered_lhalo+1:size(a,4)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,:,:,size(a,4)-staggered_uhalo+1:size(a,4)) = nrecv(:,:,:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,:,:,1:staggered_lhalo) = srecv(:,:,:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then
          a(:,:, 1:staggered_lhalo-1, :) = 0.0d0
       endif

       if (this_rank >= east) then
          a(:,:, size(a,3)-staggered_uhalo+1:size(a,3), :) = 0.0d0
       endif

       if (this_rank <= south) then
          a(:,:, :, 1:staggered_lhalo-1) = 0.0d0
       endif

       if (this_rank >= north) then
          a(:,:, :, size(a,4)-staggered_uhalo+1:size(a,4)) = 0.0d0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(:,:, 1:staggered_lhalo, :) = 0.0d0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(:,:, size(a,3)-staggered_uhalo:size(a,3), :) = 0.0d0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:,:, :, 1:staggered_lhalo) = 0.0d0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:,:, :, size(a,4)-staggered_uhalo:size(a,4)) = 0.0d0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:,:, :staggered_lhalo, :staggered_lhalo) = 0.0d0
       if (southeast_corner) a(:,:, size(a,3)-staggered_uhalo:size(a,3), :staggered_lhalo) = 0.0d0
       if (northeast_corner) a(:,:, size(a,3)-staggered_uhalo:size(a,3), &
                               size(a,4)-staggered_uhalo:size(a,4)) = 0.0d0
       if (northwest_corner) a(:,:, :staggered_lhalo, size(a,4)-staggered_uhalo:size(a,4)) = 0.0d0

    endif   ! outflow or no_ice_bc

    end associate

  end subroutine staggered_parallel_halo_real8_4d

!=======================================================================

  ! subroutines belonging to the staggered_parallel_halo_extrapolate interface

  subroutine staggered_parallel_halo_extrapolate_integer_2d(a, parallel)

    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel
    integer :: i, j

    ! begin

    associate(  &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         east   => parallel%east,    &
         west   => parallel%west,    &
         north  => parallel%north,   &
         south  => parallel%south    &
       )

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Extrapolate the staggered field into halo cells along the global boundary.
    ! Currently this is used only for kinbcmask.
    ! Note: The extrapolation region includes locally owned cells along
    !       the north and east boundaries of the global domain.

    ! First update the halos so that we are sure the interior halos are correct
    call staggered_parallel_halo(a, parallel)

    ! MJH Note: Modified code to now copy entire east and west columns rather than
    !  just the owned cells in those columns.  This avoids having the halos have
    !  potentially wrong information (i.e., a few cells in the corner don't get extrapolated into)

    if (this_rank >= east) then  ! at east edge of global domain
       ! extrapolate eastward
       do i = size(a,1)-staggered_uhalo, size(a,1)
          a(i, :) = a(size(a,1)-staggered_uhalo-1, :)
       enddo
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! extrapolate westward
       do i = 1, staggered_lhalo
          a(i, :) = a(staggered_lhalo+1, :)
       enddo
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       ! extrapolate northward
       do j = size(a,2)-staggered_uhalo, size(a,2)
          a(:, j) = a(:, size(a,2)-staggered_uhalo-1)
       enddo
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! extrapolate southward
       do j = 1, staggered_lhalo
          a(:, j) = a(:, staggered_lhalo+1)
       enddo
    endif

    end associate

  end subroutine staggered_parallel_halo_extrapolate_integer_2d


  subroutine staggered_parallel_halo_extrapolate_real8_2d(a, parallel)

    implicit none
    real(dp), dimension(:,:) :: a
    type(parallel_type) :: parallel
    integer :: i, j

    ! begin

    associate(  &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         east  => parallel%east,   &
         west  => parallel%west,   &
         north => parallel%north,  &
         south => parallel%south   &
       )

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(*,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Extrapolate the staggered field into halo cells along the global boundary.
    ! Currently this is used only for kinbcmask.
    ! Note: The extrapolation region includes locally owned cells along
    !       the north and east boundaries of the global domain.

    ! First update the halos so that we are sure the interior halos are correct
    call staggered_parallel_halo(a, parallel)

    ! MJH Note: Modified code to now copy entire east and west columns rather than
    !  just the owned cells in those columns.  This avoids having the halos have
    !  potentially wrong information (i.e., a few cells in the corner don't get extrapolated into)

! Useful for debugging small domains (the YYYY is just a tag for grepping the output, particularly if you prepend the processor number, e.g. "0YYYY")
!  do j = 1, size(a,2)
!     write(6, "(i3, 'YYYY BEFORE row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    if (this_rank >= east) then  ! at east edge of global domain
       ! extrapolate eastward
       do i = size(a,1)-staggered_uhalo, size(a,1)
          a(i, :) = a(size(a,1)-staggered_uhalo-1, :)
       enddo
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! extrapolate westward
       do i = 1, staggered_lhalo
          a(i, :) = a(staggered_lhalo+1, :)
       enddo
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       ! extrapolate northward
       do j = size(a,2)-staggered_uhalo, size(a,2)
          a(:, j) = a(:, size(a,2)-staggered_uhalo-1)
       enddo
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! extrapolate southward
       do j = 1, staggered_lhalo
          a(:, j) = a(:, staggered_lhalo+1)
       enddo
    endif

! Useful for debugging small domains
!  do j = 1, size(a,2)
!     write(6, "(i3, 'YYYY AFTER  row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo
    end associate

  end subroutine staggered_parallel_halo_extrapolate_real8_2d

!=======================================================================

!  Following routines imported from the Community Earth System Model
!  (models/utils/mct/mpeu.m_FcComms.F90)
!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gather_int - Gather an array of type integer
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em integer} 
! to the {\tt root} process. Explicit handshaking messages are used
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gather
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!

  subroutine fc_gather_int (sendbuf, sendcnt, sendtype, &
                            recvbuf, recvcnt, recvtype, &
                            root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      integer,               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer,               intent(in)  :: recvcnt
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      integer,               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   integer :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, i, count, displs
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = 1
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               if (recvcnt > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  displs = p*recvcnt
                  call mpi_irecv ( recvbuf(displs+1), recvcnt, &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         displs = mytid*recvcnt
         do i=1,sendcnt
            recvbuf(displs+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

      endif

   else
 
      call mpi_gather (sendbuf, sendcnt, sendtype, &
                       recvbuf, recvcnt, recvtype, &
                       root, comm, ier)
   endif

   return
  end subroutine fc_gather_int

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_int - Gather an array of type integer
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em integer} 
! to the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_int (sendbuf, sendcnt, sendtype, &
                              recvbuf, recvcnts, displs, recvtype, &
                              root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      integer,               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer, dimension(:), intent(in)  :: recvcnts
      integer, dimension(:), intent(in)  :: displs
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      integer,               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   integer :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = 1
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

      endif

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)

   endif

   return

  end subroutine fc_gatherv_int

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_real4 - Gather an array of type real*4
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em real*4} to
! the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_real4 (sendbuf, sendcnt, sendtype, &
                                recvbuf, recvcnts, displs, recvtype, &
                                root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      real(sp),               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer, dimension(:), intent(in)  :: recvcnts
      integer, dimension(:), intent(in)  :: displs
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      real(sp),               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   real(sp) :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = 1.0
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

      endif

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)

   endif

   return

  end subroutine fc_gatherv_real4

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_real8 - Gather an array of type real*8
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em real*8} to
! the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_real8 (sendbuf, sendcnt, sendtype, &
                                recvbuf, recvcnts, displs, recvtype, &
                                root, comm, flow_cntl)
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      real(dp),               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer, dimension(:), intent(in)  :: recvcnts
      integer, dimension(:), intent(in)  :: displs
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      real(dp),               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   real(dp) :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = 1.0
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

      endif

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)

   endif

   return

  end subroutine fc_gatherv_real8

!BOP -------------------------------------------------------------------
!
! !IROUTINE: fc_gatherv_log - Gather an array of type logical
!
! !DESCRIPTION:
! This routine gathers a {\em distributed} array of type {\em logical} 
! to the {\tt root} process. Explicit handshaking messages are uesd
! to control the number of processes communicating with the root
! at any one time.
!
! If flow_cntl optional parameter 
!    < 0 : use MPI_Gatherv
!    >= 0: use point-to-point with handshaking messages and 
!          preposting receive requests up to 
!          max(min(1,flow_cntl),max_gather_block_size) 
!          ahead if optional flow_cntl parameter is present.
!          Otherwise, fc_gather_flow_cntl is used in its place.
!    Default value is max_gather_block_size.
! !INTERFACE:
!
   subroutine fc_gatherv_log (sendbuf, sendcnt, sendtype, &
                              recvbuf, recvcnts, displs, recvtype, &
                              root, comm, flow_cntl )
!
! !USES:
!
      use mpi_mod

!
! !INPUT PARAMETERS: 
!
      logical,               intent(in)  :: sendbuf(*)
      integer,               intent(in)  :: sendcnt
      integer,               intent(in)  :: sendtype
      integer, dimension(:), intent(in)  :: recvcnts
      integer, dimension(:), intent(in)  :: displs
      integer,               intent(in)  :: recvtype
      integer,               intent(in)  :: root
      integer,               intent(in)  :: comm
      integer, optional,     intent(in)  :: flow_cntl

! !OUTPUT PARAMETERS: 
!
      logical,               intent(out) :: recvbuf(*)

!EOP ___________________________________________________________________

   logical :: signal
   logical :: fc_gather         ! use explicit flow control?
   integer :: gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer :: ier ! MPI error code

   signal = .true.
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
        fc_gather = .false.
      endif
   else
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
   endif

   if (fc_gather) then
 
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

         ! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, recvtype, p, mtag, comm, ier )
               end if
            end if
         end do

         ! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

         ! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, sendtype, root, mtag, comm, &
                            status, ier )
            call mpi_send ( sendbuf, sendcnt, sendtype, root, mtag, &
                            comm, ier )
         end if

     endif

   else
 
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)

   endif

   return

  end subroutine fc_gatherv_log

  !=======================================================================

end module cism_parallel

!=======================================================================

