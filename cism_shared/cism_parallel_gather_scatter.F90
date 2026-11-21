!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_parallel_gather_scatter.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_parallel_gather_scatter

  ! Part of the cism_parallel module:
  !  gathers to and scatters from the main task, distributed reads and writes
  !  of netCDF variables, and printing of distributed arrays
  ! Other modules use these through cism_parallel.

  use netcdf
  use glimmer_global, only : dp, sp
  use glimmer_paramets, only: iulog
  use cism_parallel_core
  use cism_parallel_broadcast

  implicit none

  ! common work space
  integer, dimension(4), save :: d_gs_mybounds
  integer, dimension(:,:), allocatable,save :: d_gs_bounds

  ! distributed gather flow control parameter
  integer,parameter :: max_gather_block_size = 64 ! max and default

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

  interface gather_var
     module procedure gather_var_integer_2d
     module procedure gather_var_logical_2d
     module procedure gather_var_real4_2d
     module procedure gather_var_real4_3d
     module procedure gather_var_real8_2d
     module procedure gather_var_real8_3d
  end interface

  interface gather_var_row
     module procedure gather_var_row_real8_2d
  end interface

  interface gather_all_var_row
     module procedure gather_all_var_row_real8_2d
  end interface

  interface gather_var_col
     module procedure gather_var_col_real8_2d
  end interface

  interface gather_all_var_col
     module procedure gather_all_var_col_real8_2d
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

  interface scatter_var
     module procedure scatter_var_integer_2d
     module procedure scatter_var_logical_2d
     module procedure scatter_var_real4_2d
     module procedure scatter_var_real4_3d
     module procedure scatter_var_real8_2d
     module procedure scatter_var_real8_3d
  end interface

  interface scatter_var_row
     module procedure scatter_var_row_real8_2d
  end interface

  interface scatter_var_col
     module procedure scatter_var_col_real8_2d
  end interface

  interface parallel_print
     module procedure parallel_print_integer_2d
     module procedure parallel_print_real8_2d
     module procedure parallel_print_real8_3d
  end interface

contains

!=======================================================================

  ! subroutines belonging to the gather_var interface

  ! WHL, July 2019:
  ! There is an issue with allocating the global_values array in the gather_var_*,
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

  subroutine gather_var_integer_2d(values, global_values, parallel)

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
       write(*,*) "gather does not currently work for"
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
       if (allocated(global_values)) deallocate(global_values)
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
       if (allocated(global_values)) deallocate(global_values)
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

  end subroutine gather_var_integer_2d


  subroutine gather_var_logical_2d(values, global_values, parallel)

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
       write(*,*) "gather does not currently work for"
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
       if (allocated(global_values)) deallocate(global_values)
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
       if (allocated(global_values)) deallocate(global_values)
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

  end subroutine gather_var_logical_2d


  subroutine gather_var_real4_2d(values, global_values, parallel)

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
       write(*,*) "gather does not currently work for"
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
       if (allocated(global_values)) deallocate(global_values)
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
       if (allocated(global_values)) deallocate(global_values)
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

  end subroutine gather_var_real4_2d


  subroutine gather_var_real4_3d(values, global_values, parallel, ld1, ud1)

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
       write(*,*) "gather does not currently work for"
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
       if (allocated(global_values)) deallocate(global_values)
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
       if (allocated(global_values)) deallocate(global_values)
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

  end subroutine gather_var_real4_3d


  subroutine gather_var_real8_2d(values, global_values, parallel)

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
       write(*,*) "gather does not currently work for"
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
       if (allocated(global_values)) deallocate(global_values)
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
       if (allocated(global_values)) deallocate(global_values)
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

  end subroutine gather_var_real8_2d


  subroutine gather_var_real8_3d(values, global_values, parallel, ld1, ud1)

    ! Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:),intent(in) :: values    ! i and j are indices 2 and 3
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

    if (uhalo==0 .and. size(values,2)==local_ewn-1) then
       ! Fixing this would require some generalization as is done for distributed_put_var
       write(*,*) "gather does not currently work for"
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
       if (allocated(global_values)) deallocate(global_values)
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
       if (allocated(global_values)) deallocate(global_values)
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

  end subroutine gather_var_real8_3d

!=======================================================================

  ! subroutines belonging to the gather_var_row interface

  subroutine gather_var_row_real8_2d(values, global_values, parallel)

    ! Gather data along a row of tasks onto the main task for that row.
    ! Based on gather_var_real8_2d.
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
       write(*,*) "ERROR: gather_var_row requires N-S array size of own_nsn"
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
!             write(iulog,*) this_rank, 'Recompute d_gs_bounds_row'
!             write(iulog,*) '  current size =', (d_gs_bounds_row(2,:) - d_gs_bounds_row(1,:) + 1)
!             write(iulog,*) '  size(values,1) =', size(values,1)
!          endif
       endif
    else
       new_bounds_row = .true.
!       if (main_task_row) write(iulog,*) this_rank, 'Allocate d_gs_bounds_row'
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

  end subroutine gather_var_row_real8_2d

!=======================================================================

  ! subroutines belonging to the gather_all_var_row interface

  subroutine gather_all_var_row_real8_2d(values, global_values, parallel)

    ! Gather global data along a row of tasks onto each task for that row.
    ! Based on gather_var_real8_2d.
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
       write(*,*) "ERROR: gather_var_row requires N-S array size of own_nsn"
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
!             write(iulog,*) this_rank, 'Recompute d_gs_bounds_row'
!             write(iulog,*) '  current size =', (d_gs_bounds_row(2,:) - d_gs_bounds_row(1,:) + 1)
!             write(iulog,*) '  size(values,1) =', size(values,1)
!          endif
       endif
    else
       new_bounds_row = .true.
!       if (main_task_row) write(iulog,*) this_rank, 'Allocate d_gs_bounds_row'
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

  end subroutine gather_all_var_row_real8_2d

!=======================================================================

  ! subroutines belonging to the gather_var_col interface

  subroutine gather_var_col_real8_2d(values, global_values, parallel)

    ! Gather data along a column of tasks onto the main task for that column.
    ! Based on gather_var_real8_2d.
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
       write(*,*) "ERROR: gather_var_row requires E-W array size of own_ewn"
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
!             write(iulog,*) this_rank, 'Recompute d_gs_bounds_col'
!             write(iulog,*) '  current size =', (d_gs_bounds_col(2,:) - d_gs_bounds_col(1,:) + 1)
!             write(iulog,*) '  size(values,1) =', size(values,1)
!          endif
       endif
    else
       new_bounds_col = .true.
!       if (main_task_col) write(iulog,*) this_rank, 'Allocate d_gs_bounds_col'
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

  end subroutine gather_var_col_real8_2d

!=======================================================================

  ! subroutines belonging to the gather_all_var_col interface

  subroutine gather_all_var_col_real8_2d(values, global_values, parallel)

    ! Gather global data along a column of tasks onto each task for that column.
    ! Based on gather_var_real8_2d.
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
       write(*,*) "ERROR: gather_var_row requires E-W array size of own_ewn"
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
!             write(iulog,*) this_rank, 'Recompute d_gs_bounds_col'
!             write(iulog,*) '  current size =', (d_gs_bounds_col(2,:) - d_gs_bounds_col(1,:) + 1)
!             write(iulog,*) '  size(values,1) =', size(values,1)
!          endif
       endif
    else
       new_bounds_col = .true.
!       if (main_task_col) write(iulog,*) this_rank, 'Allocate d_gs_bounds_col'
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

  end subroutine gather_all_var_col_real8_2d

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

    integer :: i,ierror,myn,status,x0id,y0id,x1id,y1id
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

    ! Initialize the ids to an invalid value; nf90_inq_varid does not set them
    ! when a variable is missing from the file (x0 and y0 often are)
    x0id = -1
    y0id = -1
    x1id = -1
    y1id = -1

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x0id)
    call broadcast(y0id)
    call broadcast(x1id)
    call broadcast(y1id)
    if (varid==x0id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub-1
       myn = global_ewn-1
    else if (varid==y0id) then
       mybounds(1) = nslb
       mybounds(2) = nsub-1
       myn = global_nsn-1
    else if (varid==x1id) then
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
       if (varid==x0id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub-1))
       elseif (varid==y0id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub-1))
       elseif (varid==x1id) then
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

    integer :: i,ierror,myn,status,x0id,y0id,x1id,y1id
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

    ! Initialize the ids to an invalid value; nf90_inq_varid does not set them
    ! when a variable is missing from the file (x0 and y0 often are)
    x0id = -1
    y0id = -1
    x1id = -1
    y1id = -1

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x0id)
    call broadcast(y0id)
    call broadcast(x1id)
    call broadcast(y1id)
    if (varid==x0id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub-1
       myn = global_ewn-1
    else if (varid==y0id) then
       mybounds(1) = nslb
       mybounds(2) = nsub-1
       myn = global_nsn-1
    else if (varid==x1id) then
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
       if (varid==x0id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub-1))
       elseif (varid==y0id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub-1))
       elseif (varid==x1id) then
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
       write(iulog,*) "distributed_print does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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
       write(iulog,*) "distributed_print does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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
       write(iulog,*) "distributed_print does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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

    integer :: i,ierror,myn,status,x0id,y0id,x1id,y1id
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

    ! Initialize the ids to an invalid value; nf90_inq_varid does not set them
    ! when a variable is missing from the file (x0 and y0 often are)
    x0id = -1
    y0id = -1
    x1id = -1
    y1id = -1

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x0id)
    call broadcast(y0id)
    call broadcast(x1id)
    call broadcast(y1id)
    if (varid==x0id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub-1
       myn = global_ewn-1
    else if (varid==y0id) then
       mybounds(1) = nslb
       mybounds(2) = nsub-1
       myn = global_nsn-1
    else if (varid==x1id) then
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
       if (varid==x0id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub-1))
       elseif (varid==y0id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub-1))
       elseif (varid==x1id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub))
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

    integer :: i,ierror,myn,status,x0id,y0id,x1id,y1id

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

    ! Initialize the ids to an invalid value; nf90_inq_varid does not set them
    ! when a variable is missing from the file (x0 and y0 often are)
    x0id = -1
    y0id = -1
    x1id = -1
    y1id = -1

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x0id)
    call broadcast(y0id)
    call broadcast(x1id)
    call broadcast(y1id)
    if (varid==x0id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub-1
       myn = global_ewn-1
    else if (varid==y0id) then
       mybounds(1) = nslb
       mybounds(2) = nsub-1
       myn = global_nsn-1
    else if (varid==x1id) then
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
       if (varid==x0id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub-1))
       elseif (varid==y0id) then
          allocate(global_values(global_minval_nslb:global_maxval_nsub-1))
       elseif (varid==x1id) then
          allocate(global_values(global_minval_ewlb:global_maxval_ewub))
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

  ! subroutines belonging to the scatter_var interface

  subroutine scatter_var_integer_2d(values, global_values, parallel)

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
       write(iulog,*) "scatter does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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

  end subroutine scatter_var_integer_2d


  subroutine scatter_var_logical_2d(values, global_values, parallel)

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
       write(iulog,*) "scatter does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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

  end subroutine scatter_var_logical_2d


  subroutine scatter_var_real4_2d(values, global_values, parallel)

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
       write(iulog,*) "scatter does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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

  end subroutine scatter_var_real4_2d


  subroutine scatter_var_real4_3d(values, global_values, parallel)

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
       write(iulog,*) "scatter does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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

  end subroutine scatter_var_real4_3d


  subroutine scatter_var_real8_2d(values, global_values, parallel)

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
       write(iulog,*) "scatter does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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

  end subroutine scatter_var_real8_2d


  subroutine scatter_var_real8_3d(values, global_values, parallel, deallocflag)

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
       write(iulog,*) "scatter does not currently work for"
       write(iulog,*) "variables on the staggered grid when uhalo=0"
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

  end subroutine scatter_var_real8_3d

!=======================================================================

  ! subroutines belonging to the scatter_var_row interface

  subroutine scatter_var_row_real8_2d(values, global_values, parallel)

    ! Scatter data to a row of tasks from the main task for that row.
    ! Based on scatter_var_real8_2d.
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
       write(iulog,*) "ERROR: scatter_var_row requires N-S array size of own_nsn"
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

  end subroutine scatter_var_row_real8_2d

!=======================================================================

  ! subroutines belonging to the scatter_var_col interface

  subroutine scatter_var_col_real8_2d(values, global_values, parallel)

    ! Scatter data to a column of tasks from the main task for that column
    ! Based on scatter_var_real8_2d.
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
       write(iulog,*) "ERROR: scatter_var_col requires E-W array size of own_nsn"
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

  end subroutine scatter_var_col_real8_2d

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

end module cism_parallel_gather_scatter

!=======================================================================
