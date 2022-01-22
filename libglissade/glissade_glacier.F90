!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_glacier.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glissade_glacier

    ! Subroutines for glacier tuning and tracking

    use glimmer_global 
    use glimmer_paramets, only: thk0, len0
    use glide_types
    use glimmer_log
    use cism_parallel, only: main_task, this_rank, nhalo

    implicit none

    private
    public :: glissade_glacier_init

    logical, parameter :: verbose_glacier = .true.

    ! derived type that holds info for each glaciated grid cell
    type glacier_info
       integer :: id           ! input glacier ID, usually RGI
       integer :: indxi        ! i index of cell
       integer :: indxj        ! j index of cell
    end type glacier_info

contains

!****************************************************      

  subroutine glissade_glacier_init(model)

    ! Initialize glaciers for a region
    ! If running on multiple disconnected glacier regions, this routine should be called once per region.
    !TODO: One set of logic for init, another for restart

    ! One key task is to create one-to-one maps between the input glacier_id array (typically with RGI IDs)
    !  and a local array called glacier_id_cism.  The local array assigns to each grid cell
    !  a number between 1 and nglacier where nglacier is the total number of unique glacier IDs.
    ! This allows us to loop over IDs in the range (1:nglacier), which is more efficient than
    !  looping over input glacier IDs.  The input IDs typically have large gaps.

    use cism_parallel, only: distributed_gather_var, distributed_scatter_var, &
         parallel_reduce_sum, broadcast, parallel_halo 

    type(glide_global_type),intent(inout) :: model

    ! local variables
    integer :: ewn, nsn, global_ewn, global_nsn
    integer :: itest, jtest, rtest   ! coordinates of diagnostic point

    ! temporary global arrays
    integer, dimension(:,:), allocatable :: &
         glacier_id_global,         & ! global array of the input glacier ID; maps (i,j) to RGI ID
         glacier_id_cism_global       ! global array of the CISM glacier ID; maps (i,j) to CISM glacier ID

    type(glacier_info), dimension(:), allocatable :: & 
         glacier_list                 ! sorted list of glacier IDs with i and j indices

    ! The next three arrays will have dimension (nglacier), once nglacier is computed
!!    integer, dimension(:), allocatable :: &
!!         cism_to_glacier_id           ! maps CISM ID (1:nglacier) to input glacier_id

    real(dp), dimension(:), allocatable :: &
         local_area,                & ! area per glacier (m^2)
         local_volume                 ! volume per glacier (m^3)

    integer :: &
         nglacier,                  & ! number of glaciers in global domain
         ncells_glacier,            & ! number of global grid cells occupied by glaciers at initialization
         current_id,                & ! current glacier_id from list
         gid_minval, gid_maxval       ! min and max values of glacier_id

    type(parallel_type) :: parallel   ! info for parallel communication

    integer :: i, j, nc, ng, count

    !WHL - debug
    integer, dimension(:), allocatable :: test_list
    integer ::  nlist
    real(sp) :: random

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'In glissade_glacier_init'
    endif

    parallel = model%parallel
    global_ewn = parallel%global_ewn
    global_nsn = parallel%global_nsn

    ewn = model%general%ewn
    nsn = model%general%nsn

    ! get coordinates of diagnostic point
    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ! debug - scatter test
!    if (main_task) print*, 'Scatter glacier_id_cism'
!    allocate(glacier_id_cism_global(global_ewn,global_nsn))
!    glacier_id_cism_global = 0
!    model%glacier%glacier_id_cism = 0
!    call distributed_scatter_var(model%glacier%glacier_id_cism, glacier_id_cism_global, parallel)
!    if (main_task) print*, 'Successful scatter'
!    if (allocated(glacier_id_cism_global)) deallocate(glacier_id_cism_global)

    ! Gather glacier IDs to the main task

    allocate(glacier_id_global(global_ewn, global_nsn))
    call distributed_gather_var(model%glacier%glacier_id, glacier_id_global, parallel)

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'Gathered glacier IDs to main task'
       print*, 'size(glacier_id) =', size(model%glacier%glacier_id,1), size(model%glacier%glacier_id,2)
       print*, 'size(glacier_id_global) =', size(glacier_id_global,1), size(glacier_id_global,2)
    endif

    if (verbose_glacier .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Glacier ID, rtest, itest, jtest:', rtest, itest, jtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') model%glacier%glacier_id(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Count the number of cells with glaciers

    count = 0

    ! Loop over locally owned cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          if (model%glacier%glacier_id(i,j) > 0) then
             count = count + 1
          elseif (model%glacier%glacier_id(i,j) < 0) then  ! should not happen
             print*, 'glacier_id < 0: i, j, value =', i, j, model%glacier%glacier_id(i,j)
             stop   ! TODO - exit gracefully
          endif
       enddo
    enddo

    ncells_glacier = parallel_reduce_sum(count)

    ! Allocate a global array on the main task only.
    ! On other tasks, allocate a size 0 array, since distributed_scatter_var wants arrays allocated on all tasks.
    if (main_task) then
       allocate(glacier_id_cism_global(global_ewn,global_nsn))
       glacier_id_cism_global(:,:) = 0.0d0
    else
       allocate(glacier_id_cism_global(0,0))
    endif

    if (main_task) then

       gid_minval = minval(glacier_id_global)
       gid_maxval = maxval(glacier_id_global)

       if (verbose_glacier) then
          print*, 'Total ncells   =', global_ewn * global_nsn
          print*, 'ncells_glacier =', ncells_glacier
          print*, 'glacier_id minval, maxval =', gid_minval, gid_maxval
       endif

       ! Create an unsorted list of glacier IDs, with associated i and j indices.
       ! There is one entry per glacier-covered cell.

       allocate(glacier_list(ncells_glacier))
       glacier_list(:)%id = 0
       glacier_list(:)%indxi = 0
       glacier_list(:)%indxj = 0

       count = 0

       do j = 1, global_nsn
          do i = 1, global_ewn
             if (glacier_id_global(i,j) > 0) then
                count = count + 1
                glacier_list(count)%id = glacier_id_global(i,j)
                glacier_list(count)%indxi = i
                glacier_list(count)%indxj = j
             endif
          enddo
       enddo

       deallocate(glacier_id_global)  ! no longer needed after glacier_list is built

       ! Sort the list from low to high IDs.
       ! As the IDs are sorted, the i and j indices come along for the ride.
       ! When there are multiple cells with the same glacier ID, these cells are adjacent on the list.
       ! For example, suppose the initial list is (5, 9, 7, 6, 7, 10, 4, 1, 1, 3, 1).
       ! The sorted list would be (1, 1, 1, 3, 5, 7, 7, 7, 9, 10).

       call glacier_quicksort(glacier_list, 1, ncells_glacier)

       if (verbose_glacier) then
          print*, 'Sorted glacier IDs in ascending order'
          print*, ' '
          print*, 'icell, i, j, ID for a few cells:'
          do i = 1, 10
             print*, i, glacier_list(i)%indxi, glacier_list(i)%indxj, glacier_list(i)%id
          enddo
          do i = ncells_glacier-9, ncells_glacier
             print*, i, glacier_list(i)%indxi, glacier_list(i)%indxj, glacier_list(i)%id
          enddo
       endif

!       WHL - Short list to test quicksort for integer arrays
!       print*, ' '
!       print*, 'Unsorted list:'
!       nlist = 20
!       allocate(test_list(nlist))
!       do i = 1, nlist
!          call random_number(random)
!          test_list(i) = int(random*nlist) + 1
!          print*, i, random, test_list(i)
!       enddo
!       call quicksort(test_list, 1, nlist)
!       print*, 'Sorted list:', test_list(:)

       ! Now that the glacier IDs are sorted from low to high,
       ! it is easy to count the total number of glaciers

       nglacier = 0
       current_id = 0
       do nc = 1, ncells_glacier
          if (glacier_list(nc)%id > current_id) then
             nglacier = nglacier + 1
             current_id = glacier_list(nc)%id
          endif
       enddo

       model%glacier%nglacier = nglacier

       ! Create two useful arrays:
       ! (1) The cism_to_glacier_id array maps the CISM ID (between 1 and nglacier) to the input glacier_id.
       ! (2) The glacier_id_cism array maps each glaciated grid cell (i,j) to a CISM ID.
       ! The reason to carry around i and j in the sorted glacier_list is to efficienly fill glacier_id_cism.
       ! Note: cism_to_glacier_id is part of the glacier derived type, but cannot be allocate until nglacier is known.

       allocate(model%glacier%cism_to_glacier_id(nglacier))
       model%glacier%cism_to_glacier_id(:) = 0

       if (verbose_glacier) then
          print*, ' '
          print*, 'Counted glaciers: nglacier =', nglacier
          print*, ' '
          print*, 'Pick a glacier: ng =',  nglacier/2
          print*, 'icell, i, j, glacier_id_cism_global(i,j), cism_to_glacier_id(ng)'
       endif

       ng = 0
       current_id = 0
       do nc = 1, ncells_glacier
          if (glacier_list(nc)%id > current_id) then
             ng = ng + 1
             current_id = glacier_list(nc)%id
             model%glacier%cism_to_glacier_id(ng) = glacier_list(nc)%id
          endif
          i = glacier_list(nc)%indxi
          j = glacier_list(nc)%indxj
          if (i == 0 .or. j == 0) then
             print*, 'Warning: zeroes, ng, i, j, id =', ng, i, j, glacier_list(nc)%id
             stop   ! TODO - exit gracefully
          endif
          glacier_id_cism_global(i,j) = ng
          if (ng == nglacier/2) then   ! random glacier
             print*, nc, i, j, glacier_id_cism_global(i,j), model%glacier%cism_to_glacier_id(ng)
          endif
          if (ng > nglacier) then
             print*, 'ng > nglacier, nc, i, j , ng =', nc, i, j, ng
             stop  !TODO - exit gracefully
          endif
       enddo

       deallocate(glacier_list)

       if (verbose_glacier) then
          print*, ' '
          print*, 'maxval(cism_to_glacier_id) =', maxval(model%glacier%cism_to_glacier_id) 
          print*, 'maxval(glacier_id_cism_global) =', maxval(glacier_id_cism_global)
       endif

    endif   ! main_task

    ! Communicate glacier info from the main task to all processors

    if (verbose_glacier .and. main_task) print*, 'Broadcast nglacier and cism_to_glacier_id'
    call broadcast(model%glacier%nglacier)
    nglacier = model%glacier%nglacier

    if (.not.associated(model%glacier%cism_to_glacier_id)) &
         allocate(model%glacier%cism_to_glacier_id(nglacier))
    call broadcast(model%glacier%cism_to_glacier_id)

    if (verbose_glacier .and. main_task) print*, 'Scatter glacier_id_cism'
    ! Note: glacier_id_cism_global is deallocated in the subroutine
    call distributed_scatter_var(model%glacier%glacier_id_cism, glacier_id_cism_global, parallel)
    call parallel_halo(model%glacier%glacier_id_cism, parallel)

    !TODO - Move area and volume computations to subroutines

    ! Allocate and initialize glacier area and volume

    allocate(model%glacier%area(nglacier))
    allocate(model%glacier%volume(nglacier))
    model%glacier%area(:) = 0.0d0
    model%glacier%volume(:) = 0.0d0

    allocate(local_area(nglacier))
    allocate(local_volume(nglacier))
    local_area(:) = 0.0d0
    local_volume(:) = 0.0d0

    ! Compute the initial area and volume of each glacier.
    ! We need parallel sums, since a glacier can lie on 2 or more processors.

    if (verbose_glacier .and. main_task) then
       print*, 'Compute glacier area and volume'
       print*, '   cell_area (m^3) =', model%geometry%cell_area(3,3) * len0**2
    endif

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = model%glacier%glacier_id_cism(i,j)
          if (ng >= 1) then
             local_area(ng) = local_area(ng) &
                  + model%geometry%cell_area(i,j)*len0**2
             local_volume(ng) = local_volume(ng) &
                  + model%geometry%cell_area(i,j)*len0**2 * model%geometry%thck(i,j)*thk0
          endif
       enddo
    enddo

    model%glacier%area   = parallel_reduce_sum(local_area)
    model%glacier%volume = parallel_reduce_sum(local_volume)

    if (verbose_glacier .and. main_task) then
       print*, 'Max area (km^2)   =', maxval(model%glacier%area) * 1.0d-6    ! m^2 to km^2
       print*, 'Max volume (km^3) =', maxval(model%glacier%volume) * 1.0d-9  ! m^3 to km^3
       print*, ' '
       print*, 'Selected A (km^2) and V (km^3) of large glaciers:'
       do ng = 1, nglacier
          if (model%glacier%area(ng) * 1.0d-6 > 10.0d0) then  ! 10 km^2 or more
             write(6,'(i8,2f10.3)') ng, model%glacier%area(ng)*1.0d-6, model%glacier%volume(ng)*1.0d-9
          endif
       enddo
    endif

    deallocate(local_area)
    deallocate(local_volume)

    if (main_task) print*, 'Done in glissade_glacier_init'

  end subroutine glissade_glacier_init

!****************************************************      

  recursive subroutine quicksort(A, first, last)
 
    ! Given an unsorted integer array, return an array with elements sorted from low to high.

    implicit none

    ! input/output arguments
    integer, dimension(:), intent(inout) :: A
    integer, intent(in) :: first, last
 
    ! local arguments
    integer :: temp
    integer :: pivot
    integer :: i, j

    pivot = A( (first+last)/2 )
    i = first
    j = last

    ! Partition loop
    do
       do while (A(i) < pivot)
          i = i + 1
       enddo
       do while (A(j) > pivot)
          j = j - 1
       enddo
       if (i >= j) exit
       temp = A(i)
       A(i) = A(j)
       A(j) = temp
       i = i + 1
       j = j - 1
    enddo

    if (first < i-1) call quicksort(A, first, i-1)
    if (last  > j+1) call quicksort(A, j+1, last)

!    print*, 'Done in quicksort'

  end subroutine quicksort

!****************************************************      

  recursive subroutine glacier_quicksort(A, first, last)
 
    ! Given an unsorted array of type glacier_info, return an array with
    ! glacier IDs (A%id) sorted from low to high.
    ! The logic is just like quicksort above, but tailored for the derived type.

    implicit none

    ! input/output arguments
    type(glacier_info), dimension(:), intent(inout) :: A
    integer, intent(in) :: first, last
 
    ! local arguments
    type(glacier_info) :: temp
    integer :: pivot
    integer :: i, j

    pivot = A( (first+last)/2 )%id
    i = first
    j = last

    ! Partition loop
    do
       do while (A(i)%id < pivot)
          i = i + 1
       enddo
       do while (A(j)%id > pivot)
          j = j - 1
       enddo
       if (i >= j) exit
       ! Swap A(i) with A(j). Note that A%indxi and A%indxj are swapped along with A%id.
       temp = A(i)
       A(i) = A(j)
       A(j) = temp
       i = i + 1
       j = j - 1
    enddo

    if (first < i-1) call glacier_quicksort(A, first, i-1)
    if (last  > j+1) call glacier_quicksort(A, j+1, last)

!    print*, 'Done in quicksort'

  end subroutine glacier_quicksort

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glissade_glacier

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
