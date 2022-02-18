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

!TODO:
! Set options for repeatedly reading the monthly climatological forcing
! Put a glacier section in the config file.
! Add restart logic in glissade_glacier_init.
! Decide on the list of glacier restart fields:
!    rgi_glacier_id, cism_glacier_id, glacier_area_target, glacier_volume_target,
!    glacier_mu_star, glacier_powerlaw_c
! What about nglacier?  Diagnose from size of restart arrays?
! What about ngdiag? Recompute?
! What about cism_to_rgi_glacier_id?  Recompute?
! What about array allocation?

module glissade_glacier

    ! Subroutines for glacier tuning and tracking

    use glimmer_global 
    use glimmer_paramets, only: thk0, len0
    use glimmer_physcon, only: scyr
    use glide_types
    use glimmer_log
    use cism_parallel, only: main_task, this_rank, nhalo

    implicit none

    private
    public :: verbose_glacier, glissade_glacier_init, &
         glissade_glacier_smb, glissade_glacier_inversion

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

    ! One key task is to create maps between input RGI glacier IDs (in the rgi_glacier_id array)
    !  and an array called cism_glacier_id.
    ! The cism_glacier_id array assigns to each grid cell (i,j) a number between 1 and nglacier,
    !  where nglacier is the total number of unique glacier IDs.
    ! This allows us to loop over IDs in the range (1:nglacier), which is more efficient than
    !  looping over the input glacier IDs, which often have large gaps.

    use cism_parallel, only: distributed_gather_var, distributed_scatter_var, &
         parallel_reduce_sum, broadcast, parallel_halo

    use cism_parallel, only: parallel_barrier  !WHL - debug

    type(glide_global_type),intent(inout) :: model

    ! local variables
    integer :: ewn, nsn               ! local grid dimensions
    integer :: global_ewn, global_nsn ! global grid dimensions
    integer :: itest, jtest, rtest    ! coordinates of diagnostic point
    real(dp) :: dew, dns              ! grid cell length in each direction (m)

    ! temporary global arrays
    integer, dimension(:,:), allocatable :: &
         rgi_glacier_id_global,     & ! global array of the input RGI glacier ID; maps (i,j) to RGI ID
         cism_glacier_id_global       ! global array of the CISM glacier ID; maps (i,j) to CISM glacier ID

    type(glacier_info), dimension(:), allocatable :: & 
         glacier_list                 ! sorted list of glacier IDs with i and j indices

    ! The next two arrays will have dimension (nglacier), once nglacier is computed
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
!    integer, dimension(:), allocatable :: test_list
!    integer ::  nlist
!    real(sp) :: random

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_glacier_init'
    endif

    parallel = model%parallel
    global_ewn = parallel%global_ewn
    global_nsn = parallel%global_nsn
    ewn = model%general%ewn
    nsn = model%general%nsn
    dew = model%numerics%dew
    dns = model%numerics%dns

    ! get coordinates of diagnostic point
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    if (verbose_glacier .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'RGI glacier ID, rtest, itest, jtest:', rtest, itest, jtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') model%glacier%rgi_glacier_id(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Arrays in the glacier derived type may have been allocated with dimension(1).
    ! If so, then deallocate here, and reallocate below with dimension (nglacier).
    ! Typically, nglacier is not known until after initialization.

    if (associated(model%glacier%glacierid)) deallocate(model%glacier%glacierid)
    if (associated(model%glacier%cism_to_rgi_glacier_id)) &
         deallocate(model%glacier%cism_to_rgi_glacier_id)
    if (associated(model%glacier%area)) deallocate(model%glacier%area)
    if (associated(model%glacier%volume)) deallocate(model%glacier%volume)
    if (associated(model%glacier%area_target)) deallocate(model%glacier%area_target)
    if (associated(model%glacier%volume_target)) deallocate(model%glacier%volume_target)
    if (associated(model%glacier%dvolume_dt)) deallocate(model%glacier%dvolume_dt)
    if (associated(model%glacier%mu_star)) deallocate(model%glacier%mu_star)
    if (associated(model%glacier%powerlaw_c)) deallocate(model%glacier%powerlaw_c)

    ! Count the number of cells with glaciers
    ! Loop over locally owned cells

    count = 0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          if (model%glacier%rgi_glacier_id(i,j) > 0) then
             count = count + 1
          elseif (model%glacier%rgi_glacier_id(i,j) < 0) then  ! should not happen
             print*, 'glacier_id < 0: i, j, value =', i, j, model%glacier%rgi_glacier_id(i,j)
             stop   ! TODO - exit gracefully
          endif
       enddo
    enddo

    ncells_glacier = parallel_reduce_sum(count)

    ! Gather the RGI glacier IDs to the main task
    if (main_task) allocate(rgi_glacier_id_global(global_ewn, global_nsn))
    call distributed_gather_var(model%glacier%rgi_glacier_id, rgi_glacier_id_global, parallel)

    ! Allocate a global array for the CISM glacier IDs on the main task..
    ! Allocate a size 0 array on other tasks; distributed_scatter_var wants arrays allocated on all tasks.                   

    if (main_task) then
       allocate(cism_glacier_id_global(global_ewn,global_nsn))
    else
       allocate(cism_glacier_id_global(0,0))
    endif
    cism_glacier_id_global(:,:) = 0.0d0

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'Gathered RGI glacier IDs to main task'
       print*, 'size(rgi_glacier_id) =', &
            size(model%glacier%rgi_glacier_id,1), size(model%glacier%rgi_glacier_id,2)
       print*, 'size(rgi_glacier_id_global) =', &
            size(rgi_glacier_id_global,1), size(rgi_glacier_id_global,2)
    endif

    if (main_task) then

       gid_minval = minval(rgi_glacier_id_global)
       gid_maxval = maxval(rgi_glacier_id_global)

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
             if (rgi_glacier_id_global(i,j) > 0) then
                count = count + 1
                glacier_list(count)%id = rgi_glacier_id_global(i,j)
                glacier_list(count)%indxi = i
                glacier_list(count)%indxj = j
             endif
          enddo
       enddo

       ! Deallocate the RGI global array (no longer needed after glacier_list is built)
       deallocate(rgi_glacier_id_global)  

       ! Sort the list from low to high IDs.
       ! As the IDs are sorted, the i and j indices come along for the ride.
       ! When there are multiple cells with the same glacier ID, these cells are adjacent on the list.
       ! For example, suppose the initial list is (5, 9, 7, 6, 7, 10, 4, 1, 1, 3, 1).
       ! The sorted list would be (1, 1, 1, 3, 4, 5, 6, 7, 7, 9, 10).

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

       ! Now that the glacier IDs are sorted from low to high, count the glaciers

       nglacier = 0
       current_id = 0
       do nc = 1, ncells_glacier
          if (glacier_list(nc)%id > current_id) then
             nglacier = nglacier + 1
             current_id = glacier_list(nc)%id
          endif
       enddo

       model%glacier%nglacier = nglacier

       ! Fill two useful arrays:
       ! (1) The cism_to_rgi_glacier_id array maps the CISM ID (between 1 and nglacier) to the RGI glacier_id.
       ! (2) The cism_glacier_id array maps each glaciated grid cell (i,j) to a CISM ID.
       ! By carrying i and j in the sorted glacier_list, we can efficiently fill cism_glacier_id.
       ! Note: cism_to_rgi_glacier_id cannot be allocated until nglacier is known.

       allocate(model%glacier%cism_to_rgi_glacier_id(nglacier))
       model%glacier%cism_to_rgi_glacier_id(:) = 0

       if (verbose_glacier) then
          print*, ' '
          print*, 'Counted glaciers: nglacier =', nglacier
          print*, ' '
          ng = nglacier/2
          print*, 'Random cism_glacier_id:', ng
          print*, 'icell, i, j, cism_glacier_id_global(i,j), cism_to_rgi_glacier_id(ng)'
       endif

       ng = 0
       current_id = 0
       do nc = 1, ncells_glacier
          if (glacier_list(nc)%id > current_id) then
             ng = ng + 1
             current_id = glacier_list(nc)%id
             model%glacier%cism_to_rgi_glacier_id(ng) = glacier_list(nc)%id
          endif
          i = glacier_list(nc)%indxi
          j = glacier_list(nc)%indxj
          if (i == 0 .or. j == 0) then
             print*, 'Warning: zeroes, ng, i, j, id =', ng, i, j, glacier_list(nc)%id
             stop   ! TODO - exit gracefully
          endif
          cism_glacier_id_global(i,j) = ng
          if (ng == nglacier/2) then   ! random glacier
             print*, nc, i, j, cism_glacier_id_global(i,j), model%glacier%cism_to_rgi_glacier_id(ng)
          endif
          if (ng > nglacier) then
             print*, 'ng > nglacier, nc, i, j , ng =', nc, i, j, ng
             stop  !TODO - exit gracefully
          endif
       enddo

       deallocate(glacier_list)

       if (verbose_glacier) then
          print*, ' '
          print*, 'maxval(cism_to_rgi_glacier_id) =', maxval(model%glacier%cism_to_rgi_glacier_id)
          print*, 'maxval(cism_glacier_id_global) =', maxval(cism_glacier_id_global)
       endif

    endif   ! main_task

    ! Scatter cism_glacier_id_global to all processors
    ! Note: This global array is deallocated in the distributed_scatter_var subroutine

    if (verbose_glacier .and. main_task) print*, 'Scatter cism_glacier_id'
    call distributed_scatter_var(model%glacier%cism_glacier_id, cism_glacier_id_global, parallel)
    call parallel_halo(model%glacier%cism_glacier_id, parallel)

    ! Broadcast glacier info from the main task to all processors

    if (verbose_glacier .and. main_task) print*, 'Broadcast nglacier and cism_to_rgi_glacier_id'
    call broadcast(model%glacier%nglacier)
    nglacier = model%glacier%nglacier

    if (.not.associated(model%glacier%cism_to_rgi_glacier_id)) &
         allocate(model%glacier%cism_to_rgi_glacier_id(nglacier))
    call broadcast(model%glacier%cism_to_rgi_glacier_id)

    ! Set the index of the diagnostic glacier, using the CISM glacier ID for the diagnostic point
    if (this_rank == rtest) then
       model%glacier%ngdiag = model%glacier%cism_glacier_id(itest,jtest)
    endif
    call broadcast(model%glacier%ngdiag, rtest)

    ! Allocate and fill the glacierid dimension array
    allocate(model%glacier%glacierid(nglacier))
    do ng = 1, nglacier
       model%glacier%glacierid(ng) = ng
    enddo

    ! Allocate other arrays with dimension(nglacier)
    allocate(model%glacier%area(nglacier))
    allocate(model%glacier%volume(nglacier))
    allocate(model%glacier%area_target(nglacier))
    allocate(model%glacier%volume_target(nglacier))
    allocate(model%glacier%dvolume_dt(nglacier))
    allocate(model%glacier%mu_star(nglacier))
    allocate(model%glacier%powerlaw_c(nglacier))

    ! Compute the initial area and volume of each glacier.
    ! These values will be targets for inversion.

    call glacier_area_volume(&
         ewn,           nsn,               &
         nglacier,                         &
         model%glacier%cism_glacier_id,    &
         dew*dns*len0**2,                  &
         model%geometry%thck*thk0,         &
         model%glacier%area,               &
         model%glacier%volume)

    ! Initialize the other glacier arrays

    model%glacier%area_target(:) = model%glacier%area(:)
    model%glacier%volume_target(:) = model%glacier%volume(:)
    model%glacier%dvolume_dt(:) = 0.0d0
    model%glacier%mu_star(:) = model%glacier%mu_star_const
    model%glacier%powerlaw_c(:) = model%basal_physics%powerlaw_c_const

    ! Check for zero A or V target
    if (main_task) then
       print*, ' '
       print*, 'Check for A = 0, V = 0'
       do ng = 1, nglacier
          if (model%glacier%area_target(ng) == 0.0d0 .or. &
              model%glacier%volume_target(ng) == 0.0d0) then
             print*, 'ng, A (km^2), V (km^3):', &
                  ng, model%glacier%area_target(ng)/1.0d6, model%glacier%volume_target(ng)/1.0d9
          endif
       enddo
    endif

    if (verbose_glacier .and. main_task) then
       print*, ' '
       ng = model%glacier%ngdiag
       print*, 'Glacier ID for diagnostic cell: r, i, j, ng =', rtest, itest, jtest, ng
       print*, 'area target (km^2) =', model%glacier%area_target(ng) / 1.0d6
       print*, 'volume target (km^3) =', model%glacier%volume_target(ng) / 1.0d9
!!       print*, 'dvolume_dt (km^3/yr) =', model%glacier%dvolume_dt(ng) * scyr/1.0d9
       print*, 'mu_star (mm/yr w.e./deg) =', model%glacier%mu_star(ng)
       print*, 'powerlaw_c (Pa (m/yr)^(-1/3)) =', model%glacier%powerlaw_c(ng)
       print*, 'Done in glissade_glacier_init'
    endif

  end subroutine glissade_glacier_init

!****************************************************

  subroutine glissade_glacier_smb(&
       ewn,              nsn,               &
       itest,   jtest,   rtest,             &
       nglacier,                            &
       cism_glacier_id,  mu_star,           &
       snow,             artm,              &
       glacier_smb)

    ! Compute the SMB in each grid cell using an empirical relationship
    !  based on Maussion et al. (2019):
    !
    !     SMB = snow - mu_star * max(artm - T_mlt, 0),
    !
    ! where snow = monthly mean snowfall rate,
    !       mu_star is a glacier-specific tuning parameter,
    !       atrm = monthly mean air temperature,
    !       Tmlt = monthly mean air temp above which melting occurs
    !
    ! This subroutine should be called at least once a month
    !
    ! Note: In Maussion et al., SMB and prcp are monthly mass balances in mm w.e.
    !       Not sure that mu_star should have the same units (though Fig. 3 shows
    !       units of mm w.e./yr/deg).

    use parallel, only: nhalo, main_task

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier,                    & ! total number of glaciers in the domain
         itest, jtest, rtest            ! coordinates of diagnostic point

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), dimension(nglacier), intent(in) :: &
         mu_star                        ! glacier-specific SMB tuning parameter (mm w.e./yr/deg)

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow,                        & ! monthly mean snowfall rate (mm w.e./yr)
         artm                           ! monthly mean 2m air temperature (deg C)

    real(dp), dimension(ewn,nsn), intent(out) :: &
         glacier_smb                    ! SMB in each gridcell (mm w.e./yr)

    ! local variables

    integer :: i, j, ng

    real(dp), parameter ::  &
         glacier_tmlt = -1.0d0          ! artm (deg C) above which melt occurs
                                        ! Maussion et al. suggest -1 C

    if (verbose_glacier .and. this_rank == rtest) then
       print*, 'In glissade_glacier_smb'
    endif

    ! initialize
    glacier_smb(:,:) = 0.0d0

    if (verbose_glacier .and. this_rank == rtest) then
       print*, 'Loop'
       print*, 'minval, maxval(snow) =', minval(snow), maxval(snow)
       print*, 'minval, maxval(artm) =', minval(artm), maxval(artm)
    endif

    ! compute SMB
    do j = 1, nsn
       do i = 1, ewn

          ng = cism_glacier_id(i,j)
          glacier_smb(i,j) = &
               snow(i,j) - mu_star(ng) * max(artm(i,j) - glacier_tmlt, 0.0d0)

          if (verbose_glacier .and. this_rank == rtest .and. i == itest .and. j == jtest) then
             print*, ' '
             print*, 'Glacier SMB: rank i, j =', this_rank, i, j
             print*, '   mu_star (mm/yr w.e./deg) =', mu_star(ng)
             print*, '   snow (mm/yr w.e.), artm (C) =', snow(i,j), artm(i,j)
             print*, '   SMB (mm/yr w.e.) =', glacier_smb(i,j)
          endif

       enddo
    enddo

    if (verbose_glacier .and. this_rank == rtest) then
       print*, 'Done in glissade_glacier_smb'
    endif

  end subroutine glissade_glacier_smb

!****************************************************

  subroutine glissade_glacier_inversion(&
       glacier_mu_star,                  &
       glacier_powerlaw_c,               &
       dt,                               &
       itest,   jtest,  rtest,           &
       ewn,             nsn,             &
       dew,             dns,             &
       thck,            dthck_dt,        &
       powerlaw_c_min,  powerlaw_c_max,  &
       glacier)

    use glimmer_paramets, only: len0, thk0
    use glimmer_physcon, only: scyr

    real(dp), intent(in) :: &
         dt,                          & ! time step (s)
         dew, dns                       ! grid cell dimensions (m)

    integer, intent(in) ::  &
         glacier_mu_star,             & ! flag for mu_star inversion
         glacier_powerlaw_c,          & ! flag for powerlaw_c inversion
         itest, jtest, rtest,         & ! coordinates of diagnostic cell
         ewn, nsn                       ! number of cells in each horizontal direction

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         thck,                        & ! ice thickness (m)
         dthck_dt                       ! rate of change of thickness (m/yr)

    real(dp), intent(in) :: &
         powerlaw_c_min, powerlaw_c_max ! min and max allowed values of C_p in power law (Pa (m/yr)^(-1/3))

    ! Note: The glacier type includes the following:
    ! integer ::  nglacier          ! number of glaciers in the global domain
    ! integer ::  ngdiag            ! CISM index of diagnostic glacier
    ! integer, dimension(:,:) :: cism_glacier_id  ! CISM glacier ID for each grid cell
    ! real(dp), dimension(:) :: area          ! glacier area (m^2)
    ! real(dp), dimension(:) :: volume        ! glacier volume (m^3)
    ! real(dp), dimension(:) :: dvolume_dt    ! rate of change of glacier volume (m^3/yr)
    ! real(dp), dimension(:) :: mu_star       ! SMB parameter for each glacier (mm/yr w.e./deg K)
    ! real(dp) :: mu_star_min, mu_star_max    ! min and max values allowed for mu_star
    ! real(dp), dimension(:) :: powerlaw_c    ! basal friction parameter for each glacier (Pa (m/yr)^(-1/3))

    type(glide_glacier), intent(inout) :: &
         glacier       ! glacier derived type

    ! local variables

    integer :: nglacier       ! number of glaciers
    integer :: ngdiag         ! CISM index of diagnostic glacier
    integer :: ng

    nglacier = glacier%nglacier
    ngdiag = glacier%ngdiag

    if (verbose_glacier .and. main_task) then
       print*, 'In glissade_glacier_inversion, dt (yr) =', dt
       print*, 'Diag cell (r, i, j) =', rtest, itest, jtest
       print*, '   thck (m), dthck(dt):', thck(itest, jtest), dthck_dt(itest, jtest)
       print*, 'call glacier_area_volume'
    endif

    ! Compute the current area and volume of each glacier
    ! Note: This requires global sums. For now, do the computation independently on each task.

    call glacier_area_volume(&
         ewn,           nsn,         &
         nglacier,                   &
         glacier%cism_glacier_id,    &
         dew*dns,                    &  ! m^2
         thck,                       &  ! m
         glacier%area,               &  ! m^2
         glacier%volume,             &  ! m^3
         dthck_dt,                   &  ! m/yr
         glacier%dvolume_dt)            ! m^3/yr

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'Update area (km^2) and volume (km^3) for glacier:', ngdiag
       print*, 'Current area and volume:', glacier%area(ngdiag)/1.0d6, &
            glacier%volume(ngdiag)/1.0d9
       print*, '   Target area and volume:', glacier%area_target(ngdiag)/1.0d6, &
            glacier%volume_target(ngdiag)/1.0d9
       print*, '   dV_dt (m^3/yr):', glacier%dvolume_dt(ngdiag)/1.0d9
    endif

    ! Given the current and target glacier areas, invert for mu_star

    if (glacier_mu_star == GLACIER_MU_STAR_INVERSION) then

       if (verbose_glacier .and. main_task) then
          print*, 'glacier_invert_mu_star'
       endif

       call glacier_invert_mu_star(&
            dt,                                         &
            ewn,                 nsn,                   &
            nglacier,            ngdiag,                &
            glacier%mu_star_min, glacier%mu_star_max,   &
            glacier%area,        glacier%area_target,   &
            glacier%mu_star)

    endif

    ! Given the current and target glacier volumes, invert for powerlaw_c

    if (glacier_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then

       if (verbose_glacier .and. main_task) then
          print*, 'glacier_invert_powerlaw_c'
       endif

       call glacier_invert_powerlaw_c(&
            dt,                                         &
            ewn,                 nsn,                   &
            nglacier,            ngdiag,                &
            powerlaw_c_min,      powerlaw_c_max,        &
            glacier%volume,      glacier%volume_target, &
            glacier%dvolume_dt,                         &
            glacier%powerlaw_c)

    endif

    if (verbose_glacier .and. main_task) then
       print*, 'Done in glacier_glacier_inversion'
    endif

  end subroutine glissade_glacier_inversion

!****************************************************

  subroutine glacier_invert_mu_star(&
       dt,                                  &
       ewn,              nsn,               &
       nglacier,         ngdiag,            &
       mu_star_min,      mu_star_max,       &
       area,             area_target,       &
       mu_star)

    ! Given the current glacier areas and area targets,
    ! invert for the parameter mu_star in the glacier SMB formula

    ! Note: This subroutine should be called from main_task only, since it uses
    !       glacier areas summed over all processors.

    ! input/output arguments

    real(dp), intent(in) :: &
         dt                             ! timestep (yr)

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier,                    & ! total number of glaciers in the domain
         ngdiag                         ! CISM ID of diagnostic glacier

    !TODO - Decide on max and min values.
    !       Min should be zero; don't want negative values

    real(dp), intent(in) :: &
         mu_star_min, mu_star_max       ! min and max allowed values of mu_star (mm w.e/yr/deg)

    real(dp), dimension(nglacier), intent(in) :: &
         area,                        & ! current glacier area (m^2)
         area_target                    ! area target (m^2)

    ! Note: Here, mu_star_glacier(nglacier) is the value shared by all cells in a given glacier
    ! The calling subroutine will need to map these values onto each grid cell.
    real(dp), dimension(nglacier), intent(inout) :: &
         mu_star                        ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)

    ! local variables

    integer :: ng

    real(dp), parameter :: &
         glacier_area_timescale = 100.d0  ! timescale (yr)

    real(dp) :: &
         err_area,                    & ! relative area error, (A - A_target)/A_target
         term1, term2,                & ! terms in prognostic equation for mu_star
         dmu_star                       ! change in mu_star

    character(len=100) :: message

    !TODO - Rewrite the comments below.
    ! I am going to try the inversion without a dA/dt term.
    ! This is because glacier area is going to change discontinuously
    !  as a glacier advances into or retreats from a given cell.

    ! The inversion works as follows:
    ! The change in mu_star is proportional to the current mu_star and to the relative error,
    !  err_area = (A - A_target)/A_target.
    ! If err_area > 0, we increase mu_star to make the glacier melt more and retreat.
    ! If err_area < 0, we reduce mu_star to make the glacier melt less and advance.
    ! This is done with a characteristic timescale tau.
    ! We also include a term proportional to dA/dt so that ideally, mu_star smoothly approaches
    !  the value needed to attain a steady-state A, without oscillating about the desired value.
    ! See the comments in module glissade_inversion, subroutine invert_basal_friction.
    ! We should always have mu_star >= 0.
    ! Maussion et al. (2019) suggest values of roughly 100 to 300 mm w.e./yr/deg,
    !  but with a wide range.
    ! (Wondering if values should be higher; seems like we should be able to get ~1000 mm melt
    !  in 0.1 year with (T - Tmlt) = 10 C.  This would imply mu_star = 1000 mm w.e./yr/deg.
    ! Here is the prognostic equation:
    ! dmu/dt = -mu_star * (1/tau) * (A - A_target)/A_target + (2*tau/A_target) * dA/dt

    do ng = 1, nglacier

       if (area_target(ng) > 0.0d0) then  ! this should be the case
          err_area = (area(ng) - area_target(ng)) / area_target(ng)
          term1 = -err_area / glacier_area_timescale
          dmu_star = mu_star(ng) * term1 * dt
!!          term2 = -2.0d0 * darea_dt(ng) / area_target(ng)
!!          dmu_star = mu_star(ng) * (term1 + term2) * dt

          ! Limit to prevent a large relative change in one step
          if (abs(dmu_star) > 0.05d0 * mu_star(ng)) then
             if (dmu_star > 0.0d0) then
                dmu_star =  0.05d0 * mu_star(ng)
             else
                dmu_star = -0.05d0 * mu_star(ng)
             endif
          endif

          ! Update mu_star
          mu_star(ng) = mu_star(ng) + dmu_star

          ! Limit to a physically reasonable range
          mu_star(ng) = min(mu_star(ng), mu_star_max)
          mu_star(ng) = max(mu_star(ng), mu_star_min)

          if (verbose_glacier .and. main_task .and. ng == ngdiag) then
             print*, ' '
             print*, 'Invert for mu_star: ngdiag =', ngdiag
             print*, 'A, A_target (km^2), err_area:', &
                  area(ng)/1.0d6, area_target(ng)/1.0d6, err_area
             print*, 'term1*dt:', term1*dt
             print*, 'dmu_star, new mu_star:', dmu_star, mu_star(ng)
          endif

       else   ! area_target(ng) = 0

          write(message,*) 'Error: area_target = 0 for glacier', ng
          call write_log(message, GM_FATAL)

       endif

    enddo   ! ng

  end subroutine glacier_invert_mu_star

!****************************************************

  subroutine glacier_invert_powerlaw_c(&
       dt,                                  &
       ewn,              nsn,               &
       nglacier,         ngdiag,            &
       powerlaw_c_min,   powerlaw_c_max,    &
       volume,           volume_target,     &
       dvolume_dt,       powerlaw_c)

    use glimmer_physcon, only: scyr

    ! Given the current glacier volumes and volume targets,
    ! invert for the parameter powerlaw_c in the relationship for basal sliding.

    ! Note: This subroutine should be called from main_task only, since it uses
    !       glacier volumes summed over all processors.

    ! input/output arguments

    real(dp), intent(in) :: &
         dt                             ! timestep (yr)

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier,                    & ! total number of glaciers in the domain
         ngdiag                         ! ID of diagnostic glacier

    real(dp), intent(in) :: &
         powerlaw_c_min, powerlaw_c_max ! min and max allowed values of powerlaw_c (Pa (m/yr)^(-1/3))

    real(dp), dimension(nglacier), intent(in) :: &
         volume,                      & ! current glacier volume (m^3)
         volume_target,               & ! volume target (m^3)
         dvolume_dt                     ! rate of change of volume (m^3/yr)

    ! Note: Here, powerlaw_c_glacier(nglacier) is the value shared by all cells in a given glacier
    ! The calling subroutine will need to map these values onto each grid cell.
    real(dp), dimension(nglacier), intent(inout) :: &
         powerlaw_c                     ! glacier-specific basal friction parameter (Pa (m/yr)^(-1/3))

    ! local variables

    integer :: ng

    real(dp), parameter :: &
         glacier_volume_timescale = 100.d0   ! timescale (yr)

    real(dp) :: &
         err_vol,                     & ! relative volume error, (V - V_target)/V_target
         term1, term2,                & ! terms in prognostic equation for powerlaw_c
         dpowerlaw_c                    ! change in powerlaw_c

    character(len=100) :: message

    ! The inversion works as follows:
    ! The change in C_p is proportional to the current value of C_p and to the relative error,
    !  err_vol = (V - V_target)/V_target.
    ! If err_vol > 0, we reduce C_p to make the glacier flow faster and thin.
    ! If err_vol < 0, we increase C_p to make the glacier flow slower and thicken.
    ! This is done with a characteristic timescale tau.
    ! We also include a term proportional to dV/dt so that ideally, C_p smoothly approaches
    !  the value needed to attain a steady-state V, without oscillating about the desired value.
    ! See the comments in module glissade_inversion, subroutine invert_basal_friction.
    ! Here is the prognostic equation:
    ! dC/dt = -C * (1/tau) * (V - V_target)/V_target + (2*tau/V_target) * dV/dt

    do ng = 1, nglacier

       if (volume_target(ng) > 0.0d0) then  ! this should be the case for most glaciers
          err_vol = (volume(ng) - volume_target(ng)) / volume_target(ng)
          term1 = -err_vol / glacier_volume_timescale
          term2 = -2.0d0 * dvolume_dt(ng) / volume_target(ng)
          dpowerlaw_c = powerlaw_c(ng) * (term1 + term2) * dt

          ! Limit to prevent a large relative change in one step
          if (abs(dpowerlaw_c) > 0.05d0 * powerlaw_c(ng)) then
             if (dpowerlaw_c > 0.0d0) then
                dpowerlaw_c =  0.05d0 * powerlaw_c(ng)
             else
                dpowerlaw_c = -0.05d0 * powerlaw_c(ng)
             endif
          endif

          ! Update powerlaw_c
          powerlaw_c(ng) = powerlaw_c(ng) + dpowerlaw_c

          ! Limit to a physically reasonable range
          powerlaw_c(ng) = min(powerlaw_c(ng), powerlaw_c_max)
          powerlaw_c(ng) = max(powerlaw_c(ng), powerlaw_c_min)

          if (verbose_glacier .and. main_task .and. ng == ngdiag) then
             print*, ' '
             print*, 'Invert for powerlaw_c: ngdiag =', ngdiag
             print*, 'V, V_target (km^3)', volume(ng)/1.0d9, volume_target(ng)/1.0d9
             print*, 'dV_dt (km^3/yr), relative err_vol:', dvolume_dt(ng)/1.0d9, err_vol
             print*, 'dt (yr), term1*dt, term2*dt:', dt, term1*dt, term2*dt
             print*, 'dpowerlaw_c, new powerlaw_c:', dpowerlaw_c, powerlaw_c(ng)
          endif

       else   ! volume_target(ng) = 0

          !TODO: Remove these glaciers from the inversion?
          ! For now, set C_p to the min value to minimize the thickness
          powerlaw_c(ng) = powerlaw_c_min

       endif

    enddo   ! ng

  end subroutine glacier_invert_powerlaw_c

!****************************************************

  subroutine glacier_area_volume(&
       ewn,           nsn,               &
       nglacier,      cism_glacier_id,   &
       cell_area,     thck,              &
       area,          volume,            &
       dthck_dt,      dvolume_dt)

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), intent(in) :: &
         cell_area                      ! grid cell area (m^2), assumed equal for all cells

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         thck                           ! ice thickness (m)

    real(dp), dimension(nglacier), intent(out) ::  &
         area,                & ! area of each glacier (m^2)
         volume                 ! volume of each glacier (m^3)

    real(dp), dimension(ewn,nsn), intent(in), optional ::  &
         dthck_dt               ! rate of change of ice thickness (m/yr)

    real(dp), dimension(nglacier), intent(out), optional ::  &
         dvolume_dt             ! rate of change of glacier volume (m^3/yr)

    ! local variables

    real(dp), dimension(:), allocatable :: &
         local_area, local_volume       ! area and volume on each processor, before global sum

    integer :: i, j, ng

    ! Initialize the output arrays
    area(:) = 0.0d0
    volume(:) = 0.0d0

    ! Allocate and initialize local arrays
    allocate(local_area(nglacier))
    allocate(local_volume(nglacier))
    local_area(:) = 0.0d0
    local_volume(:) = 0.0d0

    ! Compute the initial area and volume of each glacier.
    ! We need parallel sums, since a glacier can lie on two or more processors.

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'Compute glacier area and volume; cell_area (m^3) =', cell_area
    endif

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng >= 1) then
             local_area(ng) = local_area(ng) + cell_area
             local_volume(ng) = local_volume(ng) + cell_area * thck(i,j)
          endif
       enddo
    enddo

    area   = parallel_reduce_sum(local_area)
    volume = parallel_reduce_sum(local_volume)

    if (verbose_glacier .and. main_task) then
       print*, 'Max area (km^2)   =', maxval(area) * 1.0d-6    ! m^2 to km^2
       print*, 'Max volume (km^3) =', maxval(volume) * 1.0d-9  ! m^3 to km^3
       print*, ' '
       print*, 'Selected A (km^2) and  V(km^3) of large glaciers (> 3 km^3):'
       do ng = 1, nglacier
          if (volume(ng) * 1.0d-9 > 3.0d0) then  ! 3 km^3 or more
             write(6,'(i8,2f10.3)') ng, area(ng)*1.0d-6, volume(ng)*1.0d-9
          endif
       enddo
    endif

    ! Optionally, compute the rate of change of glacier volume
    if (present(dthck_dt) .and. present(dvolume_dt)) then
       ! use local_volume as a work array for dvolume_dt
       dvolume_dt(:) = 0.0d0
       local_volume(:) = 0.0d0
       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             ng = cism_glacier_id(i,j)
             if (ng >= 1) then
                local_volume(ng) = local_volume(ng) + cell_area * dthck_dt(i,j)
             endif
          enddo
       enddo
       dvolume_dt = parallel_reduce_sum(local_volume)
    endif

    deallocate(local_area)
    deallocate(local_volume)

  end subroutine glacier_area_volume

!****************************************************

  recursive subroutine quicksort(A, first, last)
 
    ! Given an unsorted integer array, return an array with elements sorted from low to high.
    ! Note: This is a template for a quicksort subroutine, but the subroutine actually called
    !       is glacier_quicksort below.

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

  end subroutine glacier_quicksort

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glissade_glacier

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
