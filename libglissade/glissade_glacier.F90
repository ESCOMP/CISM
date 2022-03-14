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
    use glimmer_paramets, only: thk0, len0, tim0, eps08
    use glimmer_physcon, only: scyr
    use glide_types
    use glimmer_log
    use cism_parallel, only: main_task, this_rank, nhalo

    implicit none

    private
    public :: verbose_glacier, glissade_glacier_init, glissade_glacier_smb, &
         glissade_glacier_advance_retreat, glissade_glacier_inversion

    logical, parameter :: verbose_glacier = .true.

    ! derived type that holds info for each glaciated grid cell
    type glacier_info
       integer :: id           ! input glacier ID, usually RGI
       integer :: indxi        ! i index of cell
       integer :: indxj        ! j index of cell
    end type glacier_info

    ! Glacier parameters used in this module
    ! Any of these could be added to the glacier derived type and set in the config file.
    ! Note: The constant, max and min values for powerlaw_c are in the basal_physics type.

    real(dp), parameter ::  &
         mu_star_const = 500.d0,            & ! uniform initial value for mu_star (mm/yr w.e/deg C)
         mu_star_min = 10.d0,               & ! min value of tunable mu_star (mm/yr w.e/deg C)
         mu_star_max = 1.0d5,               & ! max value of tunable mu_star (mm/yr w.e/deg C)
         glacier_mu_star_timescale = 1.d0,  & ! inversion timescale for mu_star (yr)
         glacier_powerlaw_c_timescale = 25.d0 ! inversion timescale for powerlaw_c (yr)

    integer, parameter :: &
         inversion_time_interval = 1          ! time interval (yr) between inversion calls; must be an integer

contains

!****************************************************      

  subroutine glissade_glacier_init(model, glacier)

    ! Initialize glaciers for an RGI region.
    ! If running with multiple disconnected glacier regions, call this subroutine once per region.
    ! Each region would be a separate instance.

    ! This subroutine creates an array called cism_glacier_id, which assigns an integer glacier ID
    !  to each CISM grid cell (i,j).  These IDs are numbered between 1 and nglacier,
    !  where nglacier is the total number of unique glacier IDs.
    ! This allows us to loop over IDs in the range (1:nglacier), which is more efficient than
    !  looping over the input RGI glacier IDs, which often have large gaps.
    ! Another array, cism_to_rgi_glacier_id, identifies the RGI ID associated with each CISM ID.
    ! The CISM input file contains the RGI IDs.

    use cism_parallel, only: distributed_gather_var, distributed_scatter_var, &
         parallel_reduce_sum, parallel_reduce_max, parallel_reduce_min, &
         broadcast, parallel_halo, parallel_globalindex

    type(glide_global_type),intent(inout) :: model

    type(glide_glacier) :: glacier    ! derived type for glacier info
                                      ! included in 'model', but passed separately to save typing

    ! local variables
    integer :: ewn, nsn               ! local grid dimensions
    integer :: global_ewn, global_nsn ! global grid dimensions
    integer :: itest, jtest, rtest    ! coordinates of diagnostic point
    real(dp) :: dew, dns              ! grid cell length in each direction (m)

    integer :: i, j, nc, ng, count
    integer :: iglobal, jglobal
    integer :: min_id, max_id
    character(len=100) :: message

    ! temporary global arrays
    integer, dimension(:,:), allocatable :: &
         rgi_glacier_id_global,     & ! global array of the input RGI glacier ID; maps (i,j) to RGI ID
         cism_glacier_id_global       ! global array of the CISM glacier ID; maps (i,j) to CISM glacier ID

    ! This type is declared at the top of the module
    type(glacier_info), dimension(:), allocatable :: & 
         glacier_list                 ! sorted list of glacier IDs with i and j indices

    integer :: &
         nglacier,                  & ! number of glaciers in global domain
         ncells_glacier,            & ! number of global grid cells occupied by glaciers at initialization
         current_id,                & ! current glacier_id from list
         gid_minval, gid_maxval       ! min and max values of glacier_id

    type(parallel_type) :: parallel   ! info for parallel communication

    !WHL - debug, for quicksort test
!    integer, dimension(:), allocatable :: test_list
!    integer ::  nlist
!    real(sp) :: random

    ! Set some local variables
    parallel = model%parallel

    global_ewn = parallel%global_ewn
    global_nsn = parallel%global_nsn
    ewn = model%general%ewn
    nsn = model%general%nsn
    dew = model%numerics%dew * len0   ! convert dew and dns to m
    dns = model%numerics%dns * len0
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_glacier_init'
       print*, ' '
       i = itest
       j = jtest
       print*, 'RGI glacier ID, rtest, itest, jtest:', rtest, itest, jtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') glacier%rgi_glacier_id(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    if (model%options%is_restart == RESTART_FALSE) then

       ! not a restart; initialize everything from the input file

       ! Note: At start-up, arrays in the glacier derived type are allocated with dimension(1),
       !  since nglacier has not yet been computed.
       ! Deallocate here, and reallocate below with dimension(nglacier).
       ! For a restart, nglacier is determined from the restart file,
       !  and these arrays should already have the correct dimensions.

       if (associated(glacier%glacierid)) deallocate(glacier%glacierid)
       if (associated(glacier%cism_to_rgi_glacier_id)) &
            deallocate(glacier%cism_to_rgi_glacier_id)
       if (associated(glacier%area)) deallocate(glacier%area)
       if (associated(glacier%volume)) deallocate(glacier%volume)
       if (associated(glacier%area_target)) deallocate(glacier%area_target)
       if (associated(glacier%volume_target)) deallocate(glacier%volume_target)
       if (associated(glacier%volume_in_init_region)) deallocate(glacier%volume_in_init_region)
       if (associated(glacier%dvolume_dt)) deallocate(glacier%dvolume_dt)
       if (associated(glacier%mu_star)) deallocate(glacier%mu_star)
       if (associated(glacier%powerlaw_c)) deallocate(glacier%powerlaw_c)

       ! Count the number of cells with glaciers
       ! Loop over locally owned cells

       count = 0
       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             if (glacier%rgi_glacier_id(i,j) > 0) then
                count = count + 1
             elseif (glacier%rgi_glacier_id(i,j) < 0) then  ! should not happen
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                write(message,*) 'RGI glacier_id < 0: i, j, value =', &
                     iglobal, jglobal, glacier%rgi_glacier_id(i,j)
                call write_log(message, GM_FATAL)
             endif
          enddo
       enddo

       ncells_glacier = parallel_reduce_sum(count)

       ! Gather the RGI glacier IDs to the main task
       if (main_task) allocate(rgi_glacier_id_global(global_ewn, global_nsn))
       call distributed_gather_var(glacier%rgi_glacier_id, rgi_glacier_id_global, parallel)

       ! Allocate a global array for the CISM glacier IDs on the main task.
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
               size(glacier%rgi_glacier_id,1), size(glacier%rgi_glacier_id,2)
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

          ! Deallocate the RGI global array (no longer needed after the glacier_list is built)
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

          glacier%nglacier = nglacier

          ! Fill two useful arrays:
          ! (1) The cism_glacier_id array maps each glaciated grid cell (i,j) to a CISM ID (between 1 and nglacier).
          ! (2) The cism_to_rgi_glacier_id array maps the CISM ID to the RGI glacier_id.
          ! By carrying i and j in the sorted glacier_list, we can efficiently fill cism_glacier_id.
          ! Note: cism_to_rgi_glacier_id cannot be allocated until nglacier is known.

          allocate(glacier%cism_to_rgi_glacier_id(nglacier))
          glacier%cism_to_rgi_glacier_id(:) = 0

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
                glacier%cism_to_rgi_glacier_id(ng) = glacier_list(nc)%id
             endif
             i = glacier_list(nc)%indxi
             j = glacier_list(nc)%indxj
             cism_glacier_id_global(i,j) = ng
             if (ng == nglacier/2) then   ! random glacier
                print*, nc, i, j, cism_glacier_id_global(i,j), glacier%cism_to_rgi_glacier_id(ng)
             endif
             if (ng > nglacier) then
                write(message,*) 'CISM glacier ID > nglacier, i, j , ng =', i, j, ng
                call write_log(message, GM_FATAL)
             endif
          enddo

          deallocate(glacier_list)

          if (verbose_glacier) then
             print*, 'maxval(cism_to_rgi_glacier_id) =', maxval(glacier%cism_to_rgi_glacier_id)
             print*, 'maxval(cism_glacier_id_global) =', maxval(cism_glacier_id_global)
          endif

       endif   ! main_task

       ! Scatter cism_glacier_id_global to all processors
       ! Note: This global array is deallocated in the distributed_scatter_var subroutine
       call distributed_scatter_var(glacier%cism_glacier_id, cism_glacier_id_global, parallel)

       ! Copy cism_glacier_id to cism_glacier_id_init, which is saved and used for mu_star inversion
       glacier%cism_glacier_id_init(:,:) = glacier%cism_glacier_id(:,:)

       ! Broadcast nglacier and cism_to_rgi_glacier_id from the main task to all processors
       call broadcast(glacier%nglacier)
       nglacier = glacier%nglacier

       if (.not.associated(glacier%cism_to_rgi_glacier_id)) &
            allocate(glacier%cism_to_rgi_glacier_id(nglacier))
       call broadcast(glacier%cism_to_rgi_glacier_id)

       ! Allocate glacier arrays with dimension(nglacier)

       allocate(glacier%glacierid(nglacier))
       allocate(glacier%area(nglacier))
       allocate(glacier%area_target(nglacier))
       allocate(glacier%volume(nglacier))
       allocate(glacier%volume_target(nglacier))
       allocate(glacier%volume_in_init_region(nglacier))
       allocate(glacier%dvolume_dt(nglacier))
       allocate(glacier%mu_star(nglacier))
       allocate(glacier%powerlaw_c(nglacier))

       ! Compute the initial area and volume of each glacier.
       ! The initial values are targets for inversion of mu_star and powerlaw_c.

       call glacier_area_volume(&
            ewn,           nsn,             &
            nglacier,                       &
            glacier%cism_glacier_id,        &
            dew*dns,                        &
            model%geometry%thck*thk0,       &
            glacier%area,                   &
            glacier%volume)

       ! Initialize other glacier arrays
       glacier%area_target(:) = glacier%area(:)
       glacier%volume_target(:) = glacier%volume(:)
       glacier%volume_in_init_region(:) = glacier%volume(:)
       glacier%dvolume_dt(:) = 0.0d0
       glacier%mu_star(:) = mu_star_const
       glacier%powerlaw_c(:) = model%basal_physics%powerlaw_c_const

       ! Check for area_target = 0 and volume_target = 0.
       ! In practice, volume_target = 0 might not be problematic;
       !  we would just lower powerlaw_c to obtain a thin glacier.
       if (main_task) then
          do ng = 1, nglacier
             if (glacier%area_target(ng) == 0.0d0) then
                write(message,*) 'Glacier area target = 0: ng =', ng
                call write_log(message, GM_FATAL)
             endif
             if (glacier%volume_target(ng) == 0.0d0) then
                write(message,*) 'Glacier volume target = 0: ng, area (km^2) =', &
                     ng, glacier%area(ng)/1.0d6
                call write_log(message)
             endif
          enddo   ! ng
       endif

    else  ! restart

       ! In this case, most glacier info has already been read from the restart file.
       ! From the restart file, nglacier is found as the length of dimension 'glacierid'.
       ! The 1D glacier arrays are then allocated with dimension(nglacier) in subroutine glide_allocarr.
       ! The following glacier arrays should be present in the restart file:
       !     rgi_glacier_id, cism_glacier_id, cism_to_rgi_glacier_id, mu_star, powerlaw_c
       ! If inverting for mu_star and powerlaw_c, the restart file will also include these arrays:
       !     area_target, volume_target, cism_glacier_id_init
       ! (Although area_target is not strictly needed for inversion, it is included as a diagnostic.)
       ! These remaining parameters are set here:
       !     glacierid, ngdiag

       nglacier = glacier%nglacier

       ! Check that the glacier arrays which are read from the restart file have nonzero values.
       ! Note: These arrays are read on all processors.

       if (maxval(glacier%mu_star) <= 0.0d0) then
          call write_log ('Error, no positive values for glacier_mu_star', GM_FATAL)
       endif

       if (maxval(glacier%powerlaw_c) <= 0.0d0) then
          call write_log ('Error, no positive values for glacier_powerlaw_c', GM_FATAL)
       endif

       max_id = maxval(glacier%cism_glacier_id_init)
       max_id = parallel_reduce_max(max_id)
       if (max_id <= 0) then
          call write_log ('Error, no positive values for cism_glacier_id_init', GM_FATAL)
       endif

       min_id = minval(glacier%cism_to_rgi_glacier_id)
       min_id = parallel_reduce_min(min_id)
       if (min_id <= 0) then
          write(message,*) 'Error, minval(cism_to_rgi_glacier_id) =', min_id
          call write_log(message, GM_FATAL)
       endif

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          if (maxval(glacier%volume_target) <= 0.0d0) then
             call write_log ('Error, no positive values for glacier_volume_target', GM_FATAL)
          endif
       endif

       ! Compute the initial area and volume of each glacier.
       ! This is not strictly necessary for a restart, but is included as a diagnostic.

       call glacier_area_volume(&
            ewn,           nsn,             &
            nglacier,                       &
            glacier%cism_glacier_id,        &
            dew*dns,                        &
            model%geometry%thck*thk0,       &
            glacier%area,                   &
            glacier%volume,                 &
            glacier%cism_glacier_id_init,   &
            glacier%volume_in_init_region)

    endif   ! not a restart

    ! The remaining code applies to both start-up and restart runs

    ! Allocate and fill the glacierid dimension array
    do ng = 1, nglacier
       glacier%glacierid(ng) = ng
    enddo

    ! Given powerlaw_c for each glacier, compute model%basal_physics%powerlaw_c,
    !  a 2D array defined at cell vertices.
    ! Set model%basal_physics%powerlaw_c = 0 at vertices that are not adjacent
    !  to any glacier cells.

    call glacier_powerlaw_c_to_2d(&
         ewn,             nsn,                 &
         nglacier,                             &
         glacier%cism_glacier_id,              &
         glacier%powerlaw_c,                   &
         model%basal_physics%powerlaw_c,       &
         parallel)

    ! Halo updates for the 2D glacier_id arrays
    call parallel_halo(glacier%rgi_glacier_id, parallel)
    call parallel_halo(glacier%cism_glacier_id, parallel)
    call parallel_halo(glacier%cism_glacier_id_init, parallel)

    ! Set the minimum thickness (m) for ice to be counted as a glacier.
    ! Choose this limit equal to the dynamics threshold (actually, slightly
    !  less in case of roundoff error).
    ! Thus, any ice that is not part of a glacier is dynamically inactive,
    !  but could receive a glacier ID and become active with thickening.

    glacier%minthck = model%numerics%thklim*thk0 - eps08

    ! Set the index of the diagnostic glacier, using the CISM glacier ID for the diagnostic point
    if (this_rank == rtest) then
       glacier%ngdiag = glacier%cism_glacier_id(itest,jtest)
    endif
    call broadcast(glacier%ngdiag, rtest)

    ! Write some values for the diagnostic glacier
    if (verbose_glacier .and. main_task) then
       print*, ' '
       ng = glacier%ngdiag
       print*, 'Glacier ID for diagnostic cell: r, i, j, ng =', rtest, itest, jtest, ng
       print*, 'area target (km^2) =', glacier%area_target(ng) / 1.0d6
       print*, 'volume target (km^3) =', glacier%volume_target(ng) / 1.0d9
       print*, 'mu_star (mm/yr w.e./deg) =', glacier%mu_star(ng)
       print*, 'powerlaw_c (Pa (m/yr)^(-1/3)) =', glacier%powerlaw_c(ng)
       print*, 'Done in glissade_glacier_init'
    endif

  end subroutine glissade_glacier_init

!****************************************************

  subroutine glissade_glacier_smb(&
       ewn,      nsn,                  &
       itest,    jtest,  rtest,        &
       nglacier,                       &
       cism_glacier_id,                &
       t_mlt,                          &
       snow,             artm,         &
       mu_star,                        &
       glacier_smb)

    ! Compute the SMB in each grid cell using an empirical relationship
    !  based on Maussion et al. (2019):
    !
    !     SMB = snow - mu_star * max(artm - T_mlt, 0),
    !
    ! where snow = monthly mean snowfall rate (mm/yr w.e.),
    !       mu_star is a glacier-specific tuning parameter (mm/yr w.e./deg C),
    !       atrm = monthly mean air temperature (deg C),
    !       T_mlt = monthly mean air temp above which ablation occurs (deg C)
    !
    ! This subroutine should be called at least once per model month.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier,                    & ! total number of glaciers in the domain
         itest, jtest, rtest            ! coordinates of diagnostic point

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), intent(in) :: &
         t_mlt                          ! min temperature (deg C) at which ablation occurs

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow,                        & ! monthly mean snowfall rate (mm w.e./yr)
         artm                           ! artm adjusted for elevation using t_lapse (deg C)

    real(dp), dimension(nglacier), intent(in) :: &
         mu_star                        ! glacier-specific SMB tuning parameter (mm w.e./yr/deg)

                                        ! defined as positive for T decreasing with height

    real(dp), dimension(ewn,nsn), intent(out) :: &
         glacier_smb                    ! SMB in each gridcell (mm w.e./yr)

    ! local variables

    integer :: i, j, ng

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_glacier_smb'
       print*, 'minval, maxval(snow) =', minval(snow), maxval(snow)
       print*, 'minval, maxval(artm) =', minval(artm), maxval(artm)
       print*, 't_mlt (deg C) =', t_mlt
    endif

    ! initialize
    glacier_smb(:,:) = 0.0d0

    ! compute SMB

    do j = 1, nsn
       do i = 1, ewn
          ng = cism_glacier_id(i,j)

          if (ng > 0) then
             glacier_smb(i,j) = snow(i,j) - mu_star(ng) * max(artm(i,j) - t_mlt, 0.0d0)
          endif

          if (verbose_glacier .and. this_rank == rtest .and. i == itest .and. j == jtest) then
             print*, ' '
             print*, 'Glacier SMB calculation: rank i, j, mu_star =', &
                  this_rank, i, j, mu_star(ng)
             print*, '   snow (mm/yr w.e.), artm (C), SMB (mm/yr w.e.) =', &
                  snow(i,j), artm(i,j), glacier_smb(i,j)
          endif

       enddo
    enddo

  end subroutine glissade_glacier_smb

!****************************************************

  subroutine glissade_glacier_advance_retreat(&
       ewn,             nsn,             &
       itest,   jtest,  rtest,           &
       usrf,            thck,            &
       acab_applied,    dt,              &
       glacier_minthck,                  &
       cism_glacier_id_init,             &
       cism_glacier_id,                  &
       parallel)

    ! Allow glaciers to advance and retreat.
    ! This subroutine should be called after the transport/SMB calculation.
    !
    ! The rules are as follows:
    ! * At start-up, glaciated cells have cism_glacier_id in the range (1, nglacier).
    !   Other cells have cism_glacier_id = 0.
    !   The initial cism_glacier_id array is saved as cism_glacier_id_init.
    ! * If a cell has H <= minthck and cism_glacier_id > 0, we set cism_glacier_id = 0.
    !   It no longer contributes to glacier area or volume.
    !   Here, minthck is a threshold for counting ice as part of a glacier.
    !   By default, minthck = model%numerics%thklim, typically 1 m.
    !   (Actually minthck is slightly less than thklim, to make sure these cells
    !   are not dynamically active.)
    ! * When a cell has H > minthck and cism_glacier_id = 0, we give it a nonzero ID:
    !   either (1) cism_glacier_id_init, if the initial ID > 0,
    !   or (2) the ID of an adjacent glaciated neighbor (the neighbor with
    !   the highest surface elevation, if there is more than one).
    !   Preference is given to (1), to preserve the original glacier outlines
    !   as much as possible.
    ! * If H > minthck in a cell with cism_glacier_id_init = 0 and no glaciated neighbors,
    !   we do not give it a glacier ID.  Instead, we set H = minthck and remove the excess ice.
    !   This ice remains dynamically inactive.
    !   Thus, there is no glacier inception; we only allow existing glaciers to advance.

    use cism_parallel, only: parallel_globalindex, parallel_halo

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest            ! coordinates of diagnostic cell

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         usrf                           ! upper surface elevation (m)

    real(dp), dimension(ewn,nsn), intent(inout) ::  &
         thck,                        & ! ice thickness (m)
         acab_applied                   ! SMB applied to ice surface (m/s)

    real(dp), intent(in) :: &
         dt,                          & ! time step (s)
         glacier_minthck                ! min ice thickness (m) counted as part of a glacier

    integer, dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init           ! cism_glacier_id at the start of the run

    integer, dimension(ewn,nsn), intent(inout) :: &
         cism_glacier_id                ! current cism glacier_id, > 0 for glaciated cells

    type(parallel_type), intent(in) :: parallel  !WHL - diagnostic only

    ! local variables

    real(dp), dimension(ewn,nsn) :: &
         cism_glacier_id_old            ! old value of cism_glacier_id

    real(dp) :: &
         usrf_max,             & ! highest elevation (m) in a neighbor cell
         dthck                   ! ice thickness loss (m)

    integer :: i, j, ii, jj, ip, jp
    integer :: iglobal, jglobal
    integer :: ng


    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_glacier_advance_retreat'
    endif

    ! Check for retreat: cells with cism_glacier_id > 0 but H = 0

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0 .and. thck(i,j) <= glacier_minthck) then
             !WHL - debug
             if (verbose_glacier .and. this_rank==rtest) then
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                print*, 'Set ID = 0: ig, jg, old ID, thck =', &
                     iglobal, jglobal, ng, thck(i,j)
             endif
             cism_glacier_id(i,j) = 0
          endif
       enddo
    enddo

    ! Check for retreat: cells with cism_glacier_id = 0 but H > H_min

    ! Save a copy of the old cism_glacier_id.
    ! This is to prevent the algorithm from depending on the loop direction.
    cism_glacier_id_old(:,:) = cism_glacier_id(:,:)

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng == 0 .and. thck(i,j) > glacier_minthck) then
             ! Assign this cell its original ID, if > 0
             if (cism_glacier_id_init(i,j) > 0) then
                cism_glacier_id(i,j) = cism_glacier_id_init(i,j)
                !WHL - debug
                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, 'Set ID = init ID: ig, jg, new ID, thck =',&
                        iglobal, jglobal, cism_glacier_id(i,j), thck(i,j)
                endif
             else  ! assign the ID of an adjacent ice-covered cell, if possible
                usrf_max = 0.0d0
                do jj = -1, 1
                   do ii = -1, 1
                      if (ii /= 0 .and. jj /= 0) then  ! one of 8 neighbors
                         ip = i + ii
                         jp = j + jj
                         if (cism_glacier_id_old(ip,jp) > 0 .and. &
                              thck(ip,jp) > glacier_minthck) then
                            if (usrf(ip,jp) > usrf_max) then
                               usrf_max = usrf(ip,jp)
                               cism_glacier_id(i,j) = cism_glacier_id(ip,jp)
                               !WHL - debug
                               if (verbose_glacier .and. this_rank == rtest) then
                                  call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                                  print*, 'Set ID = neighbor ID, ig, jg, new ID, thck =', &
                                       iglobal, jglobal, cism_glacier_id(i,j), thck(i,j)
                               endif
                            endif
                         endif
                      endif
                   enddo  ! ii
                enddo   ! jj
             endif   ! cism_glacier_id_init > 0

             ! If the cell still has cism_glacier_id = 0 and H > glacier_minthck,
             !  then cap the thickness at glacier_minthck.
             ! Note: The ice removed is used to increment acab_applied, the ice SMB in m/s.
             !       Thus, the total SMB flux will generally be more negative during time steps
             !        when this subroutine is solved.
             if (cism_glacier_id(i,j) == 0 .and. thck(i,j) > glacier_minthck) then
                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, 'Cap H = glacier_minthck, ig, jg, thck =', &
                        iglobal, jglobal, thck(i,j)
                endif
                dthck = thck(i,j) - glacier_minthck
                thck(i,j) = glacier_minthck
                acab_applied(i,j) = acab_applied(i,j) - dthck/dt   ! m/s
             endif

          endif   ! ng = 0, H > 0
       enddo   ! i
    enddo   ! j

    ! Halo updates for output arrays
    call parallel_halo(cism_glacier_id, parallel)
    call parallel_halo(thck, parallel)

  end subroutine glissade_glacier_advance_retreat

!****************************************************

  subroutine glissade_glacier_inversion(model, glacier)

    use glissade_grid_operators, only: glissade_stagger
    use cism_parallel, only: parallel_reduce_sum, staggered_parallel_halo

    ! input/output arguments

    type(glide_global_type), intent(inout) :: model

    type(glide_glacier) :: glacier  ! derived type for glacier info
                                    ! included in 'model', but passed separately to save typing

    ! local variables

    integer ::  &
         itest, jtest, rtest,     & ! coordinates of diagnostic cell
         ewn, nsn                   ! number of cells in each horizontal direction

    real(dp) :: &
         dt,                      & ! time step (yr)
         dew, dns                   ! grid cell dimensions (m)

    integer :: nglacier       ! number of glaciers
    integer :: ngdiag         ! CISM index of diagnostic glacier
    integer :: i, j, ng

    integer, dimension(model%general%ewn, model%general%nsn) ::  &
         ice_mask                   ! = 1 where ice is present (thck > thklim), else = 0


    integer, dimension(model%general%ewn, model%general%nsn) ::  &
         glacier_mask               ! = 1 where glacier ice is present (thck > thklim), else = 0

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck,                    & ! ice thickness (m)
         dthck_dt,                & ! rate of change of thickness (m/yr)
         powerlaw_c_icegrid         ! powerlaw_c on the unstaggered ice grid

    type(parallel_type) :: parallel ! info for parallel communication

    real(dp), save ::  &            ! time since the last averaging computation;
         time_since_last_avg = 0.0d0  ! set to 1 yr for now

    real(dp) ::  smb_annmean     ! annual mean SMB for a given cell

    real(dp), dimension(glacier%nglacier) :: &
         smb_init_area,        & ! SMB over initial area determined by cism_glacier_id_init
         smb_current_area        ! SMB over cufrent area determined by cism_glacier_id

    ! Note: The glacier type includes the following:
    ! integer ::  nglacier          ! number of glaciers in the global domain
    ! integer ::  ngdiag            ! CISM index of diagnostic glacier
    ! real(dp), dimension(:) :: area              ! glacier area (m^2)
    ! real(dp), dimension(:) :: area_target       ! glacier area target (m^2)
    ! real(dp), dimension(:) :: volume            ! glacier volume (m^3)
    ! real(dp), dimension(:) :: volume_target     ! glacier volume target (m^3)
    ! real(dp), dimension(:) :: dvolume_dt        ! rate of change of glacier volume (m^3/yr)
    ! real(dp), dimension(:) :: mu_star           ! SMB parameter for each glacier (mm/yr w.e./deg K)
    ! real(dp), dimension(:) :: powerlaw_c        ! basal friction parameter for each glacier (Pa (m/yr)^(-1/3))
    ! integer, dimension(:,:) :: cism_glacier_id  ! CISM glacier ID for each grid cell
    ! integer, dimension(:,:) :: cism_glacier_id_init  ! initial value of CISM glacier ID
    ! real(dp), dimension(:,:) :: snow_accum      ! snow accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: Tpos_accum      ! max(artm-T_mlt,0) accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: dthck_dt_accum  ! dthck_dt accumulated and averaged over 1 year

    ! Set some local variables

    parallel = model%parallel

    ewn = model%general%ewn
    nsn = model%general%nsn
    dew = model%numerics%dew * len0         ! convert to m
    dns = model%numerics%dns * len0         ! convert to m
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_glacier_inversion, diag cell (r, i, j) =', rtest, itest, jtest
    endif

    nglacier = glacier%nglacier
    ngdiag = glacier%ngdiag

    ! some unit conversions
    dt = model%numerics%dt * tim0/scyr          ! model units to yr
    thck = model%geometry%thck * thk0           ! model units to m
    dthck_dt = model%geometry%dthck_dt * scyr   ! m/s to m/yr

    ! Accumulate the 2D fields used for inversion: snow, Tpos and dthck_dt.

    call accumulate_glacier_fields(&
         ewn,                    nsn,                       &
         dt,                     time_since_last_avg,       &
         model%climate%snow,     glacier%snow_accum,        &  ! mm/yr w.e.
         max(model%climate%artm - glacier%t_mlt, 0.0d0),    &
         glacier%Tpos_accum,                                &  ! deg C
         dthck_dt,               glacier%dthck_dt_accum)       ! m/yr ice

    if (verbose_glacier .and. this_rank == rtest) then
       i = itest; j = jtest
       print*, 'r, i, j, time, time_since_last_avg, snow, Tpos, dthck_dt:', &
            this_rank, i, j, model%numerics%time, time_since_last_avg, &
            glacier%snow_accum(i,j), glacier%Tpos_accum(i,j), glacier%dthck_dt_accum(i,j)
    endif

    ! Check whether it is time to do the inversion.
    ! Note: model%numerics%time has units of yr.

    if (abs(time_since_last_avg - real(inversion_time_interval,dp)) < eps08) then

       if (verbose_glacier .and. this_rank == rtest) then
          print*, 'calculate_glacier_averages, time_since_last_avg =', time_since_last_avg
       endif

       ! compute annual average of glacier fields

       call calculate_glacier_averages(&
            ewn,                    nsn,   &
            time_since_last_avg,           &  ! yr
            glacier%snow_accum,            &  ! mm/yr w.e.
            glacier%Tpos_accum,            &  ! deg C
            glacier%dthck_dt_accum)           ! m/yr ice

       if (verbose_glacier .and. this_rank == rtest) then
          i = itest; j = jtest
          print*, 'Annual glacier averages, r, i, j:', rtest, itest, jtest
          print*, '   snow (mm/yr w.e.)=', glacier%snow_accum(i,j)
          print*, '   Tpos (deg C)    =', glacier%Tpos_accum(i,j)
          print*, '   dthck_dt (m/yr) =', glacier%dthck_dt_accum(i,j)
       endif

       ! Compute the current area and volume of each glacier
       ! Note: This requires global sums. For now, do the computation independently on each task.
       ! The difference between volume and volume_target is used to invert for powerlaw_c.
       ! The area is not used for inversion but is computed as a diagnostic.

       call glacier_area_volume(&
            ewn,           nsn,              &
            nglacier,                        &
            glacier%cism_glacier_id,         &
            dew*dns,                         &  ! m^2
            model%geometry%thck * thk0,      &  ! m
            glacier%area,                    &  ! m^2
            glacier%volume,                  &  ! m^3
            glacier%cism_glacier_id_init,    &
            glacier%volume_in_init_region,   &  ! m^3
            glacier%dthck_dt_accum,          &  ! m/yr
            glacier%dvolume_dt)                 ! m^3/yr

       if (verbose_glacier .and. main_task) then
          print*, ' '
          print*, 'Update area (km^2) and volume (km^3) for glacier:', ngdiag
          print*, 'Current area and volume:', glacier%area(ngdiag)/1.0d6, &
               glacier%volume(ngdiag)/1.0d9
          print*, '   Volume in init region =', glacier%volume_in_init_region(ngdiag)/1.0d9
          print*, '   Target area and volume:', glacier%area_target(ngdiag)/1.0d6, &
               glacier%volume_target(ngdiag)/1.0d9
          print*, '   dV_dt (m^3/yr):', glacier%dvolume_dt(ngdiag)/1.0d9
          print*, ' '
          print*, 'All glaciers: ng, A, A_target, Aerr, V, V_target, Verr:'
          do ng = 1, nglacier
             write(6,'(i6,3f12.4,3f14.6)') ng, glacier%area(ng)/1.0d6, glacier%area_target(ng)/1.0d6, &
                  (glacier%area(ng) - glacier%area_target(ng))/1.0d6, &
                  glacier%volume_in_init_region(ng)/1.0d9, glacier%volume_target(ng)/1.0d9, &
                  (glacier%volume_in_init_region(ng) - glacier%volume_target(ng))/1.0d9
          enddo
       endif

       ! Given the current and target glacier areas, invert for mu_star

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION) then

          call glacier_invert_mu_star(&
               ewn,                 nsn,                   &
               nglacier,            ngdiag,                &
               glacier%snow_accum,  glacier%Tpos_accum,    &
               glacier%cism_glacier_id_init,               &
               glacier%mu_star)

          smb_init_area(:) = 0.0d0
          smb_current_area(:) = 0.0d0

          !WHL - debug - compute the SMB over the original and current glacier area
          do j = nhalo+1, nsn-nhalo
             do i = nhalo+1, ewn-nhalo

                ! increment SMB over initial glacier area
                ng = glacier%cism_glacier_id_init(i,j)
                if (ng > 0) then
                   smb_annmean = glacier%snow_accum(i,j) - glacier%mu_star(ng) * glacier%Tpos_accum(i,j)
                   smb_init_area(ng) = smb_init_area(ng) + smb_annmean
                endif

                ! increment SMB over current glacier area
                ng = glacier%cism_glacier_id(i,j)
                if (ng > 0) then
                   smb_annmean = glacier%snow_accum(i,j) - glacier%mu_star(ng) * glacier%Tpos_accum(i,j)
                   smb_current_area(ng) = smb_current_area(ng) + smb_annmean
                endif

             enddo
          enddo

          ! global sums
          smb_init_area = parallel_reduce_sum(smb_init_area)
          smb_current_area = parallel_reduce_sum(smb_current_area)

          ! take area average
          where (glacier%area_target > 0.0d0) &
               smb_init_area(:) = smb_init_area(:) / glacier%area_target(:)

          where (glacier%area > 0.0d0) &
               smb_current_area(:) = smb_current_area(:) / glacier%area(:)

!          if (verbose_glacier .and. main_task) then
!             print*, ' '
!             print*, 'All glaciers: smb_init_area, smb_current_area'
!             do ng = 1, nglacier
!                write(6,'(i6,2f12.4)') ng, smb_init_area(ng), smb_current_area(ng)
!             enddo
!          endif

       endif   ! invert for mu_star

       ! Given the current and target glacier volumes, invert for powerlaw_c
       ! Note: The current volume is computed not over the entire glacier
       !       (which could be advanced or retreat compared to the initial extent),
       !       but over the initial region defined by cism_glacier_id_init.
       !       This prevents the inversion scheme from generating thickness errors
       !       to compensate for area errors.

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then

          call glacier_invert_powerlaw_c(&
               ewn,                 nsn,              &
               nglacier,            ngdiag,           &
               model%basal_physics%powerlaw_c_min,    &
               model%basal_physics%powerlaw_c_max,    &
               glacier%volume_in_init_region,         &
               glacier%volume_target,                 &
               glacier%dvolume_dt,                    &
               glacier%powerlaw_c)

       endif

       !WHL - debug
       if (verbose_glacier .and. main_task) then
!          print*, ' '
!          print*, 'All glaciers: mu_star, powerlaw_c'
!          do ng = 1, nglacier
!             write(6,*) ng, glacier%mu_star(ng), glacier%powerlaw_c(ng)
!          enddo
       endif

       ! Given powerlaw_c for each glacier, compute a 2D array of powerlaw_c,
       !  part of the basal_physics derived type.
       ! Set basal_physics%powerlaw_c = 0 at vertices that are not adjacent
       !  to any glacier cells.

       call glacier_powerlaw_c_to_2d(&
            ewn,             nsn,                 &
            nglacier,                             &
            glacier%cism_glacier_id,              &
            glacier%powerlaw_c,                   &
            model%basal_physics%powerlaw_c,       &
            parallel)

       ! Reset the accumulated fields
       call reset_glacier_fields(&
            ewn,           nsn,    &
            glacier%snow_accum,    &
            glacier%Tpos_accum,    &
            glacier%dthck_dt_accum)

    endif   ! time to do inversion

  end subroutine glissade_glacier_inversion

!****************************************************

  subroutine glacier_invert_mu_star(&
       ewn,              nsn,               &
       nglacier,         ngdiag,            &
       snow_accum,       Tpos_accum,        &
       cism_glacier_id_init,                &
       mu_star)

    ! Given the current glacier areas and area targets,
    ! invert for the parameter mu_star in the glacier SMB formula

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                 & ! number of cells in each horizontal direction
         nglacier,                 & ! total number of glaciers in the domain
         ngdiag                      ! CISM ID of diagnostic glacier

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow_accum,                  & ! time-avg snowfall for each cell (mm/yr w.e.)
         Tpos_accum                     ! time-avg of max(artm - T_mlt) for each cell (deg)

    integer, dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init           ! cism_glacier_id at the start of the run

    ! Note: Here, mu_star_glacier(nglacier) is the value shared by all cells in a given glacier
    ! The calling subroutine will need to map these values onto each grid cell.
    real(dp), dimension(nglacier), intent(inout) :: &
         mu_star                        ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)

    ! local variables
    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         glacier_snow, glacier_Tpos,  & ! global sums for each glacier
         mu_star_new                    ! new target value of mu_star, toward which we relax

    character(len=100) :: message

    ! Inversion for mu_star is more direct than inversion for powerlaw_c.
    ! Instead of solving a damped harmonic oscillator equation for mu_star,
    !  we compute mu_star for each glacier such that SMB = 0 over the initial extent.
    !
    ! The SMB for glacier ng is given by
    !      sum_ij(smb) = sum_ij(snow) - mu_star(ng) * sum_ij(Tpos),
    ! where Tpos = max(artm - T_mlt, 0),
    ! and sum_ij notes a sum over all cells (i,j) in the glacier.
    !
    ! Setting SMB = 0 and rearranging, we get
    !      mu_star(ng) = sum_ij(snow) / sum_ij(Tpos)
    !
    ! Thus, given the annual average of snow and Tpos for each grid cell in a glacier,
    ! we can find mu_star such that SMB = 0.
    !
    ! We take sums are taken over the target area of each glacier, using cism_glacier_id_init.
    ! If a glacier is too large, the net SMB will be < 0 and the glacier should shrink.
    ! Similarly, if the glacier is too small, the net SMB > 0 and the glacier should grow.
    !
    ! Optionally, by setting glacier_mu_star_timescale > inversion_time_interval,
    ! we can relax toward the computed mu_star instead of going there immediately.
    !
    ! Notes:
    !
    ! (1) This approach works only for land-based glaciers.
    !     TODO: Modify for marine-terminating glaciers.
    ! (2) If spinning up with climatological SMB, then mu_star will have the same value
    !     throughout the inversion.  This means that when the glacier advances or retreats,
    !     mu_star will not change to compensate.
    ! (3) If the glacier advances, then its net SMB should be < 0, so it should lose mass.
    !     It is possible that the steady-state glacier will have the correct total volume,
    !     but will be too advanced and too thin. An alternative is to adjust C_p
    !     based on the volume contained within the original glacier outline.
    !     TODO: Try this.  Get the volume right within the original outlines,
    !     which allows a slight advance (e.g., if the ice is too thin in the center
    !     and thick at the margins) but hopefully not far beyond those outlines.

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'In glissade_invert_mu_star'
    endif

    glacier_snow(:) = 0.0d0
    glacier_Tpos(:) = 0.0d0

    ! Compute local sums over the initial extent of each glacier
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id_init(i,j)
          if (ng > 0) then
             glacier_snow(ng) = glacier_snow(ng) + snow_accum(i,j)
             glacier_Tpos(ng) = glacier_Tpos(ng) + Tpos_accum(i,j)
          endif
       enddo
    enddo

    ! Compute global sums
    glacier_snow = parallel_reduce_sum(glacier_snow)
    glacier_Tpos = parallel_reduce_sum(glacier_Tpos)

    ! For each glacier, compute the new mu_star

    do ng = 1, nglacier

       if (glacier_Tpos(ng) > 0.0d0) then  ! ablation is nonzero

          ! Compute the value of mu_star that will give SMB = 0 over the target area
          mu_star_new(ng) = glacier_snow(ng) / glacier_Tpos(ng)

          ! Limit to a physically reasonable range
          mu_star_new(ng) = min(mu_star_new(ng), mu_star_max)
          mu_star_new(ng) = max(mu_star_new(ng), mu_star_min)

          if (verbose_glacier .and. main_task .and. ng == ngdiag) then
             print*, ' '
             print*, 'ng, sum_snow, sum_Tpos:', ng, glacier_snow(ng), glacier_Tpos(ng)
             print*, 'Old and new mu_star:', mu_star(ng), mu_star_new(ng)
          endif

          ! Relax toward the new value
          ! By default, inversion_time_interval = glacier_mu_star_timescale = 1 yr
          mu_star(ng) = mu_star(ng) + (mu_star_new(ng) - mu_star(ng)) &
               * max(inversion_time_interval/glacier_mu_star_timescale, 1.0d0)

       else   ! glacier_Tpos = 0; no ablation

          mu_star(ng) = mu_star_max

          if (verbose_glacier .and. main_task) then
             print*, 'Warning: no ablation for glacier', ng
          endif

       endif

    enddo   ! ng

  end subroutine glacier_invert_mu_star

!****************************************************

  subroutine glacier_invert_powerlaw_c(&
       ewn,              nsn,               &
       nglacier,         ngdiag,            &
       powerlaw_c_min,   powerlaw_c_max,    &
       volume,           volume_target,     &
       dvolume_dt,       powerlaw_c)

    use glimmer_physcon, only: scyr

    ! Given the current glacier volumes and volume targets,
    ! invert for the parameter powerlaw_c in the relationship for basal sliding.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier,                    & ! total number of glaciers in the domain
         ngdiag                         ! ID of diagnostic glacier

    real(dp), intent(in) :: &
         powerlaw_c_min, powerlaw_c_max ! min and max allowed values of powerlaw_c (Pa (m/yr)^(-1/3))

    real(dp), dimension(nglacier), intent(in) :: &
         volume,                      & ! current glacier volume over the target region (m^3)
         volume_target,               & ! volume target (m^3)
         dvolume_dt                     ! rate of change of volume (m^3/yr)

    real(dp), dimension(nglacier), intent(inout) :: &
         powerlaw_c                     ! glacier-specific basal friction parameter (Pa (m/yr)^(-1/3))

    ! local variables

    integer :: ng

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
    ! dC/dt = -C * (1/tau) * [(V - V_target)/V_target + (2*tau/V_target) * dV/dt]

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'In glissade_invert_powerlaw_c'
    endif

    do ng = 1, nglacier

       if (volume_target(ng) > 0.0d0) then  ! this should be the case for nearly all glaciers
          err_vol = (volume(ng) - volume_target(ng)) / volume_target(ng)
          term1 = -err_vol / glacier_powerlaw_c_timescale
          term2 = -2.0d0 * dvolume_dt(ng) / volume_target(ng)
          dpowerlaw_c = powerlaw_c(ng) * (term1 + term2) * inversion_time_interval

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
             print*, 'dt (yr), term1*dt, term2*dt:', inversion_time_interval, &
                  term1*inversion_time_interval, term2*inversion_time_interval
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

  subroutine glacier_powerlaw_c_to_2d(&
       ewn,            nsn,             &
       nglacier,                        &
       cism_glacier_id,                 &
       glacier_powerlaw_c,              &
       basal_physics_powerlaw_c,        &
       parallel)

    ! Given model%glacier%powerlaw_c(ng) for each glacier,
    ! compute basal_physics%powerlaw_c(i,j) for each vertex.

    use cism_parallel, only: staggered_parallel_halo
    use glissade_grid_operators, only: glissade_stagger

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                 & ! number of cells in each horizontal direction
         nglacier                    ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id             ! integer glacier ID in the range (1, nglacier)

    real(dp), dimension(nglacier), intent(in) :: &
         glacier_powerlaw_c          ! glacier-specific powerlaw_c from inversion

    real(dp), dimension(ewn-1,nsn-1), intent(inout) :: &
         basal_physics_powerlaw_c    ! powerlaw_c at each vertex, derived from glacier values

    !TODO - Not sure if the halo update is needed
    type(parallel_type), intent(in) :: &
         parallel                    ! info for parallel communication

    ! local variables

    integer :: i, j, ng

    real(dp), dimension(ewn,nsn) :: &
         powerlaw_c_icegrid            ! powerlaw_c at cell centers, before interpolating to vertices

    integer, dimension(ewn,nsn) :: &
         glacier_mask

    ! Copy glacier_powerlaw_c to a 2D array on the ice grid

    powerlaw_c_icegrid(:,:) = 0.0d0
    do j = 1, nsn
       do i = 1, ewn
          ng = cism_glacier_id(i,j)
          if (ng > 0) powerlaw_c_icegrid(i,j) = glacier_powerlaw_c(ng)
       enddo
    enddo

    ! Compute a mask of cells with glacier ice
    where (cism_glacier_id > 0)
       glacier_mask = 1
    elsewhere
       glacier_mask = 0
    endwhere

    ! Interpolate powerlaw_c to the velocity grid.
    ! At glacier margins, ignore powerlaw_c in cells with glacier_mask = 0
    !  (by setting stagger_margin_in = 1).
    ! Thus, powerlaw_c = 0 at vertices surrounded by cells without glaciers.
    ! Note: This could pose problems if there are dynamically active cells
    !        with cism_glacier_id = 0, but all such cells are currently inactive.

    call glissade_stagger(&
         ewn,                  nsn,    &
         powerlaw_c_icegrid,           &
         basal_physics_powerlaw_c,     &
         ice_mask = glacier_mask,      &
         stagger_margin_in = 1)

    !TODO - Is this update needed?
    call staggered_parallel_halo(basal_physics_powerlaw_c, parallel)

  end subroutine glacier_powerlaw_c_to_2d

!****************************************************

  subroutine glacier_area_volume(&
       ewn,           nsn,               &
       nglacier,      cism_glacier_id,   &
       cell_area,     thck,              &
       area,          volume,            &
       cism_glacier_id_init,             &
       volume_in_init_region,            &
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
         area,                        & ! area of each glacier (m^2)
         volume                         ! volume of each glacier (m^3)

    integer, dimension(ewn,nsn), intent(in), optional ::  &
         cism_glacier_id_init           ! initial value of cism_glacier_id

    real(dp), dimension(nglacier), intent(out), optional ::  &
         volume_in_init_region          ! volume (m^3) in the region defined by cism_glacier_id_init

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
             write(6,'(i8,2f12.6)') ng, area(ng)*1.0d-6, volume(ng)*1.0d-9
          endif
       enddo
    endif

    ! Optionally, compute the volume over the region defined by cism_glacier_id_init.
    ! The idea is that instead of choosing the current glacier volume as a target,
    !  we would match the volume over the initial glacier region.
    ! Then, CISM will not compensate for a too-far-advanced glacier by making it thin,
    !  or for a too-far-retreated glacier by making it thick.

    if (present(volume_in_init_region) .and. present(cism_glacier_id_init)) then
       local_volume(:) = 0.0d0
       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             ng = cism_glacier_id_init(i,j)
             if (ng >= 1) then
                local_volume(ng) = local_volume(ng) + cell_area * thck(i,j)
             endif
          enddo
       enddo
       volume_in_init_region = parallel_reduce_sum(local_volume)
    endif

    ! Optionally, compute the rate of change of glacier volume over the initial glacier region.
    if (present(dthck_dt) .and. present(dvolume_dt) .and. present(cism_glacier_id_init)) then
       ! use local_volume as a work array for dvolume_dt
       dvolume_dt(:) = 0.0d0
       local_volume(:) = 0.0d0
       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             ng = cism_glacier_id_init(i,j)
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

  subroutine accumulate_glacier_fields(&
       ewn,            nsn,                 &
       dt,             time_since_last_avg, &
       snow,           snow_accum,          &
       Tpos,           Tpos_accum,          &
       dthck_dt,       dthck_dt_accum)

    ! input/output variables

    integer, intent(in) ::  &
         ewn, nsn                    ! number of cells in each horizontal direction

    real(dp), intent(in) :: dt       ! time step (yr)

    real(dp), intent(inout) :: &
         time_since_last_avg         ! time (yr) since fields were last averaged

    real(dp), dimension(ewn, nsn), intent(in) ::  &
         snow,                     & ! snowfall rate (mm/yr w.e.)
         Tpos,                     & ! max(artm - T_mlt, 0) (deg C)
         dthck_dt                    ! rate of change of ice thickness (m/yr)

    real(dp), dimension(ewn, nsn), intent(inout) ::  &
         snow_accum,               & ! accumulated snow (mm/yr w.e.)
         Tpos_accum,               & ! accumulated Tpos (deg C)
         dthck_dt_accum              ! rate of change of ice thickness (m/yr)

    time_since_last_avg = time_since_last_avg + dt

    snow_accum = snow_accum + snow * dt
    Tpos_accum = Tpos_accum + Tpos * dt
    dthck_dt_accum = dthck_dt_accum + dthck_dt * dt

  end subroutine accumulate_glacier_fields

!****************************************************

  subroutine calculate_glacier_averages(&
       ewn,            nsn,      &
       time_since_last_avg,      &
       snow_accum,               &
       Tpos_accum,               &
       dthck_dt_accum)

    ! input/output variables

    integer, intent(in) ::  &
         ewn, nsn                    ! number of cells in each horizontal direction

    real(dp), intent(inout) :: &
         time_since_last_avg         ! time (yr) since fields were last averaged

    real(dp), dimension(ewn, nsn), intent(inout) ::  &
         snow_accum,               & ! snow (mm/yr w.e.)
         Tpos_accum,               & ! max(artm - T_mlt, 0) (deg C)
         dthck_dt_accum              ! rate of change of ice thickness (m/yr)

    snow_accum = snow_accum / time_since_last_avg
    Tpos_accum = Tpos_accum / time_since_last_avg
    dthck_dt_accum = dthck_dt_accum / time_since_last_avg

    time_since_last_avg = 0.0d0

  end subroutine calculate_glacier_averages

!****************************************************

  subroutine reset_glacier_fields(&
       ewn,           nsn,        &
       snow_accum,                &
       Tpos_accum,                &
       dthck_dt_accum)

    ! input/output variables

    integer, intent(in) ::  &
         ewn, nsn                    ! number of cells in each horizontal direction

    real(dp), dimension(ewn,nsn), intent(inout) ::  &
         snow_accum,               & ! snow (mm/yr w.e.)
         Tpos_accum,               & ! max(artm - T_mlt, 0) (deg C)
         dthck_dt_accum              ! rate of change of ice thickness (m/yr)

    ! Reset the accumulated fields to zero
    snow_accum = 0.0d0
    Tpos_accum = 0.0d0
    dthck_dt_accum = 0.0d0

  end subroutine reset_glacier_fields

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
