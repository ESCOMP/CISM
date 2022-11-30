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
    !TODO - Add these to the glacier derived type and make them config parameters?

    real(dp), parameter ::  &
         mu_star_const = 500.d0,                  & ! uniform initial value for mu_star (mm/yr w.e/deg C)
         mu_star_min = 2.0d1,                     & ! min value of tunable mu_star (mm/yr w.e/deg C)
         mu_star_max = 2.0d4                        ! max value of tunable mu_star (mm/yr w.e/deg C)

    real(dp), parameter ::  &
         glacier_powerlaw_c_timescale = 100.d0,   & ! inversion timescale for powerlaw_c (yr)
         glacier_powerlaw_c_thck_scale = 100.d0,  & ! inversion thickness scale for powerlaw_c (m)
         glacier_powerlaw_c_relax_factor = 0.05d0   !  controls strength of relaxation to default values (unitless)

    !TODO - Make this an input argument?
    integer, parameter :: &
         inversion_time_interval = 1                ! interval (yr) between inversion calls; must be an integer

    !WHL - Debug
    integer, parameter :: ngtot = 5

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
    real(dp) :: max_glcval

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
       if (associated(glacier%mu_star)) deallocate(glacier%mu_star)
       if (associated(glacier%smb)) deallocate(glacier%smb)
       if (associated(glacier%smb_obs)) deallocate(glacier%smb_obs)
       if (associated(glacier%delta_artm)) deallocate(glacier%delta_artm)

       ! Set the RGI ID to 0 in cells without ice.
       ! Typically, any ice-free cell should already have an RGI ID of 0,
       !  but there can be exceptions related to no-ice boundary conditions.
       where (model%geometry%thck == 0.0d0)
          glacier%rgi_glacier_id = 0
       endwhere

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

          ! Sort the list from low to high values of the RGI IDs.
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
       allocate(glacier%mu_star(nglacier))
       allocate(glacier%smb(nglacier))
       allocate(glacier%smb_obs(nglacier))
       allocate(glacier%delta_artm(nglacier))

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
       glacier%mu_star(:) = mu_star_const
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

       !WHL - debug
       ! For testing, initialize model%climate%smb_obs to something simple.
!!       model%climate%smb_obs(:,:) = 0.d0   ! mm/yr w.e.
!!       model%climate%smb_obs(:,:) = -100.d0   ! mm/yr w.e.
!!       model%climate%smb_obs(:,:) = 100.d0   ! mm/yr w.e.

       ! Given the 2D smb_obs field, compute the 1D glacier-average field.
       ! On restart, this will be read from the restart file.
       call glacier_2d_to_1d(&
            ewn,                   nsn,                           &
            nglacier,              glacier%cism_glacier_id_init,  &
            model%climate%smb_obs, glacier%smb_obs)

       ! If inverting for powerlaw_c, then initialize powerlaw_c to a constant value,
       !  and initialize the inversion target, usrf_obs.
       ! On restart, powerlaw_c and usrf_obs are read from the restart file.
       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          model%basal_physics%powerlaw_c(:,:) = model%basal_physics%powerlaw_c_const
          model%geometry%usrf_obs(:,:) = model%geometry%usrf(:,:)
       endif

       !WHL - debug - Make sure cism_glacier_id_init = 0 where (and only where) rgi_glacier_id > 0
       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             if (glacier%rgi_glacier_id(i,j) > 0) then
                if (glacier%cism_glacier_id_init(i,j) == 0) then
                   write(message,*) 'ERROR: rgi ID, cism ID =', &
                        glacier%rgi_glacier_id(i,j), glacier%cism_glacier_id_init(i,j)
                   call write_log(message, GM_FATAL)
                endif
             endif
             if (glacier%cism_glacier_id_init(i,j) > 0) then
                if (glacier%rgi_glacier_id(i,j) == 0) then
                   write(message,*) 'ERROR: rgi ID, cism ID =', &
                        glacier%rgi_glacier_id(i,j), glacier%cism_glacier_id_init(i,j)
                   call write_log(message, GM_FATAL)
                endif
             endif
          enddo
       enddo

       !WHL - debug - check for cells with thck > 0 and ng = 0
       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             if (glacier%cism_glacier_id_init(i,j) == 0 .and. &
                  model%geometry%thck(i,j)*thk0 > 1.0d0) then
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                print*, 'Warning, ng = 0 but H > 0: Init rank, i, j, ig, jg, thck:', &
                     this_rank, i, j, iglobal, jglobal, model%geometry%thck(i,j) * thk0
             endif
          enddo
       enddo

    else  ! restart

       ! In this case, most required glacier info has already been read from the restart file.
       ! Here, do some error checks and diagnostics.

       ! From the restart file, nglacier is found as the length of dimension 'glacierid'.
       ! The 1D glacier arrays are then allocated with dimension(nglacier) in subroutine glide_allocarr.
       ! The following glacier arrays should be present in the restart file:
       !     rgi_glacier_id, cism_glacier_id, cism_glacier_id_init, cism_to_rgi_glacier_id,
       !     glacier_mu_star, powerlaw_c
       ! If inverting for powerlaw_c, then usrf_obs is read from the restart file.
       ! If inverting for mu_star, then smb_obs is read from the restart file.

       nglacier = glacier%nglacier

       ! Check that the glacier arrays which are read from the restart file have nonzero values.
       ! Note: These arrays are read on all processors.

       max_id = maxval(glacier%cism_glacier_id)
       max_id = parallel_reduce_max(max_id)
       if (max_id <= 0) then
          call write_log ('Error, no positive values for cism_glacier_id', GM_FATAL)
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

       max_glcval = maxval(model%basal_physics%powerlaw_c)
       max_glcval = parallel_reduce_max(max_glcval)
       if (max_glcval <= 0.0d0) then
          call write_log ('Error, no positive values for glacier powerlaw_c', GM_FATAL)
       endif

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          max_glcval = maxval(model%geometry%usrf_obs)
          max_glcval = parallel_reduce_max(max_glcval)
          if (max_glcval <= 0.0d0) then
             call write_log ('Error, no positive values for usrf_obs', GM_FATAL)
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
            glacier%volume)

    endif   ! not a restart

    ! The remaining code applies to both start-up and restart runs

    ! Fill the glacierid dimension array
    do ng = 1, nglacier
       glacier%glacierid(ng) = ng
    enddo

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

    ! Set the relaxation value for powerlaw_c
    if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
       model%basal_physics%powerlaw_c_relax(:,:) = model%basal_physics%powerlaw_c_const
    endif

    ! If not inverting for powerlaw_c, then set delta_artm = 0.
    ! (Need delta_artm = 0 if switching from inversion to no-inversion on restart)
    if (glacier%set_powerlaw_c /= GLACIER_POWERLAW_C_INVERSION) then
       glacier%delta_artm = 0.0d0
    endif

    ! Set the index of the diagnostic glacier, using the CISM glacier ID for the diagnostic point
    if (this_rank == rtest) then
       glacier%ngdiag = glacier%cism_glacier_id(itest,jtest)
    endif
    call broadcast(glacier%ngdiag, rtest)

    ! Write some values for the diagnostic glacier
    if (verbose_glacier .and. this_rank == rtest) then
       i = itest; j = jtest
       ng = glacier%ngdiag
       print*, ' '
       print*, 'Glacier ID for diagnostic cell: r, i, j, ng =', rtest, itest, jtest, ng
       print*, 'area target (km^2) =', glacier%area_target(ng) / 1.0d6
       print*, 'volume target (km^3) =', glacier%volume_target(ng) / 1.0d9
       print*, 'powerlaw_c (Pa (m/yr)^(-1/3)) =', model%basal_physics%powerlaw_c(i,j)
       print*, 'smb_obs (mm/yr w.e.) =', glacier%smb_obs(ng)
       print*, 'Done in glissade_glacier_init'
    endif

  end subroutine glissade_glacier_init

!****************************************************

  !TODO - Pass in precip
  !       Determine whether it's snow based on artm

  subroutine glissade_glacier_smb(&
       ewn,      nsn,                  &
       itest,    jtest,  rtest,        &
       nglacier,                       &
       cism_glacier_id,                &
       t_mlt,                          &
       snow,             artm,         &
       delta_artm,       mu_star,      &
       smb,                            &
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
         delta_artm,                  & ! temperature adjustment to yield SMB ~ 0 (deg C)
         mu_star                        ! glacier-specific SMB tuning parameter (mm w.e./yr/deg)

                                        ! defined as positive for T decreasing with height

    real(dp), dimension(ewn,nsn), intent(out) :: &
         smb                            ! SMB in each gridcell (mm/yr w.e.)

    real(dp), dimension(nglacier), intent(out) :: &
         glacier_smb                    ! average SMB for each glacier (mm/yr w.e.)

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
    smb(:,:) = 0.0d0

    ! compute SMB in each glacier grid cell

    do j = 1, nsn
       do i = 1, ewn
          ng = cism_glacier_id(i,j)

          if (ng > 0) then
             smb(i,j) = snow(i,j) - mu_star(ng) * max(artm(i,j) + delta_artm(ng) - t_mlt, 0.0d0)
          endif

          if (verbose_glacier .and. this_rank == rtest .and. i == itest .and. j == jtest) then
             print*, ' '
             print*, 'Glacier SMB calculation: rank i, j, mu_star =', &
                  this_rank, i, j, mu_star(ng)
             print*, '   snow (mm/yr w.e.), artm (C), delta_artm, T - Tmlt, SMB (mm/yr w.e.) =', &
                  snow(i,j), artm(i,j), delta_artm(ng), max(artm(i,j) + delta_artm(ng) - t_mlt, 0.0d0), smb(i,j)
          endif

       enddo
    enddo

    ! Compute glacier average values

    call glacier_2d_to_1d(&
         ewn,           nsn,              &
         nglacier,      cism_glacier_id,  &
         smb,           glacier_smb)

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
    !   (Actually, minthck is slightly less than thklim, to make sure these cells
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

    integer :: i, j, ii, jj, ip, jp, ipmax, jpmax
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

    ! Check for advance: cells with cism_glacier_id = 0 but H > H_min

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
                               ipmax = ip; jpmax = jp
                            endif
                         endif
                      endif   ! neighbor cell
                   enddo   ! ii
                enddo   ! jj
                if (usrf_max > 0.0d0) then
                   cism_glacier_id(i,j) = cism_glacier_id(ipmax,jpmax)
                   if (verbose_glacier .and. this_rank == rtest) then
                      call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                      print*, 'Set ID = neighbor ID, ig, jg, new ID, thck =', &
                           iglobal, jglobal, cism_glacier_id(i,j), thck(i,j)
                   endif
                endif   ! usrf_max > 0
             endif   ! cism_glacier_id_init > 0

             ! If the cell still has cism_glacier_id = 0 and H > glacier_minthck,
             !  then cap the thickness at glacier_minthck.
             ! Note: The ice removed is used to increment acab_applied, the ice SMB in m/s.
             !       Thus, the total SMB flux can be more negative during time steps
             !        when this subroutine is called.
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
    use glissade_utils, only: glissade_usrf_to_thck
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
         thck_obs,                & ! observed ice thickness (m)
         dthck_dt,                & ! rate of change of thickness (m/yr)
         Tpos,                    & ! max(artm - T_mlt, 0.0)
         Tpos_dartm,              & ! max(artm + delta_artm - T_mlt, 0.0)
         delta_artm_2d,           & ! 2D version of glacier%artm_delta
         mu_star_2d,              & ! 2D version of glacier%mu_star
         smb_annmean                ! annual mean SMB for each glacier cell (mm/yr w.e.)

    real(dp), dimension(model%general%ewn-1, model%general%nsn-1) ::  &
         stag_thck,                   & ! ice thickness at vertices (m)
         stag_thck_obs,               & ! observed ice thickness at vertices (m)
         stag_dthck_dt                  ! rate of change of ice thickness at vertices (m/yr)

    type(parallel_type) :: parallel ! info for parallel communication

    real(dp), save ::  &              ! time since the last averaging computation;
         time_since_last_avg = 0.0d0  ! compute the average once a year

    real(dp), dimension(glacier%nglacier) :: &
         area_old,                & ! glacier%area from the previous inversion step
         darea_dt,                & ! rate of change of glacier area over the inversion interval
         smb_init_area,           & ! SMB over initial area determined by cism_glacier_id_init (should be ~smb_obs)
         smb_init_area_dartm,     & ! Same as smb_init_area, but with the corrected artm (should be ~ 0)
         smb_current_area_dartm     ! SMB over current area determined by cism_glacier_id, with the corrected artm
                                    ! (should eventually approach 0)

    ! Note: The glacier type includes the following:
    ! integer ::  nglacier          ! number of glaciers in the global domain
    ! integer ::  ngdiag            ! CISM index of diagnostic glacier
    ! real(dp), dimension(:) :: area              ! glacier area (m^2)
    ! real(dp), dimension(:) :: volume            ! glacier volume (m^3)
    ! real(dp), dimension(:) :: area_target       ! glacier area target (m^2)
    ! real(dp), dimension(:) :: volume_target     ! glacier volume target (m^3)
    ! real(dp), dimension(:) :: mu_star           ! SMB parameter for each glacier (mm/yr w.e./deg K)
    ! real(dp), dimension(:) :: smb_obs           ! observed SMB for each glacier (mm/yr w.e.)
    ! integer, dimension(:,:) :: cism_glacier_id       ! CISM glacier ID for each grid cell
    ! integer, dimension(:,:) :: cism_glacier_id_init  ! initial value of CISM glacier ID
    ! real(dp), dimension(:,:) :: dthck_dt_2d          ! snow accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: snow_2d              ! snow accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: Tpos_2d              ! max(artm-T_mlt,0) accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: Tpos_dartm_2d        ! max(artm+delta_artm-T_mlt,0) accumulated and averaged over 1 year

    ! Set some local variables

    parallel = model%parallel

    ewn = model%general%ewn
    nsn = model%general%nsn
    dew = model%numerics%dew * len0         ! convert to m
    dns = model%numerics%dns * len0         ! convert to m
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    nglacier = glacier%nglacier
    ngdiag = glacier%ngdiag

    ! some unit conversions
    dt = model%numerics%dt * tim0/scyr          ! model units to yr
    thck = model%geometry%thck * thk0           ! model units to m
    dthck_dt = model%geometry%dthck_dt * scyr   ! m/s to m/yr

    ! Compute the current area and volume of each glacier.
    ! These are not needed for inversion, but are computed as diagnostics.
    ! Note: This requires global sums. For now, do the computation independently on each task.

    call glacier_area_volume(&
         ewn,           nsn,              &
         nglacier,                        &
         glacier%cism_glacier_id,         &
         dew*dns,                         &  ! m^2
         model%geometry%thck * thk0,      &  ! m
         glacier%area,                    &  ! m^2
         glacier%volume)                     ! m^3

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'Update area (km^2) and volume (km^3) for glacier:', ngdiag
       print*, 'Current area and volume:', &
            glacier%area(ngdiag)/1.0d6, glacier%volume(ngdiag)/1.0d9
       print*, '   Target area and volume:', &
            glacier%area_target(ngdiag)/1.0d6, glacier%volume_target(ngdiag)/1.0d9
       print*, ' '
       print*, ngtot, 'glaciers: ng, A, A_target, Aerr, V, V_target, Verr:'
       do ng = 1, ngtot
          write(6,'(i6,3f12.4,3f14.6)') ng, glacier%area(ng)/1.0d6, glacier%area_target(ng)/1.0d6, &
               (glacier%area(ng) - glacier%area_target(ng))/1.0d6, &
               glacier%volume(ng)/1.0d9, glacier%volume_target(ng)/1.0d9, &
               (glacier%volume(ng) - glacier%volume_target(ng))/1.0d9
       enddo
    endif

    ! Invert for mu_star and/or powerlaw_c
    ! Note: Tpos is based on the input air temperature, artm.
    !       During the inversion, we choose mu_star such that smb = smb_obs for each glacier.
    !       Tpos_dartm is based on artm along with artm_delta, where artm_delta is an adjustment term
    !        that results in smb ~ 0. Correcting the SMB inhibits glacier advance and retreat
    !        during the spin-up, which makes it possible to invert for powerlaw_c in a quasi-steady state.

    if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .or.  &
        glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then

       ! Accumulate the 2D fields used for mu_star inversion: snow and Tpos
       ! Also accumulate dthck_dt and Tpos_dartm, which are used for powerlaw_c inversion

       if (time_since_last_avg == 0.0d0) then ! start of new averaging period

          ! Reset the accumulated fields to zero
          call reset_glacier_fields(&
               ewn,           nsn,     &
               glacier%snow_2d,        &
               glacier%Tpos_2d,        &
               glacier%Tpos_dartm_2d,  &
               glacier%dthck_dt_2d)
       endif

       Tpos(:,:) = max(model%climate%artm(:,:) - glacier%t_mlt, 0.0d0)

       ! Given delta_artm for each glacier, scatter values to the 2D CISM grid

       call glacier_1d_to_2d(&
            ewn,      nsn,   &
            nglacier,        &
            glacier%cism_glacier_id_init,  &
            glacier%delta_artm,   &
            delta_artm_2d)

       Tpos_dartm(:,:) = max(model%climate%artm(:,:) + delta_artm_2d(:,:) - glacier%t_mlt, 0.0d0)

       ! Accumulate Tpos_2d, Tpos_dartm_2d, and dthck_dt_2d over this timestep

       call accumulate_glacier_fields(&
            ewn,                    nsn,                    &
            dt,                     time_since_last_avg,    &
            model%climate%snow,     glacier%snow_2d,        &  ! mm/yr w.e.
            Tpos,                   glacier%Tpos_2d,        &  ! deg C
            Tpos_dartm,             glacier%Tpos_dartm_2d,  &  ! deg C
            dthck_dt,               glacier%dthck_dt_2d)       ! m/yr ice

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'In glissade_glacier_inversion, diag cell (r, i, j) =', rtest, itest, jtest
          i = itest; j = jtest
          print*, 'r, i, j, time, time_since_last_avg, snow, Tpos, Tpos_dartm:', &
               this_rank, i, j, model%numerics%time, time_since_last_avg, &
               glacier%snow_2d(i,j), glacier%Tpos_2d(i,j), glacier%Tpos_dartm_2d(i,j)
       endif

       ! Check whether it is time to do the inversion.
       ! Note: model%numerics%time has units of yr.
       !       inversion_time_inveral is an integer number of years.

       if (abs(time_since_last_avg - real(inversion_time_interval,dp)) < eps08) then

          if (verbose_glacier .and. this_rank == rtest) then
             print*, 'calculate_glacier_2d_to_1ds, time_since_last_avg =', time_since_last_avg
          endif

          ! Compute the average of glacier fields over the accumulation period

          call glacier_time_averages(&
               ewn,                    nsn,   &
               time_since_last_avg,           &  ! yr
               glacier%snow_2d,               &  ! mm/yr w.e.
               glacier%Tpos_2d,               &  ! deg C
               glacier%Tpos_dartm_2d,         &  ! deg C
               glacier%dthck_dt_2d)              ! m/yr ice

          if (verbose_glacier .and. this_rank == rtest) then
             i = itest; j = jtest
             print*, ' '
             print*, 'Annual averages, r, i, j:', rtest, itest, jtest
             print*, '   snow (mm/yr w.e.)  =', glacier%snow_2d(i,j)
             print*, '   Tpos (deg C)       =', glacier%Tpos_2d(i,j)
             print*, '   Tpos_dartm (deg C) =', glacier%Tpos_dartm_2d(i,j)
             print*, '   dthck_dt (m/yr)    =', glacier%dthck_dt_2d(i,j)
          endif

          ! Invert for mu_star

          if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION) then

             ! Choose mu_star for each glacier to match smb_obs over the initial glacier footprint.
             ! Note: glacier%smb_obs and glacier%mu_star are 1D, per-glacier fields.

             call glacier_invert_mu_star(&
                  ewn,                nsn,                 &
                  nglacier,           ngdiag,              &
                  glacier%cism_glacier_id_init,            &
                  glacier%smb_obs,                         &
                  glacier%snow_2d,    glacier%Tpos_2d,     &
                  glacier%mu_star)

             ! Given these values of mu_star, compute the average SMB for each glacier,
             !  based on its initial area and its current area (for diagnostic purposes only).

             ! Convert mu_star to a 2D field

             call glacier_1d_to_2d(&
                  ewn,             nsn,                           &
                  nglacier,        glacier%cism_glacier_id_init,  &
                  glacier%mu_star, mu_star_2d)

             ! Compute the SMB for each grid cell, given the appropriate mu_star

             smb_annmean(:,:) = glacier%snow_2d(:,:) - mu_star_2d(:,:) * glacier%Tpos_2d(:,:)

             ! Compute the average SMB for each glacier over the initial glacier area

             call glacier_2d_to_1d(&
                  ewn,             nsn,                          &
                  nglacier,        glacier%cism_glacier_id_init, &
                  smb_annmean,     smb_init_area)

             ! Repeat using the delta_artm correction

             smb_annmean(:,:) = glacier%snow_2d(:,:) - mu_star_2d(:,:) * glacier%Tpos_dartm_2d(:,:)

             call glacier_2d_to_1d(&
                  ewn,             nsn,                          &
                  nglacier,        glacier%cism_glacier_id_init, &
                  smb_annmean,     smb_init_area_dartm)

             ! Repeat for the current glacier area, with the delta_artm correction

             call glacier_1d_to_2d(&
                  ewn,             nsn,                      &
                  nglacier,        glacier%cism_glacier_id,  &
                  glacier%mu_star, mu_star_2d)

             smb_annmean(:,:) = glacier%snow_2d(:,:) - mu_star_2d(:,:) * glacier%Tpos_dartm_2d(:,:)

             call glacier_2d_to_1d(&
                  ewn,             nsn,                     &
                  nglacier,        glacier%cism_glacier_id, &
                  smb_annmean,     smb_current_area_dartm)


             if (verbose_glacier .and. this_rank == rtest) then
                print*, ' '
                ng = ngdiag
                print*, 'ngdiag, smb_init_area (mm/yr w.e.), smb_init_area_dartm, smb_current_area_dartm, mu_star:'
                write(6,'(i6,f20.14, 3f12.4)') ng, smb_init_area(ng), smb_init_area_dartm(ng), &
                     smb_current_area_dartm(ng), glacier%mu_star(ng)
                print*, ' '
                print*, ngtot, 'glaciers: smb_init_area (mm/yr w.e.), smb_init_area_dartm, smb_current_area, mu_star:'
                do ng = 1, ngtot
                   write(6,'(i6,f20.14, 3f12.4)') ng, smb_init_area(ng), smb_init_area_dartm(ng), &
                        smb_current_area_dartm(ng), glacier%mu_star(ng)
                enddo
             endif

          endif   ! invert for mu_star

          ! Given the current and target ice thickness, invert for powerlaw_c

          if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then

             ! When inverting for powerlaw_c, we want the glacier footprint to match the observed footprint
             !  as closely as possible.
             ! This is done by adjusting the surface temperature (artm) such that the modeled SMB is close to zero
             !  over the original glacier footprint.
             ! Here, we update delta_artm for each glacier such that SMB is close to zero.
             ! May not have SMB exactly zero because of the max term in the SMB formula.
             !
             ! If snow_2d - mu_star * Tpos_dartm_2d > 0, delta_artm will become more negative
             ! If snow_2d - mu_star * Tpos_dartm_2d < 0, delta_artm will become more positive
             !
             ! Note: The value of delta_artm computed here is not used directly for powerlaw_c inversion.
             !       Rather, it is passed into subroutine glissade_glacier_smb to minimize the change
             !        in the glacier footprint during the spin-up.

             if (verbose_glacier .and. this_rank == rtest) then
                print*, ' '
                print*, ngtot, 'glaciers: initial delta_artm'
                do ng = 1, ngtot
                   write(6,'(i6,2f12.4)') ng, glacier%delta_artm(ng)
                enddo
             endif

             call glacier_adjust_artm(&
                  ewn,              nsn,         &
                  nglacier,         ngdiag,      &
                  glacier%cism_glacier_id_init,  &
                  glacier%snow_2d,               &
                  glacier%Tpos_dartm_2d,         &
                  glacier%mu_star,               &
                  glacier%delta_artm)

             if (verbose_glacier .and. this_rank == rtest) then
                print*, ' '
                print*, ngtot, 'glaciers: new delta_artm'
                do ng = 1, ngtot
                   write(6,'(i6,f12.4)') ng, glacier%delta_artm(ng)
                enddo
             endif

             ! Given the surface elevation target, compute the thickness target.
             ! (This can change in time if the bed topography is dynamic.)
             call glissade_usrf_to_thck(&
                  model%geometry%usrf_obs * thk0,  &
                  model%geometry%topg * thk0,      &
                  model%climate%eus * thk0,        &
                  thck_obs)

             ! Interpolate thck_obs to the staggered grid
             call glissade_stagger(ewn,         nsn,              &
                                   thck_obs,    stag_thck_obs)

             ! Interpolate thck to the staggered grid
             call glissade_stagger(ewn,         nsn,              &
                                   thck,        stag_thck)

             ! Interpolate dthck_dt to the staggered grid
             call glissade_stagger(ewn,                 nsn,           &
                                   glacier%dthck_dt_2d, stag_dthck_dt)

             if (verbose_glacier .and. this_rank == rtest) then
                print*, ' '
                print*, 'call glacier_invert_powerlaw_c, time (yr) =', model%numerics%time
             endif

             call glacier_invert_powerlaw_c(&
                  ewn,              nsn,                 &
                  itest,   jtest,   rtest,               &
                  model%basal_physics%powerlaw_c_min,    &
                  model%basal_physics%powerlaw_c_max,    &
                  stag_thck,        stag_thck_obs,       &
                  stag_dthck_dt,                         &
                  model%basal_physics%powerlaw_c_relax,  &
                  model%basal_physics%powerlaw_c)

          endif   ! powerlaw_c_inversion

       endif   ! time to do inversion

    endif   ! invert for mu_star or powerlaw_c

  end subroutine glissade_glacier_inversion

!****************************************************

  subroutine glacier_invert_mu_star(&
       ewn,              nsn,           &
       nglacier,         ngdiag,        &
       cism_glacier_id_init,            &
       glacier_smb_obs,                 &
       snow_2d,          Tpos_2d,       &
       mu_star)

    ! Given an observational SMB target, invert for the parameter mu_star in the glacier SMB formula

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier,                    & ! total number of glaciers in the domain
         ngdiag                         ! CISM ID of diagnostic glacier

    integer, dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init           ! cism_glacier_id at the start of the run

    real(dp), dimension(nglacier), intent(in) :: &
         glacier_smb_obs                ! observed glacier-average SMB (mm/yr w.e.)

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow_2d,                     & ! time-avg snowfall for each cell (mm/yr w.e.)
         Tpos_2d                        ! time-avg of max(artm - T_mlt, 0) for each cell (deg)

    real(dp), dimension(nglacier), intent(inout) :: &
         mu_star                        ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)

    ! local variables
    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         glacier_snow, glacier_Tpos     ! glacier-average snowfall and Tpos

    character(len=100) :: message

    ! Inversion for mu_star is more direct than inversion for powerlaw_c.
    ! Instead of solving a damped harmonic oscillator equation for mu_star,
    !  we compute mu_star for each glacier such that SMB = smb_obs over the initial extent.
    !
    ! The SMB for glacier ng is given by
    !      sum_ij(smb) = sum_ij(snow) - mu_star(ng) * sum_ij(Tpos),
    ! where Tpos = max(artm - T_mlt, 0),
    ! and sum_ij notes a sum over all cells (i,j) in the glacier.
    !
    ! Rearranging, we get
    !      mu_star(ng) = (sum_ij(snow) - sum_ij(smb) / sum_ij(Tpos)
    !
    ! Thus, given the annual average of snow and Tpos for each grid cell in a glacier,
    ! we can find mu_star such that SMB = smb_obs.
    !
    ! We take sums are taken over the target area of each glacier, using cism_glacier_id_init.
    ! If a glacier is too large, the modeled SMB will be < 0 and the glacier should shrink.
    ! Similarly, if the glacier is too small, the modeled SMB > 0 and the glacier should grow.
    !
    ! Notes:
    !
    ! (1) This approach works only for land-based glaciers.
    !     TODO: Modify for marine-terminating glaciers.
    ! (2) Assuming climatological forcing with smb_obs prescribed, mu_star will have nearly the same value
    !     throughout the inversion.  It changes slightly as surface elevation changes, modifying the downscaled Tpos.

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'In glacier_invert_mu_star'
    endif

    ! Compute average snowfall, Tpos, and SMB over the initial extent of each glacier

    call glacier_2d_to_1d(&
         ewn,           nsn,                   &
         nglacier,      cism_glacier_id_init,  &
         snow_2d,       glacier_snow)

    call glacier_2d_to_1d(&
         ewn,           nsn,                   &
         nglacier,      cism_glacier_id_init,  &
         Tpos_2d,       glacier_Tpos)

    ! For each glacier, compute the new mu_star

    do ng = 1, nglacier

       if (glacier_Tpos(ng) > 0.0d0) then  ! ablation is nonzero

          ! Compute the value of mu_star that will give the desired SMB over the target area
          mu_star(ng) = (glacier_snow(ng) - glacier_smb_obs(ng)) / glacier_Tpos(ng)

          ! Limit to a physically reasonable range
          mu_star(ng) = min(mu_star(ng), mu_star_max)
          mu_star(ng) = max(mu_star(ng), mu_star_min)

          if (verbose_glacier .and. main_task .and. ng == ngdiag) then
             print*, ' '
             print*, 'ng, glacier-average snow, Tpos, smb_obs:', &
                  ng, glacier_snow(ng), glacier_Tpos(ng), glacier_smb_obs(ng)
             print*, 'New mu_star:', mu_star(ng)
          endif

       else   ! glacier_Tpos = 0; no ablation

          mu_star(ng) = mu_star_max

          if (verbose_glacier .and. main_task) then
             print*, 'Warning: no ablation for glacier', ng
          endif

       endif

    enddo   ! ng

  end subroutine glacier_invert_mu_star

!****************************************************

  subroutine glacier_adjust_artm(&
       ewn,              nsn,                &
       nglacier,         ngdiag,             &
       cism_glacier_id_init,                 &
       snow_2d,          Tpos_dartm_2d,      &
       mu_star,          delta_artm)

    ! Given mu_star for each glacier, compute a temperature correction delta_artm
    ! that will nudge the SMB toward zero over the initial glacier footprint.

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier,                    & ! total number of glaciers in the domain
         ngdiag                         ! CISM ID of diagnostic glacier

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow_2d,                  & ! time-avg snowfall for each cell (mm/yr w.e.)
         Tpos_dartm_2d               ! time-avg of max(artm + delta_artm - T_mlt, 0) for each cell (deg)

    integer, dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init           ! cism_glacier_id at the start of the run

    real(dp), dimension(nglacier), intent(in) :: &
         mu_star                        ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)

    real(dp), dimension(nglacier), intent(inout) :: &
         delta_artm                     ! glacier-specific temperature correction (deg)

    ! local variables
    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         glacier_snow, glacier_Tpos_dartm   ! average snow and Tpos for each glacier

    real(dp) :: artm_correction

    ! The SMB for glacier ng is given by
    !      sum_ij(smb) = sum_ij(snow) - mu_star(ng) * sum_ij(Tpos_dartm),
    ! where Tpos_dartm = max(artm + delta_artm  - T_mlt, 0),
    ! and sum_ij notes a sum over all cells (i,j) in the glacier.
    !
    ! We set SMB = 0 and replacing Tpos_dartm with Tpos_dartm + artm_correction,
    ! where we want to find artm_correction.
    !
    ! Rearranging, we get
    !
    !      artm_correction = (sum_ij(snow) - mu_star*sum_ij(Tpos_dartm)) / mu_star
    !

    ! Compute the average of snow_2d and Tpos_dartm_2d over each glacier

    call glacier_2d_to_1d(&
         ewn,              nsn,         &
         nglacier,                      &
         cism_glacier_id_init,          &
         snow_2d,                       &
         glacier_snow)

    call glacier_2d_to_1d(&
         ewn,              nsn,         &
         nglacier,                      &
         cism_glacier_id_init,          &
         Tpos_dartm_2d,                 &
         glacier_Tpos_dartm)

    ! For each glacier, compute the new delta_artm
    ! Note: Because of the threshold T > T_mlt for contributing to Tpos,
    !        not all the temperature change may be effective in increasing
    !        or decreasing ablation.
    !       So we may not end up with SMB = 0, but we will approach that target
    !        over several timesteps.

    do ng = 1, nglacier
       artm_correction = (glacier_snow(ng) - mu_star(ng) * glacier_Tpos_dartm(ng))   &
                      / mu_star(ng)
       delta_artm(ng) = delta_artm(ng) + artm_correction

       if (verbose_glacier .and. main_task .and. ng == ngdiag) then
          print*, ' '
          print*, 'glacier_adjust_artm, ng =', ng
          print*, 'glacier-average snow, Tpos_dartm, mu_star:', &
               glacier_snow(ng), glacier_Tpos_dartm(ng), mu_star(ng)
          print*, 'artm correction =', artm_correction
          print*, 'New delta_artm  =', delta_artm(ng)
       endif

    enddo

  end subroutine glacier_adjust_artm

!****************************************************

  subroutine glacier_invert_powerlaw_c(&
       ewn,              nsn,                 &
       itest,   jtest,   rtest,               &
       powerlaw_c_min,   powerlaw_c_max,      &
       stag_thck,        stag_thck_obs,       &
       stag_dthck_dt,                         &
       powerlaw_c_relax, powerlaw_c)

    ! Given the current ice thickness, rate of thickness change, and target thickness,
    ! invert for the parameter powerlaw_c in the relationship for basal sliding.
    ! Note: This subroutine is similar to subroutine invert_basal_friction
    !       in the glissade_inversion_module.  It is separate so that we can experiment
    !       with glacier inversion parameters without changing the standard ice sheet inversion.
    !       The glacier inversion parameters are currently declared at the top of this module.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest            ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         powerlaw_c_min, powerlaw_c_max ! min and max allowed values of powerlaw_c (Pa (m/yr)^(-1/3))

    real(dp), dimension(ewn-1,nsn-1), intent(in) :: &
         stag_thck,                   & ! ice thickness at vertices (m)
         stag_thck_obs,               & ! observed ice thickness at vertices (m)
         stag_dthck_dt                  ! rate of change of ice thickness at vertices (m/yr)

    real(dp), dimension(ewn-1,nsn-1), intent(in) :: &
         powerlaw_c_relax               ! powerlaw_c field to which we relax

    real(dp), dimension(ewn-1,nsn-1), intent(inout) :: &
         powerlaw_c                     ! basal friction field to be adjusted (Pa (m/yr)^(-1/3))

    ! local variables

    integer :: i, j

    real(dp), dimension(ewn-1,nsn-1) :: &
         stag_dthck                     ! stag_thck - stag_thck_obs (m)

    real(dp) :: &
         dpowerlaw_c,                 & ! change in powerlaw_c
         term_thck, term_dHdt,        & ! tendency terms for powerlaw_c based on thickness target
         term_relax                     ! tendency terms based on relaxation to default value

    ! The inversion works as follows:
    ! The change in C_p is proportional to the current value of C_p and to the relative error,
    !  err_H = (H - H_obs)/H_scale, where H is a thickness scale.
    ! If err_H > 0, we reduce C_p to make the ice flow faster and thin.
    ! If err_H < 0, we increase C_p to make the ice flow slower and thicken.
    ! This is done with a characteristic timescale tau.
    ! We also include a term proportional to dH/dt so that ideally, C_p smoothly approaches
    !  the value needed to attain a steady-state H, without oscillating about the desired value.
    ! In addition, we include a relaxation term proportional to the ratio of C_p to a default value.
    ! See the comments in module glissade_inversion, subroutine invert_basal_friction.
    !
    ! Here is the prognostic equation:
    ! dC/dt = -C * [(H - H_obs)/(H0*tau) + dH/dt * 2/H0 - r * ln(C/C_r) / tau],
    !   where tau = glacier_powerlaw_c_timescale, H0 = glacier_powerlaw_c_thck_scale,
    !         r = glacier_powerlaw_c_relax_factor, and C_r = powerlaw_c_relax.

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'In glacier_invert_powerlaw_c'
    endif

    if (glacier_powerlaw_c_thck_scale > 0.0d0 .and. glacier_powerlaw_c_timescale > 0.0d0) then

       stag_dthck(:,:) = stag_thck(:,:) - stag_thck_obs(:,:)

       ! Loop over vertices

       do j = 1, nsn-1
          do i = 1, ewn-1

             if (stag_thck(i,j) > 0.0d0) then

                term_thck = -stag_dthck(i,j) / (glacier_powerlaw_c_thck_scale * glacier_powerlaw_c_timescale)
                term_dHdt = -stag_dthck_dt(i,j) * 2.0d0 / glacier_powerlaw_c_thck_scale

                ! Add a term to relax C = powerlaw_c toward a target value, C_r = powerlaw_c_relax
                ! The log term below ensures the following:
                ! * When C /= C_r, it will relax toward C_r.
                ! * When C = C_r, there is no further relaxation.
                ! * In steady state (dC/dt = 0, dH/dt = 0), we have dthck/thck_scale = -k * ln(C/C_r),
                !    or C = C_r * exp(-dthck/(k*thck_scale)), where k is a prescribed constant

                term_relax = -glacier_powerlaw_c_relax_factor * log(powerlaw_c(i,j)/powerlaw_c_relax(i,j)) &
                     / glacier_powerlaw_c_timescale

                dpowerlaw_c = powerlaw_c(i,j) * (term_thck + term_dHdt + term_relax) * inversion_time_interval

                ! Limit to prevent a large relative change in one step
                if (abs(dpowerlaw_c) > 0.05d0 * powerlaw_c(i,j)) then
                   if (dpowerlaw_c > 0.0d0) then
                      dpowerlaw_c =  0.05d0 * powerlaw_c(i,j)
                   else
                      dpowerlaw_c = -0.05d0 * powerlaw_c(i,j)
                   endif
                endif

                ! Update powerlaw_c
                powerlaw_c(i,j) = powerlaw_c(i,j) + dpowerlaw_c

                ! Limit to a physically reasonable range
                powerlaw_c(i,j) = min(powerlaw_c(i,j), powerlaw_c_max)
                powerlaw_c(i,j) = max(powerlaw_c(i,j), powerlaw_c_min)

                if (verbose_glacier .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                   print*, ' '
                   print*, 'Invert for powerlaw_c: rank, i, j =', this_rank, i, j
                   print*, 'H, H_obs (m)', stag_thck(i,j), stag_thck_obs(i,j)
                   print*, 'dH_dt (m/yr):', stag_dthck_dt(i,j)
                   print*, 'dt (yr), term_thck*dt, term_dHdt*dt:', inversion_time_interval, &
                        term_thck*inversion_time_interval, term_dHdt*inversion_time_interval
                   print*, 'relax term:', term_relax*inversion_time_interval
                   print*, 'dpowerlaw_c, new powerlaw_c:', dpowerlaw_c, powerlaw_c(i,j)
                endif

             else   ! stag_thck = 0

                ! do nothing; keep the current value

             endif

          enddo   ! i
       enddo   ! j

    else   ! thck_scale or timescale = 0

       call write_log &
            ('Must have thck_scale and timescale > 0 for glacier powerlaw_c inversion', GM_FATAL)

    endif

    if (verbose_glacier .and. this_rank == rtest) then
       i = itest; j = jtest
       print*, ' '
       print*, 'stag_thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') stag_thck(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'stag_thck - stag_thck_obs (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') stag_dthck(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'stag_dthck_dt (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') stag_dthck_dt(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'new powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.0)',advance='no') powerlaw_c(i,j)
          enddo
          print*, ' '
       enddo
    endif   ! verbose_glacier

  end subroutine glacier_invert_powerlaw_c

!****************************************************

  subroutine glacier_2d_to_1d(&
       ewn,           nsn,              &
       nglacier,      cism_glacier_id,  &
       field_2d,      glacier_field)

    ! Given a 2D field, compute the average of the field over each glacier

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         field_2d                       ! 2D field to be averaged over glaciers

    real(dp), dimension(nglacier), intent(out) ::  &
         glacier_field                  ! field average over each glacier

    ! local variables

    integer :: i, j, ng

    integer, dimension(nglacier) :: ncells_glacier

    ncells_glacier(:) = 0
    glacier_field(:) = 0.0d0

    ! Loop over locally owned cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0) then
             ncells_glacier(ng) = ncells_glacier(ng) + 1
             glacier_field(ng) = glacier_field(ng) + field_2d(i,j)
          endif
       enddo
    enddo

    ncells_glacier = parallel_reduce_sum(ncells_glacier)
    glacier_field  = parallel_reduce_sum(glacier_field)

    where (ncells_glacier > 0)
       glacier_field = glacier_field/ncells_glacier
    endwhere

  end subroutine glacier_2d_to_1d

!****************************************************

  subroutine glacier_1d_to_2d(&
       ewn,           nsn,              &
       nglacier,      cism_glacier_id,  &
       glacier_field, field_2d)

    ! Given a 1D per-glacier field, scatter the values to the 2D grid.
    ! Each cell in a given glacier will have the same value.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), dimension(nglacier), intent(in) ::  &
         glacier_field                  ! field average over each glacier

    real(dp), dimension(ewn,nsn), intent(out) ::  &
         field_2d                       ! 2D field to be averaged over glaciers

    ! local variables

    integer :: i, j, ng

    field_2d(:,:) = 0.0d0

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0) then
             field_2d(i,j) = glacier_field(ng)
          endif
       enddo   ! i
    enddo   ! j

  end subroutine glacier_1d_to_2d

!****************************************************

  subroutine glacier_area_volume(&
       ewn,           nsn,               &
       nglacier,      cism_glacier_id,   &
       cell_area,     thck,              &
       area,          volume)

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
       print*, ' '
       print*, 'Compute glacier area and volume; cell_area (m^3) =', cell_area
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

    deallocate(local_area)
    deallocate(local_volume)

  end subroutine glacier_area_volume

!****************************************************

  subroutine accumulate_glacier_fields(&
       ewn,            nsn,                 &
       dt,             time_since_last_avg, &
       snow,           snow_2d,             &
       Tpos,           Tpos_2d,             &
       Tpos_dartm,     Tpos_dartm_2d,       &
       dthck_dt,       dthck_dt_2d)

    ! input/output variables

    integer, intent(in) ::  &
         ewn, nsn                    ! number of cells in each horizontal direction

    real(dp), intent(in) :: dt       ! time step (yr)

    real(dp), intent(inout) :: &
         time_since_last_avg         ! time (yr) since fields were last averaged

    real(dp), dimension(ewn, nsn), intent(in) ::  &
         snow,                     & ! snowfall rate (mm/yr w.e.)
         Tpos,                     & ! max(artm - T_mlt, 0) (deg C)
         Tpos_dartm,               & ! max(artm + delta_artm - T_mlt, 0) (deg C)
         dthck_dt                    ! rate of change of ice thickness (m/yr)

    real(dp), dimension(ewn, nsn), intent(inout) ::  &
         snow_2d,                  & ! accumulated snow (mm/yr w.e.)
         Tpos_2d,                  & ! accumulated Tpos (deg C)
         Tpos_dartm_2d,            & ! accumulated Tpos (deg C)
         dthck_dt_2d                 ! rate of change of ice thickness (m/yr)

    time_since_last_avg = time_since_last_avg + dt

    snow_2d = snow_2d + snow * dt
    Tpos_2d = Tpos_2d + Tpos * dt
    Tpos_dartm_2d = Tpos_dartm_2d + Tpos_dartm * dt
    dthck_dt_2d = dthck_dt_2d + dthck_dt * dt

  end subroutine accumulate_glacier_fields

!****************************************************

  subroutine glacier_time_averages(&
       ewn,            nsn,      &
       time_since_last_avg,      &
       snow_2d,                  &
       Tpos_2d,                  &
       Tpos_dartm_2d,            &
       dthck_dt_2d)

    ! input/output variables

    integer, intent(in) ::  &
         ewn, nsn                    ! number of cells in each horizontal direction

    real(dp), intent(inout) :: &
         time_since_last_avg         ! time (yr) since fields were last averaged

    real(dp), dimension(ewn, nsn), intent(inout) ::  &
         snow_2d,                  & ! snow (mm/yr w.e.)
         Tpos_2d,                  & ! max(artm - T_mlt, 0) (deg C)
         Tpos_dartm_2d,            & ! max(artm + delta_artm - T_mlt, 0) (deg C)
         dthck_dt_2d                 ! rate of change of ice thickness (m/yr)

    snow_2d = snow_2d / time_since_last_avg
    Tpos_2d = Tpos_2d / time_since_last_avg
    Tpos_dartm_2d = Tpos_dartm_2d / time_since_last_avg
    dthck_dt_2d = dthck_dt_2d / time_since_last_avg

    time_since_last_avg = 0.0d0

  end subroutine glacier_time_averages

!****************************************************

  subroutine reset_glacier_fields(&
       ewn,        nsn,        &
       snow_2d,                &
       Tpos_2d,                &
       Tpos_dartm_2d,          &
       dthck_dt_2d)

    ! input/output variables

    integer, intent(in) ::  &
         ewn, nsn                 ! number of cells in each horizontal direction

    real(dp), dimension(ewn,nsn), intent(inout) ::  &
         snow_2d,               & ! snow (mm/yr w.e.)
         Tpos_2d,               & ! max(artm - T_mlt, 0) (deg C)
         Tpos_dartm_2d,         & ! max(artm + delta_artm - T_mlt, 0) (deg C)
         dthck_dt_2d              ! rate of change of ice thickness (m/yr)

    ! Reset the accumulated fields to zero
    snow_2d = 0.0d0
    Tpos_2d = 0.0d0
    Tpos_dartm_2d = 0.0d0
    dthck_dt_2d = 0.0d0

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
