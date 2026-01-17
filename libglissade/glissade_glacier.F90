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
    use glimmer_paramets, only: iulog, eps08
    use glimmer_physcon, only: scyr, pi, rhow, rhoi
    use glide_types
    use glimmer_log
    use glimmer_utils, only: point_diag
    use cism_parallel, only: main_task, this_rank, nhalo

    implicit none

    private
    public :: verbose_glacier, glissade_glacier_init, glissade_glacier_update

    logical, parameter :: verbose_glacier = .true.

    ! derived type that holds info for each glaciated grid cell
    type glacier_info
       integer :: id           ! input glacier ID, usually RGI
       integer :: indxi        ! i index of cell
       integer :: indxj        ! j index of cell
    end type glacier_info

    ! Glacier parameters used in this module

    !TODO - Make this an input argument?
    integer, parameter :: &
         glacier_update_interval = 1    ! interval (yr) between inversion calls and other glacier updates

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

    use cism_parallel, only: gather_var, scatter_var, &
         parallel_global_sum, parallel_reduce_max, parallel_reduce_min, parallel_is_zero, &
         broadcast, parallel_halo, staggered_parallel_halo, parallel_globalindex

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
    integer :: ng_ne, ng_nw, ng_se, ng_sw
    integer :: min_id, max_id
    real(dp) :: max_glcval
    real(dp) :: theta_rad             ! latitude in radians

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

    integer, dimension(model%general%ewn,model%general%nsn) :: &
         glacier_mask                 ! = 1 for cells with glaciers (glacier_id > 0), else = 0

    type(parallel_type) :: parallel   ! info for parallel communication

    !WHL - debug, for quicksort test
!    integer, dimension(:), allocatable :: test_list
!    integer ::  nlist
!    real(sp) :: random

    ! Optional grid cell dimension correction
    ! Note: The following is an awkward way of dealing with the fact that for some of our glacier grids,
    !        the nominal grid dimensions in the input file are different from the true dimensions.
    !       For instance, we can have a 200-m input grid for glaciers at 45 N (e.g., in the Alps).
    !       The nominal cell size, 200 m, corresponds to the cell size on a projected grid.
    !       At 45 N the length correction factor is cos(45) = sqrt(2)/2, giving an actual cell length of ~140 m.
    !       The correction is as follows:
    !       (1) Set an average length correction factor, glacier%length_factor, in the config file.
    !           Multiply dew and dns by this factor so the dynamics will see the (approximately) correct length.
    !       (2) Compute a corrected cell_area(i,j) based on the latitude: cell_area -> cell_area * cos^2(lat),
    !           where cos^2(lat) is roughly equal to length_factor^2, but not exactly since lat depends on (i,j).

    ! Set some local variables
    parallel = model%parallel

    global_ewn = parallel%global_ewn
    global_nsn = parallel%global_nsn
    ewn = model%general%ewn
    nsn = model%general%nsn
    dew = model%numerics%dew
    dns = model%numerics%dns
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    if (verbose_glacier) then
       if (this_rank == rtest) write(iulog,*) 'In glissade_glacier_init'
       call point_diag(glacier%rgi_glacier_id, 'RGI glacier ID', itest, jtest, rtest, 7, 7)
    endif

    if (glacier%scale_area) then

       ! Optionally, rescale the grid cell dimensions and coordinates
       ! This is answer-changing throughout the code.
       ! Note: The global arrays model%general%x1_global, etc., which are written to output files, are not rescaled.
       !       These arrays are computed from the input file, which typically ignores the scale factor.
       if (glacier%length_scale_factor /= 1.0d0) then
          model%general%x0 = model%general%x0 * glacier%length_scale_factor
          model%general%y0 = model%general%y0 * glacier%length_scale_factor
          model%general%x1 = model%general%x1 * glacier%length_scale_factor
          model%general%y1 = model%general%y1 * glacier%length_scale_factor
          model%numerics%dew = model%numerics%dew * glacier%length_scale_factor
          model%numerics%dns = model%numerics%dns * glacier%length_scale_factor
          dew = model%numerics%dew
          dns = model%numerics%dns
       endif

       ! Rescale the grid cell areas (diagnostic only; not used for dynamic calculations).
       ! Originally computed as the (unscaled) product dew*dns; scale here by cos^2(lat).
       ! Note: These use the actual cell latitudes, as opposed to acos(length_scale_factor)
       do j = 1, nsn
          do i = 1, ewn
             theta_rad = model%general%lat(i,j) * pi/180.d0
             model%geometry%cell_area(i,j) = model%geometry%cell_area(i,j) * cos(theta_rad)**2
          enddo
       enddo
       call parallel_halo(model%geometry%cell_area, parallel)

       if (verbose_glacier .and. this_rank == rtest) then
          i = itest; j = jtest
          theta_rad = model%general%lat(i,j) * pi/180.d0
          write(iulog,*) 'Scale dew and dns: factor, new dew, dns =', &
               glacier%length_scale_factor, dew, dns
          write(iulog,*) 'Scale cell area: i, j, lat, cos(lat), cell_area =', &
               i, j, model%general%lat(i,j), cos(theta_rad), model%geometry%cell_area(i,j)
       endif

    endif   ! scale_area

    if (model%options%is_restart == NO_RESTART) then

       ! not a restart; initialize everything from the input file

       ! Note: At start-up, arrays in the glacier derived type are allocated with dimension(1),
       !  since nglacier has not yet been computed.
       ! Deallocate here, and reallocate below with dimension(nglacier).
       ! For a restart, nglacier is determined from the restart file,
       !  and these arrays should already have the correct dimensions.

       if (associated(glacier%glacierid)) deallocate(glacier%glacierid)
       if (associated(glacier%cism_to_rgi_glacier_id)) deallocate(glacier%cism_to_rgi_glacier_id)
       if (associated(glacier%area)) deallocate(glacier%area)
       if (associated(glacier%volume)) deallocate(glacier%volume)
       if (associated(glacier%area_init)) deallocate(glacier%area_init)
       if (associated(glacier%volume_init)) deallocate(glacier%volume_init)
       if (associated(glacier%area_init_extent)) deallocate(glacier%area_init_extent)
       if (associated(glacier%volume_init_extent)) deallocate(glacier%volume_init_extent)
       if (associated(glacier%area_target)) deallocate(glacier%area_target)
       if (associated(glacier%volume_target)) deallocate(glacier%volume_target)
       if (associated(glacier%smb)) deallocate(glacier%smb)
       if (associated(glacier%smb_obs)) deallocate(glacier%smb_obs)
       if (associated(glacier%mu_star)) deallocate(glacier%mu_star)
       if (associated(glacier%alpha_snow)) deallocate(glacier%alpha_snow)
       if (associated(glacier%beta_artm)) deallocate(glacier%beta_artm)

       ! Set the RGI ID to 0 in cells without ice.
       ! Typically, any ice-free cell should already have an RGI ID of 0,
       !  but there can be exceptions related to no-ice boundary conditions.
       where (model%geometry%thck == 0.0d0)
          glacier%rgi_glacier_id = 0
       endwhere

       ! Count the number of cells with glaciers
       ! Loop over locally owned cells

       glacier_mask = 0
       where (glacier%rgi_glacier_id > 0) glacier_mask = 1
       ncells_glacier = parallel_global_sum(glacier_mask, parallel)

       ! Gather the RGI glacier IDs to the main task
       if (main_task) allocate(rgi_glacier_id_global(global_ewn, global_nsn))
       call gather_var(glacier%rgi_glacier_id, rgi_glacier_id_global, parallel)

       ! Allocate a global array for the CISM glacier IDs on the main task.
       ! Allocate a size 0 array on other tasks; scatter_var wants arrays allocated on all tasks.
       if (main_task) then
          allocate(cism_glacier_id_global(global_ewn,global_nsn))
       else
          allocate(cism_glacier_id_global(0,0))
       endif
       cism_glacier_id_global(:,:) = 0.0d0

       if (verbose_glacier .and. main_task) then
          write(iulog,*) ' '
          write(iulog,*) 'Gathered RGI glacier IDs to main task'
          write(iulog,*) 'size(rgi_glacier_id) =', &
               size(glacier%rgi_glacier_id,1), size(glacier%rgi_glacier_id,2)
          write(iulog,*) 'size(rgi_glacier_id_global) =', &
               size(rgi_glacier_id_global,1), size(rgi_glacier_id_global,2)
       endif

       if (main_task) then

          gid_minval = minval(rgi_glacier_id_global)
          gid_maxval = maxval(rgi_glacier_id_global)

          if (verbose_glacier) then
             write(iulog,*) 'Total ncells   =', global_ewn * global_nsn
             write(iulog,*) 'ncells_glacier =', ncells_glacier
             write(iulog,*) 'glacier_id minval, maxval =', gid_minval, gid_maxval
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
             write(iulog,*) 'Sorted glacier IDs in ascending order'
             write(iulog,*) ' '
             write(iulog,*) 'icell, i, j, ID for a few cells:'
             do i = 1, 10
                write(iulog,*) i, glacier_list(i)%indxi, glacier_list(i)%indxj, glacier_list(i)%id
             enddo
             do i = ncells_glacier-9, ncells_glacier
                write(iulog,*) i, glacier_list(i)%indxi, glacier_list(i)%indxj, glacier_list(i)%id
             enddo
          endif

!       WHL - Short list to test quicksort for integer arrays
!       write(iulog,*) ' '
!       write(iulog,*) 'Unsorted list:'
!       nlist = 20
!       allocate(test_list(nlist))
!       do i = 1, nlist
!          call random_number(random)
!          test_list(i) = int(random*nlist) + 1
!          write(iulog,*) i, random, test_list(i)
!       enddo
!       call quicksort(test_list, 1, nlist)
!       write(iulog,*) 'Sorted list:', test_list(:)

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
             write(iulog,*) ' '
             write(iulog,*) 'Counted glaciers: nglacier =', nglacier
             write(iulog,*) ' '
             ng = nglacier/2
             write(iulog,*) 'Random cism_glacier_id:', ng
             write(iulog,*) 'icell, i, j, cism_glacier_id_global(i,j), cism_to_rgi_glacier_id(ng)'
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
                write(iulog,*) nc, i, j, cism_glacier_id_global(i,j), glacier%cism_to_rgi_glacier_id(ng)
             endif
             if (ng > nglacier) then
                write(message,*) 'CISM glacier ID > nglacier, i, j , ng =', i, j, ng
                call write_log(message, GM_FATAL)
             endif
          enddo

          deallocate(glacier_list)

          if (verbose_glacier) then
             write(iulog,*) 'maxval(cism_to_rgi_glacier_id) =', maxval(glacier%cism_to_rgi_glacier_id)
             write(iulog,*) 'maxval(cism_glacier_id_global) =', maxval(cism_glacier_id_global)
          endif

       endif   ! main_task

       ! Scatter cism_glacier_id_global to all processors
       ! Note: This global array is deallocated in the scatter_var subroutine
       call scatter_var(glacier%cism_glacier_id, cism_glacier_id_global, parallel)

       call parallel_halo(glacier%cism_glacier_id, parallel)

       ! Copy cism_glacier_id to cism_glacier_id_init
       glacier%cism_glacier_id_init(:,:) = glacier%cism_glacier_id(:,:)

       ! Broadcast nglacier and cism_to_rgi_glacier_id from the main task to all processors
       call broadcast(glacier%nglacier)
       nglacier = glacier%nglacier

       if (.not.associated(glacier%cism_to_rgi_glacier_id)) &
            allocate(glacier%cism_to_rgi_glacier_id(nglacier))
       call broadcast(glacier%cism_to_rgi_glacier_id)

       ! Set each glaciated cell to at least the minimum dynamically active thickness
       ! Adjust the upper surface accordingly
       where (glacier%cism_glacier_id > 0)
          model%geometry%thck = max(model%geometry%thck, model%numerics%thklim)
          model%geometry%usrf = &
               max(model%geometry%usrf, model%geometry%topg + model%geometry%thck)
       endwhere

       ! Allocate glacier arrays with dimension(nglacier).
       ! Note: We should avoid accessing these arrays for grid cells with cism_glacier_id = 0.
       allocate(glacier%glacierid(nglacier))
       allocate(glacier%area(nglacier))
       allocate(glacier%area_init(nglacier))
       allocate(glacier%volume(nglacier))
       allocate(glacier%volume_init(nglacier))
       allocate(glacier%area_init_extent(nglacier))
       allocate(glacier%volume_init_extent(nglacier))
       allocate(glacier%area_target(nglacier))
       allocate(glacier%volume_target(nglacier))
       allocate(glacier%smb(nglacier))
       allocate(glacier%smb_obs(nglacier))
       allocate(glacier%mu_star(nglacier))
       allocate(glacier%alpha_snow(nglacier))
       allocate(glacier%beta_artm(nglacier))

       ! Compute the initial area and volume of each glacier.
       ! These values are saved and written to the restart file.
       ! Only ice thicker than diagnostic_minthck is included in area and volume sums.

       call glacier_area_volume(&
            ewn,           nsn,               &
            parallel,                         &
            nglacier,                         &
            glacier%cism_glacier_id_init,     &
            model%geometry%cell_area,         &  ! m^2
            model%geometry%thck,              &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area_init,                &  ! m^2
            glacier%volume_init)                 ! m^3

       ! Initialize other glacier arrays
       glacier%area(:)   = glacier%area_init(:)
       glacier%volume(:) = glacier%volume_init(:)
       glacier%area_init_extent(:) = glacier%area_init(:)
       glacier%volume_init_extent(:) = glacier%volume_init(:)
       glacier%area_target(:) = glacier%area_init(:)
       glacier%volume_target(:) = glacier%volume_init(:)
       glacier%smb(:)         = 0.0d0
       glacier%mu_star(:)     = glacier%mu_star_const
       glacier%alpha_snow(:)  = glacier%alpha_snow_const
       glacier%beta_artm(:)   = 0.0d0

       ! Initially, allow nonzero SMB only in glacier-covered cells.
       ! These masks are updated at runtime.
       glacier%smb_glacier_id_init(:,:) = glacier%cism_glacier_id_init(:,:)
       glacier%smb_glacier_id(:,:) = glacier%cism_glacier_id_init(:,:)

       ! Check for area_init = 0 and volume_init = 0.
       ! In practice, volume_init = 0 might not be problematic;
       !  we would just lower powerlaw_c to obtain a thin glacier.
       ! Could have area_init = 0 if all the ice in the glacier is thinner
       !  than the diagnostic threshold.

       if (main_task) then
          do ng = 1, nglacier
             if (glacier%area_init(ng) == 0.0d0) then
                write(message,*) 'Glacier area init = 0: ng =', ng
                call write_log(message)
             endif
             if (glacier%volume_init(ng) == 0.0d0) then
                write(message,*) 'Glacier volume init = 0: ng, area (km^2) =', &
                     ng, glacier%area(ng)/1.0d6
                call write_log(message)
             endif
          enddo   ! ng
       endif

       ! Save the initial usrf to usrf_obs.
       ! This value becomes the RGI target and is read on restart
       model%geometry%usrf_obs = model%geometry%usrf

       ! If inverting for powerlaw_c, then initialize powerlaw_c to a constant value,
       !  and initialize the inversion target to the initial thickness.
       ! Note: When inverting for thickness, thck_target is the target for the baseline date,
       !       which usually is earlier than the RGI date. Thus, thck_target usually is greater than
       !       the input thickness, if the input thickness corresponds to the RGI date.
       ! On restart, powerlaw_c is read from the restart file;
       !  thck_target is not a restart field but is updated annually during the inversion.
       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          model%basal_physics%powerlaw_c(:,:) = model%basal_physics%powerlaw_c_const
          glacier%thck_target = model%geometry%thck
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
             ng = glacier%cism_glacier_id_init(i,j)
             if (ng == 0 .and. model%geometry%thck(i,j) > 1.0d0) then
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                write(iulog,*) 'Warning, ng = 0 but H > 0: Init rank, i, j, ig, jg, thck:', &
                     this_rank, i, j, iglobal, jglobal, model%geometry%thck(i,j)
             endif
          enddo
       enddo

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
           glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
          ! Make sure a nonzero smb_obs field was read in
          if (parallel_is_zero(model%climate%smb_obs)) then
             call write_log ('Error, no nonzero values for smb_obs', GM_FATAL)
          endif
       else
          ! If a nonzero smb_obs field was read in, then set to zero
          model%climate%smb_obs = 0.0d0
       endif

       ! Use the 2D smb_obs field to compute the 1D glacier-average field.
       ! On restart, this 1D field will be read from the restart file.

       call glacier_2d_to_1d(&
            ewn,                   nsn,                           &
            parallel,                                             &
            nglacier,              glacier%cism_glacier_id_init,  &
            model%climate%smb_obs, glacier%smb_obs)

    else  ! restart (either standard or hybrid)

       ! In this case, most required glacier info has already been read from the restart file.
       ! Here, do some error checks and diagnostics.

       ! From the restart file, nglacier is found as the length of dimension 'glacierid'.
       ! The 1D glacier arrays are then allocated with dimension(nglacier) in subroutine glide_allocarr.
       ! The following glacier arrays should be present in the restart file:
       !     rgi_glacier_id, cism_glacier_id, cism_glacier_id_init, cism_to_rgi_glacier_id,
       !     glacier_mu_star, and powerlaw_c.
       ! Note: Depending on the model settings, some other fields are needed too.
       !       The code below does not check for all required fields.
       ! If inverting for mu_star and alpha_snow, then usrf_obs and smb_obs should be read from the restart file.
       ! If inverting for mu_star alone, then usrf_obs should be read from the restart file.

       nglacier = glacier%nglacier

       ! Check that some glacier arrays which are read from the restart file have nonzero values.
       ! Note: These arrays are read on all processors.
       ! TODO: Use the parallel_is_zero interface.

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

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION) then
          max_glcval = maxval(model%geometry%usrf_obs)
          max_glcval = parallel_reduce_max(max_glcval)
          if (max_glcval <= 0.0d0) then
             call write_log ('Error, no positive values for usrf_obs', GM_FATAL)
          endif
          if (glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
             ! need nonzero smb_obs for inversion
             max_glcval = maxval(abs(glacier%smb_obs))
             max_glcval = parallel_reduce_max(max_glcval)
             if (max_glcval == 0.d0) then
                call write_log ('Error, no nonzero values for smb_obs', GM_FATAL)
             endif
          else   ! inverting for mu_star only; 1-equation scheme with SMB = 0
             glacier%smb_obs = 0.0d0
          endif
       endif

       ! Compute the area and volume of each glacier (diagnostic only)

       call glacier_area_volume(&
            ewn,           nsn,               &
            parallel,                         &
            nglacier,                         &
            glacier%cism_glacier_id,          &
            model%geometry%cell_area,         &  ! m^2
            model%geometry%thck,              &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area,                     &  ! m^2
            glacier%volume)                      ! m^3

       ! Repeat, summing over the initial glacier extent

       call glacier_area_volume(&
            ewn,           nsn,               &
            parallel,                         &
            nglacier,                         &
            glacier%cism_glacier_id_init,     &
            model%geometry%cell_area,         &  ! m^2
            model%geometry%thck,              &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area_init_extent,         &  ! m^2
            glacier%volume_init_extent)          ! m^3

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

    glacier%minthck = model%numerics%thklim - eps08

    ! Set the index of the diagnostic glacier, using the CISM glacier ID for the diagnostic point
    if (this_rank == rtest) then
       glacier%ngdiag = glacier%cism_glacier_id_init(itest,jtest)
       if (glacier%ngdiag == 0) then
          write(message,*) &
               'The diagnostic cell has cism_glacier_id = 0; may want to choose a different cell'
          call write_log(message, GM_WARNING)
       endif
    endif
    call broadcast(glacier%ngdiag, rtest)

    ! Define a mask whose value is 1 at vertices that border two different glaciers.
    ! At runtime, Cp is set to a large value at these vertices to reduce mass exchange between glaciers.
    !TODO: Consider removing the mask. This would allow CISM to reduce basal friction to thin the ice if needed.
    glacier%boundary_mask(:,:) = 0

    ! Loop over locally owned vertices
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng_ne = glacier%cism_glacier_id_init(i+1,j+1)
          ng_nw = glacier%cism_glacier_id_init(i,j+1)
          ng_se = glacier%cism_glacier_id_init(i+1,j)
          ng_sw = glacier%cism_glacier_id_init(i,j)
          if ( (ng_ne > 0 .and. ng_nw > 0 .and. ng_ne /= ng_nw) .or. &
               (ng_ne > 0 .and. ng_se > 0 .and. ng_ne /= ng_se) .or. &
               (ng_ne > 0 .and. ng_sw > 0 .and. ng_ne /= ng_sw) .or. &
               (ng_nw > 0 .and. ng_se > 0 .and. ng_nw /= ng_se) .or. &
               (ng_nw > 0 .and. ng_sw > 0 .and. ng_nw /= ng_sw) .or. &
               (ng_se > 0 .and. ng_sw > 0 .and. ng_se /= ng_sw) ) then
             glacier%boundary_mask(i,j) = 1
          endif
       enddo
    enddo

    call staggered_parallel_halo(glacier%boundary_mask, parallel)

    if (verbose_glacier) then
       call point_diag(glacier%cism_glacier_id_init, 'cism_glacier_id_init', itest, jtest, rtest, 7, 7)
       call point_diag(glacier%boundary_mask, 'Glacier boundary mask', itest, jtest, rtest, 7, 7)
    endif

    ! Write some values for the diagnostic glacier
    if (verbose_glacier .and. this_rank == rtest) then
       i = itest; j = jtest
       ng = glacier%ngdiag
       write(iulog,*) ' '
       write(iulog,*) 'Glacier ID for diagnostic cell: r, i, j, ng =', rtest, itest, jtest, ng
       if (ng > 0) then
          write(iulog,*) 'area_init (km^2) =', glacier%area_init(ng) / 1.0d6
          write(iulog,*) 'volume_init (km^3) =', glacier%volume_init(ng) / 1.0d9
          write(iulog,*) 'powerlaw_c (Pa (m/yr)^(-1/3)) =', model%basal_physics%powerlaw_c(i,j)
          write(iulog,*) 'smb_obs (mm/yr w.e.) =', glacier%smb_obs(ng)
          write(iulog,*) 'mu_star (mm/yr w.e./deg) =', glacier%mu_star(ng)
          write(iulog,*) 'Done in glissade_glacier_init'
       endif
    endif

  end subroutine glissade_glacier_init

!****************************************************

  subroutine glissade_glacier_update(model, glacier)

    use glissade_grid_operators, only: glissade_stagger
    use glissade_utils, only: glissade_usrf_to_thck
    use cism_parallel, only: parallel_global_sum, &
         parallel_halo, staggered_parallel_halo

    ! Do glacier inversion (if applicable), update glacier masks, and compute glacier diagnostics.
    !
    ! SMB is computed from an empirical relationship based on Maussion et al. (2019):
    !
    !     SMB = alpha_snow * snow - mu_star * max(artm - tmlt, 0),
    !
    ! where snow = monthly mean snowfall rate (mm/yr w.e.),
    !       alpha_snow is a glacier-specific tuning parameter (a scalar of order 1)
    !       mu_star is a glacier-specific tuning parameter (mm/yr w.e./deg C),
    !       atrm = monthly mean air temperature (deg C),
    !       tmlt = monthly mean air temp above which ablation occurs (deg C)

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
         ice_mask,                & ! = 1 where ice is present (thck > thklim), else = 0
         glacier_mask               ! temporary mask

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck,                    & ! ice thickness (m)
         dthck_dt,                & ! rate of change of thickness (m/yr)
         cell_area_uniform,       & ! grid cell area defined as dew*dns(m^2)
         thck_old,                & ! saved value of ice thickness (m)
         artm,                    & ! artm, baseline or current date
         snow,                    & ! snowfall, baseline or current date
         precip,                  & ! precip, baseline or current date
         Tpos,                    & ! max(artm - tmlt, 0.0)
         artm_recent,             & ! artm, recent (smb_obs) date
         snow_recent,             & ! snowfall rate (mm w.e./yr), recent date
         precip_recent,           & ! precip rate, recent date
         Tpos_recent,             & ! max(artm - tmlt, 0.0), recent date
         artm_rgi,                & ! artm, RGI date
         snow_rgi,                & ! snowfall rate, RGI date
         precip_rgi,              & ! precip rate, RGI date
         Tpos_rgi,                & ! max(artm - tmlt, 0.0), RGI date
         mu_star_2d,              & ! 2D version of glacier%mu_star
         alpha_snow_2d,           & ! 2D version of glacier%alpha_snow
         delta_smb_rgi,           & ! SMB anomaly between the baseline date and the RGI date (mm/yr w.e.)
         delta_smb_recent,        & ! SMB anomaly between the baseline date and the recent date (mm/yr w.e.)
         smb_weight_init,         & ! ratio of applied SMB to potential SMB, in range [0,1], for sums over initial area
         smb_weight_current         ! ratio of applied SMB to potential SMB, in range [0,1], for sums over current area

    real(dp), dimension(model%general%ewn-1, model%general%nsn-1) ::  &
         stag_thck,                 & ! ice thickness at vertices (m)
         stag_thck_target,          & ! target ice thickness at vertices (m)
         stag_dthck_dt                ! rate of change of ice thickness at vertices (m/yr)

    type(parallel_type) :: parallel   ! info for parallel communication

    real(dp), save ::  &              ! time since the last averaging computation (yr);
         time_since_last_avg = 0.0d0  ! compute the average once a year

    real(dp), dimension(glacier%nglacier) :: &
         area_old,                  & ! glacier%area from the previous inversion step
         darea_dt,                  & ! rate of change of glacier area over the inversion interval
         smb_init_area,             & ! SMB over initial area determined by cism_glacier_id_init
         smb_current_area,          & ! SMB over current area determined by cism_glacier_id
         smb_min, smb_max,          & ! min and max SMB for each glacier (mm/yr w.e.)
         smb_min_recent,            & ! min and max SMB for each glacier in recent climate (mm/yr w.e.)
         smb_max_recent,            & !
         aar_init, aar,             & ! accumulation area ratio for baseline climate (init and current area)
         aar_init_recent, aar_recent  ! accumulation area ratio for recent climate (init and current area)

    ! Note: The glacier type includes the following:
    ! integer ::  nglacier            ! number of glaciers in the global domain
    ! integer ::  ngdiag              ! CISM index of diagnostic glacier
    ! real(dp), dimension(:) :: area              ! glacier area (m^2)
    ! real(dp), dimension(:) :: volume            ! glacier volume (m^3)
    ! real(dp), dimension(:) :: area_init         ! initial glacier area (m^2)
    ! real(dp), dimension(:) :: volume_init       ! initial glacier volume (m^3)
    ! real(dp), dimension(:) :: area_init_extent  ! current glacier area (m^2) over initial ice extent
    ! real(dp), dimension(:) :: volume_init_extent! current glacier volume (m^3) over initial ice extent
    ! real(dp), dimension(:) :: area_target       ! target glacier area (m^2) for inversion
    ! real(dp), dimension(:) :: volume_target     ! target glacier volume (m^3) for inversion
    ! real(dp), dimension(:) :: mu_star           ! SMB parameter for each glacier (mm/yr w.e./deg K)
    ! real(dp), dimension(:) :: alpha_snow        ! snow factor for each glacier (unitless)
    ! real(dp), dimension(:) :: beta_artm         ! artm correction for each glacier (deg C)
    ! real(dp), dimension(:) :: smb_obs           ! observed SMB for each glacier (mm/yr w.e.)
    ! integer, dimension(:,:) :: cism_glacier_id       ! CISM glacier ID for each grid cell
    ! integer, dimension(:,:) :: cism_glacier_id_init  ! initial value of CISM glacier ID
    ! integer, dimension(:,:) :: smb_glacier_id        ! CISM glacier ID that determines where SMB is applied
    ! integer, dimension(:,:) :: smb_glacier_id_init   ! like smb_glacier_id, but based on cism_glacier_id_init
    ! real(dp), dimension(:,:) :: snow_annmean         ! snow accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: Tpos_annmean         ! max(artm - tmlt,0) accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: snow_recent_annmean  ! snow accumulated and averaged over 1 year, recent date
    ! real(dp), dimension(:,:) :: Tpos_recent_annmean  ! max(artm - tmlt,0) accumulated and averaged over 1 year, recent date
    ! real(dp), dimension(:,:) :: snow_rgi_annmean     ! snow accumulated and averaged over 1 year, RGI date
    ! real(dp), dimension(:,:) :: Tpos_rgi_annmean     ! max(artm - tmlt,0) accumulated and averaged over 1 year, RGI date
    ! real(dp), dimension(:,:) :: dthck_dt_annmean     ! dthck_dt accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: usrf_target          ! target surface elevation (m) for the baseline climate
    ! real(dp), dimension(:,:) :: thck_target          ! target thickness (m) for the baseline climate
    !TODO - Are any glacier fields missing?

    ! Note: The following areas are computed based on the cism_glacier_id masks, without a min thickness criterion
    real(dp), dimension(glacier%nglacier) ::  &
         area_initial, area_current,    &  ! initial and current glacier areas (m^2)
         area_advance, area_retreat        ! areas of glacier advance and retreat relative to initial mask (m^2)

    real(dp) :: area_sum
    real(dp) :: usrf_recent    ! estimated surface elevation in recent climate
    real(dp) :: usrf_rgi       ! estimated surface elevation in RGI climate
    real(dp) :: rgi_date_frac
    real(dp), parameter :: diagnostic_volume_threshold = 1.0d9  ! volume threshold for big glaciers (m^3)

    integer :: count_cgii, count_cgi
    integer :: count_sgii, count_sgi

    ! Set some local variables

    parallel = model%parallel

    ewn = model%general%ewn
    nsn = model%general%nsn
    dew = model%numerics%dew                  ! convert to m
    dns = model%numerics%dns                  ! convert to m
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    nglacier = glacier%nglacier
    ngdiag = glacier%ngdiag
    cell_area_uniform = dew*dns

    ! some unit conversions
    !       Skip these conversion and use SI units (s instead of yr) in the code.
    dt = model%numerics%dt /scyr                    ! s to yr
    dthck_dt = model%geometry%dthck_dt * scyr       ! m/s to m/yr

    ! Accumulate the 2D fields used for mu_star and alpha_snow inversion: snow and Tpos.
    ! Also accumulate dthck_dt, which is used for powerlaw_c inversion.
    ! Note: snow and Tpos are also used by subroutines glacier_advance_retreat
    !       and update_smb_glacier_id. Thus, they are accumulated and updated
    !       during forward runs with fixed mu_star and alpha_snow, not just
    !       spin-ups with inversion for mu_star and alpha_snow.

    if (time_since_last_avg == 0.0d0) then ! start of new averaging period

       ! Reset the accumulated fields to zero
       glacier%snow_annmean = 0.0d0
       glacier%Tpos_annmean = 0.0d0
       glacier%smb_applied_annmean = 0.0d0

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
           glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
          glacier%snow_recent_annmean = 0.0d0
          glacier%Tpos_recent_annmean = 0.0d0
          glacier%snow_rgi_annmean = 0.0d0
          glacier%Tpos_rgi_annmean = 0.0d0
       endif

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          glacier%dthck_dt_annmean = 0.0d0
       endif

       ! If inverting for mu_star and alpha_snow, then compute some SMB-related quantities
       ! used in the inversion.

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
           glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then

          ! Compute the SMB anomaly for the RGI and recent climates relative to the baseline climate.
          ! This is done once a year; smb, smb_rgi, and smb_recent are updated at the end of the previous year.

          where (glacier%cism_glacier_id_init > 0)
             delta_smb_rgi = glacier%smb_rgi - model%climate%smb
             glacier%delta_usrf_rgi = delta_smb_rgi*(rhow/rhoi)/1000.d0 * &
               (glacier%rgi_date - glacier%baseline_date)/2.d0
             delta_smb_recent = glacier%smb_recent - model%climate%smb
             glacier%delta_usrf_recent = delta_smb_recent*(rhow/rhoi)/1000.d0 * &
               (glacier%recent_date - glacier%baseline_date)/2.0d0  ! m ice
          elsewhere
             delta_smb_rgi = 0.0d0
             delta_smb_recent = 0.0d0
          endwhere

          ! Adjust the baseline target. The baseline target should exceed the RGI target by abs(delta_usrf_rgi),
          !  assuming the ice thins between the baseline and RGI dates.
          ! Then, provided usrf is close to usrf_target in the spin-up, usrf will be close to
          !  usrf_obs (the RGI target) when a forward run starting from the baseline date reaches the RGI date.
          !TODO - How to set usrf_target if not inverting for mu_star? Set to usrf_obs?

          glacier%usrf_target(:,:) = model%geometry%usrf_obs(:,:) - glacier%delta_usrf_rgi(:,:)

          ! Make sure the target is not below the topography
          glacier%usrf_target = &
               max(glacier%usrf_target, (model%geometry%topg + model%climate%eus))

          if (verbose_glacier .and. this_rank == rtest) then
             i = itest; j = jtest
             write(iulog,*) ' '
             write(iulog,*) 'RGI usrf correction, delta_smb:', &
                  glacier%delta_usrf_rgi(i,j), delta_smb_rgi(i,j)
             write(iulog,*)    'usrf RGI obs, new usrf_target baseline =', &
                  model%geometry%usrf_obs(i,j), glacier%usrf_target(i,j)
             write(iulog,*) 'Recent usrf correction, delta_smb:', &
                  glacier%delta_usrf_recent(i,j), delta_smb_recent(i,j)
          endif

       endif   ! set_mu_star

    endif   ! time_since_last_avg = 0

    ! Halo updates for snow and artm
    ! Note: artm_corrected, snow_corrected, and precip_corrected are the input fields.
    !       The 'corrected' suffix means that anomaly forcing, if enabled, has been included.
    !       Assuming artm_input_function = xy_lapse, a lapse rate correction has already been applied.
    ! Note: snow_calc is the snow calculation option: Either use the snowfall rate directly,
    !       or compute the snowfall rate from the precip rate and downscaled artm.
    !TODO - Not sure these are needed. Maybe can save halo updates for the annual-averaged snow and Tpos

    call parallel_halo(model%climate%artm_corrected, parallel)
    if (glacier%snow_calc == GLACIER_SNOW_CALC_SNOW) then
       call parallel_halo(model%climate%snow_corrected, parallel)
    elseif (glacier%snow_calc == GLACIER_SNOW_CALC_PRECIP_ARTM) then
       call parallel_halo(model%climate%precip_corrected, parallel)
    endif

    ! Initialize the glacier fields: artm, snow, and precip.
    ! If inverting for mu_star, then artm, snow, and precip apply to the baseline climate.
    ! For forward runs, artm and Tpos apply to the current climate.
    !
    ! The 'corrected' suffix means that anomaly forcing, if enabled, has already been included.
    ! When inverting for mu_star, the anomaly fields are used to form the 'recent' forcing fields below,
    !  but are not part of the baseline climate fields.
    !  We have enable_acab_anomaly = enable_snow_anomaly = enable_snow_anomaly = F,
    !  and thus the anomaly fields are ignored in glissade.F90.
    ! To include anomaly forcing in forward runs, we set enable_acab_anomaly = enable_snow_anomaly
    !  = enable_snow_anomaly = T. Then the anomaly fields are added to the baseline fields in glissade.F90
    !  to form the current fields.

    artm(:,:)   = model%climate%artm_corrected(:,:)
    snow(:,:)   = model%climate%snow_corrected(:,:)
    precip(:,:) = model%climate%precip_corrected(:,:)

    ! Add the beta temperature correction term for glaciers with nonzero beta_artm.
    ! Note: smb_glacier_id = smb_glacier_id_init wherever smb_glacier_id_init > 0

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = glacier%smb_glacier_id(i,j)
          if (ng > 0) then
             artm(i,j) = artm(i,j) + glacier%beta_artm(ng)
          endif
          Tpos(i,j) = max(artm(i,j) - glacier%tmlt, 0.0d0)
       enddo
    enddo

    ! Compute the snowfall rate

    if (glacier%snow_calc == GLACIER_SNOW_CALC_SNOW) then

       ! do nothing; use the input snowfall rate directly

    elseif (glacier%snow_calc == GLACIER_SNOW_CALC_PRECIP_ARTM) then

       ! compute snowfall based on precip and artm

       call glacier_calc_snow(&
            ewn,       nsn,                   &
            glacier%snow_threshold_min,       &
            glacier%snow_threshold_max,       &
            precip,                           &
            artm,                             &
            snow)

    endif   ! snow_calc

    if (verbose_glacier .and. this_rank == rtest) then
       i = itest; j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'glissade_glacier_update, diag cell (r, i, j) =', rtest, itest, jtest
       write(iulog,*) ' '
       ! Convert acab_applied from m/yr ice to mm/yr w.e.
       write(iulog,'(a32,2f10.3)') '     acab_applied, smb_applied: ', &
            model%climate%acab_applied(i,j)*scyr, &  ! m/yr ice
            model%climate%acab_applied(i,j)*scyr * 1000.d0*(rhoi/rhow)  ! mm/yr w.e.
       write(iulog,'(a32,4f10.3)') 'artm_ref, usrf_ref, usrf, diff: ', &
            model%climate%artm_ref(i,j), &
            model%climate%usrf_ref(i,j), model%geometry%usrf(i,j), &
            model%geometry%usrf(i,j) - model%climate%usrf_ref(i,j)
       write(iulog,'(a32,3f10.3)') '              artm, Tpos, snow: ', artm(i,j), Tpos(i,j), snow(i,j)
    endif   ! verbose

    ! If inverting for mu and alpha, then compute artm_ref, snow, and precip at the recent and RGI dates.
    ! Note: When inverting for mu_star and alpha, we have enable_artm_anomaly = enable_snow_anomaly =
    !       enable_precip_anomaly = F. The anomalies are used here for inversion, but are not applied
    !       in the main glissade module.

    if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
        glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then

       artm_recent(:,:)   = artm(:,:)   + model%climate%artm_anomaly(:,:)
       snow_recent(:,:)   = snow(:,:)   + model%climate%snow_anomaly(:,:)
       precip_recent(:,:) = precip(:,:) + model%climate%precip_anomaly(:,:)

       ! Compute artm and Tpos for the recent climate at the extrapolated surface elevation.
       ! We estimate usrf_recent = usrf + (dSMB/2)*dt,
       !    where dSMB = smb_recent - smb is the difference in SMB between the baseline and recent climate,
       !          (so dSMB/2 is the average SMB anomaly over that period), and dt is the number of years elapsed.
       ! In other words, assume that the entire SMB anomaly is used to melt ice, without the
       !  flow having time to adjust.

       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             artm_recent(i,j) = artm_recent(i,j) - glacier%delta_usrf_recent(i,j)*model%climate%t_lapse
             Tpos_recent(i,j) = max(artm_recent(i,j) - glacier%tmlt, 0.0d0)
          enddo
       enddo

       ! Estimate artm, Tpos, and snow or precip for the RGI climate by interpolation.

       rgi_date_frac = (glacier%rgi_date    - glacier%baseline_date) / &
                       (glacier%recent_date - glacier%baseline_date)

       artm_rgi(:,:) = &
            (1.d0 - rgi_date_frac) * artm(:,:)  &
                 +  rgi_date_frac  * artm_recent(:,:)

       Tpos_rgi(:,:) = max(artm_rgi(:,:) - glacier%tmlt, 0.0d0)

       ! Compute the snowfall rate for the RGI and recent climate.

       if (glacier%snow_calc == GLACIER_SNOW_CALC_SNOW) then

          snow_rgi(:,:) = (1.d0 - rgi_date_frac) * snow(:,:)  &
                               +  rgi_date_frac  * snow_recent(:,:)

       elseif (glacier%snow_calc == GLACIER_SNOW_CALC_PRECIP_ARTM) then

          call glacier_calc_snow(&
               ewn,       nsn,                   &
               glacier%snow_threshold_min,       &
               glacier%snow_threshold_max,       &
               precip_recent,                    &
               artm_recent,                      &
               snow_recent)

          precip_rgi(:,:) = (1.d0 - rgi_date_frac) * precip(:,:)  &
                                  + rgi_date_frac  * precip_recent(:,:)

          call glacier_calc_snow(&
               ewn,       nsn,                   &
               glacier%snow_threshold_min,       &
               glacier%snow_threshold_max,       &
               precip_rgi,                       &
               artm_rgi,                         &
               snow_rgi)

       endif   ! snow calc

       if (verbose_glacier .and. this_rank == rtest) then
          i = itest; j = jtest
          write(iulog,'(a32,3f10.3)') '         RGI artm, Tpos, snow: ', &
               artm_rgi(i,j), Tpos_rgi(i,j), snow_rgi(i,j)
          write(iulog,'(a32,3f10.3)') '      Recent artm, Tpos, snow: ', &
               artm_recent(i,j), Tpos_recent(i,j), snow_recent(i,j)
       endif

    endif   ! set_mu_star

    ! Accumulate snow_annmean, Tpos_annmean, and dthck_dt_annmean over this timestep

    time_since_last_avg = time_since_last_avg + dt

    glacier%snow_annmean = glacier%snow_annmean + snow * dt
    glacier%Tpos_annmean = glacier%Tpos_annmean + Tpos * dt
    glacier%smb_applied_annmean = glacier%smb_applied_annmean  &
         + model%climate%acab_applied*scyr * 1000.d0*(rhoi/rhow) * dt

    if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
        glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
       glacier%snow_rgi_annmean = glacier%snow_rgi_annmean + snow_rgi * dt
       glacier%Tpos_rgi_annmean = glacier%Tpos_rgi_annmean + Tpos_rgi * dt
       glacier%snow_recent_annmean = glacier%snow_recent_annmean + snow_recent * dt
       glacier%Tpos_recent_annmean = glacier%Tpos_recent_annmean + Tpos_recent * dt
    endif

    if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
       glacier%dthck_dt_annmean = glacier%dthck_dt_annmean + dthck_dt * dt
    endif

    ! Check whether it is time to do the inversion and update other glacier fields.
    ! Note: time_since_last_avg is real(dp) with units of yr;
    !       glacier_update_interval is an integer number of years.

    if (abs(time_since_last_avg - real(glacier_update_interval,dp)) < eps08) then

       ! Average the glacier fields over the accumulation period

       glacier%snow_annmean = glacier%snow_annmean / time_since_last_avg
       glacier%Tpos_annmean = glacier%Tpos_annmean / time_since_last_avg
       glacier%smb_applied_annmean = glacier%smb_applied_annmean / time_since_last_avg

    if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
        glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
          glacier%snow_rgi_annmean = glacier%snow_rgi_annmean / time_since_last_avg
          glacier%Tpos_rgi_annmean = glacier%Tpos_rgi_annmean / time_since_last_avg
          glacier%snow_recent_annmean = glacier%snow_recent_annmean / time_since_last_avg
          glacier%Tpos_recent_annmean = glacier%Tpos_recent_annmean / time_since_last_avg
       endif

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          glacier%dthck_dt_annmean = glacier%dthck_dt_annmean / time_since_last_avg
       endif

       time_since_last_avg = 0.0d0

       if (verbose_glacier .and. this_rank == rtest) then
          i = itest; j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'Annual averages, r, i, j:', rtest, itest, jtest
          write(iulog,*) '   snow (mm/yr)       =', glacier%snow_annmean(i,j)
          write(iulog,*) '   Tpos (deg C)       =', glacier%Tpos_annmean(i,j)
          write(iulog,*) '   smb_applied (mm/yr)=', glacier%smb_applied_annmean(i,j)
          if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
              glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
             write(iulog,*) '   snow_rgi (mm/yr)   =', glacier%snow_rgi_annmean(i,j)
             write(iulog,*) '   Tpos_rgi (deg C)   =', glacier%Tpos_rgi_annmean(i,j)
             write(iulog,*) '   snow_rec (mm/yr)   =', glacier%snow_recent_annmean(i,j)
             write(iulog,*) '   Tpos_rec (deg C)   =', glacier%Tpos_recent_annmean(i,j)
          endif
          if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
             write(iulog,*) '   dthck_dt (m/yr)    =', glacier%dthck_dt_annmean(i,j)
          endif
       endif

       !-------------------------------------------------------------------------
       ! Optionally, thin advanced ice in the accumulation zone to reduce spurious advance.
       ! Ice mass is redistibuted conservatively across the glacier.
       ! Note: Redistribution contributes a positive dH/dt term over the initial extent
       !        and a negative dH/dt term outside the initial extent.
       !       Need to include this contribution in dthck_dt_annmean.
       !-------------------------------------------------------------------------

       if (glacier%redistribute_advanced_ice) then

          thck_old = model%geometry%thck

          call glacier_redistribute_advanced_ice(&
            ewn,             nsn,               &
            parallel,                           &
            itest,   jtest,  rtest,             &
            nglacier,        ngdiag,            &
            real(glacier_update_interval,dp),   &  ! yr
            cell_area_uniform,                  &  ! m^2
            glacier%thinning_rate_advanced_ice, &  ! m/yr
            glacier%cism_glacier_id_init,       &
            glacier%smb_glacier_id,             &
            model%climate%smb,                  &  ! m/yr
            model%geometry%thck)                   ! m

          glacier%dthck_dt_annmean = glacier%dthck_dt_annmean + &
               (model%geometry%thck - thck_old) / real(glacier_update_interval,dp)

       endif   ! redistribute advanced ice

       ! Compute an SMB weighting factor for the inversion.
       ! Set nonzero weights for (1) initial glacier cells and (2) advanced cells in the ablation zone.
       ! Note: For advanced cells in the ablation zone, a weight of zero tends to drive spurious retreat,
       !        while a weight of 1 can allow spurious advance.
       !       An intermediate value of ~0.5 seems to work well.

       smb_weight_init(:,:) = 0.0d0

       where (glacier%cism_glacier_id_init > 0)   ! initial extent
          smb_weight_init = 1.0d0
       elsewhere (glacier%smb_glacier_id_init > 0 .and. model%climate%smb < 0.0d0)
          smb_weight_init = glacier%smb_weight_advanced_ice
       endwhere

       ! Compute the average SMB applied over the initial area of each glacier in the year just finished.
       ! During inversion for mu_star, this should be close to 0 by design.
       ! During a forward run in a warm climate, it will be negative.
       !TODO - Rename smb_init_area?

       call glacier_2d_to_1d_weighted(&
            ewn,                  nsn,       &
            parallel,                        &
            nglacier,                        &
            glacier%smb_glacier_id_init,     &
            smb_weight_init,                 &
            model%climate%smb,    smb_init_area)

       ! Repeat for the current area
       ! Note: Cells in the ablation zone where the full SMB is not applied are given partial weights.
       !       This makes the computed total SMB closer to the true SMB.
       !TODO - Compare use of smb_applied/smb to a constant smb_weight_advanced_ice
       smb_weight_current(:,:) = 0.0d0

       where (glacier%cism_glacier_id > 0)  ! current glacier cells
          smb_weight_current = 1.0d0
       elsewhere (glacier%smb_glacier_id > 0 .and. model%climate%smb < 0.0d0)
          smb_weight_current = glacier%smb_applied_annmean / model%climate%smb
       endwhere

       call glacier_2d_to_1d_weighted(&
            ewn,                  nsn,       &
            parallel,                        &
            nglacier,                        &
            glacier%smb_glacier_id,          &
            smb_weight_current,              &
            model%climate%smb,    smb_current_area)


       ! Invert for mu_star
       ! This can be done in either of two ways:
       ! (1) set_mu_star = 1, set_alpha_snow = 0 (1-parameter inversion)
       !     In this case, mu_star is chosen such that SMB ~ 0 over the initial glacier extent,
       !     given the input temperature and snow/precip fields (without the 'recent' suffix).
       ! (2) set_mu_star = 1, set_alpha_snow = 1 (2-parameter inversion)
       !     In this case, mu_star and alpha_snow are chosen jointly such that
       !     (a) SMB = 0 over the initial extent given the baseline temperature and snow/precip, and
       !     (b) SMB = smb_obs given the recent temperature and snow/precip.
       ! The code aborts at startup if set to invert for alpha_snow without inverting for mu_star.

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION) then

          if (glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then

             ! invert for both mu_star and alpha_snow, based on two SMB conditions
             ! (SMB = 0 in a balanced climate, SMB = smb_obs in an out-of-balance climate)
             ! Note: glacier%smb_obs, glacier%mu_star, and glacier%alpha_snow are 1D glacier-specific fields.

             call glacier_invert_mu_star_alpha_snow(&
                  ewn,                         nsn,                   &
                  parallel,                                           &
                  itest,       jtest,          rtest,                 &
                  nglacier,                    ngdiag,                &
                  glacier%smb_glacier_id_init,                        &
                  smb_weight_init,                                    &
                  glacier%smb_obs,                                    &
                  glacier%area_init,           glacier%volume_init,   &  ! diagnostic only
                  glacier%snow_annmean,        glacier%Tpos_annmean,       &
                  glacier%snow_recent_annmean, glacier%Tpos_recent_annmean,&
                  glacier%mu_star_const,                                   &
                  glacier%mu_star_min,         glacier%mu_star_max,        &
                  glacier%alpha_snow_const,                                &
                  glacier%alpha_snow_min,      glacier%alpha_snow_max,     &
                  glacier%beta_artm_max,       glacier%beta_artm_increment,&
                  glacier%mu_star,             glacier%alpha_snow,         &
                  glacier%beta_artm)

          else  ! not inverting for alpha_snow

             ! Invert for mu_star based on the condition SMB = 0 over the initial glacier extent,
             ! using the default value of alpha_snow (typically 1.0)
             !TODO - Make sure weights are handled OK

             call glacier_invert_mu_star(&
                  ewn,                   nsn,                  &
                  parallel,                                    &
                  itest,     jtest,      rtest,                &
                  nglacier,              ngdiag,               &
                  glacier%smb_glacier_id_init,                 &
                  smb_weight_init,                             &
                  glacier%area_init,     glacier%volume_init,  &  ! diagnostic only
                  glacier%snow_annmean,  glacier%Tpos_annmean, &
                  glacier%mu_star_const,                       &
                  glacier%mu_star_min,   glacier%mu_star_max,  &
                  glacier%beta_artm_max,                       &
                  glacier%beta_artm_increment,                 &
                  glacier%alpha_snow,                          &
                  glacier%mu_star,       glacier%beta_artm)

          endif  ! set_alpha_snow

       endif   ! set_mu_star

       ! advance/retreat diagnostics
       call glacier_area_advance_retreat(&
            ewn,           nsn,           &
            parallel,                     &
            nglacier,                     &
            glacier%cism_glacier_id_init, &
            glacier%cism_glacier_id,      &
            cell_area_uniform,            &
            area_initial,                 &
            area_current,                 &
            area_advance,                 &
            area_retreat)

       if (verbose_glacier .and. this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'Selected big glaciers:'
          write(iulog,'(a101)') &
               '  ng,   Ainit,     A,     Vinit,     V,   smb_iniA, smb_curA, mu_star, alpha_snow, beta_artm, smb_obs'
          do ng = 1, nglacier
             if (glacier%volume_init(ng) > diagnostic_volume_threshold .or. ng == ngdiag) then  ! big glacier
                write(iulog,'(i6,4f9.3,6f10.3)') ng, glacier%area_init(ng)/1.e6, glacier%area(ng)/1.e6, &
                     glacier%volume_init(ng)/1.0d9, glacier%volume(ng)/1.0d9, &
                     smb_init_area(ng), smb_current_area(ng), glacier%mu_star(ng), glacier%alpha_snow(ng), &
                     glacier%beta_artm(ng), glacier%smb_obs(ng)
             endif
          enddo
       endif

       if (verbose_glacier .and. this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'Advance/retreat diagnostics'
          write(iulog,*) '  ng  A_initial A_advance A_retreat A_current'
          do ng = 1, nglacier
             if (glacier%volume_init(ng) > 1.0d9 .or. ng == ngdiag) then  ! big glacier, > 1 km^3
                write(iulog,'(i6,6f10.3)') ng, area_initial(ng)/1.e6, area_advance(ng)/1.e6, &
                     area_retreat(ng)/1.e6, area_current(ng)/1.e6
             endif
          enddo
       endif   ! verbose_glacier

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then

          ! Given the current and target ice thickness, invert for powerlaw_c.
          ! For this to work, the SMB should be close to zero over the initial glacier footprint,
          !  to minimize thickness changes caused by the glacier being out of balance with climate.
          ! This means we must also be inverting for mu_star (and possibly also alpha_snow).
          ! The code aborts at startup if set to invert for powerlaw_c without inverting for mu_star.

          ! Given the surface elevation target, compute the thickness target.
          ! (This can change in time if the bed topography is dynamic.)

          call glissade_usrf_to_thck(&
               glacier%usrf_target,             &
               model%geometry%topg,             &
               model%climate%eus,               &
               glacier%thck_target)

          ! Interpolate thck_target to the staggered grid
          call glissade_stagger(&
               ewn,         nsn,              &
               glacier%thck_target,           &
               stag_thck_target)

          ! Interpolate thck to the staggered grid
          call glissade_stagger(&
               ewn,         nsn,              &
               model%geometry%thck,           &
               stag_thck)

          ! Interpolate dthck_dt to the staggered grid
          call glissade_stagger(&
               ewn,                      nsn,           &
               glacier%dthck_dt_annmean, stag_dthck_dt)

          ! Set stag_thck_dt = 0 at vertices that are initially ice-free.
          ! This will zero out the dH/dt term in the inversion, which inhibits oscillations
          !  in Cp and H near the terminus.
          do j = nhalo, nsn-1
             do i = nhalo, ewn-1
                if (glacier%cism_glacier_id_init(i,  j+1) == 0 .and. &
                    glacier%cism_glacier_id_init(i+1,j+1) == 0 .and. &
                    glacier%cism_glacier_id_init(i,  j)   == 0 .and. &
                    glacier%cism_glacier_id_init(i+1,j)   == 0) then
                   stag_dthck_dt(i,j) = 0.0d0
                endif
             enddo
          enddo
          call staggered_parallel_halo(stag_dthck_dt, parallel)

          ! Update powerlaw_c
          call glacier_invert_powerlaw_c(&
               ewn,                nsn,                 &
               itest,    jtest,    rtest,               &
               model%basal_physics%powerlaw_c_min,      &
               model%basal_physics%powerlaw_c_max,      &
               model%inversion%babc_timescale/scyr,     &  ! yr
               model%inversion%babc_thck_scale,         &  ! m
               model%inversion%babc_relax_factor,       &
               stag_thck,          stag_thck_target,    &
               stag_dthck_dt,                           &
               model%basal_physics%powerlaw_c_const,    &  ! relax to this value
               model%basal_physics%powerlaw_c)

          ! Set Cp to a large value at glacier boundaries, to minimize flow from one glacier to another.
          ! Flow between glaciers is often the result of failing to resolve the surface topography
          ! (e.g., a narrow ridge between two glaciers). A large Cp then substitutes for a physical barrier.
          where (glacier%boundary_mask == 1)
             model%basal_physics%powerlaw_c = model%basal_physics%powerlaw_c_max
          endwhere

       endif   ! set_powerlaw_c

       !-------------------------------------------------------------------------
       ! Update glacier IDs based on advance and retreat since the last update.
       !-------------------------------------------------------------------------

       if (verbose_glacier) then
          call point_diag(model%geometry%thck, 'Before advance_retreat, thck', itest, jtest, rtest, 7, 7)
       endif

       ! Assign nonzero IDs in grid cells where ice has reached the minimum glacier thickness.
       ! Remove IDs in grid cells where ice is now thinnier than the minimum thickness.
       ! Adjust IDs to prevent spurious advance due to SMB differences in adjacent glaciers.

       call glacier_advance_retreat(&
            ewn,             nsn,           &
            parallel,                       &
            itest,   jtest,  rtest,         &
            nglacier,                       &
            glacier%minthck,                &  ! m
            model%geometry%thck,            &  ! m
            glacier%snow_annmean,           &  ! mm/yr w.e.
            glacier%Tpos_annmean,           &  ! deg C
            glacier%mu_star,                &  ! mm/yr/deg
            glacier%alpha_snow,             &  ! unitless
            glacier%cism_glacier_id_init,   &
            glacier%cism_glacier_id)

       ! Compute smb_glacier_id, which determines where the SMB is computed. It is the union of
       !       (1) cism_glacier_id > 0
       !       (2) cism_glacier_id_init > 0
       !       (3) cells adjacent to cells with cism_glacier_id > 0
       ! Thus, a glacier ID is associated with any cell that is currently or potentially glaciated.
       ! Cells are potentially glaciated if adjacent to current glacier cells.

       call update_smb_glacier_id(&
            ewn,           nsn,             &
            itest, jtest,  rtest,           &
            glacier%nglacier,               &
            glacier%snow_annmean,           &  ! mm/yr w.e.
            glacier%Tpos_annmean,           &  ! deg C
            glacier%mu_star,                &  ! mm/yr/deg
            glacier%alpha_snow,             &  ! unitless
            glacier%cism_glacier_id_init,   &  ! initial extent
            glacier%cism_glacier_id,        &  ! current extent
            glacier%smb_glacier_id,         &
            parallel)

       ! Compute smb_glacier_id_init, as needed for inversion
       ! Note: cism_glacier_id_init is passed in twice to match the interface;
       !       the second version is redundant.

       call update_smb_glacier_id(&
            ewn,           nsn,             &
            itest, jtest,  rtest,           &
            glacier%nglacier,               &
            glacier%snow_annmean,           &  ! mm/yr w.e.
            glacier%Tpos_annmean,           &  ! deg C
            glacier%mu_star,                &  ! mm/yr/deg
            glacier%alpha_snow,             &  ! unitless
            glacier%cism_glacier_id_init,   &  ! initial extent
            glacier%cism_glacier_id_init,   &  ! treated as current extent
            glacier%smb_glacier_id_init,    &
            parallel)

       ! Where smb_glacier_id_init > 0, make sure smb_glacier_id has the same value.
       ! This piece of code requires that smb_glacier_id_init is always computed,
       !  even if not inverting.

       where (glacier%smb_glacier_id_init > 0)
          glacier%smb_glacier_id = glacier%smb_glacier_id_init
       endwhere

       ! Using the new smb_glacier_id mask, compute model%climate%smb for the next year.
       !TODO - Reduce loop size?

       do j = 1, nsn
          do i = 1, ewn
             ng = glacier%smb_glacier_id(i,j)
             if (ng > 0) then
                model%climate%smb(i,j) = &
                     glacier%alpha_snow(ng)*glacier%snow_annmean(i,j) &
                   - glacier%mu_star(ng)*glacier%Tpos_annmean(i,j)
             else
                model%climate%smb(i,j) = 0.0d0
             endif
          enddo
       enddo

       ! In advanced or potential advanced cells, zero out any positive SMB.
       ! This inhibits further advance.

       where (glacier%cism_glacier_id_init == 0 .and. glacier%smb_glacier_id > 0)
          model%climate%smb = min(model%climate%smb, 0.0d0)
       endwhere

       call parallel_halo(model%climate%smb, parallel)

       ! If inverting, then repeat for the RGI and recent SMB

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
           glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then

          do j = 1, nsn
             do i = 1, ewn
                ng = glacier%smb_glacier_id(i,j)
                if (ng > 0) then
                   glacier%smb_rgi(i,j) = &
                        glacier%alpha_snow(ng)*glacier%snow_rgi_annmean(i,j) &
                        - glacier%mu_star(ng)*glacier%Tpos_rgi_annmean(i,j)
                   glacier%smb_recent(i,j) = &
                        glacier%alpha_snow(ng)*glacier%snow_recent_annmean(i,j) &
                        - glacier%mu_star(ng)*glacier%Tpos_recent_annmean(i,j)
                else
                   glacier%smb_rgi(i,j) = 0.0d0
                   glacier%smb_recent(i,j) = 0.0d0
                endif
             enddo
          enddo

          ! In advanced or potential advanced cells, zero out any positive SMB
          where (glacier%cism_glacier_id_init == 0 .and. glacier%smb_glacier_id > 0)
             glacier%smb_rgi = min(glacier%smb_rgi, 0.0d0)
             glacier%smb_recent = min(glacier%smb_recent, 0.0d0)
          endwhere

          call parallel_halo(glacier%smb_rgi, parallel)
          call parallel_halo(glacier%smb_recent, parallel)

       endif   ! set_mu_star

       if (verbose_glacier) then
          call point_diag(model%geometry%thck, 'After advance_retreat, thck', itest, jtest, rtest, 7, 7)
          call point_diag(glacier%cism_glacier_id_init, 'cism_glacier_id_init', itest, jtest, rtest, 7, 7)
          call point_diag(glacier%smb_glacier_id_init, 'smb_glacier_id_init', itest, jtest, rtest, 7, 7)
          call point_diag(glacier%cism_glacier_id, 'New cism_glacier_id', itest, jtest, rtest, 7, 7)
          call point_diag(glacier%smb_glacier_id, 'New smb_glacier_id', itest, jtest, rtest, 7, 7)
          call point_diag(glacier%smb_applied_annmean, 'smb_applied_annmean, previous yr', itest, jtest, rtest, 7, 7)
          call point_diag(smb_weight_init, 'smb_weight_init, previous yr', itest, jtest, rtest, 7, 7)
          call point_diag(glacier%Tpos_annmean, 'Tpos_annmean', itest, jtest, rtest, 7, 7)
          call point_diag(glacier%snow_annmean, 'snow_annmean', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%smb, 'climate%smb', itest, jtest, rtest, 7, 7)
          if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
              glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
             call point_diag(glacier%smb_rgi, 'smb_rgi', itest, jtest, rtest, 7, 7)
             call point_diag(glacier%smb_recent, 'smb_recent', itest, jtest, rtest, 7, 7)
          endif   ! set_mu_star
       endif

       ! Find the minimum and maximum SMB for each glacier in the baseline climate.
       ! Note: Include only cells that are part of the initial glacier extent.

       call glacier_smb_min_max(&
            ewn,           nsn,               &
            nglacier,                         &
            glacier%cism_glacier_id_init,     &
            model%climate%smb,                &
            smb_min,       smb_max)

       ! Compute AAR for each glacier in the baseline climate.

       ! (1) Include only cells that are part of the initial glacier extent
       call glacier_accumulation_area_ratio(&
            ewn,           nsn,               &
            parallel,                         &
            nglacier,                         &
            glacier%cism_glacier_id_init,     &
            cell_area_uniform,                &
            model%climate%smb,                &
            aar_init)

       ! (2) Include all cells in the glacier
       call glacier_accumulation_area_ratio(&
            ewn,           nsn,               &
            parallel,                         &
            nglacier,                         &
            glacier%cism_glacier_id,          &
            cell_area_uniform,                &
            model%climate%smb,                &
            aar)

       if (verbose_glacier .and. this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'Glacier SMB and AAR:'
          write(iulog,*) '    ng     smb_min   smb_max   AAR_initA     AAR'
          do ng = 1, nglacier
             if (glacier%volume_init(ng) > diagnostic_volume_threshold .or. ng == ngdiag) then  ! big glacier
                write(iulog,'(i10, 2f10.1, 2f10.4 )') ng, smb_min(ng), smb_max(ng), aar_init(ng), aar(ng)
             endif
          enddo
       endif

       ! If inverting for mu_star and alpha_snow, then repeat for the recent climate

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
           glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then

          call glacier_smb_min_max(&
               ewn,            nsn,           &
               nglacier,                      &
               glacier%cism_glacier_id_init,  &
               glacier%smb_recent,            &
               smb_min_recent, smb_max_recent)

          ! (1) Include only cells that are part of the initial glacier extent
          call glacier_accumulation_area_ratio(&
               ewn,           nsn,               &
               parallel,                         &
               nglacier,                         &
               glacier%cism_glacier_id_init,     &
               cell_area_uniform,                &
               glacier%smb_recent,               &
               aar_init_recent)

          ! (2) Include all cells in the glacier
          call glacier_accumulation_area_ratio(&
               ewn,           nsn,               &
               parallel,                         &
               nglacier,                         &
               glacier%cism_glacier_id,          &
               cell_area_uniform,                &
               glacier%smb_recent,               &
               aar_recent)

          if (verbose_glacier .and. this_rank == rtest) then
             write(iulog,*) ' '
             write(iulog,*) 'Recent SMB and AAR:'
             write(iulog,*) '    ng     smb_min   smb_max   AAR_initA     AAR'
             do ng = 1, nglacier
                if (glacier%volume_init(ng) > diagnostic_volume_threshold .or. ng == ngdiag) then  ! big glacier
                   write(iulog,'(i10, 2f10.1, 2f10.4 )') ng, smb_min_recent(ng), smb_max_recent(ng), &
                        aar_init_recent(ng), aar_recent(ng)
                endif
             enddo
          endif

       endif   ! set_mu_star

       ! Compute the area and volume of each glacier

       call glacier_area_volume(&
            ewn,           nsn,               &
            parallel,                         &
            nglacier,                         &
            glacier%cism_glacier_id,          &
            model%geometry%cell_area,         &  ! m^2
            model%geometry%thck,              &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area,                     &  ! m^2
            glacier%volume)                      ! m^3

       ! Repeat, summing over the initial glacier extent (no advanced cells)
       ! Note: area_init_extent < area_init if there has been any retreat

       call glacier_area_volume(&
            ewn,           nsn,               &
            parallel,                         &
            nglacier,                         &
            glacier%cism_glacier_id_init,     &
            model%geometry%cell_area,         &  ! m^2
            model%geometry%thck,              &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area_init_extent,         &  ! m^2
            glacier%volume_init_extent)          ! m^3

       if (verbose_glacier .and. this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'Update area (km^2) and volume (km^3) for glacier:', ngdiag
          write(iulog,*) ' Initial area and volume:', &
               glacier%area_init(ngdiag)/1.0d6, glacier%volume_init(ngdiag)/1.0d9
          write(iulog,*) ' Current area and volume:', &
               glacier%area(ngdiag)/1.0d6, glacier%volume(ngdiag)/1.0d9
          write(iulog,*) 'A and V over init extent:', &
               glacier%area_init_extent(ngdiag)/1.0d6, glacier%volume_init_extent(ngdiag)/1.0d9
          write(iulog,*) 'A and V over init extent:', &
               glacier%area_init_extent(ngdiag)/1.0d6, glacier%volume_init_extent(ngdiag)/1.0d9
       endif

       ! If inverting for thickness, compute the target area and volume

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then

          call glacier_area_volume(&
               ewn,           nsn,               &
               parallel,                         &
               nglacier,                         &
               glacier%cism_glacier_id_init,     &
               model%geometry%cell_area,         &  ! m^2
               glacier%thck_target,              &  ! m
               glacier%diagnostic_minthck,       &  ! m
               glacier%area_target,              &  ! m^2
               glacier%volume_target)               ! m^3

          if (verbose_glacier .and. this_rank == rtest) then
             write(iulog,*) ' Target area and volume:', &
                  glacier%area_target(ngdiag)/1.0d6, glacier%volume_target(ngdiag)/1.0d9
          endif

       endif

       if (verbose_glacier) then

          glacier_mask = 0
          where (glacier%cism_glacier_id_init == ngdiag) glacier_mask = 1
          count_cgii = parallel_global_sum(glacier_mask, parallel)

          glacier_mask = 0
          where (glacier%cism_glacier_id == ngdiag) glacier_mask = 1
          count_cgi = parallel_global_sum(glacier_mask, parallel)

          glacier_mask = 0
          where (glacier%smb_glacier_id_init == ngdiag) glacier_mask = 1
          count_sgii = parallel_global_sum(glacier_mask, parallel)

          glacier_mask = 0
          where (glacier%smb_glacier_id == ngdiag) glacier_mask = 1
          count_sgi = parallel_global_sum(glacier_mask, parallel)

          if (this_rank == rtest) then
             write(iulog,*) ' '
             write(iulog,*) 'Mask count, ng =', ngdiag
             write(iulog,*) 'count_cgii, count_cgi =', count_cgii, count_cgi
             write(iulog,*) 'count_sgii, count_sgi =', count_sgii, count_sgi
          endif

       endif   ! verbose

    endif   ! glacier_update_inverval

  end subroutine glissade_glacier_update

!****************************************************

  subroutine glacier_invert_mu_star(&
       ewn,              nsn,           &
       parallel,                        &
       itest,   jtest,   rtest,         &
       nglacier,         ngdiag,        &
       smb_glacier_id_init,             &
       smb_weight,                      &
       glacier_area_init,glacier_volume_init, &  ! diagnostic only
       snow,             Tpos,          &
       mu_star_const,                   &
       mu_star_min,      mu_star_max,   &
       beta_artm_max,                   &
       beta_artm_increment,             &
       alpha_snow,                      &
       mu_star,          beta_artm)

    ! Given an observational SMB target, invert for the parameter mu_star in the glacier SMB formula.
    ! This assumes that the input snow field does not need to be corrected.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest,         & ! coordinates of diagnostic cell
         nglacier,                    & ! total number of glaciers in the domain
         ngdiag                         ! CISM ID of diagnostic glacier

    type(parallel_type), intent(in) :: &
         parallel                       ! info for parallel communication

    integer, dimension(ewn,nsn), intent(in) :: &
         smb_glacier_id_init            ! smb_glacier_id based on the initial glacier extent

    real(dp), dimension(nglacier), intent(in) :: &
         smb_weight                     ! weight for applying SMB; < 1 if actual melt < potential melt

    real(dp), dimension(nglacier), intent(in) :: &
         glacier_area_init,           & ! initial glacier area (m^2); diagnostic only
         glacier_volume_init            ! initial glacier volume (m^2); diagnostic only

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow,                        & ! time-avg snowfall for each cell (mm/yr w.e.)
         Tpos                           ! time-avg of max(artm - tmlt, 0) for each cell (deg)

    real(dp), intent(in) :: &
         mu_star_const,               & ! default constant value of mu_star
         mu_star_min, mu_star_max,    & ! min and max allowed values of mu_star
         beta_artm_max,               & ! max allowed magnitude of beta_artm
         beta_artm_increment            ! increment of beta_artm in each iteration

    real(dp), dimension(nglacier), intent(inout) :: &
         alpha_snow                     ! prescribed glacier-specific snow factor (unitless)

    real(dp), dimension(nglacier), intent(inout) :: &
         mu_star,                     & ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)
         beta_artm                      ! correction to artm (deg C)

    ! local variables

    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         glacier_snow, glacier_Tpos,  & ! glacier-average snowfall and Tpos
         smb_baseline                   ! SMB in baseline climate

    character(len=100) :: message

    real(dp), parameter :: Tpos_min = 0.1d0     ! deg C available for melting, min value
                                                ! values too close to zero can result in high mu_star

    integer :: count_violate_1                  ! number of glaciers violating Eq. 1
    real(dp) :: area_violate_1                  ! total area of these glaciers (m^2)
    real(dp) :: volume_violate_1                ! total volume of these glaciers (m^3)
    real(dp) :: mu_eq1

    ! Compute mu_star for each glacier such that SMB = 0 over the initial extent.
    ! The initial extent can include an ablation zone of glacier-free cells adjacent
    !  to glacier-covered cells, with weights in the range [0,1].
    !
    ! The SMB for glacier ng is given by
    !      sum_ij(smb) = alpha_snow(ng)*sum_ij(snow) - mu_star(ng) * sum_ij(Tpos),
    ! where Tpos = max(artm - tmlt, 0),
    ! and sum_ij notes a sum over all cells (i,j) in the glacier.
    !
    ! Setting sum_ij(smb) = 0 and rearranging, we get
    ! (1)   mu_star(ng) = alpha_snow(ng)*sum_ij(snow) / sum_ij(Tpos)
    !
    ! Thus, given the annual average of snow and Tpos for each grid cell in a glacier,
    !  we can find mu_star such that SMB = 0.
    ! If mu_star lies outside a prescribed range, we adjust a parameter beta_artm,
    !  which in turn changes Tpos in a way that will bring mu_star in range.
    !
    ! Notes:
    ! (1) This approach works only for land-based glaciers.
    !     TODO: Modify for marine-terminating glaciers.
    ! (2) Assuming climatological forcing with smb_obs = 0, mu_star has nearly the same value
    !     throughout the inversion.  It changes slightly as surface elevation changes.

    if (verbose_glacier .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'In glacier_invert_mu_star'
    endif

    ! Compute weighted averages of Tpos and snow over each glacier

    call glacier_2d_to_1d_weighted(&
         ewn,           nsn,                   &
         parallel,                             &
         nglacier,                             &
         smb_glacier_id_init,                  &
         smb_weight,                           &
         snow,          glacier_snow)

    call glacier_2d_to_1d_weighted(&
         ewn,           nsn,                   &
         parallel,                             &
         nglacier,                             &
         smb_glacier_id_init,                  &
         smb_weight,                           &
         Tpos,          glacier_Tpos)

    if (verbose_glacier .and. this_rank == rtest) then
       ng = ngdiag
       write(iulog,*) ' '
       write(iulog,*) 'ng, snow and Tpos with weighting =', ng, glacier_snow(ng), glacier_Tpos(ng)
    endif

    ! For each glacier, compute the new mu_star. Adjust beta_artm if necessary.

    do ng = 1, nglacier

       if (glacier_snow(ng) == 0.0d0) then

          if (verbose_glacier .and. this_rank == rtest) then
             write(iulog,*) 'WARNING: snow = 0 for glacier', ng
             !TODO - Throw a fatal error?
          endif

          mu_star(ng) = mu_star_const

       else   ! glacier_snow > 0

          if (glacier_Tpos(ng) < Tpos_min) then

             ! There is little or no ablation anywhere on the glacier.
             ! Compensate by raising artm until there is some ablation.
             ! Prescribe mu for now.

             beta_artm(ng) = beta_artm(ng) + beta_artm_increment
             mu_star(ng) = mu_star_const

          else    ! Tpos >= Tpos_min

             ! Compute the value of mu_star that will give the desired SMB = 0 over the target area
             mu_star(ng) = (alpha_snow(ng)*glacier_snow(ng)) / glacier_Tpos(ng)

             ! Note: Would use the following commented-out equation if smb_obs /= 0
             ! mu_star(ng) = (alpha_snow(ng)*glacier_snow(ng) - glacier_smb_obs(ng)) / glacier_Tpos(ng)

             ! If mu_star is out of range (based on Eq. 1), then modify beta
             if (mu_star(ng) < mu_star_min) then
                ! This could happen if Tpos is too large. Compensate by cooling.
                beta_artm(ng) = beta_artm(ng) - beta_artm_increment
                mu_star(ng) = mu_star_min
             elseif (mu_star(ng) > mu_star_max) then
                ! This could happen if Tpos is too small. Compensate by warming.
                beta_artm(ng) = beta_artm(ng) + beta_artm_increment
                mu_star(ng) = mu_star_max
             endif

          endif   ! glacier_Tpos

       endif   ! glacier_snow

    enddo   ! ng

    ! Diagnostic checks

    ! Make sure the glacier variables are now in range

    do ng = 1, nglacier

       if (mu_star(ng) < mu_star_min .or. mu_star(ng) > mu_star_max) then
          if (this_rank == rtest) then
             write(iulog,*) 'WARNING, mu out of range: ng, mu =', ng, mu_star(ng)
          endif
       endif

       beta_artm(ng) = min(beta_artm(ng),  beta_artm_max)
       beta_artm(ng) = max(beta_artm(ng), -beta_artm_max)

    enddo   ! ng

    ! Check the mass balance. The goal is that all glaciers satisfy (1).

    count_violate_1 = 0
    area_violate_1 = 0.0d0
    volume_violate_1 = 0.0d0

    do ng = 1, nglacier

       smb_baseline(ng) = alpha_snow(ng)*glacier_snow(ng) - mu_star(ng)*glacier_Tpos(ng)
       if (glacier_Tpos(ng) > 0.0d0) then
          mu_eq1 = alpha_snow(ng) * glacier_snow(ng) / glacier_Tpos(ng)
       else
          mu_eq1 = 0.0d0
       endif

       ! Check whether the glacier violates Eq. (1)
       if (verbose_glacier .and. this_rank == rtest) then
          if (abs(smb_baseline(ng)) > eps08) then
!!             write(iulog,'(a60,i6,6f10.2)') 'Eq 1 violation, ng, snow, Tpos, init mu, adj mu, beta, smb :', &
!!                  ng, glacier_snow(ng), glacier_Tpos(ng), mu_eq1, mu_star(ng), beta_artm(ng), smb_baseline(ng)
             count_violate_1 = count_violate_1 + 1
             area_violate_1 = area_violate_1 + glacier_area_init(ng)
             volume_violate_1 = volume_violate_1 + glacier_volume_init(ng)
          endif
       endif

    enddo  ! ng

    if (verbose_glacier .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'Violations of Eq. 1 (SMB = 0, baseline climate):', count_violate_1
       write(iulog,*) '   Total area, volume =', area_violate_1/1.0d6, volume_violate_1/1.0d9
       write(iulog,*) ' '
       ng = ngdiag
       write(iulog,*) 'Balance solution, ng =', ng
       write(iulog,'(a30,3f12.4)') '   mu_star, alpha_snow, beta: ', &
            mu_star(ng), alpha_snow(ng), beta_artm(ng)
       write(iulog,'(a30,3f12.4)') '   Baseline snow, Tpos, SMB : ', &
            glacier_snow(ng), glacier_Tpos(ng), smb_baseline(ng)
    endif

  end subroutine glacier_invert_mu_star

!****************************************************

  subroutine glacier_invert_mu_star_alpha_snow(&
       ewn,              nsn,            &
       parallel,                         &
       itest,   jtest,   rtest,          &
       nglacier,         ngdiag,         &
       smb_glacier_id_init,              &
       smb_weight,                       &
       glacier_smb_obs,                  &
       glacier_area_init,glacier_volume_init, &  ! diagnostic only
       snow,             Tpos,           &
       snow_recent,      Tpos_recent,    &
       mu_star_const,                    &
       mu_star_min,      mu_star_max,    &
       alpha_snow_const,                 &
       alpha_snow_min,   alpha_snow_max, &
       beta_artm_max,                    &
       beta_artm_increment,              &
       mu_star,          alpha_snow,     &
       beta_artm)

    ! Given an observational SMB target, invert for the parameters mu_star and alpha_snow.
    ! Two conditions must be satisfied:
    ! SMB = 0 given input snow and Tpos, for a period with glaciers in balance.
    ! SMB = smb_obs given input snow_recent and Tpos_recent, for a period with glaciers out of balance.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest,         & ! coordinates of diagnostic cell
         nglacier,                    & ! total number of glaciers in the domain
         ngdiag                         ! CISM ID of diagnostic glacier

    type(parallel_type), intent(in) :: &
         parallel                       ! info for parallel communication

    integer, dimension(ewn,nsn), intent(in) :: &
         smb_glacier_id_init            ! smb_glacier_id based on the initial glacier extent

    real(dp), dimension(nglacier), intent(in) :: &
         smb_weight,                  & ! weight for applying SMB; < 1 if actual melt < potential melt
         glacier_smb_obs                ! observed glacier-average SMB (mm/yr w.e.)

    real(dp), dimension(nglacier), intent(in) :: &
         glacier_area_init,           & ! initial glacier area (m^2); diagnostic only
         glacier_volume_init            ! initial glacier volume (m^2); diagnostic only

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow,                        & ! time-avg snowfall for each cell (mm/yr w.e.)
         Tpos,                        & ! time-avg of max(artm - tmlt, 0) for each cell (deg)
         snow_recent,                 & ! time-avg snowfall for each cell (mm/yr w.e.), recent date
         Tpos_recent                    ! time-avg of max(artm - tmlt, 0) for each cell (deg), recent date

    real(dp), intent(in) :: &
         mu_star_const,                  & ! default constant value of mu_star
         mu_star_min, mu_star_max,       & ! min and max allowed values of mu_star
         alpha_snow_const,               & ! default constant value of alpha_snow
         alpha_snow_min, alpha_snow_max, & ! min and max allowed values of mu_star
         beta_artm_max,                  & ! max allowed magnitude of beta_artm
         beta_artm_increment               ! increment of beta_artm in each iteration

    real(dp), dimension(nglacier), intent(inout) :: &
         mu_star,                     & ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)
         alpha_snow,                  & ! glacier-specific snow factor (unitless)
         beta_artm                      ! correction to artm (deg C)

    ! local variables

    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         glacier_snow, glacier_Tpos,                & ! glacier-average snowfall and Tpos
         glacier_snow_recent, glacier_Tpos_recent,  & ! glacier-average snowfall_recent and Tpos_recent
         smb_baseline, smb_recent,                  & ! SMB in baseline and recent climates
         smb_recent_diff,                           & ! difference between modeled and observed SMB, recent climate
         denom

    character(len=100) :: message

    real(dp), parameter :: Tpos_min = 0.1d0        ! deg C available for melting, min value
                                                   ! values too close to zero can result in high mu_star

    integer :: count_violate_1, count_violate_2    ! number of glaciers violating Eq. 1 and Eq. 2
    real(dp) :: area_violate_1, area_violate_2     ! total area of these glaciers (m^2)
    real(dp) :: volume_violate_1, volume_violate_2 ! total volume of these glaciers (m^3)
    real(dp) :: mu_eq1, deltaT

    ! Compute mu_star and alpha_snow for each glacier such that
    ! (1) snow and Tpos combine to give SMB = 0
    ! (2) snow_recent and Tpos_recent combine to give SMB = smb_obs
    ! In both cases, the SMB is computed over the initial glacier extent.
    ! Here, the initial extent includes an ablation zone of glacier-free cells adjacent
    !  to glacier-covered cells.

    ! The SMB for glacier ng is given by
    !      sum_ij(smb) = alpha_snow * sum_ij(snow) - mu_star(ng) * sum_ij(Tpos),
    ! where Tpos = max(artm - tmlt, 0),
    ! and sum_ij notes a sum over all cells (i,j) in the glacier.
    !
    ! For glaciers in balance, this becomes (dropping the sum_ij notation)
    ! (1)            0  = alpha_snow * snow - mu_star * Tpos.
    !
    ! For glaciers observed to be out of balance, this becomes
    ! (2)       smb_obs = alpha_snow * snow_recent - mu_star * Tpos_recent.
    !
    ! Rearranging and solving, we get
    !              mu_star = (-smb_obs * snow) / D,
    !           alpha_snow = (-smb_obs * Tpos) / D,
    !              where D = snow*Tpos_recent - snow_recent*Tpos
    !
    ! Ideally, both mu_star and alpha_snow fall within physically realistic ranges.
    ! If not, there is some additional logic to adjust beta_artm such that the computed mu_star
    !  moves toward a realistic range.
    !
    ! Notes:
    !     This approach works only for land-based glaciers.
    !     TODO: Modify for marine-terminating glaciers.

    if (verbose_glacier .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'In glacier_invert_mu_star_alpha_snow'
    endif

    ! Compute weighted averages of Tpos and snow over each glacier

    call glacier_2d_to_1d_weighted(&
         ewn,                  nsn,       &
         parallel,                        &
         nglacier,                        &
         smb_glacier_id_init,             &
         smb_weight,                      &
         snow,                 glacier_snow)

    call glacier_2d_to_1d_weighted(&
         ewn,                  nsn,       &
         parallel,                        &
         nglacier,                        &
         smb_glacier_id_init,             &
         smb_weight,                      &
         Tpos,                 glacier_Tpos)

    call glacier_2d_to_1d_weighted(&
         ewn,                  nsn,       &
         parallel,                        &
         nglacier,                        &
         smb_glacier_id_init,             &
         smb_weight,                      &
         snow_recent,          glacier_snow_recent)

    call glacier_2d_to_1d_weighted(&
         ewn,                  nsn,       &
         parallel,                        &
         nglacier,                        &
         smb_glacier_id_init,             &
         smb_weight,                      &
         Tpos_recent,          glacier_Tpos_recent)

    if (verbose_glacier .and. this_rank == rtest) then
       ng = ngdiag
       write(iulog,*) ' '
       write(iulog,*) 'ng, snow and Tpos with weighting =', ng, glacier_snow(ng), glacier_Tpos(ng)
       write(iulog,*) 'recent snow and Tpos with weighting =', glacier_snow_recent(ng), glacier_Tpos_recent(ng)
    endif

    ! For each glacier, compute the new mu_star and alpha_snow

    do ng = 1, nglacier

       if (glacier_snow(ng) == 0.0d0) then

          if (verbose_glacier .and. this_rank == rtest) then
             write(iulog,*) 'WARNING: snow = 0 for glacier', ng
             !TODO - Throw a fatal error?
          endif

          mu_star(ng) = mu_star_const
          alpha_snow(ng) = alpha_snow_const

       else   ! glacier_snow > 0

          ! compute D = snow*Tpos_recent - snow_recent*Tpos
          denom(ng) = glacier_snow(ng)*glacier_Tpos_recent(ng) - glacier_snow_recent(ng)*glacier_Tpos(ng)

          if (glacier_Tpos(ng) < Tpos_min) then

             ! There is little or no ablation anywhere on the glacier in the baseline climate.
             ! Compensate by raising artm (along with artm_recent) until there is some ablation.
             ! Prescribe mu and alpha for now.

             beta_artm(ng) = beta_artm(ng) + beta_artm_increment
             alpha_snow(ng) = alpha_snow_const
             mu_star(ng) = mu_star_const

          else   ! Tpos >= Tpos_min; this implies denom > 0

             if (denom(ng) * glacier_smb_obs(ng) > 0.0d0) then

                ! The glacier is either gaining mass in a warming climate or losing mass in a cooling climate.
                ! This is unrealistic and may be due to mass-balance measurement error.
                ! To keep things simple, prescribe alpha and use Eq. (1) to compute mu.

                alpha_snow(ng) = alpha_snow_const
                mu_star(ng) = alpha_snow(ng) * glacier_snow(ng) / glacier_Tpos(ng)

             else   ! usual case; compute mu and alpha using the 2-equation scheme

                mu_star(ng)    = -glacier_smb_obs(ng)*glacier_snow(ng) / denom(ng)
                alpha_snow(ng) = -glacier_smb_obs(ng)*glacier_Tpos(ng) / denom(ng)

                ! Check for mu and alpha in range.
                ! If out of range, then we can try some adjustments.
                ! One adjustment (not yet tried) is to adjust smb_obs within its stated error.
                ! Another is to prescribe alpha and use Eq. (1) to compute mu.
                ! If mu is still out of range, then try adjusting beta to change the temperature.

                if (   mu_star(ng) <    mu_star_min .or.    mu_star(ng) > mu_star_max .or. &
                    alpha_snow(ng) < alpha_snow_min .or. alpha_snow(ng) > alpha_snow_max) then

                   ! Note the discrepancy
!                   if (verbose_glacier .and. this_rank == rtest) then
!                      write(iulog,'(a46,i6,6f10.3)') 'Out of range, ng, Tp, Tp_recent, D, B, alpha, mu:', &
!                           ng, glacier_Tpos(ng), glacier_Tpos_recent(ng), denom(ng), &
!                           glacier_smb_obs(ng), alpha_snow(ng), mu_star(ng)
!                   endif

                   ! There are a number of reasons this could happen.
                   ! Assuming that Tpos and therefore D are not too small, the most likely reason
                   !  is mass-balance measurement error.
                   ! To keep things simple, cap alpha and then use Eq. (1) to compute mu.

                   alpha_snow(ng) = min(alpha_snow(ng), alpha_snow_max)
                   alpha_snow(ng) = max(alpha_snow(ng), alpha_snow_min)

                   mu_star(ng) = alpha_snow(ng) * glacier_snow(ng) / glacier_Tpos(ng)

                endif   ! mu_star and alpha in range

             endif   ! denom * smb_obs > 0

             ! If mu_star is still out of range (based on Eq. 1), then modify beta.

             if (mu_star(ng) < mu_star_min) then
                ! This could happen if Tpos is too large. Compensate by cooling.
                beta_artm(ng) = beta_artm(ng) - beta_artm_increment
                mu_star(ng) = mu_star_min
             elseif (mu_star(ng) > mu_star_max) then
                ! This could happen if Tpos is too small. Compensate by warming.
                beta_artm(ng) = beta_artm(ng) + beta_artm_increment
                mu_star(ng) = mu_star_max
             endif

          endif   ! glacier_Tpos

       endif   ! glacier_snow

    enddo   ! ng

    ! Diagnostic checks

    ! Make sure the glacier variables are now in range

    do ng = 1, nglacier

       if (mu_star(ng) < mu_star_min .or. mu_star(ng) > mu_star_max) then
          if (this_rank == rtest) then
             write(iulog,*) 'WARNING, mu out of range: ng, mu =', ng, mu_star(ng)
          endif
       endif

       if (alpha_snow(ng) < alpha_snow_min .or. alpha_snow(ng) > alpha_snow_max) then
          if (this_rank == rtest) then
             write(iulog,*) 'WARNING, alpha out of range: ng, alpha =', ng, alpha_snow(ng)
          endif
       endif

!       if (abs(beta_artm(ng)) > beta_artm_max) then
!          if (this_rank == rtest) then
!             write(iulog,*) 'WARNING, beta out of range: ng, beta =', ng, beta_artm(ng)
!          endif
!       endif

       beta_artm(ng) = min(beta_artm(ng),  beta_artm_max)
       beta_artm(ng) = max(beta_artm(ng), -beta_artm_max)

    enddo   ! ng

    ! Check the mass balance for the baseline and recent climates.
    ! The goal is that all glaciers satisfy (1), and most satisfy (2).

    count_violate_1 = 0
    count_violate_2 = 0
    area_violate_1 = 0.0d0
    area_violate_2 = 0.0d0
    volume_violate_1 = 0.0d0
    volume_violate_2 = 0.0d0

    do ng = 1, nglacier

       smb_baseline(ng) = alpha_snow(ng)*glacier_snow(ng) - mu_star(ng)*glacier_Tpos(ng)
       smb_recent(ng)   = alpha_snow(ng)*glacier_snow_recent(ng) - mu_star(ng)*glacier_Tpos_recent(ng)
       smb_recent_diff(ng) = smb_recent(ng) - glacier_smb_obs(ng)

       if (glacier_Tpos(ng) > 0.0d0) then
          mu_eq1 = alpha_snow(ng) * glacier_snow(ng) / glacier_Tpos(ng)
       else
          mu_eq1 = 0.0d0
       endif

       ! Check whether the glacier violates Eq. (1) and/or Eq. (2)

       if (verbose_glacier .and. this_rank == rtest) then
          if (abs(smb_baseline(ng)) > eps08) then
!!             write(iulog,'(a60,i6,6f10.2)') 'Eq 1 violation, ng, snow, Tpos, init mu, adj mu, beta, smb :', &
!!                  ng, glacier_snow(ng), glacier_Tpos(ng), mu_eq1, mu_star(ng), beta_artm(ng), smb_baseline(ng)
             count_violate_1 = count_violate_1 + 1
             area_violate_1 = area_violate_1 + glacier_area_init(ng)
             volume_violate_1 = volume_violate_1 + glacier_volume_init(ng)
          endif
          if (abs(smb_recent_diff(ng)) > eps08) then
!!             write(iulog,*) '   Violation of Eq. 2: ng, smb_recent_diff =', ng, smb_recent_diff(ng)
             count_violate_2 = count_violate_2 + 1
             area_violate_2 = area_violate_2 + glacier_area_init(ng)
             volume_violate_2 = volume_violate_2 + glacier_volume_init(ng)
          endif
       endif

    enddo  ! ng

    if (verbose_glacier .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'Violations of Eq. 1 (SMB = 0, baseline climate):', count_violate_1
       write(iulog,*) '   Total area, volume =', area_violate_1/1.0d6, volume_violate_1/1.0d9
       write(iulog,*) 'Violations of Eq. 2 (SMB = SMB_obs, recent climate):', count_violate_2
       write(iulog,*) '   Total area, volume =', area_violate_2/1.0d6, volume_violate_2/1.0d9
       write(iulog,*) ' '
       ng = ngdiag
       write(iulog,*) 'Balance solution, ng =', ng
       write(iulog,'(a30,3f12.4)') '   mu_star, alpha_snow, beta: ', &
            mu_star(ng), alpha_snow(ng), beta_artm(ng)
       write(iulog,'(a30,3f12.4)') '   Baseline snow, Tpos, SMB : ', &
            glacier_snow(ng), glacier_Tpos(ng), smb_baseline(ng)
       write(iulog,'(a30,3f12.4)') '     Recent snow, Tpos, SMB : ', &
            glacier_snow_recent(ng), glacier_Tpos_recent(ng), smb_recent(ng)
    endif

  end subroutine glacier_invert_mu_star_alpha_snow

!****************************************************

  subroutine glacier_invert_powerlaw_c(&
       ewn,              nsn,              &
       itest,   jtest,   rtest,            &
       powerlaw_c_min,   powerlaw_c_max,   &
       babc_timescale,   babc_thck_scale,  &
       babc_relax_factor,                  &
       stag_thck,        stag_thck_target, &
       stag_dthck_dt,                      &
       powerlaw_c_relax,                   &
       powerlaw_c)

    ! Given the current ice thickness, rate of thickness change, and target thickness,
    ! invert for the parameter powerlaw_c in the relationship for basal sliding.
    ! Note: This subroutine is similar to subroutine invert_basal_friction in glissade_inversion.F90.
    !       The main difference is that it does not include a smoothing term.
    !       In cells that become ice-free, Cp will relax back toward its default value.
    ! TODO: Call subroutine invert_basal_friction instead?

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest            ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         powerlaw_c_min, powerlaw_c_max ! min and max allowed values of powerlaw_c (Pa (m/yr)^(-1/3))

    real(dp), intent(in) :: &
         babc_timescale,              & ! inversion timescale for powerlaw_c (yr)
         babc_thck_scale,             & ! inversion thickness scale for powerlaw_c (m)
         babc_relax_factor              ! controls strength of relaxation to default values (unitless)

    real(dp), dimension(ewn-1,nsn-1), intent(in) :: &
         stag_thck,                   & ! ice thickness at vertices (m)
         stag_thck_target,            & ! target ice thickness at vertices (m)
         stag_dthck_dt                  ! rate of change of ice thickness at vertices (m/yr)

    real(dp), intent(in) :: &
         powerlaw_c_relax               ! powerlaw_c value to which we relax; must be > 0

    real(dp), dimension(ewn-1,nsn-1), intent(inout) :: &
         powerlaw_c                     ! basal friction field to be adjusted (Pa (m/yr)^(-1/3))

    ! local variables

    integer :: i, j

    real(dp), dimension(ewn-1,nsn-1) :: &
         stag_dthck                     ! stag_thck - stag_thck_target (m)

    real(dp) :: &
         dpowerlaw_c,                 & ! change in powerlaw_c
         term_thck, term_dHdt,        & ! tendency terms for powerlaw_c based on thickness target
         term_relax                     ! tendency terms based on relaxation to default value

    !WHL - debug
    real(dp), dimension(ewn-1,nsn-1) ::  &
         logC,                        & ! log_10(powerlaw_c)
         dlogC                          ! change in log_10(powerlaw_c)

    real(dp) :: logC_relax              ! log_10(powerlaw_c_relax)

    real(dp), parameter :: logmin = -99.d0   ! arbitrary negative value;
                                             ! values of log(c) below logmin are considered non-physical

    if (verbose_glacier .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'In glacier_invert_powerlaw_c'
    endif

    if (babc_thck_scale > 0.0d0 .and. babc_timescale > 0.0d0) then

       stag_dthck(:,:) = stag_thck(:,:) - stag_thck_target(:,:)

       ! Compute the log (base 10) of the current Cp.
       ! We work with log(C) instead of C itself, because the physical effects of changing C
       ! by an amount dC are much greater at low C than at high C.
       where (powerlaw_c > 0.0d0)
          logC = log10(powerlaw_c)
       elsewhere
          logC = logmin
       endwhere

       ! initialize
       dlogC = 0.0d0
       logC_relax = log10(powerlaw_c_relax)

       ! Loop over vertices
       do j = 1, nsn-1
          do i = 1, ewn-1

             ! Compute and sum the three tendency terms
             term_thck = -stag_dthck(i,j) / (babc_thck_scale*babc_timescale)
             term_dHdt = -stag_dthck_dt(i,j) * 2.0d0 / babc_thck_scale

             if (logC(i,j) > logmin) then
                term_relax = -babc_relax_factor * (logC(i,j) - logC_relax) / babc_timescale
             else
                term_relax = 0.0d0
             endif

             dlogC(i,j) = (term_thck + term_dHdt + term_relax) * glacier_update_interval

             ! Limit to prevent a large change in one step
             ! Note: glacier_update_interval has units of yr.
             if (abs(dlogC(i,j)) > 0.1d0 * glacier_update_interval) then
                if (dlogC(i,j) > 0.0d0) then
                   dlogC(i,j) =  0.1d0 * glacier_update_interval
                else
                   dlogC(i,j) = -0.1d0 * glacier_update_interval
                endif
             endif

             ! Update log(C)
             logC(i,j) = logC(i,j) + dlogC(i,j)

             ! Convert log(C) back to C
             if (logC(i,j) > logmin) then
                powerlaw_c(i,j) = 10.d0**(logC(i,j))
             else
                powerlaw_c(i,j) = 0.0d0
             endif

             ! Limit to a physically reasonable range
             powerlaw_c(i,j) = min(powerlaw_c(i,j), powerlaw_c_max)
             powerlaw_c(i,j) = max(powerlaw_c(i,j), powerlaw_c_min)

             if (verbose_glacier .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'Invert for powerlaw_c: rank, i, j =', this_rank, i, j
                write(iulog,*) 'H, H_target (m)', stag_thck(i,j), stag_thck_target(i,j)
                write(iulog,*) 'dH_dt (m/yr):', stag_dthck_dt(i,j)
                write(iulog,*) 'dt (yr), term_thck*dt, term_dHdt*dt:', glacier_update_interval, &
                     term_thck*glacier_update_interval, term_dHdt*glacier_update_interval
                write(iulog,*) 'relax term:', term_relax*glacier_update_interval
                write(iulog,*) 'dlogC, new powerlaw_c:', dlogC(i,j), powerlaw_c(i,j)
             endif

          enddo  ! i
       enddo   ! j

    else   ! thck_scale or timescale = 0

       call write_log &
            ('Must have thck_scale and timescale > 0 for glacier powerlaw_c inversion', GM_FATAL)

    endif

    if (verbose_glacier) then
       call point_diag(stag_thck, 'stag_thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(stag_thck_target, 'stag_thck_target (m)', itest, jtest, rtest, 7, 7)
       call point_diag(stag_dthck, 'stag_thck - stag_thck_target (m)', itest, jtest, rtest, 7, 7)
       call point_diag(stag_dthck_dt, 'stag_dthck_dt (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(powerlaw_c, 'new powerlaw_c', itest, jtest, rtest, 7, 7)
    endif   ! verbose_glacier

  end subroutine glacier_invert_powerlaw_c

!****************************************************

  subroutine glacier_calc_snow(&
       ewn,       nsn,       &
       snow_threshold_min,   &
       snow_threshold_max,   &
       precip,               &
       artm,                 &
       snow)

    ! Given the precip rate and surface air temperature, compute the snowfall rate.
    ! Assume that the ratio snow/precip is given by a linear ramp between two thresholds.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn                  ! number of cells in each horizontal direction

    real(dp), intent(in) :: &
         snow_threshold_min,     & ! air temperature (deg C) below which all precip falls as snow
         snow_threshold_max        ! air temperature (deg C) above which all precip falls as rain

    real(dp), dimension(ewn,nsn), intent(in) :: &
         precip,                 & ! precipitation rate (mm/yr w.e.) at reference elevation usrf_ref
         artm                      ! surface air temperature (deg C)

    real(dp), dimension(ewn,nsn), intent(out) :: &
         snow                      ! snowfall rate (mm/yr w.e.)

    ! temperature correction; precip falls as snow only at cold temperatures
    where(artm > snow_threshold_max)
       snow = 0.0d0
    elsewhere (artm < snow_threshold_min)
       snow = precip
    elsewhere
       snow = precip * (snow_threshold_max - artm) / (snow_threshold_max - snow_threshold_min)
    endwhere

  end subroutine glacier_calc_snow

!****************************************************

  subroutine glacier_redistribute_advanced_ice(&
       ewn,             nsn,          &
       parallel,                      &
       itest,   jtest,  rtest,        &
       nglacier,        ngdiag,       &
       glacier_update_interval,       & ! yr
       cell_area,                     & ! m^2
       thinning_rate_advanced_ice,    & ! m/yr
       cism_glacier_id_init,          &
       smb_glacier_id,                &
       smb,                           & ! m/yr
       thck)                            ! m

    ! Limit glacier advance in the accumulation zone.
    ! This applies to grid cells that are initially ice-free, into which ice is advected.
    ! The fix here is to thin the ice in these cells at a prescribed rate and
    !  redistribute the mass conservatively across the glacier.

    use cism_parallel, only: parallel_halo, parallel_global_sum_patch

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest,         & ! coordinates of diagnostic cell
         nglacier,                    & ! number of glaciers
         ngdiag                         ! CISM ID of diagnostic glacier

    type(parallel_type), intent(in) :: &
         parallel                       ! info for parallel communication

    real(dp), intent(in) :: &
         glacier_update_interval,     & ! time interval (yr) of the glacier update, typically 1 yr
         thinning_rate_advanced_ice     ! thinning rate (m/yr) where glaciers advance in the accumulation zone

    real(dp), dimension(ewn,nsn), intent(in) :: &
         cell_area                      ! grid cell area (m^2)

    integer, dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init,        & ! integer glacier ID at the start of the run
         smb_glacier_id                 ! integer ID for current glacier cells and adjacent glacier-free cells

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         smb                            ! surface mass balance (m/yr)

    real(dp), dimension(ewn,nsn), intent(inout) ::  &
         thck                           ! ice thickness (m)

    ! local variables

    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         glacier_area_init,           & ! glacier area based on cism_glacier_id_init
         glacier_vol_removed,         & ! total volume (m^3) removed from each advanced cells in each glacier
         glacier_dthck,               & ! thickness (m) added over the initial extent of each glacier
         glacier_vol_1,               & ! volume (m^3) of each glacier before thinning and restribution
         glacier_vol_2                  ! volume (m^3) of each glacier after thinning and restribution

    real(dp), dimension(ewn,nsn) ::  &
         dthck                          ! thickness removed (m)

    integer, dimension(ewn,nsn) :: &
         glacier_id                     ! temporary glacier ID

    glacier_id = max(cism_glacier_id_init, smb_glacier_id)

    ! Compute the total volume of each glacier before limiting advance.
    ! Note: This includes adjacent glacier-free cells that might have a small nonzero thickness
    !       (i.e., cism_glacier_id = 0 but smb_glacier_id > 0).

    glacier_vol_1 = parallel_global_sum_patch(cell_area*thck, nglacier, glacier_id, parallel)

    ! compute the area of each glacier over its initial extent

    glacier_area_init = parallel_global_sum_patch(cell_area, nglacier, cism_glacier_id_init, parallel)

    ! Compute thinning in advanced grid cells
    ! This includes potential advanced cells adjacent to current glacier cells.
    ! Note: Currently, SMB is set to 0 in advanced cells where SMB would be > 0 otherwise.
    !       The logic below (smb >= 0) ensures that ice in these cells is thinned.

    dthck = 0.0d0
    where (cism_glacier_id_init == 0 .and. smb_glacier_id > 0)   ! advanced cell
       where (smb >= 0.0d0)   ! accumulation zone
          dthck = min(thinning_rate_advanced_ice*glacier_update_interval, thck)
          thck = thck - dthck
       endwhere
    endwhere
    glacier_vol_removed = parallel_global_sum_patch(cell_area*dthck, nglacier, smb_glacier_id, parallel)

    ! Assuming conservation of volume, compute the thickness to be added to each glacier.
    ! Only cells within the initial glacier extent can thicken.
    where (glacier_area_init > 0.0d0)
       glacier_dthck = glacier_vol_removed / glacier_area_init
    elsewhere
       glacier_dthck = 0.0d0
    endwhere

    ! Redistribute the ice volume over the initial extent of each glacier
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id_init(i,j)
          if (ng > 0) then
             thck(i,j) = thck(i,j) + glacier_dthck(ng)
          endif
       enddo
    enddo

    ! Halo update
    call parallel_halo(thck, parallel)

    ! Compute the volume of each glacier after limiting advance

    glacier_vol_2 = parallel_global_sum_patch(cell_area*thck, nglacier, glacier_id, parallel)

    ! conservation check
    do ng = 1, nglacier
       if (abs(glacier_vol_2(ng) - glacier_vol_1(ng)) > eps08*glacier_vol_1(ng)) then
          write(iulog,*) 'redistribute advanced ice, conservation error: ng, vol_1, vol_2:', &
               ng, glacier_vol_1(ng)/1.d9, glacier_vol_2(ng)/1.d9
          call write_log('Volume conservation error, redistribute advanced ice', GM_FATAL)
       endif
    enddo

  end subroutine glacier_redistribute_advanced_ice

  !****************************************************

  subroutine glacier_advance_retreat(&
       ewn,             nsn,         &
       parallel,                     &
       itest,   jtest,  rtest,       &
       nglacier,                     &
       glacier_minthck,              &
       thck,                         &
       snow,                         &
       Tpos,                         &
       mu_star,                      &
       alpha_snow,                   &
       cism_glacier_id_init,         &
       cism_glacier_id)

    ! Allow glaciers to advance and retreat.
    !
    ! The rules are as follows:
    ! - At start-up, glaciated cells have cism_glacier_id in the range (1, nglacier).
    !   Other cells have cism_glacier_id = 0.
    !   The initial cism_glacier_id array is saved as cism_glacier_id_init.
    ! - If a cell has H <= glacier_minthck and cism_glacier_id > 0, we set cism_glacier_id = 0.
    !   It is no longer considered to be glaciated.
    !   Here, glacier_minthck is a threshold for counting ice as part of a glacier.
    !   By default, glacier_minthck = model%numerics%thklim, typically 1 m.
    !   (Actually, glacier_minthck is slightly less than thklim, to make sure these cells
    !   are not dynamically active.)
    ! - When a cell has H > glacier_minthck and cism_glacier_id = 0, we give it a nonzero ID:
    !   either (1) cism_glacier_id_init, if the initial ID > 0,
    !   or (2) the ID of a glaciated neighbor (the one with the most negative SMB,
    !   if there is more than one).
    ! - In rare cases, there is no glaciated neighbor. This can happen when a few cells
    !   with H close to glacier_minthck are cut off from the parent glacier.
    !   With SMB = 0, they will slowly thin dynamically, but this can take a long time.
    !   It is simpler just to set H = 0.

    use cism_parallel, only: parallel_globalindex, parallel_halo

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest,         & ! coordinates of diagnostic cell
         nglacier                       ! number of glaciers

    type(parallel_type), intent(in) :: &
         parallel                       ! info for diagnostic only

    real(dp), intent(in) :: &
         glacier_minthck                ! min ice thickness (m) counted as part of a glacier

    real(dp), dimension(ewn,nsn), intent(inout) ::  &
         thck                           ! ice thickness (m)

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         snow,                        & ! annual mean snowfall (mm/yr w.e.)
         Tpos                           ! annual mean Tpos = max(T - Tmlt, 0)

    real(dp), dimension(nglacier), intent(in) :: &
         mu_star,                     & ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)
         alpha_snow                     ! glacier-specific snow factor (unitless)

    integer, dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init           ! cism_glacier_id at the start of the run

    integer, dimension(ewn,nsn), intent(inout) :: &
         cism_glacier_id                ! current cism glacier_id, > 0 for glaciated cells

    ! local variables

    integer, dimension(ewn,nsn) :: &
         cism_glacier_id_old            ! old value of cism_glacier_id

    real(dp) :: &
         smb_min,                     & ! min SMB among a cell and its neighbors
         smb_potential                  ! SMB if the cell were in a neighbor glacier

    integer :: i, j, ii, jj, ip, jp
    integer :: iglobal, jglobal
    integer :: ng, ng_init, ng_neighbor, ng_min
    logical :: found_neighbor

    real(dp), parameter :: big_number = 1.d+20   ! arbitrary large value

    if (verbose_glacier .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'In glacier_advance_retreat'
    endif

    ! Check for retreat: cells with cism_glacier_id > 0 but H < glacier_minthck

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0 .and. thck(i,j) <= glacier_minthck) then
             if (verbose_glacier .and. this_rank==rtest) then
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                write(iulog,*) 'Set ID = 0: ig, jg, old ID, thck =', &
                     iglobal, jglobal, ng, thck(i,j)
             endif
             cism_glacier_id(i,j) = 0
          endif
       enddo
    enddo

    ! Check for advance: cells with cism_glacier_id = 0 but H > H_min

    ! Save a copy of the current cism_glacier_id.
    ! This prevents the algorithm from depending on the loop direction.
    cism_glacier_id_old(:,:) = cism_glacier_id(:,:)

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng_init = cism_glacier_id_init(i,j)
          ng = cism_glacier_id_old(i,j)

          if (ng == 0 .and. thck(i,j) > glacier_minthck) then
             ! assign this cell its original ID, if > 0
             if (ng_init > 0) then
                cism_glacier_id(i,j) = ng_init
                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   write(iulog,*) 'Set ID = init ID: ig, jg, new ID, thck =',&
                        iglobal, jglobal, cism_glacier_id(i,j), thck(i,j)
                endif
             else  ! assign the ID of an adjacent glaciated cell, if possible
                found_neighbor = .false.
                smb_min = big_number
                ng_min = 0
                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   write(iulog,*) 'Look for glaciated neighbor: ig, jg =', iglobal, jglobal
                endif
                do jj = -1, 1
                   do ii = -1, 1
                      if (ii /= 0 .or. jj /= 0) then  ! edge or diagonal neighbor
                         ip = i + ii
                         jp = j + jj
                         ng_neighbor = cism_glacier_id_old(ip,jp)
                         if (ng_neighbor > 0) then
                            found_neighbor = .true.
                            ! compute the potential SMB, assuming cell (i,j) is in glacier ng_neighbor
                            smb_potential = alpha_snow(ng_neighbor)*snow(i,j) &
                                             - mu_star(ng_neighbor)*Tpos(i,j)
                            if (smb_potential < smb_min) then
                               smb_min = smb_potential
                               ng_min = ng_neighbor
                            endif
                         endif   ! neighbor cell is glaciated
                      endif   ! neighbor cell
                   enddo   ! ii
                enddo   ! jj

                if (found_neighbor) then
                   cism_glacier_id(i,j) = ng_min  ! glacier with the most negative SMB
                   if (verbose_glacier .and. this_rank == rtest) then
                      call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                      write(iulog,*) '  Set ID = neighbor ID, ig, jg, ID, H, smb =', &
                           iglobal, jglobal, cism_glacier_id(i,j), thck(i,j), smb_min
                   endif
                else  ! no adjacent glacier cell
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   write(iulog,*) '  Warning, did not find neighbor, ig, jg =', iglobal, jglobal
                   write(iulog,*) '  Setting H = 0'
                   thck(i,j) = 0.0d0  !TODO - anything else to zero out?
                endif   ! found_neighbor

             endif   ! cism_glacier_id_init > 0
          endif   ! ng = 0, H > minthck
       enddo   ! i
    enddo   ! j

    call parallel_halo(thck, parallel)
    call parallel_halo(cism_glacier_id, parallel)

    ! Check advanced cells (beyond the initial extent) for problematic glacier IDs.
    ! This code protects against glacier 'pirating', which can occur when an advanced cell
    !  is adjacent to two different glaciers, call them A and B.
    ! Suppose the cell is fed primarily by glacier A but has the same ID as glacier B,
    !  and has a more positive SMB as a result of belonging to B rather than A.
    ! Then glacier B is pirating ice from glacier A and can advance spuriously.
    ! Here, for each advanced cell (cism_glacier_id_init = 0, cism_glacier_id > 0), we check
    !  whether the cell's SMB would be more negative if it were in a different neighbor glacier.
    ! If so, the ID is switched.

    ! Save a copy of the current cism_glacier_id.
    cism_glacier_id_old = cism_glacier_id

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng_init = cism_glacier_id_init(i,j)
          ng = cism_glacier_id_old(i,j)

          if (ng_init == 0 .and. ng > 0) then ! advanced cell
             smb_min = alpha_snow(ng)*snow(i,j) - mu_star(ng)*Tpos(i,j) ! current SMB
             ng_min = 0
             ! Identify the neighbor with the most negative SMB
             do jj = -1, 1
                do ii = -1, 1
                   if (ii /= 0 .or. jj /= 0) then  ! edge or diagonal neighbor
                      ip = i + ii
                      jp = j + jj
                      ng_neighbor = cism_glacier_id_old(ip,jp)
                      if (ng_neighbor > 0) then
                         ! compute the potential SMB, assuming cell (i,j) is in glacier ng_neighbor
                         smb_potential = alpha_snow(ng_neighbor)*snow(i,j) - mu_star(ng_neighbor)*Tpos(i,j)
                         if (smb_potential < smb_min) then
                            smb_min = smb_potential
                            ng_min = ng
                         endif
                      endif  ! neighbor is glaciated
                   endif   ! neighbor cell
                enddo   ! ii
             enddo   ! jj

             if (ng_min > 0 .and. ng_min /= ng) then
                ! Move this cell to the adjacent glacier, resulting in a more negative SMB
                cism_glacier_id(i,j) = ng_min
                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   write(iulog,*) '   Transfer to adjacent glacier, old and new IDs =', &
                        cism_glacier_id_old(i,j), cism_glacier_id(i,j)
                endif
             endif

          endif   ! advanced cell
       enddo   ! i
    enddo   ! j

    call parallel_halo(cism_glacier_id, parallel)

  end subroutine glacier_advance_retreat

!****************************************************

  subroutine update_smb_glacier_id(&
         ewn,           nsn,      &
         itest, jtest,  rtest,    &
         nglacier,                &
         snow,                    &
         Tpos,                    &
         mu_star,                 &
         alpha_snow,              &
         cism_glacier_id_init,    &
         cism_glacier_id,         &
         smb_glacier_id,          &
         parallel)

    ! Based on the current glacier extent, compute a mask of cells that can have a nonzero SMB.
    !
    ! The rules for smb_glacier_id are as follows:
    ! (1) Where cism_glacier_id > 0, set smb_glacier_id = cism_glacier_id.
    ! (2) In retreated cells (cism_glacier_id = 0, cism_glacier_id_init > 0), set smb_glacier_id = cism_glacier_id_init.
    ! (3) In potential advanced grid cells (cism_glacier_id = 0 but adjacent to cells with cism_glacier_id > 0),
    !     set smb_glacier_id to the neighboring value of cism_glacier_id.
    !     If there is more than one neighbor glacier, choose the one that would result in the most negative SMB.
    ! (4) In other cells, no SMB is needed and smb_glacier_id = 0.
    !
    ! The logic for smb_glacier_id_init is the same, except that rule (2) is redundant
    ! since the initial and 'current' extents are the same.

    use cism_parallel, only: parallel_halo, parallel_globalindex

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier,                    & ! total number of glaciers in the domain
         itest, jtest, rtest            ! coordinates of diagnostic point

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         snow,                        & ! annual mean snowfall (mm/yr w.e.)
         Tpos                           ! annual mean Tpos = max(T - Tmlt, 0)

    real(dp), dimension(nglacier), intent(in) :: &
         mu_star,                     & ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)
         alpha_snow                     ! glacier-specific snow factor (unitless)

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id_init,        & ! integer glacier ID in the range (1, nglacier); initial value
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier); current value
                                        ! = 0 in cells without glaciers

    integer, dimension(ewn,nsn), intent(out) ::  &
         smb_glacier_id                 ! integer glacier ID in the range (1, nglacier), based on input extent

    type(parallel_type), intent(in) :: parallel

    ! local variables
    integer :: i, j, ii, jj, ng, ng_min
    integer :: ip, jp
    integer :: iglobal, jglobal

    real(dp), parameter :: big_number = 1.d+20   ! arbitrary large value

    real(dp) :: &
         smb_potential,               & ! potential SMB in a given cell outside the initial footprint
         smb_min                        ! min value of SMB for a given cell with glacier-covered neighbors

    ! Initialize to cism_glacier_id
    smb_glacier_id = cism_glacier_id

    ! Set smb_glacier_id = cism_glacier_id_init in retreated cells
    where (smb_glacier_id == 0 .and. cism_glacier_id_init > 0)
       smb_glacier_id = cism_glacier_id_init
    endwhere

    ! Where cism_glacier_id = 0, look for neighbors with cism_glacier_id > 0.
    ! Extend smb_glacier_id to these cells.

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          if (cism_glacier_id_init(i,j) == 0 .and. cism_glacier_id(i,j) == 0) then ! glacier-free cell
             ! find the adjacent glacier-covered cell (if any) with the most negative SMB
             smb_min = big_number
             ng_min = 0
             do jj = -1,1
                do ii = -1,1
                   if (ii /= 0 .or. jj /= 0) then  ! edge or diagonal neighbor
                      ip = i + ii
                      jp = j + jj
                      if (cism_glacier_id(ip,jp) > 0) then  ! adjacent glacier cell
                         ng = cism_glacier_id(ip,jp)
                         ! compute the potential SMB, assuming cell (i,j) is in glacier ng
                         smb_potential = alpha_snow(ng)*snow(i,j) - mu_star(ng)*Tpos(i,j)
                         if (smb_potential < smb_min) then
                            smb_min = smb_potential
                            ng_min = ng
                         endif
                      endif   ! cism_glacier_id > 0
                   endif   ! neighbor cell
                enddo   ! ii
             enddo   ! jj
             ! If there are any adjacent glacier cells with ng > 0, add cell (i,j) to the mask
             if (ng_min > 0) then
                smb_glacier_id(i,j) = ng_min
!                if (verbose_glacier .and. this_rank == rtest) then
!                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!                   write(iulog,*) 'Set smb_glacier_id = neighbor ID: ig, jg, smb_min, neighbor ID =', &
!                        iglobal, jglobal, smb_min, smb_glacier_id(i,j)
!                endif
             endif
          endif   ! cism_glacier_id_init = cism_glacier_id = 0
       enddo   ! i
    enddo   ! j

    call parallel_halo(smb_glacier_id, parallel)

  end subroutine update_smb_glacier_id

!****************************************************

  subroutine glacier_2d_to_1d(&
       ewn,           nsn,              &
       parallel,                        &
       nglacier,      cism_glacier_id,  &
       field_2d,      glacier_field)

    ! Given a 2D field, compute the average of the field over each glacier
    !TODO - Pass in cellarea to compute an area average.

    use cism_parallel, only: parallel_global_sum_patch

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    type(parallel_type), intent(in) :: &
         parallel                       ! info for parallel communication

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         field_2d                       ! 2D field to be averaged over glaciers

    real(dp), dimension(nglacier), intent(out) ::  &
         glacier_field                  ! field average over each glacier

    ! local variables

    integer :: i, j, ng

    integer, dimension(nglacier) :: ncells_glacier

    integer, dimension(ewn,nsn) :: ones   ! matrix = 1 everywhere

    ones(:,:) = 1

    ncells_glacier = parallel_global_sum_patch(ones, nglacier, cism_glacier_id, parallel)
    glacier_field = parallel_global_sum_patch(field_2d, nglacier, cism_glacier_id, parallel)

    where (ncells_glacier > 0)
       glacier_field = glacier_field/ncells_glacier
    endwhere

  end subroutine glacier_2d_to_1d

!****************************************************

  subroutine glacier_2d_to_1d_weighted(&
       ewn,           nsn,              &
       parallel,                        &
       nglacier,                        &
       glacier_id,    weight,           &
       field_2d,      glacier_field)

    ! Given a 2D field, compute the average of the field over each glacier
    ! Certain grid cells (e.g., at the glacier periphery) can be given weights between 0 and 1.

    use cism_parallel, only: parallel_global_sum_patch

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    type(parallel_type), intent(in) :: &
         parallel                       ! info for parallel communication

    integer, dimension(ewn,nsn), intent(in) ::  &
         glacier_id                     ! integer glacier ID

    real(dp), dimension(ewn,nsn), intent(in) :: &
         weight                         ! weighting factor applied to each grid cell

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         field_2d                       ! 2D field to be averaged over glaciers

    real(dp), dimension(nglacier), intent(out) ::  &
         glacier_field                  ! field average over each glacier

    ! local variables

    real(dp), dimension(nglacier) :: sum_weights

    sum_weights = parallel_global_sum_patch(weight, nglacier, glacier_id, parallel)
    glacier_field = parallel_global_sum_patch(weight*field_2d, nglacier, glacier_id, parallel)

    where (sum_weights > 0.0d0)
       glacier_field = glacier_field/sum_weights
    endwhere

  end subroutine glacier_2d_to_1d_weighted

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
       parallel,                         &
       nglacier,      cism_glacier_id,   &
       cell_area,     thck,              &
       diagnostic_minthck,               &
       area,          volume)

    use cism_parallel, only: parallel_global_sum_patch

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    type(parallel_type), intent(in) :: &
         parallel                       ! info for parallel communication

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         cell_area,                   & ! grid cell area (m^2)
                                        ! Note: can be latitude-dependent and differ from dew*dns
         thck                           ! ice thickness (m)

    real(dp), intent(in) :: &
         diagnostic_minthck             ! minimum thickness (m) to be included in area and volume sums

    real(dp), dimension(nglacier), intent(out) ::  &
         area,                        & ! area of each glacier (m^2)
         volume                         ! volume of each glacier (m^3)

    ! local variables

    real(dp), dimension(ewn,nsn) ::  &
         diag_area, diag_volume       ! area and volume where thck >= diagnostic_minthck

    ! Compute the area and volume of each glacier.
    ! Need parallel sums, since a glacier can lie on two or more processors.

    where(thck >= diagnostic_minthck)
       diag_area = cell_area
       diag_volume = cell_area*thck
    elsewhere
       diag_area = 0.0d0
       diag_volume = 0.0d0
    endwhere

    area = parallel_global_sum_patch(diag_area, nglacier, cism_glacier_id, parallel)
    volume = parallel_global_sum_patch(diag_volume, nglacier, cism_glacier_id, parallel)

  end subroutine glacier_area_volume

!****************************************************

  subroutine glacier_area_advance_retreat(&
       ewn,           nsn,     &
       parallel,               &
       nglacier,               &
       cism_glacier_id_init,   &
       cism_glacier_id,        &
       cell_area,              &
       area_initial,           &
       area_current,           &
       area_advance,           &
       area_retreat)

    use cism_parallel, only: parallel_global_sum_patch

    ! For each glacier, compare the current glacier area (as given by cism_glacier_id)
    ! to the initial area (given by cism_glacier_id_init).
    ! Compute the area of the advanced region (ice is present now, but not at init)
    ! and the retreated region (ice was present at init, but not now).
    ! Note: For this subroutine, the area is based on the cism_glacier_id masks,
    !       so it includes cells with thck < diagnostic_min_thck.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    type(parallel_type), intent(in) :: &
         parallel                       ! info for parallel communication

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id_init,        & ! integer glacier ID in the range (1, nglacier), initial value
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier), current value

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         cell_area                      ! grid cell area (m^2)

    real(dp), dimension(nglacier), intent(out) ::  &
         area_initial,                & ! initial glacier area
         area_current,                & ! current glacier area
         area_advance,                & ! area of the region where the glacier has advanced (m^2)
         area_retreat                   ! area of the region where the glacier has retreated (m^2)

    ! local variables

    integer, dimension(ewn,nsn) :: glacier_id    ! temporary glacier ID

    integer :: ng

    ! Compute the area of each glacier over the initial and current masks.
    ! We need parallel sums, since a glacier can lie on two or more processors.

    area_initial = parallel_global_sum_patch(cell_area, nglacier, cism_glacier_id_init, parallel)

    ! current area

    area_current = parallel_global_sum_patch(cell_area, nglacier, cism_glacier_id, parallel)

    ! area where the glacier has advanced

    where (cism_glacier_id_init == 0 .and. cism_glacier_id > 0)
       glacier_id = cism_glacier_id
    elsewhere
       glacier_id = 0
    endwhere

    area_advance = parallel_global_sum_patch(cell_area, nglacier, glacier_id, parallel)

    ! area where the glacier has retreated

    where (cism_glacier_id_init > 0 .and. cism_glacier_id == 0)
       glacier_id = cism_glacier_id_init
    elsewhere
       glacier_id = 0
    endwhere

    area_retreat = parallel_global_sum_patch(cell_area, nglacier, glacier_id, parallel)

    ! bug check
    do ng = 1, nglacier
       if (area_initial(ng) + area_advance(ng) - area_retreat(ng) /= area_current(ng)) then
          write(iulog,*) ' '
          write(iulog,*) 'WARNING: area mismatch in glacier_area_advance_retreat'
          write(iulog,*) '   ng, initial, advance, retreat, current:', ng, area_initial(ng)/1.d6, &
               area_advance(ng)/1.d6, area_retreat(ng)/1.d6, area_current(ng)/1.d6
       endif
    enddo

  end subroutine glacier_area_advance_retreat

!****************************************************

  subroutine glacier_accumulation_area_ratio(&
       ewn,           nsn,     &
       parallel,               &
       nglacier,               &
       cism_glacier_id,        &
       cell_area,              &
       smb,                    &
       aar)

    ! Compute the accumulation area ratio (AAR) for each glacier.

    use cism_parallel, only: parallel_global_sum_patch

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    type(parallel_type), intent(in) :: &
         parallel                       ! info for parallel communication

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         cell_area                      ! grid cell area = dew*dns (m^2); same for all cells

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         smb                            ! surface mass balance (mm/yr w.e.)

    real(dp), dimension(nglacier), intent(out) ::  &
         aar                            ! accumulation area ratio

    ! local variables

!    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         ablat_area,                  & ! area of accumulation zone (SMB < 0)
         accum_area                     ! area of accumulation zone (SMB > 0)

    integer, dimension(ewn,nsn) :: glacier_id   ! temporary glacier ID

    ! Compute the accumulation and ablation area for each glacier
    ! Note: Grid cells with SMB = 0 are not counted in either zone.

    where (cism_glacier_id > 0 .and. smb > 0.0d0)
       glacier_id = cism_glacier_id
    elsewhere
       glacier_id = 0
    endwhere

    accum_area = parallel_global_sum_patch(cell_area, nglacier, glacier_id, parallel)

    where (cism_glacier_id > 0 .and. smb < 0.0d0)
       glacier_id = cism_glacier_id
    elsewhere
       glacier_id = 0
    endwhere

    ablat_area = parallel_global_sum_patch(cell_area, nglacier, glacier_id, parallel)

    ! Compute the AAR for each glacier

    where (accum_area + ablat_area > 0.0d0)
       aar = accum_area / (accum_area + ablat_area)
    elsewhere
       aar = 0.0d0
    endwhere

  end subroutine glacier_accumulation_area_ratio

  !****************************************************

  subroutine glacier_smb_min_max(&
       ewn,           nsn,            &
       nglacier,                      &
       cism_glacier_id,               &
       smb,                           &
       smb_min,       smb_max)

    use cism_parallel, only: parallel_reduce_min, parallel_reduce_max

    ! Find the most negative SMB in the glacier.
    ! Typically, this is the SMB in the grid cell with the lowest elevation.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in)  :: &
         cism_glacier_id                ! current cism glacier_id, > 0 for glaciated cells

    real(dp), dimension(ewn,nsn), intent(in) :: &
         smb                            ! surface mass balance (mm/yr w.e.)

    real(dp), dimension(nglacier), intent(out) :: &
         smb_min, smb_max               ! min and max SMB for each glacier (mm/yr w.e.)

    ! local variables

    integer :: i, j, ng

    smb_min(:) = 0.0d0
    smb_max(:) = 0.0d0

    ! Find the most negative SMB for each glacier on the local processor
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0) then
             if (smb(i,j) < smb_min(ng)) then
                smb_min(ng) = smb(i,j)
             endif
             if (smb(i,j) > smb_max(ng)) then
                smb_max(ng) = smb(i,j)
             endif
          endif
       enddo
    enddo

    ! global reductions
    smb_min = parallel_reduce_min(smb_min)
    smb_max = parallel_reduce_max(smb_max)

  end subroutine glacier_smb_min_max

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
