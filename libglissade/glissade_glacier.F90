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
    use glimmer_paramets, only: thk0, len0, tim0, vel0, eps08
    use glimmer_physcon, only: scyr, pi, rhow, rhoi
    use glide_types
    use glimmer_log
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

    use cism_parallel, only: distributed_gather_var, distributed_scatter_var, &
         parallel_reduce_sum, parallel_reduce_max, parallel_reduce_min, &
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
    integer :: ng_west, ng_east, ng_south, ng_north
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

    if (glacier%scale_area) then

       ! Optionally, rescale the grid cell dimensions dew and dns
       ! This is answer-changing throughout the code.
       if (glacier%length_scale_factor /= 1.0d0) then
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
          print*, 'Scale dew and dns: factor, new dew, dns =', &
               glacier%length_scale_factor, dew*len0, dns*len0
          print*, 'Scale cell area: i, j, lat, cos(lat), cell_area =', &
               i, j, model%general%lat(i,j), cos(theta_rad), model%geometry%cell_area(i,j)*len0**2
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
       if (associated(glacier%cism_to_rgi_glacier_id)) &
            deallocate(glacier%cism_to_rgi_glacier_id)
       if (associated(glacier%area)) deallocate(glacier%area)
       if (associated(glacier%volume)) deallocate(glacier%volume)
       if (associated(glacier%area_init)) deallocate(glacier%area_init)
       if (associated(glacier%volume_init)) deallocate(glacier%volume_init)
       if (associated(glacier%area_init_extent)) deallocate(glacier%area_init_extent)
       if (associated(glacier%volume_init_extent)) deallocate(glacier%volume_init_extent)
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
       allocate(glacier%smb(nglacier))
       allocate(glacier%smb_obs(nglacier))
       allocate(glacier%mu_star(nglacier))
       allocate(glacier%alpha_snow(nglacier))
       allocate(glacier%beta_artm(nglacier))

       ! Compute the initial area and volume of each glacier.
       ! Only ice thicker than diagnostic_minthck is included in area and volume sums.

       call glacier_area_volume(&
            ewn,           nsn,               &
            nglacier,                         &
            glacier%cism_glacier_id,          &
            model%geometry%cell_area*len0**2, &  ! m^2
            model%geometry%thck*thk0,         &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area,                     &  ! m^2
            glacier%volume)                      ! m^3

       ! Initialize other glacier arrays
       glacier%smb(:)         = 0.0d0
       glacier%area_init(:)   = glacier%area(:)
       glacier%volume_init(:) = glacier%volume(:)
       glacier%area_init_extent(:) = glacier%area(:)
       glacier%volume_init_extent(:) = glacier%volume(:)
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
       !  and initialize the inversion target to the initial usrf.
       ! Note: usrf_obs is the thickness (in scaled model units) at the RGI date, e.g. the
       !        Farinotti et al. consensus thickness.
       !       usrf_target_baseline is the target thickness for the baseline state, which
       !        ideally will evolve to usrf_obs between the baseline date and RGI date.
       ! On restart, powerlaw_c and usrf_obs are read from the restart file;
       !       usrf_target_baseline is not needed for exact restart.

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          model%basal_physics%powerlaw_c(:,:) = model%basal_physics%powerlaw_c_const
          glacier%usrf_target_baseline(:,:) = model%geometry%usrf(:,:)*thk0
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
             if (ng == 0 .and. model%geometry%thck(i,j)*thk0 > 1.0d0) then
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                print*, 'Warning, ng = 0 but H > 0: Init rank, i, j, ig, jg, thck:', &
                     this_rank, i, j, iglobal, jglobal, model%geometry%thck(i,j) * thk0
             endif
          enddo
       enddo

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
           glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
          ! Make sure a nonzero smb_obs field was read in
          max_glcval = maxval(abs(model%climate%smb_obs))
          max_glcval = parallel_reduce_max(max_glcval)
          if (max_glcval == 0.0d0) then
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

       ! Compute the initial area and volume of each glacier.
       ! This is not necessary for exact restart, but is included as a diagnostic.
       ! Only ice thicker than diagnostic_minthck is included in area and volume sums.

       call glacier_area_volume(&
            ewn,           nsn,               &
            nglacier,                         &
            glacier%cism_glacier_id,          &
            model%geometry%cell_area*len0**2, &  ! m^2
            model%geometry%thck*thk0,         &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area,                     &  ! m^2
            glacier%volume)                      ! m^3

       ! Compute the area and volume over the initial ice extent.

       call glacier_area_volume(&
            ewn,           nsn,               &
            nglacier,                         &
            glacier%cism_glacier_id_init,     &
            model%geometry%cell_area*len0**2, &  ! m^2
            model%geometry%thck*thk0,         &  ! m
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

    glacier%minthck = model%numerics%thklim*thk0 - eps08

    ! Set the relaxation value for powerlaw_c
    if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
       model%basal_physics%powerlaw_c_relax(:,:) = model%basal_physics%powerlaw_c_const
    endif

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

    ! Define a mask whose value is 1 at vertices along the boundary between two glaciers.
    ! At runtime, Cp is set to a large value at masked vertices to reduce flow between glaciers.
    glacier%boundary_mask(:,:) = 0

    ! Loop over locally owned cells
    do j = nhalo, nsn-nhalo
       do i = nhalo, ewn-nhalo
          ng = glacier%cism_glacier_id_init(i,j)
          if (ng > 0) then
             ng_west  = glacier%cism_glacier_id_init(i-1,j)
             ng_east  = glacier%cism_glacier_id_init(i+1,j)
             ng_south = glacier%cism_glacier_id_init(i,j-1)
             ng_north = glacier%cism_glacier_id_init(i,j+1)
             if (ng_west > 0 .and. ng_west /= ng) then
                glacier%boundary_mask(i-1,j-1) = 1
                glacier%boundary_mask(i-1,j)   = 1
             endif
             if (ng_east > 0 .and. ng_east /= ng) then
                glacier%boundary_mask(i,j-1) = 1
                glacier%boundary_mask(i,j)   = 1
             endif
             if (ng_south > 0 .and. ng_south /= ng) then
                glacier%boundary_mask(i-1,j-1) = 1
                glacier%boundary_mask(i,j-1)   = 1
             endif
             if (ng_north > 0 .and. ng_north /= ng) then
                glacier%boundary_mask(i-1,j) = 1
                glacier%boundary_mask(i,j)   = 1
             endif
          endif
       enddo
    enddo

    call staggered_parallel_halo(glacier%boundary_mask, parallel)

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'Create glacier boundary_mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') &
                  glacier%boundary_mask(i,j)
          enddo
          print*, ' '
       enddo
    endif

    ! Write some values for the diagnostic glacier
    if (verbose_glacier .and. this_rank == rtest) then
       i = itest; j = jtest
       ng = glacier%ngdiag
       print*, ' '
       print*, 'Glacier ID for diagnostic cell: r, i, j, ng =', rtest, itest, jtest, ng
       if (ng > 0) then
          print*, 'area_init (km^2) =', glacier%area_init(ng) / 1.0d6
          print*, 'volume_init (km^3) =', glacier%volume_init(ng) / 1.0d9
          print*, 'powerlaw_c (Pa (m/yr)^(-1/3)) =', model%basal_physics%powerlaw_c(i,j)
          print*, 'smb_obs (mm/yr w.e.) =', glacier%smb_obs(ng)
          print*, 'mu_star (mm/yr w.e./deg) =', glacier%mu_star(ng)
          print*, 'Done in glissade_glacier_init'
       endif
    endif

  end subroutine glissade_glacier_init

!****************************************************

  subroutine glissade_glacier_update(model, glacier)

    use glissade_grid_operators, only: glissade_stagger
    use glissade_utils, only: glissade_usrf_to_thck
    use cism_parallel, only: parallel_reduce_sum, parallel_global_sum, parallel_halo

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
         ice_mask                   ! = 1 where ice is present (thck > thklim), else = 0

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck,                    & ! ice thickness (m)
         thck_target,             & ! target ice thickness for the baseline state (m)
         dthck_dt,                & ! rate of change of thickness (m/yr)
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
    ! real(dp), dimension(:) :: volume_init_extent! current glacier volume (m^3) over initial ice extent
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
          ! Then, provided usrf is close to usrf_target_baseline in the spin-up, usrf will be close to
          !  usrf_obs (the RGI target) when a forward run starting from the baseline date reaches the RGI date.

          glacier%usrf_target_baseline(:,:) = &
               model%geometry%usrf_obs(:,:)*thk0 - glacier%delta_usrf_rgi(:,:)

          ! Make sure the target is not below the topography
          glacier%usrf_target_baseline = &
               max(glacier%usrf_target_baseline, (model%geometry%topg + model%climate%eus)*thk0)

          if (verbose_glacier .and. this_rank == rtest) then
             i = itest; j = jtest
             print*, ' '
             print*, 'RGI usrf correction, delta_smb:', &
                  glacier%delta_usrf_rgi(i,j), delta_smb_rgi(i,j)
             print*,    'usrf RGI obs, new usrf_target_baseline =', &
                  model%geometry%usrf_obs(i,j)*thk0, glacier%usrf_target_baseline(i,j)
             print*, 'Recent usrf correction, delta_smb:', &
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
       print*, ' '
       print*, 'glissade_glacier_update, diag cell (r, i, j) =', rtest, itest, jtest
       print*, ' '
       ! Convert acab_applied from m/yr ice to mm/yr w.e.
       write(6,'(a32,2f10.3)') '     acab_applied, smb_applied: ', &
            model%climate%acab_applied(i,j)*scyr*thk0/tim0, &  ! m/yr ice
            model%climate%acab_applied(i,j)*scyr*thk0/tim0 * 1000.d0*(rhoi/rhow)  ! mm/yr w.e.
       write(6,'(a32,4f10.3)') 'artm_ref, usrf_ref, usrf, diff: ', &
            model%climate%artm_ref(i,j), &
            model%climate%usrf_ref(i,j), model%geometry%usrf(i,j)*thk0, &
            model%geometry%usrf(i,j)*thk0 - model%climate%usrf_ref(i,j)
       write(6,'(a32,3f10.3)') '              artm, Tpos, snow: ', artm(i,j), Tpos(i,j), snow(i,j)
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
          write(6,'(a32,3f10.3)') '         RGI artm, Tpos, snow: ', &
               artm_rgi(i,j), Tpos_rgi(i,j), snow_rgi(i,j)
          write(6,'(a32,3f10.3)') '      Recent artm, Tpos, snow: ', &
               artm_recent(i,j), Tpos_recent(i,j), snow_recent(i,j)
       endif

    endif   ! set_mu_star

    ! Accumulate snow_annmean, Tpos_annmean, and dthck_dt_annmean over this timestep

    time_since_last_avg = time_since_last_avg + dt

    glacier%snow_annmean = glacier%snow_annmean + snow * dt
    glacier%Tpos_annmean = glacier%Tpos_annmean + Tpos * dt
    glacier%smb_applied_annmean = glacier%smb_applied_annmean  &
         + model%climate%acab_applied*(scyr*thk0/tim0) * 1000.d0*(rhoi/rhow) * dt

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
          print*, ' '
          print*, 'Annual averages, r, i, j:', rtest, itest, jtest
          print*, '   snow (mm/yr)       =', glacier%snow_annmean(i,j)
          print*, '   Tpos (deg C)       =', glacier%Tpos_annmean(i,j)
          print*, '   smb_applied (mm/yr)=', glacier%smb_applied_annmean(i,j)
          if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
              glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
             print*, '   snow_rgi (mm/yr)   =', glacier%snow_rgi_annmean(i,j)
             print*, '   Tpos_rgi (deg C)   =', glacier%Tpos_rgi_annmean(i,j)
             print*, '   snow_rec (mm/yr)   =', glacier%snow_recent_annmean(i,j)
             print*, '   Tpos_rec (deg C)   =', glacier%Tpos_recent_annmean(i,j)
          endif
          if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
             print*, '   dthck_dt (m/yr)    =', glacier%dthck_dt_annmean(i,j)
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

          thck_old = thck

          call glacier_redistribute_advanced_ice(&
            ewn,             nsn,               &
            itest,   jtest,  rtest,             &
            nglacier,        ngdiag,            &
            real(glacier_update_interval,dp),   &  ! yr
            dew*dns,                            &  ! m^2
            glacier%thinning_rate_advanced_ice, &  ! m/yr
            glacier%cism_glacier_id_init,       &
            glacier%smb_glacier_id,             &
            model%climate%smb,                  &  ! m/yr
            thck,                               &  ! m
            parallel)

          glacier%dthck_dt_annmean = glacier%dthck_dt_annmean + &
               (thck - thck_old) / real(glacier_update_interval,dp)

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
       ! Note: This subroutine assumes cell_area = dew*dns for all cells
       call glacier_area_advance_retreat(&
            ewn,           nsn,           &
            nglacier,                     &
            glacier%cism_glacier_id_init, &
            glacier%cism_glacier_id,      &
            dew*dns,                      &
            area_initial,                 &
            area_current,                 &
            area_advance,                 &
            area_retreat)

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'Selected big glaciers:'
          write(6,'(a101)') &
               '  ng,   Ainit,     A,     Vinit,     V,   smb_iniA, smb_curA, mu_star, alpha_snow, beta_artm, smb_obs'
          do ng = 1, nglacier
             if (glacier%volume_init(ng) > diagnostic_volume_threshold .or. ng == ngdiag) then  ! big glacier
                write(6,'(i6,4f9.3,6f10.3)') ng, glacier%area_init(ng)/1.e6, glacier%area(ng)/1.e6, &
                     glacier%volume_init(ng)/1.0d9, glacier%volume(ng)/1.0d9, &
                     smb_init_area(ng), smb_current_area(ng), glacier%mu_star(ng), glacier%alpha_snow(ng), &
                     glacier%beta_artm(ng), glacier%smb_obs(ng)
             endif
          enddo
       endif

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'Advance/retreat diagnostics'
          print*, '  ng  A_initial A_advance A_retreat A_current'
          do ng = 1, nglacier
             if (glacier%volume_init(ng) > 1.0d9 .or. ng == ngdiag) then  ! big glacier, > 1 km^3
                write(6,'(i6,6f10.3)') ng, area_initial(ng)/1.e6, area_advance(ng)/1.e6, &
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
               glacier%usrf_target_baseline,    &
               model%geometry%topg * thk0,      &
               model%climate%eus * thk0,        &
               thck_target)

          ! Interpolate thck_target to the staggered grid
          call glissade_stagger(&
               ewn,         nsn,              &
               thck_target, stag_thck_target)

          ! Interpolate thck to the staggered grid
          call glissade_stagger(&
               ewn,         nsn,              &
               thck,        stag_thck)

          ! Interpolate dthck_dt to the staggered grid
          call glissade_stagger(&
               ewn,                      nsn,           &
               glacier%dthck_dt_annmean, stag_dthck_dt)

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
               model%basal_physics%powerlaw_c_relax,    &
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

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'topg:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%topg(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'Before advance_retreat, thck, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i4)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       ! Assign nonzero IDs in grid cells where ice has reached the minimum glacier thickness.
       ! Remove IDs in grid cells where ice is now thinnier than the minimum thickness.
       ! Adjust IDs to prevent spurious advance due to SMB differences in adjacent glaciers.

       call glacier_advance_retreat(&
            ewn,             nsn,           &
            itest,   jtest,  rtest,         &
            nglacier,                       &
            glacier%minthck,                &  ! m
            thck,                           &  ! m
            glacier%snow_annmean,           &  ! mm/yr w.e.
            glacier%Tpos_annmean,           &  ! deg C
            glacier%mu_star,                &  ! mm/yr/deg
            glacier%alpha_snow,             &  ! unitless
            glacier%cism_glacier_id_init,   &
            glacier%cism_glacier_id,        &
            parallel)

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

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'After advance_retreat, thck, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i4)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'cism_glacier_id_init:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i11)',advance='no') glacier%cism_glacier_id_init(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'smb_glacier_id_init:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i11)',advance='no') glacier%smb_glacier_id_init(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'New cism_glacier_id:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i11)',advance='no') glacier%cism_glacier_id(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'New smb_glacier_id:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i11)',advance='no') glacier%smb_glacier_id(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'smb_applied_annmean (previous year):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f11.3)',advance='no') glacier%smb_applied_annmean(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'smb_weight_init (previous year):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f11.3)',advance='no') smb_weight_init(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'Tpos_annmean:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f11.3)',advance='no') glacier%Tpos_annmean(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'snow_annmean:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f11.3)',advance='no') glacier%snow_annmean(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'model%climate%smb:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f11.3)',advance='no') model%climate%smb(i,j)
             enddo
             print*, ' '
          enddo
          if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
              glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
             print*, ' '
             print*, 'smb_rgi:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f11.3)',advance='no') glacier%smb_rgi(i,j)
                enddo
                print*, ' '
             enddo
             print*, ' '
             print*, 'smb_recent:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f11.3)',advance='no') glacier%smb_recent(i,j)
                enddo
                print*, ' '
             enddo
          endif   ! set_mu_star
       endif   ! verbose

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
            nglacier,                         &
            glacier%cism_glacier_id_init,     &
            model%climate%smb,                &
            aar_init)

       ! (2) Include all cells in the glacier
       call glacier_accumulation_area_ratio(&
            ewn,           nsn,               &
            nglacier,                         &
            glacier%cism_glacier_id,          &
            model%climate%smb,                &
            aar)

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'Glacier SMB and AAR:'
          print*, '    ng     smb_min   smb_max   AAR_initA     AAR'
          do ng = 1, nglacier
             if (glacier%volume_init(ng) > diagnostic_volume_threshold .or. ng == ngdiag) then  ! big glacier
                write(6,'(i10, 2f10.1, 2f10.4 )') ng, smb_min(ng), smb_max(ng), aar_init(ng), aar(ng)
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
               nglacier,                         &
               glacier%cism_glacier_id_init,     &
               glacier%smb_recent,               &
               aar_init_recent)

          ! (2) Include all cells in the glacier
          call glacier_accumulation_area_ratio(&
               ewn,           nsn,               &
               nglacier,                         &
               glacier%cism_glacier_id,          &
               glacier%smb_recent,               &
               aar_recent)

          if (verbose_glacier .and. this_rank == rtest) then
             print*, ' '
             print*, 'Recent SMB and AAR:'
             print*, '    ng     smb_min   smb_max   AAR_initA     AAR'
             do ng = 1, nglacier
                if (glacier%volume_init(ng) > diagnostic_volume_threshold .or. ng == ngdiag) then  ! big glacier
                   write(6,'(i10, 2f10.1, 2f10.4 )') ng, smb_min_recent(ng), smb_max_recent(ng), &
                        aar_init_recent(ng), aar_recent(ng)
                endif
             enddo
          endif

       endif   ! set_mu_star

       ! Update the glacier area and volume (diagnostic only)

       ! Compute the new area and volume

       call glacier_area_volume(&
            ewn,           nsn,               &
            nglacier,                         &
            glacier%cism_glacier_id,          &
            model%geometry%cell_area*len0**2, &  ! m^2
            thck,                             &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area,                     &  ! m^2
            glacier%volume)                      ! m^3

       ! Compute the new area and volume over the initial ice extent
       ! Note: area_init_extent <= area_init; inequality applies if there has been any retreat

       call glacier_area_volume(&
            ewn,           nsn,               &
            nglacier,                         &
            glacier%cism_glacier_id_init,     &
            model%geometry%cell_area*len0**2, &  ! m^2
            thck,                             &  ! m
            glacier%diagnostic_minthck,       &  ! m
            glacier%area_init_extent,         &  ! m^2
            glacier%volume_init_extent)          ! m^3

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'Update area (km^2) and volume (km^3) for glacier:', ngdiag
          print*, ' Initial area and volume:', &
               glacier%area_init(ngdiag)/1.0d6, glacier%volume_init(ngdiag)/1.0d9
          print*, ' Current area and volume:', &
               glacier%area(ngdiag)/1.0d6, glacier%volume(ngdiag)/1.0d9
          print*, 'A and V over init extent:', &
               glacier%area_init_extent(ngdiag)/1.0d6, glacier%volume_init_extent(ngdiag)/1.0d9
       endif

       if (verbose_glacier) then

          ! debug - count cells in masks
          count_cgii = 0
          count_cgi = 0
          count_sgii = 0
          count_sgi = 0
          do j = nhalo+1, nsn-nhalo
             do i = nhalo+1, ewn-nhalo
                ng = glacier%cism_glacier_id_init(i,j)
                if (ng == ngdiag) count_cgii = count_cgii + 1
                ng = glacier%cism_glacier_id(i,j)
                if (ng == ngdiag) count_cgi  = count_cgi  + 1
                ng = glacier%smb_glacier_id_init(i,j)
                if (ng == ngdiag) count_sgii = count_sgii + 1
                ng = glacier%smb_glacier_id(i,j)
                if (ng == ngdiag) count_sgi  = count_sgi  + 1
             enddo
          enddo

          count_cgii = parallel_reduce_sum(count_cgii)
          count_cgi  = parallel_reduce_sum(count_cgi)
          count_sgii = parallel_reduce_sum(count_sgii)
          count_sgi  = parallel_reduce_sum(count_sgi)

          if (this_rank == rtest) then
             print*, ' '
             print*, 'Mask count, ng =', ngdiag
             print*, 'count_cgii, count_cgi =', count_cgii, count_cgi
             print*, 'count_sgii, count_sgi =', count_sgii, count_sgi
          endif

       endif   ! verbose

    endif   ! glacier_update_inverval

    ! Convert fields back to dimensionless units as needed
    model%geometry%thck = thck/thk0

  end subroutine glissade_glacier_update

!****************************************************

  subroutine glacier_invert_mu_star(&
       ewn,              nsn,           &
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
       print*, ' '
       print*, 'In glacier_invert_mu_star'
    endif

    ! Compute weighted averages of Tpos and snow over each glacier

    call glacier_2d_to_1d_weighted(&
         ewn,           nsn,                   &
         nglacier,                             &
         smb_glacier_id_init,                  &
         smb_weight,                           &
         snow,          glacier_snow)

    call glacier_2d_to_1d_weighted(&
         ewn,           nsn,                   &
         nglacier,                             &
         smb_glacier_id_init,                  &
         smb_weight,                           &
         Tpos,          glacier_Tpos)

    if (verbose_glacier .and. this_rank == rtest) then
       ng = ngdiag
       print*, ' '
       print*, 'ng, snow and Tpos with weighting =', ng, glacier_snow(ng), glacier_Tpos(ng)
    endif

    ! For each glacier, compute the new mu_star. Adjust beta_artm if necessary.

    do ng = 1, nglacier

       if (glacier_snow(ng) == 0.0d0) then

          if (verbose_glacier .and. this_rank == rtest) then
             print*, 'WARNING: snow = 0 for glacier', ng
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
             print*, 'WARNING, mu out of range: ng, mu =', ng, mu_star(ng)
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
!!             write(6,'(a60,i6,6f10.2)') 'Eq 1 violation, ng, snow, Tpos, init mu, adj mu, beta, smb :', &
!!                  ng, glacier_snow(ng), glacier_Tpos(ng), mu_eq1, mu_star(ng), beta_artm(ng), smb_baseline(ng)
             count_violate_1 = count_violate_1 + 1
             area_violate_1 = area_violate_1 + glacier_area_init(ng)
             volume_violate_1 = volume_violate_1 + glacier_volume_init(ng)
          endif
       endif

    enddo  ! ng

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'Violations of Eq. 1 (SMB = 0, baseline climate):', count_violate_1
       print*, '   Total area, volume =', area_violate_1/1.0d6, volume_violate_1/1.0d9
       print*, ' '
       ng = ngdiag
       print*, 'Balance solution, ng =', ng
       write(6,'(a30,3f12.4)') '   mu_star, alpha_snow, beta: ', &
            mu_star(ng), alpha_snow(ng), beta_artm(ng)
       write(6,'(a30,3f12.4)') '   Baseline snow, Tpos, SMB : ', &
            glacier_snow(ng), glacier_Tpos(ng), smb_baseline(ng)
    endif

  end subroutine glacier_invert_mu_star

!****************************************************

  subroutine glacier_invert_mu_star_alpha_snow(&
       ewn,              nsn,            &
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
       print*, ' '
       print*, 'In glacier_invert_mu_star_alpha_snow'
    endif

    ! Compute weighted averages of Tpos and snow over each glacier

    call glacier_2d_to_1d_weighted(&
         ewn,                  nsn,       &
         nglacier,                        &
         smb_glacier_id_init,             &
         smb_weight,                      &
         snow,                 glacier_snow)

    call glacier_2d_to_1d_weighted(&
         ewn,                  nsn,       &
         nglacier,                        &
         smb_glacier_id_init,             &
         smb_weight,                      &
         Tpos,                 glacier_Tpos)

    call glacier_2d_to_1d_weighted(&
         ewn,                  nsn,       &
         nglacier,                        &
         smb_glacier_id_init,             &
         smb_weight,                      &
         snow_recent,          glacier_snow_recent)

    call glacier_2d_to_1d_weighted(&
         ewn,                  nsn,       &
         nglacier,                        &
         smb_glacier_id_init,             &
         smb_weight,                      &
         Tpos_recent,          glacier_Tpos_recent)

    if (verbose_glacier .and. this_rank == rtest) then
       ng = ngdiag
       print*, ' '
       print*, 'ng, snow and Tpos with weighting =', ng, glacier_snow(ng), glacier_Tpos(ng)
       print*, 'recent snow and Tpos with weighting =', glacier_snow_recent(ng), glacier_Tpos_recent(ng)
    endif

    ! For each glacier, compute the new mu_star and alpha_snow

    do ng = 1, nglacier

       if (glacier_snow(ng) == 0.0d0) then

          if (verbose_glacier .and. this_rank == rtest) then
             print*, 'WARNING: snow = 0 for glacier', ng
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
!                      write(6,'(a46,i6,6f10.3)') 'Out of range, ng, Tp, Tp_recent, D, B, alpha, mu:', &
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
             print*, 'WARNING, mu out of range: ng, mu =', ng, mu_star(ng)
          endif
       endif

       if (alpha_snow(ng) < alpha_snow_min .or. alpha_snow(ng) > alpha_snow_max) then
          if (this_rank == rtest) then
             print*, 'WARNING, alpha out of range: ng, alpha =', ng, alpha_snow(ng)
          endif
       endif

!       if (abs(beta_artm(ng)) > beta_artm_max) then
!          if (this_rank == rtest) then
!             print*, 'WARNING, beta out of range: ng, beta =', ng, beta_artm(ng)
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
!!             write(6,'(a60,i6,6f10.2)') 'Eq 1 violation, ng, snow, Tpos, init mu, adj mu, beta, smb :', &
!!                  ng, glacier_snow(ng), glacier_Tpos(ng), mu_eq1, mu_star(ng), beta_artm(ng), smb_baseline(ng)
             count_violate_1 = count_violate_1 + 1
             area_violate_1 = area_violate_1 + glacier_area_init(ng)
             volume_violate_1 = volume_violate_1 + glacier_volume_init(ng)
          endif
          if (abs(smb_recent_diff(ng)) > eps08) then
!!             print*, '   Violation of Eq. 2: ng, smb_recent_diff =', ng, smb_recent_diff(ng)
             count_violate_2 = count_violate_2 + 1
             area_violate_2 = area_violate_2 + glacier_area_init(ng)
             volume_violate_2 = volume_violate_2 + glacier_volume_init(ng)
          endif
       endif

    enddo  ! ng

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'Violations of Eq. 1 (SMB = 0, baseline climate):', count_violate_1
       print*, '   Total area, volume =', area_violate_1/1.0d6, volume_violate_1/1.0d9
       print*, 'Violations of Eq. 2 (SMB = SMB_obs, recent climate):', count_violate_2
       print*, '   Total area, volume =', area_violate_2/1.0d6, volume_violate_2/1.0d9
       print*, ' '
       ng = ngdiag
       print*, 'Balance solution, ng =', ng
       write(6,'(a30,3f12.4)') '   mu_star, alpha_snow, beta: ', &
            mu_star(ng), alpha_snow(ng), beta_artm(ng)
       write(6,'(a30,3f12.4)') '   Baseline snow, Tpos, SMB : ', &
            glacier_snow(ng), glacier_Tpos(ng), smb_baseline(ng)
       write(6,'(a30,3f12.4)') '     Recent snow, Tpos, SMB : ', &
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

    real(dp), intent(in) :: &
         babc_timescale,              & ! inversion timescale for powerlaw_c (yr)
         babc_thck_scale,             & ! inversion thickness scale for powerlaw_c (m)
         babc_relax_factor              ! controls strength of relaxation to default values (unitless)

    real(dp), dimension(ewn-1,nsn-1), intent(in) :: &
         stag_thck,                   & ! ice thickness at vertices (m)
         stag_thck_target,            & ! target ice thickness at vertices (m)
         stag_dthck_dt                  ! rate of change of ice thickness at vertices (m/yr)

    real(dp), dimension(ewn-1,nsn-1), intent(in) :: &
         powerlaw_c_relax               ! powerlaw_c field to which we relax

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

    ! The inversion works as follows:
    ! The change in C_p is proportional to the current value of C_p and to the relative error,
    !  err_H = (H - H_target)/H_scale, where H is a thickness scale.
    ! If err_H > 0, we reduce C_p to make the ice flow faster and thin.
    ! If err_H < 0, we increase C_p to make the ice flow slower and thicken.
    ! This is done with a characteristic timescale tau.
    ! We also include a term proportional to dH/dt so that ideally, C_p smoothly approaches
    !  the value needed to attain a steady-state H, without oscillating about the desired value.
    ! In addition, we include a relaxation term proportional to the ratio of C_p to a default value.
    ! See the comments in module glissade_inversion, subroutine invert_basal_friction.
    !
    ! Here is the prognostic equation:
    ! dC/dt = -C * [(H - H_target)/(H0*tau) + dH/dt * 2/H0 - r * ln(C/C_r) / tau],
    !   where tau = glacier_powerlaw_c_timescale, H0 = glacier_powerlaw_c_thck_scale,
    !         r = glacier_powerlaw_c_relax_factor, and C_r = powerlaw_c_relax.

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glacier_invert_powerlaw_c'
    endif

    if (babc_thck_scale > 0.0d0 .and. babc_timescale > 0.0d0) then

       stag_dthck(:,:) = stag_thck(:,:) - stag_thck_target(:,:)

       ! Loop over vertices

       do j = 1, nsn-1
          do i = 1, ewn-1

             if (stag_thck(i,j) > 0.0d0) then

                term_thck = -stag_dthck(i,j) / (babc_thck_scale * babc_timescale)
                term_dHdt = -stag_dthck_dt(i,j) * 2.0d0 / babc_thck_scale

                ! Add a term to relax C = powerlaw_c toward a target value, C_r = powerlaw_c_relax
                ! The log term below ensures the following:
                ! * When C /= C_r, it will relax toward C_r.
                ! * When C = C_r, there is no further relaxation.
                ! * In steady state (dC/dt = 0, dH/dt = 0), we have dthck/thck_scale = -k * ln(C/C_r),
                !    or C = C_r * exp(-dthck/(k*thck_scale)), where k is a prescribed constant

                term_relax = -babc_relax_factor * log(powerlaw_c(i,j)/powerlaw_c_relax(i,j)) &
                     / babc_timescale

                dpowerlaw_c = powerlaw_c(i,j) * (term_thck + term_dHdt + term_relax) * glacier_update_interval

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
                   print*, 'H, H_target (m)', stag_thck(i,j), stag_thck_target(i,j)
                   print*, 'dH_dt (m/yr):', stag_dthck_dt(i,j)
                   print*, 'dt (yr), term_thck*dt, term_dHdt*dt:', glacier_update_interval, &
                        term_thck*glacier_update_interval, term_dHdt*glacier_update_interval
                   print*, 'relax term:', term_relax*glacier_update_interval
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
       print*, ' '
       print*, 'stag_thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') stag_thck(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'stag_thck - stag_thck_target (m):'
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
       itest,   jtest,  rtest,        &
       nglacier,        ngdiag,       &
       glacier_update_interval,       & ! yr
       cell_area,                     & ! m^2
       thinning_rate_advanced_ice,    & ! m/yr
       cism_glacier_id_init,          &
       smb_glacier_id,                &
       smb,                           & ! m/yr
       thck,                          & ! m
       parallel)

    ! Limit glacier advance in the accumulation zone.
    ! This applies to grid cells that are initially ice-free, into which ice is advected.
    ! The fix here is to thin the ice in these cells at a prescribed rate and
    !  redistribute the mass conservatively across the glacier.

    use cism_parallel, only: parallel_reduce_sum, parallel_halo

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest,         & ! coordinates of diagnostic cell
         nglacier,                    & ! number of glaciers
         ngdiag                         ! CISM ID of diagnostic glacier

    real(dp), intent(in) :: &
         glacier_update_interval,     & ! time interval (yr) of the glacier update, typically 1 yr
         cell_area,                   & ! grid cell area (m^2), assumed to be the same for each cell
         thinning_rate_advanced_ice     ! thinning rate (m/yr) where glaciers advance in the accumulation zone

    integer, dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init,        & ! integer glacier ID at the start of the run
         smb_glacier_id                 ! integer ID for current glacier cells and adjacent glacier-free cells

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         smb                            ! surface mass balance (m/yr)

    real(dp), dimension(ewn,nsn), intent(inout) ::  &
         thck                           ! ice thickness (m)

    type(parallel_type), intent(in) :: parallel   ! info for parallel communication

    ! local variables

    integer :: i, j, ng

    real(dp) :: dthck                   ! thickness change (m)

    real(dp), dimension(nglacier) :: &
         glacier_area_init,           & ! glacier area based on cism_glacier_id_init
         glacier_vol_removed,         & ! total volume (m^3) removed from each advanced cells in each glacier
         glacier_dthck,               & ! thickness (m) added over the initial extent of each glacier
         glacier_vol_1,               & ! volume (m^3) of each glacier before thinning and restribution
         glacier_vol_2                  ! volume (m^3) of each glacier after thinning and restribution

    ! Compute the total volume of each glacier before limiting advance.
    ! Note: This includes adjacent glacier-free cells that might have a small nonzero thickness
    !       (i.e., cism_glacier_id = 0 but smb_glacier_id > 0).
    !TODO: Write a sum-over-glaciers subroutine

    glacier_vol_1(:) = 0.0d0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = smb_glacier_id(i,j)
          if (ng > 0) then
             glacier_vol_1(ng) = glacier_vol_1(ng) + cell_area*thck(i,j)
          endif
       enddo
    enddo
    glacier_vol_1 = parallel_reduce_sum(glacier_vol_1)

    ! compute the area of each glacier over its initial extent
    glacier_area_init(:) = 0.0d0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id_init(i,j)
          if (ng > 0) then
             glacier_area_init(ng) = glacier_area_init(ng) + cell_area
          endif
       enddo
    enddo
    glacier_area_init = parallel_reduce_sum(glacier_area_init)

    ! Compute thinning in advanced grid cells
    ! This includes potential advanced cells adjacent to current glacier cells.
    ! Note: Currently, SMB is set to 0 in advanced cells where SMB would be > 0 otherwise.
    !       The logic below (smb >= 0) ensures that ice in these cells is thinned.

    glacier_vol_removed(:) = 0.0d0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          if (cism_glacier_id_init(i,j) == 0 .and. smb_glacier_id(i,j) > 0) then ! advanced cell
             if (smb(i,j) >= 0.d0) then   ! accumulation zone
                ng = smb_glacier_id(i,j)
                dthck = min(thinning_rate_advanced_ice*glacier_update_interval, thck(i,j))
                thck(i,j) = thck(i,j) - dthck
                glacier_vol_removed(ng) = glacier_vol_removed(ng) + cell_area*dthck
             endif
          endif
       enddo
    enddo
    glacier_vol_removed = parallel_reduce_sum(glacier_vol_removed)

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
    glacier_vol_2(:) = 0.0d0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = max(cism_glacier_id_init(i,j), smb_glacier_id(i,j))
          if (ng > 0) then
             glacier_vol_2(ng) = glacier_vol_2(ng) + cell_area*thck(i,j)
          endif
       enddo
    enddo
    glacier_vol_2 = parallel_reduce_sum(glacier_vol_2)

    ! conservation check
    do ng = 1, nglacier
       if (abs(glacier_vol_2(ng) - glacier_vol_1(ng)) > eps08*glacier_vol_1(ng)) then
          write(6,*) 'redistribute advanced ice, conservation error: ng, vol_1, vol_2:', &
               ng, glacier_vol_1(ng)/1.d9, glacier_vol_2(ng)/1.d9
          call write_log('Volume conservation error, redistribute advanced ice', GM_FATAL)
       endif
    enddo

  end subroutine glacier_redistribute_advanced_ice

  !****************************************************

  subroutine glacier_advance_retreat(&
       ewn,             nsn,         &
       itest,   jtest,  rtest,       &
       nglacier,                     &
       glacier_minthck,              &
       thck,                         &
       snow,                         &
       Tpos,                         &
       mu_star,                      &
       alpha_snow,                   &
       cism_glacier_id_init,         &
       cism_glacier_id,              &
       parallel)

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

    real(dp), intent(in) :: &
         glacier_minthck                ! min ice thickness (m) counted as part of a glacier

    real(dp), dimension(ewn,nsn), intent(inout) ::  &
         thck                           ! ice thickness (m)

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         snow,                        & ! annual mean snowfall (mm/yr w.e.)
         Tpos                           ! annual mean Tpos = min(T - Tmlt, 0)

    real(dp), dimension(nglacier), intent(in) :: &
         mu_star,                     & ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)
         alpha_snow                     ! glacier-specific snow factor (unitless)

    integer, dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init           ! cism_glacier_id at the start of the run

    integer, dimension(ewn,nsn), intent(inout) :: &
         cism_glacier_id                ! current cism glacier_id, > 0 for glaciated cells

    type(parallel_type), intent(in) :: parallel  ! diagnostic only

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
       print*, ' '
       print*, 'In glacier_advance_retreat'
    endif

    ! Check for retreat: cells with cism_glacier_id > 0 but H < glacier_minthck

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0 .and. thck(i,j) <= glacier_minthck) then
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
                   print*, 'Set ID = init ID: ig, jg, new ID, thck =',&
                        iglobal, jglobal, cism_glacier_id(i,j), thck(i,j)
                endif
             else  ! assign the ID of an adjacent glaciated cell, if possible
                found_neighbor = .false.
                smb_min = big_number
                ng_min = 0
                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, 'Look for glaciated neighbor: ig, jg =', iglobal, jglobal
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
                      print*, '  Set ID = neighbor ID, ig, jg, ID, H, smb =', &
                           iglobal, jglobal, cism_glacier_id(i,j), thck(i,j), smb_min
                   endif
                else  ! no adjacent glacier cell
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, '  Warning, did not find neighbor, ig, jg =', iglobal, jglobal
                   print*, '  Setting H = 0'
                   thck(i,j) = 0.0d0  !TODO - anything else to zero out?
                endif   ! found_neighbor

             endif   ! cism_glacier_id_init > 0
          endif   ! ng = 0, H > minthck
       enddo   ! i
    enddo   ! j

    call parallel_halo(cism_glacier_id, parallel)

    ! Check advanced cells (beyond the initial extent) for problematic glacier IDs.
    ! This code protects against glacier 'pirating', which ccan occur when an advanced cell
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
                   print*, '   Transfer to adjacent glacier, old and new IDs =', &
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
         Tpos                           ! annual mean Tpos = min(T - Tmlt, 0)

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
!                   print*, 'Set smb_glacier_id = neighbor ID: ig, jg, smb_min, neighbor ID =', &
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
       nglacier,      cism_glacier_id,  &
       field_2d,      glacier_field)

    ! Given a 2D field, compute the average of the field over each glacier
    !TODO - Pass in cellarea to compute an area average.

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

  subroutine glacier_2d_to_1d_weighted(&
       ewn,           nsn,              &
       nglacier,                        &
       glacier_id,    weight,           &
       field_2d,      glacier_field)

    ! Given a 2D field, compute the average of the field over each glacier
    ! Certain grid cells (e.g., at the glacier periphery) can be given weights between 0 and 1.

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         glacier_id                     ! integer glacier ID

    real(dp), dimension(ewn,nsn), intent(in) :: &
         weight                         ! weighting factor applied to each grid cell

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         field_2d                       ! 2D field to be averaged over glaciers

    real(dp), dimension(nglacier), intent(out) ::  &
         glacier_field                  ! field average over each glacier

    ! local variables

    integer :: i, j, ng

    real(dp), dimension(nglacier) :: sum_weights

    sum_weights(:) = 0.0d0
    glacier_field(:) = 0.0d0

    ! Loop over locally owned cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = glacier_id(i,j)
          if (ng > 0) then
             sum_weights(ng) = sum_weights(ng) + weight(i,j)
             glacier_field(ng) = glacier_field(ng) + weight(i,j) * field_2d(i,j)
          endif
       enddo
    enddo

    sum_weights = parallel_reduce_sum(sum_weights)
    glacier_field  = parallel_reduce_sum(glacier_field)
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
       nglacier,      cism_glacier_id,   &
       cell_area,     thck,              &
       diagnostic_minthck,               &
       area,          volume)

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

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

    real(dp), dimension(nglacier) ::  &
         local_area, local_volume       ! area and volume on each processor, before global sum

    integer :: i, j, ng

    ! Initialize the output arrays
    area(:) = 0.0d0
    volume(:) = 0.0d0

    ! Initialize local arrays
    local_area(:) = 0.0d0
    local_volume(:) = 0.0d0

    ! Compute the area and volume of each glacier.
    ! We need parallel sums, since a glacier can lie on two or more processors.

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0) then
             if (thck(i,j) >= diagnostic_minthck) then
                local_area(ng) = local_area(ng) + cell_area(i,j)
                local_volume(ng) = local_volume(ng) + cell_area(i,j) * thck(i,j)
             endif
          endif
       enddo
    enddo

    area   = parallel_reduce_sum(local_area)
    volume = parallel_reduce_sum(local_volume)

  end subroutine glacier_area_volume

!****************************************************

  subroutine glacier_area_advance_retreat(&
       ewn,           nsn,     &
       nglacier,               &
       cism_glacier_id_init,   &
       cism_glacier_id,        &
       cell_area,              &
       area_initial,           &
       area_current,           &
       area_advance,           &
       area_retreat)

    use cism_parallel, only: parallel_reduce_sum

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

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id_init,        & ! integer glacier ID in the range (1, nglacier), initial value
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier), current value

    real(dp), intent(in) ::  &
         cell_area                      ! grid cell area = dew*dns (m^2); same for all cells

    real(dp), dimension(nglacier), intent(out) ::  &
         area_initial,                & ! initial glacier area
         area_current,                & ! current glacier area
         area_advance,                & ! area of the region where the glacier has advanced (m^2)
         area_retreat                   ! area of the region where the glacier has retreated (m^2)

    ! local variables

    real(dp), dimension(nglacier) :: &
         local_area                     ! area on each processor, before global sum

    integer :: i, j, ng, ngi

    ! Initialize the output arrays
    area_initial(:) = 0.0d0
    area_current(:) = 0.0d0
    area_advance(:) = 0.0d0
    area_retreat(:) = 0.0d0

    ! Compute the area of each glacier over the initial and current masks.
    ! We need parallel sums, since a glacier can lie on two or more processors.

    ! init area
    local_area(:) = 0.0d0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ngi = cism_glacier_id_init(i,j)
          if (ngi > 0) then
             local_area(ngi) = local_area(ngi) + cell_area
          endif
       enddo
    enddo
    area_initial = parallel_reduce_sum(local_area)

    ! current area
    local_area(:) = 0.0d0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0) then
             local_area(ng) = local_area(ng) + cell_area
          endif
       enddo
    enddo
    area_current = parallel_reduce_sum(local_area)

    ! area where the glacier has advanced
    local_area(:) = 0.0d0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ngi = cism_glacier_id_init(i,j)
          ng  = cism_glacier_id(i,j)
          if (ngi == 0 .and. ng > 0) then
             local_area(ng) = local_area(ng) + cell_area
          endif
       enddo
    enddo
    area_advance = parallel_reduce_sum(local_area)

    ! area where the glacier has retreated
    local_area(:) = 0.0d0
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ngi = cism_glacier_id_init(i,j)
          ng  = cism_glacier_id(i,j)
          if (ngi > 0 .and. ng == 0) then
             local_area(ngi) = local_area(ngi) + cell_area
          endif
       enddo
    enddo
    area_retreat = parallel_reduce_sum(local_area)

    ! bug check
    do ng = 1, nglacier
       if (area_initial(ng) + area_advance(ng) - area_retreat(ng) /= area_current(ng)) then
          print*, ' '
          print*, 'WARNING: area mismatch in glacier_area_advance_retreat'
          print*, '   ng, initial, advance, retreat, current:', ng, area_initial(ng)/1.d6, &
               area_advance(ng)/1.d6, area_retreat(ng)/1.d6, area_current(ng)/1.d6
       endif
    enddo

  end subroutine glacier_area_advance_retreat

!****************************************************

  subroutine glacier_accumulation_area_ratio(&
       ewn,           nsn,     &
       nglacier,               &
       cism_glacier_id,        &
       smb,                    &
       aar)

    ! Compute the accumulation area ratio (AAR) for each glacier.
    ! Note: In this subroutine the grid cell area is assumed equal for all cells.

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         smb                            ! surface mass balance (mm/yr w.e.)

    real(dp), dimension(nglacier), intent(out) ::  &
         aar                            ! accumulation area ratio

    ! local variables

    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         ablat_area,                  & ! area of accumulation zone (SMB < 0)
         accum_area                     ! area of accumulation zone (SMB > 0)

    ! initialize
    ablat_area(:) = 0.0d0
    accum_area(:) = 0.0d0

    ! Compute the accumulation and ablation area for each glacier
    ! Note: Grid cells with SMB = 0 are not counted in either zone.

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id(i,j)
          if (ng > 0) then
             if (smb(i,j) > 0.0d0) then
                accum_area(ng) = accum_area(ng) + 1.0d0
             elseif (smb(i,j) < 0.0d0) then
                ablat_area(ng) = ablat_area(ng) + 1.0d0
             endif
          endif
       enddo   ! i
    enddo   ! j

    accum_area = parallel_reduce_sum(accum_area)
    ablat_area = parallel_reduce_sum(ablat_area)

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
