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
       if (associated(glacier%area_init)) deallocate(glacier%area_init)
       if (associated(glacier%volume_init)) deallocate(glacier%volume_init)
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

       ! Allocate glacier arrays with dimension(nglacier).
       ! Note: We should avoid accessing these arrays for grid cells with cism_glacier_id = 0.
       allocate(glacier%glacierid(nglacier))
       allocate(glacier%area(nglacier))
       allocate(glacier%area_init(nglacier))
       allocate(glacier%volume(nglacier))
       allocate(glacier%volume_init(nglacier))
       allocate(glacier%smb(nglacier))
       allocate(glacier%smb_obs(nglacier))
       allocate(glacier%mu_star(nglacier))
       allocate(glacier%alpha_snow(nglacier))
       allocate(glacier%beta_artm(nglacier))

       ! Compute area scale factors
       if (glacier%scale_area) then
          do j = nhalo+1, nsn-nhalo
             do i = nhalo+1, ewn-nhalo
                theta_rad = model%general%lat(i,j) * pi/180.d0
                glacier%area_factor(i,j) = cos(theta_rad)**2
             enddo
          enddo
          call parallel_halo(glacier%area_factor, parallel)
          if (verbose_glacier .and. this_rank == rtest) then
             i = itest; j = jtest
             print*, 'Scale glacier area: i, j, area_factor =', i, j, glacier%area_factor(i,j)
             print*, '   lat, theta, cos(theta) =', model%general%lat(i,j), theta_rad, cos(theta_rad)
          endif
       else
          glacier%area_factor(:,:) = 1.0d0
       endif

       ! Compute the initial area and volume of each glacier.
       ! Only ice thicker than diagnostic_minthck is included in area and volume sums.

       call glacier_area_volume(&
            ewn,           nsn,             &
            nglacier,                       &
            glacier%cism_glacier_id,        &
            dew*dns,                        &
            model%geometry%thck*thk0,       &  ! m
            glacier%diagnostic_minthck,     &  ! m
            glacier%area_factor,            &
            glacier%area,                   &  ! m^2
            glacier%volume)                    ! m^3

       ! Initialize other glacier arrays
       glacier%smb(:)         = 0.0d0
       glacier%area_init(:)   = glacier%area(:)
       glacier%volume_init(:) = glacier%volume(:)
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

       ! If inverting for powerlaw_c, then initialize powerlaw_c to a constant value,
       !  and initialize the inversion target to the initial usrf.
       ! Note: usrf_target_rgi is the thickness at the RGI date, e.g. the
       !        Farinotti et al. consensus thickness).
       !       usrf_target_baseline is the target thickness for the baseline state, which
       !        ideally will evolve to usrf_target_rgi between the baseline date and RGI date.
       ! On restart, powerlaw_c, usrf_target_baseline, and usrf_target_rgi are read from the restart file.

       if (glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          model%basal_physics%powerlaw_c(:,:) = model%basal_physics%powerlaw_c_const
          glacier%usrf_target_baseline(:,:) = model%geometry%usrf(:,:)*thk0
          glacier%usrf_target_rgi(:,:) = model%geometry%usrf(:,:)*thk0
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

    else  ! restart

       ! In this case, most required glacier info has already been read from the restart file.
       ! Here, do some error checks and diagnostics.

       ! From the restart file, nglacier is found as the length of dimension 'glacierid'.
       ! The 1D glacier arrays are then allocated with dimension(nglacier) in subroutine glide_allocarr.
       ! The following glacier arrays should be present in the restart file:
       !     rgi_glacier_id, cism_glacier_id, cism_glacier_id_init, cism_to_rgi_glacier_id,
       !     glacier_mu_star, and powerlaw_c.
       ! If inverting for powerlaw_c, then usrf_target_baseline and usrf_target_rgi are read from the restart file.
       ! If inverting for both mu_star and alpha_snow, then glacier%smb_obs is read from the restart file.

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
          max_glcval = maxval(abs(glacier%smb_rgi))
          max_glcval = parallel_reduce_max(max_glcval)
          if (max_glcval <= 0.0d0) then
             call write_log ('Error, no nonzero values for smb_rgi', GM_FATAL)
          endif
          max_glcval = maxval(glacier%usrf_target_rgi)
          max_glcval = parallel_reduce_max(max_glcval)
          if (max_glcval <= 0.0d0) then
             call write_log ('Error, no positive values for usrf_target_rgi', GM_FATAL)
          endif
       endif

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .and. &
           glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
          max_glcval = maxval(abs(glacier%smb_obs))
          max_glcval = parallel_reduce_max(max_glcval)
          if (max_glcval == 0.d0) then
             call write_log ('Error, no nonzero values for smb_obs', GM_FATAL)
          endif
       else
          ! If a nonzero smb_obs field was read in, then set to zero
          glacier%smb_obs = 0.0d0
       endif

       ! Compute the initial area and volume of each glacier.
       ! This is not necessary for exact restart, but is included as a diagnostic.
       ! Only ice thicker than diagnostic_minthck is included in area and volume sums.

       call glacier_area_volume(&
            ewn,           nsn,             &
            nglacier,                       &
            glacier%cism_glacier_id,        &
            dew*dns,                        &
            model%geometry%thck*thk0,       &  ! m
            glacier%diagnostic_minthck,     &  ! m
            glacier%area_factor,            &
            glacier%area,                   &  ! m^2
            glacier%volume)                    ! m^3

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
         tsrf,                    & ! local array for surface air temperature (deg C)
         Tpos,                    & ! max(artm - tmlt, 0.0)
         snow,                    & ! snowfall rate (mm w.e./yr)
         Tpos_aux,                & ! max(artm - tmlt, 0.0), auxiliary field
         snow_aux,                & ! snowfall rate (mm w.e./yr), auxiliary field
         artm_rgi,                & ! artm, RGI date
         precip_rgi,              & ! precip rate, RGI date
         Tpos_rgi,                & ! max(artm - tmlt, 0.0), RGI date
         snow_rgi,                & ! snowfall rate, RGI date
         mu_star_2d,              & ! 2D version of glacier%mu_star
         alpha_snow_2d,           & ! 2D version of glacier%alpha_snow
         smb_annmean_init,        & ! annual mean SMB for each glacier cell over init area (mm/yr w.e.)
         smb_annmean,             & ! annual mean SMB for each glacier cell over current area (mm/yr w.e.)
         delta_smb_rgi,           & ! SMB anomaly between the baseline date and the RGI date (mm/yr w.e.)
         delta_smb_aux              ! SMB anomaly between the baseline date and the auxiliary date (mm/yr w.e.)

    real(dp), dimension(model%general%ewn-1, model%general%nsn-1) ::  &
         stag_thck,                   & ! ice thickness at vertices (m)
         stag_thck_target,            & ! target ice thickness at vertices (m)
         stag_dthck_dt                  ! rate of change of ice thickness at vertices (m/yr)

    type(parallel_type) :: parallel ! info for parallel communication

    real(dp), save ::  &              ! time since the last averaging computation (yr);
         time_since_last_avg = 0.0d0  ! compute the average once a year

    real(dp), dimension(glacier%nglacier) :: &
         area_old,                & ! glacier%area from the previous inversion step
         darea_dt,                & ! rate of change of glacier area over the inversion interval
         smb_init_area,           & ! SMB over initial area determined by cism_glacier_id_init
         smb_new_area,            & ! SMB over new area determined by cism_glacier_id
         aar_init,                & ! accumulation area ratio over the initial area using cism_glacier_id_init
         aar                        ! accumulation area ratio over the new area using cism_glacier_id

    ! Note: The glacier type includes the following:
    ! integer ::  nglacier          ! number of glaciers in the global domain
    ! integer ::  ngdiag            ! CISM index of diagnostic glacier
    ! real(dp), dimension(:) :: area              ! glacier area (m^2)
    ! real(dp), dimension(:) :: volume            ! glacier volume (m^3)
    ! real(dp), dimension(:) :: area_init         ! initial glacier area (m^2)
    ! real(dp), dimension(:) :: volume_init       ! initial glacier volume (m^3)
    ! real(dp), dimension(:) :: mu_star           ! SMB parameter for each glacier (mm/yr w.e./deg K)
    ! real(dp), dimension(:) :: alpha_snow        ! snow factor for each glacier (unitless)
    ! real(dp), dimension(:) :: beta_artm         ! artm correction for each glacier (deg C)
    ! real(dp), dimension(:) :: smb_obs           ! observed SMB for each glacier (mm/yr w.e.)
    ! integer, dimension(:,:) :: cism_glacier_id       ! CISM glacier ID for each grid cell
    ! integer, dimension(:,:) :: cism_glacier_id_init  ! initial value of CISM glacier ID
    ! integer, dimension(:,:) :: smb_glacier_id        ! CISM glacier ID that determines where SMB is applied
    ! integer, dimension(:,:) :: smb_glacier_id_init   ! like smb_glacier_id, but based on cism_glacier_id_init
    ! real(dp), dimension(:,:) :: snow_2d              ! snow accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: Tpos_2d              ! max(artm - tmlt,0) accumulated and averaged over 1 year
    ! real(dp), dimension(:,:) :: snow_aux_2d          ! snow accumulated and averaged over 1 year, auxiliary field
    ! real(dp), dimension(:,:) :: Tpos_aux_2d          ! max(artm - tmlt,0) accumulated and averaged over 1 year, auxiliary field
    ! real(dp), dimension(:,:) :: snow_rgi_2d          ! snow accumulated and averaged over 1 year, RGI date
    ! real(dp), dimension(:,:) :: Tpos_rgi_2d          ! max(artm - tmlt,0) accumulated and averaged over 1 year, RGI date
    ! real(dp), dimension(:,:) :: dthck_dt_2d          ! dthck_dt accumulated and averaged over 1 year

    ! SMB and accumulation area diagnostics
    real(dp), dimension(:), allocatable :: &
         area_acc_init, area_abl_init, f_accum_init, &
         area_acc_new,  area_abl_new,  f_accum_new

    ! Note: The following areas are computed based on the cism_glacier_id masks, without a min thickness criterion
    real(dp), dimension(glacier%nglacier) ::  &
         area_initial, area_current,    &  ! initial and current glacier areas (m^2)
         area_advance, area_retreat        ! areas of glacier advance and retreat relative to initial mask (m^2)

    real(dp) :: area_sum
    real(dp) :: usrf_aux    ! estimated surface elevation in auxiliary climate
    real(dp) :: usrf_rgi    ! estimated surface elevation in RGI climate
    real(dp), parameter :: diagnostic_volume_threshold = 1.0d9  ! volume threshold for big glaciers (m^3)

    !TODO - Make these config parameters
    real(dp), parameter :: &
         baseline_date = 1984.d0,   & ! date of baseline climate, when glaciers are assumed to be in balance
         rgi_date = 2003.d0,        & ! RGI reference date, when we have observed glacier outlines and thickness targets
         smbobs_date = 2010.d0        ! date of recent climate data, when we have smb_obs for glaciers out of balance

    real(dp) :: rgi_date_frac

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
       !TODO - 'if' logic around the aux and rgi fields

       glacier%snow_2d = 0.0d0
       glacier%Tpos_2d = 0.0d0
       glacier%snow_aux_2d = 0.0d0
       glacier%Tpos_aux_2d = 0.0d0
       glacier%snow_rgi_2d = 0.0d0
       glacier%Tpos_rgi_2d = 0.0d0
       glacier%dthck_dt_2d = 0.0d0

       ! Compute the SMB anomaly for the RGI and auxiliary climates relative to the baseline climate.
       ! This is done once a year; smb, smb_rgi, and smb_aux are updated at the end of the previous year.

       where (glacier%smb_glacier_id_init > 0 .and. model%climate%smb /= 0.0d0 .and. glacier%smb_rgi /= 0.0d0)
          delta_smb_rgi = glacier%smb_rgi - model%climate%smb
       elsewhere
          delta_smb_rgi = 0.0d0
       endwhere
       glacier%delta_usrf_rgi(:,:) = &
            delta_smb_rgi(:,:)*(rhow/rhoi)/1000.d0 * (rgi_date - baseline_date)/2.d0

       where (glacier%smb_glacier_id_init > 0 .and. model%climate%smb /= 0.0d0 &
            .and. model%climate%smb_aux /= 0.0d0)
          delta_smb_aux = model%climate%smb_aux - model%climate%smb
       elsewhere
          delta_smb_aux = 0.0d0
       endwhere
       glacier%delta_usrf_aux(:,:) = &
            delta_smb_aux(:,:)*(rhow/rhoi)/1000.d0 * (smbobs_date - baseline_date)/2.0d0  ! m ice

       ! Adjust the baseline target. The baseline target should exceed the RGI target by abs(delta_usrf_rgi),
       !  assuming the ice thins between the baseline and RGI dates.
       ! Then, provided usrf is close to usrf_target_baseline in the spin-up, usrf will be close to
       !  usrf_target_rgi when a forward run starting from the baseline date reaches the RGI date.

       glacier%usrf_target_baseline(:,:) = &
            glacier%usrf_target_rgi(:,:) - glacier%delta_usrf_rgi(:,:)

       ! Make sure the target is not below the topography
       glacier%usrf_target_baseline = &
            max(glacier%usrf_target_baseline, (model%geometry%topg + model%climate%eus)*thk0)

       if (verbose_glacier .and. this_rank == rtest) then
          i = itest; j = jtest
          print*, ' '
          print*, 'RGI usrf correction, delta_smb:', &
               glacier%delta_usrf_rgi(i,j), delta_smb_rgi(i,j)
          print*,    'usrf_target_rgi, new usrf_target_baseline =', &
               glacier%usrf_target_rgi(i,j), glacier%usrf_target_baseline(i,j)
          print*, 'Aux usrf correction, delta_smb:', &
               glacier%delta_usrf_aux(i,j), delta_smb_aux(i,j)
       endif

    endif   ! time_since_last_avg = 0

    ! Halo updates for snow and artm
    ! Note: artm_corrected is the input artm, possibly corrected to include an anomaly term.
    ! Note: snow_calc is the snow calculation option: Either use the snowfall rate directly,
    !       or compute the snowfall rate from the precip rate and downscaled artm.
    !TODO - Not sure these are needed. Maybe can save halo updates for the annual-averaged snow and Tpos

    if (glacier%snow_calc == GLACIER_SNOW_CALC_SNOW) then
       call parallel_halo(model%climate%snow, parallel)
    elseif (glacier%snow_calc == GLACIER_SNOW_CALC_PRECIP_ARTM) then
       call parallel_halo(model%climate%precip, parallel)
    endif
    call parallel_halo(model%climate%artm_corrected, parallel)

    ! Compute artm and Tpos for the baseline climate at the current surface elevation, usrf

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = glacier%smb_glacier_id_init(i,j)
          if (ng > 0) then
             model%climate%artm(i,j) = model%climate%artm_ref(i,j)  &
                  - (model%geometry%usrf(i,j)*thk0 - model%climate%usrf_ref(i,j))*model%climate%t_lapse &
                  + glacier%beta_artm(ng)
          else
             model%climate%artm(i,j) = model%climate%artm_ref(i,j)  &
                  - (model%geometry%usrf(i,j)*thk0 - model%climate%usrf_ref(i,j))*model%climate%t_lapse
          endif
          Tpos(i,j) = max(model%climate%artm(i,j) - glacier%tmlt, 0.0d0)
       enddo
    enddo

    ! Compute artm and Tpos for the auxiliary climate at the extrapolated surface elevation, usrf_aux.
    ! We estimate usrf_aux = usrf + (dSMB/2)*dt,
    !    where dSMB = smb_aux - smb is the difference in SMB between the baseline and auxiliary climate,
    !          (so dSMB/2 is the average SMB anomaly over that period), and dt is the number of years elapsed.
    ! In other words, assume that the entire SMB anomaly is used to melt ice, without the
    !  flow having time to adjust.

    ! Note: The fields with the 'aux' suffix are used only for inversion
    !       and are needed only for cells that are initially glacier-covered.
    !       If inversion is turned off, these fields will equal 0.
    !       TODO: Add 'if inversion' logic so that only Tpos and snow are always computed?

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          usrf_aux = model%geometry%usrf(i,j)*thk0 + glacier%delta_usrf_aux(i,j)
          ng = glacier%smb_glacier_id_init(i,j)
          if (ng > 0) then
             model%climate%artm_aux(i,j) = model%climate%artm_ref_aux(i,j)  &
                  - (usrf_aux - model%climate%usrf_ref(i,j))*model%climate%t_lapse  &
                  + glacier%beta_artm(ng)
          else
             model%climate%artm_aux(i,j) = model%climate%artm_ref_aux(i,j)  &
                  - (usrf_aux - model%climate%usrf_ref(i,j))*model%climate%t_lapse
          endif
          Tpos_aux(i,j) = max(model%climate%artm_aux(i,j) - glacier%tmlt, 0.0d0)
       enddo
    enddo

    ! Estimate artm, Tpos, and snow or precip for the RGI climate by interpolation.

    rgi_date_frac = (rgi_date - baseline_date) / (smbobs_date - baseline_date)

    artm_rgi(:,:) = &
         (1.d0 - rgi_date_frac) * model%climate%artm(:,:)  &
              +  rgi_date_frac  * model%climate%artm_aux(:,:)

    Tpos_rgi(:,:) = max(artm_rgi(:,:) - glacier%tmlt, 0.0d0)

    if (glacier%snow_calc == GLACIER_SNOW_CALC_SNOW) then
    elseif (glacier%snow_calc == GLACIER_SNOW_CALC_PRECIP_ARTM) then
    endif

    ! Compute the snowfall rate for each climate.
    ! Note: Depending on glacier%snow_calc, we either use the snowfall rate directly,
    !       or compute snowfall based on the input precip and artm

    if (glacier%snow_calc == GLACIER_SNOW_CALC_SNOW) then

       snow(:,:) = model%climate%snow(:,:)
       snow_aux(:,:) = model%climate%snow_aux(:,:)

       snow_rgi(:,:) = &
            (1.d0 - rgi_date_frac) * snow(:,:)  &
                 +  rgi_date_frac  * snow_aux(:,:)

    elseif (glacier%snow_calc == GLACIER_SNOW_CALC_PRECIP_ARTM) then

       call glacier_calc_snow(&
            ewn,       nsn,                   &
            glacier%snow_threshold_min,       &
            glacier%snow_threshold_max,       &
            model%climate%precip,             &
            model%climate%artm,               &
            snow)

       call glacier_calc_snow(&
            ewn,       nsn,                   &
            glacier%snow_threshold_min,       &
            glacier%snow_threshold_max,       &
            model%climate%precip_aux,         &
            model%climate%artm_aux,           &
            snow_aux)

       precip_rgi(:,:) = &
            (1.d0 - rgi_date_frac) * model%climate%precip(:,:)  &
                 +  rgi_date_frac  * model%climate%precip_aux(:,:)

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
       print*, ' '
       print*, 'glacier lapse-rate correction, diag cell (r, i, j) =', rtest, i, j
       print*, '   usrf_ref, usrf, diff:', &
            model%climate%usrf_ref(i,j), model%geometry%usrf(i,j)*thk0, &
            model%geometry%usrf(i,j)*thk0 - model%climate%usrf_ref(i,j)
       print*, 'Baseline artm_ref, artm, Tpos, snow, smb:', &
            model%climate%artm_ref(i,j), model%climate%artm(i,j), &
            Tpos(i,j), snow(i,j), model%climate%smb(i,j)
       print*, 'RGI artm, Tpos, snow:', &
            artm_rgi(i,j), Tpos_rgi(i,j), snow_rgi(i,j)
       print*, 'Aux artm, Tpos, snow:', &
            model%climate%artm_aux(i,j), Tpos_aux(i,j), snow_aux(i,j)
       print*, ' '
    endif   ! verbose

    ! Accumulate snow_2d, Tpos_2d, and dthck_dt_2d over this timestep

    time_since_last_avg = time_since_last_avg + dt

    glacier%snow_2d = glacier%snow_2d + snow * dt
    glacier%Tpos_2d = glacier%Tpos_2d + Tpos * dt
    glacier%snow_rgi_2d = glacier%snow_rgi_2d + snow_rgi * dt
    glacier%Tpos_rgi_2d = glacier%Tpos_rgi_2d + Tpos_rgi * dt
    glacier%snow_aux_2d = glacier%snow_aux_2d + snow_aux * dt
    glacier%Tpos_aux_2d = glacier%Tpos_aux_2d + Tpos_aux * dt
    glacier%dthck_dt_2d = glacier%dthck_dt_2d + dthck_dt * dt

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_glacier_inversion, diag cell (r, i, j) =', rtest, itest, jtest
       i = itest; j = jtest
       print*, '   r, i, j, time, artm, snow, Tpos:', &
            this_rank, i, j, model%numerics%time, &
            model%climate%artm_corrected(i,j), snow(i,j), Tpos(i,j)
       print*, '   r, i, j, time, artm_aux, snow_aux, Tpos_aux:', &
            this_rank, i, j, model%numerics%time, &
            model%climate%artm_aux(i,j), snow_aux(i,j), Tpos_aux(i,j)
    endif

    ! Check whether it is time to do the inversion and update other glacier fields.
    ! Note: time_since_last_avg is real(dp) with units of yr;
    !       glacier_update_interval is an integer number of years.

    if (abs(time_since_last_avg - real(glacier_update_interval,dp)) < eps08) then

       if (verbose_glacier .and. this_rank == rtest) then
          print*, 'calculate_glacier_2d_to_1ds, time_since_last_avg =', time_since_last_avg
       endif

       ! Compute the average of glacier fields over the accumulation period

       glacier%snow_2d = glacier%snow_2d / time_since_last_avg
       glacier%Tpos_2d = glacier%Tpos_2d / time_since_last_avg
       glacier%snow_rgi_2d = glacier%snow_rgi_2d / time_since_last_avg
       glacier%Tpos_rgi_2d = glacier%Tpos_rgi_2d / time_since_last_avg
       glacier%snow_aux_2d = glacier%snow_aux_2d / time_since_last_avg
       glacier%Tpos_aux_2d = glacier%Tpos_aux_2d / time_since_last_avg
       glacier%dthck_dt_2d = glacier%dthck_dt_2d / time_since_last_avg

       time_since_last_avg = 0.0d0

       if (verbose_glacier .and. this_rank == rtest) then
          i = itest; j = jtest
          print*, ' '
          print*, 'Annual averages, r, i, j:', rtest, itest, jtest
          print*, '   snow (mm/yr)       =', glacier%snow_2d(i,j)
          print*, '   Tpos (deg C)       =', glacier%Tpos_2d(i,j)
          print*, '   snow_rgi (mm/yr)   =', glacier%snow_rgi_2d(i,j)
          print*, '   Tpos_rgi (deg C)   =', glacier%Tpos_rgi_2d(i,j)
          print*, '   snow_aux (mm/yr)   =', glacier%snow_aux_2d(i,j)
          print*, '   Tpos_aux (deg C)   =', glacier%Tpos_aux_2d(i,j)
          print*, '   dthck_dt (m/yr)    =', glacier%dthck_dt_2d(i,j)
       endif

       ! Invert for mu_star
       ! This can be done in either of two ways:
       ! (1) set_mu_star = 1, set_alpha_snow = 0 (1-parameter inversion)
       !     In this case, mu_star is chosen such that SMB ~ 0 over the initial glacier footprint, given
       !     the input temperature and snow/precip fields (without the 'aux' suffix).
       ! (2) set_mu_star = 1, set_alpha_snow = 1 (2-parameter inversion)
       !     In this case, mu_star and alpha_snow are chosen jointly such that
       !     (a) SMB = 0 over the initial footprint given the baseline temperature and snow/precip, and
       !     (b) SMB = smb_obs given the auxiliary temperature and snow/precip.
       ! The code aborts at startup if set to invert for alpha_snow without inverting for mu_star.

       if (glacier%set_mu_star == GLACIER_MU_STAR_INVERSION) then

          if (glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then

             ! invert for both mu_star and alpha_snow, based on two SMB conditions
             ! (SMB = 0 in a balanced climate, SMB = smb_obs in an out-of-balance climate)
             ! Note: glacier%smb_obs, glacier%mu_star, and glacier%alpha_snow are 1D glacier-specific fields.

             call glacier_invert_mu_star_alpha_snow(&
                  ewn,                    nsn,                   &
                  itest,     jtest,       rtest,                 &
                  nglacier,               ngdiag,                &
                  glacier%smb_glacier_id_init,                   &
                  glacier%smb_obs,                               &
                  glacier%cism_to_rgi_glacier_id,                &  ! diagnostic only
                  glacier%area_init,      glacier%volume_init,   &  ! diagnostic only
                  glacier%snow_2d,        glacier%Tpos_2d,       &
                  glacier%snow_aux_2d,    glacier%Tpos_aux_2d,   &
                  glacier%mu_star_const,                         &
                  glacier%mu_star_min,    glacier%mu_star_max,   &
                  glacier%alpha_snow_const,                      &
                  glacier%alpha_snow_min, glacier%alpha_snow_max,&
                  glacier%beta_artm_max,                         &
                  glacier%beta_artm_increment,                   &
                  glacier%mu_star,        glacier%alpha_snow,    &
                  glacier%beta_artm)

          else  ! not inverting for alpha_snow

             ! invert for mu_star based on a single SMB condition (balanced climate)
             ! Choose mu_star for each glacier to match smb = 0 over the initial glacier footprint.
             ! Use the default value of alpha_snow (typically = 1.0).

             call glacier_invert_mu_star(&
                  ewn,                  nsn,                 &
                  itest,     jtest,     rtest,               &
                  nglacier,             ngdiag,              &
                  glacier%smb_glacier_id_init,               &
                  glacier%smb_obs,                           &
                  glacier%snow_2d,      glacier%Tpos_2d,     &
                  glacier%mu_star_min,  glacier%mu_star_max, &
                  glacier%mu_star)

          endif  ! set_alpha_snow

       endif   ! invert for mu_star

       !TODO - A lot of optional diagnostic output follows.
       !       Need to consolidate and move some of it to subroutines.

       ! Given mu_star and alpha_snow, compute the average SMB for each glacier,
       !  based on its initial area and its current area (for diagnostic purposes only).

       ! Convert mu_star and alpha_snow to 2D fields, scattering over the initial glacier area

       call glacier_1d_to_2d(&
            ewn,              nsn,                         &
            nglacier,         glacier%smb_glacier_id_init, &
            glacier%mu_star,  mu_star_2d)

       call glacier_1d_to_2d(&
            ewn,                 nsn,                          &
            nglacier,            glacier%smb_glacier_id_init,  &
            glacier%alpha_snow,  alpha_snow_2d)

       ! Compute the SMB for each grid cell over the initial glacier area

       where (glacier%smb_glacier_id_init > 0)
          smb_annmean_init = alpha_snow_2d * glacier%snow_2d - mu_star_2d * glacier%Tpos_2d
       elsewhere
          smb_annmean_init = 0.0d0
       endwhere

       ! Compute the average SMB for each glacier over the initial glacier area
       ! TODO - Rename smb_init_area?

       call glacier_2d_to_1d(&
            ewn,              nsn,                         &
            nglacier,         glacier%smb_glacier_id_init, &
            smb_annmean_init, smb_init_area)

       ! Repeat for the current glacier area

       ! Convert mu_star and alpha_snow to 2D fields, scattering over the current glacier area

       call glacier_1d_to_2d(&
            ewn,              nsn,                         &
            nglacier,         glacier%smb_glacier_id,      &
            glacier%mu_star,  mu_star_2d)

       call glacier_1d_to_2d(&
            ewn,                 nsn,                      &
            nglacier,            glacier%smb_glacier_id,   &
            glacier%alpha_snow,  alpha_snow_2d)

       ! Compute the SMB for each grid cell based on the current glacier area

       where (glacier%smb_glacier_id > 0)
          smb_annmean = alpha_snow_2d * glacier%snow_2d - mu_star_2d * glacier%Tpos_2d
       elsewhere
          smb_annmean = 0.0d0
       endwhere

       call parallel_halo(smb_annmean, parallel)

       ! Compute the average SMB for each glacier over the current glacier area

       call glacier_2d_to_1d(&
            ewn,          nsn,                     &
            nglacier,     glacier%smb_glacier_id, &
            smb_annmean,  smb_new_area)

       ! some local diagnostics

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'cism_glacier_id_init:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') glacier%cism_glacier_id_init(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'cism_glacier_id:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') glacier%cism_glacier_id(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'thck:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'smb_annmean (based on initial smb_glacier_id):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') smb_annmean_init(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'smb_annmean (based on current smb_glacier_id):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') smb_annmean(i,j)
             enddo
             print*, ' '
          enddo
       endif   ! verbose

       ! accumulation and ablation area diagnostics
       !TODO - Remove since another subroutine does this?

       allocate(area_acc_init(nglacier))
       allocate(area_abl_init(nglacier))
       allocate(f_accum_init(nglacier))
       allocate(area_acc_new(nglacier))
       allocate(area_abl_new(nglacier))
       allocate(f_accum_new(nglacier))

       area_acc_init = 0.0d0
       area_abl_init = 0.0d0
       f_accum_init = 0.0d0
       area_acc_new = 0.0d0
       area_abl_new = 0.0d0
       f_accum_new = 0.0d0

       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             ! initial glacier ID
             ng = glacier%smb_glacier_id_init(i,j)
             if (ng > 0) then
                if (smb_annmean_init(i,j) >= 0.0d0) then
                   area_acc_init(ng) = area_acc_init(ng) + dew*dns
                else
                   area_abl_init(ng) = area_abl_init(ng) + dew*dns
                endif
             endif
             ! current glacier ID
             ng = glacier%smb_glacier_id(i,j)
             if (ng > 0) then
                if (smb_annmean(i,j) >= 0.0d0) then
                   area_acc_new(ng) = area_acc_new(ng) + dew*dns
                else
                   area_abl_new(ng) = area_abl_new(ng) + dew*dns
                endif
             endif
          enddo   ! i
       enddo   ! j

       area_acc_init = parallel_reduce_sum(area_acc_init)
       area_abl_init = parallel_reduce_sum(area_abl_init)
       area_acc_new = parallel_reduce_sum(area_acc_new)
       area_abl_new = parallel_reduce_sum(area_abl_new)

       do ng = 1, nglacier
          area_sum = area_acc_init(ng) + area_abl_init(ng)
          if (area_sum > 0.0d0) then
             f_accum_init(ng) = area_acc_init(ng) / area_sum
          endif
          area_sum = area_acc_new(ng) + area_abl_new(ng)
          if (area_sum > 0.0d0) then
             f_accum_new(ng) = area_acc_new(ng) / area_sum
          endif
       enddo

       ! advance/retreat diagnostics
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
          ng = ngdiag
          if (ng > 0) then
             print*, 'ngdiag, smb_init_area (mm/yr w.e.), smb_new_area, mu_star, alpha_snow, beta_artm, beta_aux:'
             write(6,'(i6,5f12.4)') ng, smb_init_area(ng), smb_new_area(ng), glacier%mu_star(ng), &
                  glacier%alpha_snow(ng), glacier%beta_artm(ng)
          endif
          print*, ' '
          print*, 'Selected big glaciers:'
          print*, 'ng,      Ainit,      A,      Vinit,      V,  smb_iniA, smb_newA, mu_star,  alpha_snow, beta_artm, smb_obs'
          do ng = 1, nglacier
             if (glacier%volume_init(ng) > diagnostic_volume_threshold .or. ng == ngdiag) then  ! big glacier
                write(6,'(i6,10f10.3)') ng, glacier%area_init(ng)/1.e6, glacier%area(ng)/1.e6, &
                     glacier%volume_init(ng)/1.0d9, glacier%volume(ng)/1.0d9, &
                     smb_init_area(ng), smb_new_area(ng), glacier%mu_star(ng), glacier%alpha_snow(ng), &
                     glacier%beta_artm(ng), glacier%smb_obs(ng)
             endif
          enddo
       endif

!!       if (verbose_glacier .and. this_rank == rtest) then
       if (verbose_glacier .and. 0 == 1) then
          print*, ' '
          print*, 'Accumulation/ablation diagnostics:'
          print*, 'ng,    A_acc_tgt, A_abl_tgt, f_acc_tgt, A_acc_new, A_abl_new, f_acc_new'
          do ng = 1, nglacier
             if (glacier%volume_init(ng) > 1.0d9 .or. ng == ngdiag) then  ! big glacier, > 1 km^3
                write(6,'(i6,6f10.3)') ng, area_acc_init(ng)/1.e6, area_abl_init(ng)/1.e6, f_accum_init(ng), &
                     area_acc_new(ng)/1.e6, area_abl_new(ng)/1.e6, f_accum_new(ng)
             endif
          enddo
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
               ewn,                 nsn,           &
               glacier%dthck_dt_2d, stag_dthck_dt)

          if (verbose_glacier .and. this_rank == rtest) then
             print*, ' '
             print*, 'call glacier_invert_powerlaw_c, time (yr) =', model%numerics%time
          endif

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

       endif   ! powerlaw_c_inversion

       !-------------------------------------------------------------------------
       ! Update glacier IDs based on advance and retreat since the last update.
       !-------------------------------------------------------------------------

       ! Assign nonzero IDs in grid cells where ice has reached the minimum glacier thickness.
       ! Remove IDs in grid cells where ice is now thinnier than the minimum thickness.
       ! Adjust IDs to prevent spurious advance due to SMB differences in adjacent glaciers.

       !TODO - Check the logic again.
       call glacier_advance_retreat(&
            ewn,             nsn,           &
            itest,   jtest,  rtest,         &
            nglacier,                       &
            glacier%minthck,                &  ! m
            thck,                           &  ! m
            smb_annmean,                    &  ! mm/yr w.e.
            glacier%snow_2d,                &  ! mm/yr w.e.
            glacier%Tpos_2d,                &  ! deg C
            glacier%mu_star,                &  ! mm/yr/deg
            glacier%alpha_snow,             &  ! unitless
            glacier%cism_glacier_id_init,   &
            glacier%cism_glacier_id,        &
            parallel)

       ! Remove snowfields, defined as isolated cells (or patches of cells) located outside
       ! the initial glacier footprint, and disconnected from the initial glacier.

       !TODO - Debug; try to avoid snowfields late in the simulation
       call remove_snowfields(&
            ewn,          nsn,              &
            parallel,                       &
            itest, jtest, rtest,            &
            thck,                           &
            glacier%cism_glacier_id_init,   &
            glacier%cism_glacier_id)

       ! Update the masks of cells where SMB can be nonzero, based on
       !  (1) initial glacier IDs, and (2) current glacier IDs.
       ! The smb_glacier_id_init mask is used for inversion.
       ! The smb_glacier_id mask determines where the SMB is applied during the next timestep.

       ! Compute smb_glacier_id as the union of
       !       (1) cgii > 0 and cgi > 0
       !       (2) cgii > 0, cgi = 0, and SMB > 0
       !       (3) cgii = 0, cgi > 0, and SMB < 0
       !       Given snow_2d, Tpos_2d, alpha, and mu, we can compute a potential SMB for each cell.
       !           Let SMB = alpha_snow * snow - mu_star * tpos, using ng corresponding to cgi, cgii, or both
       !           where alpha_snow and mu_star are per glacier, and snow and tpos are annual averages
       !       Use the potential SMB to determine smb_glacier_id in advanced and retreated cells.

       call update_smb_glacier_id(&
            ewn,           nsn,             &
            itest, jtest,  rtest,           &
            glacier%nglacier,               &
            glacier%snow_2d,                &  ! mm/yr w.e.
            glacier%Tpos_2d,                &  ! deg C
            glacier%mu_star,                &  ! mm/yr/deg
            glacier%alpha_snow,             &  ! unitless
            glacier%cism_glacier_id_init,   &
            glacier%cism_glacier_id,        &
            glacier%smb_glacier_id_init,    &
            glacier%smb_glacier_id,         &
            parallel)

       ! Using the new smb_glacier_id mask, compute model%climate%smb for the next year.
       !       Cells with smb_glacier_id = 0 have smb = 0.

       ! Use an empirical relationship based on Maussion et al. (2019):
       !
       !     SMB = alpha_snow * snow - mu_star * max(artm - tmlt, 0),
       !
       ! where snow = monthly mean snowfall rate (mm/yr w.e.),
       !       alpha_snow is a glacier-specific tuning parameter (a scalar of order 1)
       !       mu_star is a glacier-specific tuning parameter (mm/yr w.e./deg C),
       !       atrm = monthly mean air temperature (deg C),
       !       tmlt = monthly mean air temp above which ablation occurs (deg C)

       do j = 1, nsn
          do i = 1, ewn
             ng = glacier%smb_glacier_id(i,j)
             if (ng > 0) then
                model%climate%smb(i,j) = &
                     glacier%alpha_snow(ng)*glacier%snow_2d(i,j) - glacier%mu_star(ng)*glacier%Tpos_2d(i,j)
             else
                model%climate%smb(i,j) = 0.0d0
             endif
          enddo
       enddo

       do j = 1, nsn
          do i = 1, ewn
             ng = glacier%smb_glacier_id(i,j)
             if (ng > 0) then
                glacier%smb_rgi(i,j) = &
                     glacier%alpha_snow(ng)*glacier%snow_rgi_2d(i,j) &
                     - glacier%mu_star(ng)*glacier%Tpos_rgi_2d(i,j)
             else
                glacier%smb_rgi(i,j) = 0.0d0
             endif
          enddo
       enddo

       do j = 1, nsn
          do i = 1, ewn
             ng = glacier%smb_glacier_id(i,j)
             if (ng > 0) then
                model%climate%smb_aux(i,j) = &
                     glacier%alpha_snow(ng)*glacier%snow_aux_2d(i,j) &
                     - glacier%mu_star(ng)*glacier%Tpos_aux_2d(i,j)
             else
                model%climate%smb_aux(i,j) = 0.0d0
             endif
          enddo
       enddo

       call parallel_halo(model%climate%smb, parallel)
       call parallel_halo(model%climate%smb_aux, parallel)
       call parallel_halo(glacier%smb_rgi, parallel)

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'New smb_glacier_id_init:'
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
          print*, 'model%climate%smb:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f11.3)',advance='no') model%climate%smb(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'smb_rgi:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f11.3)',advance='no') glacier%smb_rgi(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'smb_aux:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f11.3)',advance='no') model%climate%smb_aux(i,j)
             enddo
             print*, ' '
          enddo
       endif

       ! Update the glacier area and volume (diagnostic only)

       call glacier_area_volume(&
            ewn,           nsn,              &
            nglacier,                        &
            glacier%cism_glacier_id,         &
            dew*dns,                         &  ! m^2
            thck,                            &  ! m
            glacier%diagnostic_minthck,      &  ! m
            glacier%area_factor,             &
            glacier%area,                    &  ! m^2
            glacier%volume)                     ! m^3

       if (verbose_glacier .and. this_rank == rtest) then
          print*, ' '
          print*, 'Update area (km^2) and volume (km^3) for glacier:', ngdiag
          print*, '   Init area and volume:', &
               glacier%area_init(ngdiag)/1.0d6, glacier%volume_init(ngdiag)/1.0d9
          print*, 'Current area and volume:', &
                  glacier%area(ngdiag)/1.0d6, glacier%volume(ngdiag)/1.0d9
          print*, ' '
       endif

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
       glacier_smb_obs,                 &
       snow_2d,          Tpos_2d,       &
       mu_star_min,      mu_star_max,   &
       mu_star)

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
         glacier_smb_obs                ! observed glacier-average SMB (mm/yr w.e.)

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow_2d,                     & ! time-avg snowfall for each cell (mm/yr w.e.)
         Tpos_2d                        ! time-avg of max(artm - tmlt, 0) for each cell (deg)

    real(dp), intent(in) :: &
         mu_star_min, mu_star_max       ! min and max allowed values of mu_star

    real(dp), dimension(nglacier), intent(inout) :: &
         mu_star                        ! glacier-specific SMB tuning parameter (mm/yr w.e./deg)

    ! local variables
    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         glacier_snow, glacier_Tpos     ! glacier-average snowfall and Tpos

    character(len=100) :: message

    ! Compute mu_star for each glacier such that SMB = smb_obs over the initial extent.
    ! Here, the initial extent includes an ablation zone of glacier-free cells adjacent
    !  to glacier-covered cells.
    !
    ! The SMB for glacier ng is given by
    !      sum_ij(smb) = sum_ij(snow) - mu_star(ng) * sum_ij(Tpos),
    ! where Tpos = max(artm - tmlt, 0),
    ! and sum_ij notes a sum over all cells (i,j) in the glacier.
    !
    ! Rearranging, we get
    !      mu_star(ng) = (sum_ij(snow) - sum_ij(smb) / sum_ij(Tpos)
    !
    ! Thus, given the annual average of snow and Tpos for each grid cell in a glacier,
    ! we can find mu_star such that SMB = smb_obs.
    !
    ! Notes:
    !
    ! (1) This approach works only for land-based glaciers.
    !     TODO: Modify for marine-terminating glaciers.
    ! (2) Assuming climatological forcing with smb_obs prescribed, mu_star has nearly the same value
    !     throughout the inversion.  It changes slightly as surface elevation changes, modifying Tpos.

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glacier_invert_mu_star'
    endif

    ! Compute average snowfall, Tpos, and SMB over the initial extent of each glacier

    call glacier_2d_to_1d(&
         ewn,           nsn,                   &
         nglacier,      smb_glacier_id_init,   &
         snow_2d,       glacier_snow)

    call glacier_2d_to_1d(&
         ewn,           nsn,                   &
         nglacier,      smb_glacier_id_init,   &
         Tpos_2d,       glacier_Tpos)

    ! For each glacier, compute the new mu_star

    do ng = 1, nglacier

       if (glacier_Tpos(ng) > 0.0d0) then  ! ablation is nonzero

          ! Compute the value of mu_star that will give the desired SMB over the target area
          mu_star(ng) = (glacier_snow(ng) - glacier_smb_obs(ng)) / glacier_Tpos(ng)

          ! Limit to a physically reasonable range
          mu_star(ng) = min(mu_star(ng), mu_star_max)
          mu_star(ng) = max(mu_star(ng), mu_star_min)

          if (verbose_glacier .and. this_rank == rtest .and. ng == ngdiag) then
             print*, ' '
             print*, 'ng, glacier-average snow, Tpos, smb_obs:', &
                  ng, glacier_snow(ng), glacier_Tpos(ng), glacier_smb_obs(ng)
             print*, 'New mu_star:', mu_star(ng)
          endif

       else   ! glacier_Tpos = 0; no ablation

          mu_star(ng) = mu_star_max

          if (verbose_glacier .and. this_rank == rtest) then
             print*, 'Warning: no ablation for glacier', ng
          endif

       endif

    enddo   ! ng

  end subroutine glacier_invert_mu_star

!****************************************************

  subroutine glacier_invert_mu_star_alpha_snow(&
       ewn,              nsn,            &
       itest,   jtest,   rtest,          &
       nglacier,         ngdiag,         &
       smb_glacier_id_init,              &
       glacier_smb_obs,                  &
       cism_to_rgi_glacier_id,           &       ! diagnostic only
       glacier_area_init,glacier_volume_init, &  ! diagnostic only
       snow_2d,          Tpos_2d,        &
       snow_aux_2d,      Tpos_aux_2d,    &
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
    ! SMB = 0 given input snow_2d and Tpos_2d, for a period with glaciers in balance.
    ! SMB = smb_obs given input snow_aux_2d and Tpos_aux_2d, for a period with glaciers out of balance.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest,         & ! coordinates of diagnostic cell
         nglacier,                    & ! total number of glaciers in the domain
         ngdiag                         ! CISM ID of diagnostic glacier

    integer, dimension(ewn,nsn), intent(in) :: &
         smb_glacier_id_init            ! smb_glacier_id based on the initial glacier extent

    real(dp), dimension(nglacier), intent(in) :: &
         glacier_smb_obs                ! observed glacier-average SMB (mm/yr w.e.)

    integer, dimension(nglacier), intent(in) :: &
         cism_to_rgi_glacier_id         ! RGI glacier ID corresponding to each CISM ID; diagnostic only

    real(dp), dimension(nglacier), intent(in) :: &
         glacier_area_init,           & ! initial glacier area (m^2); diagnostic only
         glacier_volume_init            ! initial glacier volume (m^2); diagnostic only

    real(dp), dimension(ewn,nsn), intent(in) :: &
         snow_2d,                     & ! time-avg snowfall for each cell (mm/yr w.e.)
         Tpos_2d,                     & ! time-avg of max(artm - tmlt, 0) for each cell (deg)
         snow_aux_2d,                 & ! time-avg snowfall for each cell (mm/yr w.e.), auxiliary field
         Tpos_aux_2d                    ! time-avg of max(artm - tmlt, 0) for each cell (deg), auxiliary field

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

    real(dp) :: smb_baseline, smb_aux, smb_aux_diff

    real(dp), dimension(nglacier) :: &
         glacier_snow, glacier_Tpos,          & ! glacier-average snowfall and Tpos
         glacier_snow_aux, glacier_Tpos_aux,  & ! glacier-average snowfall_aux and Tpos_aux
         denom

    character(len=100) :: message

    real(dp), parameter :: Tpos_min = 0.1d0     ! deg C available for melting, min value
                                                ! very low values can resutls in high mu_star

    integer :: count_violate_1, count_violate_2    ! number of glaciers violating Eq. 1 and Eq. 2
    real(dp) :: area_violate_1, area_violate_2     ! total area of these glaciers (m^2)
    real(dp) :: volume_violate_1, volume_violate_2 ! total volume of these glaciers (m^3)
    real(dp) :: mu_eq1, deltaT

    ! Compute mu_star and alpha_snow for each glacier such that
    ! (1) snow and Tpos combine to give SMB = 0
    ! (2) snow_aux and Tpos_aux combine to give SMB = smb_obs
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
    ! (2)       smb_obs = alpha_snow * snow_aux - mu_star * Tpos_aux.
    !
    ! Rearranging and solving, we get
    !              mu_star = (-smb_obs * snow) / D,
    !           alpha_snow = (-smb_obs * Tpos) / D,
    !              where D = snow*Tpos_aux - snow_aux*Tpos
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

    ! Compute average snowfall, Tpos, and SMB over the initial extent of each glacier

    call glacier_2d_to_1d(&
         ewn,           nsn,                   &
         nglacier,      smb_glacier_id_init,   &
         snow_2d,       glacier_snow)

    call glacier_2d_to_1d(&
         ewn,           nsn,                   &
         nglacier,      smb_glacier_id_init,   &
         Tpos_2d,       glacier_Tpos)

    call glacier_2d_to_1d(&
         ewn,           nsn,                   &
         nglacier,      smb_glacier_id_init,   &
         snow_aux_2d,   glacier_snow_aux)

    call glacier_2d_to_1d(&
         ewn,           nsn,                   &
         nglacier,      smb_glacier_id_init,   &
         Tpos_aux_2d,   glacier_Tpos_aux)

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

          ! compute D = snow*Tpos_aux - snow_aux*Tpos
          denom(ng) = glacier_snow(ng)*glacier_Tpos_aux(ng) - glacier_snow_aux(ng)*glacier_Tpos(ng)

          if (glacier_Tpos(ng) < Tpos_min) then

             ! There is little or no ablation anywhere on the glacier in the baseline climate.
             ! Compensate by raising artm (along with artm_aux) until there is some ablation.
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
!                      write(6,'(a46,i6,6f10.3)') 'Out of range, ng, Tp, Tp_aux, D, B, alpha, mu:', &
!                           ng, glacier_Tpos(ng), glacier_Tpos_aux(ng), denom(ng), &
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

       if (verbose_glacier .and. this_rank == rtest .and. ng == ngdiag) then
          print*, ' '
          print*, 'Balance solution, ng =', ng
          print*, '   New mu_star, alpha_snow, beta_artm:', &
               mu_star(ng), alpha_snow(ng), beta_artm(ng)
          print*, '   baseline snow,   Tpos,     smb:', &
               glacier_snow(ng), glacier_Tpos(ng), smb_baseline
          print*, '   recent snow_aux, Tpos_aux, smb:', &
               glacier_snow_aux(ng), glacier_Tpos_aux(ng), smb_aux
          print*, '   smb_aux_diff, smb_obs target   :', &
               smb_aux_diff, glacier_smb_obs(ng)
          print*, ' '
       endif

    enddo   ! ng

    ! Diagnostic checks

    ! Make sure the glacier variables are now in range.
    ! If they are not, there is an error in the logic above.

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

       if (abs(beta_artm(ng)) > beta_artm_max) then
          if (this_rank == rtest) then
             print*, 'WARNING, beta out of range: ng, beta =', ng, beta_artm(ng)
          endif
       endif

    enddo   ! ng

    ! Check the mass balance for the baseline and auxiliary climates.
    ! The goal is that all glaciers satisfy (1), and most satisfy (2).

    count_violate_1 = 0
    count_violate_2 = 0
    area_violate_1 = 0.0d0
    area_violate_2 = 0.0d0
    volume_violate_1 = 0.0d0
    volume_violate_2 = 0.0d0

    do ng = 1, nglacier

       smb_baseline = alpha_snow(ng)*glacier_snow(ng) - mu_star(ng)*glacier_Tpos(ng)
       smb_aux =  alpha_snow(ng)*glacier_snow_aux(ng) - mu_star(ng)*glacier_Tpos_aux(ng)
       smb_aux_diff = smb_aux - glacier_smb_obs(ng)

       if (glacier_Tpos(ng) > 0.0d0) then
          mu_eq1 = alpha_snow(ng) * glacier_snow(ng) / glacier_Tpos(ng)
       else
          mu_eq1 = 0.0d0
       endif

       ! Check whether the glacier violates Eq. (1) and/or Eq. (2)

       if (verbose_glacier .and. this_rank == rtest) then
          if (abs(smb_baseline) > eps08) then
             write(6,'(a60,i6,6f10.2)') 'Eq 1 violation, ng, snow, Tpos, init mu, adj mu, beta, smb :', &
                  ng, glacier_snow(ng), glacier_Tpos(ng), mu_eq1, mu_star(ng), beta_artm(ng), smb_baseline
             count_violate_1 = count_violate_1 + 1
             area_violate_1 = area_violate_1 + glacier_area_init(ng)
             volume_violate_1 = volume_violate_1 + glacier_volume_init(ng)
          endif
          if (abs(smb_aux_diff) > eps08) then
!!             print*, '   Violation of Eq. 2: ng, smb_aux_diff =', ng, smb_aux_diff
             count_violate_2 = count_violate_2 + 1
             area_violate_2 = area_violate_2 + glacier_area_init(ng)
             volume_violate_2 = volume_violate_2 + glacier_volume_init(ng)
          endif
       endif

    enddo  ! ng

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'Violations of Eq. 1:', count_violate_1
       print*, '   Total area, volume =', area_violate_1/1.0d6, volume_violate_1/1.0d9
       print*, 'Violations of Eq. 2:', count_violate_2
       print*, '   Total area, volume =', area_violate_2/1.0d6, volume_violate_2/1.0d9
    endif

    !WHL - debug - Make a list of glaciers with denom and smb_obs having the same sign
!!    if (verbose_glacier .and. this_rank == rtest) then
    if (verbose_glacier .and. 0 == 1) then
       print*, ' '
       print*, 'Glaciers with smb_obs inconsistent with dT = (S/S_aux)*T_aux - T'
       print*, '   ID    RGI_ID    A_init    V_init     snow   snow_aux      Tpos   Tpos_aux      dT    smb_obs'
       do ng = 1, nglacier
          deltaT = denom(ng) / glacier_snow_aux(ng)
          if (glacier_smb_obs(ng) * deltaT > 0.0d0) then
             write(6,'(i6, i10, 8f10.3)') ng, cism_to_rgi_glacier_id(ng), &
                  glacier_area_init(ng)/1.0d6, glacier_volume_init(ng)/1.0d9, &
                  glacier_snow(ng), glacier_snow_aux(ng), &
                  glacier_Tpos(ng), glacier_Tpos_aux(ng), deltaT, glacier_smb_obs(ng)
          endif
       enddo
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

  subroutine glacier_advance_retreat(&
       ewn,             nsn,             &
       itest,   jtest,  rtest,           &
       nglacier,                         &
       glacier_minthck,                  &
       thck,                             &
       smb_annmean,                      &
       snow,                             &
       Tpos,                             &
       mu_star,                          &
       alpha_snow,                       &
       cism_glacier_id_init,             &
       cism_glacier_id,                  &
       parallel)

    ! Allow glaciers to advance and retreat.
    ! This subroutine should be called after the transport/SMB calculation.
    !
    ! The rules are as follows:
    ! - At start-up, glaciated cells have cism_glacier_id in the range (1, nglacier).
    !   Other cells have cism_glacier_id = 0.
    !   The initial cism_glacier_id array is saved as cism_glacier_id_init.
    ! - If a cell has H <= minthck and cism_glacier_id > 0, we set cism_glacier_id = 0.
    !   It no longer contributes to glacier area or volume.
    !   Here, minthck is a threshold for counting ice as part of a glacier.
    !   By default, minthck = model%numerics%thklim, typically 1 m.
    !   (Actually, minthck is slightly less than thklim, to make sure these cells
    !   are not dynamically active.)
    ! - When a cell has H > minthck and cism_glacier_id = 0, we give it a nonzero ID:
    !   either (1) cism_glacier_id_init, if the initial ID > 0,
    !   or (2) the ID of an adjacent glaciated neighbor (the one where the cell would
    !   have the most negative SMB, if there is more than one).
    !   Preference is given to (1), to preserve the original glacier outlines
    !   as much as possible.

    use cism_parallel, only: parallel_globalindex, parallel_halo

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         itest, jtest, rtest,         & ! coordinates of diagnostic cell
         nglacier                       ! number of glaciers

    real(dp), intent(in) :: &
         glacier_minthck                ! min ice thickness (m) counted as part of a glacier

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         thck,                        & ! ice thickness (m)
         smb_annmean,                 & ! annual mean SMB (mm/yr w.e.)
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
         smb_min,                    &  ! min SMB possible for this cell
         smb_neighbor                   ! SMB that a cell would have in a neighbor glacier
                                        ! (due to different alpha_snow and mu_star)

    character(len=100) :: message

    integer :: i, j, ii, jj, ip, jp
    integer :: iglobal, jglobal
    integer :: ng, ng_init, ng_neighbor, ng_min
    logical :: found_neighbor

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glacier_advance_retreat'
    endif

    ! Check for retreat: cells with cism_glacier_id > 0 but H > glacier_minthck

!    do j = nhalo+1, nsn-nhalo
!       do i = nhalo+1, ewn-nhalo
!          ng = cism_glacier_id_init(i,j)
!          if (ng == 3651) then
!             call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!             print*, 'Glacier 3651: ig, jg =', iglobal, jglobal
!          endif
!       enddo
!    enddo

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


    ! Put the cell in the glacier that gives it the lowest SMB, given its own snow and Tpos.

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng = cism_glacier_id_old(i,j)
          ng_init = cism_glacier_id_init(i,j)

          if (ng == 0 .and. thck(i,j) > glacier_minthck) then
             ! assign this cell its original ID, if > 0
             if (ng_init > 0) then
                cism_glacier_id(i,j) = ng_init
                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, 'Set ID = init ID: ig, jg, new ID, thck =',&
                        iglobal, jglobal, cism_glacier_id(i,j), thck(i,j)
                endif
             else  ! assign the ID of an adjacent ice-covered cell, if possible

                smb_min = 1.0d11   ! arbitrary big number
                ng_min = 0
                found_neighbor = .false.

                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, 'Look for neighbor for cell: ig, jg, rank, i, j =', &
                        iglobal, jglobal, this_rank, i, j
                endif

                do jj = -1, 1
                   do ii = -1, 1
                      if (ii /= 0 .or. jj /= 0) then  ! one of 8 neighbors
                         ip = i + ii
                         jp = j + jj
                         ng_neighbor = cism_glacier_id_old(ip,jp)
                         !TODO - Do we need the thickness criterion?
                         if (ng_neighbor > 0 .and. thck(ip,jp) > glacier_minthck) then
                            found_neighbor = .true.
                            ! Compute the SMB this cell would have if in the neighbor glacier
                            smb_neighbor = alpha_snow(ng_neighbor) * snow(i,j) &
                                            - mu_star(ng_neighbor) * Tpos(i,j)
                            if (smb_neighbor < smb_min) then
                               smb_min = smb_neighbor
                               ng_min = ng_neighbor
                            endif
                         endif   ! neighbor cell is a glacier cell
                      endif   ! neighbor cell
                   enddo   ! ii
                enddo   ! jj
                if (found_neighbor) then
                   cism_glacier_id(i,j) = ng_min
                   if (verbose_glacier .and. this_rank == rtest) then
                      call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                      print*, 'Set ID = neighbor ID, ig, jg, new ID, thck, smb =', &
                           iglobal, jglobal, cism_glacier_id(i,j), thck(i,j), smb_min
                   endif
                else
                   !Note: This can happen if an advanced cell has a more positive SMB than its neighbor,
                   !      and the neighbor melts. We want to remove this cell from the glacier.
                   ! For now, remove ice from this cell.
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, 'WARNING, did not find neighbor, ig, jg =', iglobal, jglobal
                endif   ! found_neighbor
             endif   ! cism_glacier_id_init > 0

          endif   ! ng = 0, H > minthck
       enddo   ! i
    enddo   ! j

    call parallel_halo(cism_glacier_id, parallel)

    ! Check glacier IDs at the margin, outside the initial footprint.
    ! Switch IDs that are potentially problematic.
    !
    ! The code below protects against glacier 'pirating'.
    ! This can happen when two adjacent glaciers have both advanced: one with a large ablation rate
    !  and the other with a lower ablation rate. The SMBs favor advance of the slow-melting glacier
    !  at the expense of the fast-melting glacier. The fast-melting glacier can feed ice
    !  into the slow-melting glacier, leading to spurious advance of the slow-melting glacier.
    ! The fix here is to loop through cells where the ice has advanced (cism_glacier_id_init = 0,
    !  cism_glacier_id > 0). For each cell, check whether it has a neighbor in a different glacier.
    !  If so, compute the SMB it would have in that glacier, given a different value of alpha_snow
    !  and mu_star. If this SMB is negative and lower than the current value, make the switch.
    ! TODO - Check for unrealistic glacier expansion.
    ! Note: This should happen early in the spin-up, not as the run approaches steady state.

    ! Save a copy of the current cism_glacier_id.
    cism_glacier_id_old = cism_glacier_id

    ! Loop over local cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          ng_init = cism_glacier_id_init(i,j)
          ng = cism_glacier_id_old(i,j)
          if (ng_init == 0 .and. ng > 0) then ! advanced cell
             smb_min = min(smb_annmean(i,j), 0.0d0)
             ng_min = 0

             ! Look for edge neighbors in different glaciers
             do jj = -1, 1
                do ii = -1, 1
                   if ((abs(ii)==1 .and. jj==0) .or. (abs(jj)==1 .and. ii==0)) then  ! edge neighbor
                      ip = i + ii
                      jp = j + jj
                      ng_neighbor = cism_glacier_id_old(ip,jp)

                      if (ng_neighbor > 0 .and. ng_neighbor /= ng) then  ! different glacier

                         if (verbose_glacier .and. this_rank == rtest) then
                            call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                            print*, 'Check neighbor SMB for cell', iglobal, jglobal
                            print*, '   Local ng, neighbor ng =', ng, ng_neighbor
                         endif

                         ! compute the SMB of cell (i,j) if moved to the neighbor glacier
                         smb_neighbor = alpha_snow(ng_neighbor) * snow(i,j) &
                                         - mu_star(ng_neighbor) * Tpos(i,j)
                         if (verbose_glacier .and. this_rank == rtest) then
                            print*, '   Local SMB, SMB if in neighbor glacier =', smb_annmean(i,j), smb_neighbor
                         endif
                         if (smb_neighbor < smb_min) then
                            smb_min = smb_neighbor
                            ng_min = ng_neighbor
                         endif
                      endif
                   endif   ! neighbor cell
                enddo   ! ii
             enddo   ! jj

             if (ng_min > 0) then
                ! Move this cell to the adjacent glacier, where it will melt faster
                cism_glacier_id(i,j) = ng_min
                if (verbose_glacier .and. this_rank == rtest) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, '   Transfer to fast-melting glacier, old and new IDs =', &
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
         smb_glacier_id_init,     &
         smb_glacier_id,          &
         parallel)

    ! Based on the current glacier footprint, compute a mask of cells that can have a nonzero SMB.
    !
    ! The rules for smb_glacier_id are as follows:
    ! - Where cism_glacier_id_init > 0, set smb_glacier_id(i,j) = cism_glacier_id(i,j)
    !   and apply the SMB.
    !   Note: In ice-free retreated cells (cism_glacier_id_init > 0 but cism_glacier_id = 0),
    !   the negative SMB will be ignored.
    ! - In advanced grid cells (cism_glacier_id_init = 0 but cism_glacier_id > 0),
    !   compute a potential SMB assuming smb_glacier_id(i,j) = cism_glacier_id(i,j).
    !   Apply this SMB if negative; else set smb_glacier_id(i,j) = 0.
    ! - In other glacier-free cells (cism_glacier_id_init = cism_glacier_id = 0), check
    !   for glacier-covered edge neighbors (cism_glacier_id > 0). For each neighbor (ii,jj),
    !   compute a potential SMB assuming smb_glacier_id(i,j) = cism_glacier_id(ii,jj).
    !   Apply this SMB if negative; else set smb_glacier_id(i,j) = 0.
    !   If there are neighbors with SMB < 0 from two or more glaciers, choose the glacier ID
    !   that results in the lowest SMB.
    !
    ! The rules for smb_glacier_id_init are the same as for smb_glacier_id, except that
    !  we assume cism_glacier_id = cism_glacier_id_init, so there are no advanced
    !  or retreated cells.
    !
    ! The goal is to spin up each glacier to an extent similar to the observed extent,
    ! using a mask to limit expansion but without using fictitious SMB values.

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
         smb_glacier_id_init,         & ! integer glacier ID used for SMB calculations, based on initial extent
         smb_glacier_id                 ! integer glacier ID in the range (1, nglacier), based on current extent
                                        ! = 0 in cells where we force SMB = 0

    type(parallel_type), intent(in) :: parallel

    ! local variables
    integer :: i, j, ii, jj, ng, ng_min
    integer :: ip, jp
    integer :: iglobal, jglobal

    real(dp) :: &
         smb_potential,               & ! potential SMB in a given cell outside the initial footprint
         smb_min                        ! min value of SMB for a given cell with glacier-covered neighbors

    ! Initialize the SMB masks
    smb_glacier_id = 0

    ! Compute smb_glacier_id

    ! First, set smb_glacier_id = cism_glacier_id_init
    smb_glacier_id = cism_glacier_id_init

    ! Extend smb_glacier_id to advanced cells with SMB < 0.

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          if (cism_glacier_id_init(i,j) == 0 .and. cism_glacier_id(i,j) > 0) then ! advanced cell
             ! compute the potential SMB for this cell; apply if negative
             ng = cism_glacier_id(i,j)
             smb_potential = alpha_snow(ng)*snow(i,j) - mu_star(ng)*Tpos(i,j)
             if (smb_potential < 0.0d0) smb_glacier_id(i,j) = ng
          endif
       enddo
    enddo

    ! Where cism_glacier_id_init = cism_glacier_id = 0, look for neighbors with cism_glacier_id > 0 and SMB < 0.
    ! Extend smb_glacier_id to these cells.

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          if (cism_glacier_id_init(i,j) == 0 .and. cism_glacier_id(i,j) == 0) then ! glacier-free cell
             ! find the adjacent glacier-covered cell (if any) with the most negative SMB
             smb_min = 0.0d0
             ng_min = 0
             do jj = -1,1
                do ii = -1,1
                   if (ii /= 0 .or. jj /= 0) then  ! edge or diagonal neighbor
                      ip = i + ii
                      jp = j + jj
                      if (cism_glacier_id(ip,jp) > 0) then  ! adjacent glacier
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
             ! If there are any adjacent glacier cells with SMB < 0, add cell (i,j) to the mask
             if (ng_min > 0) then
                smb_glacier_id(i,j) = ng_min
!                if (verbose_glacier .and. this_rank == rtest) then
!                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!                   print*, 'Set smb_glacier_id = neighbor ID: ig, jg, smb_min, upstream ID =', &
!                        iglobal, jglobal, smb_min, smb_glacier_id(i,j)
!                endif
             endif
          endif   ! cism_glacier_id_init = cism_glacier_id = 0
       enddo   ! i
    enddo   ! j

    ! Compute smb_glacier_id_init

    ! First, set smb_glacier_id_init = cism_glacier_id_init
    smb_glacier_id_init = cism_glacier_id_init

    ! Where cism_glacier_id_init = 0, look for neighbors with cism_glacier_id_init > 0 and SMB < 0.
    ! Extend smb_glacier_id_init to these cells.

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          if (cism_glacier_id_init(i,j) == 0) then ! initially glacier-free cell
             ! find the adjacent glacier-covered cell (if any) with the most negative SMB
             smb_min = 0.0d0
             ng_min = 0
             do jj = -1,1
                do ii = -1,1
                   if (ii /= 0 .or. jj /= 0) then  ! edge or diagonal neighbor
                      ip = i + ii
                      jp = j + jj
                      if (cism_glacier_id_init(ip,jp) > 0) then  ! adjacent glacier
                         ng = cism_glacier_id_init(ip,jp)
                         ! compute the potential SMB, assuming cell (i,j) is in glacier ng
                         smb_potential = alpha_snow(ng)*snow(i,j) - mu_star(ng)*Tpos(i,j)
                         if (smb_potential < smb_min) then
                            smb_min = smb_potential
                            ng_min = ng
                         endif
                      endif   ! cism_glacier_id_init > 0
                   endif   ! neighbor cell
                enddo   ! ii
             enddo   ! jj
             ! If there are any adjacent glacier cells with SMB < 0, add cell (i,j) to the mask
             if (ng_min > 0) then
                smb_glacier_id_init(i,j) = ng_min
!                if (verbose_glacier .and. this_rank == rtest) then
!                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!                   print*, 'Set smb_glacier_id_init = neighbor ID: ig, jg, smb_min, upstream ID =', &
!                        iglobal, jglobal, smb_min, smb_glacier_id_init(i,j)
!                endif
             endif
          endif   ! cism_glacier_id_init = 0
       enddo   ! i
    enddo   ! j

    call parallel_halo(smb_glacier_id, parallel)
    call parallel_halo(smb_glacier_id_init, parallel)

  end subroutine update_smb_glacier_id

!****************************************************

  subroutine remove_snowfields(&
       ewn,          nsn,           &
       parallel,                    &
       itest, jtest, rtest,         &
       thck,                        &
       cism_glacier_id_init,        &
       cism_glacier_id)

    ! This subroutine is patterned after subroutine remove_icebergs in the calving module.
    ! A snowfield is defined as an isolated patch of glacier ice outside the initial glacier footprint
    !  (as defined by cism_glacier_id_init).
    ! If it becomes disconnected from the main glacier, it is removed.
    !
    ! The algorithm is as follows:
    ! (1) Mark all cells with ice (either active or inactive) with the initial color.
    !     Mark other cells with the boundary color.
    ! (2) Seed the fill by giving the fill color to active glacier cells (cism_glacier_id = 1)
    !     that are part of the initial glacier (cism_glacier_id_init = 1).
    ! (3) Recursively fill all cells that are connected to filled cells by a path
    !     that passes only through active glacier cells.
    ! (4) Repeat the recursion as necessary to spread the fill to adjacent processors.
    ! (5) Once the fill is done, any ice-covered cells that still have the initial color
    !     are considered to be isolated snowfields and are removed.
    !
    ! Notes:
    ! (1) The recursive fill applies to edge neighbors, not corner neighbors.
    !     The path back to the initial glacier must go through edges, not corners.
    ! (2) Inactive cells (thck < glacier%minthck) can be filled if adjacent to active cells, but
    !     do not further spread the fill.

    use glissade_masks, only: glissade_fill_with_buffer, initial_color, fill_color, boundary_color
    use cism_parallel, only: parallel_halo, parallel_reduce_sum, parallel_globalindex

    integer, intent(in) :: ewn, nsn                              !> horizontal grid dimensions
    type(parallel_type), intent(in) :: parallel                  !> info for parallel communication
    integer, intent(in) :: itest, jtest, rtest                   !> coordinates of diagnostic point

    real(dp), dimension(ewn,nsn), intent(inout) :: thck          !> ice thickness

    integer,  dimension(ewn,nsn), intent(in) :: &
         cism_glacier_id_init

    integer,  dimension(ewn,nsn), intent(inout) :: &
         cism_glacier_id

    ! local variables

    real(dp) :: dthck

    integer :: i, j, iglobal, jglobal

    integer :: &
         iter,                      & ! iteration counter
         max_iter,                  & ! max(ewtasks, nstasks)
         local_count,               & ! local counter for filled values
         global_count,              & ! global counter for filled values
         global_count_save            ! globalcounter for filled values from previous iteration

    integer, dimension(ewn,nsn) ::  &
         cism_glacier_mask_init,    & ! = 1 where cism_glacier_id_init > 0, else = 0
         cism_glacier_mask,         & ! = 1 where cism_glacier_id > 0, else = 0
         color                        ! integer 'color' for identifying snowfields

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'In remove_snowfields'
       print*, ' '
       print*, 'thck, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'cism_glacier_id_init:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') cism_glacier_id_init(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'cism_glacier_id:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') cism_glacier_id(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Initialize snowfield removal
    ! Note: Any cell with ice, active or inactive, receives the initial color.
    !       Inactive cells can later receive the fill color (if adjacent to active cells)
    !        but cannot further spread the fill color.
    !       This protects inactive, thickening cells at the glacier margin from being removed
    !        before they can activate.

    do j = 1, nsn
       do i = 1, ewn
          if (thck(i,j) > 0.0d0) then
             color(i,j) = initial_color
          else
             color(i,j) = boundary_color
          endif
       enddo
    enddo

    where (cism_glacier_id_init > 0)
       cism_glacier_mask_init = 1
    elsewhere
       cism_glacier_mask_init = 0
    endwhere

    where (cism_glacier_id > 0)
       cism_glacier_mask = 1
    elsewhere
       cism_glacier_mask = 0
    endwhere

    ! Loop through cells, identifying active glacier cells with cism_glacier_id_init = 1.
    ! Fill each such cell, and then recursively fill active neighbor cells (cism_glacier_id = 1).
    ! We may have to do this several times to incorporate connections between neighboring processors.

    max_iter = max(parallel%ewtasks, parallel%nstasks)
    global_count_save = 0

    do iter = 1, max_iter

       if (iter == 1) then   ! identify active glacier cells that can seed the fill

          do j = 1, nsn
             do i = 1, ewn

                ! Fill active glacier cells that are part of the initial glacier.

                if (cism_glacier_mask_init(i,j) == 1 .and. cism_glacier_mask(i,j) == 1) then

                   if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then

                      ! assign the fill color to this cell, and recursively fill neighbor cells
                      call glissade_fill_with_buffer(ewn,   nsn,   &
                                                     i,     j,     &
                                                     color, cism_glacier_mask)

                   endif

                endif
             enddo
          enddo

       else  ! count > 1

          ! Check for halo cells that were just filled on neighbor processors
          ! Note: In order for a halo cell to seed the fill on this processor, it must not only have the fill color,
          !       but also must be an active cell.

          call parallel_halo(color, parallel)

          ! west halo layer
          i = nhalo
          do j = 1, nsn
             if (color(i,j) == fill_color .and. cism_glacier_id(i,j) == 1) then
                call glissade_fill_with_buffer(ewn,   nsn,   &
                                               i+1,   j,     &
                                               color, cism_glacier_mask)
             endif
          enddo

          ! east halo layers
          i = ewn - nhalo + 1
          do j = 1, nsn
             if (color(i,j) == fill_color .and. cism_glacier_id(i,j) == 1) then
                call glissade_fill_with_buffer(ewn,   nsn,   &
                                               i-1,   j,     &
                                               color, cism_glacier_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, ewn-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. cism_glacier_id(i,j) == 1) then
                call glissade_fill_with_buffer(ewn,   nsn,   &
                                               i,     j+1,   &
                                               color, cism_glacier_mask)
             endif
          enddo

          ! north halo layer
          j = nsn-nhalo+1
          do i = nhalo+1, ewn-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. cism_glacier_id(i,j) == 1) then
                call glissade_fill_with_buffer(ewn,   nsn,   &
                                               i,     j-1,   &
                                               color, cism_glacier_mask)
             endif
          enddo

       endif  ! count = 1

       local_count = 0
       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             if (color(i,j) == fill_color) local_count = local_count + 1
          enddo
       enddo

       !WHL - If running a large problem, may want to reduce the frequency of this global sum
       global_count = parallel_reduce_sum(local_count)

       if (global_count == global_count_save) then
          if (verbose_glacier .and. main_task) &
               print*, 'Fill converged: iter, global_count =', iter, global_count
          exit
       else
          if (verbose_glacier .and. main_task) &
               print*, 'Convergence check: iter, global_count =', iter, global_count
          global_count_save = global_count
       endif

    enddo  ! count

    ! Snowfields are cells that still have the initial color and are not on land.
    ! Remove ice in these cells.
    ! TODO: How to conserve mass while doing this? Need to update acab?

    do j = 2, nsn-1
       do i = 2, ewn-1
          if (color(i,j) == initial_color) then
             if (cism_glacier_id(i,j) > 0) then
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                print*, 'Snowfield: Set cism_glacier_id = 0, ig, jg, ng, thck =', &
                     iglobal, jglobal, cism_glacier_id(i,j), thck(i,j)
             endif
             cism_glacier_id(i,j) = 0
             dthck = thck(i,j)
             thck(i,j) = 0.0d0
             !TODO - Also handle tracers?  E.g., set damage(:,i,j) = 0.d0?
          endif
       enddo
    enddo

    call parallel_halo(thck, parallel)
    call parallel_halo(cism_glacier_id, parallel)

    if (verbose_glacier .and. this_rank == rtest) then
       print*, ' '
       print*, 'Done in remove_snowfields'
       print*, ' '
       print*, 'thck, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine remove_snowfields

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
       area_factor,                      &
       area,          volume)

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier)

    real(dp), intent(in) :: &
         cell_area                      ! grid cell area (m^2), dew*dns, assumed equal for all cells

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         thck,                        & ! ice thickness (m)
         area_factor                    ! scale factor multiplying the nominal cell area, based on latitude

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
                local_area(ng) = local_area(ng) + cell_area*area_factor(i,j)
                local_volume(ng) = local_volume(ng) + cell_area*area_factor(i,j) * thck(i,j)
             endif
          endif
       enddo
    enddo

    area   = parallel_reduce_sum(local_area)
    volume = parallel_reduce_sum(local_volume)

    if (verbose_glacier .and. main_task) then
       print*, ' '
       print*, 'Compute glacier area and volume'
       print*, 'Max area (km^2)   =', maxval(area) * 1.0d-6    ! m^2 to km^2
       print*, 'Max volume (km^3) =', maxval(volume) * 1.0d-9  ! m^3 to km^3
       print*, ' '
    endif

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
    ! Note: In this subroutine the cell area is not corrected using an area scale factor.
    !       We assume all cells have equal area, cell_area = dew*dns.

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id_init,        & ! integer glacier ID in the range (1, nglacier), initial value
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier), current value

    real(dp), intent(in) :: &
         cell_area                      ! grid cell area (m^2), assumed equal for all cells

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
       cism_glacier_id_init,   &
       cism_glacier_id,        &
       cell_area,              &
       smb_annmean,            &
       aar_init,               &
       aar)

    ! Compute the accumulation area ratio (AAR) for each glacier.
    ! Note: In this subroutine the cell area is not corrected using an area scale factor.
    !       We assume all cells have equal area, cell_area = dew*dns.

    use cism_parallel, only: parallel_reduce_sum

    ! input/output arguments

    integer, intent(in) ::  &
         ewn, nsn,                    & ! number of cells in each horizontal direction
         nglacier                       ! total number of glaciers in the domain

    integer, dimension(ewn,nsn), intent(in) ::  &
         cism_glacier_id_init,        & ! integer glacier ID in the range (1, nglacier), initial value
         cism_glacier_id                ! integer glacier ID in the range (1, nglacier), current value

    real(dp), intent(in) :: &
         cell_area                      ! grid cell area (m^2), assumed equal for all cells

    real(dp), dimension(ewn,nsn), intent(in) ::  &
         smb_annmean                    ! 2D annual mean SMB (mm/yr w.e.)

    real(dp), dimension(nglacier), intent(out) ::  &
         aar_init,                    & ! AAR over the initial glacier area
         aar                            ! AAR over the current glacier area

    ! local variables

    integer :: i, j, ng

    real(dp), dimension(nglacier) :: &
         area_init, area, &
         accum_area_init, accum_area

    ! initialize
    area_init(:) = 0.0d0
    area(:) = 0.0d0
    accum_area_init(:) = 0.0d0
    accum_area(:) = 0.0d0

    ! Compute the accumulation area and total area for each glacier

    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo

          ! initial glacier ID
          ng = cism_glacier_id_init(i,j)
          if (ng > 0) then
             area_init(ng) = area_init(ng) + cell_area
             if (smb_annmean(i,j) >= 0.0d0) then
                accum_area_init(ng) = accum_area_init(ng) + cell_area
             endif
          endif

          ! current glacier ID
          ng = cism_glacier_id(i,j)
          if (ng > 0) then
             area(ng) = area(ng) + cell_area
             if (smb_annmean(i,j) >= 0.0d0) then
                accum_area(ng) = accum_area(ng) + cell_area
             endif
          endif

       enddo   ! i
    enddo   ! j

    area_init = parallel_reduce_sum(area_init)
    area = parallel_reduce_sum(area)
    accum_area_init = parallel_reduce_sum(accum_area_init)
    accum_area = parallel_reduce_sum(accum_area)

    ! Compute the AAR for each glacier

    where (area_init > 0.0d0)
       aar_init = accum_area_init / area_init
    elsewhere
       aar_init = 0.0d0
    endwhere

    where (area > 0.0d0)
       aar = accum_area / area
    elsewhere
       aar = 0.0d0
    endwhere

  end subroutine glacier_accumulation_area_ratio

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
