!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade.F90 - part of the Community Ice Sheet Model (CISM)  
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

! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

!CLEANUP - glissade.F90
!
! NOTE: MJH Lines that start with !### are ones I have identified to be deleted.
!
! This module was originally copied from glide.F90 (William Lipscomb, June 2012)
! Removed SIA-specific code, leaving only the HO code with remapping transport.
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glissade.f90 - part of the CISM ice model        + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade

  ! Driver for Glissade (parallel, higher-order) dynamical core

  use glimmer_global, only: dp
  use glimmer_log
  use glide_types
  use glide_io
  use glide_lithot
  use glimmer_config
  use glissade_test, only: &
       glissade_test_halo, glissade_test_transport, glissade_test_comm_row_col
  use glide_thck, only: glide_calclsrf  ! TODO - Make this a glissade subroutine, or inline
  use profile, only: t_startf, t_stopf
  use cism_parallel, only: this_rank, main_task, comm, nhalo

  implicit none

  integer, private, parameter :: dummyunit=99
  logical, parameter :: verbose_glissade = .false.
!!  logical, parameter :: verbose_glissade = .true.

  ! Change any of the following logical parameters to true to carry out simple tests
  logical, parameter :: test_transport = .false.    ! if true, call test_transport subroutine
  real(dp), parameter :: thk_init = 500.d0          ! initial thickness (m) for test_transport
  logical, parameter :: test_halo = .false.         ! if true, call test_halo subroutine
  logical, parameter :: test_comm_row_col = .false. ! if true, test the row and column communicators

contains

!=======================================================================

! Note: There is no glissade_config subroutine; glide_config works for all dycores.

!=======================================================================

  subroutine glissade_initialise(model, evolve_ice)

    ! initialise Glissade model instance

    use cism_parallel, only: parallel_type, distributed_gather_var,  &
         distributed_scatter_var, parallel_finalise, &
         distributed_grid, distributed_grid_active_blocks,  parallel_global_edge_mask, &
         parallel_halo, parallel_halo_extrapolate, parallel_reduce_max, &
         staggered_parallel_halo_extrapolate, staggered_no_penetration_mask, &
         parallel_create_comm_row, parallel_create_comm_col, not_parallel

    use glide_stop, only: register_model
    use glide_setup
    use glimmer_ncio, only: openall_in, openall_out, glimmer_nc_get_var, glimmer_nc_get_dimlength
    use glide_velo, only: init_velo  !TODO - Remove call to init_velo?
    use glissade_therm, only: glissade_init_therm
    use glissade_transport, only: glissade_overwrite_acab_mask, glissade_add_2d_anomaly
    use glissade_basal_water, only: glissade_basal_water_init
    use glissade_masks, only: glissade_get_masks, glissade_marine_connection_mask
    use glimmer_scales
    use glimmer_paramets, only: eps11, thk0, len0, tim0, vel0, scyr
    use glimmer_physcon, only: rhow, rhoi
    use glide_mask
    use isostasy, only: init_isostasy, isos_relaxed
    use glimmer_map_init
    use glimmer_coordinates, only: coordsystem_new
    use glissade_grid_operators, only: glissade_stagger, glissade_laplacian_smoother
    use glissade_velo_higher, only: glissade_velo_higher_init
    use glide_diagnostics, only: glide_init_diag
    use glissade_calving, only: glissade_calving_mask_init, glissade_thck_calving_threshold_init
    use glissade_inversion, only: glissade_init_inversion, verbose_inversion
    use glissade_basal_traction, only: glissade_init_effective_pressure
    use glissade_bmlt_float, only: glissade_bmlt_float_thermal_forcing_init, verbose_bmlt_float
    use glissade_grounding_line, only: glissade_grounded_fraction
    use glissade_glacier, only: glissade_glacier_init
    use glissade_utils, only: glissade_adjust_thickness, glissade_smooth_usrf, &
         glissade_smooth_topography, glissade_adjust_topography
    use glissade_utils, only: glissade_stdev
    use felix_dycore_interface, only: felix_velo_init

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    logical, intent(in), optional :: evolve_ice       ! whether ice evolution is turned on (if not present, assumed true)

    !TODO - Is glimmer_version_char still needed?
    character(len=100), external :: glimmer_version_char

    character(len=100) :: message

    real(dp) :: local_maxval, global_maxval   ! max values of a given variable; = 0 if not yet read in
    integer :: i, j, k
    logical :: l_evolve_ice  ! local version of evolve_ice

    integer, dimension(:,:), allocatable :: &
         ice_mask,          & ! = 1 where ice is present, else = 0
         floating_mask,     & ! = 1 where ice is present and floating, else = 0
         ocean_mask,        & ! = 1 if topg is below sea level and ice is absent, else = 0
         land_mask            ! = 1 if topg is at or above sea level, else = 0

    real(dp), dimension(:,:), allocatable :: &
         topg_smoothed,     & ! smoothed input topography
         thck_flotation       ! flotation thickness

    integer, dimension(:,:), allocatable :: &
         ice_domain_mask      ! = 1 where ice is potentially present and active

    logical, parameter :: &
         make_ice_domain_mask = .false.   ! set to .true. to create mask at initialization
!!         make_ice_domain_mask = .true.   ! set to .true. to create mask at initialization

    real(dp) :: usrf_max    ! max value of usrf
    real(dp) :: topg        ! model%geometry%topg - model%climate%eus
    real(dp) :: thck_flot   ! flotation thickness

    integer :: itest, jtest, rtest
    integer :: status, varid

    type(glimmer_nc_input), pointer :: infile
    type(parallel_type) :: parallel   ! info for parallel communication

    !WHL - added for optional topg_stdev calculations
    logical, parameter :: compute_topg_stdev = .false.
    real(dp), dimension(:,:), allocatable :: topg_global, topg_stdev_global
    integer, parameter :: grid_ratio = 8   ! ratio between no. of grid cells on the input grid and
                                           ! the grid on which we want to compute topg_stdev; typically 2, 4, or 8

    if (main_task) print*, 'In glissade_initialise'

    if (present(evolve_ice)) then
       l_evolve_ice = evolve_ice
    else
       l_evolve_ice = .true.
    end if

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters
    call glide_scale_params(model)

    ! Set up coordinate systems, and change the parallel values of ewn and nsn.
    ! With no_ice BCs, scalars adjacent to the global boundary (including halos) are set to zero.

    if (model%options%compute_blocks == ACTIVE_BLOCKS_ONLY .or.   &
        model%options%compute_blocks == ACTIVE_BLOCKS_INQUIRE) then

       ! Allocate memory for a global array, ice_domain_mask, on main_task.

       if (associated(model%general%ice_domain_mask)) &
            deallocate(model%general%ice_domain_mask)

       if (main_task) then
          allocate(model%general%ice_domain_mask(model%general%ewn, model%general%nsn))
       else
          allocate(model%general%ice_domain_mask(1,1))
       endif
       model%general%ice_domain_mask = 0

       ! Read ice_domain_mask from the input or restart file
       ! Note: In general, input arrays are read from subroutine glide_io_readall (called below) in glide_io.F90.
       !       However, ice_domain_mask is needed now to identify active blocks.

       infile => model%funits%in_first   ! assume ice_domain_mask is in the input or restart file

       call glimmer_nc_get_var(infile, 'ice_domain_mask', &
                               model%general%ice_domain_mask)

       if (model%options%compute_blocks == ACTIVE_BLOCKS_INQUIRE) then

          ! The subroutine will report how many tasks are needed to compute on all active blocks, and then abort.
          ! The user can then resubmit (on an optimal number of processors) with model%options%compute_blocks = ACTIVE_BLOCKS.

          call distributed_grid_active_blocks(model%general%ewn,      model%general%nsn,      &
                                              model%general%nx_block, model%general%ny_block, &
                                              model%general%ice_domain_mask,                  &
                                              model%parallel,                                 &
                                              inquire_only = .true.)

       else  ! compute_blocks = ACTIVE_BLOCKS_ONLY

          ! Set up a distributed grid with computations on active blocks only.
          ! An active block contains one or more cells with ice_domain_mask = 1.
          ! This option is supported only with outflow BCs.

          if (model%general%global_bc /= GLOBAL_BC_NO_ICE) then
             call write_log('Changing to outflow boundary conditions to support the active_blocks option')
             model%general%global_bc = GLOBAL_BC_NO_ICE
          endif

          call distributed_grid_active_blocks(model%general%ewn,      model%general%nsn,      &
                                              model%general%nx_block, model%general%ny_block, &
                                              model%general%ice_domain_mask,                  &
                                              model%parallel)

       endif   ! compute_blocks

       ! Nullify ice_domain_mask so it can be read again and scattered below.
       nullify(model%general%ice_domain_mask)

    elseif (model%general%global_bc == GLOBAL_BC_OUTFLOW) then

       call distributed_grid(model%general%ewn, model%general%nsn, &
                             model%parallel,    global_bc_in = 'outflow')

    elseif (model%general%global_bc == GLOBAL_BC_NO_ICE) then

       call distributed_grid(model%general%ewn, model%general%nsn, &
                             model%parallel,     global_bc_in = 'no_ice')

    elseif (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then

       ! Note: In this case, halo updates are the same as for periodic BC.
       !       The difference is that we also use no-penetration masks for (uvel,vvel) at the global boundary
       !       (computed by calling staggered_no_penetration_mask below).

       call distributed_grid(model%general%ewn, model%general%nsn, &
                             model%parallel,     global_bc_in = 'no_penetration')

    else  ! global_bc = GLOBAL_BC_PERIODIC

!       call distributed_grid(model%general%ewn, model%general%nsn, global_bc_in = 'periodic')

       call distributed_grid(model%general%ewn, model%general%nsn, &
                             model%parallel,  global_bc_in = 'periodic')

    endif

    if (model%options%which_ho_precond == HO_PRECOND_TRIDIAG_GLOBAL) then

       ! Set up row-based and column-based communicators (in addition to mpi_comm_world).
       ! These communicators are used to solve tridiagonal matrix problems in parallel along global rows and columns.
       ! For the row-based communicator, the task with the minimum rank in each row becomes main_task_row.
       ! For the column-based communicator, the task with the minimum rank in each column becomes main_task_column.
       call parallel_create_comm_row(comm, model%parallel)
       call parallel_create_comm_col(comm, model%parallel)

       if (test_comm_row_col) then
          call glissade_test_comm_row_col(model)
       endif

    endif  ! HO_PRECOND_TRIDIAG_GLOBAL

    ! Now that model%parallel has been set, copy to 'parallel' to save typing below
    parallel = model%parallel

    model%general%ice_grid = coordsystem_new(0.d0,               0.d0,               &
                                             model%numerics%dew, model%numerics%dns, &
                                             model%general%ewn,  model%general%nsn)

    model%general%velo_grid = coordsystem_new(model%numerics%dew/2.d0, model%numerics%dns/2.d0, &
                                              model%numerics%dew,      model%numerics%dns,      &
                                              model%general%ewn-1,     model%general%nsn-1)

    ! If the length of any dimension is unknown, then get the length now, before allocating arrays.
    ! Currently, the length of most dimensions is set in the config file.
    ! An exception is dimension glacierid, whose length (nglacier) is computed internally by CISM.
    ! On restart, we can get the length from the restart file.

    if (model%options%enable_glaciers .and. model%options%is_restart == RESTART_TRUE) then
       infile => model%funits%in_first   ! assume glacierid is a dimension in the restart file
       call glimmer_nc_get_dimlength(infile, 'glacierid', model%glacier%nglacier)
    endif

    ! allocate arrays
    call glide_allocarr(model)

    ! Compute a mask to identify cells at the edge of the global domain
    ! (Currently used only to compute bwat_mask for basal water routing)
    ! Includes a halo update for global_edge_mask
    call parallel_global_edge_mask(model%general%global_edge_mask, parallel)

    ! Set masks at global boundary for no-penetration boundary conditions
    ! Includes a halo update for the masks
    if (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then
       call staggered_no_penetration_mask(model%velocity%umask_no_penetration, &
                                          model%velocity%vmask_no_penetration, &
                                          parallel)
    endif

    ! set uniform basal heat flux (positive down)
    model%temper%bheatflx = model%paramets%geot

    ! compute sigma levels or load from external file
    ! (if not already read from config file)
    call glide_load_sigma(model,dummyunit)

    ! initialize the time step counter
    ! For restart, tstep_count will be overwritten from the restart file.
    ! Alternatively, we could initialize tstep_count as follows:
    !    model%numerics%tstep_count = nint(model%numerics%time/model%numerics%tinc)
    ! But reading from a restart file should prevent roundoff issues.
    model%numerics%tstep_count = 0

    ! open all input files
    call openall_in(model)

    ! read first time slice
    call glide_io_readall(model,model)

    ! initialize model diagnostics
    call glide_init_diag(model)

    ! Set coordinates of diagnostic point
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    ! Check that lat and lon fields were read in, if desired
    if (model%options%read_lat_lon) then
       local_maxval = maxval(abs(model%general%lat))
       global_maxval = parallel_reduce_max(local_maxval)
       if (global_maxval < eps11) then
          call write_log('Failed to read latitude (lat) field from input file', GM_FATAL)
       endif
       local_maxval = maxval(abs(model%general%lon))
       global_maxval = parallel_reduce_max(local_maxval)
       if (global_maxval < eps11) then
          call write_log('Failed to read longitude (lon) field from input file', GM_FATAL)
       endif
    endif

    ! Allocate mask arrays in case they are needed below
    allocate(ice_mask(model%general%ewn, model%general%nsn))
    allocate(floating_mask(model%general%ewn, model%general%nsn))
    allocate(land_mask(model%general%ewn, model%general%nsn))
    allocate(ocean_mask(model%general%ewn, model%general%nsn))

    ! Optionally, compute area scale factors for stereographic map projection.
    ! This should be done after reading the input file, in case the input file contains mapping info.
    ! Note: Not yet enabled for other map projections.
    ! Note: Area factors currently are used only when computing diagnostics for the log file.
    !       If area factors are computed based on projection info, then the ice area and volume
    !         computed in CISM's dycore are corrected for area distortions, giving a better estimate
    !        of the true ice area and volume.
    !       However, applying scale factors will give a mass conservation error (total dmass_dt > 0)
    !        in the diagnostics, because horizontal transport does not account for area factors.
    !        Transport conserves mass only under the assumption of rectangular grid cells.
    ! TODO - Tested only for Greenland (N. Hem.; projection origin offset from N. Pole). Test for other grids.

    if (associated(model%projection%stere)) then

       call glimmap_stere_area_factor(model%projection%stere,  &
                                      model%general%ewn,       &
                                      model%general%nsn,       &
                                      model%numerics%dew*len0, &
                                      model%numerics%dns*len0, &
                                      parallel)
    endif

    ! Write projection info to log
    call glimmap_printproj(model%projection)

    ! Optionally, adjust the input ice thickness in grid cells where there are interior lakes
    !  (usrf - thck > topg), but the ice is above flotation thickness.
    ! In these grid cells, we set thck = usrf - topg, preserving the input usrf and removing the lakes.

    if (model%options%adjust_input_thickness .and. model%options%is_restart == RESTART_FALSE) then
       call glissade_adjust_thickness(model)
    endif

    ! Optionally, smooth the input surface elevation with a Laplacian smoother.
    ! This subroutine does not change the topg, but returns thck consistent with the new usrf.
    ! If the initial usrf is rough, then multiple smoothing passes may be needed to stabilize the flow.

    if (model%options%smooth_input_usrf .and. model%options%is_restart == RESTART_FALSE) then
       call glissade_smooth_usrf(model, nsmooth = 5)
    endif   ! smooth_input_usrf

    ! Optionally, smooth the input topography with a Laplacian smoother.

    if (model%options%smooth_input_topography .and. model%options%is_restart == RESTART_FALSE) then
       call glissade_smooth_topography(model)
    endif   ! smooth_input_topography

    ! Optionally, adjust the input topography in a specified region

    if (model%options%adjust_input_topography .and. model%options%is_restart == RESTART_FALSE) then
       call glissade_adjust_topography(model)
    endif

    ! Optionally, compute the standard deviation of the topography.
    ! Ideally, this would be done when preparing the input file.
    ! But the following code can be used to compute topg_stdev from a high-resolution
    !  input file that contains topg only.
    ! For example, we can read in a 1-km file, set grid_ratio = 8, compute topg_stdev
    !  in 8x8 blocks of cells, and write it to the initial output file.
    ! In a given 8x8 block of data, each cell will contain topg_stdev for that block.
    ! If we then coarsen the output from the 1-km grid to the 8-km grid,
    !  we have a field of topg_stdev to be used in 8-km simulations.

    if (compute_topg_stdev) then

       allocate(topg_global(parallel%global_ewn, parallel%global_nsn))
       allocate(topg_stdev_global(parallel%global_ewn, parallel%global_nsn))

       ! gather topg to global grid
       call distributed_gather_var(model%geometry%topg*thk0, topg_global, parallel)

       if (main_task) then
          call glissade_stdev(parallel%global_ewn,  parallel%global_nsn,  &
                              grid_ratio,                                 &
                              model%numerics%idiag, model%numerics%jdiag, &
                              topg_global,          topg_stdev_global)
       endif

       ! scatter topg_stdev to processors, and rescale
       call distributed_scatter_var(model%geometry%topg_stdev, topg_stdev_global, parallel)
       model%geometry%topg_stdev = model%geometry%topg_stdev/thk0

       deallocate(topg_global)
       deallocate(topg_stdev_global)

    endif

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first

    if (model%options%isostasy == ISOSTASY_COMPUTE) then

       call init_isostasy(model)

    endif

    select case(model%isostasy%whichrelaxed)

    case(RELAXED_TOPO_INPUT)   ! supplied input topography is relaxed

       model%isostasy%relx = model%geometry%topg

    case(RELAXED_TOPO_COMPUTE) ! supplied topography is in equilibrium
                               !TODO - Test the case RELAXED_TOPO_COMPUTE

       call isos_relaxed(model)

    end select

    ! If a 2D bheatflx field is present in the input file, it will have been written
    !  to model%temper%bheatflx.  For the case model%options%gthf = 0, we want to use
    !  a uniform heat flux instead.
    ! If no bheatflx field is present in the input file, then we default to the
    !  prescribed uniform value, model%paramets%geot.

    if (model%options%gthf == GTHF_UNIFORM) then

       ! Check to see if this flux was present in the input file
       ! (by checking whether the flux is nonuniform over the domain)
       if (abs(maxval(model%temper%bheatflx) - minval(model%temper%bheatflx)) > 1.d-6) then  
          call write_log('Setting uniform prescribed geothermal flux')
          call write_log('(Set gthf = 1 to read geothermal flux field from input file)')
       endif

       ! set uniform basal heat flux (positive down)
       model%temper%bheatflx = model%paramets%geot

    endif

    ! Make sure the basal heat flux follows the positive-down sign convention
    if (maxval(model%temper%bheatflx) > 0.0d0) then
       write(message,*) 'Error, Input basal heat flux has positive values: '
       call write_log(trim(message))
       write(message,*) 'this_rank, maxval =', this_rank, maxval(model%temper%bheatflx)
       call write_log(trim(message))
       write(message,*) 'Basal heat flux is defined as positive down, so should be <= 0 on input'
       call write_log(trim(message), GM_FATAL)
    endif

    ! Compute the cell areas of the grid
    model%geometry%cell_area = model%numerics%dew*model%numerics%dns

    ! If running with glaciers, then process the input glacier data
    ! On start-up, this subroutine counts the glaciers.  It should be called before glide_io_createall,
    !  which needs to know nglacier to set up glacier output files with the right dimensions.
    ! On restart, most of the required glacier arrays are in the restart file, and this subroutine
    !  computes a few remaining variable.

    if (model%options%enable_glaciers) then

       !WHL - debug
       ! Glaciers are run with a no-ice BC to allow removal of inactive regions.
       ! This can be problematic when running in a sub-region that has glaciers along the global boundary.
       ! A halo update here for 'thck' will remove ice from cells along the global boundary.
       ! It is best to do this before initializing glaciers, so that ice that initially exists
       !  in these cells is removed before computing the area and thickness targets.
       !TODO - These calls are repeated a few lines below.  Try moving them up, before the call
       !       to glissade_glacier_init.  I don't think it's possible to move the glissade_glacier_init call
       !       down, because we need to compute nglacier before setting up output files.

       call parallel_halo(model%geometry%thck, parallel)
       ! calculate the lower and upper ice surface
       call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
       model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

       call glissade_glacier_init(model, model%glacier)

    endif

    ! open all output files
    call openall_out(model)

    ! create glide I/O variables
    call glide_io_createall(model, model)

    ! initialize glissade components

    ! Set some variables in halo cells
    ! Note: We need thck and artm in halo cells so that temperature will be initialized correctly
    !        (if not read from input file).
    !       We do an update here for temp in case temp is read from an input file.
    !       If temp is computed below in glissade_init_therm (based on the value of options%temp_init),
    !        then the halos will receive the correct values.

    call parallel_halo(model%geometry%thck, parallel)
    call parallel_halo(model%climate%artm,  parallel)
    call parallel_halo(model%temper%temp,   parallel)
    call parallel_halo(model%temper%tempunstag, parallel)

    ! calculate the lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    ! Note: For outflow BCs, most fields (thck, usrf, temp, etc.) are set to zero in the global halo,
    !        to create ice-free conditions. However, we might not want to set topg = 0 in the global halo,
    !        because then the global halo will be interpreted as ice-free land, whereas we may prefer to
    !        treat it as ice-free ocean. For this reason, topg is extrapolated from adjacent cells.
    !       Similarly, for no_ice BCs, we want to zero out ice state variables adjacent to the global boundary,
    !        but we do not want to zero out the topography.
    ! Note: For periodic BCs, there is an optional argument periodic_offset_ew for topg.
    !       This is for ismip-hom experiments. A positive EW offset means that
    !        the topography in west halo cells will be raised, and the topography
    !        in east halo cells will be lowered.  This ensures that the topography
    !        and upper surface elevation are continuous between halo cells
    !        and locally owned cells at the edge of the global domain.
    !       In other cases (anything but ismip-hom), periodic_offset_ew = periodic_offset_ns = 0,
    !        and this argument will have no effect.

    if (model%general%global_bc == GLOBAL_BC_OUTFLOW .or.  &
        model%general%global_bc == GLOBAL_BC_NO_ICE) then
       call parallel_halo_extrapolate(model%geometry%topg, parallel)
    else  ! other global BCs, including periodic
       call parallel_halo(model%geometry%topg, parallel, &
                          periodic_offset_ew = model%numerics%periodic_offset_ew)
    endif

    if (model%options%whichtemp == TEMP_ENTHALPY) &
         call parallel_halo(model%temper%waterfrac, parallel)

    ! halo update for kinbcmask (= 1 where uvel and vvel are prescribed, elsewhere = 0)
    ! Note: Instead of assuming that kinbcmask is periodic, we extrapolate it into the global halo
    !       (and also into the north and east rows of the global domain, which are not included 
    !       on the global staggered grid).
    call staggered_parallel_halo_extrapolate (model%velocity%kinbcmask, parallel)  ! = 1 for Dirichlet BCs

    !TODO - Remove call to init_velo in glissade_initialise?
    !       Most of what's done in init_velo is needed for SIA only, but still need velowk for call to wvelintg
    call init_velo(model)

    ! Initialize basal hydrology, if needed
    call glissade_basal_water_init(model)

    ! Initialize acab, if SMB (with different units) was read in
    if (model%options%smb_input == SMB_INPUT_MMYR_WE) then
       ! Convert units from mm/yr w.e. to m/yr ice, then convert to model units
       model%climate%acab(:,:) = (model%climate%smb(:,:) * (rhow/rhoi)/1000.d0) / scale_acab
    endif

    ! Initialize artm_corrected.  This is equal to artm, plus any prescribed temperature anomaly.
    model%climate%artm_corrected(:,:) = model%climate%artm(:,:)

    if (model%options%enable_artm_anomaly) then
       call glissade_add_2d_anomaly(model%climate%artm_corrected,          &   ! degC
                                    model%climate%artm_anomaly,            &   ! degC
                                    model%climate%artm_anomaly_timescale,  &   ! yr
                                    model%numerics%time)                       ! yr
    endif

    ! Initialize the temperature profile in each column
    call glissade_init_therm(model%options%temp_init,    model%options%is_restart,  &
                             model%general%ewn,          model%general%nsn,         &
                             model%general%upn,                                     &
                             model%numerics%idiag_local, model%numerics%jdiag_local,&
                             model%numerics%rdiag_local,                            &
                             model%numerics%sigma,       model%numerics%stagsigma,  &
                             model%numerics%dups,                                   &
                             model%geometry%thck*thk0,                              & ! m
                             model%climate%artm_corrected,                          & ! deg C
                             model%climate%acab*thk0/tim0,                          & ! m/s
                             model%temper%bheatflx,                                 & ! W/m^2, positive down
                             model%temper%pmp_offset,                               & ! deg C
                             model%temper%temp,                                     & ! deg C
                             model%temper%tempunstag)                                 ! deg C

    if (model%options%gthf == GTHF_COMPUTE) then
       call not_parallel(__FILE__,__LINE__)
       call init_lithot(model)
    end if

    ! Dycore-specific velocity solver initialization
    select case (model%options%whichdycore)

    case ( DYCORE_GLISSADE )   ! glissade finite-element

       call glissade_velo_higher_init

    case ( DYCORE_ALBANYFELIX)

       call felix_velo_init(model)

    end select

    !TODO - Add halo updates of state variables here?

    ! WHL: Set make_ice_domain_mask = .true. above to create an ice_domain_mask at initialization.
    !      False by default, because usually we will read in a mask created previously.
    !      Could make this a config option if desired.

    ! Create an ice_domain_mask which includes all cells with thck > 0 and/or topg > 0,
    !  along with one or more buffer layers of ocean cells.
    ! This could be made more complex, for instance by choosing a negative threshold for topg.
    ! For now, we just want to create a mask that can be used to identify active blocks
    !  for whole-ice-sheet experiments.

    if (make_ice_domain_mask) then

       where (model%geometry%thck > 0.0d0 .or. model%geometry%topg > 0.0d0)
!!       where (model%geometry%thck > 0.0d0)  ! uncomment for terrestrial margins
          model%general%ice_domain_mask = 1
       elsewhere
          model%general%ice_domain_mask = 0
       endwhere

       ! Extend the mask a few cells in each direction to be on the safe side.
       ! The number of buffer layers could be made a config parameter.

       allocate(ice_domain_mask(model%general%ewn,model%general%nsn))

       do k = 1, 3
          call parallel_halo(model%general%ice_domain_mask, parallel)
          ice_domain_mask = model%general%ice_domain_mask   ! temporary copy
          do j = nhalo+1, model%general%nsn - nhalo
             do i = nhalo+1, model%general%ewn - nhalo
                if (ice_domain_mask(i-1,j) == 1 .or. ice_domain_mask(i+1,j) == 1 .or. &
                     ice_domain_mask(i,j-1) == 1 .or. ice_domain_mask(i,j+1) == 1) then
                   model%general%ice_domain_mask(i,j) = 1
                endif
             enddo
          enddo
       enddo

       deallocate(ice_domain_mask)

    endif   ! make_ice_domain_mask

    ! If unstagbeta (i.e., beta on the scalar ice grid) was read from an input file,
    !  then interpolate it to beta on the staggered grid.
    ! NOTE: unstagbeta is initialized to unphys_val (a large negative number),
    !       so its maxval will be > 0 only if the field is read in.
    ! We can make an exception for ISHOM case C; for greater accuracy we set beta in 
    !  subroutine calcbeta instead of interpolating from unstagbeta (one processor only).

    if (maxval(model%velocity%unstagbeta) > 0.d0 .and.   &
               model%options%which_ho_babc /= HO_BABC_ISHOMC) then  ! interpolate to staggered grid
       call write_log('Interpolating beta from unstaggered (unstagbeta) to staggered grid (beta)')
       if (maxval(model%velocity%beta) > 0.0d0 ) then
          call write_log('Warning: the input "beta" field will be overwritten with values interpolated from the input &
               &"unstagbeta" field!')
       endif

       ! do a halo update and interpolate to the staggered grid
       ! Note: stagger_margin_in = 0 => use all values in the staggering, including where ice is absent
       call parallel_halo(model%velocity%unstagbeta,  parallel)
       call glissade_stagger(model%general%ewn,          model%general%nsn,      &
                             model%velocity%unstagbeta,  model%velocity%beta,    &
                             stagger_margin_in = 0)

    endif  ! unstagbeta > 0

    ! Note: If the beta field has been read from an external file, then model%velocity%beta should not be modified.
    !       If beta is weighted by f_ground or otherwise modified, the modified field is model%velocity%beta_internal.

    ! The MISMIP 3D test case requires reading in a spatial factor that multiplies Coulomb_C.
    ! This factor is read in on the unstaggered grid, then interpolated to the staggered grid.
    ! If this factor is not present in the input file, then set it to 1 everywhere.
    if (model%options%use_c_space_factor) then

       if (maxval(model%basal_physics%c_space_factor) > tiny(0.d0)) then  ! c_space_factor was read in

          ! do a halo update and interpolate to the staggered grid
          ! Note: stagger_margin_in = 0 => use all values in the staggering, including where ice is absent
          call parallel_halo(model%basal_physics%c_space_factor, parallel)
          call glissade_stagger(model%general%ewn,                  model%general%nsn,                       &
                                model%basal_physics%c_space_factor, model%basal_physics%c_space_factor_stag, &
                                stagger_margin_in = 0)

       else  ! c_space_factor was not read in; set to 1 everywhere, so it will be ignored when computing beta
             ! With this value, it will be ignored when computing beta

          ! Note: It would be possible here to set c_space_factor to values different from 1,
          !       if not reading it from the input file.
          model%basal_physics%c_space_factor(:,:) = 1.0d0
          model%basal_physics%c_space_factor_stag(:,:) = 1.0d0

       endif  ! maxval(c_space_factor) > 0

    else   ! use_c_space_factor = F

       ! set equal to 1 everywhere, so it will be ignored when computing beta
       model%basal_physics%c_space_factor(:,:) = 1.0d0
       model%basal_physics%c_space_factor_stag(:,:) = 1.0d0

    endif  ! use_c_space_factor

    ! Note: The basal process option is currently disabled.
    ! initialize basal process module
!!    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bmod == BAS_PROC_FASTCALC) then        
!!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!!                              model%numerics%dttem)
!!    end if      

    ! If acab is to be overwritten for some cells, then set overwrite_acab_mask = 1 for these cells.
    ! We can overwrite the input acab with a fixed value (typically negative) where
    ! (1) the input acab = 0 at initialization, or
    ! (2) the input thck <= overwrite_acab_minthck at initialization
    ! Note: This option is designed for standalone runs, and should be used only with caution for coupled runs.
    !       On restart, overwrite_acab_mask is read from the restart file.

    if (model%climate%overwrite_acab_value /= 0 .and. model%options%is_restart == RESTART_FALSE) then

       call glissade_overwrite_acab_mask(model%options%overwrite_acab,          &
                                         model%climate%acab,                    &
                                         model%geometry%thck,                   &
                                         model%climate%overwrite_acab_minthck,  &
                                         model%climate%overwrite_acab_mask)

    endif

    ! calculate mask
    ! Note: This call includes a halo update for thkmask
    !TODO - Replace with glissade masks?
    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    ! compute halo for relaxed topography
    ! Note: See comments above with regard to the halo update for topg.
    !       Outflow BCs zero out scalars beyond the global boundary, and no_ice BCs zero out scalars
    !        adjacent to or beyond the global boundary. This is an appropriate treatment for
    !        ice state variables, but not for bed topography and related fields (like relx).
    !TODO - Is this halo update necessary?
    if (model%general%global_bc == GLOBAL_BC_OUTFLOW .or. &
        model%general%global_bc == GLOBAL_BC_NO_ICE) then
       call parallel_halo_extrapolate(model%isostasy%relx, parallel)
    else
       call parallel_halo(model%isostasy%relx, parallel)
    endif

    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)

    ! optional unit tests

    if (test_halo) then
       call glissade_test_halo (model)
       call parallel_finalise
    endif
     
    if (test_transport) then
       where (model%geometry%thck > model%numerics%thklim)
          model%geometry%thck = thk_init/thk0
       elsewhere
          model%geometry%thck = 0.d0
       endwhere
    endif

    ! Compute a mask to identify cells with a marine path to ice-free ocean.
    ! Currently used for some basal melting options, to identify cells where bmlt_float can be nonzero.
    ! Note: This mask needs to be recomputed at runtime if topg or eus changes.
    !       It is currently updated at the end of glissade_isostasy_solve.
    !       Although thck is passed in, it is used only to identify ocean cells that can seed the fill.
    !       The mask should not depend on ice thickness, unless there are land-locked regions with
    !        deep topography that are alternately ice-covered and ice-free.

    call glissade_marine_connection_mask(&
         model%general%ewn,          model%general%nsn,          &
         parallel,                                               &
         model%numerics%idiag_local, model%numerics%jdiag_local, &
         model%numerics%rdiag_local,                             &
         model%geometry%thck * thk0, model%geometry%topg * thk0, &
         model%climate%eus * thk0,   0.0d0,                      &  ! thklim = 0
         model%geometry%marine_connection_mask)

    ! TODO: Move calving-related initialization to a separate subroutine.

    ! initial calving, if desired
    ! Note: Do initial calving only for a cold start with evolving ice, not for a restart
    if (l_evolve_ice .and. &
         model%options%calving_init == CALVING_INIT_ON .and. &
         model%options%is_restart == RESTART_FALSE) then

       ! ------------------------------------------------------------------------
       ! Note: The initial calving solve is treated differently from the runtime calving solve.
       !       In particular, calving-front culling is done at initialization only.
       !       Culling may also be used to remove a row of thin cells (~1 m)
       !        at the calving front, as present in some input data sets.
       !       But we do not want to remove calving_front cells every timestep.
       ! ------------------------------------------------------------------------

       call glissade_calving_solve(model, .true.)   ! init_calving = .true.

       ! The mask needs to be recalculated after calving.
       ! Note: glide_set_mask includes a halo update for thkmask.
       !TODO - Delete this call? Glissade dycore does not use thkmask.
       call glide_set_mask(model%numerics,                                &
                           model%geometry%thck,  model%geometry%topg,     &
                           model%general%ewn,    model%general%nsn,       &
                           model%climate%eus,    model%geometry%thkmask,  &
                           model%geometry%iarea, model%geometry%ivol)

    endif  ! initial calving

    ! Initialize the effective pressure calculation

    if (model%options%is_restart == RESTART_FALSE) then
       call glissade_init_effective_pressure(model%options%which_ho_effecpress,  &
                                             model%basal_physics)
    endif

    ! Initialize powerlaw_c and coulomb_c.
    ! Note: This can set powerlaw_c and coulomb_c to nonzero values when they are never used,
    !       but is simpler than checking all possible basal friction options.

    if (model%options%is_restart == RESTART_FALSE) then
       if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_CONSTANT) then
          model%basal_physics%powerlaw_c = model%basal_physics%powerlaw_c_const
       endif
       if (model%options%which_ho_coulomb_c == HO_COULOMB_C_CONSTANT) then
          model%basal_physics%coulomb_c = model%basal_physics%coulomb_c_const
       endif
    endif

    ! Optionally, do initial calculations for inversion
    ! At the start of the run (but not on restart), this might lead to further thickness adjustments,
    !  so it should be called before computing the calving mask.
    !TODO: Separate the basal friction inversion from the bmlt_basin inversion.

    if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or.  &
        model%options%which_ho_coulomb_c  == HO_COULOMB_C_INVERSION  .or.  &
        model%options%which_ho_bmlt_basin_inversion == HO_BMLT_BASIN_INVERSION_COMPUTE) then

       call glissade_init_inversion(model)

    endif  ! inversion for Cp, Cc or bmlt

    ! If using a mask to force ice retreat, then set the reference thickness (if not already read in).

    if (model%options%force_retreat /= FORCE_RETREAT_NONE) then

       ! Set reference_thck
       ! This field is loaded if present in the input file.  Otherwise, reference_thck is set to the initial thickness.
       local_maxval = maxval(model%geometry%reference_thck)
       global_maxval = parallel_reduce_max(local_maxval)
       if (global_maxval < eps11) then
          write(message,*) 'Setting reference_thck to the initial ice thickness'
          call write_log(trim(message))
          model%geometry%reference_thck = model%geometry%thck
       else
          write(message,*) 'reference_thck was read from the input/restart file'
          call write_log(trim(message))
       endif

       ! Check whether a nonzero retreat mask has been read in. If not, then write a message.
       ! Note: This is not necessarily an error.  If reading the retreat mask from a forcing file,
       !       the first nonzero mask may not be read until the model time is greater than tstart.
       local_maxval = maxval(model%geometry%ice_fraction_retreat_mask)
       global_maxval = parallel_reduce_max(local_maxval)
       if (global_maxval < eps11) then
          call write_log('Initial ice_fraction_retreat_mask = 0 everywhere')
       endif

       !WHL - debug
       if (this_rank == rtest) then
          i = itest
          j = jtest
          print*, 'force_retreat option =', model%options%force_retreat
          print*, ' '
          print*, 'reference_thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%reference_thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'ice_fraction_retreat_mask:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%ice_fraction_retreat_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    endif   ! force_retreat

    !TODO - Move to an init_calving subroutine
    if ( (model%options%whichcalving == CALVING_THCK_THRESHOLD  .or.  &
          model%options%whichcalving == EIGENCALVING  .or. model%options%whichcalving == CALVING_DAMAGE) ) then

       ! Check whether thck_calving_threshold was already read in.  If not, then compute it.
       ! Note: Do not use the restart_false logic, since it may be convenient to compute this field on restart.

       local_maxval = maxval(model%calving%thck_calving_threshold)
       global_maxval = parallel_reduce_max(local_maxval)

       write(message,*) 'thck_calving_threshold, global_maxval =', global_maxval
       call write_log(trim(message))

       if (global_maxval < eps11) then
          write(message,*) 'Set thck_calving_threshold at startup'
          call write_log(trim(message))

          ! Given the initial calving front thickness, initialize thck_calving_threshold.
          ! If calving_minthck > 0 in the config file, set thck_calving_threshold = calving_minthck everywhere.
          ! Otherwise, set thck_calving_threshold to the initial calving_front thickness, and extrapolate upstream and downstream.
          ! During the run, any ice in the calving domain with thck < thck_calving_threshold will be removed.
          ! On restart, thck_calving_threshold will be read from the restart file.
          ! Note: This is done after initial calving, which may include iceberg removal or calving-front culling.
          !       Since thck_calving_threshold = 0 initially, no ice is thinner than the threshold during the
          !        initial call above to glissade_calve_ice.  The threshold takes effect during the first timestep.

          call glissade_thck_calving_threshold_init(&
               model%general%ewn,         model%general%nsn,            &
               parallel,                                                &
               itest,   jtest,    rtest,                                &
               model%options%which_ho_calving_front,                    &
               model%geometry%thck*thk0,  model%geometry%topg*thk0,     &
               model%climate%eus*thk0,    model%numerics%thklim*thk0,   &
               model%geometry%marine_connection_mask,                   &
               model%calving%minthck,                                   &
               model%calving%thck_calving_threshold)
       endif

    endif  ! thickness-based calving

    if ((model%options%whichcalving == CALVING_GRID_MASK .or. model%options%apply_calving_mask)  &
         .and. model%options%is_restart == RESTART_FALSE) then

       ! Initialize the no-advance calving_mask
       ! Note: This is done after initial calving, which may include iceberg removal or calving-front culling.
       !       The calving front that exists after initial culling is the one that is held fixed during the simulation.
       ! Note: Typically, the calving mask is set to 1 (i.e., force calving) in all ice-free ocean cells.
       !       If usfc_obs and vsfc_obs have been read in, then the mask will be set to 0 in ice-free ocean cells
       !        where the observed velocity is nonzero.  Ice-free cells can have nonzero velocity
       !        if the input velocity comes from a different data source than the input thickness.
       ! Note: calving_front_x and calving_front_y already have units of m, so do not require multiplying by len0.
       ! On restart, calving_mask is read from the restart file.

       call glissade_calving_mask_init(&
            model%numerics%dew*len0,           model%numerics%dns*len0,           &
            parallel,                                                             &
            model%geometry%thck*thk0,          model%geometry%topg*thk0,          &  ! m
            model%climate%eus*thk0,            model%numerics%thklim*thk0,        &  ! m
            model%velocity%usfc_obs*vel0*scyr, model%velocity%vsfc_obs*vel0*scyr, &  ! m/yr
            model%calving%calving_front_x,     model%calving%calving_front_y,     &
            model%calving%calving_mask)

    endif   ! calving grid mask

    ! Note: The DIVA solver needs a halo update for effective viscosity.
    !       This is done at the end of glissade_diagnostic_variable_solve, which in most cases is sufficient.
    !       However, if we are (1) reading efvs from an input file and (2) solving for velocity before
    !        the halo update, then we need a halo update here too, to avoid symmetry issues in the velocity solver.
    !       I ran into this issue when running MISMIP+, which does cold starts (restart = 0) from files containing efvs.
    !       An update is done here regardless of code options, just to be on the safe side.
    call parallel_halo(model%stress%efvs, parallel)

    ! recalculate the lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    ! save the initial ice thickness
    model%geometry%thck_old(:,:) = model%geometry%thck(:,:)

    ! initialize ocean forcing data, if desired
    ! Currently, this is done only when using the ISMIP6 basal melting parameterization
    ! Note: Need the current value of lsrf when calling this subroutine

    if (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then

       ! update some masks
       !TODO: Move these mask updates to the thermal_forcing_init subroutine?
       !TODO: Modify glissade_get_masks so that 'parallel' is not needed
       call glissade_get_masks(model%general%ewn, model%general%nsn,       &
                               parallel,                                   &
                               model%geometry%thck, model%geometry%topg,   &
                               model%climate%eus,   0.0d0,                 &  ! thklim = 0
                               ice_mask,                                   &
                               floating_mask = floating_mask,              &
                               land_mask = land_mask)

       ! update the grounded fraction, f_ground_cell
       call glissade_grounded_fraction(model%general%ewn,             &
                                       model%general%nsn,             &
                                       parallel,                      &
                                       itest, jtest, rtest,           &  ! diagnostic only
                                       model%geometry%thck*thk0,      &
                                       model%geometry%topg*thk0,      &
                                       model%climate%eus*thk0,        &
                                       ice_mask,                      &
                                       floating_mask,                 &
                                       land_mask,                     &
                                       model%options%which_ho_ground, &
                                       model%options%which_ho_flotation_function, &
                                       model%options%which_ho_fground_no_glp,     &
                                       model%geometry%f_flotation,    &
                                       model%geometry%f_ground,       &
                                       model%geometry%f_ground_cell,  &
                                       model%geometry%topg_stdev*thk0)

       call glissade_bmlt_float_thermal_forcing_init(model, model%ocean_data)

    endif

    ! clean up
    deallocate(ice_mask)
    deallocate(floating_mask)
    deallocate(land_mask)
    deallocate(ocean_mask)

    if (main_task) print*, 'Done in glissade_initialise'

  end subroutine glissade_initialise
  
!=======================================================================

  subroutine glissade_tstep(model, time)

    ! Perform time-step of an ice model instance with the Glissade dycore

    use cism_parallel, only:  parallel_type, not_parallel

    use glimmer_paramets, only: tim0, len0, thk0
    use glimmer_physcon, only: scyr
    use glide_mask, only: glide_set_mask, calc_iareaf_iareag

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    real(dp), intent(in) :: time         ! current time in years

    ! --- Local variables ---

    integer :: i, j
    integer :: itest, jtest, rtest

    ! ========================

    ! Set debug diagnostics
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    ! Update internal clock
    model%numerics%time = time  
    model%numerics%tstep_count = model%numerics%tstep_count + 1

    if (main_task .and. verbose_glissade) print*, 'glissade_tstep, test_count =', model%numerics%tstep_count

    ! optional transport test
    ! code execution will end when this is done
    if (test_transport) then
       call glissade_test_transport (model)
       return
    endif

    ! save the old ice thickness; used for diagnostics and tendencies
    ! also used to reset thickness for the no-evolution option
    model%geometry%thck_old(:,:) = model%geometry%thck(:,:)

    ! ------------------------------------------------------------------------
    ! Calculate isostatic adjustment
    ! ------------------------------------------------------------------------
    !
    ! Note: This call used to be near the end of the glissade time step, between
    !       calving and the velocity solve. But this can be problematic, because
    !       a cell identified as grounded for calving purposes can become floating
    !       as a result of isostatic adjustment, or vice versa.
    !       It is better to compute isostasy just after the velocity solve,
    !       at the start of the next time step.
    !
    ! Matt Hoffman writes:
    ! Is this isostasy call in the right place?
    ! Consider for a forward Euler time step:
    ! With a relaxing mantle model, topg is a prognostic (time-evolving) variable:
    !      topg1 = f(topg0, thk0, ...)
    ! However, for a fluid mantle where the adjustment is instantaneous, topg is a diagnostic variable
    !(comparable to calculating floatation height of ice in the ocean):
    !      topg1 = f(thk1)
    ! In either case, the topg update should be separate from the thickness evolution (because thk1 = f(thk0, vel0=g(topg0,...)).
    ! However, if the isostasy calculation needs topg0, the icewaterload call should be made BEFORE thck is updated.
    ! If the isostasy calculation needs topg1, the icewaterload call should be made AFTER thck is updated.
    ! Also, we should think about when marinlim, usrf, lsrf, derivatives should be calculated relative to the topg update via isostasy.
    !
    ! WHL writes (May 2017):
    ! When isostasy is turned on, it is usually run with a relaxing mantle.
    ! With the call moved to the start of the time step, both the icewaterload call (if needed) and
    !  the relaxation are done before the ice thickness update. So we have
    !       topg1 = f(topg0, thk0, ...)
    !  followed by
    !       thk1  = f(thk0, vel0=g(topg0,...)
    ! I think this is what is desired.
    ! ------------------------------------------------------------------------

    call glissade_isostasy_solve(model)

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    !TODO Not sure if this is in the right place.  G1=f(G0,T0) and T1=g(G0,T0)  
    !     If we update G1 now, then we will be doing T1=g(G1,T0).
    if (model%options%gthf == GTHF_COMPUTE) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------
    ! Do the vertical thermal solve if it is time to do so.
    ! Vertical diffusion and strain heating only; no temperature advection.
    ! Note: model%numerics%tinc and model%numerics%time have units of years.
    !       dttem has scaled units, so multiply by tim0/scyr to convert to years.
    ! ------------------------------------------------------------------------ 

    if ( model%numerics%tinc > mod(model%numerics%time, model%numerics%dttem*tim0/scyr)) then

       if (model%options%which_ho_thermal_timestep == HO_THERMAL_BEFORE_TRANSPORT) then

          ! vertical thermal solve before transport
          call glissade_thermal_solve(model,  &
                                      model%numerics%dttem*tim0)   ! convert dt from model units to s

       elseif (model%options%which_ho_thermal_timestep == HO_THERMAL_SPLIT_TIMESTEP) then

          ! vertical thermal solve split into two parts, before and after transport
          call glissade_thermal_solve(model,  &
                                      model%numerics%dttem*tim0/2.0d0)

       endif

    end if

    ! ------------------------------------------------------------------------ 
    ! Compute the basal melt rate beneath floating ice.
    ! (The basal melt rate beneath grounded ice is part of the thermal solve.)
    ! ------------------------------------------------------------------------ 

    call glissade_bmlt_float_solve(model)

    ! Add bmlt_float to bmlt_ground to determine the total bmlt

    ! Note: bmlt = bmlt_ground + bmlt_float may not be equal to the applied melt rate in a given cell,
    !       if ice is thin or absent in the cell.
    ! Note: bmlt does not include bmlt_float_inversion, which is a separate input/output field.
    ! Note: bmlt_ground is computed in glissade_therm_driver.
    !       If glissade_therm_driver is called twice per time step, then bmlt_ground from
    !        the second time is not included in the transport solve, but is diagnostic only.
    !       That is, the transport scheme assumes that the bmlt_ground rate computed during the
    !        first call is applied during the entire time step.
    !       This might lead to small violations of energy conservation.
    !       TODO: Separate the bmlt_ground computation from the temperature computation?

    model%basal_melt%bmlt(:,:) = model%basal_melt%bmlt_ground(:,:) + model%basal_melt%bmlt_float(:,:)

    if (verbose_glissade .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Before thickness solver, thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bmlt_ground (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%basal_melt%bmlt_ground(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bmlt_float (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%basal_melt%bmlt_float(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

    ! ------------------------------------------------------------------------ 
    ! Calculate ice thickness and tracer evolution under horizontal transport.
    ! The surface and basal mass balances are also applied here.
    ! ------------------------------------------------------------------------ 

    call glissade_thickness_tracer_solve(model)

    ! ------------------------------------------------------------------------ 
    ! Calculate iceberg calving
    ! ------------------------------------------------------------------------ 

    call glissade_calving_solve(model, .false.)   ! init_calving = .false.

    if (verbose_glissade .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'After calving, thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') model%geometry%thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

    ! ------------------------------------------------------------------------
    ! Clean up variables in ice-free columns.
    ! This subroutine should be called after transport and calving, which may
    !  have set thck = 0 in some cells without zeroing out basal water and tracers.
    ! ------------------------------------------------------------------------ 

    call glissade_cleanup_icefree_cells(model)

    ! ------------------------------------------------------------------------ 
    ! Increment the ice age.
    ! If a cell becomes ice-free, the age is reset to zero.
    ! Note: Internally, the age has the same units as dt, but on output it will be converted to years.
    ! ------------------------------------------------------------------------ 
    
    if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
       do j = 1, model%general%nsn 
          do i = 1, model%general%ewn 
             if (model%geometry%thck(i,j) > 0.0d0) then
                model%geometry%ice_age(:,i,j) = model%geometry%ice_age(:,i,j) + model%numerics%dt
             else
                model%geometry%ice_age(:,i,j) = 0.0d0
             endif
          enddo
       enddo
    endif

    ! glissade_calve_ice adjusts thickness for calved ice.  Therefore the mask needs to be recalculated.
    ! Note: glide_set_mask includes a halo update of thkmask
    ! This time we want to calculate the optional arguments iarea and ivol because thickness 
    ! will not change further during this time step.
    !TODO - Remove this call to glide_set_mask?
    !       This subroutine is called at the beginning of glissade_velo_driver,
    !        so a call here is not needed for the velo diagnostic solve.
    !       The question is whether it is needed for the isostasy.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    ! --- Calculate global area of ice that is floating and grounded.
    !TODO  May want to calculate iareaf and iareag in glide_write_diag and remove those calculations here.  

    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%thkmask,                        &
                            model%geometry%iareaf, model%geometry%iareag)

    ! ------------------------------------------------------------------------
    ! Do the vertical thermal solve if it is time to do so.
    ! Note: A thermal solve should be done here (using option HO_THERMAL_AFTER_TRANSPORT 
    !       or HO_THERMAL_SPLIT_TIMESTEP) if it is desired to update the bed temperature 
    !       and pmp temperature after transport and before the velocity solve.
    ! ------------------------------------------------------------------------

    if ( model%numerics%tinc > mod(model%numerics%time, model%numerics%dttem*tim0/scyr)) then

       if (model%options%which_ho_thermal_timestep == HO_THERMAL_AFTER_TRANSPORT) then

          ! vertical thermal solve after transport
          call glissade_thermal_solve(model,  &
                                      model%numerics%dttem*tim0)   ! convert dt from model units to s

       elseif (model%options%which_ho_thermal_timestep == HO_THERMAL_SPLIT_TIMESTEP) then

          ! vertical thermal solve split into two parts, before and after transport
          call glissade_thermal_solve(model,  &
                                      model%numerics%dttem*tim0/2.0d0)

       endif

    end if  ! take a temperature time step

    ! ------------------------------------------------------------------------
    ! Calculate diagnostic variables, including ice velocity
    ! ------------------------------------------------------------------------

    call glissade_diagnostic_variable_solve(model)

    !TODO - Any halo updates needed at the end of glissade_tstep?

  end subroutine glissade_tstep

!=======================================================================

  subroutine glissade_bmlt_float_solve(model)

    ! Solve for basal melting beneath floating ice.

    use glimmer_paramets, only: eps08, eps11, tim0, thk0, len0
    use glimmer_physcon, only: scyr
    use glissade_bmlt_float, only: glissade_basal_melting_float, &
         glissade_bmlt_float_thermal_forcing, verbose_bmlt_float
    use glissade_transport, only: glissade_add_2d_anomaly
    use glissade_masks, only: glissade_get_masks
    use cism_parallel, only:  parallel_reduce_max

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

    integer, dimension(model%general%ewn, model%general%nsn) ::   &
         ice_mask,              & ! = 1 if ice is present (thck > 0, else = 0
         floating_mask,         & ! = 1 if ice is present (thck > 0) and floating, else = 0
         ocean_mask,            & ! = 1 if topg is below sea level and ice is absent, else = 0
         land_mask                ! = 1 if topg - eus >= 0

    real(dp), dimension(model%general%ewn, model%general%nsn) ::   &
         h_cavity                 ! ocean cavity thickness, >= 0 (m)

    ! melt rate field for ISMIP6
    real(dp), dimension(model%general%ewn, model%general%nsn) ::   &
         bmlt_float_transient     ! basal melt rate for ISMIP6 thermal forcing (m/s)

    real(dp) :: previous_time     ! time (yr) at the end of the previous timestep
    real(dp) :: time_from_start   ! time (yr) since the start of applying the anomaly
    real(dp) :: anomaly_fraction  ! fraction of full anomaly to apply
    real(dp) :: tf_anomaly        ! uniform thermal forcing anomaly (deg C), applied everywhere
    integer  :: tf_anomaly_basin  ! basin number where anomaly is applied;
                                  ! for default value of 0, apply to all basins

    real(dp) :: local_maxval, global_maxval   ! max values of a given variable
    integer :: i, j
    integer :: ewn, nsn
    real(dp) :: dew, dns
    integer :: itest, jtest, rtest

    type(parallel_type) :: parallel   ! info for parallel communication

    ! set grid dimensions
    ewn = model%general%ewn
    nsn = model%general%nsn

    dew = model%numerics%dew
    dns = model%numerics%dns

    ! set debug diagnostics
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    parallel = model%parallel

    ! ------------------------------------------------------------------------
    ! Compute the basal melt rate beneath floating ice.
    ! Note: model%basal_melt is a derived type with various fields and parameters
    ! ------------------------------------------------------------------------

    !WHL - Put other simple options in this subroutine instead of glissade_basal_melting_float?

    if (main_task .and. verbose_glissade) print*, 'Call glissade_bmlt_float_solve'

    ! Compute masks:
    ! Note: The '0.0d0' argument is thklim. Any ice with thck > 0 gets ice_mask = 1.

    !TODO: Modify glissade_get_masks so that 'parallel' is not needed
    call glissade_get_masks(ewn,                 nsn,                   &
                            parallel,                                   &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   0.0d0,                 &  ! thklim = 0
                            ice_mask,                                   &
                            floating_mask = floating_mask,              &
                            ocean_mask = ocean_mask,                    &
                            land_mask = land_mask)

    ! Compute bmlt_float depending on the whichbmlt_float option

    if (model%options%whichbmlt_float == BMLT_FLOAT_NONE) then

       model%basal_melt%bmlt_float(:,:) = 0.0d0

    elseif (model%options%whichbmlt_float == BMLT_FLOAT_EXTERNAL) then

       ! Apply the external melt rate

       model%basal_melt%bmlt_float(:,:) = model%basal_melt%bmlt_float_external(:,:)

       ! Optionally, multiply bmlt_float by a scalar adjustment factor
       if (model%basal_melt%bmlt_float_factor /= 1.0d0) then
          model%basal_melt%bmlt_float(:,:) = model%basal_melt%bmlt_float(:,:) * model%basal_melt%bmlt_float_factor
       endif

    elseif (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then

       if (this_rank == rtest .and. verbose_bmlt_float) then
          print*, ' '
          print*, 'Compute bmlt_float at runtime from current thermal forcing'
       endif

       !Note: Currently, there is no difference between ocean_data_domain = 0
       !       (compute internally) and ocean_data_domain = 1 (read from file).
       !      Thermal forcing is initialized to zero and then is loaded from
       !       the input or forcing file, if present.
       !      If ocean_data_domain = 2, then the thermal forcing is set by Glad;
       !       any values read from an input or forcing file are overwritten.
       !      CISM is not yet able to compute thermal forcing internally.
       !TODO: Add code to compute thermal forcing internally.

       ! Check for positive values of thermal forcing.
       ! If whichbmlt_float = BMLT_FLOAT_THERMAL_FORCING, but there are no positive values,
       !  something is probably wrong.

       local_maxval = maxval(model%ocean_data%thermal_forcing)
       global_maxval = parallel_reduce_max(local_maxval)
       if (global_maxval <= eps11) then
          call write_log('thermal forcing <= 0 everywhere, GM_WARNING')
       endif

       !-----------------------------------------------
       ! Optionally, apply a uniform thermal forcing anomaly everywhere.
       ! This anomaly can be phased in linearly over a prescribed timescale.
       !-----------------------------------------------

       if (model%ocean_data%thermal_forcing_anomaly /= 0.0d0) then
          time_from_start = model%numerics%time - model%ocean_data%thermal_forcing_anomaly_tstart
          if (time_from_start + eps08 > model%ocean_data%thermal_forcing_anomaly_timescale .or.  &
               model%ocean_data%thermal_forcing_anomaly_timescale == 0.0d0) then
             anomaly_fraction = 1.0d0   ! apply the full anomaly
          else
             anomaly_fraction = floor(time_from_start + eps08) &
                  / model%ocean_data%thermal_forcing_anomaly_timescale
          endif
          tf_anomaly = anomaly_fraction * model%ocean_data%thermal_forcing_anomaly
          tf_anomaly_basin = model%ocean_data%thermal_forcing_anomaly_basin
          if (this_rank == rtest .and. verbose_bmlt_float) then
             print*, 'time_from_start (yr):', time_from_start
             print*, 'ocean_data%thermal forcing anomaly  (deg):', model%ocean_data%thermal_forcing_anomaly
             print*, 'timescale (yr):', model%ocean_data%thermal_forcing_anomaly_timescale
             print*, 'fraction:', anomaly_fraction
             print*, 'current TF anomaly (deg):', tf_anomaly
             if (model%ocean_data%thermal_forcing_anomaly_timescale /= 0.0d0) then
                print*, 'anomaly applied to basin number', model%ocean_data%thermal_forcing_anomaly_basin
             endif
          endif
       else
          tf_anomaly = 0.0d0
          tf_anomaly_basin = 0
       endif

       call glissade_bmlt_float_thermal_forcing(&
            model%options%bmlt_float_thermal_forcing_param, &
            model%options%ocean_data_extrapolate,  &
            parallel,                              &
            ewn,                nsn,               &
            dew*len0,           dns*len0,          & ! m
            itest,     jtest,   rtest,             &
            ice_mask,                              &
            ocean_mask,                            &
            model%geometry%marine_connection_mask, &
            model%geometry%f_ground_cell,          &
            model%geometry%thck*thk0,              & ! m
            model%geometry%lsrf*thk0,              & ! m
            model%geometry%topg*thk0,              & ! m
            model%ocean_data,                      &
            model%basal_melt%bmlt_float,           &
            tf_anomaly_in = tf_anomaly,            & ! deg C
            tf_anomaly_basin_in = tf_anomaly_basin)

       ! Convert bmlt_float from SI units (m/s) to scaled model units
       model%basal_melt%bmlt_float(:,:) = model%basal_melt%bmlt_float(:,:) * tim0/thk0

    else  ! other options include BMLT_FLOAT_CONSTANT, BMLT_FLOAT_MISMIP, &
          !  BMLT_FLOAT_DEPTH, and BMLT_FLOAT_MISOMIP
          !TODO - Call separate subroutines for each of these options?

       call glissade_basal_melting_float(model%options%whichbmlt_float,                         &
                                         parallel,                                              &
                                         ewn,                        nsn,                       &
                                         model%numerics%dew*len0,    model%numerics%dns*len0,   &
                                         itest,                      jtest,                     &
                                         rtest,                                                 &
                                         model%general%x1,                                      & ! m
                                         model%geometry%thck*thk0,                              & ! m
                                         model%geometry%lsrf*thk0,                              & ! m
                                         model%geometry%topg*thk0,                              & ! m
                                         model%climate%eus*thk0,                                & ! m
                                         model%basal_melt,                                      & ! bmlt_float in m/s
                                         model%ocean_data)

       ! Convert bmlt_float from SI units (m/s) to scaled model units
       model%basal_melt%bmlt_float(:,:) = model%basal_melt%bmlt_float(:,:) * tim0/thk0

    endif  ! whichbmlt_float


    ! If desired, add a bmlt_anomaly field.
    ! This is done for the initMIP Greenland and Antarctic experimennts.

    if (model%options%enable_bmlt_anomaly) then

       ! Compute the previous time
       ! Note: When being ramped up, the anomaly is not incremented until after the final time step of the year.
       !       This is the reason for passing the previous time to the subroutine.
       previous_time = model%numerics%time - model%numerics%dt * tim0/scyr

       ! Add the bmlt_float anomaly where ice is present and floating
       call glissade_add_2d_anomaly(model%basal_melt%bmlt_float,              &   ! scaled model units
                                    model%basal_melt%bmlt_float_anomaly,      &   ! scaled model units
                                    model%basal_melt%bmlt_anomaly_timescale,  &   ! yr
                                    previous_time)                                ! yr

    endif

    ! Zero out bmlt_float in ice-free ocean cells.
    ! Note: Do not do this for the thermal_forcing option, because this option allows nonzero bmlt_float
    !       in ocean cells adjacent to floating cells.
    ! TODO: Look at other options and decide which ones need this logic.
    if (model%options%whichbmlt_float /= BMLT_FLOAT_THERMAL_FORCING) then
       where (ocean_mask == 1)
          model%basal_melt%bmlt_float = 0.0d0
       endwhere
    endif

    ! Reduce or zero out bmlt_float in cells with fully or partly grounded ice
    !TODO - Write a subroutine to do this calculation (in glissade_ground or glissade_bmlt_float?)
    !       The same subroutine could be called from the inversion solver.

    if (model%options%which_ho_ground == HO_GROUND_GLP_DELUXE) then

       ! Reduce bmlt_float in partly or fully grounded cells based on f_ground_cell

       if (model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_FLOATING_FRAC) then

          ! Multiply bmlt_float by the fraction of the cell that is floating.
          ! Cells that are fully grounded will have bmlt_float = 0.
          ! This option ensures smooth changes in bmlt_float as the GL migrates.
          ! However, it might allow spurious melting of grounded ice near the GL.

          where (model%geometry%f_ground_cell > 0.0d0)
             model%basal_melt%bmlt_float = model%basal_melt%bmlt_float   &
                                         * (1.0d0 - model%geometry%f_ground_cell)
          endwhere

       elseif (model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_ZERO_GROUNDED) then

          ! Where f_ground_cell > 0, set bmlt_float = 0.
          ! Cells that are even partly grounded will have bmlt_float = 0.
          ! This option ensures no spurious melting of grounded ice near the GL.
          ! However, it may underestimate melting of floating ice near the GL, especially on coarser grids.

          where (model%geometry%f_ground_cell > tiny(0.0d0))
             model%basal_melt%bmlt_float = 0.0d0
          endwhere

       elseif (model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_NO_GLP) then

          ! Zero out bmlt_float in grounded cells based on floating_mask.
          ! Note: CISM typically would not be run with this combination, but it is included for generality.

          where (floating_mask == 0)
             model%basal_melt%bmlt_float = 0.0d0
          endwhere

       endif  ! which_ho_ground_bmlt 

    else

       ! Zero out bmlt_float in grounded cells based on floating_mask
       where (floating_mask == 0)
          model%basal_melt%bmlt_float = 0.0d0
       endwhere

    endif

    ! Reduce basal melting in shallow cavities if bmlt_cavity_h0 > 0.
    ! The tanh function follows Asay-Davis et al. (2016), Eqs. 14 and 17.
    ! Note: model%basal_melt%bmlt_cavity_h0 has units of m.
    ! Note: For BMLT_FLOAT_MISMIP, this reduction is done in subroutine glissade_basal_melting_float
    !       based on model%basal_melt%bmlt_float_h0 and should not be repeated here.

    if (model%basal_melt%bmlt_cavity_h0 > 0.0d0 .and.  &
        model%options%whichbmlt_float /= BMLT_FLOAT_MISMIP) then

       ! TODO: Make sure lsrf is up to date.

       h_cavity = max(model%geometry%lsrf - model%geometry%topg, 0.0d0) * thk0  ! cavity thickness (m)

       if (this_rank==rtest .and. verbose_bmlt_float) then
          print*, ' '
          print*, 'Reduce bmlt_float in shallow cavities, bmlt_cavity_h0 =', &
               model%basal_melt%bmlt_cavity_h0
          print*, ' '
          print*, 'original bmlt_float (m/yr):'
          write(6,'(a6)',advance='no') '      '
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') i
          enddo
          print*, ' '
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%basal_melt%bmlt_float(i,j) * thk0*scyr/tim0
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'h_cavity:'
          print*, ' '
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') h_cavity(i,j)
             enddo
             print*, ' '
          enddo
          print*, ' '
          print*, 'Adjustment factor:'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') &
                     tanh(h_cavity(i,j)/model%basal_melt%bmlt_cavity_h0)
             enddo
             print*, ' '
          enddo
       endif   ! verbose_bmlt_float

       where (h_cavity > 0.0d0)
          model%basal_melt%bmlt_float = model%basal_melt%bmlt_float * &
               tanh(h_cavity/model%basal_melt%bmlt_cavity_h0)
       elsewhere
          model%basal_melt%bmlt_float = 0.0d0
       endwhere

       if (this_rank==rtest .and. verbose_bmlt_float) then
          print*, ' '
          print*, 'reduced bmlt_float (m/yr):'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%basal_melt%bmlt_float(i,j) * thk0*scyr/tim0
             enddo
             print*, ' '
          enddo
       endif

    endif   ! bmlt_cavity_h0 > 0

    !WHL - debug
    if (this_rank==rtest .and. verbose_bmlt_float) then
       print*, ' '
       print*, 'After glissade_bmlt_float_solve, which_ho_ground_bmlt =', model%options%which_ho_ground_bmlt

       if (model%options%which_ho_ground == HO_GROUND_GLP_DELUXE) then
          print*, ' '
          print*, '1 - f_ground_cell:'
          write(6,'(a6)',advance='no') '      '
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') i
          enddo
          write(6,*) ' '
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') 1.0d0 - model%geometry%f_ground_cell(i,j)
             enddo
             write(6,*) ' '
          enddo
       else
          print*, ' '
          print*, 'floating_mask:'
          write(6,'(a6)',advance='no') '      '
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') i
          enddo
          write(6,*) ' '
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') floating_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif   ! which_ho_ground

       print*, ' '
       print*, 'bmlt_float (m/yr):'
       write(6,'(a6)',advance='no') '      '
       do i = itest-3, itest+3
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') model%basal_melt%bmlt_float(i,j) * thk0*scyr/tim0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '

       print*, 'Done in glissade_bmlt_float_solve'

    endif  ! this_rank = rtest and verbose_bmlt_float

  end subroutine glissade_bmlt_float_solve

!=======================================================================

  subroutine glissade_thermal_solve(model, dt)

    ! Do the vertical thermal solve.
    ! First call a driver subroutine for vertical temperature or enthalpy evolution,
    ! and then update the basal water.
    use cism_parallel, only: parallel_type, parallel_halo

    use glimmer_paramets, only: tim0, thk0, len0
    use glimmer_physcon, only: rhow, rhoi, scyr
    use glissade_therm, only: glissade_therm_driver
    use glissade_basal_water, only: glissade_calcbwat, glissade_bwat_flux_routing
    use glissade_transport, only: glissade_add_2d_anomaly
    use glissade_grid_operators, only: glissade_vertical_interpolate
    use glissade_masks, only: glissade_get_masks

    !WHL - debug
    use cism_parallel, only: parallel_reduce_max
    use glissade_glacier, only : verbose_glacier

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    real(dp), intent(in) :: dt   ! time step (s)
    
    ! unscaled model parameters (SI units)
    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
       bmlt_ground_unscaled,   & ! basal melt rate for grounded ice (m/s)
       bwat_unscaled             ! basal water thickness (m)

    real(dp) :: previous_time       ! time (yr) at the start of this time step
                                    ! (model%numerics%time is the time at the end of the step.)

    integer :: i, j, up
    integer :: itest, jtest, rtest

    integer, dimension(model%general%ewn, model%general%nsn) ::   &
         ice_mask,              & ! = 1 if ice is present (thck > thklim_temp), else = 0
         floating_mask,         & ! = 1 if ice is present (thck > thklim_temp) and floating, else = 0
         ocean_mask,            & ! = 1 if topg is below sea level and ice is absent, else = 0
         bwat_mask                ! = 1 for cells through which basal water is routed, else = 0

    !WHL - debug
    real(dp) :: head_max

    type(parallel_type) :: parallel   ! info for parallel communication

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    parallel = model%parallel

    call t_startf('glissade_thermal_solve')

    ! Optionally, add an anomaly to the surface air temperature
    ! Typically, artm_corrected = artm, but sometimes (e.g., for ISMIP6 forcing experiments),
    !  it includes a time-dependent anomaly.
    ! Note that artm itself does not change in time.

    ! initialize
    model%climate%artm_corrected(:,:) = model%climate%artm(:,:)

    if (model%options%enable_artm_anomaly) then

       ! Note: When being ramped up, the anomaly is not incremented until after the final time step of the year.
       !       This is the reason for passing the previous time to the subroutine.
       previous_time = model%numerics%time - model%numerics%dt * tim0/scyr

       call glissade_add_2d_anomaly(model%climate%artm_corrected,          &   ! degC
                                    model%climate%artm_anomaly,            &   ! degC
                                    model%climate%artm_anomaly_timescale,  &   ! yr
                                    previous_time)                             ! yr

       if (verbose_glissade .and. this_rank==rtest) then
          i = itest
          j = jtest
          print*, 'i, j, previous_time, artm, artm anomaly, corrected artm (deg C):', &
               i, j, previous_time, model%climate%artm(i,j), model%climate%artm_anomaly(i,j), &
               model%climate%artm_corrected(i,j)
       endif

    endif

    if (main_task .and. verbose_glissade) print*, 'Call glissade_therm_driver'

    ! Downscale artm to the current surface elevation if needed.
    ! Depending on the value of artm_input_function, artm might be dependent on the upper surface elevation.
    ! The options are:
    ! (0) artm(x,y); no dependence on surface elevation
    ! (1) artm(x,y) + d(artm)/dz(x,y) * dz; artm depends on input field at reference elevation, plus vertical correction
    ! (2) artm(x,y,z); artm obtained by linear interpolation between values prescribed at adjacent vertical levels
    ! (3) artm(x,y) adjusted with a uniform lapse rate
    ! For options (1) - (3), the elevation-dependent artm is computed here.

    if (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_GRADZ) then

       ! compute artm by a lapse-rate correction to the reference value
       model%climate%artm(:,:) = model%climate%artm_ref(:,:) + &
            (model%geometry%usrf(:,:)*thk0 - model%climate%usrf_ref(:,:)) * model%climate%artm_gradz(:,:)

    elseif (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XYZ) then

       ! Note: With linear_extrapolate_in = T, the values outside the range are obtained by linear extrapolation
       !        from the top two or bottom two values.
       !       For temperature, which varies roughly linearly with elevation, this is more accurate
       !        than simply extending the top and bottom values.
       !       This call includes a halo update.

       call glissade_vertical_interpolate(model%general%ewn,      model%general%nsn,         &
                                          model%climate%nlev_smb, model%climate%smb_levels,  &
                                          model%geometry%usrf,                               &
                                          model%climate%artm_3d,                             &
                                          model%climate%artm,                                &
                                          linear_extrapolate_in = .true.)

    elseif (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_LAPSE) then

       ! compute artm by a lapse-rate correction to artm_ref
       ! T_lapse is defined as positive for T decreasing with height
       ! Note: This option is currently used for glaciers lapse rate adjustments

       model%climate%artm(:,:) = model%climate%artm_ref(:,:) - &
            (model%geometry%usrf(:,:)*thk0 - model%climate%usrf_ref(:,:)) * model%climate%t_lapse
       if (verbose_glacier .and. this_rank == rtest) then
          i = itest; j = jtest
          print*, ' '
          print*, 'rank, i, j, usrf, usrf_ref, dz:', this_rank, i, j, &
               model%geometry%usrf(i,j)*thk0, model%climate%usrf_ref(i,j), &
               model%geometry%usrf(i,j)*thk0 - model%climate%usrf_ref(i,j)
          print*, '   artm_ref, artm:', model%climate%artm_ref(i,j), model%climate%artm(i,j)
       endif

    endif   ! artm_input_function

    call parallel_halo(model%climate%artm, parallel)

    ! Note: glissade_therm_driver uses SI units
    !       Output arguments are temp, waterfrac, bpmp and bmlt_ground
    call glissade_therm_driver (model%options%whichtemp,                                      &
                                model%options%temp_init,                                      &
                                dt,                                                           & ! s
                                model%parallel,                                               &
                                model%general%ewn,          model%general%nsn,                &
                                model%general%upn,                                            &
                                model%numerics%idiag_local, model%numerics%jdiag_local,       &
                                model%numerics%rdiag_local,                                   &
                                model%numerics%sigma,       model%numerics%stagsigma,         &
                                model%numerics%dups,                                          &
                                model%numerics%thklim_temp*thk0,                              & ! m
                                model%geometry%thck*thk0,                                     & ! m
                                model%geometry%topg*thk0,                                     & ! m
                                model%geometry%lsrf*thk0,                                     & ! m
                                model%climate%eus*thk0,                                       & ! m
                                model%climate%artm_corrected,                                 & ! deg C
                                model%options%which_ho_ground,                                &
                                model%geometry%f_ground_cell,                                 & ! [0,1]
                                model%temper%bheatflx,      model%temper%bfricflx,            & ! W/m2
                                model%temper%dissip,                                          & ! deg/s
                                model%temper%pmp_threshold,                                   & ! deg C
                                model%basal_hydro%bwat*thk0,                                  & ! m
                                model%temper%temp,                                            & ! deg C
                                model%temper%waterfrac,                                       & ! unitless
                                model%temper%bpmp,                                            & ! deg C
                                model%temper%btemp_ground,                                    & ! deg C
                                model%temper%btemp_float,                                     & ! deg C
                                bmlt_ground_unscaled)                                           ! m/s

    ! Update basal hydrology, if needed
    ! Note: glissade_calcbwat uses SI units

    if (main_task .and. verbose_glissade) print*, 'Call glissade_calcbwat'

    ! convert bwat to SI units for input to glissade_calcbwat
    bwat_unscaled(:,:) = model%basal_hydro%bwat(:,:) * thk0

    !TODO - Move the following calls to a new basal hydrology solver?

    if (model%options%which_ho_bwat == HO_BWAT_FLUX_ROUTING) then

       !WHL - Temporary code for debugging: Make up a simple basal melt field.
!       model%basal_hydro%head(:,:) = &
!            model%geometry%thck(:,:)*thk0 + (rhow/rhoi)*model%geometry%topg(:,:)*thk0
!       head_max = maxval(model%basal_hydro%head)  ! max on local processor
!       head_max = parallel_reduce_max(head_max)   ! global max
!       do j = 1, model%general%nsn
!          do i = 1, model%general%ewn
!             if (head_max - model%basal_hydro%head(i,j) < 1000.d0) then
!!             if (head_max - model%basal_hydro%head(i,j) < 200.d0) then
!                bmlt_ground_unscaled(i,j) = 1.0d0/scyr    ! units are m/s
!             else
!                bmlt_ground_unscaled(i,j) = 0.0d0
!             endif
!          enddo
!       enddo

       ! Compute some masks needed below

       call glissade_get_masks(&
            model%general%ewn,    model%general%nsn,      &
            model%parallel,                               &
            model%geometry%thck,  model%geometry%topg,    &
            model%climate%eus,    model%numerics%thklim,  &
            ice_mask,                                     &
            floating_mask = floating_mask,                &
            ocean_mask = ocean_mask)

       ! Compute a mask that sets the domain for flux routing.
       ! Cells excluded from the domain are:
       ! (1) floating or ocean cells
       ! (2) cells at the edge of the global domain
       ! (3) ice-free cells in the region where the SMB is overwritten
       !     by a prescribed negative value (on the assumption that
       !     such cells are supposed to be beyond the ice margin)
       !
       ! Note: Cells with bwat_mask = 0 can have bwat_flux > 0 if they receive water
       !  from adjacent cells with bwat_mask = 1.
       ! But once the flux reaches a cell with bwat_mask = 0, it is not routed further.
       ! Thus, the total flux in cells with bwat_mask = 0 should be equal to the
       !  total input flux of basal meltwater.

       bwat_mask = 1   ! initially, include the entire domain

       where (floating_mask == 1 .or. ocean_mask == 1 .or.  &
              model%general%global_edge_mask == 1)
          bwat_mask = 0
       endwhere

       if (model%options%overwrite_acab /= OVERWRITE_ACAB_NONE .and. &
           model%climate%overwrite_acab_value < 0.0d0) then
          where (model%climate%overwrite_acab_mask == 1 .and. &
                 model%geometry%thck < model%numerics%thklim)
             bwat_mask = 0
          endwhere
       endif

       !WHL - debug - Set mask = 0 where thck = 0 for dome test
!       where (model%geometry%thck == 0)
!          bwat_mask = 0
!       endwhere

       call parallel_halo(bwat_mask, parallel)

       ! Compute the steady-state basal water flux based on a flux-routing scheme

       call glissade_bwat_flux_routing(&
            model%general%ewn,       model%general%nsn,       &
            model%numerics%dew*len0, model%numerics%dns*len0, &  ! m
            model%parallel,                                   &
            itest, jtest, rtest,                              &
            model%options%ho_flux_routing_scheme,             &
            model%geometry%thck*thk0,                         &  ! m
            model%geometry%topg*thk0,                         &  ! m
            model%numerics%thklim_temp*thk0,                  &  ! m
            bwat_mask,                                        &
            floating_mask,                                    &
            bmlt_ground_unscaled,                             &  ! m/s
            model%basal_hydro%bwatflx,                        &  ! m^3/s
            model%basal_hydro%head)                              ! m

    else  ! simpler basal water options

       call glissade_calcbwat(model%options%which_ho_bwat,      &
                              model%basal_hydro,                &
                              dt,                               &  ! s
                              model%geometry%thck*thk0,         &  ! m
                              model%numerics%thklim_temp*thk0,  &  ! m
                              bmlt_ground_unscaled,             &  ! m/s
                              bwat_unscaled)                       ! m

    endif

    ! convert bmlt and bwat from SI units (m/s and m) to scaled model units
    model%basal_melt%bmlt_ground(:,:) = bmlt_ground_unscaled(:,:) * tim0/thk0
    model%basal_hydro%bwat(:,:) = bwat_unscaled(:,:) / thk0

    ! Update tempunstag as sigma weighted interpolation from temp to layer interfaces
    do up = 2, model%general%upn-1
      model%temper%tempunstag(up,:,:) = model%temper%temp(up-1,:,:) +       &
           (model%temper%temp(up,:,:) - model%temper%temp(up-1,:,:)) *      &
           (model%numerics%sigma(up) - model%numerics%stagsigma(up-1)) /    &
           (model%numerics%stagsigma(up) - model%numerics%stagsigma(up-1))
    end do
    ! boundary conditions are identical on both grids, but temp starts at index 0
    model%temper%tempunstag(1,:,:) = model%temper%temp(0,:,:)
    model%temper%tempunstag(model%general%upn,:,:) = model%temper%temp(model%general%upn,:,:)

    !------------------------------------------------------------------------ 
    ! Halo updates
    !------------------------------------------------------------------------ 
    
    ! Note: bwat is needed in halos to compute effective pressure if which_ho_effecpress = HO_EFFECPRESS_BWAT
    call parallel_halo(model%basal_hydro%bwat, parallel)

    call t_stopf('glissade_thermal_solve')
    
  end subroutine glissade_thermal_solve

!=======================================================================

  subroutine glissade_thickness_tracer_solve(model)

    ! ------------------------------------------------------------------------ 
    ! Calculate ice thickness and tracer evolution, including horizontal transport and surface and basal mass balance.
    ! MJH: This subroutine uses velocity from the previous time step, which is appropriate for a Forward Euler time-stepping scheme.
    ! WHL: We used to have EVOL_NO_THICKNESS = -1 as a Glide option, used to hold the ice surface elevation fixed during CESM runs. 
    !      This option has been replaced by a Glint/Glad option, evolve_ice.
    !      We now have EVOL_NO_THICKESS = 5 as a glissade option.  It is used to hold the ice surface elevation fixed
    !       while allowing temperature to evolve, which can be useful for model spinup.  This option might need more testing.
    ! Note: This subroutine calls glissade_inversion_bmlt_float, because it is convenient to invert for bmlt_float
    !       after horizontal transport and before applying the surface and basal mass balance.
    ! ------------------------------------------------------------------------ 

    use cism_parallel, only: parallel_type, parallel_halo, parallel_halo_tracers,  &
         staggered_parallel_halo, parallel_reduce_max

    use glimmer_paramets, only: eps11, eps08, tim0, thk0, vel0, len0
    use glimmer_physcon, only: rhow, rhoi, scyr
    use glimmer_scales, only: scale_acab
    use glissade_therm, only: glissade_temp2enth, glissade_enth2temp
    use glissade_transport, only: glissade_mass_balance_driver, &
                                  glissade_transport_driver, &
                                  glissade_check_cfl,  &
                                  glissade_transport_setup_tracers, &
                                  glissade_transport_finish_tracers, &
                                  glissade_overwrite_acab,  &
                                  glissade_add_2d_anomaly
    use glissade_masks, only: glissade_get_masks, glissade_extend_mask, &
                              glissade_calving_front_mask
    use glissade_inversion, only: verbose_inversion
    use glissade_bmlt_float, only: verbose_bmlt_float
    use glissade_calving, only: verbose_calving
    use glissade_grid_operators, only: glissade_vertical_interpolate
    use glissade_glacier, only: verbose_glacier, glissade_glacier_smb, &
                                glissade_glacier_advance_retreat
    use glide_stop, only: glide_finalise

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! --- Local variables ---

    integer :: sc  ! subcycling index

    ! temporary arrays in SI units
    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
       thck_unscaled,     & ! ice thickness (m)
       topg_unscaled,     & ! bedrock topography (m)
       thck_new_unscaled, & ! expected new ice thickness, after mass balance (m)
       acab_unscaled,     & ! surface mass balance (m/s)
       bmlt_unscaled        ! = bmlt (m/s) if basal mass balance is included in continuity equation, else = 0

    ! masks
    integer, dimension(model%general%ewn, model%general%nsn) ::   &
       ice_mask,             & ! = 1 if thck > 0, else = 0
       floating_mask,        & ! = 1 where ice is present and floating, else = 0
       ocean_mask,           & ! = 1 if topg is below sea level and thck = 0, else = 0
       land_mask,            & ! = 1 if topg is at or above sea level, else = 0
       calving_front_mask,   & ! = 1 where ice is floating and borders an ocean cell, else = 0
       active_ice_mask,      & ! = 1 for cells that are dynamically active, else = 0
       protected_mask,       & ! = 1 for cells that are protected from later calving (HO_CALVING_FRONT_SUBGRID)
       extended_ice_sheet_mask ! extension of ice_sheet_mask to include neighbor cells

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
       thck_flotation,       & ! thickness at which ice is exactly floating
       thck_calving_front,   & ! effective thickness of ice at the calving front
       effective_areafrac      ! effective fractional area of ice at the calving front

    real(dp) :: previous_time       ! time (yr) at the start of this time step
                                    ! (model%numerics%time is the time at the end of the step.)

    real(dp) :: advective_cfl       ! advective CFL number
                                    ! If advective_cfl > 1, the model is unstable without subcycling
    real(dp) :: dt_transport        ! time step (s) for transport; = model%numerics%dt*tim0 by default

    integer :: nsubcyc              ! number of times to subcycle advection

    logical :: do_upwind_transport  ! logical for whether transport code should do upwind transport or incremental remapping
                                    ! set to true for EVOL_UPWIND, else = false

    integer :: ntracers             ! number of tracers to be transported

    integer :: i, j, k
    integer :: ewn, nsn, upn, nlev_smb
    integer :: itest, jtest, rtest

    type(parallel_type) :: parallel   ! info for parallel communication

    !WHL - debug
    real(dp) :: local_maxval, global_maxval
    character(len=100) :: message

    logical, parameter :: verbose_smb = .false.

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn
    nlev_smb = model%climate%nlev_smb

    parallel = model%parallel

    select case(model%options%whichevol)

    case(EVOL_INC_REMAP, EVOL_UPWIND, EVOL_NO_THICKNESS) 

       if (model%options%whichevol == EVOL_UPWIND) then
          do_upwind_transport = .true.
       else
          do_upwind_transport = .false.
       endif

       !-------------------------------------------------------------------------
       ! First transport ice thickness and temperature, given the horizontal velocity (u,v).
       ! Then apply the surface and basal mass balance in each grid cell.
       ! Note: The main reason to do both transport and mass balance in one subroutine is
       !  that both operations require tracer array setup and cleanup (e.g., copying the 
       !  various tracer fields into generic tracer arrays). With this arrangement,
       !  the tracer operations need to be done only once.
       !-------------------------------------------------------------------------
       ! MJH: I put the no thickness evolution option here so that it is still possible 
       !      (but not required) to use IR to advect temperature when thickness evolution is turned off.
       ! TODO  MJH If we really want to support no evolution, then we may want to implement it so that IR does not occur 
       !       at all - right now a run can fail because of a CFL violation in IR even if evolution is turned off.
       !       Do we want to support temperature evolution without thickness evolution?
       !       If so, then the current implementation may be preferred approach.

       call t_startf('inc_remap_driver')

       if (verbose_glissade .and. main_task) print*, 'Compute dH/dt'

       call t_startf('glissade_transport_driver')

       ! ------------------------------------------------------------------------
       ! Compute some masks prior to horizontal transport.
       ! ------------------------------------------------------------------------

       call glissade_get_masks(ewn,              nsn,              &
                               parallel,                           &
                               model%geometry%thck*thk0,           &   ! m
                               model%geometry%topg*thk0,           &   ! m
                               model%climate%eus*thk0,             &   ! m
                               0.0d0,                              &   ! thklim = 0
                               ice_mask,                           &
                               floating_mask = floating_mask,      &
                               ocean_mask = ocean_mask,            &
                               land_mask = land_mask,              &
                               active_ice_mask = active_ice_mask)

       ! If using a subgrid calving-front scheme, then recompute active_ice_mask.
       if (model%options%which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then

          call glissade_calving_front_mask(ewn,                    nsn,              &
                                           model%options%which_ho_calving_front,     &
                                           parallel,                                 &
                                           model%geometry%thck*thk0,                 &   ! m
                                           model%geometry%topg*thk0,                 &   ! m
                                           model%climate%eus*thk0,                   &   ! m
                                           ice_mask,               floating_mask,    &
                                           ocean_mask,             land_mask,        &
                                           calving_front_mask,     thck_calving_front, &
                                           active_ice_mask = active_ice_mask)
       endif

       ! For the enthalpy option, derive enthalpy from temperature and waterfrac.
       ! Must transport enthalpy rather than temperature/waterfrac to conserve energy.
       !TODO: Move this code to glissade_transport_setup_tracers?

       if (model%options%whichtemp == TEMP_ENTHALPY) then  ! Use IR to transport enthalpy
          ! Note: glissade_temp2enth expects SI units
          do j = 1, nsn
             do i = 1, ewn
                call glissade_temp2enth (model%numerics%stagsigma(1:upn-1),        &
                                         model%temper%temp(0:upn,i,j),     model%temper%waterfrac(1:upn-1,i,j),   &
                                         model%geometry%thck(i,j)*thk0,    model%temper%enthalpy(0:upn,i,j))
             enddo
          enddo
       endif    ! TEMP_ENTHALPY

       ! copy tracers (temp/enthalpy, etc.) into model%geometry%tracers
       call glissade_transport_setup_tracers (model)

       ! temporary in/out arrays in SI units
       thck_unscaled(:,:) = model%geometry%thck(:,:) * thk0
       topg_unscaled(:,:) = model%geometry%topg(:,:) * thk0

       ! pre-transport halo updates for thickness and tracers
       call parallel_halo(thck_unscaled, parallel)
       call parallel_halo(topg_unscaled, parallel)
       call parallel_halo_tracers(model%geometry%tracers, parallel)
       call parallel_halo_tracers(model%geometry%tracers_usrf, parallel)
       call parallel_halo_tracers(model%geometry%tracers_lsrf, parallel)

       ! pre-transport halo updates for velocity
       ! Velocity update might be needed if velo has not been updated in the halo since the previous diagnostic solve.
       !  (just to be on the safe side).

       call staggered_parallel_halo(model%velocity%uvel, parallel)
       call staggered_parallel_halo(model%velocity%vvel, parallel)

       ! --- Determine CFL limits ---
       ! Note: We are using the subcycled dt here (if subcycling is on).
       !  (See note above about the EVOL_NO_THICKNESS option and how it is affected by a CFL violation)
       !  stagthck, dusrfdew/ns and u/vvel need to be from the previous time step (and are at this point)
       ! Note: If using adaptive subcycling (with adaptive_cfl_threshold > 0), then dt_transport should
       !       be equal to dt (which is the case by default).
       !TODO - Remove the dt_transport option and simply rely on adaptive subcycling as needed?

       call glissade_check_cfl(ewn,                       nsn,                       upn-1,                    &
                               model%numerics%dew * len0, model%numerics%dns * len0, model%numerics%sigma,     &
                               parallel,                                                                       &
                               model%geomderv%stagthck * thk0,                                                 &
                               model%geomderv%dusrfdew*thk0/len0, model%geomderv%dusrfdns*thk0/len0,           &
                               model%velocity%uvel * scyr * vel0, model%velocity%vvel * scyr * vel0,           &
                               model%numerics%dt_transport * tim0 / scyr,                                      &
                               model%numerics%adaptive_cfl_threshold,                                          &
                               model%numerics%adv_cfl_dt,         model%numerics%diff_cfl_dt)

       ! Set the transport timestep.
       ! The timestep is model%numerics%dt by default, but optionally can be reduced for subcycling

       !WHL - debug
!      if (main_task) then
!         print*, 'Checked advective CFL threshold'
!         print*, 'model dt (yr) =', model%numerics%dt * tim0/scyr
!         print*, 'adv_cfl_dt    =', model%numerics%adv_cfl_dt
!      endif

       advective_cfl = model%numerics%dt*(tim0/scyr) / model%numerics%adv_cfl_dt

       if (model%numerics%adaptive_cfl_threshold > 0.0d0) then

          ! subcycle the transport when advective_cfl exceeds the threshold

          if (advective_cfl > model%numerics%adaptive_cfl_threshold) then

             ! compute the number of subcycles
             ! If advective_cfl > advective_cfl_threshold, then nsubcyc >= 2
             ! The larger the ratio, the larger the value of nsubcyc
             nsubcyc = ceiling(advective_cfl / model%numerics%adaptive_cfl_threshold)

             if (main_task) then
                print*, 'WARNING: adv_cfl_dt exceeds threshold; CFL =', advective_cfl
                print*, 'Ratio =', advective_cfl / model%numerics%adaptive_cfl_threshold
                print*, 'nsubcyc =', nsubcyc
             endif

          else
             nsubcyc = 1
          endif
          dt_transport = model%numerics%dt * tim0 / real(nsubcyc,dp)   ! convert to s

       else  ! no adaptive subcycling

          advective_cfl = model%numerics%dt*(tim0/scyr) / model%numerics%adv_cfl_dt

          ! If advective_cfl exceeds 1.0, then abort cleanly.  Otherwise, set dt_transport and proceed.
          ! Note: Usually, it would be enough to write a fatal abort message.
          !       The call to glide_finalise was added to allow CISM to finish cleanly when running
          !        a suite of automated stability tests, e.g. with the stabilitySlab.py script.
          if (advective_cfl > 1.0d0) then
             if (main_task) print*, 'advective CFL violation; call glide_finalise and exit cleanly'
             call glide_finalise(model, crash=.true.)
             stop
          else
             nsubcyc = model%numerics%subcyc
             dt_transport = model%numerics%dt_transport * tim0  ! convert to s
          endif

       endif

       !-------------------------------------------------------------------------
       ! Compute horizontal transport of mass and tracers, subcycling as needed.
       !-------------------------------------------------------------------------

       do sc = 1, nsubcyc

          if (nsubcyc > 1 .and. main_task) write(*,*) 'Subcycling transport: Cycle', sc

          ! Call the transport driver subroutine.
          ! (includes a halo update for thickness: thck_unscaled in this case)
          !
          ! Note: This subroutine assumes SI units:
          !       * dt (s)
          !       * dew, dns, thck (m)
          !       * uvel, vvel (m/s)
          !       Since thck has intent(inout), we create and pass a temporary array (thck_unscaled) with units of m.
          ! Note: tracers_ursf and tracers_lsrf are not transported, but they provide upper and lower BCs
          !       for vertical remapping. They are intent(in).

          call glissade_transport_driver(dt_transport,                                         &  ! s
                                         model%numerics%dew * len0, model%numerics%dns * len0, &
                                         ewn,          nsn,         upn-1,                     &
                                         model%numerics%sigma,                                 &
                                         parallel,                                             &
                                         model%velocity%uvel(:,:,:) * vel0,                    &  ! m/s
                                         model%velocity%vvel(:,:,:) * vel0,                    &  ! m/s
                                         thck_unscaled(:,:),                                   &  ! m
                                         model%geometry%ntracers,                              &
                                         model%geometry%tracers(:,:,:,:),                      &
                                         model%geometry%tracers_usrf(:,:,:),                   &
                                         model%geometry%tracers_lsrf(:,:,:),                   &
                                         model%options%which_ho_vertical_remap,                &
                                         upwind_transport_in = do_upwind_transport)

          ! halo updates for thickness and tracers
          !TODO: For outflow and no_ice BCs where halo routines can remove ice near the global boundary,
          !      keep track of the mass of ice removed, and incorporate it into the global mass balance.
          call parallel_halo(thck_unscaled, parallel)
          call parallel_halo_tracers(model%geometry%tracers, parallel)

       enddo     ! subcycling of transport

       if ((verbose_inversion .or. verbose_glissade .or. verbose_calving) .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'After glissade_transport_driver, thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck_unscaled(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       ! If using a subgrid calving_front scheme (but not using a no-advance calving mask),
       !  then identify thin ice that was transported beyond the CF to ice-free cells without active neighbors.
       ! In this case, model%calving%calving_mask is computed here and applied later, in glissade_calving_solve.
       !
       ! Note: We are trying to avoid spurious CF advance in the following situation:
       !  -------------
       !  |     |     |
       !  |  1  |  2  |
       !  -------------
       !  |     |     |
       !  |  3  |  4  |
       !  -------------
       !
       ! Suppose cell 1 is an active, floating interior cell; cells 2 and 3 are inactive CF cells that are filling;
       !  and cell 4 is ice-free ocean.
       ! If the transport scheme delivers ice to cell 4, then cell 4 can become an inactive CF cell,
       !  making cell 3 an active interior cell, even though it is much thinner than other interior cells.
       ! To prevent this, ice is removed from cells that do not meet at least one of the following criteria
       !  before transport:
       ! (1) ice is present (ice_mask = 1)
       ! (2) land-based
       ! (3) shares at least one edge with an active cell

       !TODO - Remove the redundant CALVING_GRID_MASK option
       if (model%options%which_ho_calving_front == HO_CALVING_FRONT_SUBGRID .and.   &
           (model%options%whichcalving /= CALVING_GRID_MASK .and. .not.model%options%apply_calving_mask)) then

          ! compute a mask of protected cells
          ! Protect cells where ice was present before advection, and protect land cells

          protected_mask(:,:) = 0
          where (ice_mask == 1 .or. land_mask == 1)
             protected_mask = 1
          endwhere

          do j = 2, nsn-1
             do i = 2, ewn-1
                if (active_ice_mask(i-1,j) == 1 .or. active_ice_mask(i+1,j) == 1 .or. &
                    active_ice_mask(i,j-1) == 1 .or. active_ice_mask(i,j+1) == 1) then
                   protected_mask(i,j) = 1
                endif
             enddo
          enddo

          ! Identify cells where thin ice should be removed.
          ! The mask is applied later, in glissade_calving_solve.
          where (protected_mask == 0 .and. thck_unscaled > 0.0d0)
             model%calving%calving_mask = 1
          elsewhere
             model%calving%calving_mask = 0
          endwhere

          call parallel_halo(model%calving%calving_mask, parallel)

          if (verbose_calving .and. this_rank == rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'After transport, calving_mask:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(i10)',advance='no') model%calving%calving_mask(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

       endif   ! subgrid calving front

       !-------------------------------------------------------------------------
       ! Prepare the surface and basal mass balance terms.
       ! Note: The basal mass balance has been computed in subroutine glissade_bmlt_float_solve.
       !-------------------------------------------------------------------------

       !TODO - Put this SMB logic in a subroutine.
       ! Downscale acab to the current surface elevation if needed.
       ! Depending on the value of smb_input_function, the SMB (i.e., acab) might be dependent on the upper surface elevation.
       ! The options are:
       ! (0) SMB(x,y); no dependence on surface elevation
       ! (1) SMB(x,y) + dSMB/dz(x,y) * dz; SMB depends on input field at reference elevation, plus vertical correction
       ! (2) SMB(x,y,z); SMB obtained by linear interpolation between values prescribed at adjacent vertical levels
       ! For options (1) and (2), the elevation-dependent SMB is computed here.
       ! Note: These conversion must be done at each time step (and not just once at initialization)
       !       if reading SMB or acab from a time-dependent forcing file.

       if (model%options%smb_input_function == SMB_INPUT_FUNCTION_XY) then

          if (model%options%smb_input == SMB_INPUT_MMYR_WE) then
             ! Convert units from mm/yr w.e. to m/yr ice, then convert to model units
             model%climate%acab(:,:) = (model%climate%smb(:,:) * (rhow/rhoi)/1000.d0) / scale_acab
          endif

       elseif (model%options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then

          if (model%options%smb_input == SMB_INPUT_MMYR_WE) then
             ! Convert smb_ref from mm/yr w.e. to model acab units.
             ! Convert smb_gradz to model acab_gradz units.
             ! Note: smb_gradz already includes a scale factor of thk0 in glide_vars.def
             model%climate%acab_ref(:,:) = (model%climate%smb_ref(:,:) * (rhow/rhoi)/1000.d0) / scale_acab
             model%climate%acab_gradz(:,:) = (model%climate%smb_gradz(:,:) * (rhow/rhoi)/1000.d0) / scale_acab
          endif

          ! compute acab by a lapse-rate correction to the reference value
          model%climate%acab(:,:) = model%climate%acab_ref(:,:) + &
               (model%geometry%usrf(:,:)*thk0 - model%climate%usrf_ref(:,:)) * model%climate%acab_gradz(:,:)

       elseif (model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then

          if (model%options%smb_input == SMB_INPUT_MMYR_WE) then
             ! Convert units from mm/yr w.e. to m/yr ice, then convert to model units
             model%climate%acab_3d(:,:,:) = (model%climate%smb_3d(:,:,:) * (rhow/rhoi)/1000.d0) / scale_acab
          endif

          ! downscale SMB to the local surface elevation (includes a halo update)
          ! Note: With linear_extrapolate_in = F, the values at top and bottom levels are simply extended upward and downward.
          !       For SMB, this is safer than linear extrapolation (especially when extrapolating upward)

          call glissade_vertical_interpolate(ewn,       nsn,                       &
                                             nlev_smb,  model%climate%smb_levels,  &
                                             model%geometry%usrf,                  &
                                             model%climate%acab_3d,                &
                                             model%climate%acab,                   &
                                             linear_extrapolate_in = .false.)

       endif   ! smb_input_function

       ! For the non-default smb_input_function options, make sure the SMB is nonzero somewhere; else abort.
       ! For the default option, do not abort, since idealized tests often have a zero SMB.

       call parallel_halo(model%climate%acab, parallel)

       if (model%options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ .or. &
           model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then
          local_maxval = maxval(abs(model%climate%acab))
          global_maxval = parallel_reduce_max(local_maxval)
          if (global_maxval < eps11) then
             write(message,*) 'Error: acab = 0 everywhere with smb_input_function =', &
                  model%options%smb_input_function
             call write_log(trim(message), GM_FATAL)
          endif
       endif

       if (verbose_smb .and. this_rank == rtest) then

          i = itest
          j = jtest
          write(6,*) 'Computing runtime acab with smb_input_function =', model%options%smb_input_function
          write(6,*) 'r, i, j =', this_rank, i, j
          write(6,*) ' '
          print*, 'usrf (m)'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%usrf(i,j) * thk0
             enddo
             write(6,*) ' '
          enddo

          if (model%options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
             write(6,*) ' '
             write(6,*) 'usrf - usrf_ref'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') &
                        (model%geometry%usrf(i,j)*thk0 - model%climate%usrf_ref(i,j))
                enddo
                write(6,*) ' '
             enddo
             write(6,*) ' '
             write(6,*) 'reference acab (m/yr ice)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%acab_ref(i,j) * scale_acab
                enddo
                write(6,*) ' '
             enddo
             write(6,*) ' '
             write(6,*) 'smb_gradz (mm/yr per m)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%smb_gradz(i,j) / thk0
                enddo
                write(6,*) ' '
             enddo
             write(6,*) ' '
             write(6,*) 'acab_gradz (m/yr per km)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%acab_gradz(i,j) * scale_acab/thk0 * 1000.d0
                enddo
                write(6,*) ' '
             enddo
             write(6,*) ' '
             write(6,*) 'adjusted acab (m/yr ice)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%acab(i,j) * scale_acab
                enddo
                write(6,*) ' '
             enddo

          elseif (model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then

             k = model%climate%nlev_smb/2 + 1  ! arbitrary k
             write(6,*) 'Diagnostic level k, usrf =', k, model%climate%smb_levels(k)*thk0
             write(6,*) ' '
             write(6,*) '3d acab(k) (m/yr ice)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%acab_3d(k,i,j) * scale_acab
                enddo
                write(6,*) ' '
             enddo
             write(6,*) ' '
             write(6,*) 'adjusted acab (m/yr ice)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%acab(i,j) * scale_acab
                enddo
                write(6,*) ' '
             enddo

          endif  ! smb_input_function

          if (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_GRADZ) then

             write(6,*) ' '
             write(6,*) 'reference artm (deg C)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%artm_ref(i,j)
                enddo
                write(6,*) ' '
             enddo
             write(6,*) ' '
             write(6,*) 'artm_gradz (deg C per km)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%artm_gradz(i,j)/thk0 * 1000.d0
                enddo
                write(6,*) ' '
             enddo
             write(6,*) ' '
             write(6,*) 'adjusted artm (deg C)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%artm(i,j)
                enddo
                write(6,*) ' '
             enddo

          elseif (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XYZ) then

             write(6,*) ' '
             write(6,*) '3d artm(k) (deg C)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%artm_3d(k,i,j)
                enddo
                write(6,*) ' '
             enddo
             write(6,*) ' '
             write(6,*) 'adjusted artm (deg C)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%climate%artm(i,j)
                enddo
                write(6,*) ' '
             enddo

          endif   ! artm_input_function

       endif  ! verbose_smb and this_rank

       ! If using a glacier-specific SMB index method, then compute the SMB and convert to acab

       if (model%options%enable_glaciers) then

          ! Halo updates for snow and artm
          ! (Not sure the artm update is needed; there is one above)
          call parallel_halo(model%climate%snow, parallel)
          call parallel_halo(model%climate%artm, parallel)

          call glissade_glacier_smb(&
               ewn,      nsn,                          &
               itest,    jtest,    rtest,              &
               model%glacier%nglacier,                 &
               model%glacier%cism_glacier_id,          &
               model%glacier%t_mlt,                    &  ! deg C
               model%climate%snow,                     &  ! mm/yr w.e.
               model%climate%artm,                     &  ! deg C
               model%glacier%mu_star,                  &  ! mm/yr w.e./deg
               model%climate%smb)                         ! mm/yr w.e.

          ! Convert SMB (mm/yr w.e.) to acab (CISM model units)
          model%climate%acab(:,:) = (model%climate%smb(:,:) * (rhow/rhoi)/1000.d0) / scale_acab
          call parallel_halo(model%climate%acab, parallel)

          if (verbose_glacier .and. this_rank == rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'Computed glacier SMB, rank, i, j =', this_rank, i, j
             print*, '   acab (m/yr ice) =', model%climate%acab(i,j)*thk0*scyr/tim0
          endif

       endif   ! enable_glaciers

       ! Compute a corrected acab field that includes any prescribed anomalies.
       ! Typically, acab_corrected = acab, but sometimes (e.g., for initMIP) it includes a time-dependent anomaly.
       ! Note that acab itself does not change in time.

       ! initialize
       model%climate%acab_corrected(:,:) = model%climate%acab(:,:)

       ! Optionally, multiply acab by a scalar adjustment factor
       if (model%climate%acab_factor /= 1.0d0) then
          model%climate%acab_corrected(:,:) = model%climate%acab_corrected(:,:) * model%climate%acab_factor
       endif

       if (model%options%enable_acab_anomaly) then

          if (model%options%smb_input == SMB_INPUT_MMYR_WE) then
             ! Convert units from mm/yr w.e. to m/yr ice, then convert to model units.
             !Note: If smb_anomaly is read from the input file, we could do this conversion once at initialization.
             !      But if read from a time-dependent forcing file, the conversion must be done repeatedly.
             model%climate%acab_anomaly(:,:) = (model%climate%smb_anomaly(:,:) * (rhow/rhoi) / 1000.d0) / scale_acab
          endif

          call parallel_halo(model%climate%acab_anomaly, parallel)

          ! Note: When being ramped up, the anomaly is not incremented until after the final time step of the year.
          !       This is the reason for passing the previous time to the subroutine.
          previous_time = model%numerics%time - model%numerics%dt * tim0/scyr

          call glissade_add_2d_anomaly(model%climate%acab_corrected,          &   ! scaled model units
                                       model%climate%acab_anomaly,            &   ! scaled model units
                                       model%climate%acab_anomaly_timescale,  &   ! yr
                                       previous_time)                             ! yr

          if (verbose_smb .and. this_rank==rtest) then
             i = itest
             j = jtest
             print*, 'i, j, previous_time, input acab, acab anomaly, corrected acab (m/yr):', &
                  i, j, previous_time, model%climate%acab(i,j)*thk0*scyr/tim0,  &
                  model%climate%acab_anomaly(i,j)*thk0*scyr/tim0,  &
                  model%climate%acab_corrected(i,j)*thk0*scyr/tim0
          endif

       endif

       ! Note: If enable_artm_anomaly = T, then artm_anomaly is computed above, before calling glissade_therm_driver.

       ! Optionally, overwrite acab_corrected where overwrite_acab_mask = 1.

       if (model%options%overwrite_acab /= 0) then

          call glissade_overwrite_acab(model%climate%overwrite_acab_mask,  &
                                       model%climate%overwrite_acab_value, &
                                       model%climate%acab_corrected)
       endif

       ! Optionally, block ice sheet inception by allowing SMB to be nonzero
       !  only in the main ice sheet and in cells adjacent to the ice sheet.
       ! Note: The ice sheet mask is computed in the diagnostic solve at the end of the previous time step
       !       (and at initialization).  Since it includes all cells that were part of the ice sheet
       !       before transport, an extended version of the mask will include all cells that potentially
       !       are part of the ice sheet after transport.

       if (model%options%block_inception) then

          ! Extend the ice sheet mask to include nearest neighbors (both edge and corner neighbors).

          call glissade_extend_mask(model%general%ewn,   model%general%nsn,     &
                                    model%geometry%ice_sheet_mask,              &
                                    extended_mask = extended_ice_sheet_mask)

          call parallel_halo(model%geometry%ice_sheet_mask, parallel)

          ! SMB is allowed to be positive only where ice_mask_smb = 1.
          ! Note: This code allows cells outside the ice sheet to melt.

          where (extended_ice_sheet_mask == 0)
             model%climate%acab_corrected = min(model%climate%acab_corrected, 0.0d0)
          endwhere

       endif

       ! Convert acab_corrected to a temporary array in SI units (m/s)
       acab_unscaled(:,:) = model%climate%acab_corrected(:,:) * thk0/tim0

       ! Convert bmlt to SI units (m/s)
       ! Note: bmlt is the sum of bmlt_ground (computed in glissade_thermal_solve) and bmlt_float
       !       (computed in glissade_bmlt_float_solve).
       ! Note: bmlt can be turned off by setting options%basal_mbal = BASAL_MBAL_NO_CONTINUITY

       if (model%options%basal_mbal == BASAL_MBAL_CONTINUITY) then    ! include bmlt in continuity equation
          bmlt_unscaled(:,:) = model%basal_melt%bmlt(:,:) * thk0/tim0
       else                                                           ! do not include bmlt in continuity equation
          bmlt_unscaled(:,:) = 0.0d0
       endif

       ! ------------------------------------------------------------------------
       ! Get masks used for the mass balance calculation.
       ! Pass thklim = 0 to identify cells with thck > 0 (not thck > thklim).
       ! Use ocean_mask to identify ocean cells where positive acab should not be applied.
       ! Use thck_calving_front to compute a fractional area for calving_front cells.
       ! TODO - Is it correct to use the old value of f_ground_cell from the start of the time step?
       !        Note that this value is used to identify CF cells where the mass balance is corrected.
       ! ------------------------------------------------------------------------

       call glissade_get_masks(ewn,              nsn,              &
                               parallel,                           &
                               thck_unscaled,                      &   ! m
                               topg_unscaled,                      &   ! m
                               model%climate%eus*thk0,             &   ! m
                               0.0d0,                              &   ! thklim = 0
                               ice_mask,                           &
                               floating_mask = floating_mask,      &
                               ocean_mask = ocean_mask,            &
                               land_mask = land_mask)

       ! Compute effective_areafrac for SMB purposes.

       call glissade_calving_front_mask(ewn,                    nsn,                &
                                        model%options%which_ho_calving_front,       &
                                        parallel,                                   &
                                        thck_unscaled,          topg_unscaled,      &   ! m
                                        model%climate%eus*thk0,                     &   ! m
                                        ice_mask,               floating_mask,      &
                                        ocean_mask,             land_mask,          &
                                        calving_front_mask,     thck_calving_front, &
                                        effective_areafrac = effective_areafrac)

       ! TODO: Zero out acab_unscaled and bmlt_unscaled in cells that are ice-free ocean after transport?
       !       Then it would not be necessary to pass ocean_mask to glissade_mass_balance_driver.

       ! Initialize the applied acab and bmlt.
       ! Note: These are smaller in magnitude than the input acab and bmlt for cells where either
       !       (1) the full column melts, and energy remains for melting, or
       !       (2) a positive mass balance is ignored, because a cell is ice-free ocean

       model%climate%acab_applied = 0.0d0
       model%basal_melt%bmlt_applied = 0.0d0

       ! ------------------------------------------------------------------------
       ! Apply the surface mass balance (acab) and basal mass balance (bmlt).
       ! Note: This subroutine assumes SI units:
       !       * dt (s)
       !       * dew, dns, thck (m)
       !       * acab, bmlt (m/s)
       ! ------------------------------------------------------------------------

       call glissade_mass_balance_driver(model%numerics%dt * tim0,                             &  ! s
                                         model%numerics%dew * len0, model%numerics%dns * len0, &  ! m
                                         ewn,         nsn,          upn-1,                     &
                                         model%numerics%sigma,                                 &
                                         parallel,                                             &
                                         thck_unscaled(:,:),                                   &  ! m
                                         acab_unscaled(:,:),                                   &  ! m/s
                                         bmlt_unscaled(:,:),                                   &  ! m/s
                                         model%climate%acab_applied(:,:),                      &  ! m/s
                                         model%basal_melt%bmlt_applied(:,:),                   &  ! m/s
                                         ocean_mask(:,:),                                      &
                                         effective_areafrac(:,:),                              &
                                         model%geometry%ntracers,                              &
                                         model%geometry%tracers(:,:,:,:),                      &
                                         model%geometry%tracers_usrf(:,:,:),                   &
                                         model%geometry%tracers_lsrf(:,:,:),                   &
                                         model%options%which_ho_vertical_remap)

       !-------------------------------------------------------------------------
       ! If running with glaciers, then adjust glacier indices based on advance and retreat.
       ! Note: This subroutine limits the ice thickness in grid cells that do not yet have
       !       a nonzero cism_glacier_id.  The acab_applied field is adjusted accordingly.
       ! Note: It would probably be OK to call this subroutine annually instead of every step.
       !       In that case, we might want to separate the special glacier acab adjustment
       !       from the rest of acab_applied.
       !-------------------------------------------------------------------------

       if (model%options%enable_glaciers) then

          call glissade_glacier_advance_retreat(&
               ewn,             nsn,                &
               itest,   jtest,  rtest,              &
               model%geometry%usrf*thk0,            &  ! m
               thck_unscaled,                       &  ! m
               model%climate%acab_applied,          &  ! m/s
               model%numerics%dt * tim0,            &  ! s
               model%glacier%minthck,               &  ! m
               model%glacier%cism_glacier_id_init,  &
               model%glacier%cism_glacier_id,       &
               parallel)

       endif   ! enable_glaciers

       !WHL - debug
       call parallel_halo(thck_unscaled, parallel)

       !-------------------------------------------------------------------------
       ! Cleanup
       !-------------------------------------------------------------------------

       ! copy tracers (temp/enthalpy, etc.) from model%geometry%tracers back to standard arrays
       call glissade_transport_finish_tracers(model)

       ! convert applied mass balance from m/s back to scaled model units
       model%climate%acab_applied(:,:) = model%climate%acab_applied(:,:)/thk0 * tim0
       model%basal_melt%bmlt_applied(:,:) = model%basal_melt%bmlt_applied(:,:)/thk0 * tim0

       ! convert thck back to scaled units
       ! (acab_unscaled is intent(in) above, so no need to scale it back)
       model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0

       ! For the enthalpy option, convert enthalpy back to temperature/waterfrac.

       if (model%options%whichtemp == TEMP_ENTHALPY) then

          ! Note: glissade_enth2temp expects SI units
          do j = 1, nsn
             do i = 1, ewn
                call glissade_enth2temp(model%numerics%stagsigma(1:upn-1),                                    &
                                        model%geometry%thck(i,j)*thk0,    model%temper%enthalpy(0:upn,i,j),   &
                                        model%temper%temp(0:upn,i,j),     model%temper%waterfrac(1:upn-1,i,j))
             enddo
          enddo
         
       endif    ! TEMP_ENTHALPY

       if (this_rank==rtest .and. verbose_glissade) then
          print*, ' '
          print*, 'After glissade_transport_driver:'
          print*, 'max, min thck (m)=', maxval(model%geometry%thck)*thk0, minval(model%geometry%thck)*thk0
          print*, 'max, min acab (m/yr) =', &
                  maxval(model%climate%acab_corrected)*scale_acab, &
                  minval(model%climate%acab_corrected)*scale_acab
          print*, 'thklim =', model%numerics%thklim * thk0
          print*, 'max, min temp =', maxval(model%temper%temp), minval(model%temper%temp)
          print*, ' '
          print*, 'thck:'
          write(6,'(a6)',advance='no') '      '
          do i = itest-5, itest+5
             write(6,'(i14)',advance='no') i
          enddo
          write(6,*) ' '
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%thck(i,j) * thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          k = upn
          print*, 'temp, k =', k
          write(6,'(a6)',advance='no') '      '
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') i
          enddo
          write(6,*) ' '
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%temper%temp(k,i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       call t_stopf('glissade_transport_driver')

       call t_stopf('inc_remap_driver')

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! restore old thickness
          model%geometry%thck(:,:) = model%geometry%thck_old(:,:)
       endif

    end select

    !------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    !TODO - Not sure this update is needed here.  It is done at the start
    !       of the diagnostic solve, but may not be needed for calving.
    !------------------------------------------------------------------------
    
    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       &
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'After mass balance, thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bmlt_applied (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%basal_melt%bmlt_applied(i,j)*(thk0/tim0)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif   ! verbose_inversion
    
  end subroutine glissade_thickness_tracer_solve

!=======================================================================

  subroutine glissade_calving_solve(model, init_calving)

    ! ------------------------------------------------------------------------ 
    ! Calculate iceberg calving
    ! ------------------------------------------------------------------------ 

    use cism_parallel, only: parallel_type, parallel_halo

    use glimmer_paramets, only: thk0, tim0, len0
    use glimmer_physcon, only: scyr
    use glissade_calving, only: glissade_calve_ice, glissade_cull_calving_front, &
         glissade_remove_icebergs, glissade_remove_isthmuses, glissade_limit_cliffs, verbose_calving
    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask,  &
         glissade_ocean_connection_mask
    use glissade_grounding_line, only: glissade_grounded_fraction
    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    logical, intent(in) :: init_calving  ! true when this subroutine is called at initialization

    ! --- Local variables ---

    real(dp), dimension(model%general%ewn, model%general%nsn) :: &
         thck_unscaled,           & ! model%geometry%thck converted to m
         thck_calving_front         ! effective ice thickness at calving front (m)

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ice_mask,                & ! = 1 if ice is present
         floating_mask,           & ! = 1 if ice is present and floating
         land_mask,               & ! = 1 if topg - eus >= 0
         ocean_mask,              & ! = 1 if ice is absent and topg - eus < 0
         active_ice_mask,         & ! = 1 if ice is present and dynamically active
         calving_front_mask         ! = 1 for calving-front cells

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ocean_connection_mask,   & ! = 1 for cells that are masked for retreat and are connected to the ocean
                                    ! through other cells that are masked for retreat
         retreat_mask               ! local version of ice_fraction_retreat_mask; excludes grounded cells

    real(dp) :: &
         maxthck,                 & ! max thickness of retreating ice
         dthck                      ! thickness loss for retreating ice

    logical :: cull_calving_front   ! true iff init_calving = T and options%cull_calving_front = T

    integer :: i, j

    integer :: nx, ny               ! horizontal grid dimensions
    integer :: itest, jtest, rtest  ! coordinates of diagnostic point

    real(dp), parameter :: &
         retreat_mask_threshold = 0.01d0  ! threshold value for removing cells based on ice_fraction_retreat_mask;
                                          !  set to a low value by default
                                          ! Could make this a config parameter

    ! variables to expand the calving mask at initialization
    logical, dimension(16) :: mask_basin  ! true for basins whose floating ice is added to the calving mask
                                          ! currently hardwired to 16 for ISMIP6
    integer :: bn    ! basin number

    type(parallel_type) :: parallel   ! info for parallel communication

    logical, parameter :: verbose_retreat = .true.

    nx = model%general%ewn
    ny = model%general%nsn

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    parallel = model%parallel

    ! Initialize

    model%calving%calving_thck = 0.0d0

    ! Thin or remove ice where retreat is forced.
    ! Note: This option is similar to apply_calving_mask.  It is different in that ice_fraction_retreat_mask
    !       is a real number in the range [0,1], allowing thinning instead of complete removal.
    !       Do not thin or remove ice if this is the initial calving call; force retreat only during runtime.
    ! There are two forced retreat options:
    ! Option 1: Thin or remove ice wherever ice_fraction_retreat_mask > 0.
    ! Option 2: Remove floating ice (but not grounded ice) where ice_fraction_retreat_mask > 0.
    !
    ! Option 1 is done now, before calling glissade_calve_ice, so that ice thinned by the retreat mask
    !        can undergo further thinning or removal by the calving scheme.
    ! Option 2 is done after the main calving solve, after thin ice at the calving front has been removed
    !  by other mechanisms.

    if (model%options%force_retreat == FORCE_RETREAT_ALL_ICE .and. .not.init_calving) then
       if (this_rank == rtest) then
          print*, 'Forcing retreat using ice_fraction_retreat_mask, time =', model%numerics%time
       endif

       if (verbose_retreat .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'Before forced retreat, thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'ice_fraction_retreat_mask:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%ice_fraction_retreat_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'maxthck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%reference_thck(i,j)*thk0 &
                     * (1.0d0 - model%geometry%ice_fraction_retreat_mask(i,j))
             enddo
             write(6,*) ' '
          enddo
       endif

       do j = 1, model%general%nsn
          do i = 1, model%general%ewn
             if (model%geometry%ice_fraction_retreat_mask(i,j) > 0.0d0) then
                maxthck = model%geometry%reference_thck(i,j) &
                     * (1.0d0 - model%geometry%ice_fraction_retreat_mask(i,j))
                dthck = model%geometry%thck(i,j) - min(maxthck, model%geometry%thck(i,j))
                model%geometry%thck(i,j) = model%geometry%thck(i,j) - dthck
                model%calving%calving_thck(i,j) = model%calving%calving_thck(i,j) + dthck
             endif
          enddo
       enddo

       if (verbose_retreat .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'After forced retreat, thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
       endif

    endif   ! force_retreat_all_ice

    !TODO - Make sure no additional halo updates are needed before glissade_calve_ice

    ! Set thickness to SI units (m)
    thck_unscaled(:,:) = model%geometry%thck(:,:)*thk0
    model%calving%calving_thck(:,:) = model%calving%calving_thck(:,:) * thk0

    ! Remove ice where forced by a calving mask.
    ! Note: whichcalving = CALVING_GRID_MASK and apply_calving_mask = T are currently redundant.
    ! TODO: Remove the CALVING_GRID_MASK option and use apply_calving_mask only (usually with marine_margin = 0).
    !       Keeping both for now to avoid breaking config files.

    if (model%options%whichcalving == CALVING_GRID_MASK .or. model%options%apply_calving_mask) then

       ! Optionally, expand the calving mask to include floating ice in select basins.
       ! Note: Currently hardwired to include 13 of the 16 ISMIP6 basins.
       !       Does not include the three largest shelves (Ross, Filchner-Ronne, Amery)

       call glissade_get_masks(nx,                       ny,                         &
                               parallel,                                             &
                               model%geometry%thck*thk0, model%geometry%topg*thk0,   &
                               model%climate%eus*thk0,   0.0d0,                      &  ! thklim = 0
                               ice_mask,                                             &
                               floating_mask = floating_mask,                        &
                               land_mask = land_mask)

       if (init_calving .and. model%options%expand_calving_mask) then

          ! Identify basins whose floating ice will be added to the calving mask
          ! Currently hardwired to the ISMIP6 basin numbers (1 to 16)
          mask_basin(:) = .true.
          mask_basin(2) = .false.   ! Amery
          mask_basin(7) = .false.   ! Ross
          mask_basin(14) = .false.  ! Filchner-Ronne

          if (verbose_calving .and. this_rank==rtest) then
             print*, 'Expanding the calving mask to ice shelves in select basins'
             print*, 'basin number, mask_basin:'
             do bn = 1, 16
                print*, bn, mask_basin(bn)
             enddo
          endif

          if (verbose_calving .and. this_rank==rtest) then
             print*, ' '
             print*, 'initial calving_mask, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(i10)',advance='no') model%calving%calving_mask(i,j)
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'floating_mask, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(i10)',advance='no') floating_mask(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

          ! For basins with mask_basin = T, add floating ice to the calving mask.
          do j = 1, model%general%nsn
             do i = 1, model%general%ewn
                bn = model%ocean_data%basin_number(i,j)
                if (mask_basin(bn) .and. floating_mask(i,j) == 1) then
                   model%calving%calving_mask(i,j) = 1
                endif
             enddo
          enddo

          call parallel_halo(model%calving%calving_mask, parallel)

       endif   ! expand_calving_mask

       if (verbose_calving .and. this_rank==rtest) then
          print*, ' '
          print*, 'Limit advance of calving front'
          print*, ' '
          print*, 'starting thck, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck_unscaled(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'floating_mask, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') floating_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'calving_mask, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') model%calving%calving_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
       endif

       ! calve ice where calving_mask = 1
       ! Optionally, if calving%timescale > 0, then there is a time scale for removal,
       !  allowing the CF to advance into masked regions.
       !TODO - Apply a time scale wherever calving%timescale > 0.
       !TODO - Move the mask logic to a subroutine.

       if (model%calving%timescale <= 1.0d0) then  ! currently have 1.0 yr in config files

          ! Remove ice in all cells with calving_mask = 1
          where (thck_unscaled > 0.0d0 .and. model%calving%calving_mask == 1)
             model%calving%calving_thck = model%calving%calving_thck + thck_unscaled
             thck_unscaled = 0.0d0
             !TODO - Reset temperature and other tracers in cells where the ice calved?
          endwhere

       else

          ! Thin the ice in floating cells where calving_mask = 1, based on a relaxation timescale

          ! In each masked floating cell, the thinning rate is max(H, H_c)/tau_c,
          !  where H_c is the calving thickness scale and tau_c the timescale.
          ! Thus the thinning rate is largest for thick ice.
          ! For thin ice, the rate has a minimum value H_c/tau_c..
          ! Note: calving%timescale has units of s (though input in yr in the config file)

          do j = 1, ny
             do i = 1, nx
                if (floating_mask(i,j) == 1 .and. model%calving%calving_mask(i,j) == 1) then
                   dthck = model%numerics%dt*tim0  &  ! dt in seconds
                        * max(thck_unscaled(i,j), model%calving%minthck) / model%calving%timescale
                   if (thck_unscaled(i,j) > dthck) then
                      model%calving%calving_thck(i,j) = model%calving%calving_thck(i,j) + dthck
                      thck_unscaled(i,j) = thck_unscaled(i,j) - dthck
                   else
                      model%calving%calving_thck(i,j) = model%calving%calving_thck(i,j) + thck_unscaled(i,j)
                      thck_unscaled(i,j) = 0.0d0
                   endif
                endif
             enddo   ! i
          enddo   ! j

          if (verbose_calving .and. this_rank==rtest) then
             print*, ' '
             print*, 'Relaxed calving, timescale (yr) =', model%calving%timescale/scyr
             print*, 'dt (yr) =', model%numerics%dt * tim0/scyr
             print*, 'calving_minthck (m) =', model%calving%minthck
             print*, ' '
             print*, 'calving_thck (m), itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%calving%calving_thck(i,j)
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'New thck (m):'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') thck_unscaled(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

      endif  ! relaxed calving

    elseif (model%options%which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then

       ! If using a subgrid calving_front scheme (but apply_calving_mask = F),
       !  remove thin ice that was transported beyond the CF to ice-free cells without active neighbors.
       ! In that case, a temporary version of model%calving%calving_mask is computed after transport and applied here.
       !TODO - Add a timescale

       where (model%calving%calving_mask == 1)
          model%calving%calving_thck = model%calving%calving_thck + thck_unscaled
          thck_unscaled = 0.0d0
          !TODO - Reset temperature and other tracers in cells where the ice calved?
       endwhere

    endif   ! apply_calving_mask

    ! ------------------------------------------------------------------------
    ! Calve ice, based on the value of whichcalving.
    ! Pass in thck, topg, etc. with units of meters.
    ! TODO: Pass in individual fields with SI units, instead of the calving derived type?
    ! ------------------------------------------------------------------------

    if (model%options%whichcalving /= CALVING_GRID_MASK) then

       call glissade_calve_ice(nx,           ny,                  &
                               model%options%whichcalving,        &
                               model%options%calving_domain,      &
                               model%options%which_ho_calving_front, &
                               parallel,                          &
                               model%calving,                     &        ! calving object; includes calving_thck (m)
                               itest, jtest, rtest,               &
                               model%numerics%dt*tim0,            &        ! s
                               model%numerics%dew*len0,           &        ! m
                               model%numerics%dns*len0,           &        ! m
                               model%numerics%sigma,              &
                               model%numerics%thklim*thk0,        &        ! m
                               thck_unscaled,                     &        ! m
                               model%isostasy%relx*thk0,          &        ! m
                               model%geometry%topg*thk0,          &        ! m
                               model%climate%eus*thk0)                     ! m
    endif

    if (init_calving .and. model%options%cull_calving_front) then

       call glissade_cull_calving_front(nx,           ny,              &
                                        parallel,                      &
                                        itest, jtest, rtest,           &
                                        thck_unscaled,                 &  ! m
                                        model%geometry%topg*thk0,      &  ! m
                                        model%climate%eus*thk0,        &  ! m
                                        model%numerics%thklim*thk0,    &  ! m
                                        model%options%which_ho_calving_front, &
                                        model%calving%ncull_calving_front,    &
                                        model%calving%calving_thck)       ! m

    endif

    if (model%options%force_retreat == FORCE_RETREAT_FLOATING_ICE) then

       ! Remove floating ice based on ice_fraction_retreat_mask.
       ! This is done after the main calving routine, to avoid complications
       !  involving thin ice near the calving front that calves after transport.
       ! The logic works as follows:
       ! * Idenfity cells with ice_fraction_retreat_mask exceeding some threshold.
       ! * Remove any such cells if they are adjacent to ocean cells, or are connected
       !   to the ocean through other identified cells.
       ! * Do not remove cells without a connection to the ocean.
       !   In other words, do not hollow out ice shelves from the interior, since
       !   this can be numerically unstable.

       ! Update masks
       call glissade_get_masks(nx,                     ny,                         &
                               parallel,                                           &
                               thck_unscaled,          model%geometry%topg*thk0,   &
                               model%climate%eus*thk0, model%numerics%thklim*thk0, &
                               ice_mask,                                           &
                               floating_mask = floating_mask,                      &
                               ocean_mask = ocean_mask)

       ! Identify floating cells with ice_fraction_retreat_mask exceeding a prescribed threshold.

       where (floating_mask == 1 .and. model%geometry%ice_fraction_retreat_mask > retreat_mask_threshold)
          retreat_mask = 1
       elsewhere
          retreat_mask = 0
       endwhere

       ! Identify cells that have retreat_mask = 1 and are either adjacent to ocean cells,
       !  or are connected to the ocean through other cells with retreat_mask = 1.

       call glissade_ocean_connection_mask(nx,            ny,           &
                                           parallel,                    &
                                           itest, jtest,  rtest,        &
                                           thck_unscaled, retreat_mask, &
                                           ocean_mask,                  &
                                           ocean_connection_mask)

       if (verbose_calving .and. this_rank==rtest) then
          print*, ' '
          print*, 'Force floating ice retreat'
          print*, ' '
          print*, 'Thickness before retreat, itest, jtest, rtest =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f7.1)',advance='no') thck_unscaled(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'floating_mask:'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i7)',advance='no') floating_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'ocean_mask:'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-7, itest+3
                write(6,'(i7)',advance='no') ocean_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'ice_fraction_retreat_mask:'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f7.2)',advance='no') model%geometry%ice_fraction_retreat_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'ocean_connection_mask:'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i7)',advance='no') ocean_connection_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       ! Remove ice from ocean-connected cells with retreat_mask = 1
       where (ocean_connection_mask == 1)
          model%calving%calving_thck = model%calving%calving_thck + thck_unscaled
          thck_unscaled = 0.0d0
          !TODO - Reset temperature and other tracers in cells where the ice calved?
       endwhere

    endif   ! force_retreat_floating_ice

    if (model%options%remove_isthmuses) then

       ! Optionally, remove isthmuses.
       ! An isthmus is defined as a floating or weakly grounded grid cell with ice-free ocean
       !  or thin floating ice on both sides.
       ! When using a calving or retreat mask derived from an ESM or other model,
       !  isthmuses may need to be removed to prevent unstable ice configurations,
       !  e.g. a shelf split into two parts connected by a bridge one cell wide.
       ! Isthmus removal should always be followed by iceberg removal.

       ! Update the masks
       call glissade_get_masks(nx,                     ny,                         &
                               parallel,                                           &
                               thck_unscaled,          model%geometry%topg*thk0,   &
                               model%climate%eus*thk0, model%numerics%thklim*thk0, &
                               ice_mask,                                           &
                               floating_mask = floating_mask,                      &
                               ocean_mask = ocean_mask,                            &
                               land_mask = land_mask)

       ! Compute f_ground_cell for isthmus removal

       call glissade_grounded_fraction(nx,          ny,               &
                                       parallel,                      &
                                       itest, jtest, rtest,           &  ! diagnostic only
                                       thck_unscaled,                 &
                                       model%geometry%topg*thk0,      &
                                       model%climate%eus*thk0,        &
                                       ice_mask,                      &
                                       floating_mask,                 &
                                       land_mask,                     &
                                       model%options%which_ho_ground, &
                                       model%options%which_ho_flotation_function, &
                                       model%options%which_ho_fground_no_glp,     &
                                       model%geometry%f_flotation,    &
                                       model%geometry%f_ground,       &
                                       model%geometry%f_ground_cell,  &
                                       model%geometry%topg_stdev*thk0)

       call glissade_remove_isthmuses(&
            nx,           ny,              &
            itest, jtest, rtest,           &
            thck_unscaled,                 &
            model%geometry%f_ground_cell,  &
            floating_mask,                 &
            ocean_mask,                    &
            model%calving%calving_thck)

    endif  ! remove isthmuses

    ! ------------------------------------------------------------------------
    ! Remove any icebergs.
    ! For the velocity solver to be robust, we require that any floating cell
    !  is connected to grounded ice along a path consisting only of active cells.
    ! Floating cells without such a connection are calved as icebergs.
    ! ------------------------------------------------------------------------

    if (model%options%remove_icebergs) then

       ! Update the basic masks

       call glissade_get_masks(nx,                     ny,                            &
                               parallel,                                              &
                               thck_unscaled,          model%geometry%topg*thk0,      &
                               model%climate%eus*thk0, model%numerics%thklim*thk0,    &
                               ice_mask,               floating_mask = floating_mask, &
                               land_mask = land_mask,  ocean_mask = ocean_mask,       &
                               active_ice_mask = active_ice_mask)

       ! If using a subgrid CF scheme, then recompute the active ice mask
       if (model%options%which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then

          call glissade_calving_front_mask(nx,                     ny,                 &
                                           model%options%which_ho_calving_front,       &
                                           parallel,                                   &
                                           thck_unscaled,                              &
                                           model%geometry%topg*thk0,                   &
                                           model%climate%eus*thk0,                     &
                                           ice_mask,               floating_mask,      &
                                           ocean_mask,             land_mask,          &
                                           calving_front_mask,     thck_calving_front, &
                                           active_ice_mask = active_ice_mask)

          if (verbose_calving .and. this_rank == rtest) then
             print*, ' '
             print*, 'Before remove_icebergs, calving_front_mask, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(i10)',advance='no') calving_front_mask(i,j)
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'thck_calving_front:'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') thck_calving_front(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

       endif   ! subgrid calving front

       ! Compute the grounded ice fraction in each grid cell

       call glissade_grounded_fraction(nx,          ny,               &
                                       parallel,                      &
                                       itest, jtest, rtest,           &  ! diagnostic only
                                       thck_unscaled,                 &
                                       model%geometry%topg*thk0,      &
                                       model%climate%eus*thk0,        &
                                       ice_mask,                      &
                                       floating_mask,                 &
                                       land_mask,                     &
                                       model%options%which_ho_ground, &
                                       model%options%which_ho_flotation_function, &
                                       model%options%which_ho_fground_no_glp,     &
                                       model%geometry%f_flotation,    &
                                       model%geometry%f_ground,       &
                                       model%geometry%f_ground_cell,  &
                                       model%geometry%topg_stdev*thk0)

       ! Remove icebergs.
       ! Icebergs are defined as floating cells that do not have a path through active cells
       !  to grounded cells (i.e., cells where f_ground_cell exceeds a threshold value).

       call glissade_remove_icebergs(nx,           ny,                     &
                                     parallel,                             &
                                     itest, jtest, rtest,                  &
                                     thck_unscaled,                        &  ! m
                                     model%geometry%f_ground_cell,         &
                                     ice_mask,                             &
                                     land_mask,                            &
                                     active_ice_mask,                      &
                                     model%calving%calving_thck)              ! m
    endif
    
    ! Optionally, impose a thickness limit on marine ice cliffs.
    ! These are defined as grounded marine-based cells adjacent to inactive calving_front cells or ice-free ocean.

    if (model%options%limit_marine_cliffs) then   ! Impose a thickness limit on marine ice cliffs

       call glissade_limit_cliffs(nx,             ny,            &
                                  parallel,                      &
                                  itest,  jtest,  rtest,         &
                                  model%numerics%dt*tim0,        &     ! s
                                  model%options%which_ho_calving_front, &
                                  model%calving%taumax_cliff,    &     ! Pa
                                  model%calving%cliff_timescale, &     ! s
                                  thck_unscaled,              &        ! m
                                  model%geometry%topg*thk0,   &        ! m
                                  model%climate%eus*thk0,     &        ! m
                                  model%numerics%thklim*thk0, &        ! m
                                  model%calving%calving_thck)          ! m
    endif

    ! Convert geometry%thck and calving%calving_thck to scaled model units
    model%geometry%thck(:,:) = thck_unscaled(:,:)/thk0
    model%calving%calving_thck(:,:) = model%calving%calving_thck(:,:)/thk0

    !TODO: Are any other halo updates needed after calving?
    ! halo updates
    call parallel_halo(model%geometry%thck, parallel)   ! Updated halo values of thck are needed below in calclsrf

    ! update the upper and lower surfaces

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       &
                        model%climate%eus,   model%geometry%lsrf)
    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

    if (verbose_calving .and. this_rank == rtest) then
       print*, ' '
       print*, 'Final calving_thck (m), itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%calving%calving_thck(i,j) * thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Final thck (m):'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%thck(i,j) * thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Final usrf (m):'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%usrf(i,j) * thk0
          enddo
          write(6,*) ' '
       enddo
    endif  ! verbose_calving

  end subroutine glissade_calving_solve

!=======================================================================

  subroutine glissade_isostasy_solve(model)

    ! ------------------------------------------------------------------------ 
    ! Calculate isostatic adjustment
    ! ------------------------------------------------------------------------ 

    use cism_parallel, only: parallel_type, parallel_halo, parallel_halo_extrapolate

    use isostasy, only: isos_compute, isos_icewaterload
    use glimmer_paramets, only: thk0
    use glissade_masks, only: glissade_marine_connection_mask

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! --- Local variables ---

    type(parallel_type) :: parallel   ! info for parallel communication

    parallel = model%parallel

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! Note: Suppose the update period is 100 years, and the time step is 1 year.
    !       Then the update will be done on the first time step of the simulation,
    !        (model%numerics%tstep_count = 1) and again on step 101, 201, etc.
    !       The update will not be done before writing output at t = 100, when
    !        model%numerics%tstep_count = 100.
    !       Thus the output file will contain the load that was applied during the
    !        preceding years, not the new load.
    !       In older code versions, the new load would have been computed on step 100.
    ! ------------------------------------------------------------------------

    if (model%options%isostasy == ISOSTASY_COMPUTE) then

       if (model%isostasy%nlith > 0) then
          if (mod(model%numerics%tstep_count-1, model%isostasy%nlith) == 0) then
             if (main_task) then
                print*, 'Update lithospheric load: tstep_count, nlith =', &
                     model%numerics%tstep_count, model%isostasy%nlith
             endif
             call isos_icewaterload(model)
             model%isostasy%new_load = .true.
          end if
       endif  ! nlith > 0

    end if
   
    ! ------------------------------------------------------------------------ 
    ! Calculate isostatic adjustment
    ! ------------------------------------------------------------------------ 

    if (model%options%isostasy == ISOSTASY_COMPUTE) then

       call isos_compute(model)

       ! update topography in halo cells
       ! Note: For outflow BCs, most fields (thck, usrf, temp, etc.) are set to zero in the global halo,
       !        to create ice-free conditions. However, we might not want to set topg = 0 in the global halo,
       !        because then the global halo will be interpreted as ice-free land, whereas we may prefer to
       !        treat it as ice-free ocean. For this reason, topg is extrapolated from adjacent cells.
       !       Similarly, for no_ice BCs, we want to zero out ice state variables adjacent to the global boundary,
       !        but we do not want to zero out the topography.
       ! Note: The topg halo update at initialization has an optional argument periodic_ew,
       !        which is needed for ismip-hom. I doubt ismip-hom will be run with active isostasy,
       !        but the argument is included to be on the safe side.
       ! TODO: Do we need similar logic for halo updates of relx?

       if (model%general%global_bc == GLOBAL_BC_OUTFLOW .or. &
           model%general%global_bc == GLOBAL_BC_NO_ICE) then
          call parallel_halo_extrapolate(model%geometry%topg, parallel)
       else  ! other global BCs, including periodic
          call parallel_halo(model%geometry%topg, parallel, &
                             periodic_offset_ew = model%numerics%periodic_offset_ew)
       endif

       ! update the marine connection mask, which depends on topg

       call glissade_marine_connection_mask(&
            model%general%ewn,          model%general%nsn,          &
            parallel,                                               &
            model%numerics%idiag_local, model%numerics%jdiag_local, &
            model%numerics%rdiag_local,                             &
            model%geometry%thck * thk0, model%geometry%topg * thk0, &
            model%climate%eus * thk0,   0.0d0,                      &  ! thklim = 0
            model%geometry%marine_connection_mask)

    end if

  end subroutine glissade_isostasy_solve

!=======================================================================

  subroutine glissade_diagnostic_variable_solve(model) 

     ! Solve diagnostic (not time-dependent) variables, in particular the ice velocity.
     ! This is needed at the end of each time step once the prognostic variables (thickness, tracers) have been updated.  
     ! It is also needed to fill out the initial state from the fields that have been read in.

    use cism_parallel, only: parallel_type, parallel_halo, &
         staggered_parallel_halo, staggered_parallel_halo_extrapolate, &
         parallel_reduce_max, parallel_reduce_min, parallel_globalindex

    use glimmer_paramets, only: eps08, tim0, len0, vel0, thk0, vis0, tau0, evs0
    use glimmer_physcon, only: rhow, rhoi, scyr
    use glimmer_scales, only: scale_acab
    use glide_thck, only: glide_calclsrf
    use glissade_velo, only: glissade_velo_driver
    use glide_velo, only: wvelintg
    use glissade_masks, only: glissade_get_masks, glissade_ice_sheet_mask, glissade_calving_front_mask
    use glissade_grid_operators, only: glissade_stagger, glissade_gradient
    use glissade_grounding_line, only: glissade_grounded_fraction, glissade_grounding_line_flux, verbose_glp
    use glissade_therm, only: glissade_interior_dissipation_sia,  &
                              glissade_interior_dissipation_first_order, &
                              glissade_flow_factor,  &
                              glissade_pressure_melting_point
    use glissade_calving, only: verbose_calving
    use felix_dycore_interface, only: felix_velo_driver
    use glissade_basal_traction, only: calc_effective_pressure
    use glissade_inversion, only: glissade_inversion_basal_friction,  &
         glissade_inversion_bmlt_basin, verbose_inversion
    use glissade_glacier, only: glissade_glacier_inversion

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

    integer :: i, j, k, n, ng
    integer :: itest, jtest, rtest

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ice_mask,           & ! = 1 where thck > thklim, else = 0
         floating_mask,      & ! = 1 where ice is present and floating, else = 0
         ocean_mask,         & ! = 1 where topg is below sea level and ice is absent
         land_mask,          & ! = 1 where topg is at or above sea level
         calving_front_mask, & ! = 1 where ice is floating and borders an ocean cell, else = 0
         marine_interior_mask  ! = 1 if ice is marine-based and borders no ocean cells, else = 0

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck_calving_front      ! effective thickness of ice at the calving front

    real(dp) :: &
         dsigma,                   & ! layer thickness in sigma coordinates
         tau_xx, tau_yy, tau_xy,   & ! stress tensor components
         eps_xx, eps_yy, eps_xy      ! strain rate tensor components
 
    real(dp) :: &
         a, b, c, root,   & ! terms in quadratic formula
         lambda1, lambda2   ! eigenvalues of horizontal strain rate tensor

    real(dp) :: nudging_factor      ! factor in range [0,1], used for inversion of powerlaw_c
    real(dp) :: weaning_time        ! time since the start of weaning (numerics%time - inversion%wean_tstart)

    type(parallel_type) :: parallel   ! info for parallel communication

    integer :: ewn, nsn, upn

    !WHL - debug
    integer :: iglobal, jglobal, ii, jj
    real(dp) :: my_max, my_min, global_max, global_min
    real(dp) :: sum_cell, sum1, sum2  ! temporary sums

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         floating_mask_old, grounded_mask_old   ! masks from previous time steps

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    parallel = model%parallel

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn

    if (verbose_glissade .and. main_task) then
       print*, 'In glissade_diagnostic_variable_solve'
    endif

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 1. First part of diagnostic solve: 
    !    Now that advection and calving are done, update geometry- and temperature-related
    !    diagnostic fields that are needed for the velocity solve.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! ------------------------------------------------------------------------
    ! Halo update for ice thickness
    ! Note: The halo update for topg is done at initialization and again (as needed)
    !       after computing isostasy.
    ! ------------------------------------------------------------------------

    call parallel_halo(model%geometry%thck, parallel)

    ! ------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    ! ------------------------------------------------------------------------
    !TODO - These are currently updated after transport. Needed for calving/isostasy, or not until here?

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       & 
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

    ! ------------------------------------------------------------------------
    ! Compute some quantities on the staggered grid.
    ! These are not used by the Glissade velocity solver, but are used for diagnostics
    !  (and optionally for SIA-based dissipation).
    ! ------------------------------------------------------------------------

    call glissade_stagger(ewn,                 nsn,  &
                          model%geometry%thck, model%geomderv%stagthck)

    call glissade_gradient(ewn,                     nsn,       &
                           model%numerics%dew,      model%numerics%dns,      &
                           model%geometry%usrf,                              &
                           model%geomderv%dusrfdew, model%geomderv%dusrfdns)

    ! ------------------------------------------------------------------------
    !TODO - Move this calculation to the calving solver?  Apply to a different flux, instead of calving?
    ! Compute masks for the ice sheet and ice caps.
    ! Ice caps are defined as ice-covered cells disconnected from the main ice sheet.
    ! Optionally, the ice sheet mask can be used to block inception outside the existing ice sheet.
    ! ------------------------------------------------------------------------

    call glissade_get_masks(ewn,                 nsn,     &
                            parallel,                                   &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   model%numerics%thklim, &
                            ice_mask)

    call glissade_ice_sheet_mask(ewn,      nsn, &
                                 parallel,                      &
                                 itest,    jtest,   rtest,      &
                                 ice_mask,                      &
                                 model%geometry%thck*thk0,      &
                                 model%geometry%ice_sheet_mask, &
                                 model%geometry%ice_cap_mask)

    call parallel_halo(model%geometry%ice_sheet_mask, parallel)
    call parallel_halo(model%geometry%ice_cap_mask, parallel)

    if (model%options%remove_ice_caps) then

       ! Remove ice caps and add them to the calving flux.
       ! If ice caps are absent in the input file, and SMB = 0 over all cells
       !  separated from the main sheet, then ice caps may never form.
       ! However, it is possible that the main ice sheet will advance under a positive SMB,
       !  and then part of that ice will melt under a negative SMB, leaving a remnant ice cap.
       ! Such remnant ice caps could flow, possibly joining the main ice sheet.
       ! Note: The ice cap mask is not updated after removal.  So if this mask is written to output,
       !       it will show where ice caps existed before they were removed.

       where (model%geometry%ice_cap_mask == 1)
          model%calving%calving_thck = model%calving%calving_thck + model%geometry%thck
          model%geometry%thck = 0.0d0
       endwhere

    endif   ! remove_ice_caps

    ! ------------------------------------------------------------------------
    ! Update some masks that are used for subsequent calculations
    ! ------------------------------------------------------------------------

    call glissade_get_masks(ewn,                 nsn,                   &
                            parallel,                                   &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   model%numerics%thklim, &
                            ice_mask,                                   &
                            floating_mask = floating_mask,              &
                            ocean_mask = ocean_mask,                    &
                            land_mask = land_mask)

    call glissade_calving_front_mask(ewn,                 nsn,     &
                                     model%options%which_ho_calving_front,       &
                                     parallel,                                   &
                                     model%geometry%thck, model%geometry%topg,   &
                                     model%climate%eus,                          &
                                     ice_mask,            floating_mask,         &
                                     ocean_mask,          land_mask,             &
                                     calving_front_mask,  thck_calving_front,    &
                                     marine_interior_mask = marine_interior_mask)

    ! ------------------------------------------------------------------------
    ! Compute the fraction of grounded ice in each cell and at each vertex.
    ! The grounded fraction at each vertex, f_ground, is used in the velocity solver
    !  to adjust the basal friction.
    ! The grounded fraction in each cell, f_ground_cell, is optionally used
    !  to adjust the melt rate beneath floating ice, and also can be used
    !  to adjust the flow enhancement factor.
    ! Note: This subroutine requires that thck and topg are up to date in halo cells.
    !
    ! See comments in subroutine glissade_grounded_fraction for details
    ! on the whichground and whichflotation_function options.
    !
    ! Computing f_ground here ensures that it is always available as a diagnostic
    ! at the end of the time step, even if the velocity solver is not called
    ! (e.g., on the first time step of a restart).
    ! ------------------------------------------------------------------------

    call glissade_grounded_fraction(ewn,          nsn,             &
                                    parallel,                      &
                                    itest, jtest, rtest,           &  ! diagnostic only
                                    model%geometry%thck*thk0,      &
                                    model%geometry%topg*thk0,      &
                                    model%climate%eus*thk0,        &
                                    ice_mask,                      &
                                    floating_mask,                 &
                                    land_mask,                     &
                                    model%options%which_ho_ground, &
                                    model%options%which_ho_flotation_function, &
                                    model%options%which_ho_fground_no_glp,     &
                                    model%geometry%f_flotation,    &
                                    model%geometry%f_ground,       &
                                    model%geometry%f_ground_cell,  &
                                    model%geometry%topg_stdev*thk0)

    !WHL - debug
    if (this_rank == rtest .and. verbose_glp) then
       print*, ' '
       print*, 'Called GLP subroutine, which_ho_ground =', model%options%which_ho_ground
       print*, ' '
       print*, 'f_flotation:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%f_flotation(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'f_ground at vertex:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.5)',advance='no') model%geometry%f_ground(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'f_ground_cell:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.5)',advance='no') model%geometry%f_ground_cell(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
    endif  ! this_rank = rtest

    ! Compute the thickness tendency dH/dt from one step to the next (m/s)
    ! This tendency is used for coulomb_c and powerlaw_c inversion.

    if ( (model%options%is_restart == RESTART_TRUE) .and. &
         (model%numerics%time == model%numerics%tstart) ) then
       ! first call after a restart; do not compute dthck_dt
    else
       model%geometry%dthck_dt(:,:) = (model%geometry%thck(:,:) - model%geometry%thck_old(:,:)) * thk0 &
                                    / (model%numerics%dt * tim0)
    endif

    ! If inverting for Cp = powerlaw_c or Cc = coulomb_c, then update it here.
    ! Note: This subroutine used to be called earlier, but now is called here
    !       in order to have f_ground_cell up to date.
    ! If running with glaciers, inversion for powerlaw_c is done elsewhere,
    !  in subroutine glissade_glacier_inversion.
    !TODO: Call when the inversion options are set, not the external options.
    !      Currently, the only thing done for the external options is to remove
    !       zero values.

    if ( model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or. &
         model%options%which_ho_powerlaw_c == HO_POWERLAW_C_EXTERNAL  .or. &
         model%options%which_ho_coulomb_c  == HO_COULOMB_C_INVERSION  .or. &
         model%options%which_ho_coulomb_c  == HO_COULOMB_C_EXTERNAL ) then

       if ( (model%options%is_restart == RESTART_TRUE) .and. &
            (model%numerics%time == model%numerics%tstart) ) then
          ! first call after a restart; do not update powerlaw_c or coulomb_c
       else
          call glissade_inversion_basal_friction(model)
       endif

    endif   ! which_ho_powerlaw_c/coulomb_c

    ! If inverting for deltaT_basin, then update it here

    if ( model%options%which_ho_bmlt_basin_inversion == HO_BMLT_BASIN_INVERSION_COMPUTE) then

       if ( (model%options%is_restart == RESTART_TRUE) .and. &
            (model%numerics%time == model%numerics%tstart) ) then
          ! first call after a restart; do not update basin-scale melting parameters

       else

          call glissade_inversion_bmlt_basin(model%numerics%dt * tim0,                  &  ! s
                                             ewn, nsn,                                  &
                                             model%numerics%dew * len0,                 &  ! m
                                             model%numerics%dns * len0,                 &  ! m
                                             itest, jtest, rtest,                       &
                                             model%ocean_data%nbasin,                   &
                                             model%ocean_data%basin_number,             &
                                             model%geometry%thck*thk0,                  &  ! m
                                             model%geometry%dthck_dt,                   &  ! m/s
                                             model%inversion%floating_thck_target*thk0, &  ! m
                                             model%inversion%dbmlt_dtemp_scale,         &  ! (m/s)/degC
                                             model%inversion%bmlt_basin_timescale,      &  ! s
                                             model%ocean_data%deltaT_basin,             &
                                             model%inversion%bmlt_basin_mass_correction,&
                                             model%inversion%bmlt_basin_number_mass_correction)

       endif  ! first call after a restart

    endif   ! which_ho_bmlt_basin_inversion

    ! If glaciers are enabled, invert for mu_star and powerlaw_c based on area and volume targets

    if (model%options%enable_glaciers .and. &
         (model%glacier%set_mu_star == GLACIER_MU_STAR_INVERSION .or.  &
          model%glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION)) then

       if (model%numerics%time == model%numerics%tstart) then

           ! first call at start-up or after a restart; do not invert

       else

          call glissade_glacier_inversion(model, model%glacier)

       endif   ! time = tstart

    endif   ! enable_glaciers with inversion

    ! ------------------------------------------------------------------------ 
    ! Calculate Glen's A
    !
    ! Notes:
    ! (1) Because flwa is not a restart variable in Glissade, no check is included 
    !      here for whether to calculate it on initial time (as is done in Glide).
    ! (2) We are passing in only vertical elements (1:upn-1) of the temp array,
    !       so that it has the same vertical dimensions as flwa.
    ! (3) The flow enhancement factor is 1 by default.
    ! (4) The waterfrac field is ignored unless whichtemp = TEMP_ENTHALPY.
    ! (5) Inputs and outputs of glissade_flow_factor should have SI units.
    ! ------------------------------------------------------------------------

    call glissade_flow_factor(model%options%whichflwa,            &
                              model%options%whichtemp,            &
                              model%numerics%stagsigma,           &
                              model%geometry%thck * thk0,         &  ! scale to m
                              model%temper%temp(1:upn-1,:,:),     &
                              model%temper%flwa,                  &  ! Pa^{-n} s^{-1}
                              model%paramets%default_flwa / scyr, &  ! scale to Pa^{-n} s^{-1}
                              model%paramets%flow_enhancement_factor,       &
                              model%paramets%flow_enhancement_factor_float, &
                              model%options%which_ho_ground,      &
                              floating_mask,                      &
                              model%geometry%f_ground_cell,       &
                              model%temper%waterfrac)

    !TODO - flwa halo update not needed?
    ! Halo update for flwa
    call parallel_halo(model%temper%flwa, parallel)

    ! ------------------------------------------------------------------------
    ! Do some additional operations if this is the first time step.
    ! The model thickness and temperature fields will have been initialized, but the
    !  thermal and transport solvers have not been called yet.
    ! Note: Some operations must be done in this subroutine when restarting; others are skipped.
    ! ------------------------------------------------------------------------

    if (model%numerics%time == model%numerics%tstart) then

       ! Compute the pressure melting point temperature, which is needed
       ! by certain basal sliding laws.

       do j = 1, nsn
          do i = 1, ewn
             call glissade_pressure_melting_point(model%geometry%thck(i,j) * thk0, &
                                                  model%temper%bpmp(i,j))
          enddo
       enddo

       ! If the velocity fields have been read in on the extended staggered mesh,
       ! then copy them to the standard staggered mesh.
       !
       ! Note: For problems with nonzero velocity along the global boundaries (e.g., MISMIP on a periodic domain),
       !        exact restart requires that the restart velocity field lies on an extended staggered mesh with
       !        an extra row and column of velocity points along the north and east boundary of the domain.
       !       In that case, uvel_extend and vvel_extend should be written to the restart file by setting
       !        restart_extend_velo = 1 in the config file. They will then be read as input fields and at this
       !        point need to be copied into uvel and vvel.
       !       (It would have been cleaner to give uvel and vvel the same dimensions as the scalar mesh,
       !        but that design decision was made many years ago and would take a lot of work to change.)

       ! If uvel_extend and vvel_extend are input fields, then copy them into uvel and vvel.
       ! Halo updates are then needed to make sure the velocities are correct along the boundaries.
 
       if  ( (maxval(abs(model%velocity%uvel_extend)) /= 0.0d0) .or. & 
             (maxval(abs(model%velocity%vvel_extend)) /= 0.0d0) ) then
          call write_log('Using uvel_extend, vvel_extend from input or restart file at initial time')
          model%velocity%uvel(:,:,:) = model%velocity%uvel_extend(:,1:ewn-1,1:nsn-1)
          model%velocity%vvel(:,:,:) = model%velocity%vvel_extend(:,1:ewn-1,1:nsn-1)
!       elseif ( (maxval(abs(model%velocity%uvel)) /= 0.0d0) .or. & 
!                (maxval(abs(model%velocity%vvel)) /= 0.0d0) ) then
!          call write_log('Using uvel, vvel from input or restart file at initial time')
       endif

       call staggered_parallel_halo(model%velocity%uvel, parallel)
       call staggered_parallel_halo(model%velocity%vvel, parallel)

       ! The DIVA solver option requires some additional fields on the staggered mesh for exact restart.
       ! If these fields were input on the extended mesh, they need to be copied to the standard staggered mesh.
       ! Halo updates are then needed to make sure they have the correct values along the boundaries.

       if (model%options%which_ho_approx == HO_APPROX_DIVA) then

          if  ( (maxval(abs(model%velocity%uvel_2d_extend)) /= 0.0d0) .or. & 
                (maxval(abs(model%velocity%vvel_2d_extend)) /= 0.0d0) ) then
             call write_log('Using uvel_2d_extend, vvel_2d_extend from input or restart file at initial time')
             model%velocity%uvel_2d(:,:) = model%velocity%uvel_2d_extend(1:ewn-1,1:nsn-1)
             model%velocity%vvel_2d(:,:) = model%velocity%vvel_2d_extend(1:ewn-1,1:nsn-1)
!          elseif ( (maxval(abs(model%velocity%uvel_2d)) /= 0.0d0) .or. & 
!                   (maxval(abs(model%velocity%vvel_2d)) /= 0.0d0) ) then
!             call write_log('Using uvel_2d, vvel_2d from input or restart file at initial time')
          endif

          if  ( (maxval(abs(model%stress%btractx_extend)) /= 0.0d0) .or. & 
                (maxval(abs(model%stress%btracty_extend)) /= 0.0d0) ) then
             model%stress%btractx(:,:) = model%stress%btractx_extend(1:ewn-1,1:nsn-1)
             model%stress%btracty(:,:) = model%stress%btracty_extend(1:ewn-1,1:nsn-1)
          endif

          call staggered_parallel_halo(model%velocity%uvel_2d, parallel)
          call staggered_parallel_halo(model%velocity%vvel_2d, parallel)
          call staggered_parallel_halo(model%stress%btractx, parallel)
          call staggered_parallel_halo(model%stress%btracty, parallel)

       endif   ! DIVA approx
             
    endif   ! time = tstart

    !------------------------------------------------------------------------------
    ! Compute the effective pressure N at the bed.
    ! Although N is not needed for all sliding options, it is computed here just in case.
    ! Note: effective pressure is part of the basal_physics derived type.
    ! Note: Ideally, bpmp and temp(nz) are computed after the transport solve,
    !       just before the velocity solve. Then they will be consistent with the
    !       current thickness field.
    !------------------------------------------------------------------------------

    !TODO - Use btemp_ground instead of temp(nz)?
    call calc_effective_pressure(model%options%which_ho_effecpress, &
                                 parallel,                          &
                                 ewn,           nsn,                &
                                 model%basal_physics,               &
                                 model%basal_hydro,                 &
                                 ice_mask,      floating_mask,      &
                                 model%geometry%thck * thk0,        &
                                 model%geometry%topg * thk0,        &
                                 model%climate%eus * thk0,          &
                                 model%temper%bpmp(:,:) - model%temper%temp(upn,:,:), &
                                 model%basal_hydro%bwat * thk0,     &   ! m
                                 model%basal_hydro%bwatflx,         &   ! m/yr
                                 model%numerics%dt * tim0/scyr,     &   ! yr
                                 itest, jtest,  rtest)

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 2. Second part of diagnostic solve: 
    !    Now that geometry- and temperature-related diagnostic fields are updated, 
    !    solve velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! Do not solve velocity for initial time on a restart because that breaks an exact restart.
    ! Note: model%numerics%tstart is the time of restart, not necessarily the value of tstart in the config file.

    if ( (model%options%is_restart == RESTART_TRUE) .and. &
         (model%numerics%time == model%numerics%tstart) ) then
  
       ! Do not solve for velocity, because this would break exact restart

    else

       ! If this is not a restart or we are not at the initial time, then proceed normally

       if ( (model%numerics%time == model%numerics%tstart) .and. &
         ( (maxval(abs(model%velocity%uvel)) /= 0.0d0) .or. & 
           (maxval(abs(model%velocity%vvel)) /= 0.0d0) ) ) then
          ! If velocity was input and this is NOT a restart, then use the input field as the first guess at the initial time.
          ! This happens automatically, but let the user know.
          ! Using this value will change the answer only within the tolerance of the nonlinear solve.  
          ! If a user already has a good guess from a previous run, they may wish to start things off with it to speed the initial solution.
          call write_log('Using uvel, vvel from input file as initial guess at initial time.')
          call write_log('If this is not desired, please remove those fields from the input file.')
       endif

       if (main_task) then
          print *, ' '
          print *, 'Compute ice velocities, time =', model%numerics%time
       endif

       !! extrapolate value of mintauf into halos to enforce periodic lateral bcs (only if field covers entire domain)
       if (model%options%which_ho_babc == HO_BABC_YIELD_PICARD) then
          call staggered_parallel_halo_extrapolate(model%basal_physics%mintauf, parallel)
       endif

       ! Call the appropriate velocity solver

       select case (model%options%whichdycore)

       case ( DYCORE_GLISSADE )   ! glissade finite-element

          call t_startf('glissade_velo_driver')
          call glissade_velo_driver(model)
          call t_stopf('glissade_velo_driver')

       case ( DYCORE_ALBANYFELIX)

          call t_startf('felix_velo_driver')
          call felix_velo_driver(model)
          call t_stopf('felix_velo_driver')

       end select

       ! Compute internal heat dissipation
       ! This is used in the prognostic temperature calculation during the next time step.
       ! Note: These glissade subroutines assume SI units on input and output

       model%temper%dissip(:,:,:) = 0.d0

       if (model%options%which_ho_disp == HO_DISP_SIA) then

          ! Compute dissipation based on the shallow-ice approximation

          call glissade_interior_dissipation_sia(ewn,  nsn,    upn,              &
                                                 model%numerics%stagsigma(:),    &
                                                 ice_mask,                       &
                                                 model%geomderv%stagthck * thk0, & ! scale to m
                                                 model%temper%flwa * vis0,       & ! scale to Pa^{-n} s^{-1}
                                                 model%geomderv%dusrfdew * thk0/len0, & ! scale to m/m
                                                 model%geomderv%dusrfdns * thk0/len0, & ! scale to m/m
                                                 model%temper%dissip)
          
       else    ! first-order dissipation                                                                                                                                                               
          call glissade_interior_dissipation_first_order(ewn,  nsn,    upn,          &
                                                         ice_mask,                   &
                                                         model%stress%tau%scalar * tau0,  &  ! scale to Pa
                                                         model%stress%efvs * evs0,   &  ! scale to Pa s
                                                         model%temper%dissip)
          
       endif    ! which_ho_disp
          
    endif     ! is_restart

    if (this_rank==rtest .and. verbose_glissade) then
       i = itest
       j = jtest
       print*, 'itest, jtest =', i, j
       print*, 'k, dissip (deg/yr):'
       do k = 1, upn-1
          print*, k, model%temper%dissip(k,i,j)*scyr
       enddo
       print*, 'ubas, vbas =', model%velocity%uvel(upn,i,j), model%velocity%vvel(upn,i,j)
       print*, 'btraction =',  model%velocity%btraction(:,i,j)
       print*, 'bfricflx =', model%temper%bfricflx(i,j)
       print*, ' '
       print*, 'After glissade velocity solve (or restart): uvel, k = 1:'
       write(6,'(a8)',advance='no') '          '
       do i = itest-5, itest+5
          write(6,'(i12)',advance='no') i
       enddo
       print*, ' '
       do j = jtest+2, jtest-2, -1
          write(6,'(i8)',advance='no') j
          do i = itest-5, itest+5
             write(6,'(f12.3)',advance='no') model%velocity%uvel(1,i,j) * (vel0*scyr)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'After glissade velocity solve (or restart): vvel, k = 1:'
       write(6,'(a8)',advance='no') '          '
       do i = itest-5, itest+5
          write(6,'(i12)',advance='no') i
       enddo
       print*, ' '
       do j = jtest+2, jtest-2, -1
          write(6,'(i8)',advance='no') j
          do i = itest-5, itest+5
             write(6,'(f12.3)',advance='no') model%velocity%vvel(1,i,j) * (vel0*scyr)
          enddo
          print*, ' '
       enddo
       
    endif  ! this_rank = rtest & verbose_glissade

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 3. Third part of diagnostic solve: 
    ! Now that velocity is solved, calculate any diagnostic fields that are
    ! a function of velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! Calculate wvel, assuming grid velocity is 0.
    ! This is calculated relative to the ice sheet base, rather than a fixed reference location.
    ! Note: This current implementation for wvel only supports whichwvel=VERTINT_STANDARD
    ! Note: For glissade, wvel is diagnostic only.
    !       Pass in the basal melt rate for grounded ice only, as in Glide.
    call wvelintg(model%velocity%uvel,                        &
                  model%velocity%vvel,                        &
                  model%geomderv,                             &
                  model%numerics,                             &
                  model%velowk,                               &
                  model%geometry%thck * 0.0d0,                &  ! Just need a 2d array of all 0's for wgrd
                  model%geometry%thck,                        &
                  model%basal_melt%bmlt,                      &
                  model%velocity%wvel)
    ! Note: halos may be wrong for wvel, but since it is currently only used as an output diagnostic variable, that is OK.

    ! compute the velocity norm (for diagnostic output)
    model%velocity%velnorm(:,:,:) = sqrt(model%velocity%uvel(:,:,:)**2 + model%velocity%vvel(:,:,:)**2)

    ! compute the mean velocity

    model%velocity%uvel_mean(:,:) = 0.0d0
    model%velocity%vvel_mean(:,:) = 0.0d0

    k = 1    ! top surface velocity associated with top half of layer 1
    model%velocity%uvel_mean(:,:) = model%velocity%uvel_mean(:,:) &
                                  + model%numerics%stagsigma(k) * model%velocity%uvel(k,:,:)
    model%velocity%vvel_mean(:,:) = model%velocity%vvel_mean(:,:) &
                                  + model%numerics%stagsigma(k) * model%velocity%vvel(k,:,:)

    do k = 2, upn-1
       model%velocity%uvel_mean(:,:) = model%velocity%uvel_mean(:,:) &
                                     + (model%numerics%stagsigma(k) - model%numerics%stagsigma(k-1)) * model%velocity%uvel(k,:,:)
       model%velocity%vvel_mean(:,:) = model%velocity%vvel_mean(:,:) &
                                     + (model%numerics%stagsigma(k) - model%numerics%stagsigma(k-1)) * model%velocity%vvel(k,:,:)
    enddo

    k = upn  ! basal velocity associated with bottom half of layer (upn-1)
    model%velocity%uvel_mean(:,:) = model%velocity%uvel_mean(:,:) &
                                  + (1.0d0 - model%numerics%stagsigma(k-1)) * model%velocity%uvel(k,:,:)
    model%velocity%vvel_mean(:,:) = model%velocity%vvel_mean(:,:) &
                                  + (1.0d0 - model%numerics%stagsigma(k-1)) * model%velocity%vvel(k,:,:)

    ! Compute the vertically integrated stress tensor (Pa) and its eigenvalues.
    ! These are used for some calving schemes.
    !TODO - Put these calculations in a utility subroutine

    if ( (model%options%is_restart == RESTART_TRUE) .and. &
         (model%numerics%time == model%numerics%tstart) ) then

       ! do nothing, since the tau eigenvalues are read from the restart file

    else  ! compute the eigenvalues of the 2D horizontal stress tensor

       model%calving%tau_eigen1(:,:) = 0.0d0
       model%calving%tau_eigen2(:,:) = 0.0d0

       do j = 1, nsn
          do i = 1, ewn

             ! compute vertically averaged stress components
             tau_xx = 0.0d0
             tau_yy = 0.0d0
             tau_xy = 0.0d0

             do k = 1, upn-1
                dsigma = model%numerics%sigma(k+1) - model%numerics%sigma(k)
                tau_xx = tau_xx + tau0 * model%stress%tau%xx(k,i,j) * dsigma
                tau_yy = tau_yy + tau0 * model%stress%tau%yy(k,i,j) * dsigma
                tau_xy = tau_xy + tau0 * model%stress%tau%xy(k,i,j) * dsigma
             enddo

             ! compute the eigenvalues of the vertically integrated stress tensor
             a = 1.0d0
             b = -(tau_xx + tau_yy)
             c = tau_xx*tau_yy - tau_xy*tau_xy
             if (b*b - 4.0d0*a*c > 0.0d0) then   ! two real eigenvalues
                root = sqrt(b*b - 4.0d0*a*c)
                lambda1 = (-b + root) / (2.0d0*a)
                lambda2 = (-b - root) / (2.0d0*a)
                if (lambda1 > lambda2) then
                   model%calving%tau_eigen1(i,j) = lambda1
                   model%calving%tau_eigen2(i,j) = lambda2
                else
                   model%calving%tau_eigen1(i,j) = lambda2
                   model%calving%tau_eigen2(i,j) = lambda1
                endif
             endif  ! b^2 - 4ac > 0

          enddo   ! i
       enddo   ! j

       ! Extrapolate tau eigenvalues to inactive CF cells where the stress tensor is not computed.

       do j = 2, nsn-1
          do i = 2, ewn-1
             if (calving_front_mask(i,j) == 1 .and. &
                  model%calving%tau_eigen1(i,j) == 0.0d0 .and. model%calving%tau_eigen2(i,j) == 0.0d0) then

                sum_cell = 0
                sum1 = 0.0d0
                sum2 = 0.0d0

                ! Set the eigenvalues to the mean values in marine interior neighbors
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (marine_interior_mask(ii,jj) == 1) then
                         sum_cell = sum_cell + 1
                         sum1 = sum1 + model%calving%tau_eigen1(ii,jj)
                         sum2 = sum2 + model%calving%tau_eigen2(ii,jj)
                      endif
                   enddo
                enddo
                if (sum_cell > 0) then
                   model%calving%tau_eigen1(i,j) = sum1/sum_cell
                   model%calving%tau_eigen2(i,j) = sum2/sum_cell
                endif

             endif  ! inactive CF cell
          enddo
       enddo

    endif   ! restart

    call parallel_halo(model%calving%tau_eigen1, parallel)
    call parallel_halo(model%calving%tau_eigen2, parallel)


    if (this_rank == rtest .and. verbose_calving .and. model%options%whichcalving == EIGENCALVING) then
       print*, ' '
       print*, 'tau eigen1 (Pa), i, j, rtest =:', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(e11.3)',advance='no') model%calving%tau_eigen1(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'tau eigen2 (Pa), i, j, rtest =:', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(e11.3)',advance='no') model%calving%tau_eigen2(i,j)
          enddo
          print*, ' '
       enddo
    endif

    ! Compute the 3D strain rate tensor (s^{-1})
    ! Note: The stress tensor tau is derived by taking strain rates at quadrature points in the velocity solve.
    !       The strain rate tensor is then diagnosed from the stress tensor.
    !       These values will be incorrect on restart, since the stress tensor is not written to the restart file.

    where (model%stress%efvs > 0.0d0) 
       model%velocity%strain_rate%scalar = tau0 * model%stress%tau%scalar / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%xz = tau0 * model%stress%tau%xz / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%yz = tau0 * model%stress%tau%yz / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%xx = tau0 * model%stress%tau%xx / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%yy = tau0 * model%stress%tau%yy / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%xy = tau0 * model%stress%tau%xy / (2.d0 * evs0 * model%stress%efvs)
    elsewhere
       model%velocity%strain_rate%scalar = 0.0d0
       model%velocity%strain_rate%xz = 0.0d0
       model%velocity%strain_rate%yz = 0.0d0
       model%velocity%strain_rate%xx = 0.0d0
       model%velocity%strain_rate%yy = 0.0d0
       model%velocity%strain_rate%xy = 0.0d0
    endwhere

    ! Compute the vertically integrated strain rate tensor (s^-1) and its eigenvalues.
    ! These could be used for eigencalving, but the current eigencalving scheme uses stresses, not strain rates.

    !TODO - Check that units are 1/s
    model%calving%eps_eigen1(:,:) = 0.0d0
    model%calving%eps_eigen2(:,:) = 0.0d0

    do j = 1, nsn
       do i = 1, ewn

          ! compute vertically averaged strain rate components
          eps_xx = 0.0d0
          eps_yy = 0.0d0
          eps_xy = 0.0d0

          do k = 1, upn-1
             dsigma = model%numerics%sigma(k+1) - model%numerics%sigma(k)
             eps_xx = eps_xx + model%velocity%strain_rate%xx(k,i,j) * dsigma
             eps_yy = eps_yy + model%velocity%strain_rate%yy(k,i,j) * dsigma
             eps_xy = eps_xy + model%velocity%strain_rate%xy(k,i,j) * dsigma
          enddo

          ! compute the eigenvalues of the vertically integrated stress tensor
          a = 1.0d0
          b = -(eps_xx + eps_yy)
          c = eps_xx*eps_yy - eps_xy*eps_xy
          if (b*b - 4.0d0*a*c > 0.0d0) then   ! two real eigenvalues
             root = sqrt(b*b - 4.0d0*a*c)
             lambda1 = (-b + root) / (2.0d0*a)
             lambda2 = (-b - root) / (2.0d0*a)
             if (lambda1 > lambda2) then
                model%calving%eps_eigen1(i,j) = lambda1
                model%calving%eps_eigen2(i,j) = lambda2
             else
                model%calving%eps_eigen1(i,j) = lambda2
                model%calving%eps_eigen2(i,j) = lambda1
             endif
          endif  ! b^2 - 4ac > 0

       enddo   ! i
    enddo   ! j

    ! Extrapolate eigenvalues to inactive CF cells where the strain rate is not computed.

    do j = 2, nsn-1
       do i = 2, ewn-1
          if (calving_front_mask(i,j) == 1 .and. &
              model%calving%eps_eigen1(i,j) == 0.0d0 .and. model%calving%eps_eigen2(i,j) == 0.0d0) then

             sum_cell = 0
             sum1 = 0.0d0
             sum2 = 0.0d0

             ! Set the eigenvalues to the mean values in marine interior neighbors
             do jj = j-1, j+1
                do ii = i-1, i+1
                   if (marine_interior_mask(ii,jj) == 1) then
                      sum_cell = sum_cell + 1
                      sum1 = sum1 + model%calving%eps_eigen1(ii,jj)
                      sum2 = sum2 + model%calving%eps_eigen2(ii,jj)
                   endif
                enddo
             enddo
             if (sum_cell > 0) then
                model%calving%eps_eigen1(i,j) = sum1/sum_cell
                model%calving%eps_eigen2(i,j) = sum2/sum_cell
             endif

          endif  ! inactive CF cell
       enddo
    enddo

    call parallel_halo(model%calving%eps_eigen1, parallel)
    call parallel_halo(model%calving%eps_eigen2, parallel)

    if (this_rank == rtest .and. verbose_calving .and. model%options%whichcalving == EIGENCALVING) then
       print*, ' '
       print*, 'strain rate eigen1 (yr^-1):'
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(e11.3)',advance='no') model%calving%eps_eigen1(i,j) * scyr
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'strain rate eigen2 (yr^-1):'
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(e11.3)',advance='no') model%calving%eps_eigen2(i,j) * scyr
          enddo
          print*, ' '
       enddo
    endif  ! this_rank = rtest

    ! Compute various vertical means.
    ! TODO - Write a utility subroutine for vertical averaging

    ! Compute the vertical mean effective viscosity
    model%stress%efvs_vertavg = 0.0d0
    do j = 1, nsn
       do i = 1, ewn
          do k = 1, upn-1
             model%stress%efvs_vertavg(i,j) = model%stress%efvs_vertavg(i,j)  &
                                            + model%stress%efvs(k,i,j) * (model%numerics%sigma(k+1) - model%numerics%sigma(k))
          enddo
       enddo
    enddo

    ! Compute the vertically integrated divergence of the horizontal velocity field.
    ! Note: Units of divu and strain_rate components are s^{-1}.
    model%velocity%divu(:,:) = 0.0d0

    do j = 1, nsn
       do i = 1, ewn
          do k = 1, upn-1
             dsigma = model%numerics%sigma(k+1) - model%numerics%sigma(k)
             model%velocity%divu(i,j) = model%velocity%divu(i,j) + &
                  (model%velocity%strain_rate%xx(k,i,j) + model%velocity%strain_rate%yy(k,i,j)) * dsigma
          enddo
       enddo
    enddo

    ! magnitude of basal traction
    model%stress%btract(:,:) = sqrt(model%stress%btractx(:,:)**2 + model%stress%btracty(:,:)**2)

    ! Copy uvel and vvel to arrays uvel_extend and vvel_extend.
    ! These arrays have horizontal dimensions (nx,ny) instead of (nx-1,ny-1).
    ! They are needed for exact restart if we have nonzero velocities along the
    !  north and east edges of the global domain, as in some test problems.
    
    model%velocity%uvel_extend(:,:,:) = 0.d0
    model%velocity%vvel_extend(:,:,:) = 0.d0

    do j = 1, nsn-1
       do i = 1, ewn-1
          model%velocity%uvel_extend(:,i,j) = model%velocity%uvel(:,i,j)
          model%velocity%vvel_extend(:,i,j) = model%velocity%vvel(:,i,j)             
       enddo
    enddo

    ! Copy some additional 2D arrays to the extended grid if using the DIVA solver.

    if (model%options%which_ho_approx == HO_APPROX_DIVA) then

       model%velocity%uvel_2d_extend(:,:) = 0.d0
       model%velocity%vvel_2d_extend(:,:) = 0.d0
       do j = 1, nsn-1
          do i = 1, ewn-1
             model%velocity%uvel_2d_extend(i,j) = model%velocity%uvel_2d(i,j)
             model%velocity%vvel_2d_extend(i,j) = model%velocity%vvel_2d(i,j)             
          enddo
       enddo

       model%stress%btractx_extend(:,:) = 0.d0
       model%stress%btracty_extend(:,:) = 0.d0
       do j = 1, nsn-1
          do i = 1, ewn-1
             model%stress%btractx_extend(i,j) = model%stress%btractx(i,j)
             model%stress%btracty_extend(i,j) = model%stress%btracty(i,j)   
          enddo
       enddo

    endif

    ! If beta is not passed in from an external file, then copy beta_internal to beta.
    ! Note: beta_internal, which is weighted by the grounded ice fraction, is the actual beta field
    !        in the glissade velocity calculation.  But if users specify 'beta' instead of
    !        'beta_internal' as an output field, this copy ensures that they get the output
    !         they expect.
    !       The copy would break exact restart, however, if beta is read from an external file.
    !        In that case, users must specify 'beta_internal' to see the weighted beta field.

    if (model%options%which_ho_babc /= HO_BABC_BETA_EXTERNAL) then
       model%velocity%beta(:,:) = model%velocity%beta_internal(:,:)
    endif

    ! DIVA needs a halo update for efvs, since the latest guess (in both local and halo cells)
    ! is used to start iterating for efvs in the next time step.
    call parallel_halo(model%stress%efvs, parallel)

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 4. Fourth part of diagnostic solve: 
    ! Diagnose some quantities that are not velocity-dependent, but may be desired for output
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------

    ! basal ice temperature
    ! This is the same as temp(upn,:,:), the lowest-level of the prognostic temperature array.
    ! However, it is set to zero for ice-free columns (unlike temp(upn) = min(artm,0.0) for ice-free columns)
    ! TODO - Make btemp a prognostic array, and limit the 3D temp array to internal layer temperatures?
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > 0.0d0) then
             model%temper%btemp(i,j) = model%temper%temp(upn,i,j)
          else
             model%temper%btemp(i,j) = 0.0d0
          endif
       enddo
    enddo

    ! surface mass balance in units of mm/yr w.e.
    ! (model%climate%acab * scale_acab) has units of m/yr of ice
    ! Note: This is not necessary (and can destroy exact restart) if the SMB was already input in units of mm/yr
    if (model%options%smb_input /= SMB_INPUT_MMYR_WE) then
       model%climate%smb(:,:) = (model%climate%acab(:,:) * scale_acab) * (1000.d0 * rhoi/rhow)
    endif

    ! surface, basal and calving mass fluxes (kg/m^2/s)
    ! positive for mass gain, negative for mass loss
    model%geometry%sfc_mbal_flux(:,:) = rhoi * model%climate%acab_applied(:,:)*thk0/tim0
    model%geometry%basal_mbal_flux(:,:) = rhoi * (-model%basal_melt%bmlt_applied(:,:)) * thk0/tim0
    model%geometry%calving_flux(:,:) = rhoi * (-model%calving%calving_thck(:,:)*thk0) / (model%numerics%dt*tim0)

    ! calving rate (m/yr ice; positive for calving)
    model%calving%calving_rate(:,:) = (model%calving%calving_thck(:,:)*thk0) / (model%numerics%dt*tim0/scyr)

    ! save old masks for diagnostics
    floating_mask_old = model%geometry%floating_mask
    grounded_mask_old = model%geometry%grounded_mask

    ! set integer masks in the geometry derived type

    ! unstaggered grid
    do j = 1, nsn
       do i = 1, ewn
          if (ice_mask(i,j) == 1) then
             model%geometry%ice_mask(i,j) = 1
             if (floating_mask(i,j) == 1) then
                model%geometry%grounded_mask(i,j) = 0
                model%geometry%floating_mask(i,j) = 1
             else
                model%geometry%grounded_mask(i,j) = 1
                model%geometry%floating_mask(i,j) = 0
             endif
          else  ! ice_mask = 0
             model%geometry%ice_mask(i,j) = 0
             model%geometry%grounded_mask(i,j) = 0
             model%geometry%floating_mask(i,j) = 0
          endif
       enddo
    enddo

    ! staggered grid
    ! set ice_mask_stag = 1 at vertices with ice_mask = 1 in any neighbor cell
    do j = 1, nsn-1
       do i = 1, ewn-1
          if (ice_mask(i,j+1)==1 .or. ice_mask(i+1,j+1)==1 .or. &
              ice_mask(i,j)  ==1 .or. ice_mask(i+1,j)  ==1) then
             model%geometry%ice_mask_stag(i,j) = 1
          else
             model%geometry%ice_mask_stag(i,j) = 0
          endif
       enddo
    enddo

    !WHL - inversion debug
    ! The goal is to spin up in a way that minimizes flipping between grounded and floating.
    if (verbose_inversion .and. model%numerics%time > model%numerics%tstart .and. &
        (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or.  &
         model%options%which_ho_coulomb_c  == HO_COULOMB_C_INVERSION) ) then
       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             if (model%geometry%floating_mask(i,j) /= floating_mask_old(i,j)) then
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                if (model%geometry%floating_mask(i,j) == 1) then
                   if (grounded_mask_old(i,j) == 1) then
                      write(6,*) 'Floating_mask flip, G to F: i, j =', iglobal, jglobal
                   else
                      write(6,*) 'Floating_mask flip, O to F: i, j =', iglobal, jglobal
                   endif
                elseif (floating_mask_old(i,j) == 1) then
                   if (model%geometry%grounded_mask(i,j) == 1) then
                      write(6,*) 'Floating_mask flip, F to G: i, j =', iglobal, jglobal
                   else
                      write(6,*) 'Floating_mask flip, F to O: i, j =', iglobal, jglobal
                   endif
                endif
             endif
          enddo
       enddo
    endif

    ! Compute grounding line fluxes
    ! Note: gl_flux_east and gl_flux_north are signed fluxes computed at cell edges;
    !       gl_flux is cell-based and is found by summing magnitudes of edge fluxes.

    call glissade_grounding_line_flux(ewn,                  nsn,                 &
                                      model%numerics%dew,   model%numerics%dns,  &
                                      model%numerics%sigma,                      &
                                      model%geometry%thck,                       &
                                      model%velocity%uvel,  model%velocity%vvel, &
                                      ice_mask,             floating_mask,       &
                                      ocean_mask,                                &
                                      model%geometry%gl_flux_east,               &
                                      model%geometry%gl_flux_north,              &
                                      model%geometry%gl_flux                      )

    !------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    !
    !
    ! TODO(wjs, 2017-05-21) I don't think we should need to update lsrf and usrf
    ! here. However, glissade_velo_higher_solve and glissade_velo_sia_solve (called from
    ! glissade_velo_driver) multiply/divide topg (and other variables) by their scale
    ! factors on entry to / exit from the routine. This can lead to roundoff-level changes
    ! in topg and other variables.
    !
    ! If we don't update usrf here, then we can get roundoff-level changes in exact
    ! restart tests when running inside a climate model: In the straight-through run
    ! (without an intervening restart), the value of usrf sent to the coupler is the one
    ! set earlier in this routine, which doesn't incorporate these roundoff-level changes
    ! to topg. The restarted run, in contrast, reads the slightly-modified topg from the
    ! restart file and recomputes usrf in initialization; thus, the values of usrf that
    ! the coupler sees in the first year differ slightly from those in the
    ! straight-through run.
    !
    ! A cleaner solution could be to avoid applying these rescalings to the fundamental
    ! model variables in glissade_velo_higher_solve and glissade_velo_sia_solve - instead,
    ! introducing temporary variables in those routines to hold the scaled
    ! quantities. Then I think it would be safe to remove the following code that updates
    ! lsrf and usrf. Or, if we completely removed these scale factors from CISM, then
    ! again I think it would be safe to remove the following code.
    ! ------------------------------------------------------------------------
    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       &
                        model%climate%eus,   model%geometry%lsrf)
    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

    if (verbose_glissade .and. main_task) then
       print*, 'Done in glissade_diagnostic_variable_solve'
    endif

  end subroutine glissade_diagnostic_variable_solve

!=======================================================================

  subroutine glissade_cleanup_icefree_cells(model)

    ! Clean up prognostic variables in ice-free cells.
    ! This means seting most tracers to zero (or min(artm,0) for the case of temperature).

    use cism_parallel, only: parallel_halo

    type(glide_global_type), intent(inout) :: model   ! model instance

    integer :: nx, ny
    integer :: i, j

    type(parallel_type) :: parallel   ! info for parallel communication

    nx = model%general%ewn
    ny = model%general%nsn

    parallel = model%parallel

    ! Make sure the ice thickness is updated in halo cells
    call parallel_halo(model%geometry%thck, parallel)

    ! Set prognostic variables in ice-free columns to default values (usually zero).
    do j = 1, ny
       do i = 1, nx

          if (model%geometry%thck_old(i,j) > 0.0d0 .and. model%geometry%thck(i,j) == 0.0d0) then

             ! basal water
             model%basal_hydro%bwat(i,j) = 0.0d0

             ! thermal variables
             if (model%options%whichtemp == TEMP_INIT_ZERO) then
                model%temper%temp(:,i,j) = 0.0d0
             else
                model%temper%temp(:,i,j) = min(model%climate%artm(i,j), 0.0d0)
             endif

             if (model%options%whichtemp == TEMP_ENTHALPY) then
                model%temper%waterfrac(:,i,j) = 0.0d0
             endif

             ! other tracers
             ! Note: Tracers should be added here as they are added to the model

             if (model%options%whichcalving == CALVING_DAMAGE) then
                model%calving%damage(:,i,j) = 0.0d0
             endif

             if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
                model%geometry%ice_age(:,i,j) = 0.0d0
             endif

          endif    ! thck = 0

       enddo
    enddo

  end subroutine glissade_cleanup_icefree_cells

!=======================================================================

end module glissade

!=======================================================================
