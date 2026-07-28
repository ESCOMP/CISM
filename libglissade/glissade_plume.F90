!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_plume.F90 - part of the Community Ice Sheet Model (CISM)  
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
!
! This module contains subroutines for a sub-ice-shelf plume model.
!
! Author: William Lipscomb
!         NSF National Center for Atmospheric Research
!         Climate and Global Dynamics Laboratory
!         Boulder, CO 80303
!         USA
!         <lipscomb@ucar.edu>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
  module glissade_plume

    use glimmer_global, only: dp
    use glimmer_physcon, only: rhoi, rhow, rhoo, grav, lhci, cpw, scyr
    use glimmer_paramets, only: iulog, eps11
    use glimmer_log
    use glimmer_utils, only: point_diag
    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo,  parallel_halo, &
         parallel_global_sum, parallel_reduce_sum, parallel_is_zero, parallel_globalindex

    implicit none
    save
    private

    public :: glissade_plume_init, glissade_plume_driver, verbose_plume

    logical, parameter :: verbose_plume = .true.

    ! prescribed MISOMIP parameters from Table 4 of Asay-Davis et al. (2016)
    ! Note: cpw and lhci are defined in glimmer_physcon
    ! Kh is from LADDIE: https://github.com/erwinlambert/laddie (accessed 7/27/26).

    !TODO - Add these to a derived type or a constants module?
    real(dp), parameter :: &
         lambda1 = -0.0573d0,        & ! liquidus slope (deg/psu)
         lambda2 =  0.0832d0,        & ! liquidus intercept (deg C)
         lambda3 = -7.53d-8,         & ! liquidus pressure coefficient (deg/Pa)
                                       ! Tb = lambda1*Sb + lambda2 + lambda3*pb
         c_drag = 2.5d-3,            & ! ocean drag coefficient (unitless)
         u_tidal = 0.01d0,           & ! tidal velocity (m/s)
         Kh = 25.d0,                 & ! horizontal diffusivity of heat and salt (m^2/s)
         eos_rho_ref = 1027.51d0,    & ! reference density for linear EOS (kg/m^3)
         eos_Tref = -1.0d0,          & ! reference temperature for linear EOS (deg C)
         eos_Sref = 34.2d0,          & ! reference salinity for linear EOS (deg C)
         eos_alpha = 3.733d-5,       & ! thermal expansion coefficient for linear EOS (deg^-1)
         eos_beta = 7.843d-4,        & ! salinity contraction coefficient for linear EOS (psu^-1)
         f_coriolis = -1.405d-4        ! Coriolis parameter (s^-1) at 75 S = 2*omega*sin(75 deg) (prescribed in text)

    ! plume parameters
    !TODO - Add to the derived type?
    real(dp), parameter :: &
         D_plume0 = 10.d0,           & ! initial plume thickness (m)
         D_plume_min = 1.0d0,        & ! min plume thickness (m) where the plume exists
         D_plume_max = 200.0d0,      & ! max plume thickness (m)
         tau_relax_entrainment = 3600. ! timescale (s) for relaxing toward D_plume_min or D_plume_max
                                       ! by imposing entrainment or detrainment

    integer, parameter :: wx = 7, wy = 8   ! block size passed to point_diag

    !TODO - Make this a config option
    logical, parameter :: &
         entrainment_gaspar = .true.

!=======================================================================

  contains

!=======================================================================

  subroutine glissade_plume_init(model, ocean_data, plume)

    ! Initialize the plume properties

    use glissade_utils, only: glissade_interpolate_3d_ocean_field_to_lsrf

    ! input/ouput arguments

    type(glide_global_type), intent(inout) :: model     !> derived type holding ice-sheet info
    type(glide_ocean_data), intent(in) :: ocean_data    !> derived type holding input ocean data
    type(glide_plume), intent(inout) :: plume           !> derived type holding plume info

    ! local variables

    integer, dimension(model%general%ewn,model%general%nsn) :: mask
    real(dp), dimension(model%general%ewn,model%general%nsn) :: &
         depth                          ! depth (m) at base of plume, negative below sea level

    integer :: ewn, nsn
    real(dp) :: dew, dns
    integer :: itest, jtest, rtest      ! coordinates of diagnostic point
    type(parallel_type) :: parallel     ! info for parallel communication

    ewn = model%general%ewn
    nsn = model%general%nsn
    dew = model%numerics%dew
    dns = model%numerics%dns
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local
    parallel = model%parallel

    if (model%options%is_restart == NO_RESTART .or. model%options%is_restart == HYBRID_RESTART) then

       ! Initialize the plume

       if (verbose_plume .and. main_task) write(iulog,*) 'Initialize the plume'

       ! Compute the ocean depth at the base of the plume
       depth = model%geometry%lsrf - model%plume%D_plume

       ! Compute T_ambient and S_ambient at the base of the plume

       if (ocean_data%misomip_profile) then

          ! use the MISOMIP profiles, Eqs. 21 and 22 in Asay-Davis et al. (2016)
          plume%T_ambient = ocean_data%T0 + (ocean_data%Tbot - ocean_data%T0)*depth/ocean_data%zb_deep
          plume%S_ambient = ocean_data%S0 + (ocean_data%Sbot - ocean_data%S0)*depth/ocean_data%zb_deep

       else

          ! interpolate from 3D ocean fields; these should have been read from an input file

          if (parallel_is_zero(model%ocean_data%thetao) .or. &
              parallel_is_zero(model%ocean_data%salinity)) then
             call write_log('Need input thetao and salinity to compute T_ambient and S_ambient', GM_FATAL)
          endif

          ! Set the mask to do the computation everywhere
          !TODO - Does this create any issues where there is no plume? Should we use plume_mask = floating_mask?
          mask = 1

          call glissade_interpolate_3d_ocean_field_to_lsrf(&
               ewn,           nsn,    &
               ocean_data%nzocn,      &
               ocean_data%zocn,       &
               mask,                  &
               depth,                 &
               ocean_data%thetao,     &
               plume%T_ambient)

          call glissade_interpolate_3d_ocean_field_to_lsrf(&
               ewn,           nsn,    &
               ocean_data%nzocn,      &
               ocean_data%zocn,       &
               mask,                  &
               depth,                 &
               ocean_data%salinity,   &
               plume%S_ambient)

       endif

       ! Spin up the plume to steady state
       ! Note: Typically, the initial spin-up takes longer than the runtime update.
       !       For an ISOMIP+ experiment with a fixed ice cavity, this is all we need to do.

       call compute_plume(&
            ewn,                 nsn,                  &
            dew,                 dns,                  &
            plume%dt_plume,      plume%tplume_spinup,  &
            itest,   jtest,      rtest,                &
            parallel,                                  &
            model%geometry%thck,                       &
            model%geometry%lsrf,                       &
            model%geometry%topg,                       &
            model%climate%eus,                         &
            plume%T_ambient,     plume%S_ambient,      &
            plume%gammaT,        plume%gammaS,         &
            ocean_data%S0,                             &   ! is this needed?
            plume%D_plume,                             &
            plume%T_plume,       plume%S_plume,        &
            plume%T_basal,       plume%S_basal,        &
            plume%u_plume_east,  plume%v_plume_north,  &
            plume%u_plume,       plume%v_plume,        &
            plume%plume_speed,   plume%drho_plume,     &
            plume%entrainment,   plume%detrainment,    &
            plume%divDu_plume,                         &
            model%basal_melt%bmlt_float)

    endif  ! not a restart

  end subroutine glissade_plume_init

!****************************************************

  subroutine glissade_plume_driver(model, ocean_data, plume)

    ! Compute melt rates using a plume model, given vertical profiles of T and S in the ambient ocean
    !
    ! The benchmark application is to the MISOMIP domain, described here:
    ! See this paper for details:
    ! X. S. Asay-Davis et al. (2016), Experimental design for three interrelated
    !    marine ice sheet and ocean model intercomparison projects:
    !    MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    !    Geosci. Model Devel., 9, 2471-2497, doi: 10.5194/gmd-9-2471-2016.

    use glissade_utils, only: glissade_interpolate_3d_ocean_field_to_lsrf

    type(glide_global_type), intent(inout) :: model     !> derived type holding ice-sheet info
    type(glide_ocean_data), intent(in) :: ocean_data    !> derived type holding input ocean data
    type(glide_plume), intent(inout) :: plume           !> derived type holding plume info

    ! local variables

    integer, dimension(model%general%ewn,model%general%nsn) :: mask
    real(dp), dimension(model%general%ewn,model%general%nsn) :: &
         depth                          ! depth (m) at base of plume, negative below sea level

    integer :: ewn, nsn
    real(dp) :: dew, dns
    integer :: itest, jtest, rtest      ! coordinates of diagnostic point
    type(parallel_type) :: parallel     ! info for parallel communication

    ewn = model%general%ewn
    nsn = model%general%nsn
    dew = model%numerics%dew
    dns = model%numerics%dns
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local
    parallel = model%parallel

    if (verbose_plume .and. main_task) write(iulog,*) 'In glissade_plume_driver'

    ! Compute the ocean depth at the base of the plume
    depth = model%geometry%lsrf - model%plume%D_plume

    ! Compute T_ambient and S_ambient at the base of the plume

    if (ocean_data%misomip_profile) then

       ! use the MISOMIP profiles, Eqs. 21 and 22 in Asay-Davis et al. (2016)
       plume%T_ambient = ocean_data%T0 + (ocean_data%Tbot - ocean_data%T0)*depth/ocean_data%zb_deep
       plume%S_ambient = ocean_data%S0 + (ocean_data%Sbot - ocean_data%S0)*depth/ocean_data%zb_deep

    else

       ! interpolate from 3D ocean fields; these should have been read from an input file

       if (parallel_is_zero(model%ocean_data%thetao) .or. &
            parallel_is_zero(model%ocean_data%salinity)) then
          call write_log('Need input thetao and salinity to compute T_ambient and S_ambient', GM_FATAL)
       endif

       ! Set the mask to do the computation everywhere
       !TODO - Does this create any issues where there is no plume? Should we use plume_mask = floating_mask?
       mask = 1

       call glissade_interpolate_3d_ocean_field_to_lsrf(&
            ewn,           nsn,    &
            ocean_data%nzocn,      &
            ocean_data%zocn,       &
            mask,                  &
            depth,                 &
            ocean_data%thetao,     &
            plume%T_ambient)

       call glissade_interpolate_3d_ocean_field_to_lsrf(&
            ewn,           nsn,    &
            ocean_data%nzocn,      &
            ocean_data%zocn,       &
            mask,                  &
            depth,                 &
            ocean_data%salinity,   &
            plume%S_ambient)

    endif

    !----------------------------------------------------------------
    ! Call the plume model to compute basal melt rates for floating ice
    !----------------------------------------------------------------

       call compute_plume(&
            ewn,                 nsn,                  &
            dew,                 dns,                  &
            plume%dt_plume,      plume%tplume_runtime, &
            itest,   jtest,      rtest,                &
            parallel,                                  &
            model%geometry%thck,                       &
            model%geometry%lsrf,                       &
            model%geometry%topg,                       &
            model%climate%eus,                         &
            plume%T_ambient,     plume%S_ambient,      &
            plume%gammaT,        plume%gammaS,         &
            ocean_data%S0,                             &   ! is this needed?
            plume%D_plume,                             &
            plume%T_plume,       plume%S_plume,        &
            plume%T_basal,       plume%S_basal,        &
            plume%u_plume_east,  plume%v_plume_north,  &
            plume%u_plume,       plume%v_plume,        &
            plume%plume_speed,   plume%drho_plume,     &
            plume%entrainment,   plume%detrainment,    &
            plume%divDu_plume,                         &
            model%basal_melt%bmlt_float)

    if (verbose_plume .and. main_task) write(iulog,*) 'Updated the plume'

  end subroutine glissade_plume_driver

!****************************************************

  subroutine compute_plume(&
       nx,               ny,               &
       dx,               dy,               &
       dt_plume,         total_time,       &
       itest,  jtest,    rtest,            &
       parallel,                           &
       thck,             lsrf,             &
       topg,             eus,              &
       T_ambient,        S_ambient,        &
       gammaT,           gammaS,           &
       S0,                                 &
       D_plume,                            &
       T_plume,          S_plume,          &
       T_basal,          S_basal,          &
       u_plume_east,     v_plume_north,    &
       u_plume,          v_plume,          &
       plume_speed,      drho_plume,       &
       entrainment,      detrainment,      &
       divDu_plume,                        &
       bmlt_float)

    !----------------------------------------------------------------------------
    ! Compute the melt rate at the ice-ocean interface from a plume model of the ocean mixed layer.
    !
    ! References:
    !
    ! P.R. Holland and D.L. Feltham, 2006: The effects of rotation and ice shelf topography
    !    on frazil-laden ice shelf water plumes. J. Phys. Oceanog., 36, 2312-2327.
    ! P.R. Holland, A. Jenkins and D.M. Holland, 2008: The response of ice shelf
    !    basal melting to variations in ocean temperature. J. Climate, 21, 2558-2572.
    ! E. Lambert, A. Juling, R.S.W. van der Wal and P.R. Holland, 2023: Modelling Antarctic
    !    ice shelf basal melt patterns using the one-layer Antarctic model for dynamical downscaling
    !    of ice–ocean exchanges (LADDIE v1.0). The Cryosphere, 17, 3203-3228.
    ! E. Lambert, F. Jesse and T. Berends, 2026: The one-Layer Antarctic model for Dynamical Downscaling
    !    of Ice–ocean Exchanges (LADDIE) version 2.0. EGUsphere, 2026 (preprint).
    !
    ! The plume model is similar to LADDIE as described in the two Lambert references.
    ! Like LADDIE, it is a 2D model of the mixed layers beneath an ice shelf, which
    !  computes sub-shelf melt rates given the ambient ocean forcing.
    ! The main differences from LADDIE are:
    ! - The model runs on a regular square mesh instead of an unstructured triangular mesh.
    ! - Plume velocity components u_plume and v_plume are computed at cell edges instead of corners.
    ! - The velocity components are diagnosed from the current geometry and forcing instead
    !   of being prognosed using a momentum advection equation.
    !
    ! The model can be applied to either idealized settings (like ISOMIP+ and MISOMIP;
    ! see Asay-Davis et al. 2016) or realistic settings (like Antarctic ice-shelf cavities).
    !----------------------------------------------------------------------------


    use glissade_masks, only: glissade_get_masks
!    use glissade_grid_operators, only: glissade_centered_gradient

    ! Input/output arguments

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    real(dp), intent(in) :: &
         dt_plume,            & ! plume timestep (s) for advection
         total_time             ! how long to run the plume model (s); the goal is to reach steady state

    integer, intent(in) :: &
         itest, jtest, rtest    ! coordinates of diagnostic point

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                & ! ice thickness (m); intent(inout) to allow calving
         lsrf,                & ! ice lower surface elevation (m, negative below sea level)
         topg                   ! bedrock elevation (m, negative below sea level)

    real(dp), intent(in) ::  &
         eus                    ! eustatic sea level (m)

    real(dp), dimension(nx,ny), intent(in) ::  &
         T_ambient,           & ! ambient ocean potential temperature at depth of ice-ocean interface (deg C)
         S_ambient              ! ambient ocean salinity at depth of ice-ocean interface (psu)

    real(dp), intent(in) :: &
         gammaT,              & ! nondimensional heat transfer coefficient
         gammaS,              & ! nondimensional salt transfer coefficient
         S0                     ! sea surface salinity (psu)

    ! Note: The following are intent(inout).
    ! D_plume, T_plume and S_plume are prognosed variables that satisfy continuity equations.
    ! T_basal and S_basal are diagnosed from scratch in plume_melt rate,
    !  but the previous values are needed to compute drho_basal.
    ! All other plume variables are diagnosed here and are intent(out).
    real(dp), dimension(nx,ny), intent(inout) :: &
         D_plume,             & ! plume thickness (m)
         T_plume,             & ! plume temperature (deg C)
         S_plume,             & ! plume salinity (psu)
         T_basal,             & ! basal ice temperature (deg C)
         S_basal                ! basal ice salinity (psu)

    !TODO - Compute divDu_plume? Not currently output
    real(dp), dimension(nx,ny), intent(out) :: &
         u_plume_east,        & ! x component of plume velocity (m/s) on east edges
         v_plume_north,       & ! y component of plume velocity (m/s) on north edges
         u_plume,             & ! x component of plume velocity (m/s) averaged to cell centers
         v_plume,             & ! y component of plume velocity (m/s) averaged to cell centers
         plume_speed,         & ! plume speed averaged to cell centers (m/s);
                                ! includes a tidal component so speed >= u_tidal
         drho_plume,          & ! density difference between ambient ocean and plume (kg/m^3)
         entrainment,         & ! entrainment rate of ambient water into plume (m/s)
         detrainment,         & ! detrainment rate of plume into ambient water (m/s)
         divDu_plume,         & ! div(Du) for plume
         bmlt_float             ! melt rate at base of floating ice (m/s)

    ! Local variables

    real(dp) :: &
         time                   ! elapsed time on the way to time_total

    integer, dimension(nx,ny) :: &
         plume_mask,          & ! = 1 for cells where scalar plume variables are computed
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed;
                                ! = 0 at closed boundaries and = 2 at open boundaries
         edge_mask_north,     & ! = 1 on north edges where plume velocity is computed;
                                ! = 0 at closed boundaries and = 2 at open boundaries
         ice_mask,            & ! = 1 if ice is present (thck > 0)
         floating_mask,       & ! = 1 where ice is present and floating, else = 0
         ocean_mask,          & ! = 1 if topg is below sea level and ice is absent, else = 0
         land_mask              ! = 1 if topg is at or above sea level, else = 0

    !TODO - Remove grav_reduced and heat transfer?
    real(dp), dimension(nx,ny) :: &
         pressure,            & ! ocean pressure at base of ice (N/m^2)
         lsrf_plume,          & ! elevation of plume-ambient interface (m, negative below sea level)
         rho_plume,           & ! plume density (kg/m^3)
         rho_ambient,         & ! ambient ocean density (kg/m^3)
         rho_basal,           & ! density of water at ice base (kg/m^3)
         drho_basal,          & ! density difference between plume and ice base (kg/m^3)
         grav_reduced,        & ! reduced gravity = grav * drho_plume/rhoo (m/s^2)
         H_cavity,            & ! thickness of ocean cavity beneath the plume (m)
         heat_transfer,       & ! rate of heat transfer from plume to ice (J/m2/s)
         theta_slope,         & ! basal slope angle (rad), used for entrainment
         ustar_plume,         & ! plume friction velocity (m/s) at cell centers
         D_plume_old            ! D_plume from previous time step

    ! plume speed on cell edges
    ! Note: u plume_east and v_plume_north (the C grid velocity components) are primary
    !        and are input/output varaibles
    !       u_plume_north and v_plume_east (the D grid components) are computed as part of
    !        the velocity solution but are not used again.
    real(dp), dimension(nx,ny) ::  &
!         u_plume_east,          & ! u_plume on east edges
         v_plume_east,          & ! v_plume on east edges
         u_plume_north            ! u_plume on north edges
!         v_plume_north            ! v_plume on north edges


    real(dp), dimension(nx,ny) ::  &
         ddrho_plume_dx_east,   & ! horizontal gradient of drho_plume on east edges
         ddrho_plume_dy_east,   & !
         ddrho_plume_dx_north,  & ! horizontal gradient of drho_plume on north edges
         ddrho_plume_dy_north,  & !
         dlsrf_plume_dx_east,   & ! horizontal gradient of lsrf_plume on east edges
         dlsrf_plume_dy_east,   & !
         dlsrf_plume_dx_north,  & ! horizontal gradient of lsrf_plume on north edges
         dlsrf_plume_dy_north

    real(dp) :: &
         dlsrf_plume_dx,        & ! lsrf gradient components at cell centers
         dlsrf_plume_dy,        &
         slope                    ! magnitude of the gradient (dlsrf_dx, dlsrf_dy)

    real(dp) ::  &
!         my_max_dt,          & ! CFL-limited time step for a given cell (s)
         L2_norm,             & ! L2 norm of residual vector from continuity equation
         L2_previous            ! L2 norm from the previous convergence check

    integer :: i, j, ig, jg
!!    integer :: iter_Dplume      ! iteration counter

    ! some variables for diagnostics
    integer :: plume_count      ! no. of plume cells

    real(dp) :: &
         D_plume_mean,        & ! mean plume thickness
         T_plume_mean,        & ! volume-weighted mean plume temperature
         S_plume_mean,        & ! volume-weighted mean plume salinity
         u_plume_mean,        & ! mean speed in plume cells
         bmlt_mean,           & ! mean melt rate in plume cells
         entrainment_mean,    & ! mean entrainment in plume cells
         detrainment_mean       ! mean detrainment in plume cells

    integer, dimension(nx,ny) :: melt_mask  ! = 1 for plume cells with lsrf < -300 m
    integer :: melt_count                   ! no. of plume cells with lsrf < -300 m
    real(dp) :: melt_sum, melt_mean         ! mean melt in cells with lsrf > -300 m

    ! parameters determining convergence of iterations
    !TODO - determine L2_target
    integer, parameter :: &
         L2_target = 0.0d0           ! convergence target for dD/dt
!!         n_check_convergence = 1,    & ! interval between convergence checks for D_plume
!!         maxiter_Dplume = 999999     ! max number of iterations of outer plume-thickness loop
                                       ! terminates when plume thickness reaches virtual steady state

    if (verbose_plume .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'In glissade_compute_plume'
    endif

    ! make sure state variables are up to date in halos
    call parallel_halo(thck, parallel)
    call parallel_halo(topg, parallel)
    call parallel_halo(lsrf, parallel)
    call parallel_halo(D_plume, parallel)
    call parallel_halo(T_plume, parallel)
    call parallel_halo(S_plume, parallel)
    call parallel_halo(T_ambient, parallel)
    call parallel_halo(S_ambient, parallel)

    !----------------------------------------------------------------
    ! compute some masks
    !----------------------------------------------------------------

    call glissade_get_masks(&
         nx,                  ny,           &
         parallel,                          &
         thck,                topg,         &
         eus,                 0.0d0,        &  ! thklim = 0
         ice_mask,                          &
         floating_mask = floating_mask,     &
         land_mask = land_mask,             &
         ocean_mask = ocean_mask)

!!    call parallel_halo(floating_mask, parallel)

    if (verbose_plume) then
       if (this_rank == rtest) write(iulog,*) 'Input ice geometry:'
       call point_diag(thck, 'thck (m)', itest, jtest, rtest, wx, wy)
       call point_diag(lsrf, 'lsrf (m)', itest, jtest, rtest, wx, wy)
       call point_diag(topg, 'topg (m)', itest, jtest, rtest, wx, wy)
       call point_diag(lsrf - topg, 'lsrf - topg (m)', itest, jtest, rtest, wx, wy)
       call point_diag(floating_mask, 'floating_mask', itest, jtest, rtest, wx, wy)
    endif

    ! Compute a mask that identifies where the plume is located
    !TODO - Refine this mask? Or cite Lambert et al. 2026 as justification?

    plume_mask = floating_mask

    !WHL - commented out for now
    ! Restrict the plume to end a few km from the calving front.
    ! Note: This may be unnecessary in CISM MISOMIP runs if the ice thickness field
    !       is smooth near the calving front.  However, the prescribed ISOMIP+ field
    !       has some strange thickness undulations near the calving front near the top
    !       and bottom domain boundaries.  Since the assumptions of the plume model
    !       may not hold near the calving front (because of lateral mixing), it may
    !       be physically justifiable anyway to cut off the model short of the front.
    !       For now I'm using a prescribed calving limit, but a thickness criterion
    !       might work too.
!    do j = 1, ny
!       do i = 1, nx
!          if  (x1(i) > plume_xmax) then
!             plume_mask(i,j) = 0
!          endif
!       enddo
!    enddo

!    call parallel_halo(plume_mask, parallel)

    ! Compute the density of the ambient ocean
    rho_ambient = eos_rho_ref * (1.d0 - eos_alpha * (T_ambient - eos_Tref)  &
                                      + eos_beta  * (S_ambient - eos_Sref) )

    ! Compute the pressure at the lower ice surface.
    pressure = -rhoo*grav*lsrf

    ! Compute the cavity thickness
    H_cavity = max(lsrf - topg, 0.0d0)

    !----------------------------------------------------------------------------
    ! Initialize D_plume, T_plume, S_plume, T_basal and S_basal as needed
    ! On the first call, the input values are zero and these fields must be initialized everywhere.
    ! On subsequent calls, these fields are initialized only if the input values are zero.
    !
    ! Note: Setting S_plume = S0 means that drho_plume = rho_ambient - rho_plume will decrease in the upslope direction,
    !        giving an upslope PGF.
    !       Setting both T_plume and S_plume to ambient values would give zero PGF, velocities, and drho_plume.
    !       The entrainment can be infinite when drho_plume = 0.
    !
    ! Note: T_basal and S_basal are diagnosed in plume_melt_rate without regard to their initial values.
    !       However, T_basal and S_basal from the previous step may be needed to compute drho_basal
    !        for entrainment.
    !
    !----------------------------------------------------------------------------

    ! loop over all cells (fields on the rhs are up to date in halos)
    do j = 1, ny
       do i = 1, nx
          if (plume_mask(i,j) == 1) then
             ! set plume to ambient temperature but with low salinity to create a positive drho_plume
             if (D_plume(i,j) == 0.0d0) D_plume(i,j) = min(D_plume0, H_cavity(i,j))
             if (T_plume(i,j) == 0.0d0) T_plume(i,j) = T_ambient(i,j)
             if (S_plume(i,j) == 0.0d0) S_plume(i,j) = S0
             ! set ice base to freezing temperature with zero salinity to create a positive drho_basal
             if (S_basal(i,j) == 0.0d0) S_basal(i,j) = 0.0d0
             if (T_basal(i,j) == 0.0d0) T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
          else   ! plume_mask = 0
             ! OK to zero out all of these?
             T_plume(i,j) = 0.0d0
             S_plume(i,j) = 0.0d0
             D_plume(i,j) = 0.0d0
             T_basal(i,j) = 0.0d0
             S_basal(i,j) = 0.0d0
          endif
       enddo
    enddo

    ! Mask out the plume in halo cells that lie outside the global domain.
    ! Also, identify global boundary cells for later use.
    ! Note: Ideally, we could zero out plume variables in the halo call by using an appropriate BC.
    ! TODO: Handle plume_mask_cell with no-penetration BCs?

    !WHL - commented out for now
!    global_bndy_west(:,:) = 0.0d0
!    global_bndy_east(:,:) = 0.0d0
!    global_bndy_south(:,:) = 0.0d0
!    global_bndy_north(:,:) = 0.0d0

!    do j = 1, ny
!       do i = 1, nx

!          call parallel_globalindex(i, j, iglobal, jglobal)

!          if (iglobal < 1 .or. iglobal > global_ewn .or. &
!              jglobal < 1 .or. jglobal > global_nsn) then
!             plume_mask_cell(i,j) = 0
!          endif

!          if (iglobal == 1) global_bndy_west(i,j) = 1
!          if (iglobal == global_ewn) global_bndy_east(i,j) = 1
!          if (jglobal == 1) global_bndy_south(i,j) = 1
!          if (jglobal == global_nsn) global_bndy_north(i,j) = 1

!       enddo
!    enddo

    !----------------------------------------------------------------------------
    ! Compute a mask to identify cell edges where plume velocities are computed.
    ! If both adjacent cells have plume_mask_cell = 1, then edge_mask = 1.
    ! At closed boundaries (one adjecent cell is grounded), edge_mask = 0.
    ! At open boundaries (one adjacent cell is open ocean), edge_mask = 2.
    !TODO - Free slip for flow parallel to edges?
    !----------------------------------------------------------------------------

    edge_mask_east = 0
    edge_mask_north = 0

    ! loop over all edges of locally owned cells
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          if (plume_mask(i,j) == 1) then
             if (plume_mask(i+1,j) == 1) then
                edge_mask_east(i,j) = 1
             elseif (lsrf(i+1,j) == topg(i+1,j)) then  ! closed boundary
                edge_mask_east(i,j) = 0
             elseif (lsrf(i+1,j) > topg(i+1,j)) then   ! open boundary
                edge_mask_east(i,j) = 2
             endif
             if (plume_mask(i,j+1) == 1) then
                edge_mask_north(i,j) = 1
             elseif (lsrf(i,j+1) == topg(i,j+1)) then  ! closed boundary
                edge_mask_north(i,j) = 0
             elseif (lsrf(i,j+1) > topg(i,j+1)) then   ! open boundary
                edge_mask_north(i,j) = 2
             endif
          endif
       enddo
    enddo

    call parallel_halo(edge_mask_east, parallel)
    call parallel_halo(edge_mask_north, parallel)

    ! Mask out edge_mask_east and edge_mask_north at edges along or outside the global domain.
    ! Note: The west and east borders have iglobal indices 0 and global_ewn, respectively.
    !       The south and north borders have jglobal indices 0 and global_nsn, respectively.
    ! TODO: Handle edge masks with no-penetration BCs?

    do j = 1, ny
       do i = 1, nx
          call parallel_globalindex(i, j, ig, jg, parallel)

          if (ig <= 0 .or. ig >= parallel%global_ewn .or. &  ! along or beyond EW boundary
              jg <= 0 .or. jg >  parallel%global_nsn) then   ! beyond NS boundary
             edge_mask_east(i,j) = 0
          endif

          if (jg <= 0 .or. jg >= parallel%global_nsn .or. &  ! along or beyond NS boundary
              ig <= 0 .or. ig >  parallel%global_ewn) then   ! beyond EW boundary
             edge_mask_north(i,j) = 0
          endif

       enddo
    enddo

    ! Compute masks for inhibiting flow toward walls of grounded ice

    !WHL - commented out for now
    ! Initialize masks to 1.0 (implying no reduction of the velocity component)
!    edge_mask_east_reduce_v(:,:) = 1.0d0
!    edge_mask_north_reduce_u(:,:) = 1.0d0

    ! Reset the masks to 0.0 or 0.5 adjacent to grounded ice
!    do j = 1, ny
!       do i = 1, nx

          ! identify east edges with a wall of grounded ice to the north or south
!          if (edge_mask_east(i,j) == 1) then
!             if ( (H_cavity(i,j+1) == 0.0d0 .and. H_cavity(i+1,j+1) == 0.0d0) .or.  &
!                  (H_cavity(i,j-1) == 0.0d0 .and. H_cavity(i+1,j-1) == 0.0d0) ) then
!                ! full wall; zero out the v component
!                edge_mask_east_reduce_v(i,j) = 0.0d0
!             elseif (H_cavity(i,j+1) == 0.0d0 .or. H_cavity(i+1,j+1) == 0.0d0 .or.  &
!                     H_cavity(i,j-1) == 0.0d0 .or. H_cavity(i+1,j-1) == 0.0d0) then
!                ! half wall; reduce the v component
!                edge_mask_east_reduce_v(i,j) = 0.5d0
!             endif
!          endif

!          ! identify north edges with a wall of grounded ice to the east or west
!          if (edge_mask_north(i,j) == 1) then
!             if ( (H_cavity(i-1,j+1) == 0.0d0 .and. H_cavity(i+1,j+1) == 0.0d0) .or.  &
!                  (H_cavity(i-1,j)   == 0.0d0 .and. H_cavity(i+1,j)   == 0.0d0) ) then
!                ! full wall; zero out the v component
!                edge_mask_north_reduce_u(i,j) = 0.0d0
!             elseif (H_cavity(i-1,j+1) == 0.0d0 .or. H_cavity(i+1,j+1) == 0.0d0 .or.  &
!                     H_cavity(i-1,j)   == 0.0d0 .or. H_cavity(i+1,j)   == 0.0d0) then
!                ! half wall; reduce the v component
!                edge_mask_north_reduce_u(i,j) = 0.5d0
!             endif
!          endif
!       enddo
!    enddo

    !TODO - Check whether these comments are still accurate
    ! Compute masks for edges with nonzero fluxes.
    ! These masks includes all edges with edge_mask_east/north = 1 (where velocity is computed).
    ! In addition, these masks include edges that have a plume on one side of the edge and
    !  open water on the other.
    ! These edges have nonzero velocity extrapolated from neighboring edges, and thus are included
    !  in computations of the divergence.

    !WHL - commented out for now
!    divu_mask_east(:,:) = edge_mask_east(:,:)
!    divu_mask_north(:,:) = edge_mask_north(:,:)

    ! east edges
!    do j = 1, ny
!       do i = 1, nx-1
!          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 0 .and. global_bndy_east(i,j) == 0) then
!             if (lsrf(i+1,j) == 0.0d0 .or. floating_mask(i+1,j) == 1) then
!                ! water in cell (i+1,j); get plume velocity from edge (i-1,j)
!                divu_mask_east(i,j) = 1
!             endif
!          elseif (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i+1,j) == 1 .and. global_bndy_west(i+1,j) == 0) then
!             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
!                ! water in cell (i,j); get plume velocity from edge (i+1,j)
!                divu_mask_east(i,j) = 1
!             endif
!          endif
!       enddo
!    enddo

    ! north edges
!    do j = 1, ny
!       do i = 1, nx
!          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 0 .and. global_bndy_north(i,j) == 0) then
!             if (lsrf(i,j+1) == 0.0d0 .or. floating_mask(i,j+1) == 1) then
!                ! water in cell (i,j+1); get plume velocity from edge (i,j-1)
!                divu_mask_north(i,j) = 1
!             endif
!          elseif (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i,j+1) == 1 .and. global_bndy_south(i,j+1) == 0) then
!             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
!                ! water in cell (i,j); get plume velocity from edge (i,j+1)
!                divu_mask_north(i,j) = 1
!             endif
!          endif
!       enddo   ! i
!    enddo   ! j

    !TODO - Should this be the lower plume surface?
    !       Modify to omit global_bndy arrays

    ! Compute the horizontal gradient of the lower ice surface.
    ! This is used to compute the pressure gradient force at velocity points, and for entrainment.
    ! Note: There are a couple of different ways to compute the PGF.
    !       (1) Jenkins et al. (1991) and HJH (2008) use grad(lsrf)
    !       (2) Holland & Feltham (2006) use grad(lsrf_plume) along with a density gradient.
    !       Method (1) is simpler and has the advantage that grad(lsrf) does not vary during plume evolution,
    !        making the PGF more stable (though possibly not as accurate).
    ! Note: The first 'lsrf' is a required argument for the subroutine.
    !       The second 'lsrf' happens to be the field whose gradient we're computing.

!    call compute_edge_gradients(&
!         nx,              ny,          &
!         dx,              dy,          &
!!         global_bndy_east,             &
!!         global_bndy_west,             &
!!         global_bndy_north,            &
!!         global_bndy_south,            &
!         plume_mask_cell,              &
!         floating_mask,                &
!         lsrf,                         &
!         lsrf,                         &
!         dlsrf_dx_east,      dlsrf_dy_east,  &
!         dlsrf_dx_north,     dlsrf_dy_north)

    !----------------------------------------------------------------------------
    ! Initialize some fields related to plume dynamics and melting.
    !----------------------------------------------------------------------------

    u_plume = 0.0d0
    v_plume = 0.0d0
    u_plume_east = 0.0d0
    v_plume_east = 0.0d0
    u_plume_north = 0.0d0
    v_plume_north = 0.0d0
    plume_speed = 0.0d0
    ustar_plume = 0.0d0
    drho_basal = 0.0d0
    drho_plume = 0.0d0
    grav_reduced = 0.0d0  !TODO - Omit?
    entrainment = 0.0d0
    detrainment = 0.0d0
    bmlt_float = 0.0d0

    divDu_plume = 0.0d0
    D_plume_old = D_plume

    if (verbose_plume) then
       if (this_rank == rtest) write(iulog,*) 'Initial plume-related fields:'
       call point_diag(plume_mask, 'plume_mask', itest, jtest, rtest, wx, wy)
       call point_diag(edge_mask_east,  'edge_mask_east', itest, jtest, rtest, wx, wy)
       call point_diag(edge_mask_north, 'edge_mask_north', itest, jtest, rtest, wx, wy)
       call point_diag(H_cavity, 'H_cavity', itest, jtest, rtest, wx, wy)
       call point_diag(T_ambient, 'T_ambient (deg C)', itest, jtest, rtest, wx, wy)
       call point_diag(S_ambient, 'S_ambient (psu)', itest, jtest, rtest, wx, wy)
       call point_diag(rho_ambient, 'rho_ambient (kg/m3)', itest, jtest, rtest, wx, wy)
       call point_diag(D_plume, 'D_plume (m)', itest, jtest, rtest, wx, wy)
       call point_diag(T_plume, 'T_plume (deg C)', itest, jtest, rtest, wx, wy)
       call point_diag(S_plume, 'S_plume (psu)', itest, jtest, rtest, wx, wy)
    endif

    time = 0.0d0

    do while(time < total_time)

       !----------------------------------------------------------------------------
       ! Iterate the plume with timestep dt_plume until we reach total_time,
       !  which ideally is long enough for the plume to reach steady state.
       ! Each timestep consists of (1) a velocity solve, (2) entrainment, detrainment,
       !  and melt rate computations, and (3) solutions of the transport equations.
       !----------------------------------------------------------------------------

       ! advance the time (units of s)
       time = time + dt_plume

       if (verbose_plume .and. this_rank == rtest) then
          write(iulog,*)
          write(iulog,*) 'Iterate plume, time (s) =', time
       endif

       !TODO - not currently computing this
       ! initialize the L2 norm to an arbitrary big number
       L2_previous = huge(0.0d0)

       ! Compute the plume density, given the current values of T_plume and S_plume.
       ! Then find the density difference between the ambient ocean and the plume.
       ! Compute the reduced gravity as a function of the density difference.

       rho_plume = eos_rho_ref * (1.d0 - eos_alpha * (T_plume - eos_Tref)  &
                                       + eos_beta  * (S_plume - eos_Sref) )

       where (plume_mask == 1)
          drho_plume = rho_ambient - rho_plume
          grav_reduced = (grav/rhoo) * drho_plume
       endwhere

       ! Compute the density at the ice base, given the current values of T_basal and S_basal.
       ! Then find the density difference between the plume and the ice base.
       ! Note: This calculation is the reason T_basal and S_basal are written to the restart file.

       rho_basal = eos_rho_ref * (1.d0 - eos_alpha * (T_basal - eos_Tref)  &
                                       + eos_beta  * (S_basal - eos_Sref) )

       where (plume_mask == 1)
          drho_basal = rho_plume - rho_basal
       endwhere

       ! Compute the elevation of the lower plume boundary

       lsrf_plume = lsrf - D_plume

       !----------------------------------------------------------------------------
       ! Compute horizontal gradients of lsrf_plume and drho_plume at each edge.
       ! Note: The subroutine includes halo updates.
       !TODO - Currently, only computes where edge_mask = 1, not edge_mask = 2
       !TODO - Put these in the velocity subroutine?
       !----------------------------------------------------------------------------

       call compute_edge_gradients(&
            nx,                   ny,                    &
            dx,                   dy,                    &
            parallel,                                    &
            edge_mask_east,       edge_mask_north,       &
            lsrf_plume,                                  &
            dlsrf_plume_dx_east,  dlsrf_plume_dy_east,   &
            dlsrf_plume_dx_north, dlsrf_plume_dy_north)

       call compute_edge_gradients(&
            nx,                   ny,                    &
            dx,                   dy,                    &
            parallel,                                    &
            edge_mask_east,       edge_mask_north,       &
            drho_plume,                                  &
            ddrho_plume_dx_east,  ddrho_plume_dy_east,   &
            ddrho_plume_dx_north, ddrho_plume_dy_north)

       !----------------------------------------------------------------------------
       ! Compute u_plume and v_plume at each edge
       ! Note: u_plume_east and v_plume_north are perpendicular to edges,
       !        whereas v_plume_north and u_plume_east are parallel to edges.
       !       Computing both u and v at each edge leads to a more graceful treatment
       !        of the Coriolis terms than computing the perpendicular components alone.
       !----------------------------------------------------------------------------

       if (verbose_plume) then
          call point_diag(rho_plume, 'rho_plume (kg/m3)', itest, jtest, rtest, wx, wy)
          call point_diag(drho_plume, 'drho_plume (kg/m3)', itest, jtest, rtest, wx, wy)
          call point_diag(rho_basal, 'rho_basal (kg/m3)', itest, jtest, rtest, wx, wy)
          call point_diag(drho_basal, 'drho_basal (kg/m3)', itest, jtest, rtest, wx, wy)
!          call point_diag(grav_reduced, 'grav_reduced (m/s2)', itest, jtest, rtest, wx, wy)
          call point_diag(dlsrf_plume_dx_east, 'dlsrf_dx_east', itest, jtest, rtest, wx, wy, '(f10.5)')
          call point_diag(dlsrf_plume_dy_north, 'dlsrf_dy_north', itest, jtest, rtest, wx, wy, '(f10.5)')
          if (this_rank == rtest) write(iulog,*) 'Compute plume velocity'
       endif

       !TODO - Pass in lsrf, drho_plume?
       call compute_plume_velocity(&
            nx,           ny,      &
            dx,           dy,      &
            itest, jtest, rtest,   &
            parallel,              &
            plume_mask,            &
            edge_mask_east,        &
            edge_mask_north,       &
            D_plume,               &
            H_cavity,              &
            drho_plume,            &
            ddrho_plume_dx_east,   &
            ddrho_plume_dy_east,   &
            ddrho_plume_dx_north,  &
            ddrho_plume_dy_north,  &
            dlsrf_plume_dx_east,   &
            dlsrf_plume_dy_east,   &
            dlsrf_plume_dx_north,  &
            dlsrf_plume_dy_north,  &
            u_plume_east,          &
            v_plume_east,          &
            u_plume_north,         &
            v_plume_north)

       !TODO - Not needed, because the update is done in the velocity solver
       ! halo updates for the velocity components used below
       ! (u_plume_north and v_plume_east are not used below)
       call parallel_halo(u_plume_east, parallel)
       call parallel_halo(v_plume_north, parallel)

       !--------------------------------------------------------------------
       ! Compute the plume speed and friction velocity at cell centers
       !--------------------------------------------------------------------

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (plume_mask(i,j) == 1) then
                u_plume(i,j) = (u_plume_east(i-1,j) + u_plume_east(i,j)) / 2.0d0
                v_plume(i,j) = (v_plume_north(i,j-1) + v_plume_north(i,j)) / 2.0d0
             endif
          enddo
       enddo

       call parallel_halo(u_plume, parallel)
       call parallel_halo(v_plume, parallel)

       do j = 1, ny
          do i = 1, nx
             plume_speed(i,j) = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)
             ustar_plume(i,j) = sqrt(c_drag) * plume_speed(i,j)
          enddo
       enddo

       !--------------------------------------------------------------------
       ! Compute the entrainment rate
       ! TODO - Make this a plume config option: 0 and 1
       ! Note: All relevant quantities for entrainment are up to date in halos.
       !--------------------------------------------------------------------

       if (entrainment_gaspar) then

          !-----------------------------------------------------------------
          ! Following Gaspar (1988), Gladish et al.(2012) and Lambert et al. (2023):
          ! Compute entrainment by relating TKE sources (friction velocity)
          !  to TKE sinks (entrainment and melt).
          ! Note: Can remove T_basal and S_basal from the restart file if
          !       not needed to compute rho_basalfor entrainment
          !-----------------------------------------------------------------

          call plume_entrainment_gaspar(&
               nx,           ny,      &
               dx,           dy,      &
               itest, jtest, rtest,   &
               parallel,              &
               plume_mask,            &
               ustar_plume,           &
               bmlt_float,            &
               drho_plume,            &
               drho_basal,            &
               H_cavity,              &
               D_plume,               &
               dt_plume,              &
               entrainment,           &
               detrainment)

       else

          !-----------------------------------------------------------------
          ! Following Bo Pederson (1980) and Jenkins (1991):
          ! Entrainment is a function of the plume speed and the basal slope angle.
          !-----------------------------------------------------------------

          ! Compute the slope angle at cell centers

          theta_slope = 0.0d0

          do j = nhalo+1, ny-nhalo
             do i = nhalo+1, nx-nhalo
                if (plume_mask(i,j) == 1) then
                   dlsrf_plume_dx = (dlsrf_plume_dx_east(i-1,j) + dlsrf_plume_dx_east(i,j)) / 2.d0
                   dlsrf_plume_dy = (dlsrf_plume_dy_north(i,j-1) + dlsrf_plume_dy_north(i,j)) / 2.d0
                   slope = sqrt(dlsrf_plume_dx**2 + dlsrf_plume_dy**2)
                   theta_slope(i,j) = atan(slope)
                endif
             enddo
          enddo

          call parallel_halo(theta_slope, parallel)

          call plume_entrainment(&
               nx,         ny,      &
               dx,         dy,      &
               itest, jtest, rtest, &
               parallel,            &
               plume_mask,          &
               theta_slope,         &
               plume_speed,         &
               H_cavity,            &
               D_plume,             &
               dt_plume,            &
               entrainment,         &
               detrainment)

       endif   ! entrainment option

       if (verbose_plume) then
          call point_diag(theta_slope, 'theta_slope (rad)', itest, jtest, rtest, wx, wy)
          call point_diag(ustar_plume, 'ustar_plume (m/s)', itest, jtest, rtest, wx, wy)
          call point_diag(plume_speed, 'plume_speed (m/s)', itest, jtest, rtest, wx, wy)
          call point_diag(entrainment*scyr, 'entrainment (m/yr)', itest, jtest, rtest, wx, wy)
          call point_diag(detrainment*scyr, 'detrainment (m/yr)', itest, jtest, rtest, wx, wy)
          if (this_rank == rtest) write(iulog,*) 'Compute melt rate: gammaT, gammaS =', gammaT, gammaS
       endif

       !--------------------------------------------------------------------
       ! Compute the basal melt rate, temperature and salinity at the plume-ice interface,
       ! Note: All relevant quantities for the melt rate are up to date in halos.
       !--------------------------------------------------------------------

       call plume_melt_rate(&
            nx,         ny,      &
            itest, jtest, rtest, &
            parallel,            &
            plume_mask,          &
            gammaT,              &
            gammaS,              &
            pressure,            &
            ustar_plume,         &
            D_plume,             &
            T_plume,             &
            S_plume,             &
            T_basal,             &
            S_basal,             &
            bmlt_float)

       ! Compute the rate of heat transfer (J/m^2/s = W/m2) from the plume to the ice base.
       ! This is equal to the melt rate (m/s) times the latent heat (J/m3) of the ice.
       !TODO - Pass bmlt_float to the transport scheme instead; save a separate array
       where (plume_mask == 1)
          heat_transfer = rhoi*lhci*bmlt_float
       elsewhere
          heat_transfer = 0.0d0
       endwhere

       if (verbose_plume) then
          if (this_rank == rtest) write(iulog,*) 'After melt calculation:'
          call point_diag(T_basal, 'T_basal (deg C)', itest, jtest, rtest, wx, wy)
          call point_diag(S_basal, 'S_basal (psu)', itest, jtest, rtest, wx, wy)
          call point_diag(bmlt_float*scyr, 'bmlt_float (m/yr)', itest, jtest, rtest, wx, wy)
          call point_diag(heat_transfer, 'heat transfer (W/m2)', itest, jtest, rtest, wx, wy)
       endif

       !TODO - Not sure if the following is needed.
       ! Determine the time step based on a CFL condition.
       ! Should be stable with a CFL number up to 1.0, but limit to 0.5 to be on the safe side.
       !WHL - Is this necessary to do for each iteration?

          !TODO - If reducing dt_plume, then it shouldn't be a parameter above
!          dt_plume = dt_plume_max
!          imax = 1
!          jmax = 1
!          do j = nhalo+1, ny-nhalo
!             do i = nhalo+1, nx-nhalo
!                if (plume_mask(i,j) == 1) then
!                   my_max_dt = 0.5d0*dx / max( abs(u_plume_east(i,j)),  abs(u_plume_east(i-1,j)), &
!                                               abs(v_plume_north(i,j)), abs(v_plume_north(i,j-1)) )
!                   if (my_max_dt < dt_plume) then
!                      dt_plume = my_max_dt
!                      imax = i
!                      jmax = j
!                   endif
!                endif
!             enddo
!          enddo

!          if (verbose_plume .and. main_task .and. dt_plume < dt_plume_max) then
!             print*, 'Limited dt_plume =', dt_plume
!          endif

       !--------------------------------------------------------------------
       ! Solve transport equations for D_plume, T_plume and S_plume,
       !  given u_plume, v_plume, entrainment, detrainment and bmlt_float.
       ! Note: Entrained water has ambient properties (T_ambient, S_ambient).
       !       Meltwater has basal properties (T basal, S basal).
       !--------------------------------------------------------------------

       call plume_transport(&
            nx,           ny,     &
            dx,           dy,     &
            itest, jtest, rtest,  &
            parallel,             &
            dt_plume,             &
            plume_mask,           &
            edge_mask_east,       &
            edge_mask_north,      &
            u_plume_east,         &
            v_plume_north,        &
            entrainment,          &
            detrainment,          &
            bmlt_float,           &
            heat_transfer,        &
            T_ambient,            &
            S_ambient,            &
            T_basal,              &
            S_basal,              &
            D_plume,              &
            T_plume,              &
            S_plume)

       ! halo updates
       call parallel_halo(D_plume, parallel)
       call parallel_halo(T_plume, parallel)
       call parallel_halo(S_plume, parallel)

       if (verbose_plume) then
          call point_diag((D_plume - D_plume_old)/(dt_plume/scyr), 'dD/dt (m/yr)', itest, jtest, rtest, wx, wy, '(f10.0)')
       endif

       !TODO - Modify this, since iter_Dplume isn't being updated
!       if (iter_Dplume >=2 .and. mod(iter_Dplume, n_check_convergence) == 0) then  ! check for convergence

          !TODO - Compute L2_norm based on dD/dt?

          ! Check for convergence: dD/dt is small everywhere
!          if (L2_norm < L2_target) then
!             if (verbose_plume .and. main_task) then
!                write(iulog,*) 'Continuity converged, time, iter, L2_norm =', time, iter_Dplume, L2_norm
!             endif
!          elseif (L2_norm < L2_previous) then ! iteration is converging; keep going
!             if (verbose_plume .and. main_task) then
!                write(iulog,*) 'Continuty not yet converged, time, iter, L2_norm =', time, iter_Dplume, L2_norm
!             endif
!          elseif (L2_norm >= L2_previous) then ! iteration is not converging
!             if (verbose_plume .and. main_task) then
!                write(iulog,*) 'Continuty not converging, time, iter, L2_norm =', time, iter_Dplume, L2_norm
!             endif
!          endif
!       endif   ! mod(iter_Dplume, n_check_convergence) = 0

       ! save variables from this iteration
       D_plume_old = D_plume
       L2_previous = L2_norm


       if (verbose_plume) then

          ! Compute mean melting at depths below 300 m.
          ! Asay-Davis et al. (2016) suggest tuning the mean melt to 30 m
          where (plume_mask == 1 .and. lsrf < -300.d0)
             melt_mask = 1
          elsewhere
             melt_mask = 0
          endwhere
          melt_count = parallel_global_sum(melt_mask, parallel)
          melt_sum = parallel_global_sum(bmlt_float, parallel, melt_mask)
          melt_mean = melt_sum/melt_count

          ! more global diagnostics
          plume_count = parallel_global_sum(plume_mask, parallel)
          D_plume_mean = parallel_global_sum(D_plume, parallel) / plume_count
          T_plume_mean = parallel_global_sum(D_plume*T_plume, parallel) / (D_plume_mean*plume_count)
          S_plume_mean = parallel_global_sum(D_plume*S_plume, parallel) / (D_plume_mean*plume_count)
          u_plume_mean = parallel_global_sum(plume_speed, parallel) / plume_count
          entrainment_mean = parallel_global_sum(entrainment, parallel) / plume_count
          detrainment_mean = parallel_global_sum(detrainment, parallel) / plume_count
          bmlt_mean = parallel_global_sum(bmlt_float, parallel) / plume_count
          if (main_task) then
             write(iulog,*) ' '
             write(iulog,*) 'Global plume diagnostics, time =', time
             write(iulog,*) 'no. of plume cells =', plume_count
             write(iulog,*) 'mean D (m) =', D_plume_mean
             write(iulog,*) 'mean T (deg C) =', T_plume_mean
             write(iulog,*) 'mean S (psu) =', S_plume_mean
             write(iulog,*) 'mean u (m/s) =', u_plume_mean
             write(iulog,*) 'mean Ent (m/yr) =', entrainment_mean*scyr
             write(iulog,*) 'mean Det (m/yr) =', detrainment_mean*scyr
             write(iulog,*) 'mean m (m/yr) =', bmlt_mean*scyr
             write(iulog,*) 'mean m below 300 m =', melt_mean*scyr
          endif
       endif

    enddo   ! time < total_time

    if (verbose_plume .and. main_task) then
       write(iulog,*) 'Plume calculation done'
    endif

  end subroutine compute_plume

!****************************************************

  subroutine compute_edge_gradients(&
       nx,               ny,              &
       dx,               dy,              &
       parallel,                          &
       edge_mask_east,   edge_mask_north, &
       field,                             &
       df_dx_east,       df_dy_east,      &
       df_dx_north,      df_dy_north)

    ! Compute horizontal gradients at the east and north edges of each cell

    ! input/output arguments

    integer, intent(in) ::  &
         nx,     ny                  ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy                  ! grid cell size (m)

    type(parallel_type), intent(in) :: &
         parallel                    ! info for parallel communication

    integer, dimension(nx,ny), intent(in) :: &
         edge_mask_east,           & ! = 1 for east edges where gradients are computed
         edge_mask_north             ! = 1 for north edges where gradients are computed

    real(dp), dimension(nx,ny), intent(in) :: &
         field                       ! input scalar field

    real(dp), dimension(nx,ny), intent(out) :: &
         df_dx_east,  df_dy_east,  & ! gradients on east edges
         df_dx_north, df_dy_north    ! gradients on north edges

    ! local variables

    integer :: i, j

    ! initialize
    df_dx_east = 0.0d0
    df_dy_east = 0.0d0
    df_dx_north = 0.0d0
    df_dy_north = 0.0d0

    ! Compute x gradients on east edges and y gradients on north edges

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (edge_mask_east(i,j) == 1) then
             df_dx_east(i,j)  = (field(i+1,j) - field(i,j)) / dx
          endif
          if (edge_mask_north(i,j) == 1) then
             df_dy_north(i,j) = (field(i,j+1) - field(i,j)) / dy
          endif
       enddo
    enddo

    call parallel_halo(df_dx_east, parallel)
    call parallel_halo(df_dy_north, parallel)

    ! Interpolate to get y gradients on east edges and x gradients on north edges

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (edge_mask_east(i,j) == 1) then
             df_dy_east(i,j)  = 0.25d0 * (df_dy_north(i,j)   + df_dy_north(i+1,j) &
                                        + df_dy_north(i,j-1) + df_dy_north(i+1,j-1))
          endif
          if (edge_mask_north(i,j) == 1) then
             df_dx_north(i,j) = 0.25d0 * (df_dx_east(i-1,j+1) + df_dx_east(i,j+1)  &
                                        + df_dx_east(i-1,j)   + df_dx_east(i,j))
          endif
       enddo
    enddo

    call parallel_halo(df_dy_east, parallel)
    call parallel_halo(df_dx_north, parallel)

    !TODO - Check the indexing above. Extend to open boundaries

  end subroutine compute_edge_gradients

!****************************************************

  subroutine compute_plume_velocity(&
       nx,           ny,       &
       dx,           dy,       &
       itest, jtest, rtest,    &
       parallel,               &
       plume_mask,             &
       edge_mask_east,         &
       edge_mask_north,        &
!       edge_mask_east_reduce_v,  &
!       edge_mask_north_reduce_u, &
       D_plume,                &
       H_cavity,               &
       drho_plume,             &
       ddrho_plume_dx_east,    &
       ddrho_plume_dy_east,    &
       ddrho_plume_dx_north,   &
       ddrho_plume_dy_north,   &
       dlsrf_plume_dx_east,    &
       dlsrf_plume_dy_east,    &
       dlsrf_plume_dx_north,   &
       dlsrf_plume_dy_north,   &
       u_plume_east,           &
       v_plume_east,           &
       u_plume_north,          &
       v_plume_north)

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask,            & ! = 1 for cells where scalar plume variables are computed
         edge_mask_east,        & ! = 1 on east edges where plume velocity is computed
         edge_mask_north          ! = 1 on north edges where plume velocity is computed

    real(dp), dimension(nx,ny), intent(in) ::  &
         D_plume,               & ! plume thickness (m)
         H_cavity,              & ! thickness of ocean cavity beneath the plume (m)
         drho_plume,            & ! density difference between plume and ambient ocean (kg/m^3)
         ddrho_plume_dx_east,   & ! horizontal gradient of drho_plume on east edges
         ddrho_plume_dy_east,   & !
         ddrho_plume_dx_north,  & ! horizontal gradient of drho_plume on north edges
         ddrho_plume_dy_north,  & !
         dlsrf_plume_dx_east,   & ! horizontal gradient of lsrf_plume on east edges
         dlsrf_plume_dy_east,   & !
         dlsrf_plume_dx_north,  & ! horizontal gradient of lsrf_plume on north edges
         dlsrf_plume_dy_north     !

!    real(dp), dimension(nx,ny), intent(in) ::  &
!         edge_mask_east_reduce_v,  & ! mask for reducing v on east edges adjacent to a wall
!         edge_mask_north_reduce_u    ! mask for reducing u on north edges adjacent to a wall


    real(dp), dimension(nx,ny), intent(out) ::  &
         u_plume_east,        & ! u_plume on east edges
         v_plume_east,        & ! v_plume on east edges
         u_plume_north,       & ! u_plume on north edges
         v_plume_north          ! v_plume on north edges

    ! local variables

    real(dp), dimension(nx,ny) :: &
         drhox,               & ! density gradient term of pgf_x
         drhoy,               & ! density gradient term of pgf_y
         dsrfx,               & ! surface gradient term of pgf_x
         dsrfy,               & ! surface gradient term of pgf_y
         pgf_x_east,          & ! x component of pressure gradient force on east edges (m^2/s^2)
         pgf_y_east,          & ! y component of pressure gradient force on east edges (m^2/s^2)
         pgf_x_north,         & ! x component of pressure gradient force on north edges (m^2/s^2)
         pgf_y_north            ! y component of pressure gradient force on north edges (m^2/s^2)

!    real(dp), dimension(nx,ny) :: &
!         latdrag_x_east,      & ! x component of lateral drag on east edges (m^2/s^2)
!         latdrag_y_east,      & ! y component of lateral drag on east edges (m^2/s^2)
!         latdrag_x_north,     & ! x component of lateral drag on north edges (m^2/s^2)
!         latdrag_y_north        ! y component of lateral drag on north edges (m^2/s^2)

    real(dp), dimension(nx,ny) :: &
         D_plume_east,        & ! D_plume averaged to east edge
         D_plume_north          ! D_plume averaged to north edge

    integer :: i, j, ig, jg

    integer :: iter_velo        ! iteration counter

    real(dp) :: grav_reduced    ! reduced gravity

    character(len=100) :: message

    logical, dimension(nx,ny) ::  &
         converged_velo_east, & ! true when velocity has converged at an east edge, else false
         converged_velo_north   ! true when velocity has converged at a north edge, else false

    integer :: &
         count_east, count_north  ! number of cells not converged on each face

    integer, parameter ::  &
         maxiter_velo = 30     ! max number of iterations of velocity loop

    ! initialize

    D_plume_east = 0.0d0
    D_plume_north = 0.0d0

    u_plume_east = 0.0d0
    v_plume_east = 0.0d0

    u_plume_north = 0.0d0
    v_plume_north = 0.0d0

    !-------------------------------------------------------------------
    ! Compute the pressure gradient force on each edge, following Lambert et al. (2023):
    ! (1) pgf_x = -(g*D^2)/(2*rhoo) d/dx(drho_plume) + g'*D d/dx(zb - D)
    ! (2) pgf_y = -(g*D^2)/(2*rhoo) d/dy(drho_plume) + g'*D d/dy(zb - D)
    !
    ! where zb = lower plume surface
    !       g' = g * drho_plume/rhoo
    !-------------------------------------------------------------------

    drhox = 0.0d0
    drhoy = 0.0d0
    dsrfx = 0.0d0
    dsrfy = 0.0d0
    pgf_x_east = 0.0d0
    pgf_y_east = 0.0d0
    pgf_x_north = 0.0d0
    pgf_y_north = 0.0d0

    ! PGF on east edges
    ! Loop over all edges of locally owned cells (includes west halo cells)
    do j = nhalo+1, ny-nhalo
       do i = nhalo, nx-nhalo
          if (edge_mask_east(i,j) == 1) then
             D_plume_east(i,j) = (D_plume(i,j) + D_plume(i+1,j)) / 2.0d0
             ! terms proportional to gradients of drho_plume
             drhox(i,j) = -0.5d0*(grav/rhoo) * D_plume_east(i,j)**2 * ddrho_plume_dx_east(i,j)
             drhoy(i,j) = -0.5d0*(grav/rhoo) * D_plume_east(i,j)**2 * ddrho_plume_dy_east(i,j)
             !WHL - check sign
!             drhox(i,j) = 0.5d0*(grav/rhoo) * D_plume_east(i,j)**2 * ddrho_plume_dx_east(i,j)
!             drhoy(i,j) = 0.5d0*(grav/rhoo) * D_plume_east(i,j)**2 * ddrho_plume_dy_east(i,j)
             ! terms proportional to gradients of lsrf_plume
             grav_reduced = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i+1,j)) / 2.0d0
             dsrfx(i,j) = grav_reduced * D_plume_east(i,j) * dlsrf_plume_dx_east(i,j)
             dsrfy(i,j) = grav_reduced * D_plume_east(i,j) * dlsrf_plume_dy_east(i,j)
             pgf_x_east(i,j) = drhox(i,j) + dsrfx(i,j)
             pgf_y_east(i,j) = drhoy(i,j) + dsrfy(i,j)
          endif
       enddo
    enddo

    if (verbose_plume) then
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'PGF components on east edges:'
       endif
       call point_diag(1.d5*drhox, '10^5*density gradient x term', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*dsrfx, '10^5*surface gradient x term', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*drhoy, '10^5*density gradient y term', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*dsrfy, '10^5*surface gradient y term', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*pgf_x_east, '10^5*pgf_x_east (m2/s2)', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*pgf_y_east, '10^5*pgf_y_east (m2/s2)', itest, jtest, rtest, wx, wy)
    endif

    ! PGF on north edges
    ! Loop over all edges of locally owned cells (includes south halo cells)
    do j = nhalo, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (edge_mask_north(i,j) == 1) then
             D_plume_north(i,j) = (D_plume(i,j) + D_plume(i,j+1)) / 2.0d0
             ! terms proportional to gradients of drho_plume
             drhox(i,j) = -0.5d0*(grav/rhoo) * D_plume_north(i,j)**2 * ddrho_plume_dx_north(i,j)
             drhoy(i,j) = -0.5d0*(grav/rhoo) * D_plume_north(i,j)**2 * ddrho_plume_dy_north(i,j)
             !WHL - check sign
!             drhox(i,j) = 0.5d0*(grav/rhoo) * D_plume_north(i,j)**2 * ddrho_plume_dx_north(i,j)
!             drhoy(i,j) = 0.5d0*(grav/rhoo) * D_plume_north(i,j)**2 * ddrho_plume_dy_north(i,j)
             ! terms proportional to gradients of lsrf_plume
             grav_reduced = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i,j+1)) / 2.0d0
             dsrfx(i,j) = grav_reduced * D_plume_north(i,j) * dlsrf_plume_dx_north(i,j)
             dsrfy(i,j) = grav_reduced * D_plume_north(i,j) * dlsrf_plume_dy_north(i,j)
             pgf_x_north(i,j) = drhox(i,j) + dsrfx(i,j)
             pgf_y_north(i,j) = drhoy(i,j) + dsrfy(i,j)
          endif   ! edge_mask_north
       enddo  ! i
    enddo  ! j

    if (verbose_plume) then
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'PGF components on north edges:'
       endif
       call point_diag(1.d5*drhox, '10^5*density gradient x term', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*dsrfx, '10^5*surface gradient x term', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*drhoy, '10^5*density gradient y term', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*dsrfy, '10^5*surface gradient y term', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*pgf_x_north, '10^5*pgf_x_north (m2/s2)', itest, jtest, rtest, wx, wy)
       call point_diag(1.d5*pgf_y_north, '10^5*pgf_y_north (m2/s2)', itest, jtest, rtest, wx, wy)
    endif

    ! initialize other fields
!    latdrag_x_east(:,:) = 0.0d0
!    latdrag_y_east(:,:) = 0.0d0
!    latdrag_x_north(:,:) = 0.0d0
!    latdrag_y_north(:,:) = 0.0d0

    converged_velo_east = .false.
    converged_velo_north = .false.

    ! Iterate as needed to compute a converged velocity at each edge
    !TODO - Keep iterating converged cells to further improve convergence?

    do iter_velo = 1, maxiter_velo

       ! Compute velocity on east edges

       if (verbose_plume .and. main_task) then
          write(iulog,*) ' '
          write(iulog,*) 'iter_velo =', iter_velo
          write(iulog,*) 'compute east edge velocities'
       endif

       call plume_velocity(&
            nx,    ny,               &
            itest, jtest, rtest,     &
            edge_mask_east,          &
            D_plume_east,            &
            pgf_x_east,              &
            pgf_y_east,              &
!            latdrag_x_east,          &
!            latdrag_y_east,          &
            u_plume_east,            &
            v_plume_east,            &
            converged_velo_east)

       ! Compute velocity on north edges

       if (verbose_plume .and. this_rank == rtest) then
          write(iulog,*) 'compute north edge velocities'
       endif

       call plume_velocity(&
            nx,    ny,                  &
            itest, jtest, rtest,        &
            edge_mask_north,            &
            D_plume_north,              &
            pgf_x_north,                &
            pgf_y_north,                &
!            latdrag_x_north,            &
!            latdrag_y_north,            &
            u_plume_north,              &
            v_plume_north,              &
            converged_velo_north)

       ! check for convergence in all cells

       count_east = 0
       count_north = 0

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (edge_mask_east(i,j) == 1 .and. .not.converged_velo_east(i,j) ) then
                count_east = count_east + 1
             endif
             if (edge_mask_north(i,j) == 1 .and. .not.converged_velo_north(i,j) ) then
                count_north = count_north + 1
                if (iter_velo > 15) then
                   call parallel_globalindex(i, j, ig, jg, parallel)
                   write(iulog,*) 'Not converged: ig, jg =', ig, jg
                endif
             endif
          enddo
       enddo

       count_east = parallel_reduce_sum(count_east)
       count_north = parallel_reduce_sum(count_north)

       if (count_east == 0 .and. count_north == 0) then
          if (verbose_plume .and. main_task) write(iulog,*) 'Plume velocity converged'
          exit   ! iter_velo loop
       elseif (iter_velo == maxiter_velo) then
          write(message,*) 'Error, glissade_plume: velocity has not converged, iter_velo =', iter_velo
          call write_log(message, GM_FATAL)
       elseif (verbose_plume .and. main_task) then
          write(iulog,*) 'Velocity not converged: count_east, count_north =', count_east, count_north
       endif

    enddo  ! iter_velo

    ! Extrapolate the final velocity to open boundaries (edge_mask = 2)

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          if (edge_mask_east(i,j) == 2) then
             if (plume_mask(i,j) == 1) then
                ! cell (i,j) is a plume cell, but cell (i+1,j) is open water;
                ! assign this edge the same velocity as the adjacent interior edge
                if (edge_mask_east(i-1,j) == 1) then
                   u_plume_east(i,j) = u_plume_east(i-1,j)
                   v_plume_east(i,j) = v_plume_east(i-1,j)
                endif
             elseif (plume_mask(i+1,j) == 1) then
                ! cell (i+1,j) is a plume cell, but cell (i,j) is open water;
                ! assign this edge the same velocity as the adjacent interior edge
                if (edge_mask_east(i+1,j) == 1) then
                   u_plume_east(i,j) = u_plume_east(i+1,j)
                   v_plume_east(i,j) = v_plume_east(i+1,j)
                endif
             endif
          endif   ! edge_mask_east = 2

          if (edge_mask_north(i,j) == 2) then
             if (plume_mask(i,j) == 1) then
                ! cell (i,j) is a plume cell, but cell (i,j+1) is open water;
                ! assign this edge the same velocity as the adjacent interior edge
                if (edge_mask_north(i,j-1) == 1) then
                   u_plume_north(i,j) = u_plume_north(i,j-1)
                   v_plume_north(i,j) = v_plume_north(i,j-1)
                endif
             elseif (plume_mask(i,j+1) == 1) then
                ! cell (i,j+1) is a plume cell, but cell (i,j) is open water;
                ! assign this edge the same velocity as the adjacent interior edge
                if (edge_mask_north(i,j+1) == 1) then
                   u_plume_north(i,j) = u_plume_north(i,j+1)
                   v_plume_north(i,j) = v_plume_north(i,j+1)
                endif
             endif
          endif   ! edge_mask_north = 2

       enddo  ! i
    enddo   ! j

    call parallel_halo(u_plume_east, parallel)
    call parallel_halo(v_plume_east, parallel)
    call parallel_halo(u_plume_north, parallel)
    call parallel_halo(v_plume_north, parallel)

    if (verbose_plume) then
       call point_diag(u_plume_east, 'u_plume_east (m/s)', itest, jtest, rtest, wx, wy)
       call point_diag(v_plume_east, 'v_plume_east (m/s)', itest, jtest, rtest, wx, wy)
       call point_diag(u_plume_north, 'u_plume_north (m/s)', itest, jtest, rtest, wx, wy)
       call point_diag(v_plume_north, 'v_plume_north (m/s)', itest, jtest, rtest, wx, wy)
    endif

    !TODO - Remove the lateral drag calculation? I don't remember why I added it

    !TODO - Now that the velocity has converged without lateral drag, try adding the lateral drag
    !       terms and recomputing the velocity. Not sure how to do this stably.

    ! Compute the lateral drag term based on the current guess for the velocity

!       call compute_lateral_drag(&
!            nx,         ny,      &
!            dx,         dy,      &
!            itest, jtest, rtest, &
!            edge_mask_east,      &  !TODO - divu_mask or edge_mask?
!            edge_mask_north,     &
!            plume_mask_cell,     &
!            D_plume,             &
!            u_plume_east,        &
!            v_plume_east,        &
!            u_plume_north,       &
!            v_plume_north,       &
!            latdrag_x_east,      &
!            latdrag_y_east,      &
!            latdrag_x_north,     &
!            latdrag_y_north)

!       if (verbose_plume .and. main_task .and. this_rank==rtest) then
!          print*, ' '
!          print*, 'Computed lateral drag terms'
!          print*, ' '
!       endif  ! verbose_plume


       !WHL - With new code, the velocity should be computed at these edges, and not extrapolated.
       !      Extrapolation can make it hard to have divergence/convergence.

       ! Extrapolate the velocity to open edges (plume on one side, open water on the other)
       !  This extrapolation is not expected to be accurate, but it prevents large convergence
       !  in cells adjacent to water.
       ! If the plume exists on neither side of the edge, the velocity remains set to zero.
       ! Also, u_plume_east = 0 on global E and W boundaries, and v_plume_north = 0 on global N and S boundaries.
       !  This prevents outflow through domain walls.
       !  Along the upper ("northern") boundary of the ISOMIP+ domain, the flow is forced to form an eastward jet.

       !TODO - Are global_bndy masks needed here? Wondering if we can avoid passing in 4 global_bndy fields.

!    do j = nhalo, ny-nhalo
!       do i = nhalo, nx-nhalo

          ! east edges
!          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 0 .and. global_bndy_east(i,j) == 0) then
!             if (lsrf(i+1,j) == 0.0d0 .or. floating_mask(i+1,j) == 1) then
                ! water in cell (i+1,j); get plume velocity from edge (i-1,j)
!                u_plume_east(i,j) = u_plume_east(i-1,j)
!             endif
!          elseif (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i+1,j) == 1 .and. global_bndy_west(i+1,j) == 0) then
!             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                ! water in cell (i,j); get plume velocity from edge (i+1,j)
!                u_plume_east(i,j) = u_plume_east(i+1,j)
!             endif
!          endif

          ! north edges
!          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 0 .and. global_bndy_north(i,j) == 0) then
!             if (lsrf(i,j+1) == 0.0d0 .or. floating_mask(i,j+1) == 1) then
                ! water in cell (i,j+1); get plume velocity from edge (i,j-1)
!                v_plume_north(i,j) = v_plume_north(i,j-1)
!             endif
!          elseif (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i,j+1) == 1 .and. global_bndy_south(i,j+1) == 0) then
!             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                ! water in cell (i,j); get plume velocity from edge (i,j+1)
!                v_plume_north(i,j) = v_plume_north(i,j+1)
!             endif
!          endif
!       enddo   ! i
!    enddo   ! j

  end subroutine compute_plume_velocity

!****************************************************

  !TODO - Remove latdrag input terms?
  subroutine plume_velocity(&
       nx,    ny,               &
       itest, jtest, rtest,     &
       edge_mask,               &
       D_plume,                 &
       pgf_x,                   &
       pgf_y,                   &
!       latdrag_x,               &
!       latdrag_y,               &
       u_plume,                 &
       v_plume,                 &
       converged_velo)
    
    ! Compute the velocity on a set of edges (either east or north)

    integer, intent(in) ::  &
         nx,  ny,           & ! number of grid cells in each dimension
         itest, jtest, rtest  ! test cell coordinates (diagnostic only)
    
    integer, dimension(nx,ny), intent(in) ::   &
         edge_mask            ! = 1 at edges where velocity is computed

    ! Note: The following variables are co-located with the velocity
    real(dp), dimension(nx,ny), intent(in) ::   &
         D_plume,           & ! plume thickness at edges (m)
         pgf_x,             & ! x component of pressure gradient force
         pgf_y                ! y component of pressure gradient force
!         latdrag_x,         & ! x component of lateral drag
!         latdrag_y            ! y component of lateral drag

    ! Note: u and v are colocated on either east edges or north edges,
    !       depending on the subroutine call; speed lies on the same edge
    real(dp), dimension(nx,ny), intent(inout) ::  &
         u_plume,           & ! x component of plume velocity (m/s) on the edge
         v_plume              ! y component of plume velocity (m/s) on the edge

!!    logical, dimension(nx,ny), intent(inout) ::  &
    logical, dimension(nx,ny), intent(out) ::  &
         converged_velo       ! true when velocity has converged at an edge, else false

    ! local variables

!    real(dp), dimension(nx,ny) ::   &
!         f_x,               &  ! pgf_x + latdrag_x
!         f_y                   ! pgf_y + latdrag_y

!    real(dp), dimension(nx,ny) ::  &
!         reduce_v,          &  ! local version of edge_mask_east_reduce_v; no reduction by default
!         reduce_u              ! local version of edge_mask_north_reduce_u; no reduction by default

    real(dp) :: &
         speed,             & ! plume speed (m/s), updated at each iteration until convergence
         x_resid, y_resid,  & ! residuals of momentum balance equations (m^2/s^2)
         denom,             & ! denominator
         a_uu, a_uv,        & ! coefficients for Newton solve
         a_vu, a_vv,        & !
         du, dv               ! change in u_plume and v_plume (m/s)
    
    character(len=128) :: message

    real(dp), parameter :: &
         maxresid_force_balance = 1.0d-8   ! max residual allowed in momentum balance equation (m^2/s^2)

    !TODO - Start with Picard, then test Newton
    logical, parameter :: &
         velo_newton = .true.  ! if true, use Newton's method; if false, use Picard method
!         velo_newton = .false.  ! if true, use Newton's method; if false, use Picard method

    integer :: i, j

    logical, parameter :: verbose_velo = .false.

    !--------------------------------------------------------------------
    ! Compute the plume velocity.
    ! Assume a balance between the pressure gradient force, basal drag and Coriolis:
    !
    ! pgf_x - c_d*|U|*u + D*f*v = 0
    ! pgf_y - c_d*|U|*v - D*f*u = 0
    !
    !          pgf = pressure gradient force
    !            D = plume thickness
    !          c_d = dimensionless ocean drag coefficient
    !            f = Coriolis coefficient
    !          |U| = sqrt(u^2 + v^2 + u_tidal^2)
    !      u_tidal = a small velocity added for regularization
    !
    ! The solution (assuming D is known) is
    !
    !                c_d*|U|*pgf_x + D*f*pgf_y
    !            u = ________________________
    !                 (D*f)^2 + (c_d*|U|)^2
    !
    !                c_d*|U|*pgf_y - D*f*pgf_x 
    !            v = ________________________
    !                 (D*f)^2 + (c_d*|U|)^2
    !
    ! Since |U| is a function of u and v, we iterate to convergence.
    !
    ! The iteration can be sped up by using Newton's method.
    ! We write   u = u0 + du
    !            v = v0 + dv
    !          |U| = U0 + d|U|/du * du + d|U|dv * dv
    ! where the partial derivatives are evaluated at (u,v) = (u0,v0).
    ! 
    ! This gives 
    !           du = (a_vv * R_x - a_uv * R_y) / det|A|
    !           dv = (a_uu * R_y - a_vu * R_x) / det|A|
    ! where    
    !          R_x = pgf_x - c_d*U0*u0 + D*f*v0 = x residual
    !          R_y = pgf_y - c_d*U0*v0 - D*f*u0 = y residual
    !
    !                | a_uu   a_uv |
    ! and        A = |             |     
    !                | a_vu   a_vv |
    !
    ! with    a_uu = c_d*(U0 + u0^2/U0)
    !         a_uv = c_d*u0*v0/U0 - D*f) 
    !         a_vu = c_d*u0*v0/U0 + D*f) 
    !         a_vv = c_d*(U0 + v0^2/U0) 
    !
    !--------------------------------------------------------------------

!    if (present(edge_mask_north_reduce_u)) then
!       reduce_u(:,:) = edge_mask_north_reduce_u(:,:)
!    else
!       reduce_u(:,:) = 1.0d0  ! no reduction
!    endif

!    if (present(edge_mask_east_reduce_v)) then
!       reduce_v(:,:) = edge_mask_east_reduce_v(:,:)
!    else
!       reduce_v(:,:) = 1.0d0  ! no reduction
!    endif

    ! Combine PGF and lateral drag into one term
!    f_x(:,:) = pgf_x(:,:) + latdrag_x(:,:)
!    f_y(:,:) = pgf_y(:,:) + latdrag_y(:,:)

    ! Loop over edges
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! Compute plume speed based on current u and v (now passed in directly)
          speed = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)

          if (edge_mask(i,j) == 1) then
       
             ! Compute residual of the momentum balance
             x_resid = pgf_x(i,j) - c_drag*speed*u_plume(i,j) + f_coriolis*D_plume(i,j)*v_plume(i,j)
             y_resid = pgf_y(i,j) - c_drag*speed*v_plume(i,j) - f_coriolis*D_plume(i,j)*u_plume(i,j)

             ! check convergence of plume velocity
             if (abs(x_resid) < maxresid_force_balance .and. abs(y_resid) < maxresid_force_balance) then
                converged_velo(i,j) = .true.
             endif

             ! Should be no harm to compute anyway, even if already converged. (Verify this)
!!             if (.not.converged_velo(i,j)) then

             if (velo_newton) then
          
                ! compute some coefficients for the Newton solve
                a_uu = c_drag * (speed + u_plume(i,j)**2/speed)
                a_vv = c_drag * (speed + v_plume(i,j)**2/speed)
                      
                a_uv = c_drag * (u_plume(i,j)*v_plume(i,j))/speed - D_plume(i,j)*f_coriolis
                a_vu = c_drag * (u_plume(i,j)*v_plume(i,j))/speed + D_plume(i,j)*f_coriolis
!                   a_uv = c_drag * (u_plume(i,j)*v_plume(i,j))/speed - reduce_v(i,j)*D_plume(i,j)*f_coriolis
!                   a_vu = c_drag * (u_plume(i,j)*v_plume(i,j))/speed + reduce_u(i,j)*D_plume(i,j)*f_coriolis
                   
                ! compute du and dv
                denom = a_uu*a_vv - a_uv*a_vu
                      
                if (abs(denom) > 0.0d0) then
                   du = (a_vv*x_resid - a_uv*y_resid) / denom
                   dv = (a_uu*y_resid - a_vu*x_resid) / denom
                         
                   u_plume(i,j) = u_plume(i,j) + du
                   v_plume(i,j) = v_plume(i,j) + dv
                      
                else  ! denom = 0.0
                   write(iulog,*) 'Error, glissade_plume: ill-posed Newton solve for velocity, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'a_uu, a_vv, a_uv, a_vu =', a_uu, a_vv, a_uv, a_vu
                   write(message,*) 'Error, glissade_plume: ill-posed Newton solve for velocity, rank, i, j:', this_rank, i, j
                   call write_log(message, GM_FATAL)
                endif
                      
             else  ! simpler Picard solve
          
                denom = (c_drag*speed)**2 + (D_plume(i,j)*f_coriolis)**2
                u_plume(i,j) = (c_drag*speed*pgf_x(i,j) + f_coriolis*D_plume(i,j)*pgf_y(i,j)) / denom
                v_plume(i,j) = (c_drag*speed*pgf_y(i,j) - f_coriolis*D_plume(i,j)*pgf_x(i,j)) / denom
          
             endif  ! Newton or Picard

!!             endif  ! .not.converged_velo

             if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) 'speed (m/s) =', speed
                write(iulog,*) 'pgf_x, pgf_y:', pgf_x(i,j), pgf_y(i,j)
!                write(iulog,*) 'latdrag_x, latdrag_y:', latdrag_x(i,j), latdrag_y(i,j)
                write(iulog,*) 'Dfv, -Dfu:', D_plume(i,j) * f_coriolis * v_plume(i,j), &
                                     -D_plume(i,j) * f_coriolis * u_plume(i,j)
                write(iulog,*) 'dragu, dragv:', c_drag * speed * u_plume(i,j), &
                                         c_drag * speed * v_plume(i,j)
                write(iulog,*) 'x/y residual:', x_resid, y_resid
                write(iulog,*) 'new u/v_plume:', u_plume(i,j), v_plume(i,j)
                write(iulog,*) 'converged =', converged_velo(i,j)
             endif

          endif  ! edge_mask
       enddo  ! i
    enddo  ! j

  end subroutine plume_velocity

!****************************************************

  subroutine plume_entrainment_gaspar(&
       nx,           ny,      &
       dx,           dy,      &
       itest, jtest, rtest,   &
       parallel,              &
       plume_mask,            &
       ustar_plume,           &
       bmlt_float,            &
       drho_plume,            &
       drho_basal,            &
       H_cavity,              &
       D_plume,               &
       dt_plume,              &
       entrainment,           &
       detrainment)

    !--------------------------------------------------------------------
    ! Compute entrainment as a function of the friction velocity and plume thickness,
    ! following Gaspar (1988), Gladish (2012) and Lambert et al. (2023).
    !
    !       (D/2)*gb'*m + (D/2)*ga'*e = mu*(u*)^3
    !
    !       where m = melt, e = entrainment, u* = friction velocity, mu = nondim parameter
    !             ga'= (grav/rhoo)*drho_plume, gb' = (grav/rhoo)*drho_basal
    !
    ! Rearrange to get e = [mu*(u*)^3 - (D/2)*gb'*m] / [(D/2)*ga']
    !
    ! Can have e < 0 for small u* and/or large m. If so, then classify as detrainment.
    !--------------------------------------------------------------------

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    type(parallel_type), intent(in) :: &
         parallel                    ! info for parallel communication

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask             ! = 1 for cells where scalar plume variables are computed

    real(dp), dimension(nx,ny), intent(in) ::  &
         ustar_plume,         & ! friction velocity (m/s)
         bmlt_float,          & ! melt rate (m/s)
         drho_plume,          & ! density difference between ambient ocean and plume (kg/m3)
         drho_basal,          & ! density difference between plume and ice base (kg/m3)
         H_cavity,            & ! cavity thickness (m)
         D_plume                ! plume thickness (m)

    real(dp), intent(in) :: &
         dt_plume               ! timestep (s)

    real(dp), dimension(nx,ny), intent(out) ::  &
         entrainment,         & ! entrainment at cell centers (m/s)
         detrainment            ! detrainment at cell centers (m/s)

    ! local variables

    real(dp) :: &
         Dmax,                & ! max plume thickness = min(D_plume_max, H_cavity)
         entrainment_min,     & ! min entrainment rate if D_plume < D_plume_min
         detrainment_min        ! min detrainment rate if D_plume > D_plume_max

    real(dp) :: numer, denom
    integer :: i, j, ig, jg

    ! entrainment parameters
    real(dp), parameter ::    &
         mu_e = 2.5d0           ! nondimensional parameter
                                 ! Gaspar (1988) and Lambert et al. (2023) set mu = 0.5;
                                 ! Gladish et al. (2012) and Lambert et al. (2026) set mu = 2.5

    entrainment = 0.0d0
    detrainment = 0.0d0

    ! loop over all cells
    do j = 1, ny
       do i = 1, nx
          if (plume_mask(i,j) == 1) then

             numer = mu_e * ustar_plume(i,j)**3 - 0.5d0*D_plume(i,j)*(grav/rhoo)*drho_basal(i,j)*bmlt_float(i,j)
             denom = 0.5d0*D_plume(i,j)*(grav/rhoo)*drho_plume(i,j)
             if (denom > 0.0d0) then
                entrainment(i,j) = numer/denom
             else   ! likely have drho_plume = 0
                entrainment(i,j) = 0.0d0
             endif

             if (entrainment(i,j) < 0.0d0) then
                detrainment(i,j) = -1.0d0*entrainment(i,j)
                entrainment(i,j) = 0.0d0
             endif

             ! Increase entrainment if D_plume < D_plume_min
             if (D_plume(i,j) < D_plume_min) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'Force entrainment: ig, jg, D_plume:', ig, jg, D_plume(i,j)
                entrainment_min = (D_plume_min - D_plume(i,j)) / tau_relax_entrainment
                entrainment(i,j) = max(entrainment(i,j), entrainment_min)
             endif

             ! Increase detrainment if D_plume > D_plume_max or H_cavity
             !WHL - Don't use the H_cavity limit, just use D_plume_max
!!             Dmax = min(H_cavity(i,j), D_plume_max)
             Dmax = D_plume_max
             if (D_plume(i,j) > Dmax) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'Force detrainment: ig, jg, D_plume:', ig, jg, D_plume(i,j)
                detrainment_min = (D_plume(i,j) - Dmax) / tau_relax_entrainment   ! < 0
                detrainment(i,j) = max(detrainment(i,j), detrainment_min)
             endif

          endif
       enddo   ! i
    enddo   ! j

  end subroutine plume_entrainment_gaspar

!****************************************************

  subroutine plume_entrainment(&
       nx,         ny,      &
       dx,         dy,      &
       itest, jtest, rtest, &
       parallel,            &
       plume_mask,          &
       theta_slope,         &
       plume_speed,         &
       H_cavity,            &
       D_plume,             &
       dt_plume,            &
       entrainment,         &
       detrainment)

    !--------------------------------------------------------------------
    ! Compute entrainment as a function of the plume speed and the slope of the
    !  plume-ambient interface, following Bo Pederson (1980) and Jenkins (1991):
    !
    ! entrainment = E0 * plume_speed * sin(theta_slope)
    !
    ! Note: plume_speed is proportional to ustar_plume, so we could replace
    !       one with the other and rescale the constant.
    !--------------------------------------------------------------------

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask             ! = 1 for cells where scalar plume variables are computed

    real(dp), dimension(nx,ny), intent(in) ::  &
         theta_slope,         & ! basal slope angle at cell centers (rad)
         plume_speed,         & ! plume speed at cell center (m/s)
         H_cavity,            & ! thickness of sub-shelf cavity (m)
         D_plume                ! plume thickness
    
    real(dp), intent(in) :: &
         dt_plume               ! timestep (s)

    !Note: Both entrainment and detrainment are >= 0 by definition
    real(dp), dimension(nx,ny), intent(out) ::  &
         entrainment,           & ! entrainment at cell centers (m/s)
         detrainment              ! detrainment at cell centers (m/s)

    ! local variables

    real(dp) :: &
         Dmax,                  & ! max plume thickness = min(D_plume_max, H_cavity)
         entrainment_min,       & ! min entrainment rate when D_plume < D_plume_min
         detrainment_min          ! min detrainment rate when D_plume > D_plume_max

    integer :: i, j, ig, jg

    ! entrainment parameters
    real(dp), parameter ::   &
!!         E0 = 0.072d0                  ! entrainment coefficient (unitless)
         E0 = 0.036d0                   ! entrainment coefficient (unitless)   ! trying a smaller value
                                       ! Bo Pederson (1980) suggests E0 = 0.072
                                       ! Jenkins (1991, JGR) suggests 0.036 to compensate for lack of Coriolis in 1D model

    entrainment = 0.0d0
    detrainment = 0.0d0

    ! loop over all cells
    do j = 1, ny
       do i = 1, nx
          if (plume_mask(i,j) == 1) then

             entrainment(i,j) = E0 * plume_speed(i,j) * sin(theta_slope(i,j))

             ! Impose a minimum entrainment rate if D_plume < D_plume_min
             if (D_plume(i,j) < D_plume_min) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'Force entrainment: ig, jg, D_plume:', ig, jg, D_plume(i,j)
                entrainment_min = (D_plume_min - D_plume(i,j)) / tau_relax_entrainment
                entrainment(i,j) = max(entrainment(i,j), entrainment_min)
             endif

             ! Impose a minimum detrainment rate if D_plume > D_plume_max
             ! Note: Arguably, D_plume should not exceed H_cavity either; test this
!!             Dmax = min(H_cavity(i,j), D_plume_max)
             Dmax = D_plume_max
             if (D_plume(i,j) > Dmax) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'Force detrainment: ig, jg, D_plume:', ig, jg, D_plume(i,j)
                detrainment_min = (D_plume(i,j) - Dmax) / tau_relax_entrainment   ! < 0
                detrainment(i,j) = max(detrainment(i,j), detrainment_min)
             endif

          endif
       enddo
    enddo

  end subroutine plume_entrainment

!****************************************************

  subroutine plume_melt_rate(&
       nx,         ny,      &
       itest, jtest, rtest, &
       parallel,            &
       plume_mask,          &
       gammaT,              &
       gammaS,              &
       pressure,            &
       ustar_plume,         &
       D_plume,             &
       T_plume,             &
       S_plume,             &
       T_basal,             &
       S_basal,             &
       bmlt_float)

    !--------------------------------------------------------------------
    ! Compute the melt rate at the ice-ocean interface.
    !
    ! Following Jenkins et al. (2010) and Asay-Davis et al. (2016),
    ! there are 3 equations for the 3 unknowns m, Tb and Sb,
    ! where m = melt rate at ice-ocean interface
    !       Tb = potential temperature at ice-ocean interface
    !       Sb = salinity at ice-ocean interface
    ! 
    ! (1) rhoi * m * L  = rhoo * cpw * u_fric * GammaT * (Tp - Tb)
    ! (2) rhoi * m * Sb = rhoo * u_fric * GammaS *(Sp - Sb)
    ! (3) Tb = lambda1*Sb + lambda2 + lambda3*pb 
    !
    ! Eqs. 1 and 2 describe heat and salt transfer at the ice-ocean interface.
    ! Eq. 3 is the linearized liquidus relation that determines the potential freezing point.
    ! Note: Asay-Davis et al. use rhow instead of rhoi on the LHS, since they define
    !       the melt rate m in units of meters of freshwater instead of meters of ice.
    !       See their Sec. 3.1.8.
    !
    ! We can rewrite these equations as
    !
    ! (1)     m = C1 * (Tp - Tb)
    ! (2)  m*Sb = C2 * (Sp - Sb)
    ! (3)    Tb = lambda1*Sb + C3
    !
    ! where C1 = (rhoo * cpw * ufric * GammaT) / (rhoi * L)
    !       C2 = (rhoo * ufric * GammaS) / rhoi
    !       C3 = lambda2 + lambda3*pb
    !
    ! Use (3) to substitute for Tb in (1): m = C1 * [Tp - lambda1*Sb - C3)
    !
    ! Then substitute for m in (2): C1*[Tp - lambda1*Sb - C3) * Sb = C2*(Sp - Sb)
    !
    ! Rearrange terms: (-lambda1*C1)*Sb^2 + [C1(Tp - C3) + C2]*Sb - C2*Sp = 0
    !
    ! Multiply by -1: (lambda1*C1)*Sb^2 + [C1(C3 - Tp) - C2]*Sb + C2*Sp = 0
    !
    ! This is a quadratic equation for Sb. Solve using the quadratic formula,
    !  then substitute to get m and Tb.
    !
    ! Note: This treatment assumes that GammaT and GammaS are spatially uniform constants.
    !       Lambert et al. (2023) have the following instead:
    !       (1) m * L = cpw * gammaT * (Tp - Tb)
    !       (2) m * Sb = gammaS * (Sp - Sb)
    !       where gammaT = ufric / [2.12d0*log(ufric*D_plume/kvw) + 12.5d0*Prandtl**(2.0d0/3.0d0) - 8.68d0]
    !             gammaS = ufric / [2.12d0*log(ufric*D_plume/kvw) + 12.5d0*Schmidt**(2.0d0/3.0d0) - 8.68d0]
    !             kvw = kinematic viscosity of seawater
    !             Prandtl and Schmidt are dimensionless numbers for turbulent transfer
    !--------------------------------------------------------------------

    ! input/output arguments
    ! Note: lambda1, lambda2, lambda2, c_drag and u_tidal are declared at the top of the module
    
    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    integer, intent(in) ::  &
         itest, jtest, rtest    ! test cell coordinates (diagnostic only)

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    integer, dimension(nx,ny), intent(in) :: &
         plume_mask             ! = 1 for cells where scalar plume variables are computed

    real(dp), intent(in) ::  &
         gammaT,              & ! nondimensional heat transfer coefficient
         gammaS                 ! nondimensional salt transfer coefficient

    real(dp), dimension(nx,ny), intent(in) :: &
         pressure,            & ! ocean pressure at base of ice (N/m^2)
         ustar_plume,         & ! plume friction velocity (m/s) on ice grid, output as a diagnostic
         D_plume,             & ! plume thickness (m)
         T_plume,             & ! plume temperature (deg C)
         S_plume                ! plume salinity (psu)

    real(dp), dimension(nx,ny), intent(out) :: &
         T_basal,             & ! basal ice temperature (deg C)
         S_basal,             & ! basal ice salinity (psu)
         bmlt_float             ! melt rate at base of floating ice (m/s)
    
    ! local variables
    
    real(dp) :: &
         C1, C2, C3,          & ! factors in the three melt-rate equations
         aa, bb, cc,          & ! factors in quadratic formula
         discriminant,        & ! (b^2 - 4ac) term in quadratic formula
         Sb1, Sb2               ! solutions of quadratic formula
    
    integer :: i, j, ig, jg

    logical :: abort            ! if true, then abort

    logical, parameter :: verbose_melt = .false.

    ! initialize
    T_basal = 0.0d0
    S_basal = 0.0d0
    bmlt_float = 0.0d0

    ! Loop over all cells
    do j = 1, ny
       do i = 1, nx
          
          if (plume_mask(i,j) == 1) then

             ! Solve a quadratic equation for S_basal
             C1 = (rhoo * cpw * ustar_plume(i,j) * gammaT) / (rhoi * lhci)
             C2 = (rhoo * ustar_plume(i,j) * gammaS) / rhoi
             C3 = lambda2 + lambda3*pressure(i,j)

             aa = lambda1*C1   ! Note: lambda1 < 0 , so aa < 0
             bb = C1*(C3 - T_plume(i,j)) - C2
             cc = C2*S_plume(i,j)

             abort = .false.
             discriminant = bb**2 - 4.d0*aa*cc
             if (discriminant >= 0.0d0) then
                Sb1 = (-bb + sqrt(discriminant)) / (2.0d0*aa)
                Sb2 = (-bb - sqrt(discriminant)) / (2.0d0*aa)
                if (Sb2 >= 0.0d0 .and. Sb1 <= 0.0d0) then
                   S_basal(i,j) = Sb2
                else
                   abort = .true.
                endif
             else
                abort = .true.
             endif

             if (abort) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'Failed to solve quadratic equation for S_plume, ig, jg =', ig, jg
                write(iulog,*) 'a, b, c =', aa, bb, cc
                write(iulog,*) 'b^2 - 4ac =', bb*bb - 4.0d0*aa*cc
                write(iulog,*) 'Sb1, Sb2 =', Sb1, Sb2
                call write_log('Failed to solve quadratic equation for S_plume', GM_FATAL)
             endif

             ! Solve for T_basal and bmlt_float
             T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
             bmlt_float(i,j) = C1 * (T_plume(i,j) - T_basal(i,j))

             if (verbose_melt .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'Melt rate calc: rank, i, j =', rtest, i, j
                write(iulog,*) 'pressure (Pa) =', pressure(i,j)
                write(iulog,*) 'C1 (m/s/deg), C2 (m/s), C3(deg C):', C1, C2, C3
                write(iulog,*) 'aa, bb, cc:=', aa, bb, cc
                write(iulog,*) 'T_basal, S_basal, bmlt_float:', T_basal(i,j), S_basal(i,j), bmlt_float(i,j)
                write(iulog,*) 'Eq. 1 LHS:', rhoi*lhci*bmlt_float(i,j)
                write(iulog,*) 'Eq. 1 RHS:', rhoo*cpw*ustar_plume(i,j)*gammaT*(T_plume(i,j) - T_basal(i,j))
                write(iulog,*) 'Eq. 2 LHS:', rhoi*bmlt_float(i,j)*S_basal(i,j)
                write(iulog,*) 'Eq. 2 RHS:', rhoo*ustar_plume(i,j)*gammaS*(S_plume(i,j) - S_basal(i,j))
                write(iulog,*) 'Eq. 3 LHS:', T_basal(i,j)
                write(iulog,*) 'Eq. 3 RHS:', lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
             endif
          endif   ! plume_mask = 1
          
       enddo   ! i
    enddo   ! j

  end subroutine plume_melt_rate

!****************************************************
    
  subroutine plume_transport(&
       nx,           ny,     &
       dx,           dy,     &
       itest, jtest, rtest,  &
       parallel,             &
       dt,                   &
       plume_mask,           &
       edge_mask_east,       &
       edge_mask_north,      &
       u_plume_east,         &
       v_plume_north,        &
       entrainment,          &
       detrainment,          &
       bmlt_float,           &
       heat_transfer,        &
       T_ambient,            &
       S_ambient,            &
       T_basal,              &
       S_basal,              &
       D_plume,              &
       T_plume,              &
       S_plume)

    !----------------------------------------------------------------------------
    ! Solve transport equations for the plume thickness, temperature and salinity.
    ! These include horizontal transport of mass, heat and salt; horizontal diffusion
    !  of heat and salt; and vertical entrainment, detrainment and melting.
    !
    ! See Eqs. 1, 2 and 4 in Lambert et al. (2023):. These describe the conservation
    !
    ! (1) dD/dt + del*(DU) = e - d + m
    !
    ! (2) d/dt(DT) + del*(DUT) = e*Ta -d*T + m*Tb - gammaT*(T - Tb) + del*(Kh*D*gradT)
    !
    ! (3) d/dt(DS) + del*(DUS) = e*Sa -d*S - del*(Kh*D*gradT)
    !
    ! where (T,S), (Tb,Sb) and (Ta,Sa) are the temperature and salinity of the plume,
    ! the ice base and the ambient ocean, respectively; D is the plume thickness;
    ! e, d and m are the rates of entraintment, detrainment and melting;
    ! gammaT is a heat transfer term, set here to (rhoi*Li*m)/(rhow*cpw);
    ! and Kh is the diffusivity of heat and salt.
    !
    ! Lambert (2023) does not specify Kh, but the LADDIE code has Kh = 25 m2/s:
    !  https://github.com/erwinlambert/laddie (7/27/26).
    ! At 8km resolution with dt_plume = 600 s, Kh = 50 is stable, but Kh = 100 is not.
    ! Explicit diffusion has a CFL limit proportional to dx^2, so the maximum
    !  stable step drops sharply with increasing resolution.
    ! To turn off diffusion, simply set Kh = 0.
    !
    ! Entrainment and detrainment are treated separately since e is linked to (Ta,Sa)
    !  and d to (Tp,Sp).
    !
    ! Salt transfer between the plume and the ice base is assumed to be negligible,
    !  so there are no terms m*Sb or gammaS*(S - Sb) in Eq. 3.
    !
    ! The horizontal transport equations are solved with a first-order upwind scheme,
    !  to reduce the cost compared to incremental remapping.
    !
    ! The diffusive terms are handled with a finite-difference scheme in flux form,
    !  using the up-gradient value of D.
    !----------------------------------------------------------------------------

    use glissade_transport, only: glissade_upwind_field

    ! input/output arguments

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    real(dp), intent(in) ::  &
         dt                     ! time step (s)

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask,          & ! = 1 for cells where the plume is present, else = 0
         edge_mask_east,      & ! = 1 for east edges with plume cells on each side
         edge_mask_north        ! = 1 for north edges with plume cells on each side

    real(dp), dimension(nx,ny), intent(in) ::  &
         u_plume_east,        & ! u_plume on east edges (m/s)
         v_plume_north,       & ! v_plume on north edges (m/s)
         entrainment,         & ! entrainment rate (m/s)
         detrainment,         & ! detrainment rate (m/s)
         bmlt_float,          & ! basal melt rate (m/s)
         heat_transfer,       & ! rate of heat transfer from plume to ice (J/m^2/s)
         T_ambient,           & ! ambient temperature (deg C)
         S_ambient,           & ! ambient salinity (psu)
         T_basal,             & ! basal temperature (deg C)
         S_basal                ! basal salinity (psu)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         D_plume,             & ! plume thickness (m)
         T_plume,             & ! plume temperature (deg C)
         S_plume                ! plume salinity (psu)

    ! local variables

    integer :: i, j, ig, jg, n
    integer :: ilo, ihi, jlo, jhi

    real(dp) :: dD, dDT, dDS           ! increments in D, D*T and D*S
    real(dp) :: gradT, gradS           ! gradients of heat and salt

    real(dp), dimension(nx,ny,3) :: &
         work                          ! work array for transport

    real(dp), dimension(nx,ny) ::  &
         diffT_east, diffT_north,    & ! diffusive fluxes of heat at cell edges (m^3*deg/s)
         diffS_east, diffS_north       ! diffusive fluxes of salt at cell edges (m^3*psu/s)

    real(dp), dimension(nx,ny) ::  &
         D_temp, T_temp, S_temp        ! temporary arrays

    character(len=100) :: message

    real(dp), parameter :: &
         T_plume_min = -3.0d0,       & ! min allowed T_plume (deg C)
         T_plume_max = 10.0d0,       & ! max allowed T_plume (deg C)
         S_plume_min = 0.0d0,        & ! min allowed S_plume (psu)
         S_plume_max = 40.0d0          ! max allowed S_plume (psu)

    call point_diag(D_plume, 'Starting D_plume (m)', itest, jtest, rtest, wx, wy)
    call point_diag(T_plume, 'T_plume (degC)', itest, jtest, rtest, wx, wy)
    call point_diag(S_plume, 'S_plume (psu)', itest, jtest, rtest, wx, wy)

    ! Make sure the input fields are in range
    !TODO - Leave this out, and only test the output?

    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             ! Note: Allow D_plume > D_plume_max; detrainment should relax toward D_plume_max
             if (D_plume(i,j) < D_plume_min) then
                call parallel_globalindex(i, j, ig, jg, parallel)
!                write(iulog,*) 'Warning, input D_plume < D_plume_min: ig, jg, D_plume =', ig, jg, D_plume(i,j)
!!                write(message,*) 'Plume transport, input D_plume < D_plume_min: ig, jg, D_plume =', &
!!                     ig, jg, D_plume(i,j)
!!                call write_log(message, GM_WARNING)
             elseif (D_plume(i,j) > D_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
!                write(iulog,*) 'Warning, input D_plume > D_plume_max: ig, jg, D_plume =', ig, jg, D_plume(i,j)
!!                write(message,*) 'Plume transport, input D_plume < D_plume_max: ig, jg, D_plume =', &
!!                     ig, jg, D_plume(i,j)
!!                call write_log(message, GM_WARNING)
             endif
             if (T_plume(i,j) < T_plume_min .or. T_plume(i,j) > T_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
!                write(iulog,*) 'Plume transport, input T_plume out of range: ig, jg, T_plume =', ig, jg, T_plume(i,j)
!                write(message,*) 'Plume transport, input T_plume out of range: ig, jg, T_plume =', &
!                     ig, jg, T_plume(i,j)
!                call write_log(message, GM_FATAL)
             endif
             if (S_plume(i,j) < S_plume_min .or. S_plume(i,j) > S_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
!                write(message,*) 'Plume transport, input S_plume out of range: ig, jg, S_plume =', &
!                     ig, jg, S_plume(i,j)
!                call write_log(message, GM_FATAL)
             endif
          endif
       enddo
    enddo

    ! Fill a work array with the fields to be transported
    work(:,:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx
          if (plume_mask(i,j) == 1) then
             work(i,j,1) = D_plume(i,j)
             work(i,j,2) = D_plume(i,j)*T_plume(i,j)
             work(i,j,3) = D_plume(i,j)*S_plume(i,j)
          endif
       enddo
    enddo

    ! Increment the work array based on vertical entrainment, detrainment, heat transfer and melting
    ! loop over locally owned cells
    do j = 1, ny
       do i = 1, nx
          if (plume_mask(i,j) == 1) then
             dD = entrainment(i,j) - detrainment(i,j) + bmlt_float(i,j)
             work(i,j,1) = work(i,j,1) + dD*dt
             ! Note: heat_transfer = rhoi*lhci*bmlt_float has units J/m^2/s, and
             !       rhow*cpw has units of (kg/m3)*J/(deg*kg) = J/(deg*m3),
             !       so heat_transfer/(rhow*cpw) has units of m*deg/s, as desired
             !TODO - Pass in bmlt_float and compute locally? Or compute in terms of (Tp - Tb)?
             dDT = entrainment(i,j)*T_ambient(i,j) - detrainment(i,j)*T_plume(i,j) + bmlt_float(i,j)*T_basal(i,j) &
                  - heat_transfer(i,j)/(rhow*cpw)
             work(i,j,2) = work(i,j,2) + dDT*dt
             ! Note: salt_transfer = 0 by assumption
             dDS = entrainment(i,j)*S_ambient(i,j) - detrainment(i,j)*S_plume(i,j) + bmlt_float(i,j)*S_basal(i,j)
             work(i,j,3) = work(i,j,3) + dDS*dt
          endif
       enddo
    enddo

    !WHL - debug
    ! Solve for D_plume, T_plume and S_plume; diagnostic only
!    D_temp = 0.0d0
!    T_temp = 0.0d0
!    S_temp = 0.0d0
!    do j = nhalo+1, ny-nhalo
!       do i = nhalo+1, nx-nhalo
!          if (plume_mask(i,j) == 1) then
!             D_temp(i,j) = work(i,j,1)
!             if (D_temp(i,j) > eps11) then
!                T_temp(i,j) = work(i,j,2)/D_temp(i,j)
!                S_temp(i,j) = work(i,j,3)/D_temp(i,j)
!             endif
!          endif
!       enddo
!    enddo
!    call point_diag(D_temp, 'After vertical calcs, D_plume (m)', itest, jtest, rtest, wx, wy)
!    call point_diag(T_temp, 'T_plume (degC)', itest, jtest, rtest, wx, wy)
!    call point_diag(S_temp, 'S_plume (psu)', itest, jtest, rtest, wx, wy)

    ! halo update before horizontal transport
    !TODO - May not be needed; I think all relevant fields are up to date
    do n = 1, 3
       call parallel_halo(work(:,:,n), parallel)
    enddo

    ! Set bounds for loops over locally owned cells and edges (inputs to glissade_upwind_field)
    ilo = nhalo + 1
    ihi = nx - nhalo
    jlo = nhalo +1
    jhi = ny - nhalo

    ! Use a first-order upwind scheme to transport D_plume

    call glissade_upwind_field(&
         nx,             ny,             &
         ilo, ihi,       jlo, jhi,       &
         dx,             dy,             &
         dt,             work(:,:,1),    &
         u_plume_east,   v_plume_north)

    ! Repeat for (D_plume*T_plume) and (D_plume*S_plume)

    call glissade_upwind_field(&
         nx,             ny,             &
         ilo, ihi,       jlo, jhi,       &
         dx,             dy,             &
         dt,             work(:,:,2),    &
         u_plume_east,   v_plume_north)

    call glissade_upwind_field(&
         nx,             ny,             &
         ilo, ihi,       jlo, jhi,       &
         dx,             dy,             &
         dt,             work(:,:,3),    &
         u_plume_east,   v_plume_north)

    ! Solve for D_plume, T_plume and S_plume

!    do j = nhalo+1, ny-nhalo
!       do i = nhalo+1, nx-nhalo
!          if (plume_mask(i,j) == 1) then
!             D_plume(i,j) = work(i,j,1)
!             if (D_plume(i,j) > eps11) then
!                T_plume(i,j) = work(i,j,2)/D_plume(i,j)
!                S_plume(i,j) = work(i,j,3)/D_plume(i,j)
!             else
!                T_plume(i,j) = 0.0d0
!                S_plume(i,j) = 0.0d0
!             endif
!          endif
!       enddo
!    enddo

    ! halo update after horizontal transport
    do n = 1, 3
       call parallel_halo(work(:,:,n), parallel)
    enddo

    ! Compute diffusive fluxes of heat and salt at each plume edge.
    ! Note: There are no diffusive fluxes at open boundaries (edge_mask = 2),
    !       since we assume gradT = gradS = 0 at open boundaries.

    diffT_east = 0.0d0
    diffS_east = 0.0d0
    diffT_north = 0.0d0
    diffS_north = 0.0d0

    ! loop over all edges of locally owned plumes
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          if (edge_mask_east(i,j) == 1) then
             gradT = (T_plume(i+1,j) - T_plume(i,j))/dx
             if (gradT > 0.0d0) then  ! heat flows from (i+1,j) to (i,j)
                diffT_east(i,j) = -1.0d0 * Kh * D_plume(i+1,j) * gradT * dy
             else  ! heat flows from (i,j) to (i+1,j)
                diffT_east(i,j) = -1.0d0 * Kh * D_plume(i,j) * gradT * dy
             endif
             gradS = (S_plume(i+1,j) - S_plume(i,j))/dx
             if (gradS > 0.0d0) then  ! salt flows from (i+1,j) to (i,j)
                diffS_east(i,j) = -1.0d0 * Kh * D_plume(i+1,j) * gradS * dy
             else  ! heat flows from (i,j) to (i+1,j)
                diffS_east(i,j) = -1.0d0 * Kh * D_plume(i,j) * gradS * dy
             endif
          endif
          if (edge_mask_north(i,j) == 1) then
             gradT = (T_plume(i,j+1) - T_plume(i,j))/dy
             if (gradT > 0.0d0) then  ! heat flows from (i,j+1) to (i,j)
                diffT_north(i,j) = -1.0d0 * Kh * D_plume(i,j+1) * gradT * dx
             else  ! heat flows from (i,j) to (i,j+1)
                diffT_north(i,j) = -1.0d0 * Kh * D_plume(i,j) * gradT * dx
             endif
             gradS = (S_plume(i+1,j) - S_plume(i,j))/dx
             if (gradS > 0.0d0) then  ! salt flows from (i,j+1) to (i,j)
                diffS_north(i,j) = -1.0d0 * Kh * D_plume(i,j+1) * gradS * dx
             else  ! heat flows from (i,j) to (i,j+1)
                diffS_north(i,j) = -1.0d0 * Kh * D_plume(i,j) * gradS * dx
             endif
          endif
       enddo
    enddo

    if (verbose_plume) then
       call point_diag(diffT_east *scyr/(dx*dy), 'T diffusion, east edges (m*deg/yr)',  itest, jtest, rtest, wx, wy)
       call point_diag(diffT_north*scyr/(dx*dy), 'T diffusion, north edges (m*deg/yr)', itest, jtest, rtest, wx, wy)
       call point_diag(diffS_east *scyr/(dx*dy), 'S diffusion, east edges (m*deg/yr)',  itest, jtest, rtest, wx, wy)
       call point_diag(diffS_north*scyr/(dx*dy), 'S_diffusion, north edges (m*deg/yr)', itest, jtest, rtest, wx, wy)
    endif

    ! Increment the work arrays for D*T and D*S based on the incoming and outgoing diffusive fluxes

    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             work(i,j,2) = work(i,j,2) + (diffT_east(i-1,j)  - diffT_east(i,j)  &
                                       +  diffT_north(i,j-1) - diffT_north(i,j)) * dt/(dx*dy)
             work(i,j,3) = work(i,j,3) + (diffS_east(i-1,j)  - diffS_east(i,j)  &
                                       +  diffS_north(i,j-1) - diffS_north(i,j)) * dt/(dx*dy)
          endif
       enddo
    enddo

    ! Solve for D_plume, T_plume and S_plume

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             D_plume(i,j) = work(i,j,1)
             if (D_plume(i,j) > eps11) then
                T_plume(i,j) = work(i,j,2)/D_plume(i,j)
                S_plume(i,j) = work(i,j,3)/D_plume(i,j)
             else
                T_plume(i,j) = 0.0d0
                S_plume(i,j) = 0.0d0
             endif
          endif
       enddo
    enddo

    ! final halo update
    call parallel_halo(D_plume, parallel)
    call parallel_halo(T_plume, parallel)
    call parallel_halo(S_plume, parallel)

    if (verbose_plume) then
       call point_diag(D_plume, 'New D_plume (m)', itest, jtest, rtest, wx, wy)
       call point_diag(T_plume, 'T_plume (degC)', itest, jtest, rtest, wx, wy)
       call point_diag(S_plume, 'S_plume (psu)', itest, jtest, rtest, wx, wy)
    endif

    ! Make sure all output fields are in range
    ! If D_plume is out of range, the entrainment or detrainment should work
    !  to bring it back in range.
    ! If T_plume or S_plume is out of range, there may be a bug.

    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             if (D_plume(i,j) < D_plume_min) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'Warning, thin plume: ig, jg, D_plume:', ig, jg, D_plume(i,j)
             endif
             if (D_plume(i,j) > D_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'Warning, thick plume: ig, jg, D_plume:', ig, jg, D_plume(i,j)
             endif
             if (T_plume(i,j) < T_plume_min .or. T_plume(i,j) > T_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'After transport, T_plume out of range: ig, jg, T_plume =', ig, jg, T_plume(i,j)
                write(message,*) 'After transport, T_plume out of range: ig, jg, T_plume =', &
                     ig, jg, T_plume(i,j)
                call write_log(message, GM_FATAL)
             endif
             if (S_plume(i,j) < S_plume_min .or. S_plume(i,j) > S_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(iulog,*) 'After transport, S_plume out of range: ig, jg, S_plume =', ig, jg, S_plume(i,j)
                write(message,*) 'After transport, S_plume out of range: ig, jg, S_plume =', &
                     ig, jg, S_plume(i,j)
                call write_log(message, GM_FATAL)
             endif
          endif
       enddo
    enddo

  end subroutine plume_transport

!****************************************************

  end module glissade_plume

!****************************************************
