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

    !TODO - Add these to a derived type or a constants module?
    real(dp), parameter :: &
         lambda1 = -0.0573d0,        & ! liquidus slope (deg/psu)
         lambda2 =  0.0832d0,        & ! liquidus intercept (deg C)
         lambda3 = -7.53d-8,         & ! liquidus pressure coefficient (deg/Pa)
                                       ! Tb = lambda1*Sb + lambda2 + lambda3*pb
         c_drag = 2.5d-3,            & ! ocean drag coefficient (unitless)
         u_tidal = 0.01d0,           & ! tidal velocity (m/s)
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
         D_plume_min = 2.0d0,        & ! min plume thickness (m) where the plume exists
         D_plume_max = 50.0d0,       & ! max plume thickness (m)
         plume_cavity_h0 = 75.d0,    & ! cavity thickness (m) below which entrainment is phased out
         tau_relax_entrainment = 3600. ! timescale (s) for relaxing toward D_plume_min or D_plume_max
                                       ! by imposing entrainment or detrainment

    integer, parameter :: wx = 7, wy = 10   ! block size passed to point_diag

    !WHL - debug
!!    logical, parameter :: verbose_velo = .true.
    logical, parameter :: verbose_velo = .false.

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
            plume%which_entrainment,                   &
            plume%gammaT,        plume%gammaS,         &
            plume%Kh,            plume%Ah,             &
            model%geometry%thck,                       &
            model%geometry%lsrf,                       &
            model%geometry%topg,                       &
            model%climate%eus,                         &
            plume%T_ambient,     plume%S_ambient,      &
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
            plume%which_entrainment,                   &
            plume%gammaT,        plume%gammaS,         &
            plume%Kh,            plume%Ah,             &
            model%geometry%thck,                       &
            model%geometry%lsrf,                       &
            model%geometry%topg,                       &
            model%climate%eus,                         &
            plume%T_ambient,     plume%S_ambient,      &
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
       which_entrainment,                  &
       gammaT,           gammaS,           &
       Kh,               Ah,               &
       thck,             lsrf,             &
       topg,             eus,              &
       T_ambient,        S_ambient,        &
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

    integer, intent(in) :: &
         which_entrainment      ! entrainment option

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
         Kh,                  & ! horizontal diffusivity (m^2/s)
         Ah                     ! horizontal viscosity (m^2/s)

    ! D_plume, T_plume and S_plume are prognosed variables that carry over from one step to the next.
    real(dp), dimension(nx,ny), intent(inout) :: &
         D_plume,             & ! plume thickness (m)
         T_plume,             & ! plume temperature (deg C)
         S_plume                ! plume salinity (psu)

    ! T_basal, S_basal and bmlt_float are diagnosed in this subroutine without any dependence
    !  on previous values. However, the Gaspar entrainment parameterization needs them as input.
    real(dp), dimension(nx,ny), intent(inout) :: &
         T_basal,             & ! basal ice temperature (deg C)
         S_basal,             & ! basal ice salinity (psu)
         bmlt_float             ! melt rate at base of floating ice (m/s)

    ! All other plume variables are diagnosed in this subroutine and are intent(out).
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
         divDu_plume            ! div(Du) for plume

    ! Local variables

    real(dp) :: &
         time                   ! elapsed time on the way to time_total

    integer, dimension(nx,ny) :: &
         plume_mask,          & ! = 1 for cells where plume variables are computed
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed;
                                ! = 0 at closed boundaries and = 2 at open boundaries
         edge_mask_north,     & ! = 1 on north edges where plume velocity is computed;
                                ! = 0 at closed boundaries and = 2 at open boundaries
         ice_mask,            & ! = 1 if ice is present (thck > 0)
         floating_mask,       & ! = 1 where ice is present and floating, else = 0
         ocean_mask,          & ! = 1 if topg is below sea level and ice is absent, else = 0
         land_mask              ! = 1 if topg is at or above sea level, else = 0

    real(dp), dimension(nx,ny) :: &
         pressure,            & ! ocean pressure at base of ice (N/m^2)
         lsrf_plume,          & ! elevation of plume-ambient interface (m, negative below sea level)
         rho_plume,           & ! plume density (kg/m^3)
         rho_ambient,         & ! ambient ocean density (kg/m^3)
         rho_basal,           & ! density of water at ice base (kg/m^3)
         drho_basal,          & ! density difference between plume and ice base (kg/m^3)
         H_cavity,            & ! thickness of ocean cavity beneath the plume (m)
         ufric_plume,         & ! plume friction velocity (m/s) at cell centers
         D_plume_old            ! D_plume from previous time step

    ! plume speed on cell edges
    ! Note: u plume_east and v_plume_north (the C grid velocity components) are primary
    !        and are input/output varaibles
    !       u_plume_north and v_plume_east (the D grid components) are computed as part of
    !        the velocity solution but are not used again.
    !TODO - Decide what to output. Maybe just the transport velocities?
    !       Maybe don't output the overall u_plume and v_plume, since these are just averages?

    real(dp), dimension(nx,ny) :: &
         v_plume_east,          & ! v_plume on east edges
         u_plume_north            ! u_plume on north edges

    real(dp), dimension(nx,ny) :: &
         u_transport_east,      & ! transport velocity u on east edges (m/s)
         v_transport_north        ! transport velocity v on north edges (m/s)

    real(dp) ::  &
!         my_max_dt,          & ! CFL-limited time step for a given cell (s)
         L2_norm,             & ! L2 norm of residual vector from continuity equation
         L2_previous            ! L2 norm from the previous convergence check

    integer :: i, j, ig, jg

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

    !WHL - debug
    real(dp) :: h1, h2, eT, mT

    ! parameters determining convergence of iterations
    !TODO - determine L2_target
    integer, parameter :: &
         L2_target = 0.0d0           ! convergence target for dD/dt
!!         n_check_convergence = 1,    & ! interval between convergence checks for D_plume
!!         maxiter_Dplume = 999999     ! max number of iterations of outer plume-thickness loop
                                       ! terminates when plume thickness reaches virtual steady state

!    logical, parameter :: average_transport_velocity = .false.
    logical, parameter :: average_transport_velocity = .true.

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

    ! Compute the density of the ambient ocean
    rho_ambient = eos_rho_ref * (1.d0 - eos_alpha * (T_ambient - eos_Tref)  &
                                      + eos_beta  * (S_ambient - eos_Sref) )

    ! Compute the pressure at the lower ice surface.
    pressure = -rhoo*grav*lsrf

    ! Compute the cavity thickness
    H_cavity = max(lsrf - topg, 0.0d0)

    !----------------------------------------------------------------------------
    ! Initialize D_plume, T_plume, S_plume, T_basal and S_basal as needed.
    ! On the first call, the input values are zero and these fields must be initialized everywhere.
    ! On subsequent calls, these fields are initialized only if the input values are zero.
    !
    ! Initialize T_plume = T_ambient and S_plume = S_ambient - 0.1 psu, following Jesse et al. (2025).
    ! Setting S_plume < S_ambient ensures stable stratification.
    ! Setting both T_plume and S_plume to ambient values would give zero PGF, velocities, and drho_plume.
    !
    ! Note: T_basal and S_basal are diagnosed in plume_melt_rate without regard to their initial values.
    !       However, T_basal and S_basal from the previous step might be needed to compute drho_basal
    !        for entrainment.
    !----------------------------------------------------------------------------

    ! loop over all cells (fields on the rhs are up to date in halos)
    do j = 1, ny
       do i = 1, nx
          if (plume_mask(i,j) == 1) then
             if (D_plume(i,j) == 0.0d0) D_plume(i,j) = min(D_plume0, H_cavity(i,j))
             if (T_plume(i,j) == 0.0d0) T_plume(i,j) = T_ambient(i,j)
             if (S_plume(i,j) == 0.0d0) S_plume(i,j) = S_ambient(i,j) - 0.1d0
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
    ! If both adjacent cells have plume_mask = 1, then edge_mask = 1.
    ! At open boundaries (one adjacent cell is open ocean), edge_mask = 2.
    ! At closed boundaries (one adjecent cell is grounded), edge_mask = 3.
    ! If neither cell is a plume cell, edge_mask = 0.
    ! For edge_mask = 1 or 2, we will compute both velocity components.
    ! For edge_mask = 3, we will compute the parallel component only,
    !  setting the perpendicular component to zero.
    !----------------------------------------------------------------------------

    edge_mask_east = 0
    edge_mask_north = 0

    ! loop over all edges of locally owned cells
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! east edges
          if (plume_mask(i,j) == 1 .and. plume_mask(i+1,j) == 1) then
             edge_mask_east(i,j) = 1
          elseif (plume_mask(i,j) == 1) then
             if (lsrf(i+1,j) > topg(i+1,j)) then   ! open boundary
                edge_mask_east(i,j) = 2
             elseif (lsrf(i+1,j) == topg(i+1,j)) then  ! closed boundary
                edge_mask_east(i,j) = 3
             endif
          elseif (plume_mask(i+1,j) == 1) then
             if (lsrf(i,j) > topg(i,j)) then   ! open boundary
                edge_mask_east(i,j) = 2
             elseif (lsrf(i,j) == topg(i,j)) then  ! closed boundary
                edge_mask_east(i,j) = 3
             endif
          endif

          ! north edges
          if (plume_mask(i,j) == 1 .and. plume_mask(i,j+1) == 1) then
             edge_mask_north(i,j) = 1
          elseif (plume_mask(i,j) == 1) then
             if (lsrf(i,j+1) > topg(i,j+1)) then   ! open boundary
                edge_mask_north(i,j) = 2
             elseif (lsrf(i,j+1) == topg(i,j+1)) then  ! closed boundary
                edge_mask_north(i,j) = 3
             endif
          elseif (plume_mask(i,j+1) == 1) then
             if (lsrf(i,j) > topg(i,j)) then   ! open boundary
                edge_mask_north(i,j) = 2
             elseif (lsrf(i,j) == topg(i,j)) then  ! closed boundary
                edge_mask_north(i,j) = 3
             endif
          endif

       enddo   ! i
    enddo   ! j

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
    ufric_plume = 0.0d0
    drho_basal = 0.0d0
    drho_plume = 0.0d0
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
       call point_diag(u_plume_east, 'u_plume_east (m/s)', itest, jtest, rtest, wx, wy)
       call point_diag(v_plume_east, 'v_plume_east (m/s)', itest, jtest, rtest, wx, wy)
       call point_diag(u_plume_north, 'u_plume_north (m/s)', itest, jtest, rtest, wx, wy)
       call point_diag(v_plume_north, 'v_plume_north (m/s)', itest, jtest, rtest, wx, wy)
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

       if (verbose_plume) then
          call point_diag(lsrf_plume, 'lsrf_plume (m)', itest, jtest, rtest, wx, wy)
          call point_diag(lsrf_plume - topg, 'lsrf_plume - topg (m)', itest, jtest, rtest, wx, wy)
          call point_diag(rho_ambient, 'rho_ambient (kg/m3)', itest, jtest, rtest, wx, wy)
          call point_diag(rho_plume, 'rho_plume (kg/m3)', itest, jtest, rtest, wx, wy)
          call point_diag(drho_plume, 'drho_plume (kg/m3)', itest, jtest, rtest, wx, wy)
          call point_diag(rho_basal, 'rho_basal (kg/m3)', itest, jtest, rtest, wx, wy)
          call point_diag(drho_basal, 'drho_basal (kg/m3)', itest, jtest, rtest, wx, wy)
          if (this_rank == rtest) write(iulog,*) 'Compute plume velocity'
       endif

       !----------------------------------------------------------------------------
       ! Compute u_plume and v_plume at each edge
       ! Note: u_plume_east and v_plume_north are perpendicular to edges,
       !        whereas v_plume_north and u_plume_east are parallel to edges.
       !       Computing both u and v at each edge leads to a more graceful treatment
       !        of the Coriolis terms than computing the perpendicular components alone.
       ! Note: The output velocities are correct in halos.
       !----------------------------------------------------------------------------

       call compute_plume_velocity(&
            nx,           ny,       &
            dx,           dy,       &
            itest, jtest, rtest,    &
            parallel,               &
            Ah,                     &
            plume_mask,             &
            edge_mask_east,         &
            edge_mask_north,        &
            D_plume,                &
            H_cavity,               &
            drho_plume,             &
            lsrf_plume,             &
            u_plume_east,           &
            v_plume_east,           &
            u_plume_north,          &
            v_plume_north)

       !--------------------------------------------------------------------
       ! Compute the plume speed and friction velocity at cell centers
       !--------------------------------------------------------------------

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (plume_mask(i,j) == 1) then
                u_plume(i,j) = 0.25d0*(u_plume_east(i-1,j) + u_plume_east(i,j) &
                                     + u_plume_north(i,j-1) + u_plume_north(i,j))
                v_plume(i,j) = 0.25d0*(v_plume_east(i-1,j) + v_plume_east(i,j) &
                                     + v_plume_north(i,j-1) + v_plume_north(i,j))
             endif
          enddo
       enddo

       call parallel_halo(u_plume, parallel)
       call parallel_halo(v_plume, parallel)

       do j = 1, ny
          do i = 1, nx
             plume_speed(i,j) = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)
             ufric_plume(i,j) = sqrt(c_drag) * plume_speed(i,j)
          enddo
       enddo

       !--------------------------------------------------------------------
       ! Compute the entrainment rate by one of two methods:
       ! - Jenkins(1991): Entrainment is a function of slope and plume speed.
       ! - Gaspar (1998): Entrainment is a function of sources and sinks of
       !                  turbulent kinetic energy.
       ! Note: At this point in the code, all relevant quantities for entrainment
       !       are up to date in halos.
       !--------------------------------------------------------------------

       if (this_rank == rtest) then
          write(iulog,*) 'call plume_entrainment, which_entrainment =', which_entrainment
       endif

       call plume_entrainment(&
            nx,         ny,      &
            dx,         dy,      &
            itest, jtest, rtest, &
            parallel,            &
            which_entrainment,   &
            plume_mask,          &
            edge_mask_east,      &
            edge_mask_north,     &
            lsrf_plume,          &
            plume_speed,         &
            bmlt_float,          &
            drho_plume,          &
            drho_basal,          &
            H_cavity,            &
            D_plume,             &
            dt_plume,            &
            entrainment,         &
            detrainment)

       if (verbose_plume) then
          if (this_rank == rtest) write(iulog,*) 'Compute melt rate'
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
            ufric_plume,         &
            D_plume,             &
            T_plume,             &
            S_plume,             &
            T_basal,             &
            S_basal,             &
            bmlt_float)

       if (verbose_plume) then
          if (this_rank == rtest) write(iulog,*) 'After melt calculation:'
          call point_diag(T_basal, 'T_basal (deg C)', itest, jtest, rtest, wx, wy)
          call point_diag(S_basal, 'S_basal (psu)', itest, jtest, rtest, wx, wy)
          call point_diag(bmlt_float*scyr, 'bmlt_float (m/yr)', itest, jtest, rtest, wx, wy)
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

       !WHL - debug - write expressions for heat loss.
       if (this_rank == rtest) then
          i = itest
          j = jtest
          h1 = gammaT * ufric_plume(i,j) * (T_plume(i,j) - T_basal(i,j))
          h2 = rhoi*lhci*bmlt_float(i,j)/(rhoo*cpw)
          eT = entrainment(i,j) * (T_ambient(i,j) - T_plume(i,j))
          mT = bmlt_float(i,j) * (T_plume(i,j) - T_basal(i,j))
          write(iulog,*) 'Transfer of heat to ice (deg*m/yr):', h1*scyr, h2*scyr
          write(iulog,*) 'E, Ta - Tp, heat gain from ent (deg*m/yr):', entrainment(i,j)*scyr, T_ambient(i,j) - T_plume(i,j), et*scyr
          write(iulog,*) 'm, Tp - Tb, heat loss from bmlt (deg*m/yr):', bmlt_float(i,j)*scyr, T_plume(i,j) - T_basal(i,j), mT*scyr

       endif

       !--------------------------------------------------------------------
       ! Optionally, compute a weighted average of the local velocity components
       ! on each edge with velocity components from neighboring edges,
       ! and use this weighted velocity for transport.
       !--------------------------------------------------------------------

       !WHL - debug - Try averaging the plume velocities in a different way
       if (average_transport_velocity) then

          call compute_transport_velocity(&
               nx,               ny,               &
               itest,   jtest,   rtest,            &
               parallel,                           &
               edge_mask_east,   edge_mask_north,  &
               u_plume_east,     v_plume_east,     &
               u_plume_north,    v_plume_north,    &
               u_transport_east, v_transport_north)

       else

          u_transport_east = u_plume_east
          v_transport_north = v_plume_north

       endif

       call point_diag(u_transport_east, 'u_transport_east (m/s)', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(v_transport_north, 'v_transport_north (m/s)', itest, jtest, rtest, wx, wy, '(f10.5)')

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
            Kh,                   &
            plume_mask,           &
            edge_mask_east,       &
            edge_mask_north,      &
!            u_plume_east,         &
!            v_plume_north,        &
            u_transport_east,         &
            v_transport_north,        &
            entrainment,          &
            detrainment,          &
            bmlt_float,           &
            T_ambient,            &
            S_ambient,            &
            T_basal,              &
            S_basal,              &
            D_plume,              &
            T_plume,              &
            S_plume,              &
            divDu_plume)

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

  subroutine compute_plume_velocity(&
       nx,           ny,       &
       dx,           dy,       &
       itest, jtest, rtest,    &
       parallel,               &
       Ah,                     &
       plume_mask,             &
       edge_mask_east,         &
       edge_mask_north,        &
       D_plume,                &
       H_cavity,               &
       drho_plume,             &
       lsrf_plume,             &
       u_plume_east,           &
       v_plume_east,           &
       u_plume_north,          &
       v_plume_north)

    !-------------------------------------------------------------------------
    ! Compute the u and v velocity components at the edges of each plume cell
    !-------------------------------------------------------------------------

    use glissade_grid_operators, only: glissade_stagger

    ! input/output arguments

    integer, intent(in) ::  &
         nx,     ny               ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy               ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest      ! diagnostic indices

    type(parallel_type), intent(in) :: &
         parallel                 ! info for parallel communication

    real(dp), intent(in) :: &
         Ah                       ! horizontal viscosity (m^2/s)

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask,            & ! = 1 for cells where plume variables are computed
         edge_mask_east,        & ! = 1 on east edges where plume velocity is computed
         edge_mask_north          ! = 1 on north edges where plume velocity is computed

    real(dp), dimension(nx,ny), intent(in) ::  &
         D_plume,               & ! plume thickness (m)
         H_cavity,              & ! thickness of ocean cavity beneath the plume (m)
         drho_plume,            & ! density difference between plume and ambient ocean (kg/m^3)
         lsrf_plume               ! elevation of plume-ambient interface (m, negative below sea level)

    ! Note: These are intent(inout) since we use the input values as the initial guess
    real(dp), dimension(nx,ny), intent(inout) ::  &
         u_plume_east,          & ! u_plume on east edges
         v_plume_east,          & ! v_plume on east edges
         u_plume_north,         & ! u_plume on north edges
         v_plume_north            ! v_plume on north edges

    ! local variables

    real(dp), dimension(nx,ny) ::  &
         dlsrf_plume_dx_east,   & ! horizontal gradient of lsrf_plume on east edges
         dlsrf_plume_dy_east,   & !
         dlsrf_plume_dx_north,  & ! horizontal gradient of lsrf_plume on north edges
         dlsrf_plume_dy_north,  & !
         ddrho_plume_dx_east,   & ! horizontal gradient of drho_plume on east edges
         ddrho_plume_dy_east,   & !
         ddrho_plume_dx_north,  & ! horizontal gradient of drho_plume on north edges
         ddrho_plume_dy_north

    real(dp), dimension(nx,ny) :: &
         drhox,                 & ! density gradient term of pgf_x
         drhoy,                 & ! density gradient term of pgf_y
         dsrfx,                 & ! surface gradient term of pgf_x
         dsrfy,                 & ! surface gradient term of pgf_y
         pgf_x_east,            & ! x component of pressure gradient force on east edges (m^2/s^2)
         pgf_y_east,            & ! y component of pressure gradient force on east edges (m^2/s^2)
         pgf_x_north,           & ! x component of pressure gradient force on north edges (m^2/s^2)
         pgf_y_north              ! y component of pressure gradient force on north edges (m^2/s^2)

    real(dp), dimension(nx,ny) :: &
         wall_factor_east,      & ! factor to reduce northward or southward flow on east edges (unitless)
         wall_factor_north        ! factor to reduce eastward or westward flow on north edges (unitless)

    real(dp), dimension(nx,ny) :: &
         uctr, vctr,            & ! u_plume and v_plume averaged to cell centers
         stagu, stagv,          & ! u_plume and v_plume averaged to cell corners
         Avisc_east_lhs,        & ! LHS lateral viscosity coefficents on east edges (m/s)
         Avisc_east_rhsu,       & ! RHS lateral viscosity coefficents on east edges, u equation (m^2/s^2)
         Avisc_east_rhsv,       & ! RHS lateral viscosity coefficents on east edges, v equation (m^2/s^2)
         Avisc_north_lhs,       & ! LHS lateral viscosity coefficents on north edges (m/s)
         Avisc_north_rhsu,      & ! RHS lateral viscosity coefficents on north edges, u equation (m^2/s^2)
         Avisc_north_rhsv         ! RHS lateral viscosity coefficents on north edges, v equation (m^2/s^2)

    real(dp), dimension(nx,ny) :: &
         D_plume_east,          & ! D_plume averaged to east edge
         D_plume_north            ! D_plume averaged to north edge

    real(dp), dimension(nx-1,ny-1) :: &
         stagD_plume              ! D_plume averaged to cell corners

    integer :: i, j, ig, jg

    integer :: iter_velo          ! iteration counter

    real(dp) :: &
         grav_reduced,          & ! reduced gravity (m/s^2)
         xterm, xterm1, xterm2, & ! components of viscosity terms
         yterm, yterm1, yterm2

    character(len=100) :: message

    logical, dimension(nx,ny) ::  &
         converged_velo_east,   & ! true when velocity has converged at an east edge, else false
         converged_velo_north     ! true when velocity has converged at a north edge, else false

    integer :: &
         count_east, count_north  ! number of cells not converged on each face

    integer, parameter ::  &
         maxiter_velo = 50        ! max number of iterations of velocity loop

    !----------------------------------------------------------------------------
    ! Compute horizontal gradients of lsrf_plume and drho_plume at each edge.
    ! Use Neumann BC for drho_plume; this gives zero gradients at all boundaries.
    ! Use the actual value of lsrf_plume at boundaries; this potentially gives
    !  a large gradient at open boundaries.
    ! Note: The subroutine includes halo updates.
    !----------------------------------------------------------------------------

    call compute_edge_gradients(&
         nx,                   ny,                    &
         dx,                   dy,                    &
         parallel,                                    &
         plume_mask,                                  &
         edge_mask_east,       edge_mask_north,       &
         lsrf_plume,                                  &
         dlsrf_plume_dx_east,  dlsrf_plume_dy_east,   &
         dlsrf_plume_dx_north, dlsrf_plume_dy_north)

    call compute_edge_gradients(&
         nx,                   ny,                    &
         dx,                   dy,                    &
         parallel,                                    &
         plume_mask,                                  &
         edge_mask_east,       edge_mask_north,       &
         drho_plume,                                  &
         ddrho_plume_dx_east,  ddrho_plume_dy_east,   &
         ddrho_plume_dx_north, ddrho_plume_dy_north,  &
         neumann_bc = .true.)

    !-------------------------------------------------------------------
    ! Compute the pressure gradient force on each edge, following Lambert et al. (2023):
    ! (1) pgf_x = -(g*D^2)/(2*rhoo) d/dx(drho_plume) + g'*D d/dx(zb - D)
    ! (2) pgf_y = -(g*D^2)/(2*rhoo) d/dy(drho_plume) + g'*D d/dy(zb - D)
    !
    ! where zb = lower plume surface
    !       g' = g * drho_plume/rhoo = reduced gravity
    !-------------------------------------------------------------------

    D_plume_east = 0.0d0
    D_plume_north = 0.0d0
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
          if (edge_mask_east(i,j) == 1) then   ! plume cell on each side
             D_plume_east(i,j) = (D_plume(i,j) + D_plume(i+1,j)) / 2.0d0
             ! terms proportional to gradients of drho_plume
             drhox(i,j) = -0.5d0*(grav/rhoo) * D_plume_east(i,j)**2 * ddrho_plume_dx_east(i,j)
             drhoy(i,j) = -0.5d0*(grav/rhoo) * D_plume_east(i,j)**2 * ddrho_plume_dy_east(i,j)
             grav_reduced = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i+1,j)) / 2.0d0
             dsrfx(i,j) = grav_reduced * D_plume_east(i,j) * dlsrf_plume_dx_east(i,j)
             dsrfy(i,j) = grav_reduced * D_plume_east(i,j) * dlsrf_plume_dy_east(i,j)
          elseif (edge_mask_east(i,j) > 1) then  ! boundary cell; plume cell on only one side
             if (plume_mask(i,j) == 1) then  ! plume cell to the west
                D_plume_east(i,j) = D_plume(i,j)
                grav_reduced = (grav/rhoo) * drho_plume(i,j)/2.0d0
             elseif (plume_mask(i+1,j) == 1) then  ! plume cell to the east
                D_plume_east(i,j) = D_plume(i+1,j)
                grav_reduced = (grav/rhoo) * drho_plume(i+1,j)/2.0d0
             endif
             ! PGF force at boundaries includes the surface gradient term but not the density term
             dsrfx(i,j) = grav_reduced * D_plume_east(i,j) * dlsrf_plume_dx_east(i,j)
             dsrfy(i,j) = grav_reduced * D_plume_east(i,j) * dlsrf_plume_dy_east(i,j)
          endif   ! edge_mask_east
          pgf_x_east(i,j) = drhox(i,j) + dsrfx(i,j)
          pgf_y_east(i,j) = drhoy(i,j) + dsrfy(i,j)
       enddo
    enddo

    if (verbose_plume) then
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'PGF components on east edges:'
       endif
       call point_diag(1.d3*drhox, '10^3*density gradient x term', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*drhoy, '10^3*density gradient y term', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*dsrfx, '10^3*surface gradient x term', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*dsrfy, '10^3*surface gradient y term', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*pgf_x_east, '10^3*pgf_x_east (m2/s2)', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*pgf_y_east, '10^3*pgf_y_east (m2/s2)', itest, jtest, rtest, wx, wy, '(f10.5)')
    endif

    ! Reset the PGF components
    drhox = 0.0d0
    drhoy = 0.0d0
    dsrfx = 0.0d0
    dsrfy = 0.0d0

    ! PGF on north edges
    ! Loop over all edges of locally owned cells (includes south halo cells)
    do j = nhalo, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (edge_mask_north(i,j) == 1) then
             D_plume_north(i,j) = (D_plume(i,j) + D_plume(i,j+1)) / 2.0d0
             ! terms proportional to gradients of drho_plume
             drhox(i,j) = -0.5d0*(grav/rhoo) * D_plume_north(i,j)**2 * ddrho_plume_dx_north(i,j)
             drhoy(i,j) = -0.5d0*(grav/rhoo) * D_plume_north(i,j)**2 * ddrho_plume_dy_north(i,j)
             grav_reduced = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i,j+1)) / 2.0d0
             dsrfx(i,j) = grav_reduced * D_plume_north(i,j) * dlsrf_plume_dx_north(i,j)
             dsrfy(i,j) = grav_reduced * D_plume_north(i,j) * dlsrf_plume_dy_north(i,j)
          elseif (edge_mask_north(i,j) > 1) then  ! boundary cell; plume cell on only one side
             if (plume_mask(i,j) == 1) then  ! plume cell to the south
                D_plume_north(i,j) = D_plume(i,j)
                grav_reduced = (grav/rhoo) * drho_plume(i,j)/2.0d0
             elseif (plume_mask(i,j+1) == 1) then  ! plume cell to the north
                D_plume_north(i,j) = D_plume(i,j+1)
                grav_reduced = (grav/rhoo) * drho_plume(i,j+1)/2.0d0
             endif
             ! PGF force at boundaries includes the surface gradient term but not the density term
             dsrfx(i,j) = grav_reduced * D_plume_north(i,j) * dlsrf_plume_dx_north(i,j)
             dsrfy(i,j) = grav_reduced * D_plume_north(i,j) * dlsrf_plume_dy_north(i,j)
          endif   ! edge_mask_north
          pgf_x_north(i,j) = drhox(i,j) + dsrfx(i,j)
          pgf_y_north(i,j) = drhoy(i,j) + dsrfy(i,j)
       enddo  ! i
    enddo  ! j

    if (verbose_plume) then
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'PGF components on north edges:'
       endif
       call point_diag(1.d3*drhox, '10^3*density gradient x term', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*drhoy, '10^3*density gradient y term', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*dsrfx, '10^3*surface gradient x term', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*dsrfy, '10^3*surface gradient y term', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*pgf_x_north, '10^3*pgf_x_north (m2/s2)', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(1.d3*pgf_y_north, '10^3*pgf_y_north (m2/s2)', itest, jtest, rtest, wx, wy, '(f10.5)')
    endif

    !--------------------------------------------------------------------
    ! Compute a field that will reduce or eliminate the Coriolis term on edges
    !  adjacent to closed boundaries.
    ! This allows a jet of strong flow along the boundary, instead of trapping
    !  and thickening the plume in cells next to the boundary.
    !--------------------------------------------------------------------

    !TODO - Count it as closed if edge_mask_north or edge_mask_east = 0 also
    wall_factor_east = 1.0d0
    wall_factor_north = 1.0d0

    do j = nhalo+1, ny-nhalo
       do i = nhalo, nx-nhalo
          if (plume_mask(i,j) == 1) then
             ! check for a closed boundary to the north
             ! if present, compute a factor that will reduce the Coriolos force on east edges
             if (edge_mask_north(i,j) == 3 .and. edge_mask_north(i+1,j) == 3) then
                wall_factor_east(i,j) = 0.0d0
             elseif (edge_mask_north(i,j) == 3 .or. edge_mask_north(i+1,j) == 3) then
                wall_factor_east(i,j) = 0.5d0
             endif
             ! check for a closed boundary to the south
             ! if present, compute a factor that will reduce the Coriolis force on east edges
             if (edge_mask_north(i,j-1) == 3 .and. edge_mask_north(i+1,j-1) == 3) then
                wall_factor_east(i,j) = 0.0d0
             elseif (edge_mask_north(i,j-1) == 3 .or. edge_mask_north(i+1,j-1) == 3) then
                wall_factor_east(i,j) = 0.5d0
             endif
             ! check for a closed boundary to the west
             ! if present, compute a factor that will reducethe Coriolis force on north edges
             if (edge_mask_east(i-1,j-1) == 3 .and. edge_mask_east(i-1,j) == 3) then
                wall_factor_north(i,j) = 0.0d0
             elseif (edge_mask_east(i-1,j) == 3 .or. edge_mask_east(i-1,j) == 3) then
                wall_factor_north(i,j) = 0.5d0
             endif
             ! check for a closed boundary to the east
             ! if present, compute a factor that will reducethe Coriolis force on north edges
             if (edge_mask_east(i,j-1) == 3 .and. edge_mask_east(i,j) == 3) then
                wall_factor_north(i,j) = 0.0d0
             elseif (edge_mask_east(i,j-1) == 3 .or. edge_mask_east(i,j) == 3) then
                wall_factor_north(i,j) = 0.5d0
             endif
          endif   ! plume_mask
       enddo   ! i
    enddo   ! j

    call parallel_halo(wall_factor_east, parallel)
    call parallel_halo(wall_factor_north, parallel)

    if (verbose_plume) then
       call point_diag(wall_factor_east, 'wall_factor_east', itest, jtest, rtest, wx, wy, '(f10.1)')
       call point_diag(wall_factor_north, 'wall_factor_north', itest, jtest, rtest, wx, wy, '(f10.1)')
    endif

    !--------------------------------------------------------------------
    ! Compute the LHS lateral viscosity terms on each edge, excluding boundaries.
    ! These coefficients multiply the current velocity components on the edge
    !  and do not change during the velocity iteration below.
    ! For details, see the comments below in subroutine plume_velocity.
    !--------------------------------------------------------------------

    ! Interpolate D_plume to the staggered grid; these values appear in some viscosity terms.
    ! Include D_plume in the average only for cells with plume_mask = 1.

    call glissade_stagger(&
         nx,           ny,            &
         D_plume,      stagD_plume,   &
         plume_mask,   stagger_margin_in = 1)

    Avisc_east_lhs = 0.0d0
    Avisc_north_lhs = 0.0d0

    if (Ah > 0.0d0) then
       ! Loop over all edges of locally owned cells (including south and west halo cells)
       do j = nhalo, ny-nhalo
          do i = nhalo, nx-nhalo
             if (edge_mask_east(i,j) == 1) then   ! plume cell on each side
                xterm = (2.0d0*Ah/dx**2) * (D_plume(i,j) + D_plume(i+1,j))
                yterm = (2.0d0*Ah/dy**2) * (stagD_plume(i,j-1) + stagD_plume(i,j))
                Avisc_east_lhs(i,j) = xterm + yterm
             endif
             if (edge_mask_north(i,j) == 1) then
                xterm = (2.0d0*Ah/dx**2) * (stagD_plume(i-1,j) + stagD_plume(i,j))
                yterm = (2.0d0*Ah/dy**2) * (D_plume(i,j) + D_plume(i,j+1))
                Avisc_north_lhs(i,j) = xterm + yterm
             endif
          enddo
       enddo
       call parallel_halo(Avisc_east_lhs, parallel)
       call parallel_halo(Avisc_north_lhs, parallel)
    endif

    converged_velo_east = .false.
    converged_velo_north = .false.

    !--------------------------------------------------------------------
    ! Iterate as needed to compute a converged velocity at each edge
    !--------------------------------------------------------------------

    do iter_velo = 1, maxiter_velo

       !--------------------------------------------------------------------
       ! Update the RHS viscosity terms using the four neighboring velocities
       !  from the previous iteration.
       ! For details, see the comments below in subroutine plume_velocity.
       !--------------------------------------------------------------------

       Avisc_east_rhsu = 0.0d0
       Avisc_east_rhsv = 0.0d0
       Avisc_north_rhsu = 0.0d0
       Avisc_north_rhsv = 0.0d0

       if (Ah > 0.0d0) then

          !--------------------------------------------------------------------
          ! Compute RHS terms on east edges, based on u_plume_north and v_plume_north.
          !--------------------------------------------------------------------

          ! Average uvel_north and vvel_north to cell centers and corners
          uctr = 0.0d0
          vctr = 0.0d0
          stagu = 0.0d0
          stagv = 0.0d0

          do j = nhalo+1, ny-nhalo
             do i = nhalo+1, nx-nhalo
                if (plume_mask(i,j) == 1) then
                   uctr(i,j) = 0.5d0*(u_plume_north(i,j-1) + u_plume_north(i,j))
                   vctr(i,j) = 0.5d0*(v_plume_north(i,j-1) + v_plume_north(i,j))
                   stagu(i,j) = 0.5d0*(u_plume_north(i,j) + u_plume_north(i+1,j))
                   stagv(i,j) = 0.5d0*(v_plume_north(i,j) + v_plume_north(i+1,j))
                endif
             enddo
          enddo

          call parallel_halo(uctr, parallel)
          call parallel_halo(vctr, parallel)
          call parallel_halo(stagu, parallel)
          call parallel_halo(stagv, parallel)

          ! Loop over all edges of locally owned cells (including south and west halo cells).
          do j = nhalo, ny-nhalo
             do i = nhalo, nx-nhalo
                if (edge_mask_east(i,j) == 1) then   ! exclude edges on boundaries
                   xterm = (2.0d0*Ah/dx**2) * (D_plume(i,j)*uctr(i,j) + D_plume(i+1,j)*uctr(i+1,j))
                   yterm = (2.0d0*Ah/dy**2) * (stagD_plume(i,j-1)*stagu(i,j-1) + stagD_plume(i,j)*stagu(i,j))
                   Avisc_east_rhsu(i,j) = xterm + yterm
                   xterm = (2.0d0*Ah/dx**2) * (D_plume(i,j)*vctr(i,j) + D_plume(i+1,j)*vctr(i+1,j))
                   yterm = (2.0d0*Ah/dy**2) * (stagD_plume(i,j-1)*stagv(i,j-1) + stagD_plume(i,j)*stagv(i,j))
                   Avisc_east_rhsv(i,j) = xterm + yterm

                   !WHL - debug
                   if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                      write(iulog,*) ' '
                      write(iulog,*) 'updating viscosity terms, east edge, iter =', iter_velo
                      write(iulog,*) '-Avisc_east_lhs:', -Avisc_east_lhs(i,j)
                      write(iulog,*) '   -Avisc_east_lhs*u:', -Avisc_east_lhs(i,j)*u_plume_east(i,j)
                      write(iulog,*) '   -Avisc_east_lhs*v:', -Avisc_east_lhs(i,j)*v_plume_east(i,j)
                      write(iulog,*) 'Avisc_east_rhsu:', Avisc_east_rhsu(i,j)
                      write(iulog,*) '   W term:', (2.0d0*Ah/dx**2)*D_plume(i,j)*uctr(i,j)
                      write(iulog,*) '   E term:', (2.0d0*Ah/dx**2)*D_plume(i+1,j)*uctr(i+1,j)
                      write(iulog,*) '   S term:', (2.0d0*Ah/dy**2)*stagD_plume(i,j-1)*stagu(i,j-1)
                      write(iulog,*) '   N term:', (2.0d0*Ah/dy**2)*stagD_plume(i,j)*stagu(i,j)
                      write(iulog,*) 'Avisc_east_rhsv:', Avisc_east_rhsv(i,j)
                      write(iulog,*) '   W term:', (2.0d0*Ah/dx**2)*D_plume(i,j)*vctr(i,j)
                      write(iulog,*) '   E term:', (2.0d0*Ah/dx**2)*D_plume(i+1,j)*vctr(i+1,j)
                      write(iulog,*) '   S term:', (2.0d0*Ah/dy**2)*stagD_plume(i,j-1)*stagv(i,j-1)
                      write(iulog,*) '   N term:', (2.0d0*Ah/dy**2)*stagD_plume(i,j)*stagv(i,j)
                   endif

                endif
             enddo
          enddo

          call parallel_halo(Avisc_east_rhsu, parallel)
          call parallel_halo(Avisc_east_rhsv, parallel)

          !--------------------------------------------------------------------
          ! Compute RHS terms on north edges, based on u_plume_east and v_plume_east.
          !--------------------------------------------------------------------

          ! Average uvel_east and vvel_east to cell centers and corners
          uctr = 0.0d0
          vctr = 0.0d0
          stagu = 0.0d0
          stagv = 0.0d0

          do j = nhalo+1, ny-nhalo
             do i = nhalo+1, nx-nhalo
                if (plume_mask(i,j) == 1) then
                   uctr(i,j) = 0.5d0*(u_plume_east(i-1,j) + u_plume_east(i,j))
                   vctr(i,j) = 0.5d0*(v_plume_east(i-1,j) + v_plume_east(i,j))
                   stagu(i,j) = 0.5d0*(u_plume_east(i,j) + u_plume_east(i,j+1))
                   stagv(i,j) = 0.5d0*(v_plume_east(i,j) + v_plume_east(i,j+1))
                endif
             enddo
          enddo

          call parallel_halo(uctr, parallel)
          call parallel_halo(vctr, parallel)
          call parallel_halo(stagu, parallel)
          call parallel_halo(stagv, parallel)

          ! Loop over all edges of locally owned cells (including south and west halo cells).
          do j = nhalo, ny-nhalo
             do i = nhalo, nx-nhalo
                if (edge_mask_north(i,j) == 1) then   ! exclude edges on boundaries
                   xterm = (2.0d0*Ah/dx**2) * (stagD_plume(i-1,j)*stagu(i-1,j) + stagD_plume(i,j)*stagu(i,j))
                   yterm = (2.0d0*Ah/dy**2) * (D_plume(i,j)*uctr(i,j) + D_plume(i,j+1)*uctr(i,j+1))
                   Avisc_north_rhsu(i,j) = xterm + yterm
                   xterm = (2.0d0*Ah/dx**2) * (stagD_plume(i-1,j)*stagv(i-1,j) + stagD_plume(i,j)*stagv(i,j))
                   yterm = (2.0d0*Ah/dy**2) * (D_plume(i,j)*vctr(i,j) + D_plume(i,j+1)*vctr(i,j+1))
                   Avisc_north_rhsv(i,j) = xterm + yterm

                   !WHL - debug
                   if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                      write(iulog,*) ' '
                      write(iulog,*) 'updating viscosity terms, north edge, iter =', iter_velo
                      write(iulog,*) '-Avisc_north_lhs:', -Avisc_north_lhs(i,j)
                      write(iulog,*) '   -Avisc_north_lhs*u:', -Avisc_north_lhs(i,j)*u_plume_north(i,j)
                      write(iulog,*) '   -Avisc_north_lhs*v:', -Avisc_north_lhs(i,j)*v_plume_north(i,j)
                      write(iulog,*) 'Avisc_north_rhsu:', Avisc_north_rhsu(i,j)
                      write(iulog,*) '   W term:', (2.0d0*Ah/dx**2)*stagD_plume(i-1,j)*stagu(i-1,j)
                      write(iulog,*) '   E term:', (2.0d0*Ah/dx**2)*stagD_plume(i,j)*stagu(i,j)
                      write(iulog,*) '   S term:', (2.0d0*Ah/dy**2)*D_plume(i,j)*uctr(i,j)
                      write(iulog,*) '   N term:', (2.0d0*Ah/dy**2)*D_plume(i,j+1)*uctr(i,j+1)
                      write(iulog,*) 'Avisc_north_rhsv:', Avisc_north_rhsv(i,j)
                      write(iulog,*) '   W term:', (2.0d0*Ah/dx**2)*stagD_plume(i-1,j)*stagv(i-1,j)
                      write(iulog,*) '   E term:', (2.0d0*Ah/dx**2)*stagD_plume(i,j)*stagv(i,j)
                      write(iulog,*) '   S term:', (2.0d0*Ah/dy**2)*D_plume(i,j)*vctr(i,j)
                      write(iulog,*) '   N term:', (2.0d0*Ah/dy**2)*D_plume(i,j+1)*vctr(i,j+1)
                   endif
                endif
             enddo
          enddo

          call parallel_halo(Avisc_north_rhsu, parallel)
          call parallel_halo(Avisc_north_rhsv, parallel)

       endif   ! Ah > 0

       !--------------------------------------------------------------------
       ! Compute velocity on east edges
       !--------------------------------------------------------------------

       if (verbose_plume .and. main_task) then
          write(iulog,*) ' '
          write(iulog,*) 'iter_velo =', iter_velo
          write(iulog,*) 'compute east edge velocities'
       endif

       call plume_velocity(&
            nx,    ny,               &
            itest, jtest, rtest,     &
            'east',                  &
            edge_mask_east,          &
            wall_factor_east,        &
            D_plume_east,            &
            pgf_x_east,              &
            pgf_y_east,              &
            Avisc_east_lhs,          &
            Avisc_east_rhsu,         &
            Avisc_east_rhsv,         &
            u_plume_east,            &
            v_plume_east,            &
            converged_velo_east)

       !--------------------------------------------------------------------
       ! Compute velocity on north edges
       !--------------------------------------------------------------------

       if (verbose_plume .and. this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'compute north edge velocities'
       endif

       call plume_velocity(&
            nx,    ny,               &
            itest, jtest, rtest,     &
            'north',                 &
            edge_mask_north,         &
            wall_factor_north,       &
            D_plume_north,           &
            pgf_x_north,             &
            pgf_y_north,             &
            Avisc_north_lhs,         &
            Avisc_north_rhsu,        &
            Avisc_north_rhsv,        &
            u_plume_north,           &
            v_plume_north,           &
            converged_velo_north)

       ! check for convergence in all cells

       count_east = 0
       count_north = 0

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (edge_mask_east(i,j) > 0 .and. .not.converged_velo_east(i,j) ) then
                count_east = count_east + 1
             endif
             if (edge_mask_north(i,j) > 0 .and. .not.converged_velo_north(i,j) ) then
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

    call parallel_halo(u_plume_east, parallel)
    call parallel_halo(v_plume_east, parallel)
    call parallel_halo(u_plume_north, parallel)
    call parallel_halo(v_plume_north, parallel)

    if (verbose_plume) then
       call point_diag(u_plume_east, 'u_plume_east (m/s)', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(v_plume_east, 'v_plume_east (m/s)', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(u_plume_north, 'u_plume_north (m/s)', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(v_plume_north, 'v_plume_north (m/s)', itest, jtest, rtest, wx, wy, '(f10.5)')
    endif

  end subroutine compute_plume_velocity

!****************************************************

  subroutine compute_edge_gradients(&
       nx,               ny,              &
       dx,               dy,              &
       parallel,                          &
       plume_mask,                        &
       edge_mask_east,   edge_mask_north, &
       field,                             &
       df_dx_east,       df_dy_east,      &
       df_dx_north,      df_dy_north,     &
       neumann_bc)

    !----------------------------------------------------------------------------
    ! Compute horizontal gradients at the east and north edges of each cell.
    ! There are two ways to compute gradients at boundaries (edges with plume_mask = 1
    !  on one side and plume_mask = 0 on the other):
    ! (1) Compute gradients in the usual way, assuming that field values outside
    !     the plume are valid. This may be done for lsrf_plume.
    ! (2) Neumann BC: Set gradients to zero at the boundary. This is done for drho_plume.
    !----------------------------------------------------------------------------

    ! input/output arguments

    integer, intent(in) ::  &
         nx,     ny                  ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy                  ! grid cell size (m)

    type(parallel_type), intent(in) :: &
         parallel                    ! info for parallel communication

    integer, dimension(nx,ny), intent(in) :: &
         plume_mask,               & ! = 1 for cells where plume gradients are computed
         edge_mask_east,           & ! boundary type for east edges (regular, open or closed)
         edge_mask_north             ! boundary type for north edges (regular, open or closed)

    real(dp), dimension(nx,ny), intent(in) :: &
         field                       ! input scalar field

    real(dp), dimension(nx,ny), intent(out) :: &
         df_dx_east,  df_dy_east,  & ! gradients on east edges
         df_dx_north, df_dy_north    ! gradients on north edges

    logical, intent(in), optional :: &
         neumann_bc                  ! if true, apply a Neumann BC (zero gradient at boundaries)

    ! local variables

    integer :: i, j

    logical :: neumann               ! local version of neumann_bc

    integer :: wt_nw, wt_ne, wt_sw, wt_se   ! binary weights for neighboring edges
    integer :: count

    if (present(neumann_bc)) then
       neumann = neumann_bc
    else
       neumann = .false.
    endif

    ! initialize
    df_dx_east = 0.0d0
    df_dy_east = 0.0d0
    df_dx_north = 0.0d0
    df_dy_north = 0.0d0

    !----------------------------------------------------------------------------
    ! Compute x gradients on east edges and y gradients on north edges
    !----------------------------------------------------------------------------

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          if (edge_mask_east(i,j) == 1) then  ! regular boundary
             df_dx_east(i,j) = (field(i+1,j) - field(i,j)) / dx
          elseif (edge_mask_east(i,j) == 2 .or. edge_mask_east(i,j) == 3) then   ! open or closed boundary
             if (neumann) then
                df_dx_east(i,j) = 0.0d0
             else
                df_dx_east(i,j) = (field(i+1,j) - field(i,j)) / dx
             endif
          endif

          if (edge_mask_north(i,j) == 1) then
             df_dy_north(i,j) = (field(i,j+1) - field(i,j)) / dy
          elseif (edge_mask_north(i,j) == 2 .or. edge_mask_north(i,j) == 3) then   ! open or closed boundary
             if (neumann) then
                df_dy_north(i,j) = 0.0d0
             else
                df_dy_north(i,j) = (field(i,j+1) - field(i,j)) / dy
             endif
          endif

       enddo
    enddo

    call parallel_halo(df_dx_east, parallel)
    call parallel_halo(df_dy_north, parallel)

    !----------------------------------------------------------------------------
    ! Average neighboring values to get y gradients on east edges and x gradients on north edges.
    ! That is, df/dy on east edges is the average of df/dy on the neighboring north edges,
    !  and df/dx on north edges is the average of df/dx on the neighboring east edges.
    ! The average excludes edges along boundaries.
    !
    !           N     |     N                |             |
    !        ---|-----------|---           E---           ---E
    !                 |                      |             |
    !                 |                      |             |
    !               E---                  ---|------|------|---
    !                 |                      |      N      |
    !                 |                      |             |
    !        ---|-----------|---           E---           ---E
    !           N     |     N                |             |
    !----------------------------------------------------------------------------

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          if (edge_mask_east(i,j) /= 0) then   ! plume cell on at least one side
             ! assign a weight of 0 or 1 to each of the four neighboring north edges
             wt_nw = 0; wt_ne = 0; wt_sw = 0; wt_se = 0
             if (edge_mask_north(i,j) == 1) wt_nw = 1
             if (edge_mask_north(i+1,j) == 1) wt_ne = 1
             if (edge_mask_north(i,j-1) == 1) wt_sw = 1
             if (edge_mask_north(i+1,j-1) == 1) wt_se = 1
             count = wt_nw + wt_ne + wt_sw + wt_se
             if (count > 0) then
                df_dy_east(i,j) = &
                     (wt_nw*df_dy_north(i,j)   + wt_ne*df_dy_north(i+1,j) &
                    + wt_sw*df_dy_north(i,j-1) + wt_se*df_dy_north(i+1,j-1)) / count
             endif
          endif

          if (edge_mask_north(i,j) /= 0) then   ! plume cell on at least one side
             ! assign a weight of 0 or 1 to each of the four neighboring east edges
             wt_nw = 0; wt_ne = 0; wt_sw = 0; wt_se = 0
             if (edge_mask_east(i-1,j+1) == 1) wt_nw = 1
             if (edge_mask_east(i,j+1) == 1) wt_ne = 1
             if (edge_mask_east(i-1,j) == 1) wt_sw = 1
             if (edge_mask_east(i,j) == 1) wt_se = 1
             count = wt_nw + wt_ne + wt_sw + wt_se
             if (count > 0) then
                df_dx_north(i,j) = &
                     (wt_nw*df_dy_east(i-1,j+1) + wt_ne*df_dy_east(i,j+1) &
                    + wt_sw*df_dy_east(i-1,j)   + wt_se*df_dy_east(i,j)) / count
             endif
          endif

       enddo   ! i
    enddo   ! j

    call parallel_halo(df_dy_east, parallel)
    call parallel_halo(df_dx_north, parallel)

  end subroutine compute_edge_gradients

!****************************************************

  subroutine plume_velocity(&
       nx,    ny,               &
       itest, jtest, rtest,     &
       edge,                    &
       edge_mask,               &
       wall_factor,             &
       D_plume_edge,            &
       pgf_x,                   &
       pgf_y,                   &
       Avisc_lhs,               &
       Avisc_rhsu,              &
       Avisc_rhsv,              &
       u_plume,                 &
       v_plume,                 &
       converged_velo)

    !TODO - Compute Avisc_rhs locally based on input velocities?

    !--------------------------------------------------------------------
    ! Compute the velocity on a set of edges (either east or north)
    ! The calculation is independent of whether these are east or north edges
    !  except at closed boundaries, where the perpendicular component is set to zero.
    !--------------------------------------------------------------------

    ! input/output arguments

    integer, intent(in) ::  &
         nx,  ny,           & ! number of grid cells in each dimension
         itest, jtest, rtest  ! test cell coordinates (diagnostic only)

    character(len=*), intent(in) ::   &
         edge                 ! 'east' or 'north'

    integer, dimension(nx,ny), intent(in) ::   &
         edge_mask            ! = 1 at edges where velocity is computed

    real(dp), dimension(nx,ny), intent(in) ::   &
         wall_factor,       & ! factor to reduce flow toward the boundary on edges (unitless)
         D_plume_edge,      & ! plume thickness at edges (m)
         pgf_x,             & ! x component of pressure gradient force, at edges (m^2/s^2)
         pgf_y,             & ! y component of pressure gradient force, at edges (m^2/s^2)
         Avisc_lhs,         & ! viscous coefficient on LHS; multiplies the edge velocity (m/s)
         Avisc_rhsu,        & ! viscous term on RHS of u equation (m^2/s^2)
         Avisc_rhsv           ! viscous term on RHS of v equation (m^2/s^2)

    ! Note: u and v are colocated on either east edges or north edges,
    !       depending on the subroutine call
    real(dp), dimension(nx,ny), intent(inout) ::  &
         u_plume,           & ! x component of plume velocity (m/s) on the edge
         v_plume              ! y component of plume velocity (m/s) on the edge

    logical, dimension(nx,ny), intent(out) ::  &
         converged_velo       ! true when velocity has converged at an edge, else false

    ! local variables

    real(dp) :: &
         speed,             & ! plume speed (m/s), updated at each iteration until convergence
         f_x, f_y,          & ! combined PGF and viscosity terms (m^2/s^2)
         f_cor,             & ! product of f_coriolis and wall_factor; goes to 0 along closed boundaries
         cUA,               & ! c_drag*speed + Avisc_lhs (m/s); multiplies the local velocity
         x_resid, y_resid,  & ! residuals of momentum balance equations (m^2/s^2)
         resid,             &
         u_new, v_new,      & ! guesses for new u and v
         denom,             & ! denominator
         a_uu, a_uv,        & ! coefficients for Newton solve
         a_vu, a_vv,        & !
         du, dv               ! change in u and v (m/s)

    character(len=128) :: message

    real(dp), parameter :: &
         maxresid_force_balance = 1.0d-8   ! max residual allowed in momentum balance equation (m^2/s^2)

    !TODO - Start with Picard, then test Newton
    logical, parameter :: &
!         velo_newton = .false.  ! if true, use Newton's method; if false, use Picard method
         velo_newton = .true.  ! if true, use Newton's method; if false, use Picard method

    integer :: i, j

    !--------------------------------------------------------------------
    ! Compute the plume velocity, assuming a balance between the pressure gradient force,
    !  basal drag, Coriolis force and horizontal viscosity:
    !
    ! (1) pgf_x - c*|U|*u + f*D*v + del*(Ah*D*grad(u)) = 0
    ! (2) pgf_y - c*|U|*v - f*D*u + del*(Ah*D*grad(v)) = 0
    !
    !            D = plume thickness (m)
    !          pgf = pressure gradient force (m^2/s^2)
    !            c = dimensionless ocean drag coefficient
    !          |U| = sqrt(u^2 + v^2 + u_tidal^2)
    !      u_tidal = uniform small velocity added for regularization
    !            f = Coriolis coefficient (1/s)
    !           Ah = uniform horizontal viscosity coefficient (m^2/s)
    !
    ! On edges adjacent to boundaries, the Coriolis terms are multiplied by a factor of 0.0 or 0.5
    ! (depending on whether there are two adjacent closed boundaries or just one).
    ! This reduces the Coriolis-driven flow toward the boundary, allowing a PGF-driven jet to form
    !  along the boundary.
    !
    ! The viscosity is small in much of the domain but can be large near closed boundaries.
    ! The viscosity terms at east edges are discretized as
    !
    ! del*(Ah*D*grad(u)) = Ah*[d/dx(D du/dx) + d/dy(D du/dy)]
    !                    = (Ah/dx^2) * [D(i+1,j)*(uctr(i+1,j) - ueast(i,j))/(dx/2) - D(i,j)*(ueast(i,j) - uctr(i-1,j))/(dx/2)]
    !                    + (Ah/dy^2) * [stagD(i,j)*(stagu(i,j) - ueast(i,j))/(dy/2) - stagD(i,j-1)*(ueast(i,j) - stagu(i,j-1))/(dy/2)]
    ! where stagD denotes D on the staggered grid.
    ! The v expressions are analogous.
    !
    ! The viscosity terms at north edges are discretized as
    ! del*(Ah*D*grad(u)) = (Ah/dx^2) * [stagD(i,j)*(stagu(i,j) - unorth(i,j))/(dx/2) - stagD(i-1,j)*(unorth(i,j) - stagu(i-1,j))/(dx/2)]
    !                    + (Ah/dy^2) * [D(i,j+1)*(uctr(i,j+1) - unorth(i,j))/(dy/2) - D(i,j)*(unorth(i,j) - uctr(i,j))/(dy/2)]
    !
    ! We can rewrite (1) and (2) as
    !
    ! (1) pgf_x - c|U|*u + f*D*v - A_lhs*u = -A_rhsu
    ! (2) pgf_y - c|U|*v - f*D*u - A_lhs*v = -A_rhsv
    !
    ! where A_lhs is a viscosity term that multiplies the local edge velocity,
    !  and A_rhsu and A_rhsv are viscosity terms that include neighboring velocities.
    ! The velocities in the RHS terms are lagged by one iteration.
    !
    ! With Ah = 0, the u and v solutions are
    !
    !                c|U|*pgf_x + f*D*pgf_y
    !            u = ________________________
    !                 (f*D)^2 + (c|U|)^2
    !
    !                c|U|*pgf_y - f*D*pgf_x
    !            v = ________________________
    !                 (f*D)^2 + (c|U|)^2
    !
    ! With nonzero Ah, the solutions are
    !
    !                (c|U| + A_lhs)*(pgf_x + A_rhsu) + f*D*(pgf_y + A_rhsv)
    !            u = ______________________________________________________
    !                            (f*D)^2 + (c|U| + A_lhs)^2
    !
    !                (c|U| + A_lhs)*(pgf_y + A_rhsv) - f*D*(pgf_x + A_rhsu)
    !            v = ______________________________________________________
    !                            (f*D)^2 + (c|U| + A_lhs)^2
    !
    ! where the A_rhs terms include the dependence on neighboring velocity components.
    ! We iterate the solution to convergence, updating U, A_rhsu and A_rhsv after each iteration.
    ! The iterative loop is in the subroutine above.
    ! This subroutine computes the solution for the current iteration.
    !
    ! At closed boundaries on east edges, the parallel (v) component can be nonzero,
    !  but the perpendicular (u) component is zero. The solution for v is
    !
    !                  c|U| + A_lhs
    !            v =  ______________
    !                 pgf_y + A_rhsv
    !
    ! At closed boundaries on north edges, the parallel (u) component can be nonzero,
    !  but the perpendicular (v) component is zero. The solution for u is
    !
    !                  c|U| + A_lhs
    !            u =  ______________
    !                 pgf_x + A_rhsu
    !
    ! Convergence can be sped up using Newton's method.
    ! We write   u = u0 + du
    !            v = v0 + dv
    !          |U| = U0 + (d|U|/du)*du + (d|U|dv)*dv
    ! where '0' denotes the current guess for the solution,
    ! and the partial derivatives are evaluated at (u,v) = (u0,v0).
    ! 
    ! By using a Taylor expansion and dropping higher-order terms, it can be shown that
    !           du = (a_vv*R_x - a_uv*R_y) / det|M|
    !           dv = (a_uu*R_y - a_vu*R_x) / det|M|
    ! where    
    !          R_x = pgf_x - (c*U0 + A_lhs)*u0 + f*D*v0 + A_rhsu = x residual
    !          R_y = pgf_y - (c*U0 + A_lhs)*v0 - f*D*u0 + A_rhsv = y residual
    !
    !                | m_uu   m_uv |
    ! and        M = |             |
    !                | m_vu   m_vv |
    !
    ! with    m_uu = c*(U0 + u0^2/U0) + A_lhs
    !         m_uv = c*u0*v0/U0 - f*D
    !         m_vu = c*u0*v0/U0 + f*D
    !         m_vv = c*(U0 + v0^2/U0) + A_lhs
    !--------------------------------------------------------------------

    ! Compute the u and v velocity components at each edge,
    ! with the plume speed and the RHS viscosity terms lagged by one iteration.

    ! Loop over edges
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
             write(iulog,*) 'edge, starting u/v plume (m/s) =', trim(edge), u_plume(i,j), v_plume(i,j)
          endif

          ! Compute the plume speed based on the input u and v
          speed = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)

          ! Set the Coriolis term, reduced as needed near closed boundaries
          f_cor = f_coriolis * wall_factor(i,j)

          ! Combine the PGF and RHS viscosity terms
          f_x = pgf_x(i,j) + Avisc_rhsu(i,j)
          f_y = pgf_y(i,j) + Avisc_rhsv(i,j)

          ! Combine the drag and LHS viscosity terms
          cUA = c_drag*speed + Avisc_lhs(i,j)

          if (edge_mask(i,j) == 3) then   ! closed boundary; solve for one component only

             if (trim(edge) == 'east') then   ! u = 0; solve for v
                y_resid = f_y - cUA*v_plume(i,j)
                if (abs(y_resid) < maxresid_force_balance) converged_velo(i,j) = .true.
                if (velo_newton) then
                   a_vv = c_drag*(speed + v_plume(i,j)**2/speed) + Avisc_lhs(i,j)
                   dv = y_resid / a_vv
                   v_plume(i,j) = v_plume(i,j) + dv
                else  ! Picard solve
                   ! actually not a Picard solve; a standard Picard solve is very slow to converge
                   v_new = f_y / cUA
                   v_plume(i,j) = 0.5d0 * (v_plume(i,j) + v_new)
                endif
             elseif (trim(edge) == 'north') then   ! v = 0; solve for u
                x_resid = f_x - cUA*u_plume(i,j)
                if (abs(x_resid) < maxresid_force_balance) converged_velo(i,j) = .true.
                if (velo_newton) then
                   a_uu = c_drag*(speed + u_plume(i,j)**2/speed) + Avisc_lhs(i,j)
                   du = x_resid / a_uu
                   u_plume(i,j) = u_plume(i,j) + du
                else  ! Picard solve
                   ! actually not a Picard solve; a standard Picard solve is very slow to converge
                   u_new = f_x / cUA
                   u_plume(i,j) = 0.5d0 * (u_plume(i,j) + u_new)
                endif
             endif   ! east or north edge

          elseif (edge_mask(i,j) > 0) then  ! regular or open boundary; solve for both components

             ! Compute the residual of the u and v equations and check convergence
             x_resid = f_x - cUa*u_plume(i,j) + f_cor*D_plume_edge(i,j)*v_plume(i,j)
             y_resid = f_y - cUa*v_plume(i,j) - f_cor*D_plume_edge(i,j)*u_plume(i,j)
             resid = sqrt(x_resid**2 + y_resid**2)
             if (resid < maxresid_force_balance) converged_velo(i,j) = .true.

             if (velo_newton) then

                ! compute some coefficients for the Newton solve
                a_uu = c_drag*(speed + u_plume(i,j)**2/speed) + Avisc_lhs(i,j)
                a_vv = c_drag*(speed + v_plume(i,j)**2/speed) + Avisc_lhs(i,j)
                a_uv = c_drag*(u_plume(i,j)*v_plume(i,j))/speed - f_cor*D_plume_edge(i,j)
                a_vu = c_drag*(u_plume(i,j)*v_plume(i,j))/speed + f_cor*D_plume_edge(i,j)

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

             else  ! Picard solve

                denom = (f_cor*D_plume_edge(i,j))**2 + cUA**2
                u_plume(i,j) = (cUA*f_x + f_cor*D_plume_edge(i,j)*f_y) / denom
                v_plume(i,j) = (cUA*f_y - f_cor*D_plume_edge(i,j)*f_x) / denom

             endif  ! Newton or Picard

          endif  ! edge_mask

          if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
             write(iulog,*) 'edge, i, j:', trim(edge), i, j
             write(iulog,*) 'speed (m/s) =', speed
             write(iulog,*) 'edgeD:', D_plume_edge(i,j)
             write(iulog,*) 'pgf_x, pgf_y(m2/s2):', pgf_x(i,j), pgf_y(i,j)
             write(iulog,*) 'f*D (m/s):', f_cor*D_plume_edge(i,j)
             write(iulog,*) 'fDv, fDu (m2/s2):', f_cor*D_plume_edge(i,j)*v_plume(i,j), &
                  f_cor*D_plume_edge(i,j)*u_plume(i,j)
             write(iulog,*) 'c|U| (m/s):', c_drag*speed
             write(iulog,*) 'Fdu, Fdv (m2/s2):', -c_drag*speed*u_plume(i,j), -c_drag*speed*v_plume(i,j)
             write(iulog,*) 'Avisc_lhs (m/s):', Avisc_lhs(i,j)
             write(iulog,*) 'LHS Fviscu, Fviscv (m2/s2):', -Avisc_lhs(i,j)*u_plume(i,j), -Avisc_lhs(i,j)*v_plume(i,j)
             write(iulog,*) 'RHS Fviscu, Fviscv (m2/s2):', Avisc_rhsu(i,j), Avisc_rhsv(i,j)
             write(iulog,*) 'x/y residual (m2/s2):', x_resid, y_resid
             write(iulog,*) 'new u/v_plume (m/s):', u_plume(i,j), v_plume(i,j)
             write(iulog,*) 'converged =', converged_velo(i,j)
          endif

       enddo  ! i
    enddo  ! j

  end subroutine plume_velocity

!****************************************************

  subroutine plume_entrainment(&
       nx,         ny,      &
       dx,         dy,      &
       itest, jtest, rtest, &
       parallel,            &
       which_entrainment,   &
       plume_mask,          &
       edge_mask_east,      &
       edge_mask_north,     &
       lsrf_plume,          &
       plume_speed,         &
       bmlt_float,          &
       drho_plume,          &
       drho_basal,          &
       H_cavity,            &
       D_plume,             &
       dt_plume,            &
       entrainment,         &
       detrainment)

    !--------------------------------------------------------------------
    ! Compute plume entrainment and detrainment by one of several methods
    !--------------------------------------------------------------------

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    integer, intent(in) :: &
         which_entrainment      ! entrainment option

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask,          & ! = 1 for cells where scalar plume variables are computed
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed
         edge_mask_north        ! = 1 on north edges where plume velocity is computed

    real(dp), dimension(nx,ny), intent(in) ::  &
         lsrf_plume,          & ! elevation of plume-ambient interface (m)
         plume_speed,         & ! plume speed at cell center (m/s)
         bmlt_float,          & ! melt rate (m/s)
         drho_plume,          & ! density difference between ambient ocean and plume (kg/m3)
         drho_basal,          & ! density difference between plume and ice base (kg/m3)
         H_cavity,            & ! thickness of sub-shelf cavity (m)
         D_plume                ! plume thickness
    
    real(dp), intent(in) :: &
         dt_plume               ! timestep (s)

    !Note: Both entrainment and detrainment are >= 0 by definition
    real(dp), dimension(nx,ny), intent(out) ::  &
         entrainment,           & ! entrainment at cell centers (m/s)
         detrainment              ! detrainment at cell centers (m/s)

    ! local variables

    real(dp), dimension(nx,ny) ::  &
         ufric_plume,           & ! friction velocity (m/s)
         dlsrf_plume_dx_east,   & ! horizontal gradient of lsrf_plume on east edges
         dlsrf_plume_dy_east,   & !
         dlsrf_plume_dx_north,  & ! horizontal gradient of lsrf_plume on north edges
         dlsrf_plume_dy_north,  & !
         theta_slope              ! basal slope angle at cell centers (rad)

    real(dp) :: &
         dlsrf_plume_dx,        & ! lsrf gradient components at cell centers
         dlsrf_plume_dy,        &
         slope,                 & ! magnitude of the gradient (dlsrf_dx, dlsrf_dy)
         entrainment_min,       & ! min entrainment rate when D_plume < D_plume_min
         detrainment_min          ! min detrainment rate when D_plume > D_plume_max

    integer :: i, j, ig, jg
    real(dp) :: numer, denom

    ! entrainment parameters
    real(dp), parameter ::   &
!!         E0 = 0.072d0             ! entrainment coefficient (unitless)
         E0 = 0.036d0             ! entrainment coefficient (unitless) for the Jenkins (1991) scheme
                                  ! Bo Pederson (1980) suggests E0 = 0.072
                                  ! Jenkins (1991, JGR) suggests 0.036 to compensate for lack of Coriolis in 1D model
    real(dp), parameter ::    &
         mu_e = 2.5d0             ! nondimensional parameter for the Gaspar (1988) scheme
                                  ! Gaspar (1988) and Lambert et al. (2023) set mu = 0.5;
                                  ! Gladish et al. (2012) and Lambert et al. (2026) set mu = 2.5

    logical, parameter :: verbose_entrainment = .true.

    entrainment = 0.0d0
    detrainment = 0.0d0

    ! Given plume_speed, compute ufric_plume
    ufric_plume = sqrt(c_drag)*plume_speed

    if (which_entrainment == PLUME_ENTRAINMENT_JENKINS) then

       !--------------------------------------------------------------------
       ! Compute entrainment as a function of the plume speed and the slope of the
       !  plume-ambient interface, following Bo Pederson (1980) and Jenkins (1991).
       !
       ! entrainment = E0 * plume_speed * sin(theta_slope)
       !
       ! Note: plume_speed is proportional to ufric_plume, so we could replace
       !       one with the other and rescale the constant.
       !--------------------------------------------------------------------

       ! Compute the slope angle at cell centers.
       ! This requires computing gradients of lsrf_plume at cell edges
       !  and averaging the gradients to cell centers.
       ! Set gradients to zero at boundaries to avoid large entrainment near boundaries.

       call compute_edge_gradients(&
            nx,                   ny,                    &
            dx,                   dy,                    &
            parallel,                                    &
            plume_mask,                                  &
            edge_mask_east,       edge_mask_north,       &
            lsrf_plume,                                  &
            dlsrf_plume_dx_east,  dlsrf_plume_dy_east,   &
            dlsrf_plume_dx_north, dlsrf_plume_dy_north,  &
            neumann_bc = .true.)

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

       if (verbose_entrainment) then
          call point_diag(theta_slope, 'theta_slope (rad)', itest, jtest, rtest, wx, wy)
       endif

       ! Eq. 6 from Jenkins (1991)
       entrainment = E0 * plume_speed * sin(theta_slope)

       if (verbose_entrainment .and. this_rank == rtest) then
          i = itest
          j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'Jenkins entrainment, rank, i, j =', this_rank, i, j
          write(iulog,*) 'D (m), sin(theta), speed (m/s):', D_plume(i,j), sin(theta_slope(i,j)), plume_speed(i,j)
          write(iulog,*) 'entrainment =', entrainment(i,j)
       endif

    elseif (which_entrainment == PLUME_ENTRAINMENT_GASPAR) then

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

       ! loop over all cells
       do j = 1, ny
          do i = 1, nx
             if (plume_mask(i,j) == 1) then
                numer = mu_e * ufric_plume(i,j)**3 - 0.5d0*D_plume(i,j)*(grav/rhoo)*drho_basal(i,j)*bmlt_float(i,j)
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

                if (verbose_entrainment .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) ' '
                   write(iulog,*) 'Gaspar entrainment, rank, i, j =', this_rank, i, j
                   write(iulog,*) 'D (m), u_fric (m/s):', D_plume(i,j), ufric_plume(i,j)
                   write(iulog,*) 'drho_plume (kg/m3), drho_basal, bmlt (m/yr):', &
                        drho_plume(i,j), drho_basal(i,j), bmlt_float(i,j)
                   write(iulog,*) 'ufric term (m3^s3):', mu_e * ufric_plume(i,j)**3
                   write(iulog,*) 'bmlt term  (m3/s3):', 0.5d0*D_plume(i,j)*(grav/rhoo)*drho_basal(i,j)*bmlt_float(i,j)
                   if (entrainment(i,j) > 0.0d0) then
                      write(iulog,*) 'entrainment =', entrainment(i,j)
                   else
                      write(iulog,*) 'detrainment =', detrainment(i,j)
                   endif
                endif

             endif
          enddo   ! i
       enddo   ! j

    endif  ! which_entrainment

    if (verbose_plume) then
       call point_diag(ufric_plume, 'ufric_plume (m/s)', itest, jtest, rtest, wx, wy, '(f10.5)')
       call point_diag(plume_speed, 'plume_speed (m/s)', itest, jtest, rtest, wx, wy)
       call point_diag(entrainment*scyr, 'Before adjusting, entrainment (m/yr)', itest, jtest, rtest, wx, wy)
       call point_diag(detrainment*scyr, 'Before adjusting, detrainment (m/yr)', itest, jtest, rtest, wx, wy)
    endif

    !--------------------------------------------------------------------
    ! Reduce entrainment in thin cavities.
    ! Entrainment = 0 for D_plume >= H_cavity
    !--------------------------------------------------------------------

    if (plume_cavity_h0 > 0.0d0) then
       do j = 1, ny
          do i = 1, nx
             if (plume_mask(i,j) == 1 .and. H_cavity(i,j) - D_plume(i,j) < plume_cavity_h0) then
                entrainment(i,j) = entrainment (i,j) * max(0.0d0, (H_cavity(i,j) - D_plume(i,j))/plume_cavity_h0)
             endif
          enddo
       enddo
    endif

    !--------------------------------------------------------------------
    ! Adjust entrainment or detrainment if D_plume is outside a desired range.
    !--------------------------------------------------------------------

    do j = 1, ny
       do i = 1, nx
          if (plume_mask(i,j) == 1) then

             ! Increase entrainment if D_plume < D_plume_min
             if (D_plume(i,j) < D_plume_min) then
                if (verbose_entrainment) then
                   call parallel_globalindex(i, j, ig, jg, parallel)
                   write(iulog,*) 'Force entrainment: ig, jg, D_plume:', ig, jg, D_plume(i,j)
                endif
                entrainment_min = (D_plume_min - D_plume(i,j)) / tau_relax_entrainment
                if (detrainment(i,j) > 0.0d0) then
                   entrainment_min = entrainment_min - detrainment(i,j)
                   detrainment(i,j) = 0.0d0
                endif
                entrainment(i,j) = max(entrainment(i,j), entrainment_min)
             endif

             ! Increase detrainment if D_plume > D_plume_max
             if (D_plume(i,j) > D_plume_max) then
                if (verbose_entrainment) then
                   call parallel_globalindex(i, j, ig, jg, parallel)
                   write(iulog,*) 'Force detrainment: ig, jg, D_plume:', ig, jg, D_plume(i,j)
                endif
                detrainment_min = (D_plume(i,j) - D_plume_max) / tau_relax_entrainment
                detrainment(i,j) = max(detrainment(i,j), detrainment_min)
             endif

          endif
       enddo   ! i
    enddo   ! j

    if (verbose_plume) then
       call point_diag(H_cavity - D_plume, 'H_cavity - D_plume (m)', itest, jtest, rtest, wx, wy)
       call point_diag(entrainment*scyr, 'After adjusting, entrainment (m/yr)', itest, jtest, rtest, wx, wy)
       call point_diag(detrainment*scyr, 'After adjusting, detrainment (m/yr)', itest, jtest, rtest, wx, wy)
    endif

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
       ufric_plume,         &
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
         ufric_plume,         & ! plume friction velocity (m/s) on ice grid, output as a diagnostic
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
             C1 = (rhoo * cpw * ufric_plume(i,j) * gammaT) / (rhoi * lhci)
             C2 = (rhoo * ufric_plume(i,j) * gammaS) / rhoi
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
                write(iulog,*) 'Eq. 1 RHS:', rhoo*cpw*ufric_plume(i,j)*gammaT*(T_plume(i,j) - T_basal(i,j))
                write(iulog,*) 'Eq. 2 LHS:', rhoi*bmlt_float(i,j)*S_basal(i,j)
                write(iulog,*) 'Eq. 2 RHS:', rhoo*ufric_plume(i,j)*gammaS*(S_plume(i,j) - S_basal(i,j))
                write(iulog,*) 'Eq. 3 LHS:', T_basal(i,j)
                write(iulog,*) 'Eq. 3 RHS:', lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
             endif
          endif   ! plume_mask = 1
          
       enddo   ! i
    enddo   ! j

  end subroutine plume_melt_rate

!****************************************************

  subroutine compute_transport_velocity(&
       nx,               ny,               &
       itest,   jtest,   rtest,            &
       parallel,                           &
       edge_mask_east,   edge_mask_north,  &
       u_plume_east,     v_plume_east,     &
       u_plume_north,    v_plume_north,    &
       u_transport_east, v_transport_north)

    !----------------------------------------------------------------------------
    ! Given u and v computed at each east and north edge, compute the C-grid velocities
    !  to be used for transport.
    ! This is done by taking a weighted average of the local velocity at the edge,
    !  combined with the velocities at up to four neighboring edges.
    ! The goal is to obtain average velocities that represent the overall flow
    !  more accurately than the local velocities.
    !----------------------------------------------------------------------------

    ! input/output arguments

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    integer, dimension(nx,ny), intent(in) ::  &
!         plume_mask,          & ! = 1 for cells where the plume is present, else = 0
         edge_mask_east,       & ! = 1 for east edges with plume cells on each side
         edge_mask_north         ! = 1 for north edges with plume cells on each side
    
    real(dp), dimension(nx,ny), intent(in) ::  &
         u_plume_east,         & ! u_plume on east edges (m/s)
         v_plume_east,         & ! v_plume on east edges (m/s)
         u_plume_north,        & ! u_plume on north edges (m/s)
         v_plume_north           ! v_plume on north edges (m/s)

    real(dp), dimension(nx,ny), intent(out) ::  &
         u_transport_east,     & ! transport velocity u on east edges (m/s)
         v_transport_north       ! transport velocity v on north edges (m/s)

    ! local variables

    integer :: i, j
    integer :: wt_nw, wt_ne, wt_sw, wt_se   ! binary weights for neighboring edges
    integer :: count

    real(dp) :: neighbor_average    ! average velocity over four neighboring edges

    u_transport_east = 0.0d0
    v_transport_north = 0.0d0

    !----------------------------------------------------------------------------
    ! Compute the average velocities over neighboring edges and combine with
    !  the velocity at the local edge.
    ! Neighboring edges outside the plume domain are excluded.
    !----------------------------------------------------------------------------

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          ! east edges
          if (edge_mask_east(i,j) > 0) then   ! plume cell on at least one side
             ! assign a weight of 0 or 1 to each of the four neighboring north edges
             wt_nw = 0; wt_ne = 0; wt_sw = 0; wt_se = 0
!             if (edge_mask_north(i,j) == 1) wt_nw = 1
!             if (edge_mask_north(i+1,j) == 1) wt_ne = 1
!             if (edge_mask_north(i,j-1) == 1) wt_sw = 1
!             if (edge_mask_north(i+1,j-1) == 1) wt_se = 1
             if (edge_mask_north(i,j) > 0) wt_nw = 1
             if (edge_mask_north(i+1,j) > 0) wt_ne = 1
             if (edge_mask_north(i,j-1) > 0) wt_sw = 1
             if (edge_mask_north(i+1,j-1) > 0) wt_se = 1
             count = wt_nw + wt_ne + wt_sw + wt_se
             if (count > 0) then
                ! compute the weighted average velocity
                neighbor_average = &
                     (wt_nw*u_plume_north(i,j)   + wt_ne*u_plume_north(i+1,j) &
                    + wt_sw*u_plume_north(i,j-1) + wt_se*u_plume_north(i+1,j-1)) / count
                u_transport_east(i,j) = &
                     (4*u_plume_east(i,j) + count*neighbor_average) / (4 + count)
             else
                ! use the local velocity
                u_transport_east(i,j) = u_plume_east(i,j)
             endif
          endif

          ! north edges
          if (edge_mask_north(i,j) > 0) then   ! plume cell on at least one side
             ! assign a weight of 0 or 1 to each of the four neighboring east edges
             wt_nw = 0; wt_ne = 0; wt_sw = 0; wt_se = 0
!             if (edge_mask_east(i-1,j+1) == 1) wt_nw = 1
!             if (edge_mask_east(i,j+1) == 1) wt_ne = 1
!             if (edge_mask_east(i-1,j) == 1) wt_sw = 1
!             if (edge_mask_east(i,j) == 1) wt_se = 1
             if (edge_mask_east(i-1,j+1) > 0) wt_nw = 1
             if (edge_mask_east(i,j+1) > 0) wt_ne = 1
             if (edge_mask_east(i-1,j) > 0) wt_sw = 1
             if (edge_mask_east(i,j) > 0) wt_se = 1
             count = wt_nw + wt_ne + wt_sw + wt_se
             if (count > 0) then
                ! compute the weighted average velocity
                neighbor_average = &
                     (wt_nw*v_plume_east(i-1,j+1) + wt_ne*v_plume_east(i,j+1) &
                    + wt_sw*v_plume_east(i-1,j)   + wt_se*v_plume_east(i,j)) / count
                v_transport_north(i,j) = &
                     (4*v_plume_north(i,j) + count*neighbor_average) / (4 + count)
             else
                ! use the local velocity
                v_transport_north(i,j) = v_plume_north(i,j)
             endif
          endif

       enddo   ! i
    enddo   ! j

    ! Zero out perpendicular velocities at closed boundaries.
    ! The input perpendicular components (u_plume_east on east edges
    ! and v_plume_north on north edges) are zero by construction,
    ! but the averaging above can introduce nonzero values.

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          if (edge_mask_east(i,j) == 3) then   ! closed east boundary
             u_transport_east(i,j) = 0.0d0
          endif
          if (edge_mask_north(i,j) == 3) then   ! closed north boundary
             v_transport_north(i,j) = 0.0d0
          endif
       enddo
    enddo

    call parallel_halo(u_transport_east, parallel)
    call parallel_halo(v_transport_north, parallel)

  end subroutine compute_transport_velocity

!****************************************************

  subroutine plume_transport(&
       nx,           ny,     &
       dx,           dy,     &
       itest, jtest, rtest,  &
       parallel,             &
       dt_plume,             &
       Kh,                   &
       plume_mask,           &
       edge_mask_east,       &
       edge_mask_north,      &
       u_east,               &
       v_north,              &
       entrainment,          &
       detrainment,          &
       bmlt_float,           &
       T_ambient,            &
       S_ambient,            &
       T_basal,              &
       S_basal,              &
       D_plume,              &
       T_plume,              &
       S_plume,              &
       divDu_plume)

    !----------------------------------------------------------------------------
    ! Solve transport equations for the plume thickness, temperature and salinity.
    ! These include horizontal transport of mass, heat and salt; horizontal diffusion
    !  of heat and salt; and vertical entrainment, detrainment and melting.
    !
    ! See Eqs. 1, 2 and 4 in Lambert et al. (2023):. These describe the conservation
    !
    ! (1) dD/dt+ del*(DU) = e - d + m
    !
    ! (2) d/dt(DT) + del*(DUT) = e*Ta -d*T + m*Tb - gammaT*(T - Tb) + del*(Kh*D*gradT)
    !
    ! (3) d/dt(DS) + del*(DUS) = e*Sa -d*S - del*(Kh*D*gradT)
    !
    ! where (T,S), (Tb,Sb) and (Ta,Sa) are the temperature and salinity of the plume,
    ! the ice base and the ambient ocean, respectively; D is the plume thickness;
    ! U = (u,v) is the velocity, with components perpendicular to the east and north edges;
    ! e, d and m are the rates of entraintment, detrainment and melting;
    ! gammaT is a heat transfer term; and Kh is the diffusivity of heat and salt.
    !
    ! Note that u and v are not necessarily to u_plume_east and v_plume_north
    !  as computed locally on each edge. For transport, we can compute u and v
    !  as weighted averages over several edges.
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
         dt_plume,            & ! plume time step (s)
         Kh                     ! horizontal plume_diffusivity (m^2/s);
                                ! assumed to be equal for heat and salt

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask,          & ! = 1 for cells where the plume is present, else = 0
         edge_mask_east,      & ! = 1 for east edges with plume cells on each side
         edge_mask_north        ! = 1 for north edges with plume cells on each side

    !TODO - Change names of u_east and v_north?
    real(dp), dimension(nx,ny), intent(in) ::  &
!         u_plume_east,        & ! u_plume on east edges (m/s)
!         v_plume_north,       & ! v_plume on north edges (m/s)
         u_east,              & ! transport velocity u on east edges (m/s)
         v_north,             & ! transport velocity v on north edges (m/s)
         entrainment,         & ! entrainment rate (m/s)
         detrainment,         & ! detrainment rate (m/s)
         bmlt_float,          & ! basal melt rate (m/s)
         T_ambient,           & ! ambient temperature (deg C)
         S_ambient,           & ! ambient salinity (psu)
         T_basal,             & ! basal temperature (deg C)
         S_basal                ! basal salinity (psu)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         D_plume,             & ! plume thickness (m)
         T_plume,             & ! plume temperature (deg C)
         S_plume,             & ! plume salinity (psu)
         divDu_plume            ! divergence of D*u (m/s)

    !WHL - debug
    real(dp), dimension(nx,ny) ::  &
         divDuT_plume

    ! local variables

    integer :: i, j, ig, jg, n
    integer :: ilo, ihi, jlo, jhi

    real(dp) :: dD, dDT, dDS           ! increments in D, D*T and D*S
    real(dp) :: gradT, gradS           ! gradients of heat and salt
    real(dp) :: heat_loss              ! rate of heat loss from the plume to the ice (m*deg/s)

    real(dp), dimension(nx,ny,3) :: &
         work                          ! work array for transport

    real(dp), dimension(nx,ny) ::  &
         diffT_east, diffT_north,    & ! diffusive fluxes of heat at cell edges (m^3*deg/s)
         diffS_east, diffS_north,    & ! diffusive fluxes of salt at cell edges (m^3*psu/s)
         flux_east,  flux_north        ! fluxes of D*u at each edge (m^2/s)

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

    !----------------------------------------------------------------------------
    ! Make sure the input fields are in range
    !----------------------------------------------------------------------------
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

    !----------------------------------------------------------------------------
    ! Fill a work array with the fields to be transported
    !----------------------------------------------------------------------------

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

    !----------------------------------------------------------------------------
    ! Increment the work array based on vertical entrainment, detrainment, melting
    ! and heat loss
    !----------------------------------------------------------------------------

    ! loop over locally owned cells
    do j = 1, ny
       do i = 1, nx
          if (plume_mask(i,j) == 1) then
             dD = entrainment(i,j) - detrainment(i,j) + bmlt_float(i,j)
             work(i,j,1) = work(i,j,1) + dD*dt_plume
             ! Compute the rate of heat transfer (J/m^2/s = W/m2) from the plume to the ice base.
             ! Note: rhoi*lhci*bmlt_float has units J/m^2/s, and rhow*cpw has units of (kg/m3)*J/(deg*kg) = J/(deg*m3),
             !       so heat_loss has units of m*deg/s, as desired
             !TODO - Compute in terms of (Tp - Tb)?
             heat_loss = rhoi*lhci*bmlt_float(i,j)/(rhoo*cpw)
             dDT = entrainment(i,j)*T_ambient(i,j) - detrainment(i,j)*T_plume(i,j) + bmlt_float(i,j)*T_basal(i,j) &
                  - heat_loss
             work(i,j,2) = work(i,j,2) + dDT*dt_plume
             ! Note: salt_transfer from the ice = 0
             dDS = entrainment(i,j)*S_ambient(i,j) - detrainment(i,j)*S_plume(i,j)
             work(i,j,3) = work(i,j,3) + dDS*dt_plume
          endif
       enddo
    enddo

    !----------------------------------------------------------------------------
    ! Compute horizontal transport using a first-order upwind scheme
    !----------------------------------------------------------------------------

    ! Set bounds for loops over locally owned cells and edges (inputs to glissade_upwind_field)
    ilo = nhalo + 1
    ihi = nx - nhalo
    jlo = nhalo +1
    jhi = ny - nhalo


    call glissade_upwind_field(&
         nx,             ny,             &
         ilo, ihi,       jlo, jhi,       &
         dx,             dy,             &
         dt_plume,       work(:,:,1),    &
         u_east,         v_north)

    ! Repeat for (D_plume*T_plume) and (D_plume*S_plume)

    call glissade_upwind_field(&
         nx,             ny,             &
         ilo, ihi,       jlo, jhi,       &
         dx,             dy,             &
         dt_plume,       work(:,:,2),    &
         u_east,         v_north)

    call glissade_upwind_field(&
         nx,             ny,             &
         ilo, ihi,       jlo, jhi,       &
         dx,             dy,             &
         dt_plume,       work(:,:,3),    &
         u_east,         v_north)

    ! halo update after horizontal transport
    do n = 1, 3
       call parallel_halo(work(:,:,n), parallel)
    enddo

    !----------------------------------------------------------------------------
    ! Compute diffusive fluxes of heat and salt at each plume edge.
    ! Note: There are no diffusive fluxes at open boundaries (edge_mask = 2),
    !       since we assume gradT = gradS = 0 at open boundaries.
    !----------------------------------------------------------------------------

    if (Kh > 0.0d0) then

       diffT_east = 0.0d0
       diffS_east = 0.0d0
       diffT_north = 0.0d0
       diffS_north = 0.0d0

       ! loop over all edges of locally owned plume cells
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
                                          +  diffT_north(i,j-1) - diffT_north(i,j)) * dt_plume/(dx*dy)
                work(i,j,3) = work(i,j,3) + (diffS_east(i-1,j)  - diffS_east(i,j)  &
                                          +  diffS_north(i,j-1) - diffS_north(i,j)) * dt_plume/(dx*dy)
             endif
          enddo
       enddo

    endif   ! Kh > 0

    !----------------------------------------------------------------------------
    ! Solve for D_plume, T_plume and S_plume
    !----------------------------------------------------------------------------

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

    !----------------------------------------------------------------------------
    ! Compute the mass divergence in each cell
    ! In steady state, this should balance the sum of entrainment and melting
    !----------------------------------------------------------------------------

    ! First find the flux at each east and north edge, using the upstream value of D_plume

    flux_east = 0.0d0
    flux_north = 0.0d0

    ! loop over all edges of locally owned cells
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          if (u_east(i,j) > 0.0d0) then
             flux_east(i,j) = D_plume(i,j)*u_east(i,j)
          else
             flux_east(i,j) = D_plume(i+1,j)*u_east(i,j)
          endif
          if (v_north(i,j) > 0.0d0) then
             flux_north(i,j) = D_plume(i,j)*v_north(i,j)
          else
             flux_north(i,j) = D_plume(i,j+1)*v_north(i,j)
          endif
       enddo
    enddo

    call parallel_halo(flux_east, parallel)
    call parallel_halo(flux_north, parallel)

    ! Sum the fluxes to compute the divergence
    divDu_plume = 0.0d0

    ! loop over locally owned cells
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          if (plume_mask(i,j) == 1) then
             divDu_plume(i,j) = (flux_east(i,j)  - flux_east(i-1,j)) / dx &
                              + (flux_north(i,j) - flux_north(i,j-1))/ dy
          endif
       enddo
    enddo


    ! final halo update
    call parallel_halo(D_plume, parallel)
    call parallel_halo(T_plume, parallel)
    call parallel_halo(S_plume, parallel)
    call parallel_halo(divDu_plume, parallel)

    if (verbose_plume) then
       call point_diag(D_plume, 'New D_plume (m)', itest, jtest, rtest, wx, wy)
       call point_diag(T_plume, 'T_plume (degC)', itest, jtest, rtest, wx, wy)
       call point_diag(S_plume, 'S_plume (psu)', itest, jtest, rtest, wx, wy)
       call point_diag(divDu_plume*scyr, 'divDu_plume (deg*m/yr)', itest, jtest, rtest, wx, wy)
    endif

    ! WHL - debug - Repeat for DuT
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          if (u_east(i,j) > 0.0d0) then
             flux_east(i,j) = D_plume(i,j)*u_east(i,j)*T_plume(i,j)
          else
             flux_east(i,j) = D_plume(i+1,j)*u_east(i,j)*T_plume(i+1,j)
          endif
          if (v_north(i,j) > 0.0d0) then
             flux_north(i,j) = D_plume(i,j)*v_north(i,j)*T_plume(i,j)
          else
             flux_north(i,j) = D_plume(i,j+1)*v_north(i,j)*T_plume(i,j+1)
          endif
       enddo
    enddo

    call parallel_halo(flux_east, parallel)
    call parallel_halo(flux_north, parallel)

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          if (plume_mask(i,j) == 1) then
             divDuT_plume(i,j) = (flux_east(i,j)  - flux_east(i-1,j)) / dx &
                               + (flux_north(i,j) - flux_north(i,j-1))/ dy
          endif
       enddo
    enddo
    if (verbose_plume) then
       i = itest
       j = jtest
       write(iulog,*) 'i, j, div(DuT) (deg*m/yr):', divDuT_plume(i,j)*scyr
       write(iulog,*) '   flux W, E:', flux_east(i-1,j)*scyr, flux_east(i,j)*scyr
       write(iulog,*) '   flux S, N:', flux_north(i,j-1)*scyr, flux_north(i,j)*scyr
    endif

    !----------------------------------------------------------------------------
    ! Make sure all output fields are in range
    ! If D_plume is out of range, the entrainment or detrainment should work
    !  to bring it back in range.
    ! If T_plume or S_plume is out of range, there may be a bug.
    !----------------------------------------------------------------------------

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
