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
    use cism_parallel, only: this_rank, main_task, nhalo, lhalo, uhalo, &
         parallel_halo, parallel_reduce_max, parallel_global_sum, &
         parallel_is_zero, parallel_globalindex

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
         ! Note: eos_Tref = -1 C and eos_Sref = 34.2 psu are not used
         f_coriolis = -1.405d-4        ! Coriolis parameter (s^-1) at 75 S = 2*omega*sin(75 deg) (prescribed in text)

    ! plume parameters
    !TODO - Add to the derived type?
    real(dp), parameter :: &
         D_plume0 = 5.d0,            & ! initial plume thickness (m)
         D_plume_min = 1.0d0,        & ! min plume thickness (m) where the plume exists
         D_plume_max = 50.0d0          ! max plume thickness (m)

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
            plume%T_plume,       plume%S_plume,        &
            plume%D_plume,                             &
            plume%T_basal,       plume%S_basal,        &
            plume%u_plume,       plume%v_plume,        &
            plume%u_plume_Cgrid, plume%v_plume_Cgrid,  &   ! is this needed?
            plume%ustar_plume,   plume%drho_plume,     &
            plume%entrainment,   plume%detrainment,    &
            plume%divDu_plume,                         &
            model%basal_melt%bmlt_float)

    endif  ! not a restart

    if (verbose_plume .and. main_task) write(iulog,*) 'Spun up the plume'

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
            plume%T_plume,       plume%S_plume,        &
            plume%D_plume,                             &
            plume%T_basal,       plume%S_basal,        &
            plume%u_plume,       plume%v_plume,        &
            plume%u_plume_Cgrid, plume%v_plume_Cgrid,  &   ! is this needed?
            plume%ustar_plume,   plume%drho_plume,     &
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
       T_plume,          S_plume,          &
       D_plume,                            &
       T_basal,          S_basal,          &
       u_plume,          v_plume,          &
       u_plume_Cgrid,    v_plume_Cgrid,    &
       ustar_plume,      drho_plume,       &
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

    !TODO - Add logic to stop when we reach total_time
    real(dp), intent(in) :: &
         dt_plume,            & ! plume timestep (s) for advection
         total_time             ! how long to run the plume model (s)

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

    real(dp), dimension(nx,ny), intent(inout) :: &
         T_plume,             & ! plume temperature (deg C)
         S_plume,             & ! plume salinity (psu)
         D_plume                ! plume thickness (m)

    ! Note: Plume velocities are computed at cell edges, and then are interpolated
    !       to cell centers as a diagnostic.
    !TODO - Is this a C grid or a CD grid?

    real(dp), dimension(nx,ny), intent(out) :: &
         u_plume,             & ! x component of plume velocity (m/s) at cell corners
         v_plume,             & ! y component of plume velocity (m/s) at cell corners
         u_plume_Cgrid,       & ! x component of plume velocity (m/s) on C grid (east edges)
         v_plume_Cgrid,       & ! y component of plume velocity (m/s) on C grid (north edges)
         ustar_plume,         & ! plume friction velocity (m/s) at cell centers
         drho_plume,          & ! density difference between plume and ambient ocean (kg/m^3)
         T_basal,             & ! basal ice temperature (deg C)
         S_basal,             & ! basal ice salinity (psu)
         entrainment,         & ! entrainment rate of ambient water into plume (m/s)
         detrainment,         & ! detrainment rate of plume into ambient water (m/s)
         divDu_plume,         & ! div(Du) for plume
         bmlt_float             ! melt rate at base of floating ice (m/s)

    ! Local variables

    integer, dimension(nx,ny) :: &
         plume_mask,          & ! = 1 for cells where scalar plume variables are computed
         ice_mask,            & ! = 1 if ice is present (thck > 0)
         floating_mask,       & ! = 1 where ice is present and floating, else = 0
         ocean_mask,          & ! = 1 if topg is below sea level and ice is absent, else = 0
         land_mask              ! = 1 if topg is at or above sea level, else = 0

    !TODO - Remove unused variables
    real(dp), dimension(nx,ny) :: &
         pressure,            & ! ocean pressure at base of ice (N/m^2)
         lsrf_plume,          & ! elevation of plume-ambient interface (m, negative below sea level)
         rho_plume,           & ! plume density (kg/m^3)
         rho_ambient,         & ! ambient ocean density (kg/m^3)
         H_cavity,            & ! thickness of ocean cavity beneath the plume (m)
         heat_transfer,       & ! rate of heat transfer from plume to ice (J/m2/s)
         dD_plume,            & ! change in D_plume (m)
         D_plume_old           ! D_plume from previous time step

    real(dp), dimension(nx,ny) ::  &
         u_plume_east,          & ! u_plume on east edges
         v_plume_east,          & ! v_plume on east edges
         u_plume_north,         & ! u_plume on north edges
         v_plume_north,         & ! v_plume on north edges
         plume_speed_east,      & ! plume speed on east edges (m/s)
         plume_speed_north        ! plume speed on north edges (m/s)

    ! Note: edge_mask = 0 for closed boundaries, = 2 for open boundaries (at least for now)
    integer, dimension(nx,ny) :: &
         edge_mask_east,        & ! = 1 on east edges where plume velocity is computed
         edge_mask_north          ! = 1 on north edges where plume velocity is computed

    real(dp), dimension(nx,ny) ::  &
         dlsrf_plume_dx_east,   & ! horizontal gradient of lsrf on east edges
         dlsrf_plume_dy_east,   & !
         dlsrf_plume_dx_north,  & ! horizontal gradient of lsrf on north edges
         dlsrf_plume_dy_north

    real(dp) :: &
         dlsrf_plume_dx,        & ! lsrf gradient components at cell centers
         dlsrf_plume_dy,        &
         slope                    ! magnitude of the gradient (dlsrf_dx, dlsrf_dy)

    real(dp), dimension(nx,ny) ::  &
         theta_slope            ! basal slope angle (rad), used for entrainment

    real(dp), dimension(nx-1,ny-1) :: &
         plume_speed            ! plume speed at vertices (m/s)

    real(dp) :: &
         time,                & ! elapsed time during the relaxation of the plume thickness (s)
         my_max_dt              ! CFL-limited time step for a given cell (s)

    real(dp) ::  &
         L2_norm,             & ! L2 norm of residual vector from continuity equation
         L2_previous            ! L2 norm from the previous convergence check

    integer :: i, j
    integer :: iglobal, jglobal        ! global i and j indices
    integer :: iter_Dplume             ! iteration counter

    ! parameters determining convergence of iterations
    !TODO - determine L2_target
    integer, parameter :: &
         n_check_convergence = 1,    & ! interval between convergence checks for D_plume
         L2_target = 0.0d0,          & ! convergence target for dD/dt
         maxiter_Dplume = 999999       ! max number of iterations of outer plume-thickness loop
                                       ! terminates when plume thickness reaches virtual steady state

    if (verbose_plume .and. main_task) then
       write(iulog,*) ' '
       write(iulog,*) 'In glissade_compute_plume'
    endif

    ! compute some masks

    call glissade_get_masks(&
         nx,                  ny,           &
         parallel,                          &
         thck,                topg,         &
         eus,                 0.0d0,        &  ! thklim = 0
         ice_mask,                          &
         floating_mask = floating_mask,     &
         land_mask = land_mask,             &
         ocean_mask = ocean_mask)

    call parallel_halo(floating_mask, parallel)

    ! Compute a mask that identifies where the plume is located
    !TODO - Refine this mask? Or cite Lambert 2026 as justification

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
    !TODO - Do we need pressure at base of plume?
    pressure = -rhoo*grav*lsrf

    ! Compute the cavity thickness
    H_cavity = max(lsrf - topg, 0.0d0)

    ! Set T_plume, S_plume and D_plume as needed
    ! On the first call, the input values are zero and these fields must be initialized everywhere.
    ! On subsequent calls, these fields are initialized only if the input values are zero.
    !TODO - Why not S_plume = S_ambient?
    !
    !    Earlier code has this comment:
    !      Set S_plume = S0 everywhere.
    !      This means that drho_plume = rho_ambient - rho_plume will decrease in the upslope direction.
    !      Setting both T_plume and S_plume to ambient values would give zero velocities and melt rates.

    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             if (T_plume(i,j) == 0.0d0) T_plume(i,j) = T_ambient(i,j)
             if (S_plume(i,j) == 0.0d0) S_plume(i,j) = S0
             if (D_plume(i,j) == 0.0d0) D_plume(i,j) = min(D_plume0, H_cavity(i,j))
          else   ! plume_mask = 0
             ! Zero out T, S and D
             T_plume(i,j) = 0.0d0
             S_plume(i,j) = 0.0d0
             D_plume(i,j) = 0.0d0
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

    ! Compute masks on cell edges, where plume velocities are computed.
    ! The mask = 1 if both adjacent cells have plume_mask_cell = 1.
    ! At closed boundaries (adjecent cell is grounded), set mask = 0.
    ! At open boundaries (adjecent cell is floating or open ocean), set mask = 2.
    !TODO - Free slip for flow parallel to edges?

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

    !WHL - commented out for now
!    do j = 1, ny
!       do i = 1, nx

!          call parallel_globalindex(i, j, iglobal, jglobal)

!          if (iglobal <= 0 .or. iglobal >= global_ewn .or. &  ! along or beyond EW boundary
!              jglobal <= 0 .or. jglobal >  global_nsn) then   ! beyond NS boundary
!             edge_mask_east(i,j) = 0
!          endif

!          if (jglobal <= 0 .or. jglobal >= global_nsn .or. &  ! along or beyond NS boundary
!              iglobal <= 0 .or. iglobal >  global_ewn) then   ! beyond EW boundary
!             edge_mask_north(i,j) = 0
!          endif

!       enddo
!    enddo

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

    !----------------------------------------------------------------
    ! Initialize some fields that are updated during the iteration.
    ! Note: T_plume, S_plume and D_plume are intent(inout) and already have initial values.
    !----------------------------------------------------------------

    !TODO - Not sure this is needed; do we solve for these?
    ! Initialize T and S at the base of the ice.
    ! Start with the same salinity as the underlying water, with T at the freezing point
    where (plume_mask == 1)
       S_basal = S_plume
       T_basal = lambda1*S_plume + lambda2 + lambda3*pressure
    elsewhere
       S_basal = S_ambient
       T_basal = lambda1*S_ambient + lambda2 + lambda3*pressure
    endwhere

    ! Initialize other fields
    u_plume = 0.0d0
    v_plume = 0.0d0
    plume_speed = 0.0d0

    u_plume_east = 0.0d0
    v_plume_east = 0.0d0
    u_plume_north = 0.0d0
    v_plume_north = 0.0d0
    plume_speed_east = 0.0d0
    plume_speed_north = 0.0d0

    ustar_plume = 0.0d0
    divDu_plume = 0.0d0
    entrainment = 0.0d0
    detrainment = 0.0d0
    bmlt_float = 0.0d0

    if (verbose_plume) then
       if (main_task) then
          write(iulog,*) ' '
          write(iulog,*) 'Initial fields:'
       endif
       call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(lsrf, 'lsrf (m)', itest, jtest, rtest, 7, 7)
       call point_diag(topg, 'topg (m)', itest, jtest, rtest, 7, 7)
       call point_diag(lsrf - topg, 'lsrf - topg (m)', itest, jtest, rtest, 7, 7)
       call point_diag(floating_mask, 'floating_mask', itest, jtest, rtest, 7, 7)
       call point_diag(plume_mask, 'plume_mask', itest, jtest, rtest, 7, 7)
       call point_diag(edge_mask_east,  'edge_mask_east', itest, jtest, rtest, 7, 7)
       call point_diag(edge_mask_north, 'edge_mask_north', itest, jtest, rtest, 7, 7)
       call point_diag(H_cavity, 'H_cavity', itest, jtest, rtest, 7, 7)
       call point_diag(T_ambient, 'T_ambient (deg C)', itest, jtest, rtest, 7, 7)
       call point_diag(S_ambient, 'S_ambient (psu)', itest, jtest, rtest, 7, 7)
       call point_diag(D_plume, 'D_plume (m)', itest, jtest, rtest, 7, 7)
       call point_diag(T_plume, 'T_plume (deg C)', itest, jtest, rtest, 7, 7)
       call point_diag(S_plume, 'S_plume (psu)', itest, jtest, rtest, 7, 7)
    endif

    time = 0.0d0

    !--------------------------------------------------------------------
    ! Iterate the plume to steady state. The solution method is:
    ! (1) Given the current ice geometry and D_plume, compute the plume velocity,
    !     entrainment, detrainment and melt rate.
    ! (2) Using the continuity equation, advance D_plume in time.
    ! (3) Repeat until the plume reaches a steady state.
    !--------------------------------------------------------------------

    ! initialize the L2 norm to an arbitrary big number
    L2_previous = huge(0.0d0)

    do iter_Dplume = 1, maxiter_Dplume   ! plume_thickness iteration

       ! advance the time (units of s)
       !TODO - Do we need to keep track of this, or just iter_Dplume?
       time = time + dt_plume

       if (verbose_plume) then
          if (main_task) write(iulog,*) 'iter_D_plume, time =', iter_Dplume, time
       endif

       ! Compute the plume density, given the current estimates of T_plume and S_plume.
       ! Then find the density difference between the ambient ocean and the plume.

       rho_plume = eos_rho_ref * (1.d0 - eos_alpha * (T_plume - eos_Tref)  &
                                       + eos_beta  * (S_plume - eos_Sref) )

       where (plume_mask == 1)
          drho_plume = rho_ambient - rho_plume
       elsewhere
          drho_plume = 0.0d0
       endwhere

       ! Compute the elevation of the lower plume boundary

       lsrf_plume = lsrf - D_plume

       ! Compute edge gradients of lsrf_plume
       ! These are used to compute the pressure gradient force

       dlsrf_plume_dx_east = 0.0d0
       dlsrf_plume_dy_east = 0.0d0
       dlsrf_plume_dx_north = 0.0d0
       dlsrf_plume_dy_north = 0.0d0

       ! Compute x gradients on east edges and y gradients on north edges

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (edge_mask_east(i,j) == 1) then
                dlsrf_plume_dx_east(i,j) = (lsrf_plume(i+1,j) - lsrf_plume(i,j)) / dx
             endif
             if (edge_mask_north(i,j) == 1) then
                dlsrf_plume_dy_north(i,j) = (lsrf_plume(i,j+1) - lsrf_plume(i,j)) / dy
             endif
          enddo
       enddo

       call parallel_halo(dlsrf_plume_dx_east, parallel)
       call parallel_halo(dlsrf_plume_dy_north, parallel)

       ! Interpolate to get y gradients on east edges and x gradients on north edges

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (edge_mask_east(i,j) == 1) then
                dlsrf_plume_dy_east(i,j) = 0.25d0 * &
                     (dlsrf_plume_dy_north(i,j)   + dlsrf_plume_dy_north(i+1,j)  &
                    + dlsrf_plume_dy_north(i,j-1) + dlsrf_plume_dy_north(i+1,j-1))
             endif
             if (edge_mask_north(i,j) == 1) then
                dlsrf_plume_dx_north(i,j) = 0.25d0 * &
                     (dlsrf_plume_dx_east(i-1,j+1) + dlsrf_plume_dx_east(i,j+1)  &
                    + dlsrf_plume_dx_east(i-1,j)   + dlsrf_plume_dx_east(i,j))
             endif
          enddo
       enddo

       call parallel_halo(dlsrf_plume_dy_east, parallel)
       call parallel_halo(dlsrf_plume_dx_north, parallel)

       !TODO - Check the indexing above. Extend to open boundaries

       ! Compute the slope angle at cell centers.  This is used to compute entrainment.
       !TODO - Compare to the method used in the nonlocal-slope scheme.

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

       ! Compute u_plume and v_plume at each edge
       ! Note: v_plume_north and u_plume_east are parallel to edges
       !       Computing both u and v at each edge leads to a more graceful treatment
       !        of the Coriolis terms than computing the perpendicular components alone.

       call compute_plume_velocity(&
            nx,           ny,      &
            dx,           dy,      &
            itest, jtest, rtest,   &
            parallel,              &
            plume_mask,            &
            edge_mask_east,        &
            edge_mask_north,       &
            dlsrf_plume_dx_east,   &
            dlsrf_plume_dy_east,   &
            dlsrf_plume_dx_north,  &
            dlsrf_plume_dy_north,  &
            drho_plume,            &
            D_plume,               &
            H_cavity,              &
            u_plume_east,          &
            v_plume_east,          &
            u_plume_north,         &
            v_plume_north,         &
            plume_speed_east,      &
            plume_speed_north)

       ! Compute the entrainment rate, given u_plume, v_plume and theta_slope

       call plume_entrainment(&
            nx,         ny,      &
            dx,         dy,      &
            itest, jtest, rtest, &
            plume_mask,          &
            theta_slope,         &
            u_plume_east,        &
            v_plume_north,       &
            entrainment)

       ! Compute the detrainment rate where D_plume exceeds D_plume_max

       call plume_detrainment(&
            nx,           ny,     &
            itest, jtest, rtest,  &
            H_cavity,             &
            D_plume,              &
            detrainment)

       ! Compute the basal melt rate, temperature and salinity at the plume-ice interface,
       ! given the plume properties.

       call plume_melt_rate(&
            nx,         ny,      &
            itest, jtest, rtest, &
            parallel,            &
            plume_mask,          &
            gammaT,              &
            gammaS,              &
            pressure,            &
            u_plume_east,        &
            v_plume_north,       &
            D_plume,             &
            T_plume,             &
            S_plume,             &
            ustar_plume,         &
            T_basal,             &
            S_basal,             &
            bmlt_float)

       ! halo updates
       call parallel_halo(T_plume, parallel)
       call parallel_halo(S_plume, parallel)

       if (verbose_plume) then
          if (main_task) write(iulog,*) 'Plume properties before advancing D_plume:'
          call point_diag(drho_plume, 'drho_plume (kg/m3)', itest, jtest, rtest, 7, 7)
          call point_diag(u_plume_east, 'u_plume_east (m/s)', itest, jtest, rtest, 7, 7)
          call point_diag(u_plume_north, 'u_plume_north (m/s)', itest, jtest, rtest, 7, 7)
          call point_diag(v_plume_east, 'v_plume_east (m/s)', itest, jtest, rtest, 7, 7)
          call point_diag(v_plume_north, 'v_plume_north (m/s)', itest, jtest, rtest, 7, 7)
          call point_diag(plume_speed, 'plume_speed (m/s)', itest, jtest, rtest, 7, 7)
          call point_diag(entrainment, 'entrainment (m/s)', itest, jtest, rtest, 7, 7)
          call point_diag(detrainment, 'detrainment (m/s)', itest, jtest, rtest, 7, 7)
          call point_diag(T_plume, 'T_plume (deg C)', itest, jtest, rtest, 7, 7)
          call point_diag(S_plume, 'S_plume (psu)', itest, jtest, rtest, 7, 7)
          call point_diag(T_basal, 'T_basal (deg C)', itest, jtest, rtest, 7, 7)
          call point_diag(S_basal, 'S_basal (psu)', itest, jtest, rtest, 7, 7)
          call point_diag(bmlt_float*scyr, 'bmlt_float (m/yr)', itest, jtest, rtest, 7, 7)
       endif

       ! Compute the rate of heat transfer (J/m^2/s) from the plume to the ice base
       where (plume_mask == 1)
          heat_transfer = rhow*cpw*ustar_plume*gammaT*(T_plume - T_basal)
       elsewhere
          heat_transfer = 0.0d0
       endwhere

       if (verbose_plume .and. main_task) then
          write(iulog,*) 'Advance the plume thickness, dt_plume, time (s) =', dt_plume, time
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

       ! Solve transport equations for D_plume, T_plume and S_plume,
       !  given u_plume, v_plume, entrainment, detrainment and bmlt_float.
       ! Note: Entrained water has ambient properties (T_ambient, S_ambient).
       !       Meltwater has basal properties (T basal, S basal).

       call plume_transport(&
            nx,           ny,     &
            dx,           dy,     &
            itest, jtest, rtest,  &
            parallel,             &
            dt_plume,             &
            plume_mask,           &
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

       if (verbose_plume) then
          call point_diag(D_plume_old, 'Old D_plume (m)', itest, jtest, rtest, 7, 7)
          call point_diag(D_plume, 'New D_plume (m)', itest, jtest, rtest, 7, 7)
          call point_diag((D_plume - D_plume_old)/dt_plume, 'dD/dt (m/s)', itest, jtest, rtest, 7, 7)
          call point_diag(divDu_plume, 'divergence (m/s)', itest, jtest, rtest, 7, 7)
       endif

       if (iter_Dplume >=2 .and. mod(iter_Dplume, n_check_convergence) == 0) then  ! check for convergence

          !TODO - Compute L2_norm based on dD/dt

          ! Check for convergence: dD/dt is small everywhere
          if (L2_norm < L2_target) then
             if (verbose_plume .and. main_task) then
                write(iulog,*) 'Continuity converged, time, iter, L2_norm =', time, iter_Dplume, L2_norm
             endif
          elseif (L2_norm < L2_previous) then ! iteration is converging; keep going
             if (verbose_plume .and. main_task) then
                write(iulog,*) 'Continuty not yet converged, time, iter, L2_norm =', time, iter_Dplume, L2_norm
             endif
          elseif (L2_norm >= L2_previous) then ! iteration is not converging
             if (verbose_plume .and. main_task) then
                write(iulog,*) 'Continuty not converging, time, iter, L2_norm =', time, iter_Dplume, L2_norm
             endif
          endif

          ! save variables from this iteration
          D_plume_old = D_plume
          L2_previous = L2_norm

       endif   ! mod(iter_Dplume, n_check_convergence) = 0

    enddo   ! iter_Dplume

    if (verbose_plume) then
       if (main_task) then
          write(iulog,*) ' '
          write(iulog,*) 'Final plume properties:'
       endif
       call point_diag(D_plume, 'D_plume (m)', itest, jtest, rtest, 7, 7)
       call point_diag(T_plume, 'T_plume (m)', itest, jtest, rtest, 7, 7)
       call point_diag(S_plume, 'S_plume (m)', itest, jtest, rtest, 7, 7)
       call point_diag(T_basal, 'T_basal (m)', itest, jtest, rtest, 7, 7)
       call point_diag(S_basal, 'S_basal (m)', itest, jtest, rtest, 7, 7)
       call point_diag(u_plume, 'u_plume (m)', itest, jtest, rtest, 7, 7)
       call point_diag(v_plume, 'v_plume (m)', itest, jtest, rtest, 7, 7)
       call point_diag(entrainment, 'entrainment (m/s)', itest, jtest, rtest, 7, 7)
       call point_diag(detrainment, 'detrainment (m/s)', itest, jtest, rtest, 7, 7)
       call point_diag(bmlt_float*scyr, 'bmlt_float (m/yr)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine compute_plume

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
       dlsrf_plume_dx_east,    &
       dlsrf_plume_dy_east,    &
       dlsrf_plume_dx_north,   &
       dlsrf_plume_dy_north,   &
       drho_plume,             &
       D_plume,                &
       H_cavity,               &
       u_plume_east,           &
       v_plume_east,           &
       u_plume_north,          &
       v_plume_north,          &
       plume_speed_east,       &
       plume_speed_north)

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
         dlsrf_plume_dx_east,   & ! horizontal gradient of lsrf_plume on east edges
         dlsrf_plume_dy_east,   & !
         dlsrf_plume_dx_north,  & ! horizontal gradient of lsrf_plume on north edges
         dlsrf_plume_dy_north,  & !
         drho_plume,            & ! density difference between plume and ambient ocean (kg/m^3)
         D_plume,               & ! plume thickness (m)
         H_cavity                 ! thickness of ocean cavity beneath the plume (m)

!    real(dp), dimension(nx,ny), intent(in) ::  &
!         edge_mask_east_reduce_v,  & ! mask for reducing v on east edges adjacent to a wall
!         edge_mask_north_reduce_u    ! mask for reducing u on north edges adjacent to a wall


    real(dp), dimension(nx,ny), intent(inout) ::  &
         plume_speed_east,    & ! plume speed on east edges
         plume_speed_north      ! plume speed on north edges

    real(dp), dimension(nx,ny), intent(out) ::  &
         u_plume_east,        & ! u_plume on east edges
         v_plume_east,        & ! v_plume on east edges
         u_plume_north,       & ! u_plume on north edges
         v_plume_north          ! v_plume on north edges

    ! local variables

    real(dp), dimension(nx,ny) :: &
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
         D_plume_north,       & ! D_plume averaged to north edge
         grav_reduced_east,   & ! reduced gravity on east edge
         grav_reduced_north     ! reduced gravity on north edge

    integer :: i, j

    integer :: iter_velo        ! iteration counter

    character(len=100) :: message

    logical, dimension(nx,ny) ::  &
         converged_velo_east, & ! true when velocity has converged at an east edge, else false
         converged_velo_north   ! true when velocity has converged at a north edge, else false

    logical :: &
         converged_all_velo     ! true when velocity has converged at all edges, else false

    integer, parameter ::  &
         maxiter_velo = 100     ! max number of iterations of velocity loop

    ! initialize

    u_plume_east = 0.0d0
    v_plume_east = 0.0d0

    u_plume_north = 0.0d0
    v_plume_north = 0.0d0

    D_plume_east = 0.0d0
    D_plume_north = 0.0d0

    grav_reduced_east = 0.0d0
    grav_reduced_north = 0.0d0

    pgf_x_east = 0.0d0
    pgf_y_east = 0.0d0
    pgf_x_north = 0.0d0
    pgf_y_north = 0.0d0

    !TODO - Use method (2)?
    ! Note: There are a couple of different ways to compute the PGF.
    !       (1) Jenkins et al. (1991) and HJH (2008) use grad(lsrf)
    !       (2) Holland & Feltham (2006) use grad(lsrf_plume) along with a density gradient.
    !       Method (1) is simpler and has the advantage that grad(lsrf) does not vary during plume evolution,
    !        making the PGF more stable (though possibly not as accurate).

    ! Compute the pressure gradient force on each edge
    !TODO - Add the terms proportional to d/dx and d/dy(drho_plume)

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! PGF on east edge
          if (edge_mask_east(i,j) == 1) then
             D_plume_east(i,j) = (D_plume(i,j) + D_plume(i+1,j)) / 2.0d0
             grav_reduced_east(i,j) = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i+1,j)) / 2.0d0
             pgf_x_east(i,j) = grav_reduced_east(i,j) * D_plume_east(i,j) * dlsrf_plume_dx_east(i,j)
             pgf_y_east(i,j) = grav_reduced_east(i,j) * D_plume_east(i,j) * dlsrf_plume_dy_east(i,j)
          endif

          ! PGF on north edge
          if (edge_mask_north(i,j) == 1) then
             D_plume_north(i,j) = (D_plume(i,j) + D_plume(i,j+1)) / 2.0d0
             grav_reduced_north(i,j) = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i,j+1)) / 2.0d0
             pgf_x_north(i,j) = grav_reduced_north(i,j) * D_plume_north(i,j) * dlsrf_plume_dx_north(i,j)
             pgf_y_north(i,j) = grav_reduced_north(i,j) * D_plume_north(i,j) * dlsrf_plume_dy_north(i,j)
          endif   ! edge_mask_north

       enddo  ! i
    enddo  ! j

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

       if (verbose_plume .and. main_task) then
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

       converged_all_velo = .true.

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (edge_mask_east(i,j) == 1 .and. .not.converged_velo_east(i,j) ) then
                converged_all_velo = .false.
             endif
             if (edge_mask_north(i,j) == 1 .and. .not.converged_velo_north(i,j) ) then
                converged_all_velo = .false.
             endif
          enddo
       enddo

       if (converged_all_velo) then
          if (verbose_plume .and. main_task) write(iulog,*) 'Plume velocity converged'
          exit   ! iter_velo loop
       elseif (iter_velo == maxiter_velo) then
          write(message,*) 'Error, glissade_plume: velocity has not converged, iter_velo =', iter_velo
          call write_log(message, GM_FATAL)
       endif

    enddo  ! iter_velo

    ! Extrapolate velocity to open boundaries (edge_mask = 2)

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

       enddo
    enddo

    call parallel_halo(u_plume_east, parallel)
    call parallel_halo(v_plume_east, parallel)
    call parallel_halo(u_plume_north, parallel)
    call parallel_halo(v_plume_north, parallel)

    if (verbose_plume) then
       call point_diag(pgf_x_east, 'pgf_x_east', itest, jtest, rtest, 7, 7)
       call point_diag(pgf_y_east, 'pgf_y_east', itest, jtest, rtest, 7, 7)
       call point_diag(pgf_x_north, 'pgf_x_north', itest, jtest, rtest, 7, 7)
       call point_diag(pgf_y_north, 'pgf_y_north', itest, jtest, rtest, 7, 7)
       call point_diag(u_plume_east, 'u_plume_east', itest, jtest, rtest, 7, 7)
       call point_diag(v_plume_east, 'v_plume_east', itest, jtest, rtest, 7, 7)
       call point_diag(u_plume_north, 'u_plume_north', itest, jtest, rtest, 7, 7)
       call point_diag(v_plume_north, 'v_plume_north', itest, jtest, rtest, 7, 7)
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

    ! Update the plume speed at the edges (including the u_tidal term)
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          plume_speed_east(i,j)  = sqrt(u_plume_east(i,j)**2 + v_plume_east(i,j)**2 + u_tidal**2)
          plume_speed_north(i,j) = sqrt(u_plume_north(i,j)**2 + v_plume_north(i,j)**2 + u_tidal**2)
       enddo   ! i
    enddo   ! j

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
    
    ! Used to be intent(in), but now are module variables
!    real(dp), intent(in) ::   &
!         u_tidal,           & ! tidal velocity (m/s)
!         c_drag,            & ! ocean drag coefficient (unitless)   
!         f_coriolis           ! Coriolis parameter (s^-1)
    
    integer, dimension(nx,ny), intent(in) ::   &
         edge_mask            ! = 1 at edges where velocity is computed

    ! Note: The following variables are co-located with the velocity
    real(dp), dimension(nx,ny), intent(in) ::   &
         D_plume,           & ! plume thickness at edges (m)
         pgf_x,             & ! x component of pressure gradient force
         pgf_y                ! y component of pressure gradient force
!         latdrag_x,         & ! x component of lateral drag
!         latdrag_y            ! y component of lateral drag
    

    real(dp), dimension(nx,ny), intent(inout) ::  &
         u_plume,           & ! x component of plume velocity (m/s)
         v_plume              ! x component of plume velocity (m/s)

!!    logical, dimension(nx,ny), intent(inout) ::  &
    logical, dimension(nx,ny), intent(out) ::  &
         converged_velo        ! true when velocity has converged at an edge, else false

    ! local variables

!    real(dp), dimension(nx,ny) ::   &
!         f_x,               &  ! pgf_x + latdrag_x
!         f_y                   ! pgf_y + latdrag_y

!    real(dp), dimension(nx,ny) ::  &
!         reduce_v,          &  ! local version of edge_mask_east_reduce_v; no reduction by default
!         reduce_u              ! local version of edge_mask_north_reduce_u; no reduction by default

    real(dp) :: &
         plume_speed,       & ! plume speed (m/s) based on input (u_plume, v_plume)
         x_resid, y_resid,  & ! residuals of momentum balance equations (m^2/s^2)
         denom,             & ! denominator
         a_uu, a_uv,        & ! coefficients for Newton solve
         a_vu, a_vv,        & !
         du, dv               ! change in u_plume and v_plume (m/s)
    
    character(len=128) :: message

    real(dp), parameter :: &
         maxresid_force_balance = 1.0d-8 ! max residual allowed in momentum balance equation (m^2/s^2)

    !TODO - Start with Picard, then test Newton
    logical, parameter :: &
!         velo_newton = .true.  ! if true, use Newton's method; if false, use Picard method
         velo_newton = .false.  ! if true, use Newton's method; if false, use Picard method

    integer :: i, j

    !--------------------------------------------------------------------
    ! Compute the plume velocity.
    ! Assume a balance between the pressure gradient force, basal drag and Coriolis:
    !
    ! pgf_x - c_d*|U|*u + D*f*v = 0
    ! pgf_y - c_d*|U|*v - D*f*u = 0
    !
    !  where pgf_x = g' * D * db/dx (m^2/s^2) 
    !        pgf_y = g' * D * db/dy (m^2/s^2) 
    !            D = plume boundary-layer thickness
    !           g' = reduced gravity = g*(rhoa - rhop)/rhoo
    !         rhoa = ambient ocean density
    !         rhop = plume density
    !         rhoo = reference ocean density
    !            b = elevation of shelf base   !TODO - plume base?
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
    ! If reduce_u < 1 or reduce_v < 1, then the Coriolis term in these equations
    ! is reduced proportionately, so as to inhibit flow into walls.
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

          ! Compute plume speed based on current u and v
          plume_speed = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)

!!          if (edge_mask(i,j) == 1 .and. .not.converged_velo(i,j) ) then
          if (edge_mask(i,j) == 1) then
       
             ! Compute residual of the momentum balance
             x_resid = pgf_x(i,j) - c_drag*plume_speed*u_plume(i,j) + f_coriolis*D_plume(i,j)*v_plume(i,j)
             y_resid = pgf_y(i,j) - c_drag*plume_speed*v_plume(i,j) - f_coriolis*D_plume(i,j)*u_plume(i,j)
!             x_resid = f_x(i,j) - c_drag*plume_speed*u_plume(i,j) + reduce_v(i,j)*f_coriolis*D_plume(i,j)*v_plume(i,j)
!             y_resid = f_y(i,j) - c_drag*plume_speed*v_plume(i,j) - reduce_u(i,j)*f_coriolis*D_plume(i,j)*u_plume(i,j)

             ! check convergence of plume velocity
             if (abs(x_resid) < maxresid_force_balance .and. abs(y_resid) < maxresid_force_balance) then
                converged_velo(i,j) = .true.
             endif

             ! Should be no harm to compute anyway, even if already converged. (Verify this)
!!             if (.not.converged_velo(i,j)) then

             if (velo_newton) then
          
                ! compute some coefficients for the Newton solve
                a_uu = c_drag * (plume_speed + u_plume(i,j)**2/plume_speed)
                a_vv = c_drag * (plume_speed + v_plume(i,j)**2/plume_speed)
                      
                a_uv = c_drag * (u_plume(i,j)*v_plume(i,j))/plume_speed - D_plume(i,j)*f_coriolis
                a_vu = c_drag * (u_plume(i,j)*v_plume(i,j))/plume_speed + D_plume(i,j)*f_coriolis
!                   a_uv = c_drag * (u_plume(i,j)*v_plume(i,j))/plume_speed - reduce_v(i,j)*D_plume(i,j)*f_coriolis
!                   a_vu = c_drag * (u_plume(i,j)*v_plume(i,j))/plume_speed + reduce_u(i,j)*D_plume(i,j)*f_coriolis
                   
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
          
                denom = (c_drag*plume_speed)**2 + (D_plume(i,j)*f_coriolis)**2
!                   u_plume = (c_drag*plume_speed*f_x(i,j) + reduce_v(i,j)*D_plume(i,j)*f_coriolis*f_y(i,j)) / denom
!                   v_plume = (c_drag*plume_speed*f_y(i,j) - reduce_u(i,j)*D_plume(i,j)*f_coriolis*f_x(i,j)) / denom
                u_plume(i,j) = (c_drag*plume_speed*pgf_x(i,j) + f_coriolis*D_plume(i,j)*pgf_y(i,j)) / denom
                v_plume(i,j) = (c_drag*plume_speed*pgf_y(i,j) - f_coriolis*D_plume(i,j)*pgf_x(i,j)) / denom
          
             endif  ! Newton or Picard

!!             endif  ! .not.converged_velo

             if (verbose_plume .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'plume_speed (m/s) =', plume_speed
                write(iulog,*) 'pgf_x, pgf_y:', pgf_x(i,j), pgf_y(i,j)
!                write(iulog,*) 'latdrag_x, latdrag_y:', latdrag_x(i,j), latdrag_y(i,j)
                write(iulog,*) 'Dfv, -Dfu:', D_plume(i,j) * f_coriolis * v_plume(i,j), &
                                     -D_plume(i,j) * f_coriolis * u_plume(i,j)
                write(iulog,*) 'dragu, dragv:', c_drag * plume_speed * u_plume(i,j), &
                                         c_drag * plume_speed * v_plume(i,j)
                write(iulog,*) 'x/y residual:', x_resid, y_resid
                write(iulog,*) 'new u/v_plume:', u_plume(i,j), v_plume(i,j)
                write(iulog,*) 'converged =', converged_velo(i,j)
             endif

          endif  ! edge_mask
       enddo  ! i
    enddo  ! j

  end subroutine plume_velocity

!****************************************************

  subroutine plume_entrainment(&
       nx,         ny,      &
       dx,         dy,      &
       itest, jtest, rtest, &
       plume_mask,          &
       theta_slope,         &
       u_plume_east,        &
       v_plume_north,       &
       entrainment)

    !--------------------------------------------------------------------
    ! Compute entrainment as a function of the plume speed and the slope of the
    !  plume-ambient interface, following Bo Pederson (1980) and Jenkins (1991).
    !--------------------------------------------------------------------

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask             ! = 1 for cells where scalar plume variables are computed

    !TODO - Also pass in u_plume_north and v_plume_east?
    real(dp), dimension(nx,ny), intent(in) ::  &
         u_plume_east,          & ! u component of plume velocity on east edges (m/s)
         v_plume_north,         & ! v component of plume velocity on north edges (m/s)
         theta_slope            ! basal slope angle at cell centers (rad)

    real(dp), dimension(nx,ny), intent(out) ::  &
         entrainment              ! entrainment at cell centers (m/s)

    ! local variables

    real(dp) :: &
         u_plume_cell,           & ! u_plume averaged to cell center (m/s)
         v_plume_cell,           & ! v_plume averaged to cell center (m/s)
         plume_speed_cell          ! plume speed at cell center (m/s)

    integer :: i, j

    ! entrainment parameters
    real(dp), parameter ::   &
!!         H0_cavity = 10.d0,          & ! cavity thickness (m) below which the entrainment gradually approaches zero
         E0 = 0.072d0                  ! entrainment coefficient (unitless)
                                       ! Bo Pederson (1980) suggests E0 = 0.072
                                       ! Jenkins (1991, JGR) uses 0.036 to compensate for lack of Coriolis in 1D model

    entrainment = 0.0d0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             u_plume_cell = 0.5d0 * (u_plume_east(i-1,j) + u_plume_east(i,j))
             v_plume_cell = 0.5d0 * (v_plume_north(i,j-1) + v_plume_north(i,j))
             plume_speed_cell = sqrt(u_plume_cell**2 + v_plume_cell**2)
             entrainment(i,j) = E0 * plume_speed_cell * sin(theta_slope(i,j))
          endif
       enddo
    enddo

  end subroutine plume_entrainment

!****************************************************

  subroutine plume_detrainment(&
       nx,       ny,   &
       itest, jtest, rtest, &
       H_cavity,       &
       D_plume,        &
       detrainment)

    ! Compute detrainment.
    ! This is not a physically based mechanism, just a regularization to prevent very thick plumes.
    ! Ideally, detrainment = 0 almost everywhere.

    integer, intent(in) ::  &
         nx,     ny           ! number of grid cells in each dimension

    integer, intent(in) ::  &
         itest, jtest, rtest  ! test cell coordinates (diagnostic only)

    real(dp), dimension(nx,ny), intent(in) ::  &
         H_cavity,          & ! cavity thickness (m), lsrf - topg
         D_plume              ! plume thickess (m)

    real(dp), dimension(nx,ny), intent(out) ::  &
         detrainment          ! plume detrainment rate (m/s)

    ! local variables

    integer :: i, j

    ! detrainment parameters
    real(dp), parameter ::  &
         tau_detrainment = 3600.d0      ! detrainment time scale (s)

    detrainment = 0.0d0

    do j = 1, ny
       do i = 1, nx
          !TODO - Strictly limit D_plume <= H_cavity?
          if (D_plume(i,j) > H_cavity(i,j)) then
             detrainment(i,j) = (D_plume(i,j) - H_cavity(i,j)) / tau_detrainment
          elseif (D_plume(i,j) > D_plume_max) then
             detrainment(i,j) = (D_plume(i,j) - D_plume_max) / tau_detrainment
          endif
       enddo
    enddo

  end subroutine plume_detrainment

!****************************************************

  subroutine plume_melt_rate(&
       nx,         ny,      &
       itest, jtest, rtest, &
       parallel,            &
       plume_mask,          &
       gammaT,              &
       gammaS,              &
       pressure,            &
       u_plume_east,        &
       v_plume_north,       &
       D_plume,             &
       T_plume,             &
       S_plume,             &
       ustar_plume,         &
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
    ! (1) rhoi * m * L  = rhoo * cpw * u_fric * gammaT * (Tp - Tb)
    ! (2) rhoi * m * Sb = rhoo * u_fric * gammaS *(Sp - Sb)
    ! (3) Tb = lambda1*Sb + lambda2 + lambda3*pb 
    !
    ! Eqs. 1 and 2 describe heat and salt transfer at the ice-ocean interface.
    ! Eq. 3 is the linearized liquidus relation that determines the potential freezing point.
    ! Note: Asay-Davis et al. use rhow instead of rhoo on the LHS, since they define
    !       the melt rate m in units of meters of freshwater instead of meters of ice.
    !       See their Sec. 3.1.8.
    !
    ! We can rewrite these equations as
    !
    ! (1)     m = C1 * (Tp - Tb)
    ! (2)  m*Sb = C2 * (Sp - Sb)
    ! (3)    Tb = lambda1*Sb + C3
    !
    ! where C1 = (rhoo * cpw * ufric * gammaT) / (rhoi * L)
    !       C2 = (rhoo * ufric * gammaS) / rhoi
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
    ! Note: This treatment assumes that gammaT and gammaS are spatially uniform constants.
    !       Lambert et al. (2023) have the following instead:
    !       (1) m * L = cpw * gammaT * (Tp - Tb)
    !       (2) m * Sb = gammaS * (Sp - Sb)
    !       where gammaT = ustar_plume / 2.12d0*log(ustar_plume*D_plume/kvw) + 12.5d0*Prandtl**(2.0d0/3.0d0) - 8.68d0
    !             gammaS = ustar_plume / 2.12d0*log(ustar_plume*D_plume/kvw) + 12.5d0*Schmidt**(2.0d0/3.0d0) - 8.68d0
    !             kvw = kinematic viscosity of seawater
    !             Prandtl and Schmidt are dimensionless numbers for turbulent transfer
    !--------------------------------------------------------------------

    ! input/output variables
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
         u_plume_east,        & ! u_plume on east edges (m/s)
         v_plume_north,       & ! v_plume on north edges (m/s)
         D_plume,             & ! plume thickness (m)
         T_plume,             & ! plume temperature (deg C)
         S_plume                ! plume salinity (psu)

    real(dp), dimension(nx,ny), intent(out) :: &
         ustar_plume,         & ! plume friction velocity (m/s) on ice grid, output as a diagnostic
         T_basal,             & ! basal ice temperature (deg C)
         S_basal,             & ! basal ice salinity (psu)
         bmlt_float             ! melt rate at base of floating ice (m/s)
    
    ! local variables
    
    real(dp) :: &
         u_plume, v_plume,    & ! plume velocity components at cell center (m/s)
         C1, C2, C3,          & ! factors in melt-rate equations
         aa, bb, cc,          & ! factors in quadratic formula
         discriminant,        & ! (b^2 - 4ac) term in quadratic formula
         Sb1, Sb2               ! solutions of quadratic formula
    
    integer :: i, j, ig, jg

    logical :: abort            ! if true, then abort

    ! initialize
    ustar_plume = 0.0d0
    T_basal = 0.0d0
    S_basal = 0.0d0
    bmlt_float = 0.0d0

    ! Loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          
          if (plume_mask(i,j) == 1) then

             ! Interpolate the plume speed to the cell center, and compute the friction velocity ustar.
             u_plume = (u_plume_east(i-1,j) + u_plume_east(i,j)) / 2.0d0
             v_plume = (v_plume_north(i,j-1) + v_plume_north(i,j)) / 2.0d0
             ustar_plume(i,j) = sqrt(c_drag*(u_plume**2 + v_plume**2 + u_tidal**2))

             ! Solve a quadratic equation for S_basal
             C1 = (rhoo * cpw * ustar_plume(i,j) * gammaT) / (rhoi * lhci)
             C2 = (rhoo * ustar_plume(i,j) * gammaS) / rhoi
             C3 = lambda2 + lambda3*pressure(i,j)

             aa = lambda1*C1
             bb = C1*(C3 - T_plume(i,j)) - C2
             cc = C2*S_plume(i,j)

             abort = .false.
             discriminant = bb**2 - 4.d0*aa*cc
             if (discriminant >= 0.0d0) then
                Sb1 = (-bb + sqrt(discriminant)) / (2.0d0*aa)
                Sb2 = (-bb - sqrt(discriminant)) / (2.0d0*aa)
                if (Sb1 >= 0.0d0 .and. Sb2 <= 0.0d0) then
                   S_basal(i,j) = Sb1
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
                call write_log('Failed to solve quadratic equation for S_plume', GM_FATAL)
             endif

             ! Solve for T_basal and bmlt_float
             T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
             bmlt_float(i,j) = C1 * (T_plume(i,j) - lambda1*S_basal(i,j) - C3)

             if (verbose_plume .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'Melt rate calc: rank, i, j =', rtest, i, j
                write(iulog,*) 'pressure (Pa) =', pressure(i,j)
                write(iulog,*) 'C1 (m/s/deg), C2 (m/s), C3(deg C):', C1, C2, C3
                write(iulog,*) 'aa, bb, cc:=', aa, bb, cc
                write(iulog,*) 'T_basal, S_basal, bmlt_float:', T_basal(i,j), S_basal(i,j), bmlt_float(i,j)
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

    ! Solve transport equations for the plume thickness, temperature and salinity.
    ! Includes upwind-weighted horizontal transport as well as local entrainment,
    !  detrainment and melting.
    ! TODO: Replace upwind with remapping transport?

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
         plume_mask             ! = 1 for cells where the plume is present, else = 0

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

    integer :: i, j, ig, jg
    integer :: ilo, ihi, jlo, jhi

    real(dp) :: dD, dDT, dDS    ! increments in D, D*T and d*S
    real(dp), dimension(nx,ny,3) :: work   ! work array for transport

    character(len=100) :: message

    real(dp), parameter :: &
         T_plume_min = -3.0d0,       & ! min allowed T_plume (deg C)
         T_plume_max = 10.0d0,       & ! max allowed T_plume (deg C)
         S_plume_min = 0.0d0,        & ! min allowed S_plume (psu)
         S_plume_max = 40.0d0          ! max allowed S_plume (psu)

    ! Compute local column adjustments from entrainment, detrainment and melting

    ! Make sure all input fields are in range
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             ! Note: Allow D_plume > D_plume_max; detrainment should relax toward D_plume_max
             if (D_plume(i,j) < D_plume_min .or. D_plume(i,j) < 2.0d0*D_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(message,*) 'Plume transport, input D_plume out of range: ig, jg, D_plume =', &
                     ig, jg, D_plume(i,j)
                call write_log(message, GM_FATAL)
             endif
             ! Note: Det
             if (T_plume(i,j) < T_plume_min .or. T_plume(i,j) > T_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(message,*) 'Plume transport, input T_plume out of range: ig, jg, T_plume =', &
                     ig, jg, T_plume(i,j)
                call write_log(message, GM_FATAL)
             endif
             if (S_plume(i,j) < S_plume_min .or. S_plume(i,j) > S_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(message,*) 'Plume transport, input S_plume out of range: ig, jg, S_plume =', &
                     ig, jg, S_plume(i,j)
                call write_log(message, GM_FATAL)
             endif
          endif
       enddo
    enddo

    ! Fill a work array with fields to be incremented
    work(:,:,1) = D_plume
    work(:,:,2) = D_plume*T_plume
    work(:,:,3) = D_plume*S_plume

    ! Increment the work array
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             dD = entrainment(i,j) - detrainment(i,j) + bmlt_float(i,j)
             work(i,j,1) = work(i,j,1) + dD*dt
             ! Make sure the adjusted D_plume >= D_plume_min. If not, then the detrainment is excessive.
             if (work(i,j,1) < D_plume_min) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(message,*) 'Error, adjusted D_plume < D_plume_min: ig, jg, D_plume =', &
                     ig, jg, work(i,j,1)
                call write_log(message, GM_FATAL)
             endif
             ! Note: heat_transfer has units J/m^2/s, so heat_transfer/(rhow*cpw) has units of m*deg/s, as desired
             dDT = entrainment(i,j)*T_ambient(i,j) - detrainment(i,j)*T_plume(i,j) + bmlt_float(i,j)*T_basal(i,j) &
                  - heat_transfer(i,j)/(rhow*cpw)
             work(i,j,2) = work(i,j,2) + dDT*dt
             ! Note: salt_transfer = 0 by assumption
             dDS = entrainment(i,j)*S_ambient(i,j) - detrainment(i,j)*S_plume(i,j) + bmlt_float(i,j)*S_basal(i,j)
             work(i,j,3) = work(i,j,3) + dDS*dt
          endif
       enddo
    enddo

    call parallel_halo(work, parallel)

    ! Set bounds for loops over locally owned cells and edges
    ilo = nhalo + 1
    ihi = nx - nhalo
    jlo = nhalo + 1
    jhi = ny - nhalo

    ! Use a first-order accurate upwind scheme to transport D_plume

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

    ! Back out D_plume, T_plume and S_plume

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

    call parallel_halo(D_plume, parallel)
    call parallel_halo(T_plume, parallel)
    call parallel_halo(S_plume, parallel)

    ! Require D_plume >= D_plume_min.
    ! If an adjustment is needed, then keep T_plume and S_plume unchanged for simplicity.

    where (plume_mask == 1)
       D_plume = max(D_plume, D_plume_min)
    endwhere

    ! Make sure all output fields are in range
    ! If not, then (assuming the input fields were in range) there is a bug in the transport.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask(i,j) == 1) then
             ! Note: Allow D_plume > D_plume_max; detrainment should relax toward D_plume_max
             if (D_plume(i,j) < D_plume_min .or. D_plume(i,j) < 2.0d0*D_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(message,*) 'Plume transport, input D_plume out of range: ig, jg, D_plume =', &
                     ig, jg, D_plume(i,j)
                call write_log(message, GM_FATAL)
             endif
             ! Note: Det
             if (T_plume(i,j) < T_plume_min .or. T_plume(i,j) > T_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(message,*) 'Plume transport, input T_plume out of range: ig, jg, T_plume =', &
                     ig, jg, T_plume(i,j)
                call write_log(message, GM_FATAL)
             endif
             if (S_plume(i,j) < S_plume_min .or. S_plume(i,j) > S_plume_max) then
                call parallel_globalindex(i, j, ig, jg, parallel)
                write(message,*) 'Plume transport, input S_plume out of range: ig, jg, S_plume =', &
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
