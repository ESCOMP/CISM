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
    use glimmer_physcon, only: rhoi, rhow, rhoo, grav, lhci, cpw, pi, scyr
    use glimmer_paramets, only: iulog, eps11
    use glimmer_log
    use glimmer_utils, only: point_diag
    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo, lhalo, uhalo, &
         parallel_halo, parallel_reduce_max, parallel_global_sum, parallel_globalindex

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

  subroutine glissade_plume_init(model, plume)

    ! Initialize the plume properties

    ! input/ouput arguments

    type(glide_global_type), intent(inout) :: model   !> derived type holding ice-sheet info
    type(glide_plume), intent(inout) :: plume    !> derived type holding plume info

    ! local variables

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

       if (plume%misomip_domain) then

          ! MISOMIP+ profiles, Eqs. 21 and 22
          plume%T_ambient = plume%T0 + (plume%Tbot - plume%T0) * (model%geometry%lsrf / plume%zbed_deep)
          plume%S_ambient = plume%S0 + (plume%Sbot - plume%S0) * (model%geometry%lsrf / plume%zbed_deep)

       else
          !TODO - Work out how to initialize T_ambient and S_ambient
       endif   ! misomip_domain

       ! Spin up the plume to steady state
       ! Note: Typically, the initial spin-up takes longer than the runtime update.
       !       For an ISOMIP+ experiment with a fixed ice cavity, this is all we need to do.

       call compute_plume(&
            ewn,                 nsn,                &
            dew,                 dns,                &
            itest,   jtest,      rtest,              &
            parallel,                                &
            model%geometry%thck,                     &
            model%geometry%lsrf,                     &
            model%geometry%topg,                     &
            model%climate%eus,                       &
            plume%T_ambient,     plume%S_ambient,    &
            plume%gammaT,        plume%gammaS,       &
            plume%S0,                                &
            plume%T_plume,       plume%S_plume,      &
            plume%D_plume,                           &
            plume%T_basal,       plume%S_basal,      &
            plume%u_plume,       plume%v_plume,      &
            plume%u_plume_Cgrid, plume%v_plume_Cgrid,&
            plume%ustar_plume,                       &
            plume%drho_plume,                        &
            plume%entrainment,   plume%detrainment,  &
            plume%divDu_plume,                       &
            model%basal_melt%bmlt_float)

    endif  ! not a restart

    if (verbose_plume .and. main_task) write(iulog,*) 'Spun up the plume'

  end subroutine glissade_plume_init

!****************************************************

  subroutine glissade_plume_driver(model, plume)

    ! Compute melt rates using a plume model, given vertical profiles of T and S in the ambient ocean
    !
    ! The benchmark application is to the MISOMIP domain, described here:
    ! See this paper for details:
    ! X. S. Asay-Davis et al. (2016), Experimental design for three interrelated
    !    marine ice sheet and ocean model intercomparison projects:
    !    MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    !    Geosci. Model Devel., 9, 2471-2497, doi: 10.5194/gmd-9-2471-2016.

    type(glide_global_type), intent(inout) :: model   !> derived type holding ice-sheet info
    type(glide_plume), intent(inout) :: plume    !> derived type holding plume info

    ! local variables

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

    ! Update the ambient T and S

    if (plume%misomip_domain) then

       ! MISOMIP+ profiles, Eqs. 21 and 22
       plume%T_ambient = plume%T0 + (plume%Tbot - plume%T0) * (model%geometry%lsrf / plume%zbed_deep)
       plume%S_ambient = plume%S0 + (plume%Sbot - plume%S0) * (model%geometry%lsrf / plume%zbed_deep)

    else
       !TODO - Figure out how to do this for other domains
    endif

    !----------------------------------------------------------------
    ! Call the plume model to compute basal melt rates for floating ice
    !----------------------------------------------------------------

    call compute_plume(&
         ewn,                 nsn,                &
         dew,                 dns,                &
         itest,   jtest,      rtest,              &
         parallel,                                &
         model%geometry%thck,                     &
         model%geometry%lsrf,                     &
         model%geometry%topg,                     &
         model%climate%eus,                       &
         plume%T_ambient,     plume%S_ambient,    &
         plume%gammaT,        plume%gammaS,       &
         plume%S0,                                &  ! is this needed?
         plume%T_plume,       plume%S_plume,      &
         plume%D_plume,                           &
         plume%T_basal,       plume%S_basal,      &
         plume%u_plume,       plume%v_plume,      &
         plume%u_plume_Cgrid, plume%v_plume_Cgrid,&  ! is this needed?
         plume%ustar_plume,                       &
         plume%drho_plume,                        &
         plume%entrainment,   plume%detrainment,  &
         plume%divDu_plume,                       &
         model%basal_melt%bmlt_float)

    if (verbose_plume .and. main_task) write(iulog,*) 'Updated the plume'

  end subroutine glissade_plume_driver

!****************************************************

  subroutine compute_plume(&
       nx,               ny,               &
       dx,               dy,               &
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
       ustar_plume,                        &
       drho_plume,                         &
       entrainment,      detrainment,      &
       divDu_plume,                        &
       bmlt_float)

    ! Compute the melt rate at the ice-ocean interface from a steady-state plume model
    !
    ! References:
    !
    ! P.R. Holland and D.L. Feltham, 2006: The effects of rotation and ice shelf topography
    !    on frazil-laden ice shelf water plumes. J. Phys. Oceanog., 36, 2312-2327.
    ! P.R. Holland, A. Jenkins and D.M. Holland, 2008: The response of ice shelf
    !    basal melting to variations in ocean temperature. J. Climate, 21, 2558-2572.
    !
    ! TODO - Add Lambert (2023) and other references

    use glissade_masks, only: glissade_get_masks
!    use glissade_grid_operators, only: glissade_centered_gradient

    ! Input/output arguments

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

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
!         D_plume_cap,         & ! min(D_plume, H_cavity)
         dD_plume,            & ! change in D_plume (m)
         T_plume_old,         & ! T_plume from previous time step
         S_plume_old,         & ! S_plume from previous time step
         T_basal_old,         & ! T_basal from previous time step
         S_basal_old,         & ! S_basal from previous time step
         D_plume_old,         & ! D_plume from previous time step
         drho_plume_old,      & ! drho_plume from previous time step
         bmlt_float_old         ! melt rate from previous time step (m/s)

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

    ! plume model parameters
    ! Stable explicit time step TBD: try 10 minutes for now
    real(dp), parameter :: &
         dt_plume = 600.d0             ! time step (s) for continuity equation

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
    !TODO - Refine this mask?

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

       ! advance the time
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

       call compute_entrainment(&
            nx,         ny,      &
            dx,         dy,      &
            itest, jtest, rtest, &
            plume_mask,          &
            theta_slope,         &
            u_plume_east,        &
            v_plume_north,       &
            entrainment)

       ! Compute the detrainment rate where D_plume exceeds its max value

       call compute_detrainment(&
            nx,           ny,     &
            itest, jtest, rtest,  &
            H_cavity,             &
            D_plume,              &
            detrainment)

       ! Compute the basal melt rate, temperature and salinity at the plume-ice interface,
       ! given the plume velocity and entrainment rate.
       !Note: This subroutine currently updates T_plume and S_plume without any time lag.
       !TODO: Make T_plume and S_plume evolve incrementally.

       call compute_melt_rate(&
            nx,         ny,      &
            gammaT,              &
            gammaS,              &
            plume_mask,          &
            pressure,            &
            entrainment,         &
            u_plume_east,        &
            v_plume_north,       &
            T_ambient,           &
            S_ambient,           &
            T_basal,             &
            S_basal,             &
            T_plume,             &
            S_plume,             &
            itest, jtest, rtest, &
            ustar_plume,         &
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

       call compute_plume_transport(&
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
          bmlt_float_old = bmlt_float
          S_plume_old = S_plume
          T_plume_old = T_plume
          S_basal_old = S_basal
          T_basal_old = T_basal
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
!!       floating_mask,       &
!!       global_bndy_east,    &
!!       global_bndy_west,    &
!!       global_bndy_north,   &
!!       global_bndy_south,   &
!!       divu_mask_east,      &
!!       divu_mask_north,     &
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
!!         global_bndy_east,      & ! = 1 along east global boundary, else = 0
!!         global_bndy_west,      & ! = 1 along west global boundary, else = 0
!!         global_bndy_north,     & ! = 1 along north global boundary, else = 0
!!         global_bndy_south,     & ! = 1 along south global boundary, else = 0
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

  subroutine compute_entrainment(&
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

  end subroutine compute_entrainment

!****************************************************

  subroutine compute_detrainment(&
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

  end subroutine compute_detrainment

!****************************************************

  subroutine compute_melt_rate(&
       nx,         ny,      &
       gammaT,              &
       gammaS,              &
       plume_mask,          &
       pressure,            &
       entrainment,         &
       u_plume_east,        &
       v_plume_north,       &
       T_ambient,           &
       S_ambient,           &
       T_basal,             &
       S_basal,             &
       T_plume,             &
       S_plume,             &
       itest, jtest, rtest, &
       ustar_plume,         &
       bmlt_float)

    !TODO - Change to a 3-equation scheme, with Tp and Sp computed differently?
    !       Not sure how to handle advection for Tp and Sp.
    !       One option might be to pass in starting values (based on advection)
    !        and then solve below for dTp, with dTp inversely proportional to D_plume.
    !        I.e., the plume has a heat capacity.
    !
    !--------------------------------------------------------------------
    ! Compute the melt rate at the ice-ocean interface.
    !
    ! There are 5 equations for 5 unknowns: m, Tb, Sb, Tp and Sp
    ! where m = melt rate at ice-ocean interface
    !       Tb = potential temperature at ice-ocean interface
    !       Sb = salinity at ice-ocean interface
    !       Tp = potential temperature of boundary-layer plume
    !       Sp = salinity of boundary-layer plume
    ! 
    ! (1) rhow * m * L  = rhoo * cw * u_fric * gammaT * (Tp - Tb)
    ! (2) rhow * m * Sb = rhoo * u_fric * gammaS *(Sp - Sb)
    ! (3) Tb = lambda1*Sb + lambda2 + lambda3*pb 
    ! (4) L * m = -cw * e * (Tp - Ta)
    ! (5) Sp * m = -e * (Sp - Sa)
    !
    ! Eq. 1 and 2 describe heat and salt transfer at the ice-ocean interface.
    ! Eq. 3 is the linearized liquidus relation that determines the potential freezing point.
    ! Eq. 4 and 5 describe heat and salt entrainment from the ambient ocean to the boundary-layer plume,
    !  where Ta and Sa are the potential temperature and salinity of the ambient ocean.
    !
    ! We can rewrite (1) and (2) as
    !
    ! (1)     m = T_factor * (Tp - Tb)
    ! (2)  Sb*m = S_factor * (Sp - Sb)
    !
    ! where T_factor = (rhoo * cw * ufric * gammaT) / (rhow * L)
    !       S_factor = (rhoo * ufric * gammaS) / rhow
    !
    ! Rearrange (4):  Tp = Ta - (L/(cw*e)) * m
    ! 
    ! Use (3) and (4) to replace Tp and Tb in (1):
    !
    ! (1')    m = m1*Sb + m2
    ! where  m1 = -T_factor*lambda1/denom
    !        m2 =  T_factor*(Ta - lambda2 - lambda3*p)/denom
    !     denom = 1 + T_factor*L/(cw*e)
    ! 
    ! Use (5) to replace S in (2):
    !
    ! (2')   Sb = S_factor*e*Sa / ((m+S_factor)*(m+e))
    !
    ! Use (2') to replace Sb in (1') to form a cubic equation for m:
    !
    ! (1'')  a*m^3 + b*m^2 + c*m + d = 0
    !
    !   where a = 1
    !         b = S_factor + e - m2
    !         c = S_factor*e - m2*(S_factor + e)
    !         d = -S_factor*e*(m1*Sa + m2)
    !
    ! Use the cubic_solver subroutine to find m.
    !
    ! Given m, back out the other 4 unknowns.
    !--------------------------------------------------------------------
    
    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    integer, intent(in) ::  &
         itest, jtest, rtest    ! test cell coordinates (diagnostic only)

    ! Note: gammaS and gammaT are config parameters and are passed in as arguments.
    !       Other MISOMIP parameters are declared at the top of the module.
    
    real(dp), intent(in) ::  &
         gammaT,              & ! nondimensional heat transfer coefficient
         gammaS                 ! nondimensional salt transfer coefficient
    
    integer, dimension(nx,ny), intent(in) :: &
         plume_mask             ! = 1 for cells where scalar plume variables are computed

    real(dp), dimension(nx,ny), intent(in) :: &
         pressure,            & ! ocean pressure at base of ice (N/m^2)
         entrainment,         & ! entrainment rate of ambient water into plume (m/s)
         u_plume_east,        & ! u_plume on east edges (m/s)
         v_plume_north,       & ! v_plume on north edges (m/s)
         T_ambient,           & ! ambient ocean potential temperature at depth of ice-ocean interface (deg C)
         S_ambient              ! ambient ocean salinity at depth of ice-ocean interface (psu)
    
    real(dp), dimension(nx,ny), intent(out) :: &
         ustar_plume,         & ! plume friction velocity (m/s) on ice grid, output as a diagnostic
         T_basal,             & ! basal ice temperature (deg C)
         S_basal,             & ! basal ice salinity (psu)
         T_plume,             & ! plume temperature (deg C)
         S_plume,             & ! plume salinity (psu)
         bmlt_float             ! melt rate at base of floating ice (m/s)
    
    ! local variables
    
    real(dp) :: &
         u_plume, v_plume,    & ! plume velocity components at cell center (m/s)
         plume_speed,         & ! plume speed at cell center (m/s)
         T_factor, S_factor,  & ! factors in melt-rate equations
         denom,               & ! denominator
         m1, m2,              & ! factors in relation between m and Sb
         ma, mb, mc, md,      & ! coefficients in cubic equation for m
         bmlt_float_avg         ! average value of bmlt_float in main cavity
    
    integer :: i, j
    
    !WHL - debug -  Test cubic solver
!       ma =    2.d0
!       mb =  -30.d0
!       mc =  162.d0
!       md = -350.d0
!       call cubic_solver(ma, mb, mc, md, solution)
!       write(iulog,*) 'Trial cubic solution =', solution
!       write(iulog,*) 'True solution =', (10.d0 + sqrt(108.d0))**(1.d0/3.d0) - (-10.d0 + sqrt(108.d0))**(1.d0/3.d0) + 5.d0


    ! Loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          
          if (plume_mask(i,j) == 1 .and. entrainment(i,j) > 0.0d0) then
             
             ! Interpolate the plume speed to the cell center, and compute the friction velocity ustar.
             
             u_plume = (u_plume_east(i,j) + u_plume_east(i-1,j)) / 2.0d0
             v_plume = (v_plume_north(i,j) + v_plume_north(i,j-1)) / 2.0d0
             plume_speed = sqrt(u_plume**2 + v_plume**2 + u_tidal**2)
             ustar_plume(i,j) = sqrt(c_drag) * plume_speed

             T_factor = (rhoo * cpw * ustar_plume(i,j) * gammaT) / (rhow * lhci)
             S_factor = (rhoo * ustar_plume(i,j) * gammaS) / rhow
             
             denom = 1.d0 + (T_factor*lhci)/(cpw*entrainment(i,j))
             m1 = -lambda1 * T_factor / denom
             m2 = T_factor * (T_ambient(i,j) - lambda2 - lambda3*pressure(i,j)) / denom
             
             ma = 1.d0
             mb = S_factor + entrainment(i,j) - m2
             mc = S_factor*entrainment(i,j) - m2*(S_factor + entrainment(i,j))
             md = -S_factor*entrainment(i,j)*(m1*S_ambient(i,j) + m2)
             
             ! Solve the cubic equation
             call cubic_solver(&
                  ma, mb, mc, md, &
                  bmlt_float(i,j))

             if (verbose_plume .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'Melt rate calc: rank, i, j =', rtest, i, j
                write(iulog,*) 'pressure (Pa) =', pressure(i,j)
                write(iulog,*) 'T_factor (m/s/deg), S_factor (m/s)=', T_factor, S_factor
                write(iulog,*) 'entrainment (m/s) =', entrainment(i,j)
                write(iulog,*) 'm1 (m/s/psu) =', m1
                write(iulog,*) 'm2 (m/s) =', m2
                write(iulog,*) 'denom =', denom
                write(iulog,*) 'a, b, c, d =', ma, mb, mc, md
                write(iulog,*) 'residual of cubic solve =', ma*bmlt_float(i,j)**3 + mb*bmlt_float(i,j)**2 + mc*bmlt_float(i,j) + md
             endif
             
             ! Given the melt rate, compute T_basal and S_basal
!               S_basal(i,j) = (S_factor * entrainment(i,j) * S_ambient(i,j)) /  &
!                               ( (bmlt_float(i,j) + S_factor) * (bmlt_float(i,j) + entrainment(i,j)) )
             S_basal(i,j) = (bmlt_float(i,j) - m2) / m1
             T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)

             ! Given m, compute T_plume and S_plume
             T_plume(i,j) = T_ambient(i,j) - (lhci/(cpw*entrainment(i,j))) * bmlt_float(i,j)
             S_plume(i,j) = S_ambient(i,j) * entrainment(i,j) / (bmlt_float(i,j) + entrainment(i,j))

             !WHL - debug - check for NaNs
             if (T_plume(i,j) /= T_plume(i,j) .or. S_plume(i,j) /= S_plume(i,j) .or. &
                 T_basal(i,j) /= T_basal(i,j) .or. S_basal(i,j) /= S_basal(i,j) .or. &
                 bmlt_float(i,j) /= bmlt_float(i,j)) then
                write(iulog,*) 'Bad values, i, j =', i, j
                write(iulog,*) 'T_plume, S_plume:', T_plume(i,j), S_plume(i,j)
                write(iulog,*) 'T_basal, S_basal:', T_basal(i,j), S_basal(i,j)
                write(iulog,*) 'bmlt_float:', bmlt_float(i,j)
                stop
             endif

          else    ! plume_mask = 0
             
             bmlt_float(i,j) = 0.0d0
             
             S_plume(i,j) = S_ambient(i,j)
             T_plume(i,j) = T_ambient(i,j)
             
             S_basal(i,j) = S_ambient(i,j)
             T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
             
          endif   ! plume_mask and entrainment > 0
          
       enddo   ! i
    enddo   ! j

  end subroutine compute_melt_rate

!****************************************************
    
  !TODO - Move this subroutine to a utility module?
  !TODO - Pass 3 complex roots in and out.
  subroutine cubic_solver(&
       a, b, c, d, &
       x1,         &
       x2_r, x2_i, &
       x3_r, x3_i)

    !------------------------------------------------
    ! Find the real root of a cubic equation:
    !
    !    ax^3 + bx^2 + cx = d = 0
    !
    ! Do this by making the substitution
    !
    !    x = y - b/(3a)
    !
    ! to convert to a depressed cubic:
    !
    !    y^3 + py + q = 0
    !
    ! where p = (1/a) * (c - b^2/(3a))
    !       q = (1/a) * (d + 2b^3/(27a^2) - bc/(3a))
    !
    !------------------------------------------------

    real(dp), intent(in) ::  &
         a, b, c, d       ! coefficients of cubic equation
                          ! assumed to be real

    real(dp), intent(out) ::  &
         x1               ! real solution of cubic equation

    real(dp), intent(out), optional ::  &
         x2_r, x2_i,    & ! other solutions of cubic equation
         x3_r, x3_i       ! could be either real or complex

    real(dp) :: &
         p, q             ! coefficients of depressed cubic

    real(dp) :: &
         Delta            ! discriminant

    real(dp) :: &
         y1,            & ! solutions of depressed cubic
         y2_r, y2_i,    & !
         y3_r, y3_i

    real(dp) :: &
         u, v,          & ! some intermediate factors
         fu, fv,        &
         phi

    real(dp), parameter :: &
         p333 = 1.d0/3.d0

    !WHL - debug
    logical, parameter :: verbose_cubic = .false.

    ! compute coefficients of depressed cubic, y^3 + py + q = 0

    p = (3.d0*c/a - (b/a)**2) / 3.d0
    q = (2.d0*(b/a)**3 - 9.d0*b*c/(a*a) + 27.d0*d/a) / 27.d0

    ! compute the discriminant
    Delta = (p/3.d0)**3 + (q/2.d0)**2

    if (verbose_cubic) then
       write(iulog,*) 'Delta =', Delta
       if (Delta > 0.d0) then
          write(iulog,*) 'One real root, 2 complex conjugate'
       elseif (Delta == 0.d0) then
          write(iulog,*) 'Three real roots of which at least two are equal'
       elseif (Delta < 0.d0) then
          write(iulog,*) 'Three distinct real roots'
       endif
    endif

    if (Delta >= 0.d0) then   

       if (Delta > 0.d0) then    ! one real root, two complex roots
          fu = -q/2.d0 + sqrt(Delta)
          fv = -q/2.d0 - sqrt(Delta)
       else  ! Delta = 0; three real roots of which at least two are equal
          fu = -q/2.d0
          fv = fu
       endif
 
       ! some logic to avoid taking cube roots of negative numbers
       if (fu >= 0.d0) then
          u = fu**p333
       else
          u = -(-fu)**p333
       endif

       if (fv >= 0.d0) then
          v = fv**p333
       else
          v = -(-fv)**p333
       endif

       ! form solutions of depressed cubic
       y1 = u + v       ! real
       y2_r = -(u+v)/2.d0
       y2_i =  (u-v)*sqrt(3.d0)/2.d0
       y3_r = -(u+v)/2.d0
       y3_r = -(u-v)*sqrt(3.d0)/2.d0

       if (verbose_cubic) then
          write(iulog,*) 'a, b, c, d:', a, b, c, d
          write(iulog,*) 'p, q:', p, q
          write(iulog,*) 'y1 =', y1
          write(iulog,*) 'x1 =', x1
       endif

    else  ! Delta < 0; three distinct real roots
          ! use a trigonometric formulation

       phi = acos(-q/(2.d0*sqrt(abs(p)**3/27.d0)))

       y1 =    2.d0 * sqrt(abs(p)/3.d0) * cos(phi/3.d0)
       y2_r = -2.d0 * sqrt(abs(p)/3.d0) * cos((phi+pi)/3.d0)
       y2_i =  0.d0
       y3_r = -2.d0 * sqrt(abs(p)/3.d0) * cos((phi-pi)/3.d0)
       y3_i =  0.d0

       if (verbose_cubic) then
          write(iulog,*) 'a, b, c, d:', a, b, c, d
          write(iulog,*) 'p, q:', p, q
          write(iulog,*) 'y1, y2, y3 =', y1, y2_r, y3_r
          write(iulog,*) 'b/3a =', b/(3.d0*a)
          write(iulog,*) 'x1 =', y1 - b/(3.d0*a)
          write(iulog,*) 'x2 =', y2_r - b/(3.d0*a)
          write(iulog,*) 'x3 =', y3_r - b/(3.d0*a)
       endif

    endif

    ! Recover the solutions
    ! Mostly likely we are only interested in x1, but compute the others if requested

    x1 = y1 - b/(3.d0*a)

    if (present(x2_r) .and. present(x2_i) .and. present(x3_r) .and. present(x3_i)) then
       x2_r = y2_r - b/(3.d0*a)
       x2_i = y2_i
       x3_r = y3_r - b/(3.d0*a)
       x3_i = y3_i
    endif

  end subroutine cubic_solver

!****************************************************

  subroutine compute_plume_transport(&
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
             dDT = entrainment(i,j)*T_ambient(i,j) - detrainment(i,j)*T_plume(i,j) + bmlt_float(i,j)*T_basal(i,j)
             work(i,j,2) = work(i,j,2) + dDT*dt
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

  end subroutine compute_plume_transport

!****************************************************

  end module glissade_plume

!****************************************************
