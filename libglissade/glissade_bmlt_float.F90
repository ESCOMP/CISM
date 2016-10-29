!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_bmlt_float.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
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
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

module glissade_bmlt_float

  use glimmer_global, only: dp
  use glimmer_physcon, only: rhoo, rhow, grav, lhci, scyr, pi
  use glimmer_log
  use glide_types
  use parallel

  implicit none
  
  private
  public :: glissade_basal_melting_float

  logical :: verbose_velo = .false.
  logical :: verbose_continuity = .true.
  logical :: verbose_melt = .false.

contains

!****************************************************

  subroutine glissade_basal_melting_float(whichbmlt_float,              &
                                          ewn,         nsn,             &
                                          dew,         dns,             &
                                          itest,       jtest,    rtest, &
                                          x1,                           &
                                          thck,        lsrf,            &
                                          topg,        eus,             &
                                          bmltfloat)

    use glissade_masks, only: glissade_get_masks
    use glimmer_paramets, only: tim0, thk0

    ! Compute the rate of basal melting for floating ice by one of several methods.

    !-----------------------------------------------------------------
    ! Input/output arguments
    !-----------------------------------------------------------------

    integer, intent(in) :: whichbmlt_float            ! method for computing melt rate of floating ice

    integer, intent(in) ::  &
         ewn, nsn,             & ! grid dimensions
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dew, dns                ! grid spacing in x, y (m)

    real(dp), dimension(:), intent(in) :: &
         x1                      ! x1 grid coordinates (m), ice grid
                                 ! used with bmlt_float_xlim for MISMIP+ Ice2r

    real(dp), dimension(:,:), intent(in) :: &
         thck,                 & ! ice thickness (m)
         lsrf,                 & ! elevation of lower ice surface (m)
         topg                    ! elevation of bed topography (m)

    real(dp), intent(in) :: &
         eus                     ! eustatic sea level (m), = 0. by default

    type(glide_bmltfloat), intent(inout) :: &
         bmltfloat               ! derived type with fields and parameters for MISMIP+ and MISOMIP

    !-----------------------------------------------------------------
    ! Note: The bmltfloat derived type includes the 2D output field bmlt_float,
    !       along with several 2D fields used in the plume model.
    ! 
    ! Also, bmltfloat includes a number of prescribed parameters for MISMIP+ and MISOMIP:
    !
    !       MISMIP+ Ice1
    !       - bmlt_float_omega     ! time scale for basal melting (s-1), default = 0.2/yr
    !       - bmlt_float_h0        ! scale for sub-shelf cavity thickness (m), default = 75 m
    !       - bmlt_float_z0        ! scale for ice draft (m), default = -100 m
    !
    !       MISMIP+ Ice2
    !       - bmlt_float_rate      ! constant melt rate (m/s), default = 100 m/yr
    !       - bmlt_float_xlim      ! melt rate = 0 for abs(x) < bmlt_float_xlim (m), default = 480000 m
    !
    !       MISOMIP
    !       - T0                   ! sea surface temperature (deg C), default = -1.9 C
    !       - Tbot                 ! temperature at the sea floor (deg C), default = 1.0 C (warm), -1.9 C (cold)
    !       - S0                   ! sea surface salinity (psu), default = 33.8 psu
    !       - Sbot                 ! salinity at the sea floor (psu), default = 34.7 psu (warm), 34.55 psu (cold)
    !       - zbed_deep            ! min sea floor elevation (m), default = -720 m
    !       - gammaT               ! nondimensional heat transfer coefficient, default = 5.0e-2
    !       - gammaS               ! nondimensional salt transfer coefficient, default = gammaT/35
    !
    ! Note: Basal melt rates are > 0 for melting, < 0 for freeze-on
    !-----------------------------------------------------------------
 
    !----------------------------------------------------------------
    ! Local variables and pointers set to components of bmltfloat derived type
    !----------------------------------------------------------------      

    real(dp), dimension(:,:), pointer :: &
         bmlt_float,          & ! basal melt rate for floating ice (m/s) (> 0 for melt, < 0 for freeze-on)
                                ! Note: converted from m/s to scaled model units on output
         T_basal,             & ! basal ice temperature; at freezing point (deg C)
         S_basal,             & ! basal ice salinity; at freezing point (psu)
         u_plume,             & ! x component of plume velocity (m/s) at cell centers
         v_plume,             & ! y component of plume velocity (m/s) at cell centers
         u_plume_Cgrid,       & ! x component of plume velocity (m/s) on C grid (east edges)
         v_plume_Cgrid,       & ! y component of plume velocity (m/s) on C grid (east edges)
         ustar_plume,         & ! plume friction velocity (m/s)
         drho_plume,          & ! density difference between plume and ambient ocean (kg/m3)
         T_plume,             & ! plume temperature (deg C)
         S_plume,             & ! plume salinity (psu)
         D_plume,             & ! plume thickness (m)
         entrainment,         & ! entrainment rate of ambient water into plume (m/s)
         detrainment,         & ! detrainment rate of plume into ambient water (m/s)
                                ! Note: entrainment/detrainment rates are converted from m/s to scaled model units on output
         divDu_plume,         & ! divergence of D_plume*u_plume (m/s)
         T_ambient,           & ! ambient ocean temperature below ice and plume (deg C)
         S_ambient              ! ambient ocean salinity below ice and plume (psu)

    ! parameters for MISMIP+ Ice1 experiments
    real(dp) :: &
         bmlt_float_omega,  & ! time scale for basal melting (s-1)
         bmlt_float_h0,     & ! scale for sub-shelf cavity thickness (m)
         bmlt_float_z0        ! scale for ice draft (m)

    ! parameters for MISMIP+ Ice2 experiments
    real(dp) :: &
         bmlt_float_rate,   & ! constant melt rate (m/s)
         bmlt_float_xlim      ! melt rate = 0 for abs(x) < bmlt_float_xlim (m)

    ! parameters for MISOMIP
    real(dp) ::  &
         T0,                & ! sea surface temperature (deg C)
         Tbot,              & ! temperature at the sea floor (deg C)
         S0,                & ! sea surface salinity (psu)
         Sbot,              & ! salinity at the sea floor  (psu)
         zbed_deep,         & ! min sea floor elevation (m)
         gammaT,            & ! nondimensional heat transfer coefficient
         gammaS               ! nondimensional salt transfer coefficient

    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------

    integer, dimension(ewn,nsn) ::  &
         ice_mask,            & ! = 1 where ice temperature is computed (thck > thklim), else = 0
         floating_mask,       & ! = 1 where ice is present and floating, else = 0
         ocean_mask             ! = 1 where topg is below sea level and ice is absent

    integer :: i, j
    real(dp) :: h_cavity        ! depth of ice cavity beneath floating ice (m)
    real(dp) :: z_draft         ! draft of floating ice (m below sea level)

!!    logical, parameter :: verbose_bmlt = .false.
    logical, parameter :: verbose_bmlt = .true.

!TODO - Make first_call depend on whether we are restarting
!!    logical :: first_call = .false.
    logical :: first_call = .true.

    !----------------------------------------------------------------      
    ! Assign local pointers and variables to components of bmltfloat derived type
    !----------------------------------------------------------------      

    ! MISMIP+
    bmlt_float_omega = bmltfloat%bmlt_float_omega
    bmlt_float_h0 = bmltfloat%bmlt_float_h0
    bmlt_float_z0 = bmltfloat%bmlt_float_z0
    bmlt_float_rate = bmltfloat%bmlt_float_rate
    bmlt_float_xlim = bmltfloat%bmlt_float_xlim

    ! MISOMIP
    T0 = bmltfloat%T0 
    Tbot = bmltfloat%Tbot 
    S0 = bmltfloat%S0 
    Sbot = bmltfloat%Sbot 
    zbed_deep = bmltfloat%zbed_deep
    gammaT = bmltfloat%gammaT
    gammaS = bmltfloat%gammaS

    bmlt_float  => bmltfloat%bmlt_float

    if (whichbmlt_float == BMLT_FLOAT_MISOMIP) then
       ! the following fields are used or computed by the plume model
       T_basal       => bmltfloat%T_basal
       S_basal       => bmltfloat%S_basal
       u_plume       => bmltfloat%u_plume
       v_plume       => bmltfloat%v_plume
       u_plume_Cgrid => bmltfloat%u_plume_Cgrid
       v_plume_Cgrid => bmltfloat%v_plume_Cgrid
       ustar_plume   => bmltfloat%ustar_plume
       drho_plume    => bmltfloat%drho_plume
       T_plume       => bmltfloat%T_plume
       S_plume       => bmltfloat%S_plume
       D_plume       => bmltfloat%D_plume
       entrainment   => bmltfloat%entrainment
       detrainment   => bmltfloat%detrainment
       divDu_plume   => bmltfloat%divDu_plume
       T_ambient     => bmltfloat%T_ambient
       S_ambient     => bmltfloat%S_ambient
    endif

    if (main_task .and. verbose_bmlt) print*, 'Computing bmltfloat, whichbmlt_float =', whichbmlt_float

    ! Compute masks: 
    ! - ice_mask = 1 where thck > thklim
    ! - floating_mask = 1 where ice is present and floating;
    ! - ocean_mask = 1 where topg is below sea level and ice is absent
    !Note: The '0.0d0' argument is thklim. Here, any ice with thck > 0 gets ice_mask = 1.

    call glissade_get_masks(ewn,           nsn,           &
                            thck,          topg,          &
                            eus,           0.0d0,         &
                            ice_mask,      floating_mask, &
                            ocean_mask)

    ! Compute the basal melt rate for floating ice

    bmlt_float(:,:) = 0.0d0

    if (whichbmlt_float == BMLT_FLOAT_NONE) then

       ! nothing to do; bmlt_float already set to zero

    elseif (whichbmlt_float == BMLT_FLOAT_CONSTANT) then

       ! Set melt rate to a constant value for floating ice.
       ! This includes ice-free ocean cells, in case ice is advected to those cells by the transport scheme.

       do j = 1, nsn
          do i = 1, ewn

             if (floating_mask(i,j) == 1 .or. ocean_mask(i,j) == 1) then   ! ice is present and floating, or ice-free ocean

                ! Note: For MISMIP+ experiment Ice2r, melting is masked out where x < 480 km

                if (abs(x1(i)) >= bmlt_float_xlim) then   ! melting is allowed
                   bmlt_float(i,j) = bmlt_float_rate
                endif

                !WHL - debug
                if (j==jtest .and. this_rank==rtest) then
!!                if (i==itest .and. j==jtest .and. this_rank==rtest) then
!!                   print*, 'rank, i, j, bmlt_float:', this_rank, i, j, bmlt_float(i,j)
                endif
                   
             endif   ! ice is present and floating

          enddo
       enddo

    elseif (whichbmlt_float == BMLT_FLOAT_MISMIP) then

       ! compute melt rate based on bed depth and cavity thickness
       !
       ! The MISMIP+ formula is as follows:
       !
       ! bmlt_float = omega * tanh(H_c/H_0) * max(z_0 - z_d, 0)
       !
       ! where H_c = lsrf - topg is the cavity thickness
       !       z_d = lsrf - eus is the ice draft
       !       omega = a time scale = 0.2 yr^{-1} by default
       !       H_0 = 75 m by default
       !       z_0 = 100 m by default

       do j = 1, nsn
          do i = 1, ewn

             ! allow basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme

             if (floating_mask(i,j) == 1 .or. ocean_mask(i,j) == 1) then   ! ice is present and floating, or ice-free ocean

                h_cavity = lsrf(i,j) - topg(i,j)
                z_draft = lsrf(i,j) - eus
                bmlt_float(i,j) = bmlt_float_omega * tanh(h_cavity/bmlt_float_h0) * max(bmlt_float_z0 - z_draft, 0.0d0)

                   !WHL - debug
!!                   if (ns == 5) then
!!                      print*, 'cavity, tanh, draft, d_draft, melt rate (m/yr):', i, j, h_cavity, tanh(h_cavity/bmlt_float_h0), &
!!                           z_draft, max(bmlt_float_z0 - z_draft, 0.d0), bmlt_float(i,j)*31536000.d0
!!                   endif

             endif   ! ice is present and floating

          enddo
       enddo

    elseif (whichbmlt_float == BMLT_FLOAT_MISOMIP) then

       ! Compute melt rates using a plume model, given vertical profiles of T and S in the ambient ocean
       !
       ! See this paper for details:
       ! X. S. Asay-Davis et al. (2016), Experimental design for three interrelated 
       !    marine ice sheet and ocean model intercomparison projects: 
       !    MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
       !    Geosci. Model Devel., 9, 2471-2497, doi: 10.5194/gmd-9-2471-2016.

       !WHL - debug
       print*, 'itest, jtest =', itest, jtest

       !WHL - debug
       print*, ' '
       print*, 'lsrf field, rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f12.5)',advance='no') lsrf(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'topg field, rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f12.5)',advance='no') topg(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'lsrf - topg, rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f12.5)',advance='no') lsrf(i,j) - topg(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'floating_mask, rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i8)',advance='no') floating_mask(i,j)
          enddo
          write(6,*) ' '
       enddo

       ! Given the ice draft in each floating grid cell, compute the ambient ocean T and S
       !  using the prescribed MISOMIP profile.

       where (floating_mask == 1) 

          ! MISOMIP+ profiles, Eqs. 21 and 22 
          T_ambient(:,:) = T0 + (Tbot - T0) * (lsrf(:,:) / zbed_deep)
          S_ambient(:,:) = S0 + (Sbot - S0) * (lsrf(:,:) / zbed_deep)

       elsewhere

          T_ambient = T0
          S_ambient = S0

       endwhere

       ! Note: The plume model expects floating_mask, T_ambient and S_ambient to be correct in halo cells.
       !       This is likely the case already, but do halo updates just in case.
       
       call parallel_halo(floating_mask)
       call parallel_halo(T_ambient)
       call parallel_halo(S_ambient)

       ! If D_plume has already been computed, then convert from scaled units to meters
       if (.not. first_call) then
          D_plume(:,:) = D_plume(:,:)*thk0
       endif

       !----------------------------------------------------------------
       ! Call the plume model to compute basal melt rates for floating ice
       !----------------------------------------------------------------
       
       call glissade_plume_melt_rate(&
            first_call,                         &
            ewn,              nsn,              &
            dew,              dns,              &
            x1,                                 &
            lsrf,             topg,             &
            floating_mask,                      &
            itest, jtest,     rtest,            &
            T_ambient,        S_ambient,        &
            gammaT,           gammaS,           &
            S0,                                 &
            T_basal,          S_basal,          &
            u_plume,          v_plume,          &
            u_plume_Cgrid,    v_plume_Cgrid,    &
            D_plume,                            &
            ustar_plume,      drho_plume,       &
            T_plume,          S_plume,          &
            entrainment,      detrainment,      &
            divDu_plume,      bmlt_float)

    endif   ! whichbmlt_float

    ! Set first_call to false. 
    ! Next time, the plume variables just computed (T_plume, S_plume, D_plume)
    !  will be taken as initial conditions.
    first_call = .false.

    ! Convert output parameters from SI units to scaled model parameters
    bmlt_float(:,:) = bmlt_float(:,:) * tim0/thk0

    if (whichbmlt_float == BMLT_FLOAT_MISOMIP) then
       entrainment(:,:) = entrainment(:,:) * tim0/thk0
       detrainment(:,:) = detrainment(:,:) * tim0/thk0
       divDu_plume(:,:) = divDu_plume(:,:) * tim0/thk0
       D_plume(:,:) = D_plume(:,:)/thk0
    endif

    !WHL - debug
    if (this_rank == rtest .and. verbose_bmlt) then
       i = itest
       j = jtest
       print*, 'After rescaling, bmlt_float =', bmlt_float(i,j)
    endif

  end subroutine glissade_basal_melting_float

!****************************************************

  subroutine glissade_plume_melt_rate(&
       first_call,                      &
       nx,               ny,            &
       dx,               dy,            &
       x1,                              &
       lsrf,             topg,          &
       floating_mask,                   &
       itest, jtest,     rtest,         &
       T_ambient,        S_ambient,     &
       gammaT,           gammaS,        &
       S0,                              &
       T_basal,          S_basal,       &
       u_plume,          v_plume,       &
       u_plume_Cgrid,    v_plume_Cgrid, &
       D_plume,                         &
       ustar_plume,      drho_plume,    &
       T_plume,          S_plume,       &
       entrainment,      detrainment,   &
       divDu_plume,      bmlt_float)

    use glissade_grid_operators, only: &
         glissade_stagger, glissade_unstagger, glissade_centered_gradient

    ! Compute the melt rate at the ice-ocean interface based on a steady-state plume model

    ! References:
    !
    ! P.R. Holland and D.L. Feltham, 2006: The effects of rotation and ice shelf topography 
    !    on frazil-laden ice shelf water plumes. J. Phys. Oceanog., 36, 2312-2327.
    !    
    ! P.R. Holland, A. Jenkins and D.M. Holland, 2008: The response of ice shelf 
    !    basal melting to variations in ocean temperature. J. Climate, 21, 2558-2572. 

    ! Input/output arguments

    logical, intent (in) :: &
         first_call             ! if true, then use simple initial conditions to start the plume calculation
                                ! if false, then start from the input values of T_plume, S_plume and D_plume

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    real(dp), dimension(:), intent(in) :: &
         x1                     ! x1 grid coordinates (m), ice grid

    real(dp), dimension(nx,ny), intent(in) ::  &
         lsrf,                & ! ice lower surface elevation (m, negative below sea level)
         topg                   ! bedrock elevation (m, negative below sea level)

    integer, dimension(nx,ny), intent(in) :: &
         floating_mask          ! = 1 where ice is floating and melt rates are computed, else = 0

    integer, intent(in) :: &
         itest, jtest, rtest

    real(dp), dimension(nx,ny), intent(in) ::  &
         T_ambient,           & ! ambient ocean potential temperature at depth of ice-ocean interface (deg C)
         S_ambient              ! ambient ocean salinity at depth of ice-ocean interface (psu)

    real(dp), intent(in) :: &
         gammaT,              & ! nondimensional heat transfer coefficient
         gammaS,              & ! nondimensional salt transfer coefficient
         S0                     ! sea surface salinity (psu)

    ! Note: T_plume, S_plume and D_plume can either be initialized below (if first_call = F)
    !       or passed in (if first_call = T).
    real(dp), dimension(nx,ny), intent(inout) :: &
         T_plume,             & ! plume temperature (deg C)
         S_plume,             & ! plume salinity (psu)
         D_plume                ! plume thickness (m)

    !TODO - Switch u_plume and v_plume to the ice grid for diagnostics
    ! Note: Plume velocities are computed on the C grid, and then are interpolated
    !       to cell centers as a diagnostic.

    real(dp), dimension(nx,ny), intent(out) :: &
         u_plume,             & ! x component of plume velocity (m/s) at cell centers
         v_plume,             & ! y component of plume velocity (m/s) at cell centers
         u_plume_Cgrid,       & ! x component of plume velocity (m/s) on C grid (east edges)
         v_plume_Cgrid,       & ! y component of plume velocity (m/s) on C grid (north edges)
         ustar_plume,         & ! plume friction velocity (m/s) on ice grid
         drho_plume,          & ! density difference between plume and ambient ocean (kg/m^3)
         T_basal,             & ! basal ice temperature (deg C)
         S_basal,             & ! basal ice salinity (psu)
         entrainment,         & ! entrainment rate of ambient water into plume (m/s)
         detrainment,         & ! detrainment rate of plume into ambient water (m/s)
         divDu_plume,         & ! div(Du) for plume
         bmlt_float             ! melt rate at base of floating ice (m/s)

    ! Local variables

    real(dp), dimension(nx,ny) :: &
         pressure,            & ! ocean pressure at base of ice (N/m^2)
         lsrf_plume,          & ! elevation of plume-ambient interface (m, negative below sea level)
         rho_plume,           & ! plume density (kg/m^3)
         rho_ambient,         & ! ambient ocean density (kg/m^3)
         eta_plume,           & ! displacement of plume surface, D_plume - H_cavity (m)
         dD_plume,            & ! change in D_plume (m)
         T_plume_old,         & ! T_plume from previous iteration
         S_plume_old,         & ! S_plume from previous iteration
         T_basal_old,         & ! T_basal from previous iteration
         S_basal_old,         & ! S_basal from previous iteration
         D_plume_old,         & ! D_plume from previous iteration
         eta_plume_old,       & ! eta_plume from previous iteration
         entrainment_old,     & ! entrainment from previous iteration
         divDu_plume_old,     & ! divDu_plume from previous iteration
         bmlt_float_old         ! melt rate from previous iteration (m/s)

    real(dp), dimension(nx,ny) ::  &
         u_plume_east,          & ! u_plume on east edges
         v_plume_east,          & ! v_plume on east edges
         u_plume_east_old,      & ! old values of u_plume on east edges
         v_plume_east_old,      & ! old values of v_plume on east edges
         plume_speed_east,      & ! plume speed on east edges
         entrainment_east,      & ! entrainment on east edges
         u_plume_north,         & ! u_plume on north edges
         v_plume_north,         & ! v_plume on north edges
         u_plume_north_old,     & ! old values of u_plume on north edges
         v_plume_north_old,     & ! old values of v_plume on north edges
         plume_speed_north,     & ! plume speed on north edges
         entrainment_north        ! entrainment on north edges

    real(dp), dimension(nx-1,ny-1) :: &
         plume_speed,         & ! plume speed (m/s)
         dlsrf_dx,            & ! horizontal gradient of lsrf
         dlsrf_dy,            & !
         dlsrf_plume_dx,      & ! horizontal gradient of lsrf_plume
         dlsrf_plume_dy,      & !
         deta_plume_dx,       & ! horizontal gradient of eta_plume
         deta_plume_dy

    integer, dimension(nx,ny) :: &
         plume_mask_cell,     & ! = 1 at vertices where scalar plume variables are computed
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed
         edge_mask_north        ! = 1 on north edges where plume velocity is computed

    integer, dimension(nx-1,ny-1) :: &
         plume_mask_velo        ! = 1 at vertices where plume velocity is computed, else = 0

    ! Note: The north and south global_bndy masks are used for ISOMIP+.
    !       Not sure if they would be needed in a more realistic run.
    integer, dimension(nx,ny) :: &
         global_bndy_east,    & ! = 1 along east global boundary, else = 0
         global_bndy_west,    & ! = 1 along west global boundary, else = 0
         global_bndy_north,   & ! = 1 along north global boundary, else = 0
         global_bndy_south      ! = 1 along south global boundary, else = 0

    !TODO - Are these work masks needed?
    integer, dimension(nx,ny) :: &
         ice_mask_work          ! work mask on ice grid

    integer, dimension(nx-1,ny-1) :: &
         velo_mask_work         ! work mask on velocity grid

    real(dp) :: &
         lsrf_min,                      & ! global min value of lsrf (m)
         H_cavity,                      & ! thickness of ocean cavity beneath the plume (m)
         pgf_x, pgf_y,                  & ! components of pressure gradient force (N/m^2)
         slope,                         & ! slope of plume-ambient interface (unitless)
         theta_slope,                   & ! atan of slope of plume-ambient interface (rad)
         plume_tendency,                & ! tendency of plume thickness (m/s)
         D_plume_east, D_plume_west,    & ! terms in discretization of plume divergence
         D_plume_north, D_plume_south,  &
         dDu_dx, dDv_dy,                &
         eta_plume_unrelaxed              ! value of eta_plume without relaxation

    real(dp) ::   &
         D_plume_edge,          & ! D_plume averaged to edge
         grav_reduced_edge,     & ! reduced gravity on edge
         dlsrf_dx_edge,         & ! dlsrf_dx on edge
         dlsrf_dy_edge,         & ! dlsrf_dy on edge
         dlsrf_plume_dx_edge,   & ! dlsrf_plume_dx on edge
         dlsrf_plume_dy_edge,   & ! dlsrf_plume_dy on edge
         deta_plume_dx_edge,    & ! deta_plume_dx on edge
         deta_plume_dy_edge,    & ! deta_plume_dx on edge
         h_cavity_edge            ! cavity thickness (lsrf - topg) on edge (m)

    real(dp) :: &
         T_factor, S_factor,  & ! factors in melt-rate equations
         denom,               & ! denominator
         m1, m2,              & ! factors in relation between m and Sb
         ma, mb, mc, md,      & ! coefficients in cubic equation for m
         bmlt_float_avg         ! average value of bmlt_float in main cavity

    real(dp) :: &
         time,                & ! elapsed time during the relaxation of the plume thickness (s)            
         dt_plume,            & ! time step (s)
         my_max_dt              ! CFL-limited time step for a given cell (s) 

    integer :: i, j
    integer :: iglobal, jglobal    ! global i and j indices
    integer :: iter_melt           ! iteration counter
    integer :: ncells_sub300       ! number of cells below 300 m depth (ISOMIP+ diagnostic)

    integer :: imax, jmax          ! i and j indices of cells with extreme values

    ! max of various quantities for a given iteration
    real(dp) :: &
         Dmax, etamax, detamax, speedmax, entrainmax, bmltmax

    real(dp) :: &
         err_continuity,      & ! max plume tendency; measure of convergence of continuity equation
         err_melt               ! max difference in bmlt_float between iterations

    logical :: &
         converged_continuity,   & ! true if continuity equation has converged in all cells, else = false
         converged_melt            ! true if melt rate has converged in all cells, else = false

    real(dp), parameter :: &
!!         dt_plume_max = 900.d0,    & ! max time step for plume thickness iteration (s)
         dt_plume_max = 30.d0,    & ! max time step for plume thickness iteration (s)
                                    ! Shortened as needed to satisfy CFL
!!         time_max = scyr            ! max time (s) before giving up on convergence 
         time_max = 1000000.d0      ! max time (s) before giving up on convergence 

!!    logical, parameter :: free_surface = .false.
    logical, parameter :: free_surface = .true.

    ! prescribed MISOMIP parameters (from Table 4 of Asay-Davis et al.)
    real(dp), parameter :: &
         spec_heat_water = 3974.d0,  & ! specific heat of seawater (J/kg/deg)
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
         !WHL - Zero Coriolis to solve an easier problem
!!         f_coriolis = 0.0d0            ! Coriolis parameter (s^-1) at 75 S = 2*omega*sin(75 deg) (prescribed in text)

    ! other parameters in the plume model
    real(dp), parameter :: &    
         plume_xmax = 630000.d0,     & ! limit of the plume in the x direction (m), determined by the calving front location
         D_plume0 = 1.d0,            & ! initial plume thickness at lowest elevation, lsrf_min (m)
         D_plume_dz = 0.02d0,        & ! rate of change of initial plume thickness with increasing z (m/m)
         H0_cavity = 25.d0,          & ! cavity thickness (m) below which the entrainment gradually approaches zero
         E0 = 0.072d0                  ! entrainment coefficient (unitless)
                                       ! Bo Pederson (1980) suggests E0 = 0.072
                                       ! Jenkins (1991, JGR) uses 0.036 to compensate for lack of Coriolis in 1D model

    ! detrainment parameters
    !TODO - Increase D_plume_max to eliminate detrainment in the main cavity?
    real(dp), parameter ::  &
         D_plume_max = 100.d0,        & ! max allowed plume thickness (m) before relaxation kicks in
         tau_detrainment = 3600.d0     ! detrainment time scale (s)

    ! parameters determining convergence of iterations
    integer, parameter :: &
         maxiter_melt = 999999         ! max number of iterations of outer melt-rate loop

    real(dp), parameter :: &
         maxerr_continuity = 1.0d-6,  & ! max plume tendency (m/s) allowed for steady state
         maxerr_melt = 5.0d-2/scyr      ! max err_melt allowed for steady state (m/yr converted to m/s)

    ! relaxation parameters
    ! Value of 1 means to use the new value.  Lower values give a greater contribution from the old value.
    real(dp), parameter :: &
         relax_m = 0.1d0,    &
         relax_u = 0.1d0,    &
         relax_E = 1.0d0,    &
         relax_D = 0.1d0,    &
         relax_eta = 0.01d0

    !WHL - debug and diagnostics
    real(dp) :: solution
    integer :: count_neg, count_pos   ! no. of cells with negative and positive div(Du)

    if (main_task) then
       print*, ' '
       print*, 'In glissade_plume_melt_rate, first_call =', first_call
       print*, 'Test point: r, i, j =', rtest, itest, jtest
    endif

    !----------------------------------------------------------------
    ! Initialize some fields that are held fixed during the iteration
    !----------------------------------------------------------------      

    ! set work masks to 1 everywhere
    ! These are inputs to glissade_stagger and glissade_unstagger
    ice_mask_work(:,:) = 1
    velo_mask_work(:,:) = 1

    ! Compute the density of the ambient ocean

    rho_ambient(:,:) = eos_rho_ref * (1.d0 - eos_alpha * (T_ambient(:,:) - eos_Tref)  &
                                           + eos_beta  * (S_ambient(:,:) - eos_Sref) )

    ! Compute the pressure at the lower ice surface.
    pressure(:,:) = -rhoo*grav*lsrf(:,:)

    ! Compute a mask for where the plume thickness is computed.

    plume_mask_cell(:,:) = floating_mask(:,:)

    ! Compute plume variables for locally owned cells only.
    !TODO - Can skip if loops below are only over locally owned cells.
    do j = 1, ny
       if (j <= nhalo .or. j > ny-nhalo) then
          plume_mask_cell(:,j) = 0
       endif
    enddo

    do i = 1, nx
       if (i <= nhalo .or. i > nx-nhalo) then
          plume_mask_cell(i,:) = 0
       endif
    enddo

    ! Restrict the plume to end a few km from the calving front.
    ! Note: This may be unnecessary in CISM MISOMIP runs if the ice thickness field
    !       is smooth near the calving front.  However, the prescribed ISOMIP+ field
    !       has some strange thickness undulations near the calving front near the top
    !       and bottom domain boundaries.  Since the assumptions of the plume model 
    !       may not hold near the calving front (because of lateral mixing), it may
    !       be physically justifiable anyway to cut off the model short of the front.
    !       For now I'm using a prescribed calving limit, but a thickness criterion
    !       might work too.
    
    do j = 1, ny
       do i = 1, nx
          if  (x1(i) > plume_xmax) then
             plume_mask_cell(i,j) = 0
          endif
       enddo
    enddo

    call parallel_halo(plume_mask_cell)

    ! Mask out the plume in halo cells that lie outside the global domain.
    ! Also, identify global boundary cells for later use.
    ! Note: Ideally, we could zero out plume variables in the halo call by using an appropriate BC.  
    ! TODO: Handle plume_mask_cell with no-penetration BCs?
 
    global_bndy_west(:,:) = 0.0d0
    global_bndy_east(:,:) = 0.0d0
    global_bndy_south(:,:) = 0.0d0
    global_bndy_north(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx

          call parallel_globalindex(i, j, iglobal, jglobal)

          if (iglobal < 1 .or. iglobal > global_ewn .or. &
              jglobal < 1 .or. jglobal > global_nsn) then
             plume_mask_cell(i,j) = 0
          endif

          if (iglobal == 1) global_bndy_west(i,j) = 1
          if (iglobal == global_ewn) global_bndy_east(i,j) = 1
          if (jglobal == 1) global_bndy_south(i,j) = 1
          if (jglobal == global_nsn) global_bndy_north(i,j) = 1

       enddo
    enddo

    ! Zero out T_plume, S_plume and D_plume outside of the plume
    where (plume_mask_cell == 0)
       T_plume = 0.0d0
       S_plume = 0.0d0
       D_plume = 0.0d0
    endwhere

    ! Compute a mask on the velocity grid. Plume velocity is computed only where plume_mask_velo = 1.
    ! A vertex must be surrounded by four cells with plume_mask_cell = 1 to compute velocity.

    plume_mask_velo(:,:) = 0

    ! loop over all vertices that border locally owned cells
    do j = nhalo, ny-nhalo+1   !TODO - only to ny-nhalo?
       do i = nhalo, nx-nhalo+1
          if (plume_mask_cell(i,j+1)==1 .and. plume_mask_cell(i+1,j+1)==1 .and.  &
              plume_mask_cell(i,j)  ==1 .and. plume_mask_cell(i+1,j)  ==1) then
             plume_mask_velo(i,j) = 1
          endif
       enddo
    enddo

    call staggered_parallel_halo(plume_mask_velo)

    ! Mask out plume_mask_velo at vertices along or outside the global domain.
    ! Note: The west and east borders have iglobal indices 0 and global_ewn, respectively.
    !       The south and north borders have jglobal indices 0 and global_nsn, respectively.
    ! TODO: Handle plume_mask_velo with no-penetration BCs?
 
    do j = 1, ny-1
       do i = 1, nx-1

          call parallel_globalindex(i, j, iglobal, jglobal)

          if (iglobal <= 0 .or. iglobal >= global_ewn .or. &
              jglobal <= 0 .or. jglobal >= global_nsn) then
             plume_mask_velo(i,j) = 0
          endif

       enddo
    enddo

    ! Compute masks on cell edges, where C-grid velocities are computed.
    ! The mask is true if both adjacent cells have plume_mask_cell = 1.
    ! Note: If one cell has plume_mask_cell = 1 and the other is open water or
    !        floating ice, then the velocity will be extrapolated from a neighbor.
    !       If one cell has plume_mask_cell = 1 and the other is grounded ice
    !        or is outside the global domain, then the velocity will be set to zero.

    edge_mask_east(:,:) = 0
    edge_mask_north(:,:) = 0

    ! loop over all edges of locally owned cells
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 1) then
             edge_mask_east(i,j) = 1
          endif
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 1) then
             edge_mask_north(i,j) = 1
          endif
       enddo
    enddo

    call parallel_halo(edge_mask_east)
    call parallel_halo(edge_mask_north)

    ! Mask out edge_mask_east and edge_mask_north at edges along or outside the global domain.
    ! Note: The west and east borders have iglobal indices 0 and global_ewn, respectively.
    !       The south and north borders have jglobal indices 0 and global_nsn, respectively.
    ! TODO: Handle edge masks with no-penetration BCs?

    do j = 1, ny
       do i = 1, nx

          call parallel_globalindex(i, j, iglobal, jglobal)

          if (iglobal <= 0 .or. iglobal >= global_ewn .or. &  ! along or beyond EW boundary
              jglobal <= 0 .or. jglobal >  global_nsn) then   ! beyond NS boundary
             edge_mask_east(i,j) = 0
          endif

          if (jglobal <= 0 .or. jglobal >= global_nsn .or. &  ! along or beyond NS boundary
              iglobal <= 0 .or. iglobal >  global_ewn) then   ! beyond EW boundary
             edge_mask_north(i,j) = 0
          endif

       enddo
    enddo

    print*, ' '
    print*, 'plume_mask_cell, rank =', rtest
    do j = jtest+3, jtest-3, -1
       do i = itest-3, itest+3
          write(6,'(i8)',advance='no') plume_mask_cell(i,j)
       enddo
       write(6,*) ' '
    enddo
    
    print*, ' '
    print*, 'plume_mask_velo, rank =', rtest
    do j = jtest+3, jtest-3, -1
       write(6,'(a6)',advance='no') '      '
       do i = itest-3, itest+3
          write(6,'(i8)',advance='no') plume_mask_velo(i,j)
       enddo
       write(6,*) ' '
    enddo

    print*, ' '
    print*, 'edge_mask_north, rank =', rtest
    do j = jtest+3, jtest-3, -1
       do i = itest-3, itest+3
          write(6,'(i8)',advance='no') edge_mask_north(i,j)
       enddo
       write(6,*) ' '
    enddo
    
    print*, ' '
    print*, 'edge_mask_east, rank =', rtest
    do j = jtest+3, jtest-3, -1
       write(6,'(a6)',advance='no') '      '
       do i = itest-3, itest+3
          write(6,'(i8)',advance='no') edge_mask_east(i,j)
       enddo
       write(6,*) ' '
    enddo

    ! Compute the horizontal gradient of the lower ice surface.
    ! This is used in the computation of the pressure gradient force at velocity points.
    ! Note: With gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY combined with plume_mask_cell, 
    !       gradients are ignored at edges with a plume cell on one side and a non-plume cell on the other.
    !       This prevents spurious large gradients near the plume edge.
    ! Note: There are a couple of different ways to compute the PGF.
    !       (1) Jenkins et al. (1991) and HJH (2008) use grad(lsrf)
    !       (2) Holland & Feltham (2006) use grad(lsrf_plume) along with a density gradient.
    !       Method (1) is simpler and has the advantage that grad(lsrf) does not vary during plume evolution,
    !        making the PGF more stable (though possibly not as accurate).
    
    call glissade_centered_gradient(nx,               ny,            &
                                    dx,               dy,            &
                                    lsrf,                            &
                                    dlsrf_dx,         dlsrf_dy,      &
                                    plume_mask_cell,                 &
                                    gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY)

    !----------------------------------------------------------------
    ! Initialize some fields that are updated during the iteration.
    ! Note: T_plume, S_plume and D_plume are intent(inout) and already have initial values.
    !----------------------------------------------------------------

    if (first_call) then

       print*, ' '
       print*, 'First call: creating simple initial conditions for T_plume, S_plume and D_plume'

       ! Initialize the plume temperature and salinity.
       !WHL - For now, setting S_plume = S0 everywhere. 
       ! This means that drho_plume = rho_ambient - rho_plume will decrease in the upslope direction.
       where (plume_mask_cell == 1)
          T_plume = T_ambient
          S_plume = S0
       endwhere

       ! Initialize the plume thickness.
       ! This is tricky. Since entrainment is non-negative and is equal to div*(Du),
       !  we ideally want div*(Du) > 0 in most cells. If the flow is upslope and D increases
       !  upslope, then div*(Du) will generally be positive in most of the domain.
       ! D_plume is constrained not to be thicker than the sub-shelf ocean cavity at
       !  the start of the simulation. If convergence results in D_plume > H_cavity,
       !  there will be a pressure gradient force tending to reduce D_plume.
       ! Note: Units for D_plume here are meters.
       !       Would have to change if D_plume is input in scaled model units

       lsrf_min = minval(lsrf)
       lsrf_min = parallel_reduce_min(lsrf_min)
       !WHL - Use an absolute level instead of the global min?

       do j = 1, ny
          do i = 1, nx
             if (plume_mask_cell(i,j) == 1) then
                D_plume(i,j) = D_plume0 + D_plume_dz * (lsrf(i,j) - lsrf_min)
                H_cavity = lsrf(i,j) - topg(i,j)
                D_plume(i,j) = min(D_plume(i,j), H_cavity)
             endif
          enddo
       enddo

    else

       print*, ' '
       print*, 'Using input values of T_plume, S_plume and D_plume'

    endif   ! first_call

    if (verbose_melt) then

       print*, ' '
       print*, 'T_plume, rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f8.3)',advance='no') T_plume(i,j)
          enddo
          write(6,*) ' '
       enddo
       
       print*, ' '
       print*, 'S_plume, rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f8.3)',advance='no') S_plume(i,j)
          enddo
          write(6,*) ' '
       enddo
       
    endif

    if (verbose_continuity) then

       print*, ' '
       print*, 'D_plume (m), rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f12.5)',advance='no') D_plume(i,j)
          enddo
          write(6,*) ' '
       enddo

    endif

    ! Initialize T and S at the base of the ice.
    ! Start with the same salinity as the underlying water, with T at the freezing point
    where (plume_mask_cell == 1)
       S_basal = S_plume
       T_basal = lambda1*S_plume + lambda2 + lambda3*pressure
    elsewhere
       S_basal = S_ambient
       T_basal = lambda1*S_ambient + lambda2 + lambda3*pressure
    endwhere

    ! Initialize the surface displacement eta_plume = max(D_plume - H_cavity, 0),
    !  where H_cavity = lsrf - topg.
    ! eta_plume = 0 in most of the domain, but eta_plume > 0 in shallow regions
    !  where the flow is convergent and the plume extends down to the seafloor.
    ! The result is a strong outflow from cells where eta > 0.

    where (plume_mask_cell == 1)
       eta_plume = max(D_plume - (lsrf - topg), 0.0d0)
    elsewhere
       eta_plume = 0.0d0
    endwhere

    ! Initialize other fields
    u_plume(:,:) = 0.0d0
    v_plume(:,:) = 0.0d0
    plume_speed(:,:) = 0.0d0

    u_plume_east(:,:) = 0.0d0
    v_plume_east(:,:) = 0.0d0
    plume_speed_east(:,:) = 0.0d0

    u_plume_north(:,:) = 0.0d0
    v_plume_north(:,:) = 0.0d0
    plume_speed_north(:,:) = 0.0d0

    ustar_plume(:,:) = 0.0d0
    divDu_plume(:,:) = 0.0d0
    entrainment(:,:) = 0.0d0
    detrainment(:,:) = 0.0d0
    bmlt_float(:,:) = 0.0d0

    !WHL - debug
    if (verbose_melt) then

       print*, ' '
       print*, 'T_ambient (deg C), rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f8.3)',advance='no') T_ambient(i,j)
          enddo
          write(6,*) ' '
       enddo
       
       print*, ' '
       print*, 'T_ambient - T_plume, rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f8.3)',advance='no') T_ambient(i,j) - T_plume(i,j)
          enddo
          write(6,*) ' '
       enddo
       
       print*, ' '
       print*, 'S_ambient (psu), rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f8.3)',advance='no') S_ambient(i,j)
          enddo
          write(6,*) ' '
       enddo
       
       print*, ' '
       print*, 'S_ambient - S_plume, rank =', rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f8.3)',advance='no') S_ambient(i,j) - S_plume(i,j)
          enddo
          write(6,*) ' '
       enddo
       
    endif   ! verbose melt

    !WHL - debug
    if (main_task) then
       print*, ' '
       print*, 'Start melt-rate iteration'
    endif

    time = 0.0d0

    !--------------------------------------------------------------------
    ! Iterate to compute the melt rate at the ice-ocean interface
    ! The solution method is:
    ! (1) Given D, eta, and rho for the plume, compute u and v from the momentum balance.
    ! (2) Given u and v, compute the entrainment, and then compute D from continuity.
    ! (3) Given the entrainment and friction velocity, compute the basal melt rate
    !     and new values of T_b, S_b, T_p and S_p.
    ! Iterate to convergence.
    !--------------------------------------------------------------------

    do iter_melt = 1, maxiter_melt   ! outer melt-rate iteration

       if (main_task) then 
          print*, ' '
          print*, 'iter_melt =', iter_melt
       endif

       ! halo updates
       call parallel_halo(T_plume)
       call parallel_halo(S_plume)
       call parallel_halo(D_plume)

       ! save variables from previous iteration
       entrainment_old(:,:) = entrainment(:,:)
       divDu_plume_old(:,:) = divDu_plume(:,:)
       D_plume_old(:,:) = D_plume(:,:)
       eta_plume_old(:,:) = eta_plume(:,:)
       S_plume_old(:,:) = S_plume(:,:)
       T_plume_old(:,:) = T_plume(:,:)
       S_basal_old(:,:) = S_basal(:,:)
       T_basal_old(:,:) = T_basal(:,:)
       bmlt_float_old(:,:) = bmlt_float(:,:)

       u_plume_east_old(:,:) = u_plume_east(:,:)
       v_plume_east_old(:,:) = v_plume_east(:,:)
       u_plume_north_old(:,:) = u_plume_north(:,:)
       v_plume_north_old(:,:) = v_plume_north(:,:)

       ! Compute the plume density
       rho_plume(:,:) = eos_rho_ref * (1.d0 - eos_alpha * (T_plume(:,:) - eos_Tref)  &
                                            + eos_beta  * (S_plume(:,:) - eos_Sref) )

       ! Note: With gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY combined with plume_mask_cell, 
       !       gradients are ignored at edges with a plume cell on one side and a non-plume cell on the other.
       !       This prevents spurious large gradients near the plume edge.

       ! Find the density difference between the ambient ocean and the plume.
       ! Compute on the ice grid, then interpolate to the velocity grid.

       drho_plume(:,:) = rho_ambient(:,:) - rho_plume(:,:)

       ! TODO - What to do where drho_plume < 0? Set plume_mask_cell = 0?
       ! TODO = Print where drho_plume < 0.
 
       ! Compute the elevation of the plume-ambient interface
       ! Note: lsrf and lsrf_plume are negative below sea level, and D_plume >= 0 by definition

       lsrf_plume(:,:) = lsrf(:,:) - D_plume(:,:)

       ! Compute the gradient of the plume-ambient interface.
       ! This gradient appears in the entrainment term.
       ! Note: With gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY combined with plume_mask_cell, 
       !       gradients are ignored at edges with a plume cell on one side and a non-plume cell on the other.
       !       This prevents spurious large gradients near the plume edge.
       
       call glissade_centered_gradient(nx,               ny,             &
                                       dx,               dy,             &
                                       lsrf_plume,                       &
                                       dlsrf_plume_dx,   dlsrf_plume_dy, &
                                       plume_mask_cell,                  &
                                       gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY)

       ! Compute the gradient of eta.
       ! This gradient appears in the free-surface part of the pressure gradient force.

       call glissade_centered_gradient(nx,               ny,             &
                                       dx,               dy,             &
                                       eta_plume,                        &
                                       deta_plume_dx,    deta_plume_dy,  &
                                       plume_mask_cell,                  &
                                       gradient_margin_in = HO_GRADIENT_MARGIN_ICE_ONLY)

       if (free_surface) then
          print*, ' '
          print*, 'Free surface calculation is ON'
       else
          print*, ' '
          print*, 'Free surface calculation is OFF'
       endif   ! free_surface

       if (verbose_melt) then

          print*, ' '
          print*, 'T_plume, rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f8.3)',advance='no') T_plume(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'S_plume, rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f8.3)',advance='no') S_plume(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'drho_plume (kg/m^3), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') drho_plume(i,j)
             enddo
             write(6,*) ' '
          enddo
          
       endif

       if (verbose_velo) then

          print*, ' '
          print*, 'dlsrf_dx, rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(e12.3)',advance='no') dlsrf_dx(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'dlsrf_dy, rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(e12.3)',advance='no') dlsrf_dy(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'eta_plume, rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(e12.3)',advance='no') eta_plume(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'deta_plume_dx, rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(e12.3)',advance='no') deta_plume_dx(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'deta_plume_dy, rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(e12.3)',advance='no') deta_plume_dy(i,j)
             enddo
             write(6,*) ' '
          enddo
          
       endif

       if (this_rank == rtest) then
          print*, ' '
          print*, 'rtest, itest, jtest:', rtest, itest, jtest
       endif

       ! compute u and v on east cell edges
       if (main_task) then
          print*, ' '
          print*, 'compute east edge velocities'
       endif

       do j = nhalo, ny-nhalo
          do i = nhalo, nx-nhalo
             if (edge_mask_east(i,j) == 1) then

                ! Compute horizontal pressure gradient force, not including the free-surface term.
                ! Based on HJH 2008
                ! On east edges, the x derivatives are based on values in the two adjacent cells,
                !  to preserve the advantages of a C grid.
                ! The y derivatives are averaged from the neighboring vertices.
                
                D_plume_edge = (D_plume(i,j) + D_plume(i+1,j)) / 2.0d0
                grav_reduced_edge = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i+1,j)) / 2.0d0

                dlsrf_dx_edge = (lsrf(i+1,j) - lsrf(i,j)) / dx
                dlsrf_dy_edge = (dlsrf_dy(i,j) + dlsrf_dy(i,j-1)) / 2.0d0

                pgf_x = grav_reduced_edge * D_plume_edge * dlsrf_dx_edge
                pgf_y = grav_reduced_edge * D_plume_edge * dlsrf_dy_edge

                !WHL - debug
                if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, 'xterm 1, yterm 1:', grav_reduced_edge * D_plume_edge * dlsrf_dx_edge, &
                                                grav_reduced_edge * D_plume_edge * dlsrf_dy_edge
                endif

                ! Optionally, add the free-surface term
                if (free_surface) then

                   deta_plume_dx_edge = (eta_plume(i+1,j) - eta_plume(i,j)) / dx
                   deta_plume_dy_edge = (deta_plume_dy(i,j) + deta_plume_dy(i,j-1)) / 2.0d0

                   if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                      print*, 'xterm 2, yterm 2:', -grav * D_plume_edge * deta_plume_dx_edge, &
                                                   -grav * D_plume_edge * deta_plume_dy_edge
                   endif

                   pgf_x = pgf_x - grav*D_plume_edge*deta_plume_dx_edge
                   pgf_y = pgf_y - grav*D_plume_edge*deta_plume_dy_edge

                endif

                !TODO - remove plume_speed as local variable; update below
                call compute_local_velocity(&
                     i, j,                    &   ! diagnostic only
                     itest, jtest, rtest,     &   ! diagnostic only
                     u_tidal,                 &
                     c_drag,                  &
                     f_coriolis,              &
                     D_plume_edge,            &
                     pgf_x,                   &
                     pgf_y,                   &
                     u_plume_east(i,j),       &
                     v_plume_east(i,j),       &
                     plume_speed_east(i,j))

             endif   ! edge_mask_east
          enddo  ! i
       enddo  ! j

       ! relax toward the velocity just computed
       u_plume_east(:,:) = (1.0d0 - relax_u)*u_plume_east_old(:,:) + relax_u*u_plume_east(:,:)
       v_plume_east(:,:) = (1.0d0 - relax_u)*v_plume_east_old(:,:) + relax_u*v_plume_east(:,:)

       ! compute u and v on north cell edges
       if (main_task) then
          print*, ' '
          print*, 'compute north edge velocities'
       endif

       do j = nhalo, ny-nhalo
          do i = nhalo, nx-nhalo
             if (edge_mask_north(i,j) == 1) then

                ! Compute horizontal pressure gradient force, not including the free-surface term
                ! Based on HJH 2008
                ! On north edges, the y derivatives are based on values in the two adjacent cells,
                !  to preserve the advantages of a C grid.
                ! The x derivatives are averaged from the neighboring vertices.
                
                D_plume_edge = (D_plume(i,j) + D_plume(i,j+1)) / 2.0d0
                grav_reduced_edge = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i,j+1)) / 2.0d0

                dlsrf_dx_edge = (dlsrf_dx(i,j) + dlsrf_dx(i-1,j)) / 2.0d0
                dlsrf_dy_edge = (lsrf(i,j+1) - lsrf(i,j)) / dy

                pgf_x = grav_reduced_edge * D_plume_edge * dlsrf_dx_edge
                pgf_y = grav_reduced_edge * D_plume_edge * dlsrf_dy_edge

                !WHL - debug
                if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, 'xterm 1, yterm 1:', grav_reduced_edge * D_plume_edge * dlsrf_dx_edge, &
                                                grav_reduced_edge * D_plume_edge * dlsrf_dy_edge
                endif

                ! Optionally, add the free-surface term

                if (free_surface) then

                   deta_plume_dx_edge = (deta_plume_dx(i,j) + deta_plume_dx(i-1,j)) / 2.0d0
                   deta_plume_dy_edge = (eta_plume(i,j+1) - eta_plume(i,j)) / dy

                   !WHL - debug
                   if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                      print*, 'xterm 2, yterm 2:', -grav * D_plume_edge * deta_plume_dx_edge, &
                                                   -grav * D_plume_edge * deta_plume_dy_edge
                   endif

                   pgf_x = pgf_x - grav*D_plume_edge*deta_plume_dx_edge
                   pgf_y = pgf_y - grav*D_plume_edge*deta_plume_dy_edge

                endif

                call compute_local_velocity(&
                     i, j,                    &   ! diagnostic only
                     itest, jtest, rtest,     &   ! diagnostic only
                     u_tidal,                 &
                     c_drag,                  &
                     f_coriolis,              &
                     D_plume_edge,            &
                     pgf_x,                   &
                     pgf_y,                   &
                     u_plume_north(i,j),      &
                     v_plume_north(i,j),      &
                     plume_speed_north(i,j))

             endif   ! edge_mask_north
          enddo  ! i
       enddo  ! j

       ! relax toward the velocity just computed
       u_plume_north(:,:) = (1.0d0 - relax_u)*u_plume_north_old(:,:) + relax_u*u_plume_north(:,:)
       v_plume_north(:,:) = (1.0d0 - relax_u)*v_plume_north_old(:,:) + relax_u*v_plume_north(:,:)

       ! Extrapolate the velocity to open edges (plume on one side, open water on the other)
       !  This extrapolation is not expected to be accurate, but it prevents large convergence
       !  in cells adjacent to water.
       ! If the plume exists on neither side of the edge, the velocity remains set to zero.
       ! Also, u_plume_east = 0 on global E and W boundaries, and v_plume_north = 0 on global N and S boundaries.
       !  This prevents outflow through domain walls.
       !  Along the upper ("northern") boundary of the ISOMIP+ domain, the flow is forced to form an eastward jet. 

       do j = nhalo, ny-nhalo
          do i = nhalo, nx-nhalo

             ! east edges
             if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 0 .and. global_bndy_east(i,j) == 0) then
                if (lsrf(i+1,j) == 0.0d0 .or. floating_mask(i+1,j) == 1) then
                   ! water in cell (i+1,j); get plume velocity from edge (i-1,j)
                   u_plume_east(i,j) = u_plume_east(i-1,j)
                endif
             elseif (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i+1,j) == 1 .and. global_bndy_west(i+1,j) == 0) then
                if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                   ! water in cell (i,j); get plume velocity from edge (i+1,j)
                   u_plume_east(i,j) = u_plume_east(i+1,j)
                endif
             endif

             ! north edges
             if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 0 .and. global_bndy_north(i,j) == 0) then
                if (lsrf(i,j+1) == 0.0d0 .or. floating_mask(i,j+1) == 1) then
                   ! water in cell (i,j+1); get plume velocity from edge (i,j-1)
                   v_plume_north(i,j) = v_plume_north(i,j-1)
                endif
             elseif (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i,j+1) == 1 .and. global_bndy_south(i,j+1) == 0) then
                if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                   ! water in cell (i,j); get plume velocity from edge (i,j+1)
                   v_plume_north(i,j) = v_plume_north(i,j+1)
                endif
             endif
             
          enddo   ! i
       enddo   ! j

       ! update the plume velocity on edges
       plume_speed_east(:,:) = sqrt(u_plume_east(:,:)**2 + v_plume_east(:,:)**2 + u_tidal**2)
       plume_speed_north(:,:) = sqrt(u_plume_north(:,:)**2 + v_plume_north(:,:)**2 + u_tidal**2)

       ! Interpolate the plume velocity to cell centers, and compute the friction velocity.

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             u_plume(i,j) = (u_plume_east(i,j) + u_plume_east(i-1,j)) / 2.0d0
             v_plume(i,j) = (v_plume_north(i,j) + v_plume_north(i,j-1)) / 2.0d0
             plume_speed(i,j) = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)
             ustar_plume(i,j) = sqrt(c_drag) * plume_speed(i,j)
          enddo
       enddo

!!       ! Compute the friction velocity on the velocity grid
!!       stag_ustar_plume(:,:) = sqrt(c_drag) * plume_speed(:,:)

       ! Interpolate to the ice grid
!!       call glissade_unstagger(&
!!            nx,                     &
!!            ny,                     &
!!            stag_ustar_plume,       &
!!            ustar_plume,            &
!!            velo_mask_work,         &
!!            stagger_margin_in = 0)

       !--------------------------------------------------------------------
       ! Compute entrainment as a function of the plume speed and the slope of the
       !  plume-ambient interface, following Bo Pederson (1980) and Jenkins (1991).
       ! Entrainment is computed at cell edges (where the slope is computed
       !  most naturally) and then interpolated to the ice grid.
       !--------------------------------------------------------------------

       entrainment_east(:,:) = 0.0d0

       do j = nhalo, ny-nhalo
          do i = nhalo, nx-nhalo

             D_plume_edge = (D_plume(i,j) + D_plume(i+1,j)) / 2.0d0

             dlsrf_plume_dx_edge = (lsrf(i+1,j) - lsrf(i,j)) / dx
             dlsrf_plume_dy_edge = (dlsrf_plume_dy(i,j) + dlsrf_plume_dy(i,j-1)) / 2.0d0
             slope = sqrt(dlsrf_plume_dx(i,j)**2 + dlsrf_plume_dy(i,j)**2)
             theta_slope = atan(slope)

             ! Note: Here the cavity thickness is defined as the thickness between the lower plume surface
             !       and the topograpy.  Entrainment -> 0 as h_cavity_edge -> 0
             h_cavity_edge = ( (lsrf_plume(i,j) - topg(i,j)) + (lsrf_plume(i+1,j) - topg(i+1,j)) ) / 2.0d0
             h_cavity_edge = max(h_cavity_edge, 0.0d0)

             entrainment_east(i,j) = E0 * plume_speed_east(i,j) * sin(theta_slope) * tanh(h_cavity_edge/H0_cavity)

          enddo
       enddo

       ! compute entrainment on north edges

       entrainment_north(:,:) = 0.0d0

       do j = nhalo, ny-nhalo
          do i = nhalo, nx-nhalo

             D_plume_edge = (D_plume(i,j) + D_plume(i,j+1)) / 2.0d0

             dlsrf_plume_dx_edge = (dlsrf_plume_dx(i,j) + dlsrf_plume_dx(i-1,j)) / 2.0d0
             dlsrf_plume_dy_edge = (lsrf(i,j+1) - lsrf(i,j)) / dy
             slope = sqrt(dlsrf_plume_dx_edge**2 + dlsrf_plume_dy_edge**2)
             theta_slope = atan(slope)

             h_cavity_edge = ( (lsrf_plume(i,j) - topg(i,j)) + (lsrf_plume(i,j+1) - topg(i,j+1)) ) / 2.0d0
             h_cavity_edge = max(h_cavity_edge, 0.0d0)

             entrainment_north(i,j) = E0 * plume_speed_north(i,j) * sin(theta_slope) * tanh(h_cavity_edge/H0_cavity)

          enddo
       enddo

       ! interpolate entrainment from edges to cell centers

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo

             entrainment(i,j) = ( entrainment_east(i-1,j) + entrainment_east(i,j) +   &
                                  entrainment_north(i,j-1) + entrainment_north(i,j) ) / 4.0d0 
          enddo
       enddo

       !WHL - Relaxing.  May not be needed for entrainment, if already done for velocity.
       entrainment(:,:) = (1.0d0 - relax_E)*entrainment_old(:,:) + relax_E*entrainment(:,:)

       ! Compute detrainment.
       ! This is not a physically based mechanism, just a regularization to prevent
       !  very thick plumes in regions where the flow is not well behaved.
       ! Ideally, detrainment = 0 almost everywhere.

       if (free_surface) then
          where (D_plume > D_plume_max)
             detrainment = (D_plume - D_plume_max) / tau_detrainment
          endwhere
       else    ! not a free surface; no grad(eta) term in PGF
          where (eta_plume > 0.0d0)
             detrainment(:,:) = eta_plume(:,:) / tau_detrainment
          endwhere
       endif

       !--------------------------------------------------------------------
       ! Compute div*(Du) on the ice grid (baroclinic component of divergence).
       ! Use the upstream value for D. Thus as D increases, the outflow from
       !  a cell increases but the inflow does not.
       !--------------------------------------------------------------------

       !WHL - debug
       count_neg = 0
       count_pos = 0

       ! loop over locally owned cells
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (plume_mask_cell(i,j) == 1) then

                ! d/dx(Du)

                if (u_plume_east(i,j) > 0.0d0) then
                   D_plume_east = D_plume(i,j)
                else
                   D_plume_east = D_plume(i+1,j)
                endif
                
                if (u_plume_east(i-1,j) > 0.0d0) then
                   D_plume_west = D_plume(i-1,j)
                else
                   D_plume_west = D_plume(i,j)
                endif

                dDu_dx = (D_plume_east*u_plume_east(i,j) - D_plume_west*u_plume_east(i-1,j)) / dx

                ! d/dy(Dv)

                if (v_plume_north(i,j) > 0.0d0) then
                   D_plume_north = D_plume(i,j)
                else
                   D_plume_north = D_plume(i,j+1)
                endif

                if (v_plume_north(i,j-1) > 0.0d0) then
                   D_plume_south = D_plume(i,j-1)
                else
                   D_plume_south = D_plume(i,j)
                endif
                
                dDv_dy = (D_plume_north*v_plume_north(i,j) - D_plume_south*v_plume_north(i,j-1)) / dy
                   
                divDu_plume(i,j) = dDu_dx + dDv_dy

                if (verbose_continuity .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'old D_plume =', D_plume(i,j)
                   print*, 'D_west, D_east:', D_plume_west, D_plume_east
                   print*, 'D_north, D_south:', D_plume_north, D_plume_south
                   print*, 'u_west, u_east:', u_plume_east(i-1,j), u_plume_east(i,j)
                   print*, 'v_north, v_south:', v_plume_north(i,j-1), v_plume_north(i,j)
                   print*, 'dDu_dx, dDv_dy:', dDu_dx, dDv_dy
                   print*, 'div(Du) (m/s):', divDu_plume(i,j)
                endif
                
                !WHL - debug
                if (divDu_plume(i,j) < 0.0d0) then
                   count_neg = count_neg + 1
                else
                   count_pos = count_pos + 1
                endif
                
             endif  ! plume_mask_cell
          enddo   ! i
       enddo   ! j
       
       if (verbose_continuity) then
          print*, ' '
          print*, 'count_neg, count_pos:', count_neg, count_pos
       endif

       !WHL - Relaxing
       divDu_plume(:,:) = (1.0d0 - relax_E)*divDu_plume_old(:,:) + relax_E*divDu_plume(:,:)

       ! Determine the time step for this iteration.
       ! Must satisfy a CFL condition. 

       dt_plume = dt_plume_max
       imax = 1
       jmax = 1
       
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (plume_mask_cell(i,j) == 1) then
                my_max_dt = 0.1d0*dx / max( abs(u_plume_east(i,j)),  abs(u_plume_east(i-1,j)), &
                                            abs(v_plume_north(i,j)), abs(v_plume_north(i,j-1)) )
                if (my_max_dt < dt_plume) then
                   dt_plume = my_max_dt
                   imax = i
                   jmax = j
                endif
             endif
          enddo
       enddo
       
       !--------------------------------------------------------------------
       ! Adjust D_plume to satisfy the steady-state continuity equation:
       !
       !    del*(Du) = e - d
       !
       ! where e is entrainment and d is detrainment.
       !  
       ! This is a work in progress. The best method might be a matrix solve.
       ! For now, try relaxing to the solution using
       ! 
       !       dD/dt = (e - d) - del*(Du), where u is the full plume velocity
       !
       !--------------------------------------------------------------------

       ! Increment the time and the plume thickness
       time = time + dt_plume

       !WHL - debug
       if (main_task) then
          print*, ' '
          print*, 'Solve continuity equation'
          print*, 'Current dt (s) =', dt_plume
          print*, 'Total time (s) =', time
       endif
       
       err_continuity = 0.0d0

       ! Increment the plume thickness
       ! loop over locally owned cells
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (plume_mask_cell(i,j) == 1) then
                
                plume_tendency = entrainment(i,j) - detrainment(i,j) - divDu_plume(i,j)
                dD_plume(i,j) = plume_tendency*dt_plume
                D_plume(i,j) = D_plume(i,j) + dD_plume(i,j)
                
                if (this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'i, j, oldD, dD, new D:', i, j, D_plume_old(i,j), dD_plume(i,j), D_plume(i,j)
                endif

                ! Relax toward new D_plume
                D_plume(i,j) = (1.0d0 - relax_D)*D_plume_old(i,j) + relax_D*D_plume(i,j)

                if (this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, 'relaxed D:', D_plume(i,j)
                endif

                ! This should not happen with an adaptive time step
                if (D_plume(i,j) < 0.0d0) then
                   print*, 'ERROR: Exceeded CFL for plume adjustment:', i, j
                   print*, 'rank, i, j, D_plume, correction:', this_rank, i, j, D_plume(i,j), dD_plume(i,j)
                   stop
                endif
                
                ! Compute the new value of eta_plume, given the new (relaxed) D_plume
                eta_plume_unrelaxed = max(D_plume(i,j) - (lsrf(i,j) - topg(i,j)), 0.0d0)

                ! Relax toward eta_plume
                ! May need slow adjustment to suppress gravity waves
                eta_plume(i,j) = (1.0d0 - relax_eta)*eta_plume_old(i,j) + relax_eta*eta_plume_unrelaxed
                
                ! Correct D_plume to account for relaxation of eta
                D_plume(i,j) = D_plume(i,j) + (eta_plume(i,j) - eta_plume_unrelaxed)

                if (this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, 'unrelaxed eta:', eta_plume_unrelaxed
                   print*, 'relaxed eta:', eta_plume(i,j)
                   print*, 'eta correction:', eta_plume_unrelaxed - eta_plume(i,j)
                   print*, 'Corrected D_plume:', D_plume(i,j)
                endif

                if (abs(plume_tendency) > err_continuity) then
                   err_continuity = abs(plume_tendency)
                   imax = i
                   jmax = j
                endif
                
             endif  ! plume_mask_cell
          enddo   ! i
       enddo   ! j
       
       !TODO - Add global maxval for err_continuity

       if (verbose_continuity) then

          print*, ' '
          print*, 'v_plume_north (m/s), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') v_plume_north(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'u_plume_east (m/s), rank =', rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(a6)',advance='no') '      '
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') u_plume_east(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'plume_speed (m/s), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(e12.3)',advance='no') plume_speed(i,j)
             enddo
             write(6,*) ' '
          enddo

!          print*, ' '
!          print*, 'detrainment (m/s), rank =', rtest
!          do j = jtest+3, jtest-3, -1
!             do i = itest-3, itest+3
!                write(6,'(e12.3)',advance='no') detrainment(i,j)
!             enddo
!             write(6,*) ' '
!          enddo

          print*, ' '
          print*, 'entrainment (m/s), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(e12.3)',advance='no') entrainment(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'divDu_plume (m/s), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(e12.3)',advance='no') divDu_plume(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'dD_plume (m), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') dD_plume(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'New D_plume (m), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') D_plume(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'New eta_plume (m), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') eta_plume(i,j)
             enddo
             write(6,*) ' '
          enddo

       endif  ! verbose_continuity

       ! print location of max continuity error
       print*, ' '
       print*, 'i, j, D_plume, dD, max tendency (m/s):', imax, jmax, D_plume(imax,jmax), dD_plume(imax,jmax), err_continuity

       !WHL - debug - Find location of max plume speed
       ! loop over locally owned cells
       speedmax = 0.0d0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (plume_speed(i,j) > speedmax) then
                speedmax = plume_speed(i,j)
                imax = i
                jmax = j
             endif
          enddo
       enddo
       print*, 'i, j, max plume_speed (m/s):', imax, jmax, plume_speed(imax,jmax)

       !WHL - debug - Find location of max entrainment rate
       ! loop over locally owned cells
       entrainmax = 0.0d0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (entrainment(i,j) > entrainmax) then
                entrainmax = entrainment(i,j)
                imax = i
                jmax = j
             endif
          enddo
       enddo
       print*, 'i, j, max entrainment (m/yr):', imax, jmax, entrainment(imax,jmax)*scyr

       !WHL - debug - Find location of max plume thickness
       ! loop over locally owned cells
       Dmax = 0.0d0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (D_plume(i,j) > Dmax) then
                Dmax = D_plume(i,j)
                imax = i
                jmax = j
             endif
          enddo
       enddo
       print*, 'i, j, max D_plume (m):', imax, jmax, D_plume(imax,jmax)

       !WHL - debug - Find location of max eta
       ! loop over locally owned cells
       etamax = 0.0d0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (eta_plume(i,j) > etamax) then
                etamax = eta_plume(i,j)
                imax = i
                jmax = j
             endif
          enddo
       enddo
       print*, 'i, j, max(eta_plume):', imax, jmax, eta_plume(imax,jmax)

       !WHL - debug - Find location of max change in eta
       ! loop over locally owned cells
       detamax = 0.0d0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (abs(eta_plume(i,j) - eta_plume_old(i,j)) > detamax) then
                detamax = abs(eta_plume(i,j) - eta_plume_old(i,j))
                imax = i
                jmax = j
             endif
          enddo
       enddo
       print*, 'i, j, old eta, new eta, d_eta:', imax, jmax, eta_plume_old(imax,jmax), &
            eta_plume(imax,jmax), eta_plume(imax,jmax) - eta_plume_old(imax,jmax)

       ! check convergence of the continuity equation in all cells
       !TODO - Not sure this is needed.  Maybe just need to go back and check error in momentum equation with new D.
       !TODO - Add global sum here for parallel runs

       if (err_continuity > maxerr_continuity) then  ! not converged
          print*, ' '
          print*, 'Continuity loop calculation not converged: iter, time, max tendency (m/s) =', &
               iter_melt, time, err_continuity
          converged_continuity = .false.
          
       else   ! converged; exit the continuity loop

          print*, ' '
          print*, 'Continuity calculation CONVERGED, iter, time, max tendency (m/s) =', &
               iter_melt, time, err_continuity
          converged_continuity = .true.

       endif

       !TODO - Put the melt rate calculation in a subroutine?
       !--------------------------------------------------------------------
       ! Compute the melt rate at the ice-ocean interface.
       !
       ! There are 5 equations for 5 unknowns: m, Tb, Sb, T and S
       ! where m = melt rate at ice-ocean interface
       !       Tb = potential temperature at ice-ocean interface
       !       Sb = salinity at ice-ocean interface
       !       T  = potential temperature of boundary-layer plume
       !       S  = salinity of boundary-layer plume
       ! 
       ! (1) rhow * m * L  = rhoo * cw * u_fric * gammaT * (T - Tb)
       ! (2) rhow * m * Sb = rhoo * u_fric * gammaS *(S - Sb)
       ! (3) Tb = lambda1*Sb + lambda2 + lambda3*pb 
       ! (4) L * m = -cw * e * (T - Ta)
       ! (5) S * m = -e * (S - Sa)
       !
       ! Eq. 1 and 2 describe heat and salt transfer at the ice-ocean interface.
       ! Eq. 3 is the linearized liquidus relation that determines the potential freezing point.
       ! Eq. 4 and 5 describe heat and salt entrainment from the ambient ocean to the boundary-layer plume,
       !  where Ta and Sa are the potential temperature and salinity of the ambient ocean.
       !
       ! We can rewrite (1) and (2) as
       !
       ! (1)     m = T_factor * (T - Tb)
       ! (2)  Sb*m = S_factor * (S - Sb) 
       !
       ! where T_factor = (rhoo * cw * ufric * gammaT) / (rhow * L)
       !       S_factor = (rhoo * ufric * gammaS) / rhow
       !
       ! Rearrange (4):  T = Ta - (L/(cw*e)) * m
       ! 
       ! Use (3) and (4) to replace T and Tb in (1):
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

       !WHL - debug -  Test cubic solver
!       ma =    2.d0
!       mb =  -30.d0
!       mc =  162.d0
!       md = -350.d0
!       call cubic_solver(ma, mb, mc, md, solution)
!       print*, 'Trial cubic solution =', solution
!       print*, 'True solution =', (10.d0 + sqrt(108.d0))**(1.d0/3.d0) - (-10.d0 + sqrt(108.d0))**(1.d0/3.d0) + 5.d0

       ! Loop over locally owned cells
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             
             if (plume_mask_cell(i,j) == 1 .and. entrainment(i,j) > 0.0d0) then

                T_factor = (rhoo * spec_heat_water * ustar_plume(i,j) * gammaT) / (rhow * lhci)
                S_factor = (rhoo * ustar_plume(i,j) * gammaS) / rhow
                
                denom = 1.d0 + (T_factor*lhci)/(spec_heat_water*entrainment(i,j))
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

                if (verbose_melt .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Melt rate calc: rank, i, j =', rtest, i, j
                   print*, 'lsrf (m), pressure (Pa) =', lsrf(i,j), pressure(i,j)
                   print*, 'T_factor (m/s/deg), S_factor (m/s)=', T_factor, S_factor
                   print*, 'entrainment (m/s) =', entrainment(i,j)
                   print*, 'm1 (m/s/psu) =', m1
                   print*, 'm2 (m/s) =', m2
                   print*, 'denom =', denom
                   print*, 'a, b, c, d =', ma, mb, mc, md
                   print*, 'residual of cubic solve =', ma*bmlt_float(i,j)**3 + mb*bmlt_float(i,j)**2 + mc*bmlt_float(i,j) + md
                endif
                
                ! Given the melt rate, compute Sb and Tb
!                S_basal(i,j) = (S_factor * entrainment(i,j) * S_ambient(i,j)) /  &
!                               ( (bmlt_float(i,j) + S_factor) * (bmlt_float(i,j) + entrainment(i,j)) )
                S_basal(i,j) = (bmlt_float(i,j) - m2) / m1
                T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)

                ! Given m, compute S and T for plume
                T_plume(i,j) = T_ambient(i,j) - (lhci/(spec_heat_water*entrainment(i,j))) * bmlt_float(i,j)
                S_plume(i,j) = S_ambient(i,j) * entrainment(i,j) / (bmlt_float(i,j) + entrainment(i,j))

                ! Relax toward solution
                S_plume(i,j) = (1.0d0 - relax_m)*S_plume_old(i,j) + relax_m*S_plume(i,j)
                T_plume(i,j) = (1.0d0 - relax_m)*T_plume_old(i,j) + relax_m*T_plume(i,j)
                S_basal(i,j) = (1.0d0 - relax_m)*S_basal_old(i,j) + relax_m*S_basal(i,j)
                T_basal(i,j) = (1.0d0 - relax_m)*T_basal_old(i,j) + relax_m*T_basal(i,j)
                !TODO - Use freezing relation instead of relaxation parameter?
!!                T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)

                !WHL - debug - check for NaNs
                if (T_plume(i,j) /= T_plume(i,j) .or. S_plume(i,j) /= S_plume(i,j) .or. &
                    T_basal(i,j) /= T_basal(i,j) .or. S_basal(i,j) /= S_basal(i,j) .or. &
                    bmlt_float(i,j) /= bmlt_float(i,j)) then
                   print*, 'Bad values, i, j =', i, j
                   print*, 'T_plume, S_plume:', T_plume(i,j), S_plume(i,j)
                   print*, 'T_basal, S_basal:', T_basal(i,j), S_basal(i,j)
                   print*, 'bmlt_float:', bmlt_float(i,j)
                   stop
                endif

                if (this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Computed melt, rank, i, j =', this_rank, i, j
                   print*, 'Entrainment, bmlt_float (m/yr) =', i, j, entrainment(i,j)*scyr, bmlt_float(i,j)*scyr
                   print*, 'T_b, S_b =', T_basal(i,j), S_basal(i,j)
                   print*, 'T_p, S_p =', T_plume(i,j), S_plume(i,j)
                endif

             else    ! plume_mask_cell = 0

                bmlt_float(i,j) = 0.0d0

                S_plume(i,j) = S_ambient(i,j)
                T_plume(i,j) = T_ambient(i,j)

                S_basal(i,j) = S_ambient(i,j)
                T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
                   
             endif   ! plume_mask_cell and entrainment > 0

          enddo   ! i
       enddo   ! j
          
       !WHL - debug - Find location of max melt rate
       ! loop over locally owned cells
       bmltmax = 0.0d0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (bmlt_float(i,j) > bmltmax) then
                bmltmax = bmlt_float(i,j)
                imax = i
                jmax = j
             endif
          enddo
       enddo
       print*, ' '
       print*, 'i, j, max(bmlt_float):', imax, jmax, bmlt_float(imax,jmax)*scyr
       
       ! check convergence of melt rate in all grid cells
       ! if converged, then exit the outer loop

       err_melt = 0.d0

       do j = 1, ny
          do i = 1, nx
             if (abs(bmlt_float(i,j) - bmlt_float_old(i,j)) > err_melt) then
                err_melt = abs(bmlt_float(i,j) - bmlt_float_old(i,j))
                imax = i
                jmax = j
             endif
          enddo
       enddo

       if (err_melt > maxerr_melt) then

          print*, ' '
          print*, 'Melt rate not converged: iter, rank, i, j, m_old, m, err target, errmax (m/yr) =', &
               iter_melt, this_rank, imax, jmax, bmlt_float_old(imax,jmax)*scyr, bmlt_float(imax,jmax)*scyr, maxerr_melt*scyr, err_melt*scyr
          converged_melt = .false.

       else  ! converged; exit the melt rate loop

          print*, ' '
          print*, 'Melt rate CONVERGED: iter, time(s), rank, i, j, m_old, m, err target, errmax (m/yr) =', &
               iter_melt, time, this_rank, imax, jmax, &
               bmlt_float_old(imax,jmax)*scyr, bmlt_float(imax,jmax)*scyr, maxerr_melt*scyr, err_melt*scyr
          converged_melt = .true.

       endif

       if (converged_continuity .and. converged_melt) then
          print*, ' '
          print*, 'Both continuity and melt have CONVERGED'
          exit
       elseif (time >= time_max) then
!!       elseif (iter_melt == maxiter_melt) then
          print*, ' '
          print*, 'Continuity and melt NOT CONVERGED after time =', time
          print*, 'Moving on'
          exit
       else   
          ! iterate again
       endif
       
    enddo   ! iter_melt
       
    ! Compute the final value of eta.
    ! The reason for this is that once eta > 0, it can never quite get back to 0,
    !  but only asymptotically approach 0 from above (assuming we are relaxing eta).
    ! This calculation sets eta = 0 wherever D_plume < lsrf - topg.
    !TODO - Check whether this correction leads to oscillations in D_plume and eta.
    where (plume_mask_cell == 1)
       eta_plume = max(D_plume - (lsrf - topg), 0.0d0)
    elsewhere
       eta_plume = 0.0d0
    endwhere

    ! Copy u_plume_east and v_plume_north into u_plume_Cgrid and v_plume_Cgrid for output.
    ! Note: v_plume_east and u_plume_north are used internally but are not part of output.
    u_plume_Cgrid(:,:) = u_plume_east(:,:)
    v_plume_Cgrid(:,:) = v_plume_north(:,:)
  
    !WHL - Tuning diagnostic
    ! Compute the mean melt rate in cells with lsrf < -300 m
    ! The goal for ISOMIP+ is to be close to 30 m/yr

    bmlt_float_avg = 0.d0
    ncells_sub300 = 0

    do j = 1, ny
       do i = 1, nx
          if (plume_mask_cell(i,j)==1 .and. lsrf(i,j) < -300) then
             ncells_sub300 = ncells_sub300 + 1
             bmlt_float_avg = bmlt_float_avg + bmlt_float(i,j)*scyr
          endif
       enddo
    enddo

    bmlt_float_avg = bmlt_float_avg/ncells_sub300
    print*, ' '
    print*, 'ncells_sub300, bmlt_float_avg:', ncells_sub300, bmlt_float_avg
    print*, ' '
    print*, 'Done in glissade_plume_melt_rate'
    print*, ' '

  end subroutine glissade_plume_melt_rate

!****************************************************

  subroutine compute_local_velocity(&
       i, j,                    &
       itest, jtest, rtest,     &
       u_tidal,                 &
       c_drag,                  &
       f_coriolis,              &
       D_plume,                 &
       pgf_x,                   &
       pgf_y,                   &
       u_plume,                 &
       v_plume,                 &
       plume_speed)
    
    integer, intent(in) ::  &
         i, j,              & ! grid cell coordinates (diagnostic only)
         itest, jtest, rtest  ! test cell coordinates (diagnostic only)
    
    real(dp), intent(in) ::   &
         u_tidal,           & ! tidal velocity (m/s)
         c_drag,            & ! ocean drag coefficient (unitless)   
         f_coriolis           ! Coriolis parameter (s^-1)
    
    ! Note: The following variables are co-located with the velocity
    real(dp), intent(in) ::   &
         D_plume,           & ! plume thickness (m)
         pgf_x,             & ! x component of pressure gradient force
         pgf_y                ! y component of pressure gradient force
    
    real(dp), intent(inout) ::  &
         u_plume,           & ! x component of plume velocity (m/s)
         v_plume              ! x component of plume velocity (m/s)

    real(dp), intent(out) ::  &
         plume_speed          ! plume speed (m/s)
    
    ! local variables

    real(dp) :: &
         x_resid, y_resid,  & ! residuals of momentum balance equations (m^2/s^2)
         denom,             & ! denominator
         a_uu, a_uv,        & ! coefficients for Newton solve
         a_vu, a_vv,        & !
         du, dv               ! change in u_plume and v_plume (m/s)
    
    integer ::  &
         iter_velo            ! iteration counter

    character(len=100) :: message

    integer, parameter ::  &
         maxiter_velo = 50    ! max number of iterations of velocity loop
    
    real(dp), parameter :: &
         maxresid_force_balance = 1.0d-8 ! max residual allowed in momentum balance equation (m^2/s^2)
    
    logical, parameter :: &
         velo_newton = .true.   ! if true, use Newton's method; if false, use Picard method

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
    !            b = elevation of shelf base
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
    ! The iteration is sped up by using Newton's method.
    ! We write   u = u0 + du
    !            v = v0 + dv
    !          |U| = U0 + d|U|/du * du + d|U|dv * dv
    ! where the partial derivatives are evaluated at (u,v) = (u_0,v_0).
    ! 
    ! This gives 
    !           du = (a_vv * R_x - a_uv * R_y) / det|A|
    !           dv = (a_uu * R_x - a_vu * R_x) / det|A|
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
    !--------------------------------------------------------------------

    do iter_velo = 1, maxiter_velo
       
       ! Compute plume speed based on current u and v
       plume_speed = sqrt(u_plume**2 + v_plume**2 + u_tidal**2)
       
       ! Compute residual of the momentum balance equations
       x_resid = pgf_x - c_drag*plume_speed*u_plume + f_coriolis*D_plume*v_plume
       y_resid = pgf_y - c_drag*plume_speed*v_plume - f_coriolis*D_plume*u_plume
       
       ! check convergence of plume velocity
       if (abs(x_resid) < maxresid_force_balance .and. abs(y_resid) < maxresid_force_balance) then
          if (this_rank == rtest .and. i==itest .and. j==jtest) then
             print*, 'Velocity converged: iter_velo, u/v_plume (m/s):', iter_velo, u_plume, v_plume
          endif
          exit
       elseif (iter_velo == maxiter_velo) then
          write(message,*) 'Error, glissade_plume: velocity has not converged, i, j =', i, j
          call write_log(message, GM_FATAL)
       endif
       
       if (velo_newton) then
          
          ! compute some coefficients for the Newton solve
          a_uu = c_drag * (plume_speed + u_plume**2/plume_speed)
          a_uv = c_drag * (u_plume*v_plume)/plume_speed - D_plume*f_coriolis
          a_vu = c_drag * (u_plume*v_plume)/plume_speed + D_plume*f_coriolis
          a_vv = c_drag * (plume_speed + v_plume**2/plume_speed)
          
          ! compute du and dv
          denom = a_uu*a_vv - a_uv*a_vu
          
          if (abs(denom) > 0.0d0) then
             du = (a_vv*x_resid - a_uv*y_resid) / denom
             dv = (a_uu*y_resid - a_vu*x_resid) / denom
             
             u_plume = u_plume + du
             v_plume = v_plume + dv
             
          else  ! denom = 0.0
             write(6,*) 'Error, glissade_plume: ill-posed Newton solve for velocity, rank, i, j:', this_rank, i, j
             write(6,*) 'a_uu, a_vv, a_uv, a_vu =', a_uu, a_vv, a_uv, a_vu
             write(message,*) 'Error, glissade_plume: ill-posed Newton solve for velocity, rank, i, j:', this_rank, i, j
             call write_log(message, GM_FATAL)
          endif
          
       else  ! simpler Picard solve
          
          denom = (c_drag*plume_speed)**2 + (D_plume*f_coriolis)**2
          u_plume = (c_drag*plume_speed*pgf_x + D_plume*f_coriolis*pgf_y) / denom
          v_plume = (c_drag*plume_speed*pgf_y - D_plume*f_coriolis*pgf_x) / denom
          
       endif

       if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
          print*, ' '
          print*, 'iter_velo, plume_speed (m/s) =', iter_velo, plume_speed
          print*, 'pgf_x, pgf_y:', pgf_x, pgf_y
          print*, 'Dfv, -Dfu:', D_plume * f_coriolis * v_plume, &
                               -D_plume * f_coriolis * u_plume
          print*, 'dragu, dragv:', c_drag * plume_speed * u_plume, &
                                   c_drag * plume_speed * v_plume
          print*, 'x/y residual:', x_resid, y_resid
          print*, 'new u/v_plume:', u_plume, v_plume
       endif
       
    enddo  ! iter_velo
    
  end subroutine compute_local_velocity

!****************************************************
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
    logical, parameter :: verbose = .false.

    ! compute coefficients of depressed cubic, y^3 + py + q = 0

    p = (3.d0*c/a - (b/a)**2) / 3.d0
    q = (2.d0*(b/a)**3 - 9.d0*b*c/(a*a) + 27.d0*d/a) / 27.d0

    ! compute the discriminant
    Delta = (p/3.d0)**3 + (q/2.d0)**2

    if (verbose) then
       print*, 'Delta =', Delta
       if (Delta > 0.d0) then
          print*, 'One real root, 2 complex conjugate'
       elseif (Delta == 0.d0) then
          print*, 'Three real roots of which at least two are equal'
       elseif (Delta < 0.d0) then
          print*, 'Three distinct real roots'
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

       if (verbose) then
          print*, 'a, b, c, d:', a, b, c, d
          print*, 'p, q:', p, q
          print*, 'y1 =', y1
          print*, 'x1 =', x1
       endif

    else  ! Delta < 0; three distinct real roots
          ! use a trigonometric formulation

       phi = acos(-q/(2.d0*sqrt(abs(p)**3/27.d0)))

       y1 =    2.d0 * sqrt(abs(p)/3.d0) * cos(phi/3.d0)
       y2_r = -2.d0 * sqrt(abs(p)/3.d0) * cos((phi+pi)/3.d0)
       y2_i =  0.d0
       y3_r = -2.d0 * sqrt(abs(p)/3.d0) * cos((phi-pi)/3.d0)
       y3_i =  0.d0

       if (verbose) then
          print*, 'a, b, c, d:', a, b, c, d
          print*, 'p, q:', p, q
          print*, 'y1, y2, y3 =', y1, y2_r, y3_r
          print*, 'b/3a =', b/(3.d0*a)
          print*, 'x1 =', y1 - b/(3.d0*a)
          print*, 'x2 =', y2_r - b/(3.d0*a)
          print*, 'x3 =', y3_r - b/(3.d0*a)
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

end module glissade_bmlt_float

!****************************************************
