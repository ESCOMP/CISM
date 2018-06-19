!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_bmlt_float.F90 - part of the Community Ice Sheet Model (CISM)  
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

  logical :: verbose_velo = .true.
  logical :: verbose_continuity = .true.
  logical :: verbose_melt = .true.

    !WHL - Should the MISOMIP parameters go elsewhere?
    !      Note: gammaS and gammaT are namelist parameters

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
!!        f_coriolis = 0.0d0            ! Coriolis parameter (s^-1) at 75 S = 2*omega*sin(75 deg) (prescribed in text)

    !TODO - Put each parameter in the appropriate subroutine, or remove it.
    ! relaxation parameters
    ! Value of 1 means to use the new value.  Lower values give a greater contribution from the old value.
    real(dp), parameter :: &
         relax_u = 1.0d0,    &
!!         relax_u = 0.5d0,    &
         relax_E = 1.0d0,    &
!!         relax_D = 0.5d0,    &   ! relax_D and relax_eta moved to thickness solver
!!         relax_eta = 0.01d0, &
!!         relax_m = 1.0d0
!!         relax_TS = 0.5d0
         relax_TS = 1.0d0

    !WHL - debug - not sure this matters much
    logical, parameter :: cap_Dplume = .true.
!!    logical, parameter :: cap_Dplume = .false.

contains

!****************************************************

  subroutine glissade_basal_melting_float(whichbmlt_float,              &
                                          ewn,         nsn,             &
                                          dew,         dns,             &
                                          itest,       jtest,    rtest, &
                                          x1,                           &
                                          thck,        lsrf,            &
                                          topg,        eus,             &
                                          basal_melt,  plume)

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
         lsrf,                 & ! elevation of lower ice surface (m)
         thck                    ! ice thickness (m)

    real(dp), dimension(:,:), intent(in) :: &
         topg                    ! elevation of bed topography (m)

    real(dp), intent(in) :: &
         eus                     ! eustatic sea level (m), = 0. by default

    type(glide_basal_melt), intent(inout) :: &
         basal_melt              ! derived type with fields and parameters related to basal melting

    type(glide_plume), intent(inout) :: &
         plume                   ! derived type with fields and parameters for plume model

    !-----------------------------------------------------------------
    ! Note: The basal_melt derived type includes the 2D output field bmlt_float,
    !        along with a number of prescribed parameters for MISMIP+:
    !
    !       MISMIP+ Ice1
    !       - bmlt_float_omega     ! time scale for basal melting (s-1), default = 0.2/yr
    !       - bmlt_float_h0        ! scale for sub-shelf cavity thickness (m), default = 75 m
    !       - bmlt_float_z0        ! scale for ice draft (m), default = -100 m
    !
    !       MISMIP+ Ice2
    !       - bmlt_float_const     ! constant melt rate (m/s), default = 100 m/yr
    !       - bmlt_float_xlim      ! melt rate = 0 for abs(x) < bmlt_float_xlim (m), default = 480000 m
    !
    ! Note: The plume derived type includes plume-related 2D output fields,
    !        along with a number of prescribed parameters for MISOMIP:
    !
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
    ! Local variables and pointers set to components of basal_melt and plume derived types
    !----------------------------------------------------------------      

    real(dp), dimension(:,:), pointer :: &
         bmlt_float             ! basal melt rate for floating ice (m/s) (> 0 for melt, < 0 for freeze-on)

    real(dp), dimension(:,:), pointer :: &
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

    real(dp) :: frz_ramp_factor    ! multiplying factor for linear ramp at depths with basal freezing
    real(dp) :: melt_ramp_factor   ! multiplying factor for linear ramp at depths with basal melting

    logical, parameter :: verbose_bmlt = .false.

!TODO - Make first_call depend on whether we are restarting
!!    logical :: first_call = .false.
    logical :: first_call = .true.

    !-----------------------------------------------------------------
    ! Compute the basal melt rate for floating ice
    !-----------------------------------------------------------------

    if (main_task .and. verbose_bmlt) print*, 'Computing bmlt_float, whichbmlt_float =', whichbmlt_float

    ! Set bmlt_float pointer and initialize
    bmlt_float  => basal_melt%bmlt_float
    bmlt_float(:,:) = 0.0d0

    ! Compute masks:
    ! - ice_mask = 1 where thck > 0
    ! - floating_mask = 1 where thck > 0 and ice is floating;
    ! - ocean_mask = 1 where topg is below sea level and ice is absent
    !Note: The '0.0d0' argument is thklim. Here, any ice with thck > 0 gets ice_mask = 1.

    call glissade_get_masks(ewn,           nsn,            &
                            thck,          topg,           &
                            eus,           0.0d0,          &
                            ice_mask,                      &
                            floating_mask = floating_mask, &
                            ocean_mask = ocean_mask)

    if (whichbmlt_float == BMLT_FLOAT_CONSTANT) then

       ! Set melt rate to a constant value for floating ice.
       ! This includes ice-free ocean cells, in case ice is advected to those cells by the transport scheme.

       do j = 1, nsn
          do i = 1, ewn

             if (floating_mask(i,j) == 1 .or. ocean_mask(i,j) == 1) then   ! ice is present and floating, or ice-free ocean

                ! Note: For MISMIP+ experiment Ice2r, melting is masked out where x < 480 km

                if (abs(x1(i)) >= basal_melt%bmlt_float_xlim) then   ! melting is allowed
                   bmlt_float(i,j) = basal_melt%bmlt_float_const
                endif

                !WHL - debug
                if (j==jtest .and. this_rank==rtest) then
!!                if (i==itest .and. j==jtest .and. this_rank==rtest) then
!!                   print*, 'rank, i, j, bmlt_float:', this_rank, i, j, bmlt_float(i,j)
                endif
                   
             endif   ! ice is present and floating

          enddo
       enddo

    elseif (whichbmlt_float == BMLT_FLOAT_MISMIP) then   ! MISMIP+

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

             ! compute basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme

             if (floating_mask(i,j) == 1 .or. ocean_mask(i,j) == 1) then   ! ice is present and floating, or ice-free ocean

                h_cavity = lsrf(i,j) - topg(i,j)
                z_draft = lsrf(i,j) - eus
                bmlt_float(i,j) = basal_melt%bmlt_float_omega * tanh(h_cavity/basal_melt%bmlt_float_h0) &
                                * max(basal_melt%bmlt_float_z0 - z_draft, 0.0d0)

                !debug
!                if (j == jtest .and. verbose_bmlt) then
!                   print*, 'cavity, tanh, thck, draft, melt rate (m/yr):', i, j, h_cavity, &
!                         tanh(h_cavity/basal_melt%bmlt_float_h0), thck(i,j), z_draft, bmlt_float(i,j)*scyr
!                endif

             endif   ! ice is present and floating

          enddo
       enddo

       !WHL - debug
       if (verbose_bmlt .and. this_rank == rtest) then
          print*, 'itest, jtest, rtest =', itest, jtest, rtest
          print*, ' '
          print*, 'thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'bmlt_float (m/yr), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') bmlt_float(i,j)*scyr
             enddo
             write(6,*) ' '
          enddo
       endif

    elseif (whichbmlt_float == BMLT_FLOAT_DEPTH) then

       ! Compute melt rates as a piecewise linear function of depth, generally with greater melting at depth.
       ! This scheme is similar to the MISMIP scheme, with the additional option of near-surface freezing.
       ! The maximum melting and freezing rates are set independently, with melting usually of greater magnitude.
       ! The melting/freezing rates fall linearly from their max values to zero over ranges defined by
       !  zmeltmax, zmelt0 and zfrzmax.
       ! The melt rate is set to a maximum value where z_draft <= zmeltmax,
       !  then decreases linearly to 0 as z_draft increases from zmeltmax to zmelt0.
       ! The freezing rate is set to a maximum value where z_draft >= zfrzmax,
       !  then decreases linearly to 0 as z_draft decreases from zfrzmax to zmelt0.
       ! (Here, z_draft < 0 by definition.)

       ! Compute ramp factors
       ! These factors are set to avoid divzero whe zmelt0 = zmeltmax, or zmelt0 = zfrzmax

       if (basal_melt%bmlt_float_depth_zfrzmax > basal_melt%bmlt_float_depth_zmelt0) then
          frz_ramp_factor = 1.0d0 / (basal_melt%bmlt_float_depth_zfrzmax - basal_melt%bmlt_float_depth_zmelt0)
       else
          frz_ramp_factor = 0.0d0
       endif

       if (basal_melt%bmlt_float_depth_zmelt0 > basal_melt%bmlt_float_depth_zmeltmax) then
          melt_ramp_factor = 1.0d0 / (basal_melt%bmlt_float_depth_zmelt0 - basal_melt%bmlt_float_depth_zmeltmax)
       else
          melt_ramp_factor = 0.0d0
       endif

       do j = 1, nsn
          do i = 1, ewn

             ! compute basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme
             if (floating_mask(i,j) == 1 .or. ocean_mask(i,j) == 1) then   ! ice is present and floating, or ice-free ocean

                z_draft = lsrf(i,j) - eus

                if (basal_melt%warm_ocean_mask(i,j) == 1) then  ! warm ocean profile; enforce minimum melt rate

                   if (z_draft > basal_melt%bmlt_float_depth_zmeltmin) then
                      bmlt_float(i,j) = basal_melt%bmlt_float_depth_meltmin
                   elseif (z_draft > basal_melt%bmlt_float_depth_zmeltmax) then
                      ! melting with a linear ramp from meltmin to meltmax
                      bmlt_float(i,j) = basal_melt%bmlt_float_depth_meltmin  &
                           + (basal_melt%bmlt_float_depth_meltmax - basal_melt%bmlt_float_depth_meltmin) *  &
                           (z_draft - basal_melt%bmlt_float_depth_zmeltmin) &
                           / (basal_melt%bmlt_float_depth_zmeltmax - basal_melt%bmlt_float_depth_zmeltmin)
                   elseif (z_draft <= basal_melt%bmlt_float_depth_meltmax) then
                      ! max melting
                      bmlt_float(i,j) = basal_melt%bmlt_float_depth_meltmax
                   endif

                else   ! standard depth-dependent profile

                   if (z_draft > basal_melt%bmlt_float_depth_zfrzmax) then
                      ! max freezing
                      bmlt_float(i,j) = -basal_melt%bmlt_float_depth_frzmax   ! frzmax >=0 by definition
                   elseif (z_draft > basal_melt%bmlt_float_depth_zmelt0) then
                      ! freezing with a linear taper from frzmax to zero
                      bmlt_float(i,j) = -basal_melt%bmlt_float_depth_frzmax * &
                           frz_ramp_factor * (z_draft - basal_melt%bmlt_float_depth_zmelt0)
                   elseif (z_draft > basal_melt%bmlt_float_depth_zmeltmax) then
                      ! melting with a linear taper from meltmax to zero
                      bmlt_float(i,j) = basal_melt%bmlt_float_depth_meltmax * &
                           melt_ramp_factor * (basal_melt%bmlt_float_depth_zmelt0 - z_draft)
                   elseif (z_draft <= basal_melt%bmlt_float_depth_meltmax) then
                      ! max melting
                      bmlt_float(i,j) = basal_melt%bmlt_float_depth_meltmax
                   endif

                endif   ! warm_ocean_mask


                ! As with the MISMIP scheme, reduce the melting as the cavity thickness approaches zero.
                ! A small value of bmlt_float_h0 allows more melting in very thin cavities.
                h_cavity = max(lsrf(i,j) - topg(i,j), 0.0d0)
                bmlt_float(i,j) = bmlt_float(i,j) * tanh(h_cavity/basal_melt%bmlt_float_h0)

             endif   ! ice is present and floating

          enddo
       enddo

       !debug
       if (verbose_bmlt .and. this_rank == rtest) then
          print*, 'itest, jtest, rtest =', itest, jtest, rtest
          print*, ' '
          print*, 'topg (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') topg(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'z_draft (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') min(lsrf(i,j) - eus, 0.0d0)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'h_cavity (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') max(lsrf(i,j) - topg(i,j), 0.0d0)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'bmlt_float (m/yr), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') bmlt_float(i,j)*scyr
             enddo
             write(6,*) ' '
          enddo
       endif

    elseif (whichbmlt_float == BMLT_FLOAT_MISOMIP) then

       ! Compute melt rates using a plume model, given vertical profiles of T and S in the ambient ocean
       !
       ! See this paper for details:
       ! X. S. Asay-Davis et al. (2016), Experimental design for three interrelated 
       !    marine ice sheet and ocean model intercomparison projects: 
       !    MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
       !    Geosci. Model Devel., 9, 2471-2497, doi: 10.5194/gmd-9-2471-2016.

       ! Assign local pointers and variables to components of the plume derived type

       T0 = plume%T0
       Tbot = plume%Tbot
       S0 = plume%S0
       Sbot = plume%Sbot
       zbed_deep = plume%zbed_deep
       gammaT = plume%gammaT
       gammaS = plume%gammaS

       ! the following fields are used or computed by the plume model
       T_basal       => plume%T_basal
       S_basal       => plume%S_basal
       u_plume       => plume%u_plume
       v_plume       => plume%v_plume
       u_plume_Cgrid => plume%u_plume_Cgrid
       v_plume_Cgrid => plume%v_plume_Cgrid
       ustar_plume   => plume%ustar_plume
       drho_plume    => plume%drho_plume
       T_plume       => plume%T_plume
       S_plume       => plume%S_plume
       D_plume       => plume%D_plume
       entrainment   => plume%entrainment
       detrainment   => plume%detrainment
       divDu_plume   => plume%divDu_plume
       T_ambient     => plume%T_ambient
       S_ambient     => plume%S_ambient

       if (verbose_bmlt .and. this_rank == rtest) then
          print*, 'itest, jtest, rtest =', itest, jtest, rtest
          print*, ' '

          print*, ' '
          print*, 'thck:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, 'lsrf:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') lsrf(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'topg:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') topg(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'lsrf - topg:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.5)',advance='no') lsrf(i,j) - topg(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'floating_mask:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i8)',advance='no') floating_mask(i,j)
             enddo
             write(6,*) ' '
          enddo

       endif  ! verbose_bmlt

       ! Given the ice draft in each floating grid cell, compute the ambient ocean T and S
       !  using the prescribed MISOMIP profile.

       where (floating_mask == 1) 

          ! MISOMIP+ profiles, Eqs. 21 and 22 
!!          T_ambient(:,:) = T0 + (Tbot - T0) * (lsrf(:,:) / zbed_deep)
!!          S_ambient(:,:) = S0 + (Sbot - S0) * (lsrf(:,:) / zbed_deep)
          T_ambient = T0 + (Tbot - T0) * (lsrf / zbed_deep)
          S_ambient = S0 + (Sbot - S0) * (lsrf / zbed_deep)

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
            thck,                               &  ! temporary, for calving
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

       ! convert plume fields to scaled units
       entrainment(:,:) = entrainment(:,:) * tim0/thk0
       detrainment(:,:) = detrainment(:,:) * tim0/thk0
       divDu_plume(:,:) = divDu_plume(:,:) * tim0/thk0
       D_plume(:,:) = D_plume(:,:)/thk0

    endif   ! whichbmlt_float

    ! Set first_call to false. 
    ! Next time, the plume variables just computed (T_plume, S_plume, D_plume)
    !  will be taken as initial conditions.
    first_call = .false.

  end subroutine glissade_basal_melting_float

!****************************************************

  subroutine glissade_plume_melt_rate(&
       first_call,                      &
       nx,               ny,            &
       dx,               dy,            &
       x1,                              &
       thck,                            &  ! temporary, for calving
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

!!    real(dp), dimension(nx,ny), intent(inout) ::  &
    real(dp), dimension(nx,ny), intent(in) ::  &
         thck                   ! ice thickness (m); intent(inout) to allow calving

!!    real(dp), dimension(nx,ny), intent(inout) ::  &
    real(dp), dimension(nx,ny), intent(in) ::  &
         lsrf                   ! ice lower surface elevation (m, negative below sea level)
                                ! intent(inout) to allow calving

    real(dp), dimension(nx,ny), intent(in) ::  &
         topg                   ! bedrock elevation (m, negative below sea level)

    !WHL - Change to intent(inout) to allow calving of floating ice
!!    integer, dimension(nx,ny), intent(inout) :: &
    integer, dimension(nx,ny), intent(in) :: &
         floating_mask          ! = 1 where ice is present and floating, else = 0

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

    !TODO - Change 'old' to 'latest'.  'old' is at the start of the time step, 'latest' at the start of the iteration.
    real(dp), dimension(nx,ny) :: &
         pressure,            & ! ocean pressure at base of ice (N/m^2)
         lsrf_plume,          & ! elevation of plume-ambient interface (m, negative below sea level)
         rho_plume,           & ! plume density (kg/m^3)
         rho_ambient,         & ! ambient ocean density (kg/m^3)
         H_cavity,            & ! thickness of ocean cavity beneath the plume (m)
         D_plume_cap,         & ! min(D_plume, H_cavity)
         eta_plume,           & ! displacement of plume surface, D_plume - H_cavity (m)
         dD_plume,            & ! change in D_plume (m)
         T_plume_latest,      & ! T_plume from latest iteration
         S_plume_latest,      & ! S_plume from latest iteration
         T_basal_latest,      & ! T_basal from latest iteration
         S_basal_latest,      & ! S_basal from latest iteration
         bmlt_float_latest,   & ! melt rate from latest iteration (m/s)
         T_plume_old,         & ! T_plume from previous time step
         S_plume_old,         & ! S_plume from previous time step
         T_basal_old,         & ! T_basal from previous time step
         S_basal_old,         & ! S_basal from previous time step
         D_plume_old,         & ! D_plume from previous time step
         eta_plume_old,       & ! eta_plume from previous time step
         drho_plume_old,      & ! drho_plume from previous time step
         bmlt_float_old         ! melt rate from previous time step (m/s)

    real(dp), dimension(nx,ny) ::  &
         u_plume_east,          & ! u_plume on east edges
         v_plume_east,          & ! v_plume on east edges
         u_plume_north,         & ! u_plume on north edges
         v_plume_north,         & ! v_plume on north edges
         plume_speed_east,      & ! plume speed on east edges (m/s)
         plume_speed_north        ! plume speed on north edges (m/s)

    !TODO - Old values might not be needed  
    real(dp), dimension(nx,ny) ::  &
         entrainment_latest,    & ! entrainment from previous iteration
         u_plume_east_latest,   & ! latest values of u_plume on east edges, from previous iteration
         v_plume_east_latest,   & ! latest values of v_plume on east edges
         u_plume_north_latest,  & ! old values of u_plume on north edges
         v_plume_north_latest     ! old values of v_plume on north edges

    real(dp), dimension(nx,ny) ::  &
         dlsrf_dx_east,       & ! horizontal gradient of lsrf on east edges
         dlsrf_dy_east,       & !
         dlsrf_dx_north,      & ! horizontal gradient of lsrf on north edges
         dlsrf_dy_north

    real(dp) :: &
         dlsrf_dx, dlsrf_dy,  & ! lsrf gradient components at cell centers
         slope                  ! magnitude of the gradient (dlsrf_dx, dlsrf_dy)

    real(dp), dimension(nx,ny) ::  &
         theta_slope            ! basal slope angle (rad), used for entrainment
 
    real(dp), dimension(nx-1,ny-1) :: &
         plume_speed            ! plume speed at vertices (m/s)
!!         dlsrf_plume_dx,      & ! horizontal gradient of lsrf_plume
!!         dlsrf_plume_dy

    ! Note: We have edge_mask_east/north = 1 only if the plume exists on each side of an edge.
    !       We can have divu_mask_east/north = 1 if the plume exists on only one side of the edge,
    !        with floating ice or open water on the other. At these edges, velocity is not computed
    !        directly but is extrapolated from a neighbor.
    !       Thus, cells with edge_mask_east/north are a subset of cells with divu_mask_east/north = 1. 

    integer, dimension(nx,ny) :: &
         plume_mask_cell,     & ! = 1 for cells where scalar plume variables are computed
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed
         edge_mask_north,     & ! = 1 on north edges where plume velocity is computed
         divu_mask_east,      & ! = 1 on east edges where divergence terms are computed
         divu_mask_north        ! = 1 on north edges where divergence terms are computed

    !TODO - New fields for testing
    ! Note: These masks are real-valued because they can have a value of 0, 0.5 or 1
    real(dp), dimension(nx,ny) :: &
         edge_mask_east_reduce_v,  & ! mask for reducing v on east edges adjacent to a wall
         edge_mask_north_reduce_u    ! mask for reducing u on north edges adjacent to a wall

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
         lsrf_min                         ! global min value of lsrf (m)

    real(dp) :: &
         plume_tendency,                & ! tendency of plume thickness (m/s)
         D_plume_east, D_plume_west,    & ! terms in discretization of plume divergence
         D_plume_north, D_plume_south,  &
         dDu_dx, dDv_dy,                &
         eta_plume_unrelaxed              ! value of eta_plume without relaxation

    real(dp) :: &
         bmlt_float_avg         ! average value of bmlt_float in main cavity

    real(dp) :: &
         time,                & ! elapsed time during the relaxation of the plume thickness (s)            
         dt_plume,            & ! time step (s)
         my_max_dt              ! CFL-limited time step for a given cell (s) 

    real(dp) ::  &
         L2_norm,             & ! L2 norm of residual vector from continuity equation
         L2_previous            ! L2 norm from the previous convergence check 

    integer, parameter :: &
         n_check_residual = 1     ! how often to compute the residual and check for convergence

    integer :: i, j
    integer :: iglobal, jglobal        ! global i and j indices
    integer :: iter_Dplume, iter_melt  ! iteration counters
    integer :: ncells_sub300           ! number of cells below 300 m depth (ISOMIP+ diagnostic)

    integer :: imax, jmax              ! i and j indices of cells with extreme values

    ! max of various quantities for a given iteration
    real(dp) :: &
         Dmax, etamax, detamax, speedmax, entrainmax, bmltmax

    real(dp) :: &
         max_tendency,        & ! max plume tendency; measure of convergence of continuity equation
         err_melt               ! max difference in bmlt_float between iterations

    logical :: &
         converged_continuity,   & ! true if continuity equation has converged in all cells, else = false
         converged_melt            ! true if melt rate has converged in all cells, else = false

    ! Parameters in the plume model
    real(dp), parameter :: &    
         plume_xmax = 630000.d0,     & ! limit of the plume in the x direction (m), determined by the calving front location
         D_plume0 = 1.d0,            & ! initial plume thickness at lowest elevation, lsrf_min (m)
         D_plume_dz = 0.02d0           ! rate of change of initial plume thickness with increasing z (m/m)

    real(dp), parameter :: &
         dt_plume_max = 300.d0,    & ! max time step for plume thickness iteration (s)
                                    ! Shortened as needed to satisfy CFL
!!         time_max = scyr            ! max time (s) before giving up on convergence 
         time_max = 10000.d0      ! max time (s) before giving up on melt rate convergence 

    real(dp), parameter ::  &
         L2_target = 1.0d-6           ! small target value for L2 norm of continuity equation residual

    !TODO - Currently not working with free_surface = false.  eta inflates anyway.
    logical, parameter :: free_surface = .false.
!!    logical, parameter :: free_surface = .true.   ! true if computing PG force due to slope in free surface

    ! parameters determining convergence of iterations
    integer, parameter :: &
         maxiter_melt = 100,         & ! max number of iterations of inner melt-rate loop
                                       ! terminates when bmlt_float, u_plume, v_plume and dD_plume/dt are consistent
         maxiter_Dplume = 999999       ! max number of iterations of outer plume-thickness loop
                                       ! terminates when plume thickness reaches virtual steady state
    real(dp), parameter :: &
         maxerr_melt = 1.0d-3/scyr     ! max err_melt allowed for steady state (m/yr converted to m/s)

    real(dp), parameter :: &
         eta_plume_min = 1.0d-8       ! threshold thickness (m) for eta_plume
                                      ! Set eta_plume = 0 when eta_plume < eta_plume_min

    real(dp), parameter :: &
!!         H_cavity_min = 1.0d0        ! threshold cavity thickness (m) for plume to exist
         H_cavity_min = 0.0d0        ! threshold cavity thickness (m) for plume to exist

    !WHL - Calving in this subroutine not currently supported
!!    real(dp), parameter :: &
!!         thck_min = 100.d0           ! threshold thickness (m) for floating ice; thinner ice calves

    character(len=6) ::  &
         nonlinear_method             ! 'Picard' or 'Newton'

    !WHL - debug and diagnostics
    real(dp) :: dD_dt
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

    !WHL - debug
    ! Calve thin floating ice if necessary.  Generally, this should be done by CISM's calving solver.
    !TODO - If uncommenting these lines, then thck and lsrf must be intent(inout)
!!    do j = 1, ny
!!       do i = 1, nx
!!          if (floating_mask(i,j)==1 .and. thck(i,j) > 0.0d0 .and. thck(i,j) < thck_min) then
!!             print*, 'Calve thin ice: i, j, thck, topg =', i, j, thck(i,j), topg(i,j)
!!             thck(i,j) = 0.0d0
!!             lsrf(i,j) = 0.0d0
!!             floating_mask(i,j) = 0
!!          endif
!!       enddo
!!    enddo

    ! Compute the density of the ambient ocean

    rho_ambient(:,:) = eos_rho_ref * (1.d0 - eos_alpha * (T_ambient(:,:) - eos_Tref)  &
                                           + eos_beta  * (S_ambient(:,:) - eos_Sref) )

    ! Compute the pressure at the lower ice surface.
    pressure(:,:) = -rhoo*grav*lsrf(:,:)

    ! Compute a mask for where the plume thickness is computed.
    ! Initialize to agree with floating_mask.
    plume_mask_cell(:,:) = floating_mask(:,:)

    ! Mask out cells that are not locally owned.
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

    !WHL - debug
    !TODO - Support cavities of any thickness, no matter how small?
    ! Compute the ocean cavity thickness.
    ! Optionally, mask out cells with very narrow cavities.
    ! Since entrainment goes to zero in small cavities, basal melt in these cells should be small. 

    H_cavity(:,:) = max(lsrf(:,:) - topg(:,:), 0.0d0)
    do j = 1, ny
       do i = 1, nx
          if (H_cavity(i,j) > 0.0d0 .and. H_cavity(i,j) < H_cavity_min) then
             plume_mask_cell(i,j) = 0
             print*, 'Mask out thin cavity: i, j, H_cavity =', i, j, H_cavity(i,j)
          endif
       enddo
    enddo

    !WHL - debug
    print*, ' '
    print*, 'H_cavity:'
    do j = jtest+3, jtest-3, -1
       do i = itest-3, itest+3
          write(6,'(f12.5)',advance='no') H_cavity(i,j)
       enddo
       write(6,*) ' '
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

    ! Compute masks for inhibiting flow toward walls of grounded ice

    ! Initialize masks to 1.0 (implying no reduction of the velocity component)
    edge_mask_east_reduce_v(:,:) = 1.0d0
    edge_mask_north_reduce_u(:,:) = 1.0d0

    ! Reset the masks to 0.0 or 0.5 adjacent to grounded ice
    do j = 1, ny
       do i = 1, nx

          ! identify east edges with a wall of grounded ice to the north or south
          if (edge_mask_east(i,j) == 1) then
             if ( (H_cavity(i,j+1) == 0.0d0 .and. H_cavity(i+1,j+1) == 0.0d0) .or.  &
                  (H_cavity(i,j-1) == 0.0d0 .and. H_cavity(i+1,j-1) == 0.0d0) ) then
                ! full wall; zero out the v component
                edge_mask_east_reduce_v(i,j) = 0.0d0
             elseif (H_cavity(i,j+1) == 0.0d0 .or. H_cavity(i+1,j+1) == 0.0d0 .or.  &
                     H_cavity(i,j-1) == 0.0d0 .or. H_cavity(i+1,j-1) == 0.0d0) then
                ! half wall; reduce the v component
                edge_mask_east_reduce_v(i,j) = 0.5d0
             endif
          endif

          ! identify north edges with a wall of grounded ice to the east or west
          if (edge_mask_north(i,j) == 1) then
             if ( (H_cavity(i-1,j+1) == 0.0d0 .and. H_cavity(i+1,j+1) == 0.0d0) .or.  &
                  (H_cavity(i-1,j)   == 0.0d0 .and. H_cavity(i+1,j)   == 0.0d0) ) then
                ! full wall; zero out the v component
                edge_mask_north_reduce_u(i,j) = 0.0d0
             elseif (H_cavity(i-1,j+1) == 0.0d0 .or. H_cavity(i+1,j+1) == 0.0d0 .or.  &
                     H_cavity(i-1,j)   == 0.0d0 .or. H_cavity(i+1,j)   == 0.0d0) then
                ! half wall; reduce the v component
                edge_mask_north_reduce_u(i,j) = 0.5d0

             endif
          endif
          
       enddo
    enddo

    !TODO - Check whether these comments are still acurate
    ! Compute masks for edges with nonzero fluxes.
    ! These masks includes all cells with edge_mask_east/north = 1 (where velocity is computed).
    ! In addition, these masks include cells that have a plume on one side of the edge and
    !  open water or floating ice on the other.
    ! These edges have nonzero velocity extrapolated from neighboring edges, and thus are included
    !  in computations of the divergence.

    divu_mask_east(:,:) = edge_mask_east(:,:)
    divu_mask_north(:,:) = edge_mask_north(:,:)

    ! east edges
    do j = 1, ny
       do i = 1, nx-1
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 0 .and. global_bndy_east(i,j) == 0) then
             if (lsrf(i+1,j) == 0.0d0 .or. floating_mask(i+1,j) == 1) then
                ! water in cell (i+1,j); get plume velocity from edge (i-1,j)
                divu_mask_east(i,j) = 1
             endif
          elseif (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i+1,j) == 1 .and. global_bndy_west(i+1,j) == 0) then
             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                ! water in cell (i,j); get plume velocity from edge (i+1,j)
                divu_mask_east(i,j) = 1
             endif
          endif
       enddo
    enddo

    ! north edges
    do j = 1, ny
       do i = 1, nx
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 0 .and. global_bndy_north(i,j) == 0) then
             if (lsrf(i,j+1) == 0.0d0 .or. floating_mask(i,j+1) == 1) then
                ! water in cell (i,j+1); get plume velocity from edge (i,j-1)
                divu_mask_north(i,j) = 1
             endif
          elseif (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i,j+1) == 1 .and. global_bndy_south(i,j+1) == 0) then
             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                ! water in cell (i,j); get plume velocity from edge (i,j+1)
                divu_mask_north(i,j) = 1
             endif
          endif
       enddo   ! i
    enddo   ! j

    ! Compute the horizontal gradient of the lower ice surface.
    ! This is used to compute the pressure gradient force at velocity points, and for entrainment.
    ! Note: There are a couple of different ways to compute the PGF.
    !       (1) Jenkins et al. (1991) and HJH (2008) use grad(lsrf)
    !       (2) Holland & Feltham (2006) use grad(lsrf_plume) along with a density gradient.
    !       Method (1) is simpler and has the advantage that grad(lsrf) does not vary during plume evolution,
    !        making the PGF more stable (though possibly not as accurate).
    ! Note: The first 'lsrf' is a required argument for the subroutine.
    !       The second 'lsrf' happens to be the field whose gradient we're computing.

    call compute_edge_gradients(&
         nx,              ny,          &
         dx,              dy,          &
         global_bndy_east,             &
         global_bndy_west,             &
         global_bndy_north,            &
         global_bndy_south,            &
         plume_mask_cell,              &
         floating_mask,                &
         lsrf,                         &
         lsrf,                         &
         dlsrf_dx_east,      dlsrf_dy_east,  &
         dlsrf_dx_north,     dlsrf_dy_north)
    
    ! Compute the slope angle at cell centers.  This is used to compute entrainment.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo 
          dlsrf_dx = (dlsrf_dx_east(i-1,j) + dlsrf_dx_east(i,j)) / 2.d0
          dlsrf_dy = (dlsrf_dy_north(i,j-1) + dlsrf_dy_north(i,j)) / 2.d0
          slope = sqrt(dlsrf_dx**2 + dlsrf_dy**2)
          theta_slope(i,j) = atan(slope)
       enddo
    enddo

    ! is this call needed?
    call parallel_halo(theta_slope)

    print*, ' '
    print*, 'plume_mask_cell, rank =', rtest
    do j = jtest+3, jtest-3, -1
       do i = itest-3, itest+3
          write(6,'(i8)',advance='no') plume_mask_cell(i,j)
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

    print*, ' '
    print*, 'edge_mask_north_reduce_u, rank =', rtest
    do j = jtest+3, jtest-3, -1
       do i = itest-3, itest+3
          write(6,'(f8.1)',advance='no') edge_mask_north_reduce_u(i,j)
       enddo
       write(6,*) ' '
    enddo
    
    print*, ' '
    print*, 'edge_mask_east_reduce_v, rank =', rtest
    do j = jtest+3, jtest-3, -1
       write(6,'(a6)',advance='no') '      '
       do i = itest-3, itest+3
          write(6,'(f8.1)',advance='no') edge_mask_east_reduce_v(i,j)
       enddo
       write(6,*) ' '
    enddo

    print*, ' '
    print*, 'divu_mask_north, rank =', rtest
    do j = jtest+3, jtest-3, -1
       do i = itest-3, itest+3
          write(6,'(i8)',advance='no') divu_mask_north(i,j)
       enddo
       write(6,*) ' '
    enddo
    
    print*, ' '
    print*, 'divu_mask_east, rank =', rtest
    do j = jtest+3, jtest-3, -1
       write(6,'(a6)',advance='no') '      '
       do i = itest-3, itest+3
          write(6,'(i8)',advance='no') divu_mask_east(i,j)
       enddo
       write(6,*) ' '
    enddo

    !----------------------------------------------------------------
    ! Initialize some fields that are updated during the iteration.
    ! Note: T_plume, S_plume and D_plume are intent(inout) and already have initial values.
    !----------------------------------------------------------------

    if (first_call) then

       print*, ' '
       print*, 'First call: creating simple initial conditions for T_plume, S_plume and D_plume'

       ! Initialize the plume temperature and salinity.
       ! Set S_plume = S0 everywhere.
       ! This means that drho_plume = rho_ambient - rho_plume will decrease in the upslope direction.
       ! Setting both T_plume and S_plume to ambient values would give zero velocities and melt rates. 

       where (plume_mask_cell == 1)
          T_plume = T_ambient
          S_plume = S0
       endwhere

       ! Initialize the plume thickness.
       ! This is tricky. Since entrainment is non-negative and is equal to div*(Du),
       !  we ideally want div*(Du) > 0 in most cells. If the flow is upslope and D increases
       !  upslope, then div*(Du) will generally be positive in most of the domain.
!!!!       ! Initally, D_plume is constrained not to be thicker than the sub-shelf ocean cavity.
!!!!       ! If convergence results in D_plume > H_cavity, there will be a pressure gradient force 
!!!!       !  tending to reduce D_plume.
       ! Note: Units for D_plume here are meters.
       !       Would have to change if D_plume is input in scaled model units

       lsrf_min = minval(lsrf)
       lsrf_min = parallel_reduce_min(lsrf_min)
       !WHL - Use an absolute level instead of the global min?

       do j = 1, ny
          do i = 1, nx
             if (plume_mask_cell(i,j) == 1) then
                D_plume(i,j) = D_plume0 + D_plume_dz * (lsrf(i,j) - lsrf_min)
!!                D_plume(i,j) = min(D_plume(i,j), H_cavity(i,j))
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


    if (verbose_continuity) then

       if (free_surface) then
          print*, ' '
          print*, 'Free surface calculation is ON'
       else
          print*, ' '
          print*, 'Free surface calculation is OFF'
       endif   ! free_surface

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

    ! Initialize other fields
    u_plume(:,:) = 0.0d0
    v_plume(:,:) = 0.0d0
    plume_speed(:,:) = 0.0d0

    u_plume_east(:,:) = 0.0d0
    v_plume_east(:,:) = 0.0d0
    u_plume_north(:,:) = 0.0d0
    v_plume_north(:,:) = 0.0d0
    plume_speed_east(:,:) = 0.0d0
    plume_speed_north(:,:) = 0.0d0

    ustar_plume(:,:) = 0.0d0
    divDu_plume(:,:) = 0.0d0
    entrainment(:,:) = 0.0d0
    detrainment(:,:) = 0.0d0
    bmlt_float(:,:) = 0.0d0

    eta_plume(:,:) = 0.0d0

    !WHL - debug
    if (main_task) then
       print*, ' '
       print*, 'Start melt-rate iteration'
    endif

    time = 0.0d0

    !--------------------------------------------------------------------
    ! Iterate to compute the melt rate at the ice-ocean interface
    ! The solution method is:
    ! (1) Given the current guesses for entrainment and friction velocity, compute bmlt_float
    !     and new values of T_basal, S_basal, T_plume and S_plume.
    ! (2) Given the current guesses for D_plume and rho_plume, compute u_plume and v_plume
    !      from the momentum balance.
    ! (3) Given u_plume and v_plume, compute the entrainment.
    ! (4) Given u_plume, v_plume and entrainment, compute the change in D_plume from continuity.
    ! There is an inner loop that iterates until u_plume, v_plume, dDplume/dt and bmlt_float
    !  have converged, meaning that the dynamic and thermodynamic fields are mutually consistent.
    ! This loop is wrapped by an outer loop that continues until D_plume is sufficiently
    !  close to steady state for all cells.
    !--------------------------------------------------------------------

    do iter_Dplume = 1, maxiter_Dplume   ! outer plume_thickness iteration

       if (main_task) then 
          print*, ' '
          print*, 'iter_D_plume =', iter_Dplume
       endif

       !TODO - Remove some of these if not needed for relaxation
       ! save variables from previous time step
       D_plume_old(:,:) = D_plume(:,:)

       if (free_surface) eta_plume_old(:,:) = eta_plume(:,:)

       bmlt_float_old(:,:) = bmlt_float(:,:)
       S_plume_old(:,:) = S_plume(:,:)
       T_plume_old(:,:) = T_plume(:,:)
       S_basal_old(:,:) = S_basal(:,:)
       T_basal_old(:,:) = T_basal(:,:)

       ! initialize the nonlinear method ('Newton' or 'Picard')
       !TODO - Test with an initial Picard
       nonlinear_method = 'Newton'
!!!       nonlinear_method = 'Picard'

       ! initialize the L2 norm to an arbitrary big number
       L2_previous = huge(0.0d0)

       ! Compute the plume velocity, solve the continuity equation for the plume thickness, and compute the melt rate.
       ! This is done iteratively until convergence.

       do iter_melt = 1, maxiter_melt

          if (main_task) then 
             print*, ' '
             print*, 'iter_melt = ', iter_melt
             print*, 'nonlinear_method = ', trim(nonlinear_method)
          endif

          ! save values from latest iteration (used for relaxation)
          bmlt_float_latest(:,:) = bmlt_float(:,:)
          S_plume_latest(:,:) = S_plume(:,:)
          T_plume_latest(:,:) = T_plume(:,:)
          S_basal_latest(:,:) = S_basal(:,:)
          T_basal_latest(:,:) = T_basal(:,:)

          ! Compute the plume density, given the current estimates of T_plume and S_plume.
          ! Then find the density difference between the ambient ocean and the plume. 

          rho_plume(:,:) = eos_rho_ref * (1.d0 - eos_alpha * (T_plume(:,:) - eos_Tref)  &
                                               + eos_beta  * (S_plume(:,:) - eos_Sref) )

          where (plume_mask_cell == 1)
             drho_plume = rho_ambient - rho_plume
          elsewhere
             drho_plume = 0.0d0
          endwhere

          ! TODO - What to do where drho_plume < 0? Set plume_mask_cell = 0?
          ! TODO = Print where drho_plume < 0.

          if (verbose_melt) then

             print*, ' '
             print*, 'New drho_plume (kg/m^3), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f12.4)',advance='no') drho_plume(i,j)
                enddo
                write(6,*) ' '
             enddo
          
          endif  ! verbose_melt

          ! TODO - Currently, T_plume and S_plume must be relaxed to avoid oscillations as a result
          !         of changes in drho_plume.  Think about incorporating du/dT and du/dS terms into
          !         the velocity solve to suppress these oscillations and improve convergence.
          !         Dependence of drho_plume on T_plume and S_plume is straightforward. 
          ! But for now, try passing in the value of drho_plume at the old time.

          ! Compute the velocity (u_plume, v_plume), given the current estimates of drho_plume and D_plume.

          i = itest
          j = jtest
          print*, ' '
          print*, 'Before velocity: D_plume, eta_plume =', D_plume(i,j), eta_plume(i,j)

          call compute_plume_velocity(&
               nx,       ny,        &
               dx,       dy,        &
               itest, jtest, rtest, &
               plume_mask_cell,     &
!!               floating_mask,       &
!!               global_bndy_east,    &
!!               global_bndy_west,    &
!!               global_bndy_north,   &
!!               global_bndy_south,   &
               divu_mask_east,      &
               divu_mask_north,     &
               edge_mask_east,      &
               edge_mask_north,     &
               edge_mask_east_reduce_v,  &
               edge_mask_north_reduce_u, &
               free_surface,        &
!!               lsrf,                &
               dlsrf_dx_east,      dlsrf_dy_east,  &
               dlsrf_dx_north,     dlsrf_dy_north, &
               drho_plume,          &
               D_plume,             &
               eta_plume,           &
               H_cavity,            &
               u_plume_east,        &
               v_plume_east,        &
               u_plume_north,       &
               v_plume_north,       &
               plume_speed_east,    &
               plume_speed_north)

          ! Determine the time step based on a CFL condition.
          ! Should be stable with a CFL number up to 1.0, but limit to 0.5 to be on the safe side.

          !WHL - Is this necessary to do for each iteration?
          dt_plume = dt_plume_max
          imax = 1
          jmax = 1
          
          do j = nhalo+1, ny-nhalo
             do i = nhalo+1, nx-nhalo
                if (plume_mask_cell(i,j) == 1) then
                   my_max_dt = 0.5d0*dx / max( abs(u_plume_east(i,j)),  abs(u_plume_east(i-1,j)), &
                                               abs(v_plume_north(i,j)), abs(v_plume_north(i,j-1)) )
                   if (my_max_dt < dt_plume) then
                      dt_plume = my_max_dt
                      imax = i
                      jmax = j
                   endif
                endif
             enddo
          enddo

          if (main_task .and. dt_plume < dt_plume_max) then
             print*, 'Limited dt_plume =', dt_plume
          endif

          ! Compute the entrainment rate, given u_plume and v_plume.
          
          call compute_entrainment(&
               nx,         ny,      &
               dx,         dy,      &
               itest, jtest, rtest, &
               divu_mask_east,      &
               divu_mask_north,     &
               plume_mask_cell,     &
!!               floating_mask,       &
!!               global_bndy_east,    &
!!               global_bndy_west,    &
!!               global_bndy_north,   &
!!               global_bndy_south,   &
!!               lsrf,                &
!!               H_cavity,            &
!!               D_plume,             &
               theta_slope,         &
               u_plume_east,        &
               v_plume_north,       &
!!               plume_speed_east,    &
!!               plume_speed_north,   &
!!               dlsrf_dx_east,      dlsrf_dy_east,  &
!!               dlsrf_dx_north,     dlsrf_dy_north, &
               entrainment)

          !WHL - Relaxing.  May not be needed.
!       entrainment(:,:) = (1.0d0 - relax_E)*entrainment_latest(:,:) + relax_E*entrainment(:,:)

          ! Compute the detrainment rate, given D_plume.

          call compute_detrainment(&
               nx,       ny,   &
               itest, jtest, rtest, &
               free_surface,   &
               dt_plume,       &
               H_cavity,       &
               D_plume,        &
               eta_plume,      &
               detrainment)
          
          ! Check for convergence of the continuity equation, dD/dt = e - d - del*(Du)

          !WHL - debug
          ! Here, I'd like to try a flexible strategy for convergence.  Maybe start with a few Picard solves,
          !  then switch to Newton as long as it's working, but come back to Picard if Newton is failing.
          ! That didn't work.  Try adjusting dt instead.

          if (iter_melt >=2 .and. mod(iter_melt, n_check_residual) == 0) then  ! time to check for convergence

             call compute_dynamic_residual(&
                  nx,       ny,        &
                  dx,       dy,        &
                  dt_plume,            &
                  itest, jtest, rtest, &
                  plume_mask_cell,     &
                  edge_mask_east,      &
                  edge_mask_north,     &
                  divu_mask_east,      &
                  divu_mask_north,     &
                  H_cavity,            &
                  entrainment,         &
                  detrainment,         &
                  D_plume_old,         &  ! value from previous time step
                  D_plume,             &
                  u_plume_east,        &
                  v_plume_north,       &
                  u_plume_north,       &  ! diagnostic only
                  v_plume_east,        &  ! diagnostic only
                  divDu_plume,         &
                  L2_norm)
             
             print*, ' '
             print*, 'COMPUTED RESIDUAL: iter_melt, L2_norm, L2_previous, L2_target:', &
                  iter_melt, L2_norm, L2_previous, L2_target

             if (L2_norm < L2_target) then
                print*, 'CONTINUITY CONVERGED, time, iter, L2_norm =', time, iter_melt, L2_norm
!!                exit
                if (converged_melt) then
                   print*, 'Melt rate has also converged; exit'
                   exit
                else
                   print*, 'Melt rate has not converged; continue'
                endif

             elseif (iter_melt == maxiter_melt) then
                print*, 'CONTINUITY FAILED TO CONVERGE, time, iter, L2_norm =', time, iter_melt, L2_norm
                print*, ' '
                exit
!!                stop
                !WHL - debug - Try a shorter time step next time?
             elseif (L2_norm < L2_previous) then ! iteration is converging; keep trying
                print*, 'CONTINUITY NOT YET CONVERGED, time, iter, L2_norm =', time, iter_melt, L2_norm
                if (trim(nonlinear_method) == 'Picard') then
!!!                   nonlinear_method = 'Newton'
!!!                   print*, 'Switching from Picard to Newton solve'
                endif
             elseif (L2_norm >= L2_previous) then ! iteration is not converging
                print*, 'CONTINUITY NOT CONVERGING, time, iter, L2_norm =', time, iter_melt, L2_norm
                ! if Newton, then try switching to Picard
                ! if Picard, then punt
                if (trim(nonlinear_method) == 'Newton') then
!!!                   print*, 'CONVERGENCE is failing with Newton; switch to Picard'
!!!                   nonlinear_method = 'Picard'
                else  ! already Picard
!!!                   print*, 'CONVERGENCE is failing with Picard; something is wrong'
!!!                   stop
                endif
             endif

             L2_previous = L2_norm

          endif   ! mod(iter_melt, n_check_convergence) = 0

          !WHL TODO - Another possibility: Adjust dt based on the iteration count.
          ! If no convergence after a certain number of iterations, then try dt -> dt/2, and start over.
          ! Keep reducing dt until we converge.
          ! Then the question is whether we can increase dt again the next time.
 
          ! The solution has not yet converged; solve the continuity equation for the new D_plume.

          print*, ' '
          print*, 'Compute plume thickness, dt_plume, time (s) =', dt_plume, time

          ! Solve the continuity equation for D_plume, given u_plume, v_plume, entrainment and detrainment.

          call compute_plume_thickness(&
               nx,       ny,     &
               dx,       dy,     &
               dt_plume,         &
               plume_mask_cell,  &
               u_plume_east,     &
               v_plume_east,     &
               u_plume_north,    &
               v_plume_north,    &
               plume_speed_east, &
               plume_speed_north,&
               edge_mask_east,   &
               edge_mask_north,  &
               divu_mask_east,   &
               divu_mask_north,  &
               H_cavity,         &
               entrainment,      &
               detrainment,      &
               itest,    jtest,  &
               rtest,            &
               iter_melt,        &  !WHL - debug
               D_plume_old,      &
               D_plume)

          ! halo updates
          call parallel_halo(D_plume)

          !WHL - some temporary diagnostics
          if (verbose_continuity) then

             print*, ' '
             print*, 'entrainment (m/s), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(e12.3)',advance='no') entrainment(i,j)
                enddo
                write(6,*) ' '
             enddo
             
             print*, ' '
             print*, 'detrainment (m/s), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(e12.3)',advance='no') detrainment(i,j)
                enddo
                write(6,*) ' '
             enddo
             
             print*, ' '
             print*, 'divergence (m/s), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(e12.3)',advance='no') divDu_plume(i,j)
                enddo
                write(6,*) ' '
             enddo
             
             print*, ' '
             print*, 'dD_dt (m/s), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(e12.3)',advance='no') (D_plume(i,j) - D_plume_old(i,j)) / dt_plume
                enddo
                write(6,*) ' '
             enddo
             
             print*, ' '
             print*, 'residual (m/s), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   if (plume_mask_cell(i,j) == 1) then
                      dD_dt = (D_plume(i,j) - D_plume_old(i,j)) / dt_plume
                      write(6,'(e12.3)',advance='no') dD_dt - entrainment(i,j) + detrainment(i,j) + divDu_plume(i,j)
                   else
                      write(6,'(e12.3)',advance='no') 0.0d0
                   endif
                enddo
                write(6,*) ' '
             enddo
             
!          print*, ' '
!          print*, 'plume_speed (m/s), rank =', rtest
!          do j = jtest+3, jtest-3, -1
!             do i = itest-3, itest+3
!                write(6,'(e12.3)',advance='no') plume_speed(i,j)
!             enddo
!             write(6,*) ' '
!          enddo

             print*, ' '
             print*, 'u_plume_east (m/s), rank =', rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(a6)',advance='no') '      '
                do i = itest-3, itest+3
                   write(6,'(f12.7)',advance='no') u_plume_east(i,j)
                enddo
                write(6,*) ' '
             enddo
             
             print*, ' '
             print*, 'v_plume_north (m/s), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f12.7)',advance='no') v_plume_north(i,j)
                enddo
                write(6,*) ' '
             enddo
             
             print*, ' '
             print*, 'H_cavity (m), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f12.6)',advance='no') H_cavity(i,j)
                enddo
                write(6,*) ' '
             enddo

             print*, ' '
             print*, 'Old D_plume (m), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f12.6)',advance='no') D_plume_old(i,j)
                enddo
                write(6,*) ' '
             enddo
             
             print*, ' '
             print*, 'New D_plume (m), rank =', rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f12.6)',advance='no') D_plume(i,j)
                enddo
                write(6,*) ' '
             enddo
             
             if (free_surface) then
                print*, ' '
                print*, 'New eta_plume (m), rank =', rtest
                do j = jtest+3, jtest-3, -1
                   do i = itest-3, itest+3
!!                   write(6,'(e14.7)',advance='no') eta_plume(i,j)
                      write(6,'(f14.9)',advance='no') eta_plume(i,j)
                   enddo
                   write(6,*) ' '
                enddo
             endif

          endif  ! verbose_continuity
          
          ! Compute the basal melt rate, temperature and salinity at the plume-ice interface,
          ! given the plume velocity and entrainment rate.

          call compute_plume_melt_rate(&
               nx,         ny,      &
               gammaT,              &
               gammaS,              &
               plume_mask_cell,     &
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
          call parallel_halo(T_plume)
          call parallel_halo(S_plume)

          if (verbose_melt) then

             if (this_rank == rtest) then
                i = itest
                j = jtest
                print*, ' '
                print*, 'Computed melt, rank, i, j, bmlt_float (m/yr) =', this_rank, i, j, bmlt_float(i,j)*scyr
                print*, 'T_b, S_b =', T_basal(i,j), S_basal(i,j)
                print*, 'T_p, S_p =', T_plume(i,j), S_plume(i,j)
             endif

             print*, ' '
             print*, 'New bmlt_float (m/yr), i, j, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f12.6)',advance='no') bmlt_float(i,j)*scyr
                enddo
                write(6,*) ' '
             enddo

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
             
          endif  ! verbose_melt
          
          ! Relax T and S toward solution
          !TODO - Not necessary to relax T_basal and S_basal?  Only T_plume and S_plume needed for drho_plume.
          bmlt_float(:,:) = (1.0d0 - relax_TS)*bmlt_float_latest(:,:) + relax_TS*bmlt_float(:,:)
          S_plume(:,:) = (1.0d0 - relax_TS)*S_plume_latest(:,:) + relax_TS*S_plume(:,:)
          T_plume(:,:) = (1.0d0 - relax_TS)*T_plume_latest(:,:) + relax_TS*T_plume(:,:)
          S_basal(:,:) = (1.0d0 - relax_TS)*S_basal_latest(:,:) + relax_TS*S_basal(:,:)
          T_basal(:,:) = (1.0d0 - relax_TS)*T_basal_latest(:,:) + relax_TS*T_basal(:,:)
!      !TODO - Use freezing relation instead of relaxation parameter?
       !!     T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)

          ! check convergence of melt rate in all grid cells

          converged_melt = .false.

          if (iter_melt > 1) then

             err_melt = 0.d0
          
             do j = 1, ny
                do i = 1, nx
                   if (abs(bmlt_float(i,j) - bmlt_float_latest(i,j)) > err_melt) then
                      err_melt = abs(bmlt_float(i,j) - bmlt_float_latest(i,j))
                      imax = i
                      jmax = j
                   endif
                enddo
             enddo
             
             if (err_melt > maxerr_melt) then
          
                print*, ' '
                print*, 'Melt rate has NOT CONVERGED:'
                print*, '   iter, time(s), rank, i, j, m_latest, m, err target, errmax (m/yr) =', &
                     iter_melt, time, this_rank, imax, jmax,&
                     bmlt_float_latest(imax,jmax)*scyr, bmlt_float(imax,jmax)*scyr, maxerr_melt*scyr, err_melt*scyr
                
             else  ! converged
          
                print*, ' '
                print*, 'Melt rate has CONVERGED everywhere'
                print*, '  iter, time(s), rank, i, j, m_latest, m, err target, errmax (m/yr) =', &
                     iter_melt, time, this_rank, imax, jmax, &
                     bmlt_float_latest(imax,jmax)*scyr, bmlt_float(imax,jmax)*scyr, maxerr_melt*scyr, err_melt*scyr

                converged_melt = .true.

             endif

          endif  ! iter_melt > 1

       enddo  ! iter_melt

       ! Increment the time
       time = time + dt_plume

       !WHL - debug
       if (main_task) then
          print*, ' '
          print*, 'Plume velocity, thickness and melt rate CONVERGED:'
          print*, 'Total time (s) =', time
       endif
       
       ! Interpolate the plume speed to cell centers (for diagnostics).

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             u_plume(i,j) = (u_plume_east(i,j) + u_plume_east(i-1,j)) / 2.0d0
             v_plume(i,j) = (v_plume_north(i,j) + v_plume_north(i,j-1)) / 2.0d0
             plume_speed(i,j) = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)
          enddo
       enddo
       
       ! Plume diagnostics

       if (verbose_continuity) then

          ! Find location of maximum plume thickness tendency

          max_tendency = 0.0d0

          do j = nhalo+1, ny-nhalo
             do i = nhalo+1, nx-nhalo
                if (plume_mask_cell(i,j) == 1) then

                   ! Check for negative D_plume. This should not happen with an adaptive time step.
                   if (D_plume(i,j) < 0.0d0) then
                      print*, 'ERROR: Exceeded CFL for plume adjustment:', i, j
                      print*, 'rank, i, j, D_plume, correction:', this_rank, i, j, D_plume(i,j), D_plume(i,j) - D_plume_old(i,j)
                      stop
                   endif
                   
                   ! Keep track of the maximum tendency.
                   ! We are trying to drive the tendency to a small value everywhere in the cavity.
                   plume_tendency = (D_plume(i,j) - D_plume_old(i,j)) / dt_plume
                   if (abs(plume_tendency) > max_tendency) then
                      max_tendency = abs(plume_tendency)
                      imax = i
                      jmax = j
                   endif

                   if (this_rank == rtest .and. i==itest .and. j==jtest) then
                      print*, ' '
                      print*, 'i, j, oldD, new D:', i, j, D_plume_old(i,j), D_plume(i,j)
                   endif

                endif  ! plume_mask_cell
             enddo  ! i
          enddo  ! j
          
          !TODO - Add global maxval for max_tendency

          ! print location of max change in D_plume
          print*, ' '
          print*, 'i, j, D_plume, dD, max tendency (m/s):', imax, jmax, D_plume(imax,jmax),&
               D_plume(imax,jmax) - D_plume_old(imax,jmax), max_tendency

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
          print*, 'i, j, max plume speed (m/s):', imax, jmax, plume_speed(imax,jmax)

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

          if (free_surface) then

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

          endif   ! free_surface

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
       
       endif   ! verbose_continuity


    enddo   ! iter_Dplume
       
       
    ! Compute the final value of eta.
!    where (plume_mask_cell == 1)
!       eta_plume = max(D_plume - (lsrf - topg), 0.0d0)
!    elsewhere
!       eta_plume = 0.0d0
!    endwhere

    !--------------------------------------------------------------------
    ! Compute various diagnostic quantities.
    !--------------------------------------------------------------------

    ! Copy u_plume_east and v_plume_north into u_plume_Cgrid and v_plume_Cgrid for output.
    ! Note: v_plume_east and u_plume_north are used internally but are not part of output.
    !TODO - Eliminate the duplicate arrays?
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

  subroutine compute_edge_gradients(&
       nx,              ny,          &
       dx,              dy,          &
       global_bndy_east,             &
       global_bndy_west,             &
       global_bndy_north,            &
       global_bndy_south,            &
       plume_mask_cell,              &
       floating_mask,                &
       lsrf,                         &
       field,                        &
       df_dx_east,      df_dy_east,  &
       df_dx_north,     df_dy_north)
   
    ! Compute the gradients of a scalar field on east and north cell edges.
    ! The procedure for east edges as follows:
    ! (1) Initialize all gradients to zero.
    ! (2) If the plume exists on both sides of an east edge, compute df/dx in the standard way.
    !     Similarly, if the plume exists on both sides of a north edge, compute df/dy in the standard way.
    ! (3) If the edge has a plume cell on one side and floating ice or open water on the other,
    !     and it is not a global boundary edge, then extrapolate the gradient from an adjacent edge.
    ! (4) Compute df/dy on east edges by averaging from adjacent north edges, and compute
    !     df/dx on north edges by extrapolating from adjacent east edges.
    
    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)
    
    integer, dimension(nx,ny), intent(in) ::  &
         global_bndy_east,   & ! = 1 for edges at global boundaries, else = 0
         global_bndy_west,   &
         global_bndy_north,  &
         global_bndy_south,  &
         plume_mask_cell,    & ! = 1 for cells where scalar plume variables are computed
         floating_mask         ! = 1 where ice is present and floating, else = 0
    
    real(dp), dimension(nx,ny), intent(in) ::  &
         lsrf                  ! lower ice surface (m); used to diagnose open ocean
    
    
    real(dp), dimension(nx,ny), intent(in) :: &
         field                 ! scalar field
    
    real(dp), dimension(nx,ny), intent(out) :: &
         df_dx_east,  df_dy_east,   &  ! gradient components on east edges
         df_dx_north, df_dy_north      ! gradient component on north edges
    
    ! local variables

    integer :: i, j

    ! initialize
    df_dx_east(:,:) = 0.0d0
    df_dy_east(:,:) = 0.0d0
   
    df_dx_north(:,:) = 0.0d0
    df_dy_north(:,:) = 0.0d0
   
    ! Compute gradients at edges with plume cells on each side

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! east edges
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 1) then
             df_dx_east(i,j) = (field(i+1,j) - field(i,j)) / dx
          endif

          ! north edges
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 1) then
             df_dy_north(i,j) = (field(i,j+1) - field(i,j)) / dy
          endif

       enddo
    enddo

    ! Set gradients at edges that have a plume cell on one side and floating ice or water on the other.
    ! Extrapolate the gradient from the nearest neighbor edge.
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! east edges
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 0 .and. global_bndy_east(i,j) == 0) then
             if (lsrf(i+1,j) == 0.0d0 .or. floating_mask(i+1,j) == 1) then
                df_dx_east(i,j) = df_dx_east(i-1,j)
             endif
          endif
          if (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i+1,j) == 1 .and. global_bndy_west(i,j) == 0) then
             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                df_dx_east(i,j) = df_dx_east(i+1,j)
             endif
          endif

          ! north edges
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 0 .and. global_bndy_north(i,j) == 0) then
             if (lsrf(i,j+1) == 0.0d0 .or. floating_mask(i,j+1) == 1) then
                df_dy_north(i,j) = df_dy_north(i,j-1)
             endif
          endif
          if (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i,j+1) == 1 .and. global_bndy_south(i,j) == 0) then
             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                df_dy_north(i,j) = df_dy_north(i,j+1)
             endif
          endif

       enddo
    enddo

    ! Average over 4 neighboring edges to estimate the y derivative on east edges and the x derivative on north edges.

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! y derivative on east edges
          df_dy_east(i,j) = 0.25d0 * (df_dy_north(i,j)   + df_dy_north(i+1,j)  &
                                    + df_dy_north(i,j-1) + df_dy_north(i+1,j-1))

          ! x derivative on north edges
          df_dx_north(i,j) = 0.25d0 * (df_dx_east(i-1,j+1) + df_dx_east(i,j+1)  &
                                     + df_dx_east(i-1,j)   + df_dx_east(i,j))

       enddo
    enddo

    !TODO - Add a halo update for parallel runs

  end subroutine compute_edge_gradients

!****************************************************

  subroutine compute_plume_velocity(&
       nx,       ny,        &
       dx,       dy,        &
       itest, jtest, rtest, &
       plume_mask_cell,     &
!!       floating_mask,       &
!!       global_bndy_east,    &
!!       global_bndy_west,    &
!!       global_bndy_north,   &
!!       global_bndy_south,   &
       divu_mask_east,      &
       divu_mask_north,     &
       edge_mask_east,      &
       edge_mask_north,     &
       edge_mask_east_reduce_v,  &
       edge_mask_north_reduce_u, &
       free_surface,        &
!!       lsrf,                &
       dlsrf_dx_east,      dlsrf_dy_east,  &
       dlsrf_dx_north,     dlsrf_dy_north, &
       drho_plume,          &
       D_plume,             &
       eta_plume,           &
       H_cavity,            &
       u_plume_east,        &
       v_plume_east,        &
       u_plume_north,       &
       v_plume_north,       &
       plume_speed_east,    &
       plume_speed_north)

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask_cell,       & ! = 1 for cells where scalar plume variables are computed
!!         floating_mask,         & ! = 1 where ice is present and floating, else = 0
!!         global_bndy_east,      & ! = 1 along east global boundary, else = 0
!!         global_bndy_west,      & ! = 1 along west global boundary, else = 0
!!         global_bndy_north,     & ! = 1 along north global boundary, else = 0
!!         global_bndy_south,     & ! = 1 along south global boundary, else = 0
         edge_mask_east,        & ! = 1 on east edges where plume velocity is computed
         edge_mask_north,       & ! = 1 on north edges where plume velocity is computed
         divu_mask_east,        & ! = 1 on east edges where divergence terms are computed
         divu_mask_north          ! = 1 on north edges where divergence terms are computed
                           
    real(dp), dimension(nx,ny), intent(in) ::  &
         edge_mask_east_reduce_v,  & ! mask for reducing v on east edges adjacent to a wall
         edge_mask_north_reduce_u    ! mask for reducing u on north edges adjacent to a wall

    logical, intent(in) ::  &
         free_surface           ! true if computing PG force due to slope in free surface

    real(dp), dimension(nx,ny), intent(in) ::  &
!!         lsrf,                & ! ice lower surface elevation (m, negative below sea level)
         dlsrf_dx_east,       & ! horizontal gradient of lsrf on east edges
         dlsrf_dy_east,       & !
         dlsrf_dx_north,      & ! horizontal gradient of lsrf on north edges
         dlsrf_dy_north

    !WHL - intent(inout) to allow temporary perturbations
    real(dp), dimension(nx,ny), intent(inout) ::  &
         D_plume                ! plume thickness (m)

    !WHL - Pass in eta or H_cavity but not both?
    real(dp), dimension(nx,ny), intent(in) ::  &
         drho_plume,          & ! density difference between plume and ambient ocean (kg/m^3)
         H_cavity               ! thickness of ocean cavity beneath the plume (m)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         eta_plume              ! displacement of plume surface, D_plume - H_cavity (m)
                                ! intent(inout) to allow perturbations
                                ! TODO - compute locally from D_plume - H_cavity

    real(dp), dimension(nx,ny), intent(inout) ::  &
         u_plume_east,        & ! u_plume on east edges
         v_plume_east,        & ! v_plume on east edges
         u_plume_north,       & ! u_plume on north edges
         v_plume_north,       & ! v_plume on north edges
         plume_speed_east,    & ! plume speed on east edges
         plume_speed_north      ! plume speed on north edges

    ! local variables

    real(dp) :: &
         deta_plume_dx,       & ! horizontal gradient of eta_plume
         deta_plume_dy

    real(dp), dimension(nx,ny) :: &
         pgf_x_east,          & ! x component of pressure gradient force on east edges (m^2/s^2)
         pgf_y_east,          & ! y component of pressure gradient force on east edges (m^2/s^2)
         pgf_x_north,         & ! x component of pressure gradient force on north edges (m^2/s^2)
         pgf_y_north            ! y component of pressure gradient force on north edges (m^2/s^2)

    real(dp), dimension(nx,ny) :: &
         latdrag_x_east,      & ! x component of lateral drag on east edges (m^2/s^2)
         latdrag_y_east,      & ! y component of lateral drag on east edges (m^2/s^2)
         latdrag_x_north,     & ! x component of lateral drag on north edges (m^2/s^2)
         latdrag_y_north        ! y component of lateral drag on north edges (m^2/s^2)

    real(dp), dimension(nx,ny) :: &
         D_plume_east,        & ! D_plume averaged to east edge
         D_plume_north,       & ! D_plume averaged to north edge
         grav_reduced_east,   & ! reduced gravity on east edge
         grav_reduced_north     ! reduced gravity on north edge

    integer :: i, j

    integer ::  &
         iter_velo              ! iteration counter

    character(len=100) :: message

    logical, dimension(nx,ny) ::  &
         converged_velo_east, & ! true when velocity has converged at an east edge, else false
         converged_velo_north   ! true when velocity has converged at a north edge, else false

    logical :: &
         converged_all_velo     ! true when velocity has converged at all edges, else false

    integer, parameter ::  &
         maxiter_velo = 100     ! max number of iterations of velocity loop

    ! initialize
    D_plume_east(:,:) = 0.0d0
    D_plume_north(:,:) = 0.0d0

    grav_reduced_east(:,:) = 0.0d0
    grav_reduced_north(:,:) = 0.0d0

    ! Note: There are a couple of different ways to compute the PGF.
    !       (1) Jenkins et al. (1991) and HJH (2008) use grad(lsrf)
    !       (2) Holland & Feltham (2006) use grad(lsrf_plume) along with a density gradient.
    !       Method (1) is simpler and has the advantage that grad(lsrf) does not vary during plume evolution,
    !        making the PGF more stable (though possibly not as accurate).

    ! Compute the pressure gradient force on each edge

    pgf_x_east(:,:) = 0.0d0
    pgf_y_east(:,:) = 0.0d0
    pgf_x_north(:,:) = 0.0d0
    pgf_y_north(:,:) = 0.0d0

    !TODO - Use edge_mask_east instead?
    !       Maybe divu_mask_east is correct, since I stopped extrapolating, but should be called edge_mask_east.
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! PGF on east edge
          if (divu_mask_east(i,j) == 1) then
             
             ! Compute horizontal pressure gradient force, not including the free-surface term.
             ! Based on HJH 2008
             ! On east edges, the x derivatives are based on values in the two adjacent cells,
             !  to preserve the advantages of a C grid.
             ! The y derivatives are averaged from the neighboring vertices.
             ! Note: D_plume = 0 and drho_plume = 0 where plume_mask_cell = 0.

             D_plume_east(i,j) = (D_plume(i,j) + D_plume(i+1,j)) / 2.0d0
             grav_reduced_east(i,j) = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i+1,j)) / 2.0d0
             
             pgf_x_east(i,j) = grav_reduced_east(i,j) * D_plume_east(i,j) * dlsrf_dx_east(i,j)
             pgf_y_east(i,j) = grav_reduced_east(i,j) * D_plume_east(i,j) * dlsrf_dy_east(i,j)
             
             !WHL - debug
             if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'grav_reduced_east = ', grav_reduced_east(i,j)
                print*, 'D_plume_east = ', D_plume_east(i,j)
                print*, 'xterm 1, yterm 1:', grav_reduced_east(i,j) * D_plume_east(i,j) * dlsrf_dx_east(i,j), &
                                             grav_reduced_east(i,j) * D_plume_east(i,j) * dlsrf_dy_east(i,j)
             endif
             
             ! Optionally, add the free-surface term
             ! TODO - Comment on treatment of deta/dy terms at plume boundary
             !TODO - Treat eta_plume gradients as lsrf gradients?  With an edge_gradient subroutine?
             !       Might not be necessary, since the subroutine does a special treatment where the plume
             !        borders open water or plume-free cells, to avoid sharp gradients.
             !        Should not have sharp gradients in eta, because eta = 0 outside the plume.
             if (free_surface) then
                
                deta_plume_dx = (eta_plume(i+1,j) - eta_plume(i,j)) / dx

                deta_plume_dy = 0.25d0/dy * ( (eta_plume(i,j+1)   - eta_plume(i,j))     * divu_mask_north(i,j)    &
                                            + (eta_plume(i,j)     - eta_plume(i,j-1))   * divu_mask_north(i,j-1)  &
                                            + (eta_plume(i+1,j+1) - eta_plume(i+1,j))   * divu_mask_north(i+1,j)  &
                                            + (eta_plume(i+1,j)   - eta_plume(i+1,j-1)) * divu_mask_north(i+1,j-1) )
!                deta_plume_dy = 0.25d0/dy * ( (eta_plume(i,j+1)   - eta_plume(i,j))     * edge_mask_north(i,j)    &
!                                            + (eta_plume(i,j)     - eta_plume(i,j-1))   * edge_mask_north(i,j-1)  &
!                                            + (eta_plume(i+1,j+1) - eta_plume(i+1,j))   * edge_mask_north(i+1,j)  &
!                                            + (eta_plume(i+1,j)   - eta_plume(i+1,j-1)) * edge_mask_north(i+1,j-1) )
                
                pgf_x_east(i,j) = pgf_x_east(i,j) - grav * D_plume_east(i,j) * deta_plume_dx
                pgf_y_east(i,j) = pgf_y_east(i,j) - grav * D_plume_east(i,j) * deta_plume_dy

                if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, 'deta_dx, deta_dy =', deta_plume_dx, deta_plume_dy
                   print*, 'xterm 2, yterm 2:', -grav * D_plume_east(i,j) * deta_plume_dx, &
                                                -grav * D_plume_east(i,j) * deta_plume_dy
                   print*, 'pgf x/y at east edge:', pgf_x_east(i,j), pgf_y_east(i,j)
                endif

             endif  ! free_surface
             
          endif   ! divu_mask_east


          ! PGF on north edge
          if (divu_mask_north(i,j) == 1) then
             
             ! Compute horizontal pressure gradient force, not including the free-surface term
             ! Based on HJH 2008
             ! On north edges, the y derivatives are based on values in the two adjacent cells,
             !  to preserve the advantages of a C grid.
             ! The x derivatives are averaged from the neighboring vertices.
             
             D_plume_north(i,j) = (D_plume(i,j) + D_plume(i,j+1)) / 2.0d0
             grav_reduced_north(i,j) = (grav/rhoo) * (drho_plume(i,j) + drho_plume(i,j+1)) / 2.0d0
             
             pgf_x_north(i,j) = grav_reduced_north(i,j) * D_plume_north(i,j) * dlsrf_dx_north(i,j)
             pgf_y_north(i,j) = grav_reduced_north(i,j) * D_plume_north(i,j) * dlsrf_dy_north(i,j)

             !WHL - debug
             if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'grav_reduced_north(i,j) = ', grav_reduced_north(i,j)
                print*, 'D_plume_north(i,j) = ', D_plume_north(i,j)
                print*, 'xterm 1, yterm 1:', grav_reduced_north(i,j) * D_plume_north(i,j) * dlsrf_dx_north(i,j), &
                                             grav_reduced_north(i,j) * D_plume_north(i,j) * dlsrf_dy_north(i,j)
             endif

             ! Optionally, add the free-surface term

             if (free_surface) then
                
                deta_plume_dx = 0.25d0/dx * ( (eta_plume(i,j+1)   - eta_plume(i-1,j+1)) * divu_mask_east(i-1,j+1) &
                                            + (eta_plume(i+1,j+1) - eta_plume(i,j+1))   * divu_mask_east(i,j+1)   &
                                            + (eta_plume(i,j)     - eta_plume(i,j-1))   * divu_mask_east(i-1,j)   &
                                            + (eta_plume(i+1,j)   - eta_plume(i,j))     * divu_mask_east(i,j) )
!                deta_plume_dx = 0.25d0/dx * ( (eta_plume(i,j+1)   - eta_plume(i-1,j+1)) * edge_mask_east(i-1,j+1) &
!                                            + (eta_plume(i+1,j+1) - eta_plume(i,j+1))   * edge_mask_east(i,j+1)   &
!                                            + (eta_plume(i,j)     - eta_plume(i,j-1))   * edge_mask_east(i-1,j)   &
!                                            + (eta_plume(i+1,j)   - eta_plume(i,j))     * edge_mask_east(i,j) )

                deta_plume_dy = (eta_plume(i,j+1) - eta_plume(i,j)) / dy
                
                pgf_x_north(i,j) = pgf_x_north(i,j) - grav * D_plume_north(i,j) * deta_plume_dx
                pgf_y_north(i,j) = pgf_y_north(i,j) - grav * D_plume_north(i,j) * deta_plume_dy

                if (verbose_velo .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'deta_dx, deta_dy:', deta_plume_dx, deta_plume_dy
                   print*, 'xterm 2, yterm 2:', -grav * D_plume_north(i,j) * deta_plume_dx, &
                                                -grav * D_plume_north(i,j) * deta_plume_dy
                   print*, 'pgf x/y at north edge:', pgf_x_north(i,j), pgf_y_north(i,j)
                endif
                
             endif
             
          endif   ! divu_mask_north

       enddo  ! i
    enddo  ! j

    ! initialize other fields
    latdrag_x_east(:,:) = 0.0d0
    latdrag_y_east(:,:) = 0.0d0
    latdrag_x_north(:,:) = 0.0d0
    latdrag_y_north(:,:) = 0.0d0

    converged_velo_east(:,:) = .false.
    converged_velo_north(:,:) = .false.

    ! Iterate as needed to compute a converged velocity at each edge

    do iter_velo = 1, maxiter_velo

       !WHL - debug
       print*, ' '
       print*, 'iter_velo =', iter_velo

       ! Compute velocity on east edges

       if (main_task) then
          print*, ' '
          print*, 'compute east edge velocities: r, i, j =', rtest, itest, jtest
       endif

       call compute_velocity(&
            nx,    ny,               &
            itest, jtest, rtest,     &   ! diagnostic only
            edge_mask_east,          &
            D_plume_east,            &
            pgf_x_east,              &
            pgf_y_east,              &
            latdrag_x_east,          &
            latdrag_y_east,          &
            u_plume_east,            &
            v_plume_east,            &
            converged_velo_east,     &
            edge_mask_east_reduce_v = edge_mask_east_reduce_v)

       ! Compute velocity on north edges

       if (main_task) then
          print*, ' '
          print*, 'compute north edge velocities'
       endif
    
       call compute_velocity(&
            nx,    ny,                  &
            itest, jtest, rtest,        &   ! diagnostic only
            edge_mask_north,            &
            D_plume_north,              &
            pgf_x_north,                &
            pgf_y_north,                &
            latdrag_x_north,            &
            latdrag_y_north,            &
            u_plume_north,              &
            v_plume_north,              &
            converged_velo_north,       &
            edge_mask_north_reduce_u = edge_mask_north_reduce_u)

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
          exit   ! iter_velo loop
       elseif (iter_velo == maxiter_velo) then
          write(message,*) 'Error, glissade_plume: velocity has not converged, iter_velo =', iter_velo
          call write_log(message, GM_FATAL)
       endif
       
    enddo  ! iter_velo
    
    !TODO - Now that the velocity has converged without lateral drag, try adding the lateral drag
    !       terms and recomputing the velocity. Not sure how to do this stably.

       !WHL - debug
       if (verbose_velo .and. main_task .and. this_rank==rtest) then
          print*, ' '
          print*, 'Computed new velocity'
          print*, ' '

          print*, ' '
          print*, 'pgf_x_east:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') pgf_x_east(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'pgf_y_east:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') pgf_y_east(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'pgf_x_north:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') pgf_x_north(i,j)
             enddo
             write(6,*) ' '
          enddo
       
          print*, ' '
          print*, 'pgf_y_north:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') pgf_y_north(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'u_plume_east:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') u_plume_east(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'v_plume_east:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') v_plume_east(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'u_plume_north:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') u_plume_north(i,j)
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'v_plume_north:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') v_plume_north(i,j)
             enddo
             write(6,*) ' '
          enddo

       endif          

       ! Compute the lateral drag term based on the current guess for the velocity

       call compute_lateral_drag(&
            nx,         ny,      &
            dx,         dy,      &
            itest, jtest, rtest, &
            edge_mask_east,      &  !TODO - divu_mask or edge_mask?
            edge_mask_north,     &
            plume_mask_cell,     &
            D_plume,             &
            u_plume_east,        &
            v_plume_east,        &
            u_plume_north,       &
            v_plume_north,       &
            latdrag_x_east,      &
            latdrag_y_east,      &
            latdrag_x_north,     &
            latdrag_y_north)

       if (verbose_velo .and. main_task .and. this_rank==rtest) then

          print*, ' '
          print*, 'Computed lateral drag terms'
          print*, ' '
          
          print*, ' '
          print*, 'latdrag_x_east:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') latdrag_x_east(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'latdrag_y_east:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') latdrag_y_east(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'latdrag_x_north:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') latdrag_x_north(i,j)
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'latdrag_y_north:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f12.6)',advance='no') latdrag_y_north(i,j)
             enddo
             write(6,*) ' '
          enddo
          
       endif  ! verbose_velo


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

    ! Compute the plume speed at the edges (including the u_tidal term)
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          plume_speed_east(i,j) = sqrt(u_plume_east(i,j)**2 + v_plume_east(i,j)**2 + u_tidal**2)
          plume_speed_north(i,j) = sqrt(u_plume_north(i,j)**2 + v_plume_north(i,j)**2 + u_tidal**2)
       enddo   ! i
    enddo   ! j

  end subroutine compute_plume_velocity

!****************************************************

  subroutine compute_velocity(&
       nx,    ny,               &
       itest, jtest, rtest,     &
       edge_mask,               &
       D_plume,                 &
       pgf_x,                   &
       pgf_y,                   &
       latdrag_x,               &
       latdrag_y,               &
       u_plume,                 &
       v_plume,                 &
       converged_velo,          &
       edge_mask_east_reduce_v, &
       edge_mask_north_reduce_u)
    
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
         pgf_y,             & ! y component of pressure gradient force
         latdrag_x,         & ! x component of lateral drag
         latdrag_y            ! y component of lateral drag
    
    real(dp), dimension(nx,ny), intent(inout) ::  &
         u_plume,           & ! x component of plume velocity (m/s)
         v_plume              ! x component of plume velocity (m/s)

    logical, dimension(nx,ny), intent(inout) ::  &
         converged_velo        ! true when velocity has converged at an edge, else false

    !TODO - Remove these terms if lateral drag works
    real(dp), dimension(nx,ny), intent(in), optional :: &
         edge_mask_east_reduce_v,  & ! mask for reducing v on east edges adjacent to a wall
         edge_mask_north_reduce_u    ! mask for reducing u on north edges adjacent to a wall

    ! local variables

    real(dp), dimension(nx,ny) ::   &
         f_x,               &  ! pgf_x + latdrag_x 
         f_y                   ! pgf_y + latdrag_y 

    real(dp), dimension(nx,ny) ::  &
         reduce_v,          &  ! local version of edge_mask_east_reduce_v; no reduction by default
         reduce_u              ! local version of edge_mask_north_reduce_u; no reduction by default

    real(dp) :: &
         plume_speed,       & ! plume speed (m/s)
         x_resid, y_resid,  & ! residuals of momentum balance equations (m^2/s^2)
         denom,             & ! denominator
         a_uu, a_uv,        & ! coefficients for Newton solve
         a_vu, a_vv,        & !
         du, dv               ! change in u_plume and v_plume (m/s)
    
    character(len=100) :: message

    real(dp), parameter :: &
         maxresid_force_balance = 1.0d-8 ! max residual allowed in momentum balance equation (m^2/s^2)
    
    logical, parameter :: &
         velo_newton = .true.  ! if true, use Newton's method; if false, use Picard method

    integer :: i, j

    !TODO - Add lateral drag to the equations
    !       Can be handled numerically by combining with pgf in a single force term

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

    if (present(edge_mask_north_reduce_u)) then
       reduce_u(:,:) = edge_mask_north_reduce_u(:,:)
    else
       reduce_u(:,:) = 1.0d0  ! no reduction
    endif

    if (present(edge_mask_east_reduce_v)) then
       reduce_v(:,:) = edge_mask_east_reduce_v(:,:)
    else
       reduce_v(:,:) = 1.0d0  ! no reduction
    endif

    ! Combine PGF and lateral drag into one term
    f_x(:,:) = pgf_x(:,:) + latdrag_x(:,:)
    f_y(:,:) = pgf_y(:,:) + latdrag_y(:,:)

    ! Loop over edges of locally owned cells
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          if (edge_mask(i,j) == 1 .and. .not.converged_velo(i,j) ) then
       
             ! Compute plume speed based on current u and v
             plume_speed = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)
       
             ! Compute residual of the momentum balance
!               x_resid = pgf_x - c_drag*plume_speed*u_plume + f_coriolis*D_plume*v_plume
!               y_resid = pgf_y - c_drag*plume_speed*v_plume - f_coriolis*D_plume*u_plume
             x_resid = f_x(i,j) - c_drag*plume_speed*u_plume(i,j) + reduce_v(i,j)*f_coriolis*D_plume(i,j)*v_plume(i,j)
             y_resid = f_y(i,j) - c_drag*plume_speed*v_plume(i,j) - reduce_u(i,j)*f_coriolis*D_plume(i,j)*u_plume(i,j)

             ! check convergence of plume velocity

             if (abs(x_resid) < maxresid_force_balance .and. abs(y_resid) < maxresid_force_balance) then

                converged_velo(i,j) = .true.

                ! diagnostic print
                if (this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Velocity converged: u/v_plume (m/s):', u_plume(i,j), v_plume(i,j)
                endif

             endif

             if (.not.converged_velo(i,j)) then

                if (velo_newton) then
          
                   ! compute some coefficients for the Newton solve
                   a_uu = c_drag * (plume_speed + u_plume(i,j)**2/plume_speed)
                   a_vv = c_drag * (plume_speed + v_plume(i,j)**2/plume_speed)
                      
                   a_uv = c_drag * (u_plume(i,j)*v_plume(i,j))/plume_speed - reduce_v(i,j)*D_plume(i,j)*f_coriolis
                   a_vu = c_drag * (u_plume(i,j)*v_plume(i,j))/plume_speed + reduce_u(i,j)*D_plume(i,j)*f_coriolis
                   
                   ! compute du and dv
                   denom = a_uu*a_vv - a_uv*a_vu
                      
                   if (abs(denom) > 0.0d0) then
                      du = (a_vv*x_resid - a_uv*y_resid) / denom
                      dv = (a_uu*y_resid - a_vu*x_resid) / denom
                         
                      u_plume(i,j) = u_plume(i,j) + du
                      v_plume(i,j) = v_plume(i,j) + dv
                      
                   else  ! denom = 0.0
                      write(6,*) 'Error, glissade_plume: ill-posed Newton solve for velocity, rank, i, j:', this_rank, i, j
                      write(6,*) 'a_uu, a_vv, a_uv, a_vu =', a_uu, a_vv, a_uv, a_vu
                      write(message,*) 'Error, glissade_plume: ill-posed Newton solve for velocity, rank, i, j:', this_rank, i, j
                      call write_log(message, GM_FATAL)
                   endif
                      
                else  ! simpler Picard solve
          
                   denom = (c_drag*plume_speed)**2 + (D_plume(i,j)*f_coriolis)**2
                   u_plume = (c_drag*plume_speed*f_x(i,j) + reduce_v(i,j)*D_plume(i,j)*f_coriolis*f_y(i,j)) / denom
                   v_plume = (c_drag*plume_speed*f_y(i,j) - reduce_u(i,j)*D_plume(i,j)*f_coriolis*f_x(i,j)) / denom
          
                endif  ! Newton or Picard

             endif  ! .not.converged_velo

             if (verbose_velo .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'plume_speed (m/s) =', plume_speed
                print*, 'pgf_x, pgf_y:', pgf_x(i,j), pgf_y(i,j)
                print*, 'latdrag_x, latdrag_y:', latdrag_x(i,j), latdrag_y(i,j)
                print*, 'Dfv, -Dfu:', D_plume(i,j) * f_coriolis * v_plume(i,j), &
                                     -D_plume(i,j) * f_coriolis * u_plume(i,j)
                print*, 'dragu, dragv:', c_drag * plume_speed * u_plume(i,j), &
                                         c_drag * plume_speed * v_plume(i,j)
                print*, 'x/y residual:', x_resid, y_resid
                print*, 'new u/v_plume:', u_plume(i,j), v_plume(i,j)
             endif

          endif  ! edge_mask
       enddo  ! i
    enddo  ! j
    
  end subroutine compute_velocity

!****************************************************

  subroutine compute_lateral_drag(&
       nx,         ny,      &
       dx,         dy,      &
       itest, jtest, rtest, &
       edge_mask_east,      &
       edge_mask_north,     &
       plume_mask_cell,     &
       D_plume,             &
       u_plume_east,        &
       v_plume_east,        &
       u_plume_north,       &
       v_plume_north,       &
       latdrag_x_east,      &
       latdrag_y_east,      &
       latdrag_x_north,     &
       latdrag_y_north)

    ! Compute lateral drag using a simple Laplacian formulation.
    !
    ! The drag terms in the x and y momentum balance equations, respectively, are
    ! 
    !      d/dx(A*D*du/dx) + d/dy(A*D*du/dy)
    !      d/dx(A*D*dv/dx) + d/dy(A*D*dv/dy)
    !
    ! where A is a spatially uniform drag coefficient.
    !
    ! Assume free-slip and no-penetration BC.
    ! In other words, velocity components parallel to the wall are assumed to have zero gradient,
    !  whereas components perpendicular to the wall are set to zero.

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    integer, dimension(nx,ny), intent(in) ::  &
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed
         edge_mask_north,     & ! = 1 on north edges where plume velocity is computed
         plume_mask_cell        ! = 1 for cells where scalar plume variables are computed

    real(dp), dimension(nx,ny), intent(in) ::  &
         D_plume,             & ! plume thickness (m)
         u_plume_east,        & ! u component of plume velocity on east edges (m/s)
         v_plume_east,        & ! v component of plume velocity on east edges (m/s)
         u_plume_north,       & ! u component of plume velocity on north edges (m/s)
         v_plume_north          ! v component of plume velocity on north edges (m/s)

    real(dp), dimension(nx,ny), intent(out) ::  &
         latdrag_x_east,      & ! lateral drag in x velocity component on east edges (m^2/s^2)
         latdrag_y_east,      & ! lateral drag in y velocity component on east edges (m^2/s^2)
         latdrag_x_north,     & ! lateral drag in x velocity component on north edges (m^2/s^2)
         latdrag_y_north        ! lateral drag in y velocity component on north edges (m^2/s^2)

    ! local variables

    real(dp), dimension(nx-1,ny-1) ::  &
         D_plume_vertex         ! D_plume averaged to vertices

    integer :: i, j

    real(dp), parameter ::  &
         A_latdrag = 10.0d0    ! lateral drag coefficient (m^2/s)

    real(dp) :: term1, term2, dx_term, dy_term   ! various terms in discretization

    ! initialize
    latdrag_x_east(:,:) = 0.0d0
    latdrag_y_east(:,:) = 0.0d0
    latdrag_x_north(:,:) = 0.0d0
    latdrag_y_north(:,:) = 0.0d0
    
    ! Average D_plume to vertices

    do j = 1, ny-1
       do i = 1, nx-1
          D_plume_vertex(i,j) = 0.25d0 * &
               (D_plume(i,j+1)* plume_mask_cell(i,j+1) + D_plume(i+1,j+1) * plume_mask_cell(i+1,j+1)  &
              + D_plume(i,j)  * plume_mask_cell(i,j)   + D_plume(i+1,j)   * plume_mask_cell(i+1,j))
       enddo
    enddo

    ! Compute lateral drag terms

    ! Loop over edges of locally owned cells
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! lateral drag on east edges

          if (edge_mask_east(i,j) == 1) then

             ! first latdrag_x

             ! d/dx(D*du/dx)
             term1 = D_plume(i+1,j) * (u_plume_east(i+1,j) - u_plume_east(i,j)) / dx
             term2 = D_plume(i,j)   * (u_plume_east(i,j) - u_plume_east(i-1,j)) / dx
             dx_term = (term1 - term2) / dx

             ! d/dy(D*du/dy)
             ! Enforce free-slip BC by zeroing the velocity gradient at plume boundaries
             if (edge_mask_east(i,j+1) == 1) then
                term1 = D_plume_vertex(i,j)   * (u_plume_east(i,j+1) - u_plume_east(i,j)) / dy
             else
                term1 = 0.0d0
             endif
             if (edge_mask_east(i,j-1) == 1) then
                term2 = D_plume_vertex(i,j-1) * (u_plume_east(i,j) - u_plume_east(i,j-1)) / dy
             else
                term2 = 0.0d0
             endif
             dy_term = (term1 - term2) / dy

             latdrag_x_east(i,j) = A_latdrag * (dx_term + dy_term)

             ! then latdrag_y

             ! d/dx(D*dv/dx)
             ! Enforce free-slip BC by zeroing the velocity gradient at plume boundaries
             if (edge_mask_east(i+1,j) == 1) then
                term1 = D_plume(i+1,j) * (v_plume_east(i+1,j) - v_plume_east(i,j)) / dx
             else
                term1 = 0.0d0
             endif
             if (edge_mask_east(i-1,j) == 1) then
                term2 = D_plume(i,j)   * (v_plume_east(i,j) - v_plume_east(i-1,j)) / dx
             else
                term2 = 0.0d0
             endif
             dx_term = (term1 - term2) / dx

             ! d/dy(D*du/dy)
             term1 = D_plume_vertex(i,j)   * (v_plume_east(i,j+1) - v_plume_east(i,j)) / dy
             term2 = D_plume_vertex(i,j-1) * (v_plume_east(i,j) - v_plume_east(i,j-1)) / dy
             dy_term = (term1 - term2) / dy

             latdrag_y_east(i,j) = A_latdrag * (dx_term + dy_term)

          endif  ! edge_mask_east

          ! lateral drag on north edge

          if (edge_mask_north(i,j) == 1) then

             ! first latdrag_x

             ! d/dx(D*dv/dx)
             term1 = D_plume_vertex(i,j)   * (u_plume_north(i+1,j) - u_plume_north(i,j)) / dx
             term2 = D_plume_vertex(i-1,j) * (u_plume_north(i,j) - u_plume_north(i-1,j)) / dx
             dx_term = (term1 - term2) / dx

             ! d/dy(D*dv/dy)
             ! Enforce free-slip BC by zeroing the velocity gradient at plume boundaries
             if (edge_mask_north(i,j+1) == 1) then
                term1 = D_plume(i,j+1) * (u_plume_north(i,j+1) - u_plume_north(i,j)) / dy
             else
                term1 = 0.0d0
             endif
             if (edge_mask_north(i,j-1) == 1) then
                term2 = D_plume(i,j)   * (u_plume_north(i,j) - u_plume_north(i,j-1)) / dy
             else
                term2 = 0.0d0
             endif
             dy_term = (term1 - term2) / dy

             latdrag_x_north(i,j) = A_latdrag * (dx_term + dy_term)

             ! then latdrag_y

             ! d/dx(D*dv/dx)
             ! Enforce free-slip BC by zeroing the velocity gradient at plume boundaries
             if (edge_mask_north(i+1,j) == 1) then
                term1 = D_plume_vertex(i,j)   * (v_plume_north(i+1,j) - v_plume_north(i,j)) / dx
             else
                term1 = 0.0d0
             endif
             if (edge_mask_north(i-1,j) == 1) then
                term2 = D_plume_vertex(i-1,j) * (v_plume_north(i,j) - v_plume_north(i-1,j)) / dx
             else
                term2 = 0.0d0
             endif
             dx_term = (term1 - term2) / dx

             ! d/dy(D*dv/dy)
             term1 = D_plume(i,j+1) * (v_plume_north(i,j+1) - v_plume_north(i,j)) / dy
             term2 = D_plume(i,j)   * (v_plume_north(i,j) - v_plume_north(i,j-1)) / dy
             dy_term = (term1 - term2) / dy

             latdrag_y_north(i,j) = A_latdrag * (dx_term + dy_term)

          endif  ! edge_mask_north

       enddo  ! i
    enddo  ! j

  end subroutine compute_lateral_drag

!****************************************************

  subroutine compute_entrainment(&
       nx,         ny,      &
       dx,         dy,      &
       itest, jtest, rtest, &
       divu_mask_east,      &
       divu_mask_north,     &
       plume_mask_cell,     &
!!       floating_mask,       &
!!       global_bndy_east,    &
!!       global_bndy_west,    &
!!       global_bndy_north,   &
!!       global_bndy_south,   &
!!       lsrf,                &
!!       H_cavity,            &
!!       D_plume,             &
       theta_slope,         &
       u_plume_east,        &
       v_plume_north,       &
!!       plume_speed_east,    &
!!       plume_speed_north,   &
!!       dlsrf_dx_east,      dlsrf_dy_east,  &
!!       dlsrf_dx_north,     dlsrf_dy_north, &
       entrainment)

    !--------------------------------------------------------------------
    ! Compute entrainment as a function of the plume speed and the slope of the
    !  plume-ambient interface, following Bo Pederson (1980) and Jenkins (1991).
    ! Entrainment is computed at cell edges (where the slope is computed
    !  most naturally) and then interpolated to the ice grid.
    !--------------------------------------------------------------------

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    integer, dimension(nx,ny), intent(in) ::  &
         divu_mask_east,      & ! = 1 on east edges where plume velocity is computed
         divu_mask_north,     & ! = 1 on north edges where plume velocity is computed
         plume_mask_cell        ! = 1 for cells where scalar plume variables are computed
!!         plume_mask_cell,     & ! = 1 for cells where scalar plume variables are computed
!!         floating_mask,       & ! = 1 where ice is present and floating, else = 0
!!         global_bndy_east,    & ! = 1 along east global boundary, else = 0
!!         global_bndy_west,    & ! = 1 along west global boundary, else = 0
!!         global_bndy_north,   & ! = 1 along north global boundary, else = 0
!!         global_bndy_south      ! = 1 along south global boundary, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
!!         lsrf,                & ! elevation of lower ice surface (m, negative below sea level)
!!         H_cavity,            & ! ocean cavity thickness (m), lsrf - topg
!!         D_plume,             & ! plume thickness (m)
!!         plume_speed_east,    & ! plume speed on east edges
!!         plume_speed_north,   & ! plume speed on north edges
         u_plume_east,          & ! u component of plume velocity on east edges (m/s)
         v_plume_north,         & ! v component of plume velocity on north edges (m/s)
         theta_slope            ! basal slope angle at cell centers (rad)
!!         dlsrf_dx_east,       & ! horizontal gradient of lsrf on east edges
!!         dlsrf_dy_east,       & !
!!         dlsrf_dx_north,      & ! horizontal gradient of lsrf on north edges
!!         dlsrf_dy_north         !

    real(dp), dimension(nx,ny), intent(out) ::  &
         entrainment              ! entrainment at cell centers (m/s)

    ! local variables

!    real(dp), dimension(nx,ny) :: &
!         lsrf_plume,            & ! elevation of plume-ambient interface (m, negative below sea level)
!         D_plume_cap,           & ! min(D_plume, H_cavity)
!         dlsrf_plume_dx_east,   & ! horizontal gradient of lsrf_plume on east edges
!         dlsrf_plume_dy_east,   & !
!         dlsrf_plume_dx_north,  & ! horizontal gradient of lsrf_plume on north edges
!         dlsrf_plume_dy_north,  & !
!         entrainment_east,      & ! entrainment on east edges
!         entrainment_north        ! entrainment on north edges

    real(dp) :: &
         u_plume_cell,           & ! u_plume averaged to cell center (m/s)
         v_plume_cell,           & ! v_plume averaged to cell center (m/s)
         plume_speed_cell          ! plume speed at cell center (m/s)

!    real(dp) ::   &
!         slope,                 & ! slope of plume-ambient interface (unitless)
!         theta_slope              ! atan of slope of basal ice interface (rad)

    integer :: i, j
    
    ! entrainment parameters
    real(dp), parameter ::   &
!!         H0_cavity = 10.d0,          & ! cavity thickness (m) below which the entrainment gradually approaches zero
         E0 = 0.072d0                  ! entrainment coefficient (unitless)
                                       ! Bo Pederson (1980) suggests E0 = 0.072
                                       ! Jenkins (1991, JGR) uses 0.036 to compensate for lack of Coriolis in 1D model

    ! Compute the elevation of the plume-ambient interface
    ! Note: lsrf and lsrf_plume are negative below sea level.
    !       D_plume is capped such that lsrf - D_plume >= topg
    !WHL -  Still works if  D_plume is not capped?


!!    lsrf_plume(:,:) = lsrf(:,:) - D_plume(:,:)

    ! Compute the horizontal gradient of lsrf_plume on east and north edges.
    ! This gradient appears in the entrainment term.

!    call compute_edge_gradients(&
!       nx,                   ny,     &
!       dx,                   dy,     &
!       global_bndy_east,             &
!       global_bndy_west,             &
!       global_bndy_north,            &
!       global_bndy_south,            &
!       plume_mask_cell,              &
!       floating_mask,                &
!       lsrf,                         &
!       lsrf_plume,                   &
!       dlsrf_plume_dx_east,  dlsrf_plume_dy_east,  &
!       dlsrf_plume_dx_north, dlsrf_plume_dy_north)

    ! Compute entrainment at cell centers as a function of basal slope and plume speed

    entrainment(:,:) = 0.0d0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then
             u_plume_cell = 0.5d0 * (u_plume_east(i-1,j) + u_plume_east(i,j))
             v_plume_cell = 0.5d0 * (v_plume_north(i,j-1) + v_plume_north(i,j))
             plume_speed_cell = sqrt(u_plume_cell**2 + v_plume_cell**2)
             entrainment(i,j) = E0 * plume_speed_cell * sin(theta_slope(i,j))

             !WHL - debug
             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, 'i, j, theta, plume speed, entrainment:', i, j, theta_slope(i,j), plume_speed_cell, entrainment(i,j)
             endif
             
          endif
       enddo
    enddo

    !TODO - delete the following

    ! Compute entrainment on east edges

!    entrainment_east(:,:) = 0.0d0

!    do j = nhalo, ny-nhalo
!       do i = nhalo, nx-nhalo

!          if (divu_mask_east(i,j) == 1) then

!!!             slope = sqrt(dlsrf_plume_dx_east(i,j)**2 + dlsrf_plume_dy_east(i,j)**2)
!             slope = sqrt(dlsrf_dx_east(i,j)**2 + dlsrf_dy_east(i,j)**2)
!             theta_slope = atan(slope)
!!!             entrainment_east(i,j) = E0 * plume_speed_east(i,j) * sin(theta_slope) * tanh(H_cavity(i,j)/H0_cavity)
!             entrainment_east(i,j) = E0 * plume_speed_east(i,j) * sin(theta_slope)

             !WHL - debug
!             if (i==itest .and. j==jtest) then
!                print*, ' '
!                print*, 'i, j, entrainment_east:', i, j, entrainment_east(i,j)
!                print*, 'plume_speed_east:', plume_speed_east(i,j)
!!!                print*, 'lsrf(i,j), D_plume(i,j):', lsrf(i,j), D_plume(i,j)
!!!                print*, 'lsrf(i+1,j), D_plume(i+1,j):', lsrf(i+1,j), D_plume(i+1,j)
!!!                print*, 'lsrf_plume(i,j), lsrf_plume(i+1,j):', lsrf_plume(i,j), lsrf_plume(i+1,j)
!                print*, 'dlsrf_dx, dlsrf_dy, slope:', dlsrf_dx_east(i,j), dlsrf_dy_east(i,j), slope
!             endif

!          endif  ! divu_mask_east
!       enddo  ! i
!    enddo  ! j
    
    ! Compute entrainment on north edges

!    entrainment_north(:,:) = 0.0d0
    
!    do j = nhalo, ny-nhalo
!       do i = nhalo, nx-nhalo
!          if (divu_mask_north(i,j) == 1) then

!             slope = sqrt(dlsrf_dx_north(i,j)**2 + dlsrf_dy_north(i,j)**2)
!             theta_slope = atan(slope)
!!!             entrainment_north(i,j) = E0 * plume_speed_north(i,j) * sin(theta_slope) * tanh(H_cavity(i,j)/H0_cavity)
!             entrainment_north(i,j) = E0 * plume_speed_north(i,j) * sin(theta_slope)
          
             !WHL - debug
!             if (i==itest .and. j==jtest) then
!                print*, ' '
!                print*, 'i, j, entrainment_north:', i, j, entrainment_north(i,j)
!                print*, 'plume_speed_north:', plume_speed_north(i,j)
!!!                print*, 'lsrf(i,j), D_plume(i,j):', lsrf(i,j), D_plume(i,j)
!!!                print*, 'lsrf(i,j+1), D_plume(i,j+1):', lsrf(i,j+1), D_plume(i,j+1)
!!!                print*, 'lsrf_plume(i,j), lsrf_plume(i,j+1):', lsrf_plume(i,j), lsrf_plume(i,j+1)
!                print*, 'dlsrf_dx, dlsrf_dy, slope:', dlsrf_dx_north(i,j), dlsrf_dy_north(i,j), slope
!             endif

!          endif  ! divu_mask_north
!       enddo  ! i
!    enddo  ! j
    
    ! interpolate entrainment from edges to cell centers
    ! Note: I tried setting e = 0 when D_plume -> H_cavity.
    !       However, this leads to oscillations in bmlt_float when D_plume is close to H_cavity.

!    entrainment(:,:) = 0.0d0

!    do j = nhalo+1, ny-nhalo
!       do i = nhalo+1, nx-nhalo
!          if (plume_mask_cell(i,j) == 1) then
!             entrainment(i,j) = ( entrainment_east(i-1,j) + entrainment_east(i,j) +   &
!                                  entrainment_north(i,j-1) + entrainment_north(i,j) ) / 4.0d0 
!          endif
!       enddo
!    enddo
    
  end subroutine compute_entrainment

!****************************************************

  subroutine compute_detrainment(&
       nx,       ny,   &
       itest, jtest, rtest, &
       free_surface,   &
       dt_plume,       &
       H_cavity,       &
       D_plume,        &
       eta_plume,      &
       detrainment)

    ! Compute detrainment.
    ! This is not a physically based mechanism, just a regularization to prevent very thick plumes.
    ! Ideally, detrainment = 0 almost everywhere.

    integer, intent(in) ::  &
         nx,     ny           ! number of grid cells in each dimension

    integer, intent(in) ::  &
         itest, jtest, rtest  ! test cell coordinates (diagnostic only)

    logical, intent(in) :: &
         free_surface         ! true if computing PG force due to slope in free surface

    real(dp), intent(in) :: &
         dt_plume             ! time step (s)

    !WHL - Remove H_cavity?
    real(dp), dimension(nx,ny), intent(in) ::  &
         H_cavity,          & ! cavity thickness (m), lsrf - topg
         D_plume,           & ! plume thickess (m)
         eta_plume            ! displacement of plume surface, D_plume - H_cavity (m)

    real(dp), dimension(nx,ny), intent(out) ::  &
         detrainment          ! plume detrainment rate (m/s)

    ! local variables

    integer :: i, j

    ! detrainment parameters
    real(dp), parameter ::  &
         D_plume_max = 50.d0,        & ! plume thickness threshold (m) where detrainment begins
         tau_detrainment = 3600.d0      ! detrainment time scale (s)

    detrainment(:,:) = 0.0d0

    if (free_surface) then

       do j = 1, ny
          do i = 1, nx
             if (D_plume(i,j) > D_plume_max) then
                detrainment(i,j) = (D_plume(i,j) - D_plume_max) / tau_detrainment
             endif
          enddo
       enddo

    else    ! not a free surface; no grad(eta) term in PGF
    
       do j = 1, ny
          do i = 1, nx

                !WHL - Testing different options here.  Set detrainment to 0 if trying to converge the matrix.
!!             if (D_plume(i,j) > H_cavity(i,j)) then
!!                detrainment(i,j) = eta_plume(i,j) / tau_detrainment
!!                detrainment(i,j) = (D_plume(i,j) - H_cavity(i,j)) / dt_plume
!!                detrainment(i,j) = 0.0d0
!!             endif

             if (D_plume(i,j) > D_plume_max) then
                detrainment(i,j) = (D_plume(i,j) - D_plume_max) / tau_detrainment
             else
                detrainment(i,j) = 0.0d0
             endif

             if (i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, ' '
                print*, 'i, j, D, H_cavity, detrainment:', i, j, D_plume(i,j), H_cavity(i,j), detrainment(i,j)
             endif

          enddo
       enddo

    endif
    
  end subroutine compute_detrainment

!****************************************************

  subroutine compute_dynamic_residual(&
               nx,       ny,        &
               dx,       dy,        &
               dt_plume,            &
               itest, jtest, rtest, &
               plume_mask_cell,     &
               edge_mask_east,      &
               edge_mask_north,     &
               divu_mask_east,      &
               divu_mask_north,     &
               H_cavity,            &
               entrainment,         &
               detrainment,         &
               D_plume_old,         &
               D_plume,             &
               u_plume_east,        &
               v_plume_north,       &
               u_plume_north,       &  ! diagnostic only
               v_plume_east,        &  ! diagnostic only
               divDu_plume,         &
               L2_norm)

    ! Check for convergence of the continuity equation, dD/dt = e - d - del*(Du)

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    real(dp), intent(in) :: &
         dt_plume               ! time step (s)

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask_cell,     & ! = 1 for cells where scalar plume variables are computed
         edge_mask_east,      & ! = 1 on east edges with plume cells on each side
         edge_mask_north,     & ! = 1 on north edges with plume cels on each side
         divu_mask_east,      & ! = 1 on east edges where divergence terms are computed, else = 0
         divu_mask_north        ! = 1 on north edges where divergence terms are computed, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         H_cavity,            & ! ocean cavity thickness (m), lsrf - topg
         entrainment,         & ! entrainment at cell centers (m/s)
         detrainment,         & ! detrainment at cell centers (m/s)
         D_plume_old,         & ! old plume thickness from previous time step (m)
         D_plume,             & ! latest guess for plume thickness (m)
         u_plume_east,        & ! u_plume on east edges (m/s)
         v_plume_north          ! v_plume on north edges (m/s)

    ! diagnostic only
    real(dp), dimension(nx,ny), intent(in) ::  &
         u_plume_north,      & ! u_plume on north edges (m/s)
         v_plume_east          ! v_plume on east edges (m/s)

    real(dp), dimension(nx,ny), intent(out) ::  &
         divDu_plume            ! plume divergence, div(Du)
                                ! computed here and output as a diagnostic

    real(dp), intent(out) ::  &
         L2_norm                ! L2 norm of residual vector

    ! local variables

    real(dp), dimension(nx,ny) :: &
         D_plume_east_up,     & ! upstream plume thickness at each east edge (m)
         D_plume_north_up       ! upstream plume thickness at each north edge (m)

    real(dp), dimension(nx,ny) :: &
         D_plume_cap                      ! min(D_plume, H_cavity)

    real(dp) :: &
         D_plume_east, D_plume_west,    & ! terms in discretization of plume divergence
         D_plume_north, D_plume_south,  &
         dDu_dx, dDv_dy,                &
         dD_dt,                         & ! rate of change of D_plume (m/s)
         resid                            ! local residual (m/s)

    integer :: i, j

    !WHL - debug
    real(dp) :: max_resid
    integer :: imax, jmax
    max_resid = 0.0d0

    L2_norm = 0.0d0

    !WHLcap - Use capped value of D_plume

!    if (cap_Dplume) then
!       D_plume_cap(:,:) = min(D_plume(:,:), H_cavity(:,:))
!    else
!       !WHL - Or maybe not?
!       D_plume_cap(:,:) = D_plume(:,:)
!    endif

    !TODO - Inline the code for computing upstream-biased divergence?

    D_plume_east_up(:,:) = 0.0d0
    D_plume_north_up(:,:) = 0.0d0

    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             ! find the upstream value of D_plume on the east edge

             if (edge_mask_east(i,j) == 1) then ! plume exists in both neighbor cells

                if (u_plume_east(i,j) > 0.0d0) then
!                   upos_mask(i,j) = 1
                   D_plume_east_up(i,j) = D_plume(i,j)
                elseif (u_plume_east(i,j) < 0.0d0) then
!                   uneg_mask(i,j) = 1
                   D_plume_east_up(i,j) = D_plume(i+1,j)
                endif

             elseif (divu_mask_east(i,j) == 1) then  ! plume exists on only one side of the edge

                if (plume_mask_cell(i,j) == 1 .and. u_plume_east(i,j) > 0.0d0) then
!                   upos_mask(i,j) = 1
                   D_plume_east_up(i,j) = D_plume(i,j)
                elseif (plume_mask_cell(i+1,j) == 1 .and. u_plume_east(i,j) < 0.0d0) then
!                   uneg_mask(i,j) = 1
                   D_plume_east_up(i,j) = D_plume(i+1,j)
                endif
                
             endif


             ! find the upstream value of D_plume on the north edge

             if (edge_mask_north(i,j) == 1) then  ! plume exists in both neighbor cells

!!             D_plume_north(i,j) = min(D_plume(i,j), D_plume(i,j+1))

                if (v_plume_north(i,j) > 0.0d0) then
!                   vpos_mask(i,j) = 1
                   D_plume_north_up(i,j) = D_plume(i,j)
                elseif (v_plume_north(i,j) < 0.0d0) then
!                   vneg_mask(i,j) = 1
                   D_plume_north_up(i,j) = D_plume(i,j+1)
                endif
                
             elseif (divu_mask_north(i,j) == 1) then  ! plume exists on only one side of the edge
                
                if (plume_mask_cell(i,j) == 1 .and. v_plume_north(i,j) > 0.0d0) then
!                   vpos_mask(i,j) = 1
                   D_plume_north_up(i,j) = D_plume(i,j)
                elseif (plume_mask_cell(i,j+1) == 1 .and. v_plume_north(i,j) < 0.0d0) then
!                   vneg_mask(i,j) = 1
                   D_plume_north_up(i,j) = D_plume(i,j+1)
                endif
                
             endif

          endif  ! plume_mask_cell
       enddo  ! i
    enddo  ! j


    ! Compute the divergence in each cell
    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             dDu_dx = (D_plume_east_up(i,j)*u_plume_east(i,j) - D_plume_east_up(i-1,j)*u_plume_east(i-1,j)) / dx
             dDv_dy = (D_plume_north_up(i,j)*v_plume_north(i,j) - D_plume_north_up(i,j-1)*v_plume_north(i,j-1)) / dy

!             if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 1) then
!                D_plume_north = min(D_plume_cap(i,j), D_plume_cap(i,j+1))
!             elseif (plume_mask_cell(i,j) == 1) then
!                D_plume_north = D_plume_cap(i,j)
!             elseif (plume_mask_cell(i,j+1) == 1) then
!                D_plume_north = D_plume_cap(i,j+1)
!             endif
             
!             if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j-1) == 1) then
!                D_plume_south = min(D_plume_cap(i,j), D_plume_cap(i,j-1))
!             elseif (plume_mask_cell(i,j) == 1) then
!                D_plume_south = D_plume_cap(i,j)
!             elseif (plume_mask_cell(i,j-1) == 1) then
!                D_plume_south = D_plume_cap(i,j-1)
!             endif

!             if (divu_mask_north(i,j) == 1) then
!                D_plume_north = 0.5d0 * (D_plume_cap(i,j) + D_plume_cap(i,j+1))

!                if (v_plume_north(i,j) > 0.0d0) then
!!                   D_plume_north = D_plume(i,j)
!                   D_plume_north = D_plume_cap(i,j)
!                else
!!                   D_plume_north = D_plume(i,j+1)
!                   D_plume_north = D_plume_cap(i,j+1)
!                endif
!             else
!                D_plume_north = 0.0d0
!             endif

!             if (divu_mask_north(i,j-1) == 1) then
!                D_plume_south = 0.5d0 * (D_plume_cap(i,j) + D_plume_cap(i,j-1))

!                if (v_plume_north(i,j-1) > 0.0d0) then
!!                   D_plume_south = D_plume(i,j-1)
!                   D_plume_south = D_plume_cap(i,j-1)
!                else
!!                   D_plume_south = D_plume(i,j)
!                   D_plume_south = D_plume_cap(i,j)
!                endif
!             else
!                D_plume_south = 0.0d0
!             endif

!             dDv_dy = (D_plume_north*v_plume_north(i,j) - D_plume_south*v_plume_north(i,j-1)) / dy
             
             divDu_plume(i,j) = dDu_dx + dDv_dy
             
!!             if (D_plume(i,j) < H_cavity(i,j)) then
                dD_dt = (D_plume(i,j) - D_plume_old(i,j)) / dt_plume
                resid = dD_dt - entrainment(i,j) + detrainment(i,j) + divDu_plume(i,j)
!!             else  ! cavity is filled; solving for steady-state eta_plume instead of dD/dt
!!                dD_dt = 0.0d0
!!                resid = -entrainment(i,j) + detrainment(i,j) + divDu_plume(i,j)
!!             endif

             !WHL - debug
             if (i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, ' '
                print*, 'Divergence, i, j =', i, j
                print*, 'dD/dt =', dD_dt
                print*, 'u_plume_west/east =', u_plume_east(i-1,j), u_plume_east(i,j)
!!                print*, 'u_plume_south/north =', u_plume_north(i,j-1), u_plume_north(i,j)
!!                print*, 'v_plume_west/east =', v_plume_east(i-1,j), v_plume_east(i,j)
                print*, 'v_plume_south/north =', v_plume_north(i,j-1), v_plume_north(i,j)
!!                print*, 'D_plume(i-1,j),(i,j) =', D_plume(i-1,j), D_plume(i,j)
!!                print*, 'D_plume(i,j),(i,j+1) =', D_plume(i,j), D_plume(i,j+1)
                print*, 'D_plume_west/east =', D_plume_east_up(i-1,j), D_plume_east_up(i,j)
                print*, 'D_plume_south/north =', D_plume_north_up(i,j-1), D_plume_north_up(i,j)
                print*, 'dDu_dx, dDv_dy =', dDu_dx, dDv_dy
                print*, 'divDu =', divDu_plume(i,j)
                print*, 'local residual =', resid
             endif
             
             L2_norm = L2_norm + resid*resid
             
             ! WHL - debug
             if (abs(resid) > max_resid) then
                max_resid = abs(resid)
                imax = i
                jmax = j
             endif

          endif  ! plume_mask_cell = 1
       enddo   ! i
    enddo   ! j
    
    !WHL - debug
    print*, 'i, j max_resid:', imax, jmax, max_resid

    !TODO - Add a global sum for parallel code
    L2_norm = sqrt(L2_norm)

  end subroutine compute_dynamic_residual

   
!****************************************************

  subroutine compute_plume_thickness(&
       nx,       ny,     &
       dx,       dy,     &
       dt_plume,         &
       plume_mask_cell,  &
       u_plume_east,     &
       v_plume_east,     &
       u_plume_north,    &
       v_plume_north,    &
       plume_speed_east, &
       plume_speed_north,&
       edge_mask_east,   &  !WHL - Do we need both pairs of masks?
       edge_mask_north,  &
       divu_mask_east,   &
       divu_mask_north,  &
       H_cavity,         &
       entrainment,      &
       detrainment,      &
       itest,    jtest,  &
       rtest,            &
       iter_melt,        &  !WHL - debug
       D_plume_old,      &
       D_plume)

    !--------------------------------------------------------------------
    ! Solve the continuity equation for D_plume:
    ! 
    !    dD/dt = e - d - del*(Du)
    !
    ! where e is entrainment and d is detrainment.
    !  
    ! The equation is solved in delta form:
    !
    !    delta_D = (D_old - D_cur) + dt*(e-d) - dt*div(Du)
    !
   ! where D_cur is the current guess for D, passed into the subroutine as D_plume.
    ! where D_old is the value of D from the previous time step.
    ! 
    ! TODO - Update the method description.
    ! The current guess for the velocity (u_plume,v_plume) is also passed into the subroutine.
    ! Given u_plume and D_plume, the divergence term is expanded to first order as
    ! 
    !    div(Du) = div[(D_old + delta_D) * u_plume]
    !
    ! The terms containing delta_D are moved to the LHS and inserted in a matrix,
    ! giving a problem of the form
    !
    !    A*delta_D = rhs
    !
    ! This subroutine is called repeatedly until the residual is sufficiently small.
    !
    ! For now, I am using SLAP to solve the matrix.
    ! Later, I plan to use a homegrown parallel solver.
    !--------------------------------------------------------------------

    ! for sparse_easy_solve
    use glimmer_sparse_type
    use glimmer_sparse

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    real(dp), intent(in) ::  &
         dt_plume               ! time step for plume solver (s)

    !TODO - Only one pair of edge masks?
    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask_cell,     & ! = 1 for cells where scalar plume variables are computed, else = 0
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed, else = 0
         edge_mask_north,     & ! = 1 on north edges where plume velocity is computed, else = 0
         divu_mask_east,      & ! = 1 on east edges where divergence terms are computed, else = 0
         divu_mask_north        ! = 1 on north edges where divergence terms are  computed, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         u_plume_east,        & ! u_plume on east edges (m/s)
         v_plume_east,        & ! v_plume on east edges (m/s)
         u_plume_north,       & ! u_plume on north edges (m/s)
         v_plume_north,       & ! v_plume on north edges (m/s)
         plume_speed_east,    & ! plume speed on east edges (m/s)
         plume_speed_north,   & ! plume speed on north edges (m/s)
         H_cavity,            & ! ocean cavity thickness (m), lsrf - topg
         entrainment,         & ! entrainment at cell centers (m/s)
         detrainment            ! detrainment at cell centers (m/s)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    integer, intent(in) :: iter_melt  !WHL - debug

    real(dp), dimension(nx,ny), intent(in) ::  &
         D_plume_old            ! old plume thickness (m) at the previous time step

    real(dp), dimension(nx,ny), intent(inout) ::  &
         D_plume                ! on input, the latest guess for the plume thickness
                                ! on output, the new guess for the plume thickness

    ! local variables

    real(dp), dimension(nx,ny) ::  &
!!         D_plume_latest,      & ! D_plume from the most recent iteration (m)
!!         D_plume_cap,         & ! min(D_plume, H_cavity)
!!         D_plume_east,        & ! plume thickness at each east edge (m)
!!         D_plume_north,       & ! plume thickness at each north edge (m)
         D_plume_east_up,     & ! upstream plume thickness at each east edge (m)
         D_plume_north_up       ! upstream plume thickness at each north edge (m)

    integer, dimension(nx,ny) ::  &
         upos_mask,   & ! = 1 at edges where u_plume > 0, else = 0
         uneg_mask,   & ! = 1 at edges where u_plume < 0, else = 0
         vpos_mask,   & ! = 1 at edges where v_plume > 0, else = 0
         vneg_mask      ! = 1 at edges where v_plume < 0, else = 0

    real(dp), dimension(-1:1,-1:1,nx,ny) ::  &
         A_plume        ! array holding nonzero matrix elements on the structured mesh
                        ! up to 9 nonzero elements per row of the matrix

    type(sparse_matrix_type) ::  &
         matrix         ! sparse matrix for SLAP solver, defined in glimmer_sparse_types
                        ! includes nonzeroes, order, col, row, val

    real(dp), dimension(:), allocatable ::  &
         rhs,         & ! right-hand-side vector, passed to solver
         answer         ! answer vector, returned from solver

    real(dp) ::  &
         err            ! solution error the solver

    integer, dimension(nx,ny) ::  &
         cellID         ! integer ID for each cell

    integer, dimension(nx*ny) ::  &
         iCellIndex, jCellIndex   ! indices for mapping cellID back to i and j

!!    real(dp) :: &
!!         D_east, D_north,      & ! current estimate of D_plume, averaged to east and north edges
!!         denom                   ! denominator in the expression for velocity

    integer :: niters   ! iteration counter
    integer :: i, j     ! horizontal indices
    integer :: iA, jA   ! horizontal index shifts, in range -1:1
    integer :: n        ! matrix index
    integer :: count    ! counter
    integer :: matrix_order  ! size of square matrix
    integer :: max_nonzeros  ! max number of nonzero elements in matrix

    ! SLAP linear solver (BICG or GMRES)
    ! For ISOMIP+, BICG is a bit faster, requiring 3 or 4 linear iterations compared to 6 or 7 for GMRES.
    !TODO - Replace with homegrown solver

    integer, parameter :: &
         whichsparse = HO_SPARSE_BICG
!!         whichsparse = HO_SPARSE_GMRES

!!    real(dp), parameter :: &
!!         relax_D = 0.5d0       !WHL - Remove relax_D?

    real(dp) :: dDu_dx, dDv_dy

    !WHL - debug
!!    real(dp) :: diag_east, diag_west, diag_north, diag_south
!!    real(dp) :: offdiag_east, offdiag_west, offdiag_north, offdiag_south
!!    real(dp), dimension(nx,ny) :: eta_plume_latest

    print*, ' '
    print*, 'In plume_thickness_solver: itest, jtest =', itest, jtest

    ! count plume cells in matrix solve
    ! loop over locally owned cells
    count = 0
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j)==1) then
             count = count + 1
             cellID(i,j) = count
             iCellIndex(count) = i
             jCellIndex(count) = j
          endif
       enddo
    enddo
   
    ! initialize and allocate
    matrix_order = count
    max_nonzeros = matrix_order * 9

    allocate(matrix%row(max_nonzeros), matrix%col(max_nonzeros), matrix%val(max_nonzeros))
    allocate(rhs(matrix_order), answer(matrix_order))

    A_plume(:,:,:,:) = 0.0d0
    rhs(:) = 0.0d0
    answer(:) = 0.0d0

    ! Save the latest iterate for D_plume.
    ! Note: This is different from D_plume_old, the value at the start of the time step.
!!    D_plume_latest(:,:) = D_plume(:,:)

    !WHL - debug
    i = itest
    j = jtest
    print*, 'i, j, latest D =', i, j, D_plume(i,j)

    ! Compute D_plume edge values used in the finite difference expression for divergence.
    ! No solution is ideal here:
    ! With a centered difference for D_plume on edges, thin plumes upstream from thick plumes can overempty.
    ! With an upstream difference, there is no overemptying, but there can be oscillations as a result
    !  of u or v changing sign between one iteration and the next.
    ! Choosing the min value is not very accurate (only first order, and too restrictive of outflow when
    !  a thick plume is upstream of a thin plume), but it inhibits overemptying and oscillations.

    ! Store for each edge the plume thickness in the upstream cell.

    ! Note: The use of edge_mask and divu_mask is a bit subtle.
    !       If edge_mask = 1, the plume exists on both sides of the edge, and one cell is clearly upstream.
    !       If edge_mask = 0 but divu_mask = 1, the plume exists on only one side of the edge.
    !        We identify an upstream plume thickness only if the plume exists in the upstream cell
    !        and the flow is out of the plume domain.
    !        (This generally is the case, but we check here to be sure.)

    upos_mask(:,:) = 0
    uneg_mask(:,:) = 0
    vpos_mask(:,:) = 0
    vneg_mask(:,:) = 0

!!    D_plume_east(:,:) = 0.0d0
!!    D_plume_north(:,:) = 0.0d0

    D_plume_east_up(:,:) = 0.0d0
    D_plume_north_up(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx

          ! mark the u component as positive or negative

          if (edge_mask_east(i,j) == 1) then ! plume exists in both neighbor cells

!!             D_plume_east(i,j) = min(D_plume(i,j), D_plume(i+1,j))

             if (u_plume_east(i,j) > 0.0d0) then
                upos_mask(i,j) = 1
                D_plume_east_up(i,j) = D_plume(i,j)
             elseif (u_plume_east(i,j) < 0.0d0) then
                uneg_mask(i,j) = 1
                D_plume_east_up(i,j) = D_plume(i+1,j)
             endif

          elseif (divu_mask_east(i,j) == 1) then  ! plume exists on only one side of the edge

             if (plume_mask_cell(i,j) == 1 .and. u_plume_east(i,j) > 0.0d0) then
                upos_mask(i,j) = 1
                D_plume_east_up(i,j) = D_plume(i,j)
             elseif (plume_mask_cell(i+1,j) == 1 .and. u_plume_east(i,j) < 0.0d0) then
                uneg_mask(i,j) = 1
                D_plume_east_up(i,j) = D_plume(i+1,j)
             endif

          endif

          ! mark the v component as positive or negative

          if (edge_mask_north(i,j) == 1) then  ! plume exists in both neighbor cells

!!             D_plume_north(i,j) = min(D_plume(i,j), D_plume(i,j+1))

             if (v_plume_north(i,j) > 0.0d0) then
                vpos_mask(i,j) = 1
                D_plume_north_up(i,j) = D_plume(i,j)
             elseif (v_plume_north(i,j) < 0.0d0) then
                vneg_mask(i,j) = 1
                D_plume_north_up(i,j) = D_plume(i,j+1)
             endif

          elseif (divu_mask_north(i,j) == 1) then  ! plume exists on only one side of the edge

             if (plume_mask_cell(i,j) == 1 .and. v_plume_north(i,j) > 0.0d0) then
                vpos_mask(i,j) = 1
                D_plume_north_up(i,j) = D_plume(i,j)
             elseif (plume_mask_cell(i,j+1) == 1 .and. v_plume_north(i,j) < 0.0d0) then
                vneg_mask(i,j) = 1
                D_plume_north_up(i,j) = D_plume(i,j+1)
             endif

          endif

       enddo  ! i
    enddo  ! j

    !--------------------------------------------------------------------
    ! Solve the equation dD/dt = e - d - div(Du)
    ! This is done in delta form:
    ! delta_D = (D_old - D_cur) + dt*(e-d) - dt*div(Du)
    !    where D_old is the old value of D and D_cur is the current guess.
    !--------------------------------------------------------------------
 
    ! compute nonzero matrix elements
    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             ! right-hand side
             ! This term includes
             ! (1) D_old - D_plume
             ! (2) (e - d)*dt
             ! (3) -dt * div(D_plume*u_plume), a divergence term based on the current guesses for D and u
             !
             n = cellID(i,j)

             rhs(n) = D_plume_old(i,j) - D_plume(i,j)

             rhs(n) = rhs(n) + dt_plume * (entrainment(i,j) - detrainment(i,j)) &
                             - (dt_plume/dx) * (D_plume_east_up(i,j) * u_plume_east(i,j) * divu_mask_east(i,j) &
                                              - D_plume_east_up(i-1,j) * u_plume_east(i-1,j) * divu_mask_east(i-1,j)) &
                             - (dt_plume/dy) * (D_plume_north_up(i,j) * v_plume_north(i,j) * divu_mask_north(i,j) &
                                              - D_plume_north_up(i,j-1) * v_plume_north(i,j-1) * divu_mask_north(i,j-1))

             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, '(D - D_old)/dt =', (D_plume(i,j) - D_plume_old(i,j)) / dt_plume
                print*, '(e - d) =', entrainment(i,j) - detrainment(i,j)

                dDu_dx = (1.0d0/dx) * (D_plume_east_up(i,j)*u_plume_east(i,j)*divu_mask_east(i,j) &
                                     - D_plume_east_up(i-1,j)*u_plume_east(i-1,j)*divu_mask_east(i-1,j))
                dDv_dy = (1.0d0/dy) * (D_plume_north_up(i,j)*v_plume_north(i,j)*divu_mask_north(i,j) &
                                     - D_plume_north_up(i,j-1)*v_plume_north(i,j-1)*divu_mask_north(i,j-1))
                print*, '-div(Du) =', -dDu_dx - dDv_dy
                print*, 'u_plume_west/east =', u_plume_east(i-1,j), u_plume_east(i,j)
                print*, 'v_plume_south/north =', v_plume_north(i,j-1), v_plume_north(i,j)
                print*, 'D_plume_west/east_up =', D_plume_east_up(i-1,j), D_plume_east_up(i,j)
                print*, 'D_plume_south/north_up =', D_plume_north_up(i,j-1), D_plume_north_up(i,j)
                print*, 'dDu_dx, dDv_dy =', dDu_dx, dDv_dy
                print*, 'residual (m/s) =', rhs(n)/dt_plume
             endif

             ! initialize the matrix diagonal
             A_plume(0,0,i,j) = 1.0d0

             ! Add matrix terms associated with u_plume*delta_D
             ! Note: The upos, uneg, vpos and vneg masks are constructed such that 
             !       if a mask = 1 at an edge, then the plume exists in the upstream cell.

             ! diagonal element
             A_plume(0,0,i,j) = A_plume(0,0,i,j)  &
                              + (dt_plume/dx) * u_plume_east(i,j) * upos_mask(i,j)  &
                              - (dt_plume/dx) * u_plume_east(i-1,j) * uneg_mask(i-1,j)  &
                              + (dt_plume/dy) * v_plume_north(i,j) * vpos_mask(i,j)  &
                              - (dt_plume/dy) * v_plume_north(i,j-1) * vneg_mask(i,j-1)

!             if (D_plume(i,j) == D_plume_east(i,j)) then
!                A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dx) * u_plume_east(i,j)
!             endif

!             if (D_plume(i,j) == D_plume_east(i-1,j)) then
!                A_plume(0,0,i,j) = A_plume(0,0,i,j) - (dt_plume/dx) * u_plume_east(i-1,j)
!             endif

!             if (D_plume(i,j) == D_plume_north(i,j)) then
!                A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dy) * v_plume_north(i,j)
!             endif

!             if (D_plume(i,j) == D_plume_north(i,j-1)) then
!                A_plume(0,0,i,j) = A_plume(0,0,i,j) - (dt_plume/dy) * v_plume_north(i,j-1)
!             endif

             ! off-diagonal elements
             A_plume(1,0,i,j)  = A_plume(1,0,i,j)  &
                               + (dt_plume/dx) * u_plume_east(i,j) * uneg_mask(i,j)
             A_plume(-1,0,i,j) = A_plume(-1,0,i,j)  &
                               - (dt_plume/dx) * u_plume_east(i-1,j) * upos_mask(i-1,j)
             A_plume(0,1,i,j)  = A_plume(0,1,i,j)  &
                               + (dt_plume/dy) * v_plume_north(i,j) * vneg_mask(i,j)
             A_plume(0,-1,i,j) = A_plume(0,-1,i,j)  &
                               - (dt_plume/dy) * v_plume_north(i,j-1) * vpos_mask(i,j-1)

             ! Note: If plume_mask_cell in the neighboring cell = 0, then D_plume in that cell is fixed at 0.
             !       In that case, there is no matrix element associated with the cell.

             ! TODO - May not need divu_mask_east if u_plume_east = 0 wherever divu_mask_east = 0
!             if (D_plume(i+1,j) == D_plume_east(i,j) .and. plume_mask_cell(i+1,j) == 1) then
!                A_plume(1,0,i,j)  = A_plume(1,0,i,j) + (dt_plume/dx) * u_plume_east(i,j)
!             endif

!             if (D_plume(i-1,j) == D_plume_east(i-1,j) .and. plume_mask_cell(i-1,j) == 1) then
!                A_plume(-1,0,i,j) = A_plume(-1,0,i,j) - (dt_plume/dx) * u_plume_east(i-1,j)
!             endif

!             if (D_plume(i,j+1) == D_plume_north(i,j) .and. plume_mask_cell(i,j+1) == 1) then
!                A_plume(0,1,i,j)  = A_plume(0,1,i,j) + (dt_plume/dy) * v_plume_north(i,j)
!             endif

!             if (D_plume(i,j-1) == D_plume_north(i,j-1) .and. plume_mask_cell(i,j-1) == 1) then
!                A_plume(0,-1,i,j) = A_plume(0,-1,i,j) - (dt_plume/dy) * v_plume_north(i,j-1)
!             endif

             if (i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, ' '
                print*, 'i, j, A_plume, rhs:', &
                     i, j, A_plume(0,0,i,j), A_plume(1,0,i,j), A_plume(-1,0,i,j), A_plume(0,1,i,j), A_plume(0,-1,i,j), &
                     rhs(cellID(i,j))
             endif

          endif  ! plume_mask_cell
          
       enddo  ! i
    enddo  ! j

    !WHL - Put a halo update here when running in parallel

!    print*, 'min, max A:', minval(A_plume), maxval(A_plume)
!    print*, 'min, max rhs:', minval(rhs), maxval(rhs)
!    print*, 'SLAP format'

    ! place nonzero elements in SLAP matrix format
    count = 0

    do n = 1, matrix_order

       i = iCellIndex(n)
       j = jCellIndex(n)

       if (plume_mask_cell(i,j) == 1) then

          ! loop over neighbor cells that can contribute terms to this matrix row

          do jA = -1,1
             do iA = -1,1

                if (A_plume(iA,jA,i,j) /= 0.0d0) then
                   count = count + 1
                   matrix%row(count) = n
                   matrix%col(count) = cellID(i+iA,j+jA)
                   matrix%val(count) = A_plume(iA,jA,i,j)
                      
                   if (matrix%col(count) == 0) then
                      print*, 'Bad matrix column: i, j, iA, jA =', i, j, iA, jA
                      stop
                   endif
                   
                   if (j==jtest) then
!!                         print*, 'i, j, iA, jA, row, col, val, rhs:', &
!!                              i, j, iA, jA, matrix%row(count), matrix%col(count), matrix%val(count), rhs(cellID(i,j))
                   endif
                      
                endif

             enddo   ! iA
          enddo  ! jA
          
       endif  ! plume_mask_cell

    enddo  ! n

    ! Set other matrix parameters
    matrix%order = matrix_order
    matrix%nonzeros = count
    matrix%symmetric = .false.

    ! call the SLAP solver

    call sparse_easy_solve(matrix,  rhs,      answer,   &
                           err,     niters,   whichsparse)
       
    print*, 'Called sparse_easy_solve: niters, err =', niters, err

    ! Update D_plume, given the answer vector from the solver.

    do n = 1, matrix_order
       i = iCellIndex(n)
       j = jCellIndex(n)

       if (i==itest .and. j==jtest) then
          print*, ' '
          print*, 'After solve, i, j, D_plume_old, delta_D, new D_plume:', &
                   i, j, D_plume_old(i,j), answer(n), D_plume(i,j) + answer(n)
       endif

       D_plume(i,j) = D_plume(i,j) + answer(n)

    enddo

    !TODO - Is this necessary?
    ! Check for D_plume > H_cavity.

    do j = 1, ny
       do i = 1, nx
          if (D_plume(i,j) > H_cavity(i,j)) then
!!!             print*, 'D_plume > H_cavity: i, j, D_plume, H_cavity =', i, j, D_plume(i,j), H_cavity(i,j)
          endif
       enddo
    enddo


  end subroutine compute_plume_thickness

!****************************************************

  subroutine compute_plume_thickness_with_eta(&
       nx,       ny,     &
       dx,       dy,     &
       dt_plume,         &
       free_surface,     &
       nonlinear_method, &
       plume_mask_cell,  &
       u_plume_east,     &
       v_plume_east,     &
       u_plume_north,    &
       v_plume_north,    &
       du_deta_west,     &
       du_deta_east,     &
       du_deta_northwest,&
       du_deta_northeast,&
       du_deta_southwest,&
       du_deta_southeast,&
       dv_deta_south,    &
       dv_deta_north,    &
       dv_deta_northwest,&
       dv_deta_northeast,&
       dv_deta_southwest,&
       dv_deta_southeast,&
       plume_speed_east, &
       plume_speed_north,&
       edge_mask_east,   &  !WHL - Do we need both pairs of masks?
       edge_mask_north,  &
       divu_mask_east,   &
       divu_mask_north,  &
       H_cavity,         &
       entrainment,      &
       detrainment,      &
       itest,    jtest,  &
       rtest,            &
       iter_melt,        &  !WHL - debug
       D_plume_old,      &
       D_plume,          &
       eta_plume)

    !--------------------------------------------------------------------
    ! Solve the continuity equation for D_plume:
    ! 
    !    dD/dt = e - d - del*(Du)
    !
    ! where e is entrainment and d is detrainment.
    !  
    ! The equation is solved in delta form:
    !
    !    delta_D = (D_old - D_cur) + dt*(e-d) - dt*div(Du)
    !
    ! where D_cur is the current guess for D, passed into the subroutine as D_plume.
    ! 
    ! TODO - Update the method description.
    ! The current velocity u_cur = (u_plume,v_plume) is also passed into the subroutine.
    ! Given u_cur and D_cur, the divergence term is expanded to first order as
    ! 
    !    div(Du) = div(D_cur*u_cur + u_cur*delta_D + D_cur*delta_u)
    !
    ! where 
    !    delta_u = (du/d_eta) * delta_eta
    !    delta_eta = delta_D where eta > 0, and delta_eta = 0 otherwise.
    !
    ! The terms containing delta_D are moved to the LHS and inserted in a matrix,
    ! giving a problem of the form
    !
    !    A*delta_D = rhs
    !
    ! This subroutine is called repeatedly until the residual is sufficiently small.
    !
    ! For now, I am using SLAP to solve the  matrix.
    ! Later, I plan to use a homegrown parallel solver.
    !--------------------------------------------------------------------

    ! for sparse_easy_solve
    use glimmer_sparse_type
    use glimmer_sparse

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    real(dp), intent(in) ::  &
         dt_plume               ! time step for plume solver (s)

    logical, intent(in) ::  &
         free_surface           ! true if computing PG force due to slope in free surface

    character(len=6), intent(in) ::  &
         nonlinear_method       ! method for solving nonlinear equations, 'Picard' or 'Newton' 

    !TODO - Only one pair of edge masks?
    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask_cell,     & ! = 1 for cells where scalar plume variables are computed, else = 0
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed, else = 0
         edge_mask_north,     & ! = 1 on north edges where plume velocity is computed, else = 0
         divu_mask_east,      & ! = 1 on east edges where divergence terms are computed, else = 0
         divu_mask_north        ! = 1 on north edges where divergence terms are  computed, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         u_plume_east,        & ! u_plume on east edges (m/s)
         v_plume_east,        & ! v_plume on east edges (m/s)
         u_plume_north,       & ! u_plume on north edges (m/s)
         v_plume_north,       & ! v_plume on north edges (m/s)
         plume_speed_east,    & ! plume speed on east edges (m/s)
         plume_speed_north,   & ! plume speed on north edges (m/s)
         H_cavity,            & ! ocean cavity thickness (m), lsrf - topg
         entrainment,         & ! entrainment at cell centers (m/s)
         detrainment            ! detrainment at cell centers (m/s)

   real(dp), dimension(nx,ny), intent(in) :: &
         du_deta_west,         & ! dependence of u_east on eta(i,j), west of the velocity point
         du_deta_east,         & ! dependence of u_east on eta(i+1,j), east of the velocity point
         du_deta_northwest,    & ! dependence of u_east on eta(i,j+1), northwest of the velocity point
         du_deta_northeast,    & ! dependence of u_east on eta(i+1,j+1), northeast of the velocity point
         du_deta_southwest,    & ! dependence of u_east on eta(i,j-1), southwest of the velocity point
         du_deta_southeast,    & ! dependence of u_east on eta(i+1,j-1), southeast of the velocity point
         dv_deta_south,        & ! dependence of v_north on eta(i,j), south of the velocity point
         dv_deta_north,        & ! dependence of v_north on eta(i,j+1), north of the velocity point
         dv_deta_northwest,    & ! dependence of v_north on eta(i-1,j+1), northwest of the velocity point
         dv_deta_northeast,    & ! dependence of v_north on eta(i+1,j+1), northeast of the velocity point
         dv_deta_southwest,    & ! dependence of v_north on eta(i-1,j), southwest of the velocity point
         dv_deta_southeast       ! dependence of v_north on eta(i+1,j), southeast of the velocity point

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    integer, intent(in) :: iter_melt  !WHL - debug

    real(dp), dimension(nx,ny), intent(in) ::  &
         D_plume_old            ! old plume thickness (m) at the previous time step

    real(dp), dimension(nx,ny), intent(inout) ::  &
         D_plume                ! on input, the current guess for the plume thickness (m)
                                ! on output, the new guess for the plume thickness
                                ! Note: D_plume is capped at H_cavity

    real(dp), dimension(nx,ny), intent(out) ::  &
         eta_plume              ! displacement of plume surface (m), once D_plume = H_cavity

    ! local variables

    !WHL - Remove _up variables?
    real(dp), dimension(nx,ny) ::  &
         D_plume_latest,      & ! D_plume from the most recent iteration (m)
         D_plume_cap,         & ! min(D_plume, H_cavity)
         D_plume_east,        & ! plume thickness at each east edge (m)
         D_plume_north,       & ! plume thickness at each north edge (m)
         D_plume_east_up,     & ! upstream plume thickness at each east edge (m)
         D_plume_north_up       ! upstream plume thickness at each north edge (m)

    real(dp), dimension(nx,ny) ::  &
         ux,                  & ! pre-multiplier for deta/dx on east edges
         uy,                  & ! pre-multiplier for deta/dy on east edges
         vx,                  & ! pre-multiplier for deta/dx on north edges
         vy                     ! pre-multiplier for deta/dy on north edges

    integer, dimension(nx,ny) ::  &
         upos_mask,   & ! = 1 at edges where u_plume > 0, else = 0
         uneg_mask,   & ! = 1 at edges where u_plume < 0, else = 0
         vpos_mask,   & ! = 1 at edges where v_plume > 0, else = 0
         vneg_mask      ! = 1 at edges where v_plume < 0, else = 0

    integer, dimension(nx,ny) ::  &
         eta_mask       ! = 1 where D_plume = H_cavity and eta_plume can be > 0, else = 0

    real(dp), dimension(-1:1,-1:1,nx,ny) ::  &
         A_plume        ! array holding nonzero matrix elements on the structured mesh
                        ! up to 9 nonzero elements per row of the matrix

    type(sparse_matrix_type) ::  &
         matrix         ! sparse matrix for SLAP solver, defined in glimmer_sparse_types
                        ! includes nonzeroes, order, col, row, val

    real(dp), dimension(:), allocatable ::  &
         rhs,         & ! right-hand-side vector, passed to solver
         answer         ! answer vector, returned from solver

    real(dp) ::  &
         err            ! solution error the solver

    integer, dimension(nx,ny) ::  &
         cellID         ! integer ID for each cell

    integer, dimension(nx*ny) ::  &
         iCellIndex, jCellIndex   ! indices for mapping cellID back to i and j

    real(dp) :: &
         D_east, D_north,      & ! current estimate of D_plume, averaged to east and north edges
         denom                   ! denominator in the expression for velocity

    integer :: niters   ! iteration counter
    integer :: i, j     ! horizontal indices
    integer :: iA, jA   ! horizontal index shifts, in range -1:1
    integer :: n        ! matrix index
    integer :: count    ! counter
    integer :: matrix_order  ! size of square matrix
    integer :: max_nonzeros  ! max number of nonzero elements in matrix

    ! SLAP linear solver (BICG or GMRES)
    ! For ISOMIP+, BICG is a bit faster, requiring 3 or 4 linear iterations compared to 6 or 7 for GMRES.
    !TODO - Replace with homegrown solver

    integer, parameter :: &
         whichsparse = HO_SPARSE_BICG
!!         whichsparse = HO_SPARSE_GMRES

    real(dp), parameter :: &
         relax_D = 0.5d0,       &  !WHL - Remove relax_D?
!!         relax_eta = 1.0d0  ! relaxation parameter to prevent D from oscillating about H_cavity
         relax_eta = 0.50d0  ! relaxation parameter to prevent D from oscillating about H_cavity
                             ! 0 < relax_eta <=1; a smaller value gives stronger overrelaxation.
                             ! In practice, a value of ~0.5 seems to work well.
                             ! With smaller values, convergence is slower but may be more robust.

    !WHL - debug
    real(dp) :: dDu_dx, dDv_dy
    real(dp) :: diag_east, diag_west, diag_north, diag_south
    real(dp) :: offdiag_east, offdiag_west, offdiag_north, offdiag_south
    real(dp), dimension(nx,ny) :: eta_plume_latest

    logical, parameter :: apply_jacobian = .true.
!!    logical, parameter :: apply_jacobian = .false.

    print*, ' '
    print*, 'In plume_thickness_solver: itest, jtest =', itest, jtest

    ! count plume cells in matrix solve
    ! loop over locally owned cells
    count = 0
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j)==1) then
             count = count + 1
             cellID(i,j) = count
             iCellIndex(count) = i
             jCellIndex(count) = j
          endif
       enddo
    enddo

    ! initialize and allocate
    matrix_order = count
    max_nonzeros = matrix_order * 9

    allocate(matrix%row(max_nonzeros), matrix%col(max_nonzeros), matrix%val(max_nonzeros))
    allocate(rhs(matrix_order), answer(matrix_order))

    A_plume(:,:,:,:) = 0.0d0
    rhs(:) = 0.0d0
    answer(:) = 0.0d0

    ! Save the latest iterates for D_plume and eta_plume.
    ! Note: These are different from D_plume_old and eta_plume_old, the values at the start of the time step.
    D_plume_latest(:,:) = D_plume(:,:)
    eta_plume_latest(:,:) = eta_plume(:,:)

    ! Create a mask: eta_mask = 1 where D_plume = H_cavity, else = 0

    where (plume_mask_cell == 1 .and. D_plume >= H_cavity)
       eta_mask = 1
    elsewhere
       eta_mask = 0
    endwhere

    !WHL - debug
    print*, ' '
    print*, 'New eta_mask, rank =', rtest
    do j = jtest+3, jtest-3, -1
       do i = itest-3, itest+3
          write(6,'(i12)',advance='no') eta_mask(i,j)
       enddo
       write(6,*) ' '
    enddo

    !WHL - debug
    i = itest
    j = jtest
    print*, 'i, j, latest D, eta =', i, j, D_plume(i,j), eta_plume(i,j)

    ! Compute D_plume edge values used in the finite difference expression for divergence.
    ! No solution is ideal here:
    ! With a centered difference for D_plume on edges, thin plumes upstream from thick plumes can overempty.
    ! With an upstream difference, there is no overemptying, but there can be oscillations as a result
    !  of u or v changing sign between one iteration and the next.
    ! Choosing the min value is not very accurate (only first order, and too restrictive of outflow when
    !  a thick plume is upstream of a thin plume), but it inhibits overemptying and oscillations.

    ! Store for each edge the plume thickness in the upstream cell.

    ! Note: The use of edge_mask and divu_mask is a bit subtle.
    !       If edge_mask = 1, the plume exists on both sides of the edge, and one cell is clearly upstream.
    !       If edge_mask = 0 but divu_mask = 1, the plume exists on only one side of the edge.
    !        We identify an upstream plume thickness only if the plume exists in the upstream cell
    !        and the flow is out of the plume domain.
    !        (This generally is the case, but we check here to be sure.)

    upos_mask(:,:) = 0
    uneg_mask(:,:) = 0
    vpos_mask(:,:) = 0
    vneg_mask(:,:) = 0

    D_plume_east(:,:) = 0.0d0
    D_plume_north(:,:) = 0.0d0

    !WHL - Remove _up variables?
    D_plume_east_up(:,:) = 0.0d0
    D_plume_north_up(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx

          ! mark the u component as positive or negative
          !WHL - Remove these calculations?
          if (edge_mask_east(i,j) == 1) then ! plume exists in both neighbor cells

             D_plume_east(i,j) = min(D_plume(i,j), D_plume(i+1,j))

!             if (u_plume_east(i,j) > 0.0d0) then
!                upos_mask(i,j) = 1
!                D_plume_east_up(i,j) = D_plume_cap(i,j)
!             elseif (u_plume_east(i,j) < 0.0d0) then
!                uneg_mask(i,j) = 1
!                D_plume_east_up(i,j) = D_plume_cap(i+1,j)
!             endif

          elseif (divu_mask_east(i,j) == 1) then  ! plume exists on only one side of the edge

             if (plume_mask_cell(i,j) == 1) then
                D_plume_east(i,j) = D_plume(i,j)
             else
                D_plume_east(i,j) = D_plume(i+1,j)
             endif

!             if (plume_mask_cell(i,j) == 1 .and. u_plume_east(i,j) > 0.0d0) then
!                upos_mask(i,j) = 1
!                D_plume_east_up(i,j) = D_plume_cap(i,j)
!             elseif (plume_mask_cell(i+1,j) == 1 .and. u_plume_east(i,j) < 0.0d0) then
!                uneg_mask(i,j) = 1
!                D_plume_east_up(i,j) = D_plume_cap(i+1,j)
!             endif

          endif

          ! mark the v component as positive or negative
          if (edge_mask_north(i,j) == 1) then

             D_plume_north(i,j) = min(D_plume(i,j), D_plume(i,j+1))

!             if (v_plume_north(i,j) > 0.0d0) then
!                vpos_mask(i,j) = 1
!                D_plume_north_up(i,j) = D_plume_cap(i,j)
!             elseif (v_plume_north(i,j) < 0.0d0) then
!                vneg_mask(i,j) = 1
!                D_plume_north_up(i,j) = D_plume_cap(i,j+1)
!             endif

          elseif (divu_mask_north(i,j) == 1) then  ! plume exists on only one side of the edge

             if (plume_mask_cell(i,j) == 1) then
                D_plume_north(i,j) = D_plume(i,j)
             else
                D_plume_north(i,j) = D_plume(i,j+1)
             endif

!             if (plume_mask_cell(i,j) == 1 .and. v_plume_north(i,j) > 0.0d0) then
!                vpos_mask(i,j) = 1
!                D_plume_north_up(i,j) = D_plume_cap(i,j)
!             elseif (plume_mask_cell(i,j+1) == 1 .and. v_plume_north(i,j) < 0.0d0) then
!                vneg_mask(i,j) = 1
!                D_plume_north_up(i,j) = D_plume_cap(i,j+1)
!             endif

          endif

       enddo  ! i
    enddo  ! j

    ! Compute terms that are multiplied by eta on the LHS.
    ! At east edges we expand u = u0 + ux * deta/dx + uy * deta/dy.
    ! At north edges we expand v = v0 + vx * delta_eta + vy * deta/dy.

    ux(:,:) = 0.0d0
    uy(:,:) = 0.0d0
    vx(:,:) = 0.0d0
    vy(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx
          
          if (divu_mask_east(i,j) == 1) then
             D_east = 0.5d0 * (D_plume(i,j) + D_plume(i+1,j))
             denom = (c_drag*plume_speed_east(i,j))**2 + (f_coriolis*D_east)**2
!!             ux(i,j) = -c_drag * plume_speed_east(i,j) * grav * D_east / (dx * denom)
             ux(i,j) = -c_drag * plume_speed_east(i,j) * grav * D_east / denom
!!             f_east(i,j) = f_coriolis * grav * D_east**2 / (4.0d0 * dy * denom)
             uy(i,j) = -f_coriolis * grav * D_east**2 / denom

          endif

          if (divu_mask_north(i,j) == 1) then
             D_north = 0.5d0 * (D_plume(i,j) + D_plume(i,j+1))
             denom = (c_drag*plume_speed_north(i,j))**2 + (f_coriolis*D_north)**2
!!                c_north(i,j) = c_drag * plume_speed_north(i,j) * grav * D_north / (dy * denom)
             vy(i,j) = -c_drag * plume_speed_north(i,j) * grav * D_north / denom
!!                f_north(i,j) = f_coriolis * grav * D_north**2 / (4.0d0 * dx * denom)
             vx(i,j) = f_coriolis * grav * D_north**2 / denom
          endif
             
          if (i==itest .and. j==jtest) then
             print*, 'i, j, denom:', i, j, denom
             print*, 'ux_west/dx, ux_east/dx:',   ux(i-1,j)/dx, ux(i,j)/dx
             print*, '      plume_speed_west:',   plume_speed_east(i-1,j)
             print*, '      plume_speed_east:',   plume_speed_east(i,j)
             print*, 'uy_west/dx, uy_east/dx:',   uy(i-1,j)/dx, uy(i,j)/dx
             print*, 'vx_south/dy, vx_north/dt:', vx(i,j-1)/dy, vx(i,j)/dy
             print*, 'vy_south/dy, vy_north/dy:', vy(i,j-1)/dy, vy(i,j)/dy
          endif
             
       enddo
    enddo

    !--------------------------------------------------------------------
    ! Solve the equation dD/dt = e - d - div(Du).
    ! This is done in delta form:
    ! delta_D = (D_old - D_cur) + dt*(e-d) - dt*div(Du)
    !    where D_cur is the current guess for D.
    ! Terms involving eta are placed on the LHS and solved implicitly.
    !--------------------------------------------------------------------
 
    ! compute nonzero matrix elements
    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             ! right-hand side
             ! This term includes
             ! (1) D_old - D_cur, the difference between D_plume_old and the current guess
             ! (2) (e - d)*dt
             ! (3) -dt * div(D_cur*u_cur), a divergence term based on the current guesses for D and u
             !
             ! Note: Term (1) appears when we are solving for D, but not when solving for eta.
             !       When solving for eta, we seek a balance between entrainment/detrainment and divergence, with no change in D.

             !TODO - May not need divu_mask_east, if u_plume_east = 0 at edges with divu_mask_east = 0
             n = cellID(i,j)

             if (D_plume(i,j) < H_cavity(i,j)) then
                rhs(n) = D_plume_old(i,j) - D_plume(i,j)
                if (i==itest .and. j==jtest) print*, 'put delta_D on rhs, i, j, =', itest, jtest
             endif

             rhs(n) = rhs(n) + dt_plume * (entrainment(i,j) - detrainment(i,j)) &
                             - (dt_plume/dx) * (D_plume_east(i,j) * u_plume_east(i,j) * divu_mask_east(i,j) &
                                              - D_plume_east(i-1,j) * u_plume_east(i-1,j) * divu_mask_east(i-1,j)) &
                             - (dt_plume/dy) * (D_plume_north(i,j) * v_plume_north(i,j) * divu_mask_north(i,j) &
                                              - D_plume_north(i,j-1) * v_plume_north(i,j-1) * divu_mask_north(i,j-1))

             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, '(D - D_old)/dt =', (D_plume(i,j) - D_plume_old(i,j)) / dt_plume
                print*, '(e - d) =', entrainment(i,j) - detrainment(i,j)

                dDu_dx = (1.0d0/dx) * (D_plume_east(i,j)*u_plume_east(i,j)*divu_mask_east(i,j) &
                                     - D_plume_east(i-1,j)*u_plume_east(i-1,j)*divu_mask_east(i-1,j))
                dDv_dy = (1.0d0/dy) * (D_plume_north(i,j)*v_plume_north(i,j)*divu_mask_north(i,j) &
                                     - D_plume_north(i,j-1)*v_plume_north(i,j-1)*divu_mask_north(i,j-1))
                print*, '-div(Du) =', -dDu_dx - dDv_dy
                print*, 'u_plume_west/east =', u_plume_east(i-1,j), u_plume_east(i,j)
                print*, 'v_plume_south/north =', v_plume_north(i,j-1), v_plume_north(i,j)
                print*, 'D_plume_west/east =', D_plume_east(i-1,j), D_plume_east(i,j)
                print*, 'D_plume_south/north =', D_plume_north(i,j-1), D_plume_north(i,j)
                print*, 'dDu_dx, dDv_dy =', dDu_dx, dDv_dy
                print*, 'residual (m/s) =', rhs(n)/dt_plume
             endif

             ! initialize the matrix diagonal
!!             A_plume(0,0,i,j) = 1.0d0

             ! Add matrix terms associated with delta_eta
             ! Note: We solve for delta_D in cells with D < H_cavity.  Elsewhere, we solve for delta_eta.

             if (D_plume(i,j) < H_cavity(i,j)) then
                A_plume(0,0,i,j) = 1.0d0
             endif

             !WHL - With apply_jacobian = T, this would become (check dx term in denom):
!!             A_plume( 1,0,i,j) = A_plume(1,0,i,j) + eta_mask(i+1,j) * (dt_plume/dx) * D_plume_east(i,j) * du_deta_east(i,j)
!!             A_plume( 0,0,i,j) = A_plume(0,0,i,j) - eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i,j) * du_deta_west(i,j)


             ! Add terms associated with deta/dx on east edge
             A_plume( 1,0,i,j) = A_plume(1,0,i,j) + &
                  eta_mask(i+1,j) * (dt_plume/dx) * D_plume_east(i,j) * ux(i,j) / dx
             A_plume( 0,0,i,j) = A_plume(0,0,i,j) - &
                  eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i,j) * ux(i,j) / dx

             ! Add terms associated with deta/dx on west edge
             A_plume( 0,0,i,j) = A_plume( 0,0,i,j) - &
                  eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i-1,j) * ux(i-1,j) / dx
             A_plume(-1,0,i,j) = A_plume(-1,0,i,j) + &
                  eta_mask(i-1,j) * (dt_plume/dx) * D_plume_east(i-1,j) * ux(i-1,j) / dx

             !WHL - debug
             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, 'u_east, deta/dx terms:'
                print*, 'A(-1,0)', eta_mask(i-1,j) * (dt_plume/dx) * D_plume_east(i-1,j) * ux(i-1,j) / dx
                print*, 'A(0,0)', -eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i-1,j) * ux(i-1,j) / dx
                print*, 'A(0,0)', -eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i,j) * ux(i,j) / dx
                print*, 'A(1,0)', eta_mask(i+1,j) * (dt_plume/dx) * D_plume_east(i,j) * ux(i,j) / dx
             endif
             
             ! Add terms associated with deta/dy on east edge
             if (divu_mask_north(i,j) == 1) then
                A_plume(0,1,i,j) = A_plume(0,1,i,j) + &
                     eta_mask(i,j+1) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)
                A_plume(0,0,i,j) = A_plume(0,0,i,j) - &
                     eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)

                !WHL - debug
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Increment north(i,j)'
                   print*, 'dA(0,1)',  eta_mask(i,j+1) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)
                   print*, 'dA(0,0)', -eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy) 
                endif

             endif

             if (divu_mask_north(i+1,j) == 1) then
                A_plume(1,1,i,j) = A_plume(1,1,i,j) + &
                     eta_mask(i+1,j+1) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)
                A_plume(1,0,i,j) = A_plume(1,0,i,j) - &
                     eta_mask(i+1,j)   * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)

                !WHL - debug
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Increment north(i+1,j)'
                   print*, 'dA(1,1)',  eta_mask(i+1,j+1) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy) 
                   print*, 'dA(1,0)', -eta_mask(i+1,j) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy) 
                endif

             endif

             if (divu_mask_north(i,j-1) == 1) then
                A_plume(0, 0,i,j) = A_plume(0, 0,i,j) + &
                     eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)
                A_plume(0,-1,i,j) = A_plume(0,-1,i,j) - &
                     eta_mask(i,j-1) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)

                !WHL - debug
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Increment north(i,j-1)'
                   print*, 'dA(0,0)',   eta_mask(i,j) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy) 
                   print*, 'dA(0,-1)', -eta_mask(i,j-1) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy) 
                endif

             endif

             if (divu_mask_north(i+1,j-1) == 1) then
                A_plume(1, 0,i,j) = A_plume(1, 0,i,j) + &
                     eta_mask(i+1,j)   * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)
                A_plume(1,-1,i,j) = A_plume(1,-1,i,j) - &
                     eta_mask(i+1,j-1) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy)

                !WHL - debug
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Increment north(i+1,j-1)'
                   print*, 'dA(1,0)',   eta_mask(i+1,j) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy) 
                   print*, 'dA(1,-1)', -eta_mask(i+1,j-1) * (dt_plume/dx) * D_plume_east(i,j) * uy(i,j) / (4.0d0 * dy) 
                endif

             endif

             ! Add terms associated with deta/dy on west edge
             if (divu_mask_north(i-1,j) == 1) then
                A_plume(-1,1,i,j) = A_plume(-1,1,i,j) - &
                     eta_mask(i-1,j+1) * (dt_plume/dx) * D_plume_east(i-1,j) * uy(i-1,j) / (4.0d0 * dy)
                A_plume(-1,0,i,j) = A_plume(-1,0,i,j) + &
                     eta_mask(i-1,j)   * (dt_plume/dx) * D_plume_east(i-1,j) * uy(i-1,j) / (4.0d0 * dy)
             endif

             if (divu_mask_north(i,j) == 1) then
                A_plume(0,1,i,j) = A_plume(0,1,i,j) - &
                     eta_mask(i,j+1) * (dt_plume/dx) * D_plume_east(i-1,j) * uy(i-1,j) / (4.0d0 * dy)
                A_plume(0,0,i,j) = A_plume(0,0,i,j) + &
                     eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i-1,j) * uy(i-1,j) / (4.0d0 * dy)
             endif

             if (divu_mask_north(i-1,j-1) == 1) then
                A_plume(-1, 0,i,j) = A_plume(-1, 0,i,j) - &
                     eta_mask(i-1,j)   * (dt_plume/dx) * D_plume_east(i-1,j) * uy(i-1,j) / (4.0d0 * dy)
                A_plume(-1,-1,i,j) = A_plume(-1,-1,i,j) + &
                     eta_mask(i-1,j-1) * (dt_plume/dx) * D_plume_east(i-1,j) * uy(i-1,j) / (4.0d0 * dy)
             endif

             if (divu_mask_north(i,j-1) == 1) then
                A_plume(0, 0,i,j) = A_plume(0, 0,i,j) - &
                     eta_mask(i,j)   * (dt_plume/dx) * D_plume_east(i-1,j) * uy(i-1,j) / (4.0d0 * dy)
                A_plume(0,-1,i,j) = A_plume(0,-1,i,j) + &
                     eta_mask(i,j-1) * (dt_plume/dx) * D_plume_east(i-1,j) * uy(i-1,j) / (4.0d0 * dy)
             endif

             ! Add terms associated with deta/dy on north edge
             A_plume(0, 1,i,j) = A_plume(0,1,i,j) + &
                  eta_mask(i,j+1)  * (dt_plume/dy) * D_plume_north(i,j) * vy(i,j) / dy
             A_plume(0, 0,i,j) = A_plume(0,0,i,j) - &
                  eta_mask(i,j)    * (dt_plume/dy) * D_plume_north(i,j) * vy(i,j) / dy

             ! Add terms associated with deta/dy on south edge
             A_plume(0, 0,i,j) = A_plume(0, 0,i,j) - &
                  eta_mask(i,j)   * (dt_plume/dy) * D_plume_north(i,j-1) * vy(i,j-1) / dy
             A_plume(0,-1,i,j) = A_plume(0,-1,i,j) + &
                  eta_mask(i,j-1) * (dt_plume/dy) * D_plume_north(i,j-1) * vy(i,j-1) / dy

             !WHL - debug
             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, 'v_south: deta/dy terms:'
                print*, 'A(0,1)',  eta_mask(i,j+1)  * (dt_plume/dy) * D_plume_north(i,j) * vy(i,j) / dy
                print*, 'A(0,0)', -eta_mask(i,j) * (dt_plume/dy) * D_plume_north(i,j) * vy(i,j) / dy
                print*, 'A(0,0)', -eta_mask(i,j) * (dt_plume/dy) * D_plume_north(i,j-1) * vy(i,j-1) / dy
                print*, 'A(0,-1)', eta_mask(i,j-1) * (dt_plume/dy) * D_plume_north(i,j-1) * vy(i,j-1) / dy
             endif
             

             ! Add terms associated with deta/dx on north edge
             if (divu_mask_east(i,j+1) == 1) then
                A_plume(1,1,i,j) = A_plume(1,1,i,j) + &
                     eta_mask(i+1,j+1) * (dt_plume/dy) * D_plume_north(i,j) * vx(i,j) / (4.0d0 * dy)
                A_plume(0,1,i,j) = A_plume(0,1,i,j) - &
                     eta_mask(i,j+1)   * (dt_plume/dy) * D_plume_north(i,j) * vx(i,j) / (4.0d0 * dy)
             endif

             if (divu_mask_east(i,j) == 1) then
                A_plume(1,0,i,j) = A_plume(1,0,i,j) + &
                     eta_mask(i+1,j) * (dt_plume/dy) * D_plume_north(i,j) * vx(i,j) / (4.0d0 * dy)
                A_plume(0,0,i,j) = A_plume(0,0,i,j) - &
                     eta_mask(i,j)   * (dt_plume/dy) * D_plume_north(i,j) * vx(i,j) / (4.0d0 * dy)
             endif

             if (divu_mask_east(i-1,j+1) == 1) then
                A_plume( 0,1,i,j) = A_plume( 0,1,i,j) + &
                     eta_mask(i,j+1)   * (dt_plume/dy) * D_plume_north(i,j) * vx(i,j) / (4.0d0 * dy)
                A_plume(-1,1,i,j) = A_plume(-1,1,i,j) - &
                     eta_mask(i-1,j+1) * (dt_plume/dy) * D_plume_north(i,j) * vx(i,j) / (4.0d0 * dy)
             endif

             if (divu_mask_east(i-1,j) == 1) then
                A_plume( 0,0,i,j) = A_plume( 0,0,i,j) + &
                     eta_mask(i,j)   * (dt_plume/dy) * D_plume_north(i,j) * vx(i,j) / (4.0d0 * dy)
                A_plume(-1,0,i,j) = A_plume(-1,0,i,j) - &
                     eta_mask(i-1,j) * (dt_plume/dy) * D_plume_north(i,j) * vx(i,j) / (4.0d0 * dy)
             endif

             ! Add terms associated with deta/dx on south edge
             if (divu_mask_east(i,j) == 1) then
                A_plume(1,0,i,j) = A_plume(1,0,i,j) - &
                     eta_mask(i+1,j) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)
                A_plume(0,0,i,j) = A_plume(0,0,i,j) + &
                     eta_mask(i,j)   * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)

                !WHL - debug
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Increment east(i,j)'
                   print*, 'dA(1,0)', -eta_mask(i+1,j) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy) 
                   print*, 'dA(0,0)',  eta_mask(i,j) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)
                endif

             endif

             if (divu_mask_east(i,j-1) == 1) then
                A_plume(1,-1,i,j) = A_plume(1,-1,i,j) - &
                     eta_mask(i+1,j-1) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)
                A_plume(0,-1,i,j) = A_plume(0,-1,i,j) + &
                     eta_mask(i,j-1)   * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)

                !WHL - debug
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Increment east(i,j-1)'
                   print*, 'dA(1,-1)', -eta_mask(i+1,j-1) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy) 
                   print*, 'dA(0,-1)',  eta_mask(i,j-1) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)
                endif

             endif

             if (divu_mask_east(i-1,j) == 1) then
                A_plume( 0,0,i,j) = A_plume( 0,0,i,j) - &
                     eta_mask(i,j)   * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)
                A_plume(-1,0,i,j) = A_plume(-1,0,i,j) + &
                     eta_mask(i-1,j) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)

                !WHL - debug
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Increment east(i-1,j)'
                   print*, 'dA(0,0)', -eta_mask(i,j) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy) 
                   print*, 'dA(-1,0)',  eta_mask(i-1,j) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)
                endif

             endif

             if (divu_mask_east(i-1,j-1) == 1) then
                A_plume( 0,-1,i,j) = A_plume( 0,-1,i,j) - &
                     eta_mask(i,j-1)   * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)
                A_plume(-1,-1,i,j) = A_plume(-1,-1,i,j) + &
                     eta_mask(i-1,j-1) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)

                !WHL - debug
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Increment east(i-1,j-1)'
                   print*, 'dA(0,-1)', -eta_mask(i,j-1) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy) 
                   print*, 'dA(-1,-1)',  eta_mask(i-1,j-1) * (dt_plume/dy) * D_plume_north(i,j-1) * vx(i,j-1) / (4.0d0 * dy)
                endif

             endif

             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, 'After adding eta terms:'
                print*, 'i, j, rhs:', i, j, rhs(cellID(i,j))
                print*, 'A(-1:1, 1):', A_plume(-1:1, 1,i,j)
                print*, 'A(-1:1, 0):', A_plume(-1:1, 0,i,j)
                print*, 'A(-1:1,-1):', A_plume(-1:1,-1,i,j)
                print*, ' '
             endif

          endif  ! plume_mask_cell
          
       enddo  ! i
    enddo  ! j

    !WHL - Put a halo update here when running in parallel

!    print*, 'min, max A:', minval(A_plume), maxval(A_plume)
!    print*, 'min, max rhs:', minval(rhs), maxval(rhs)
!    print*, 'SLAP format'

    ! place nonzero elements in SLAP matrix format
    count = 0

    do n = 1, matrix_order

       i = iCellIndex(n)
       j = jCellIndex(n)

       if (plume_mask_cell(i,j) == 1) then

          ! loop over neighbor cells that can contribute terms to this matrix row

          do jA = -1,1
             do iA = -1,1

                if (A_plume(iA,jA,i,j) /= 0.0d0) then
                   count = count + 1
                   matrix%row(count) = n
                   matrix%col(count) = cellID(i+iA,j+jA)
                   matrix%val(count) = A_plume(iA,jA,i,j)
                      
                   if (matrix%col(count) == 0) then
                      print*, 'Bad matrix column: i, j, iA, jA =', i, j, iA, jA
                      stop
                   endif
                   
                   if (j==jtest) then
!!                         print*, 'i, j, iA, jA, row, col, val, rhs:', &
!!                              i, j, iA, jA, matrix%row(count), matrix%col(count), matrix%val(count), rhs(cellID(i,j))
                   endif
                      
                endif

             enddo   ! iA
          enddo  ! jA
          
       endif  ! plume_mask_cell

    enddo  ! n

    ! Set other matrix parameters
    matrix%order = matrix_order
    matrix%nonzeros = count
    matrix%symmetric = .false.

    ! call the SLAP solver

    call sparse_easy_solve(matrix,  rhs,      answer,   &
                           err,     niters,   whichsparse)
       
    print*, 'Called sparse_easy_solve: niters, err =', niters, err

    ! Update D_plume and eta_plume, given the answer vector from the solver.
    ! Note: Where D_plume < H_cavity, the answer vector contains delta_D_plume.
    !       Where D_plume = H_cavity, the answer vector contains delta_eta_plume.

    do n = 1, matrix_order
       i = iCellIndex(n)
       j = jCellIndex(n)

       if (D_plume(i,j) < H_cavity(i,j)) then
          D_plume(i,j) = D_plume(i,j) + answer(n)
       else
          eta_plume(i,j) = eta_plume(i,j) + answer(n)
       endif

       if (i==itest .and. j==jtest) then
          print*, ' '
          print*, 'After solve, i, j =:', itest, jtest
          print*, 'H_cavity, D_plume_latest, D_plume, eta_latest, eta:', &
                   H_cavity(i,j), D_plume_latest(i,j), D_plume(i,j), eta_plume_latest(i,j), eta_plume(i,j)
       endif

    enddo

    !WHL - debug
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Predictions without Coriolis:'
                print*, ' '
                print*, 'old u_plume_east =', u_plume_east(i,j)
!!                print*, 'predicted change =', ux(i,j)/dx * ( (eta_plume(i+1,j) - eta_plume_latest(i+1,j))  &
!!                                                           - (eta_plume(i,j) - eta_plume_latest(i,j)) ) 
                print*, 'predicted change =', ux(i,j)/dx * ( (eta_plume(i+1,j) - eta_plume_latest(i+1,j))  &
                                                           - (eta_plume(i,j) - eta_plume_latest(i,j)) ) 
                print*, ' '
                print*, 'old u_plume_west =', u_plume_east(i-1,j)
                print*, 'predicted change =', ux(i-1,j)/dx * ( (eta_plume(i,j) - eta_plume_latest(i,j))  &
                                                             - (eta_plume(i-1,j) - eta_plume_latest(i-1,j)) )
                print*, ' '
                print*, 'old v_plume_north =', v_plume_north(i,j)
                print*, 'predicted change =', vy(i,j)/dy  * ( (eta_plume(i,j+1) - eta_plume_latest(i,j+1))  &
                                                            - (eta_plume(i,j) - eta_plume_latest(i,j)) )
                print*, ' '
                print*, 'old v_plume_south =', v_plume_north(i,j-1)
                print*, 'predicted change =', vy(i,j-1)/dy * ( (eta_plume(i,j) - eta_plume_latest(i,j))  &
                                                             - (eta_plume(i,j-1) - eta_plume_latest(i,j-1)) )
             endif

          endif
       enddo
    enddo

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             ! Where D_plume > H_cavity, set D_plume = H_cavity.
             ! We will solve for eta_plume in this cell on the next iteration.

             if (D_plume_latest(i,j) < H_cavity(i,j) .and. D_plume(i,j) >= H_cavity(i,j)) then
                !WHL - debug
!!                eta_plume(i,j) = relax_eta * (D_plume(i,j) - H_cavity(i,j))
                D_plume(i,j) = H_cavity(i,j)
                print*, 'INFLATE: i, j, new D_plume, eta_plume:', i, j, D_plume(i,j), eta_plume(i,j)
             endif

             if (eta_plume_latest(i,j) >= 0.0d0) then
                eta_plume(i,j) = eta_plume_latest(i,j) + relax_eta*(eta_plume(i,j) - eta_plume_latest(i,j))
             endif

             ! Where eta_plume_latest > 0 but eta_plume <= 0, set eta_plume = 0.
             ! We will solve for delta_D_plume in this cell on the next iteration.
             if (D_plume_latest(i,j) >= H_cavity(i,j) .and. eta_plume(i,j) <= 0.0d0) then
!!                D_plume(i,j) = H_cavity(i,j)*(1.d0 - 1.d-11)  ! set D_plume to slightly less than H_cavity
                D_plume(i,j) = H_cavity(i,j) + relax_eta*eta_plume(i,j)
                eta_plume(i,j) = 0.0d0
                print*, 'DEFLATE: i, j, D_plume, H_cavity - D_plume:', i, j, H_cavity(i,j), H_cavity(i,j) - D_plume(i,j)
             endif

             if (i==itest .and. j==jtest) then
                print*, 'After correction, i, j =:', itest, jtest
                print*, 'H_cavity, D_plume_latest, D_plume, eta_latest, eta:', &
                         H_cavity(i,j), D_plume_latest(i,j), D_plume(i,j), eta_plume_latest(i,j), eta_plume(i,j)
             endif

          endif  ! plume_mask_cell
       enddo  ! i
    enddo  ! j

    !TODO - Is this necessary?
    ! Check for D_plume > H_cavity.

    do j = 1, ny
       do i = 1, nx
          if (D_plume(i,j) > H_cavity(i,j)) then
             print*, 'Error, D_plume > H_cavity'
             print*, 'i, j, D_plume, H_cavity:', i, j, D_plume(i,j), H_cavity(i,j)
             stop
          endif
       enddo
    enddo

  end subroutine compute_plume_thickness_with_eta

!****************************************************

  subroutine compute_plume_thickness_old(&
       nx,       ny,     &
       dx,       dy,     &
       dt_plume,         &
       free_surface,     &
       nonlinear_method, &
       plume_mask_cell,  &
       u_plume_east,     &
       v_plume_east,     &
       u_plume_north,    &
       v_plume_north,    &
       du_deta_west,     &
       du_deta_east,     &
       du_deta_northwest,&
       du_deta_northeast,&
       du_deta_southwest,&
       du_deta_southeast,&
       dv_deta_south,    &
       dv_deta_north,    &
       dv_deta_northwest,&
       dv_deta_northeast,&
       dv_deta_southwest,&
       dv_deta_southeast,&
       plume_speed_east, &
       plume_speed_north,&
       edge_mask_east,   &  !WHL - Do we need both pairs of masks?
       edge_mask_north,  &
       divu_mask_east,   &
       divu_mask_north,  &
       H_cavity,         &
       entrainment,      &
       detrainment,      &
       itest,    jtest,  &
       rtest,            &
       iter_melt,        &  !WHL - debug
       D_plume_old,      &
       D_plume,          &
       eta_plume)

    !--------------------------------------------------------------------
    ! Solve the continuity equation for D_plume:
    ! 
    !    dD/dt = e - d - del*(Du)
    !
    ! where e is entrainment and d is detrainment.
    !  
    ! The equation is solved in delta form:
    !
    !    delta_D = (D_old - D_cur) + dt*(e-d) - dt*div(Du)
    !
    ! where D_cur is the current guess for D, passed into the subroutine as D_plume.
    ! 
    ! The current velocity u_cur = (u_plume,v_plume) is also passed into the subroutine.
    ! Given u_cur and D_cur, the divergence term is expanded to first order as
    ! 
    !    div(Du) = div(D_cur*u_cur + u_cur*delta_D + D_cur*delta_u)
    !
    ! where 
    !    delta_u = (du/d_eta) * delta_eta
    !    delta_eta = delta_D where eta > 0, and delta_eta = 0 otherwise.
    !
    ! The terms containing delta_D are moved to the LHS and inserted in a matrix,
    ! giving a problem of the form
    !
    !    A*delta_D = rhs
    !
    ! This subroutine is called repeatedly until the residual is sufficiently small.
    !
    ! For now, I am using SLAP to solve the  matrix.
    ! Later, I plan to use a homegrown parallel solver.
    !--------------------------------------------------------------------

    ! for sparse_easy_solve
    use glimmer_sparse_type
    use glimmer_sparse

    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)

    real(dp), intent(in) ::  &
         dt_plume               ! time step for plume solver (s)

    logical, intent(in) ::  &
         free_surface           ! true if computing PG force due to slope in free surface

    character(len=6), intent(in) ::  &
         nonlinear_method       ! method for solving nonlinear equations, 'Picard' or 'Newton' 

    integer, dimension(nx,ny), intent(in) ::  &
         plume_mask_cell,     & ! = 1 for cells where scalar plume variables are computed, else = 0
         edge_mask_east,      & ! = 1 on east edges where plume velocity is computed, else = 0
         edge_mask_north,     & ! = 1 on north edges where plume velocity is computed, else = 0
         divu_mask_east,      & ! = 1 on east edges where divergence terms are computed, else = 0
         divu_mask_north        ! = 1 on north edges where divergence terms are  computed, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         u_plume_east,        & ! u_plume on east edges (m/s)
         v_plume_east,        & ! v_plume on east edges (m/s)
         u_plume_north,       & ! u_plume on north edges (m/s)
         v_plume_north,       & ! v_plume on north edges (m/s)
         du_deta_west,        & ! dependence of u on eta(i,j), west of the velocity point
         du_deta_east,        & ! dependence of u on eta(i+1,j), east of the velocity point
         du_deta_northwest,   & ! dependence of u on eta(i,j+1), northwest of the velocity point
         du_deta_northeast,   & ! dependence of u on eta(i+1,j+1), northeast of the velocity point
         du_deta_southwest,   & ! dependence of u on eta(i,j-1), southwest of the velocity point
         du_deta_southeast,   & ! dependence of u on eta(i+1,j-1), southeast of the velocity point
         dv_deta_south,       & ! dependence of v on eta(i,j), south of the velocity point
         dv_deta_north,       & ! dependence of v on eta(i,j+1), north of the velocity point
         dv_deta_northwest,   & ! dependence of v on eta(i-1,j+1), northwest of the velocity point
         dv_deta_northeast,   & ! dependence of v on eta(i+1,j+1), northeast of the velocity point
         dv_deta_southwest,   & ! dependence of v on eta(i-1,j), southwest of the velocity point
         dv_deta_southeast,   & ! dependence of v on eta(i+1,j), southeast of the velocity point
         plume_speed_east,    & ! plume speed on east edges (m/s)
         plume_speed_north,   & ! plume speed on north edges (m/s)
         H_cavity,            & ! ocean cavity thickness (m), lsrf - topg
         entrainment,         & ! entrainment at cell centers (m/s)
         detrainment            ! detrainment at cell centers (m/s)

    integer, intent(in) :: &
         itest, jtest, rtest    ! diagnostic indices

    integer, intent(in) :: iter_melt  !WHL - debug

    real(dp), dimension(nx,ny), intent(in) ::  &
         D_plume_old            ! old plume thickness (m) at the previous time step

    real(dp), dimension(nx,ny), intent(inout) ::  &
         D_plume                ! on input, the current guess for the plume thickness (m)
                                ! on output, the new guess for the plume thickness

    real(dp), dimension(nx,ny), intent(out) ::  &
         eta_plume              ! displacement of plume surface, max(D_plume - H_cavity, 0.0)

    ! local variables

    !WHL - Remove _up variables?
    real(dp), dimension(nx,ny) ::  &
         D_plume_latest,      & ! D_plume from the most recent iteration (m)
         D_plume_cap,         & ! min(D_plume, H_cavity)
         D_plume_east,        & ! plume thickness at each east edge (m)
         D_plume_north,       & ! plume thickness at each north edge (m)
         D_plume_east_up,     & ! upstream plume thickness at each east edge (m)
         D_plume_north_up       ! upstream plume thickness at each north edge (m)

    real(dp), dimension(nx,ny) ::  &
         c_east,              & ! term in du/deta at east edges, proportional to c_drag
         f_east,              & ! term in du/deta at east edges, proportional to f_coriolis
         c_north,             & ! term in du/deta at north edges, proportional to c_drag
         f_north                ! term in du/deta at north edges, proportional to f_coriolis

    integer, dimension(nx,ny) ::  &
         upos_mask,   & ! = 1 at edges where u_plume > 0, else = 0
         uneg_mask,   & ! = 1 at edges where u_plume < 0, else = 0
         vpos_mask,   & ! = 1 at edges where v_plume > 0, else = 0
         vneg_mask      ! = 1 at edges where v_plume < 0, else = 0

    real(dp), dimension(-1:1,-1:1,nx,ny) ::  &
         A_plume        ! array holding nonzero matrix elements on the structured mesh
                        ! up to 9 nonzero elements per row of the matrix

    type(sparse_matrix_type) ::  &
         matrix         ! sparse matrix for SLAP solver, defined in glimmer_sparse_types
                        ! includes nonzeroes, order, col, row, val

    real(dp), dimension(:), allocatable ::  &
         rhs,         & ! right-hand-side vector, passed to solver
         answer         ! answer vector, returned from solver

    real(dp) ::  &
         err            ! solution error the solver

    integer, dimension(nx,ny) ::  &
         cellID         ! integer ID for each cell

    integer, dimension(nx*ny) ::  &
         iCellIndex, jCellIndex   ! indices for mapping cellID back to i and j

    real(dp) :: &
         D_east, D_north,      & ! current estimate of D_plume, averaged to east and north edges
         denom                   ! denominator in the expression for velocity

    integer :: niters   ! iteration counter
    integer :: i, j     ! horizontal indices
    integer :: iA, jA   ! horizontal index shifts, in range -1:1
    integer :: n        ! matrix index
    integer :: count    ! counter
    integer :: matrix_order  ! size of square matrix
    integer :: max_nonzeros  ! max number of nonzero elements in matrix

    real(dp) :: &
         uu_term, vv_term, uv_term

    ! SLAP linear solver (BICG or GMRES)
    ! For ISOMIP+, BICG is a bit faster, requiring 3 or 4 linear iterations compared to 6 or 7 for GMRES.
    !TODO - Replace with homegrown solver

    integer, parameter :: &
         whichsparse = HO_SPARSE_BICG
!!         whichsparse = HO_SPARSE_GMRES

    real(dp), parameter :: &
         relax_D = 0.5d0,       &  !WHL - Remove relax_D?
!!         relax_eta = 0.20d0
         relax_eta = 0.50d0  ! relaxation parameter to prevent D from oscillating about H_cavity
                             ! 0 < relax_eta <=1; a smaller value gives stronger overrelaxation.
                             ! In practice, a value of ~0.5 seems to work well.
                             ! With smaller values, convergence is slower but may be more robust.

    !WHL - debug
    real(dp) :: dDu_dx, dDv_dy
    real(dp) :: diag_east, diag_west, diag_north, diag_south
    real(dp) :: offdiag_east, offdiag_west, offdiag_north, offdiag_south
    real(dp) :: diag_north_mod, offdiag_north_mod
    real(dp) :: diag_south_mod, offdiag_south_mod
    real(dp) :: diag_north_mod1a, diag_north_mod1b, diag_north_mod2a, diag_north_mod2b
    real(dp) :: diag_north_mod3, diag_north_mod4, diag_north_mod5, diag_north_mod6
    real(dp) :: diag_south_mod1a, diag_south_mod1b, diag_south_mod2a, diag_south_mod2b
    real(dp) :: diag_south_mod3, diag_south_mod4, diag_south_mod5, diag_south_mod6
    real(dp), dimension(nx,ny) :: eta_plume_latest
    real(dp) :: uu_north, vv_north, uv_north, uu_south, vv_south, uv_south

    logical, parameter :: apply_jacobian = .true.
!!    logical, parameter :: apply_jacobian = .false.

    print*, ' '
    print*, 'In plume_thickness_solver: itest, jtest =', itest, jtest

    ! count plume cells in matrix solve
    ! loop over locally owned cells
    count = 0
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j)==1) then
             count = count + 1
             cellID(i,j) = count
             iCellIndex(count) = i
             jCellIndex(count) = j
          endif
       enddo
    enddo
   
    ! initialize and allocate
    matrix_order = count
    max_nonzeros = matrix_order * 9

    allocate(matrix%row(max_nonzeros), matrix%col(max_nonzeros), matrix%val(max_nonzeros))
    allocate(rhs(matrix_order), answer(matrix_order))

    A_plume(:,:,:,:) = 0.0d0
    rhs(:) = 0.0d0
    answer(:) = 0.0d0

    ! Save the latest iterate for D_plume.
    ! Note: This is different from D_plume_old, which is the plume thickness at the start of the time step.
    D_plume_latest(:,:) = D_plume(:,:)

    ! Given the latest iterate for D_plume, compute eta_plume
    eta_plume(:,:) = 0.0d0
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             eta_plume(i,j) = max(D_plume(i,j) - H_cavity(i,j), 0.0d0)

             if (i==itest .and. j==jtest) then
                print*, 'Current D, eta =', D_plume(i,j), eta_plume(i,j)
             endif

          endif
       enddo
    enddo

    !WHL - debug
    eta_plume_latest(:,:) = eta_plume(:,:)

    ! Set the capped value of D_plume, for purposes of computing the divergence.
    if (cap_Dplume) then
       D_plume_cap(:,:) = min(D_plume(:,:), H_cavity(:,:))
    else
       D_plume_cap(:,:) = D_plume(:,:)
    endif

    ! Compute D_plume edge values used in the finite difference expression for divergence.
    ! No solution is ideal here:
    ! With a centered difference for D_plume on edges, thin plumes upstream from thick plumes can overempty.
    ! With an upstream difference, there is no overemptying, but there can be oscillations as a result
    !  of u or v changing sign between one iteration and the next.
    ! Choosing the min value is not very accurate (only first order, and too restrictive of outflow when
    !  a thick plume is upstream of a thin plume), but it inhibits overemptying and oscillations.

!!    ! Store for each edge the plume thickness in the upstream cell.

    ! Note: The use of edge_mask and divu_mask is a bit subtle.
    !       If edge_mask = 1, the plume exists on both sides of the edge, and one cell is clearly upstream.
    !       If edge_mask = 0 but divu_mask = 1, the plume exists on only one side of the edge.
    !        We identify an upstream plume thickness only if the plume exists in the upstream cell
    !        and the flow is out of the plume domain.
    !        (This generally is the case, but we check here to be sure.)

    upos_mask(:,:) = 0
    uneg_mask(:,:) = 0
    vpos_mask(:,:) = 0
    vneg_mask(:,:) = 0

    D_plume_east(:,:) = 0.0d0
    D_plume_north(:,:) = 0.0d0

    !WHL - Remove _up variables?
    D_plume_east_up(:,:) = 0.0d0
    D_plume_north_up(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx

!!          D_plume_east(i,j) = 0.5d0 * (D_plume(i,j) + D_plume(i+1,j))
!!          D_plume_north(i,j) = 0.5d0 * (D_plume(i,j) + D_plume(i,j+1))
!          if (divu_mask_east(i,j) == 1) then
!             D_plume_east(i,j) = 0.5d0 * (D_plume_cap(i,j) + D_plume_cap(i+1,j))
!          else
!             D_plume_east(i,j) = 0.0d0
!          endif

!          if (divu_mask_north(i,j) == 1) then
!             D_plume_north(i,j) = 0.5d0 * (D_plume_cap(i,j) + D_plume_cap(i,j+1))
!          else
!             D_plume_north(i,j) = 0.0d0
!          endif

          ! mark the u component as positive or negative
          !WHL - Remove these calculations?
          if (edge_mask_east(i,j) == 1) then ! plume exists in both neighbor cells

             D_plume_east(i,j) = min(D_plume_cap(i,j), D_plume_cap(i+1,j))

!             if (u_plume_east(i,j) > 0.0d0) then
!                upos_mask(i,j) = 1
!                D_plume_east_up(i,j) = D_plume_cap(i,j)
!             elseif (u_plume_east(i,j) < 0.0d0) then
!                uneg_mask(i,j) = 1
!                D_plume_east_up(i,j) = D_plume_cap(i+1,j)
!             endif

          elseif (divu_mask_east(i,j) == 1) then  ! plume exists on only one side of the edge

             if (plume_mask_cell(i,j) == 1) then
                D_plume_east(i,j) = 0.5d0 * D_plume_cap(i,j)
             else
                D_plume_east(i,j) = 0.5d0 * D_plume_cap(i+1,j)
             endif

!             if (plume_mask_cell(i,j) == 1 .and. u_plume_east(i,j) > 0.0d0) then
!                upos_mask(i,j) = 1
!                D_plume_east_up(i,j) = D_plume_cap(i,j)
!             elseif (plume_mask_cell(i+1,j) == 1 .and. u_plume_east(i,j) < 0.0d0) then
!                uneg_mask(i,j) = 1
!                D_plume_east_up(i,j) = D_plume_cap(i+1,j)
!             endif

          endif

          ! mark the v component as positive or negative
          if (edge_mask_north(i,j) == 1) then

             D_plume_north(i,j) = min(D_plume_cap(i,j), D_plume_cap(i,j+1))

!             if (v_plume_north(i,j) > 0.0d0) then
!                vpos_mask(i,j) = 1
!                D_plume_north_up(i,j) = D_plume_cap(i,j)
!             elseif (v_plume_north(i,j) < 0.0d0) then
!                vneg_mask(i,j) = 1
!                D_plume_north_up(i,j) = D_plume_cap(i,j+1)
!             endif

          elseif (divu_mask_north(i,j) == 1) then  ! plume exists on only one side of the edge

             if (plume_mask_cell(i,j) == 1) then
                D_plume_north(i,j) = 0.5d0 * D_plume_cap(i,j)
             else
                D_plume_north(i,j) = 0.5d0 * D_plume_cap(i,j+1)
             endif

!             if (plume_mask_cell(i,j) == 1 .and. v_plume_north(i,j) > 0.0d0) then
!                vpos_mask(i,j) = 1
!                D_plume_north_up(i,j) = D_plume_cap(i,j)
!             elseif (plume_mask_cell(i,j+1) == 1 .and. v_plume_north(i,j) < 0.0d0) then
!                vneg_mask(i,j) = 1
!                D_plume_north_up(i,j) = D_plume_cap(i,j+1)
!             endif

          endif

       enddo  ! i
    enddo  !j

    ! Compute terms in the expression for du/d_eta.  These terms contain D_plume but are treated as constant,
    ! the dependence of u on D (where eta = 0) is much weaker than the dependence on eta (where eta > 0).
    ! Note: These terms need to be set to 0.0 for masked-out edges with zero velocity.

    c_east(:,:) = 0.0d0
    f_east(:,:) = 0.0d0
    c_north(:,:) = 0.0d0
    f_north(:,:) = 0.0d0

    print*, 'trim(nonlinear_method) = ', trim(nonlinear_method)

    if (trim(nonlinear_method) == 'Newton') then

       ! compute some matrix terms associated with the dependence of plume speed on eta_plume.

       !TODO -  Check the plume_speed values at the edge of the plume?
       !TODO - Modify for divu_mask_east/north
       !    Note: edge_mask_east = 1 implies that the plume exists in cells (i,j) and (i+1,j).
       !          edge_mask_north = 1 implies that the plume exists in cells (i,j) and (i,j+1).
       !TODO - Think about whether we should extrapolate c and f terms at such edges
       !       from the neighboring edges with edge_mask = 1.
       !TODO - Think about whether to use D_plume_cap for D_east and D_north.
       !       Maybe so, provided we use D_plume_cap in the velocity calculation.
       
       do j = 1, ny
          do i = 1, nx
          
             !TODO - Could replace D_east by D_plume_east(i,j)
             if (divu_mask_east(i,j) == 1) then
!!                D_east = 0.5d0 * (D_plume(i,j) + D_plume(i+1,j))
                D_east = 0.5d0 * (D_plume_cap(i,j) + D_plume_cap(i+1,j))
                denom = (c_drag*plume_speed_east(i,j))**2 + (f_coriolis*D_east)**2
                c_east(i,j) = c_drag * plume_speed_east(i,j) * grav * D_east / (dx * denom)
                f_east(i,j) = f_coriolis * grav * D_east**2 / (4.0d0 * dy * denom)
             endif
             
             !TODO - Could replace D_north by D_plume_north(i,j)
             if (divu_mask_north(i,j) == 1) then
!!                D_north = 0.5d0 * (D_plume(i,j) + D_plume(i,j+1))
                D_north = 0.5d0 * (D_plume_cap(i,j) + D_plume_cap(i,j+1))
                denom = ((c_drag*plume_speed_north(i,j))**2 + (f_coriolis*D_north)**2)
                c_north(i,j) = c_drag * plume_speed_north(i,j) * grav * D_north / (dy * denom)
                f_north(i,j) = f_coriolis * grav * D_north**2 / (4.0d0 * dx * denom)
             endif
             
             if (i==itest .and. j==jtest) then
                print*, 'i, j, denom:', i, j, denom
                print*, 'c_west, c_east:', c_east(i-1,j), c_east(i,j)
                print*, 'c_south, c_north:', c_north(i,j-1), c_north(i,j)
                print*, 'f_east, f_north:', f_east(i,j), f_north(i,j)
             endif
             
          enddo
       enddo

    endif  ! nonlinear method = Newton

    !--------------------------------------------------------------------
    ! Solve the equation dD/dt = e - d - div(Du).
    ! This is done in delta form:
    ! delta_D = (D_old - D_cur) + dt*(e-d) - dt*div(Du)
    !    where D_cur is the current guess for D.
    ! This subroutine is called repeatedly until the residual is sufficiently small.
    !
    !--------------------------------------------------------------------
 
    ! compute nonzero matrix elements
    ! loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             ! right-hand side
             ! This term includes
             ! (1) D_old - D_cur, the difference between D_plume_old and the current guess
             ! (2) (e - d)*dt
             ! (3) -dt * div(D_cur*u_cur), a divergence term based on the current guesses for D and u
             ! Note: The finite-difference expression for the divergence uses upstream plume thicknesses
             !       to avoid overemptying cells.

             !TODO - May not need divu_mask_east, if u_plume_east = 0 at edges with divu_mask_east = 0
             n = cellID(i,j)
!             rhs(n) = D_plume_old(i,j) - D_plume(i,j)  &
!                    + (entrainment(i,j) - detrainment(i,j))*dt_plume  &
!                    - (dt_plume/dx) * (D_plume_east_up(i,j)*u_plume_east(i,j)*divu_mask_east(i,j) &
!                                     - D_plume_east_up(i-1,j)*u_plume_east(i-1,j)*divu_mask_east(i-1,j)) &
!                    - (dt_plume/dy) * (D_plume_north_up(i,j)*v_plume_north(i,j)*divu_mask_north(i,j) &
!                                     - D_plume_north_up(i,j-1)*v_plume_north(i,j-1)*divu_mask_north(i,j-1))
             rhs(n) = D_plume_old(i,j) - D_plume(i,j)  &
                    + (entrainment(i,j) - detrainment(i,j))*dt_plume  &
                    - (dt_plume/dx) * (D_plume_east(i,j)*u_plume_east(i,j)*divu_mask_east(i,j) &
                                     - D_plume_east(i-1,j)*u_plume_east(i-1,j)*divu_mask_east(i-1,j)) &
                    - (dt_plume/dy) * (D_plume_north(i,j)*v_plume_north(i,j)*divu_mask_north(i,j) &
                                     - D_plume_north(i,j-1)*v_plume_north(i,j-1)*divu_mask_north(i,j-1))

             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, '(D - D_old)/dt =', (D_plume(i,j) - D_plume_old(i,j)) / dt_plume
                print*, '(e - d) =', entrainment(i,j) - detrainment(i,j)

                dDu_dx = (1.0d0/dx) * (D_plume_east(i,j)*u_plume_east(i,j)*divu_mask_east(i,j) &
                                     - D_plume_east(i-1,j)*u_plume_east(i-1,j)*divu_mask_east(i-1,j))
                dDv_dy = (1.0d0/dy) * (D_plume_north(i,j)*v_plume_north(i,j)*divu_mask_north(i,j) &
                                     - D_plume_north(i,j-1)*v_plume_north(i,j-1)*divu_mask_north(i,j-1))
                print*, '-div(Du) =', -dDu_dx - dDv_dy
                print*, 'u_plume_west/east =', u_plume_east(i-1,j), u_plume_east(i,j)
!!                print*, 'u_plume_south/north =', u_plume_north(i,j-1), u_plume_north(i,j)
!!                print*, 'v_plume_west/east =', v_plume_east(i-1,j), v_plume_east(i,j)
                print*, 'v_plume_south/north =', v_plume_north(i,j-1), v_plume_north(i,j)
                print*, 'D_plume(i-1,j),(i,j) =', D_plume(i-1,j), D_plume(i,j)
                print*, 'D_plume(i,j),(i,j+1) =', D_plume(i,j), D_plume(i,j+1)
                print*, 'D_plume_west/east =', D_plume_east(i-1,j), D_plume_east(i,j)
                print*, 'D_plume_south/north =', D_plume_north(i,j-1), D_plume_north(i,j)
                print*, 'dDu_dx, dDv_dy =', dDu_dx, dDv_dy

                print*, 'residual (m/s) =', rhs(n)/dt_plume
             endif

             ! initialize the matrix diagonal
             A_plume(0,0,i,j) = 1.0d0

             ! Add matrix terms associated with u_latest*delta_D
             ! Note: The upos, uneg, vpos and vneg masks are constructed such that 
             !       if a mask = 1 at an edge, then the plume exists in the upstream cell.

             ! diagonal element
!             A_plume(0,0,i,j) = A_plume(0,0,i,j)  &
!                              + (dt_plume/dx) * u_plume_east(i,j) * upos_mask(i,j)  &
!                              - (dt_plume/dx) * u_plume_east(i-1,j) * uneg_mask(i-1,j)  &
!                              + (dt_plume/dy) * v_plume_north(i,j) * vpos_mask(i,j)  &
!                              - (dt_plume/dy) * v_plume_north(i,j-1) * vneg_mask(i,j-1)
!             A_plume(0,0,i,j) = A_plume(0,0,i,j)  &
!                              + 0.5d0 * (dt_plume/dx) * u_plume_east(i,j) * divu_mask_east(i,j) &
!                              - 0.5d0 * (dt_plume/dx) * u_plume_east(i-1,j) * divu_mask_east(i-1,j) &
!                              + 0.5d0 * (dt_plume/dy) * v_plume_north(i,j) * divu_mask_north(i,j) &
!                              - 0.5d0 * (dt_plume/dy) * v_plume_north(i,j-1) * divu_mask_north(i,j-1)
             if (D_plume(i,j) == D_plume_east(i,j)) then
                A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dx) * u_plume_east(i,j)
             endif

             if (D_plume(i,j) == D_plume_east(i-1,j)) then
                A_plume(0,0,i,j) = A_plume(0,0,i,j) - (dt_plume/dx) * u_plume_east(i-1,j)
             endif

             if (D_plume(i,j) == D_plume_north(i,j)) then
                A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dy) * v_plume_north(i,j)
             endif

             if (D_plume(i,j) == D_plume_north(i,j-1)) then
                A_plume(0,0,i,j) = A_plume(0,0,i,j) - (dt_plume/dy) * v_plume_north(i,j-1)
             endif

             ! off-diagonal elements
!             A_plume(1,0,i,j)  = A_plume(1,0,i,j)  &
!                               + (dt_plume/dx) * u_plume_east(i,j) * uneg_mask(i,j)
!             A_plume(-1,0,i,j) = A_plume(-1,0,i,j)  &
!                               - (dt_plume/dx) * u_plume_east(i-1,j) * upos_mask(i-1,j)
!             A_plume(0,1,i,j)  = A_plume(0,1,i,j)  &
!                               + (dt_plume/dy) * v_plume_north(i,j) * vneg_mask(i,j)
!             A_plume(0,-1,i,j) = A_plume(0,-1,i,j)  &
!                               - (dt_plume/dy) * v_plume_north(i,j-1) * vpos_mask(i,j-1)

             ! Note: If plume_mask_cell in the neighboring cell = 0, then D_plume in that cell is fixed at 0.
             !       In that case, there is no matrix element associated with the cell.
             ! TODO - May not need divu_mask_east if u_plume_east = 0 wherever divu_mask_east = 0
             if (D_plume(i+1,j) == D_plume_east(i,j) .and. plume_mask_cell(i+1,j) == 1) then
                A_plume(1,0,i,j)  = A_plume(1,0,i,j) + (dt_plume/dx) * u_plume_east(i,j)
             endif

             if (D_plume(i-1,j) == D_plume_east(i-1,j) .and. plume_mask_cell(i-1,j) == 1) then
                A_plume(-1,0,i,j) = A_plume(-1,0,i,j) - (dt_plume/dx) * u_plume_east(i-1,j)
             endif

             if (D_plume(i,j+1) == D_plume_north(i,j) .and. plume_mask_cell(i,j+1) == 1) then
                A_plume(0,1,i,j)  = A_plume(0,1,i,j) + (dt_plume/dy) * v_plume_north(i,j)
             endif

             if (D_plume(i,j-1) == D_plume_north(i,j-1) .and. plume_mask_cell(i,j-1) == 1) then
                A_plume(0,-1,i,j) = A_plume(0,-1,i,j) - (dt_plume/dy) * v_plume_north(i,j-1)
             endif

             if (i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, ' '
                print*, 'i, j, A_plume, rhs:', &
                     i, j, A_plume(0,0,i,j), A_plume(1,0,i,j), A_plume(-1,0,i,j), A_plume(0,1,i,j), A_plume(0,-1,i,j), &
                     rhs(cellID(i,j))
             endif

             if (trim(nonlinear_method) == 'Newton') then

                ! Compute terms associated with D_latest*delta_u
                ! These terms take into account the change of u with changes in eta.
                ! I.e., we let delta_u = du/deta * delta_eta
                ! Here, delta_eta = delta_D only where eta > 0.  Elsewhere, delta_eta = 0.
                !
                ! Note: Some factors of dx and dy in the denominator have already been incorporated
                !        in the c and f terms (including the factor of 1/4 in the f terms).
                !       The dx and dy below are associated with taking the divergence.
                ! Note: c_east and f_east are nonzero only if divu_mask_east is nonzero.
                !       c_north and f_north are nonzero only if divu_mask_north is nonzero.
                !
                ! Note: Because of grav in the numerator, these matrix elements are much larger
                !       than those computed above. This means that once D >= H_cavity, changes
                !       in eta are small. If a layer with D < H_cavity inflates to have eta > 0,
                !       we set eta = 0 on this iteration to prevent unrealistic large value of eta.
                !       Inflation can continue at a slower rate on subsequent iterations.
                
                ! TODO: Deal with oscillations of D in the neighborhood of H_cav?

                ! add c terms
                ! WHL - These help a lot with convergence

!!              go to 400

                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'Add du/deta and dv/deta terms:'
                endif

                ! East edge

                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'i, j, c_east, du_deta_west, du_deta_east:', &
                        i, j, c_east(i,j), du_deta_west(i,j), du_deta_east(i,j)
                endif
                
                if (apply_jacobian) then

                   if (D_plume(i,j) >= H_cavity(i,j)) then
!!                   if (plume_mask_cell(i,j) == 1) then
                      A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dx) * D_plume_east(i,j) * du_deta_west(i,j)
                   endif
                  
                   if (D_plume(i+1,j) >= H_cavity(i+1,j)) then
!!                   if (plume_mask_cell(i+1,j) == 1) then
                      A_plume(1,0,i,j) = A_plume(1,0,i,j) + (dt_plume/dx) * D_plume_east(i,j) * du_deta_east(i,j)
                   endif

                   if (D_plume(i,j+1) >= H_cavity(i,j+1)) then
                      A_plume(0,1,i,j) = A_plume(0,1,i,j) + (dt_plume/dx) * D_plume_east(i,j) * du_deta_northwest(i,j)
                   endif

                   if (D_plume(i+1,j+1) >= H_cavity(i+1,j+1)) then
                      A_plume(1,1,i,j) = A_plume(1,1,i,j) + (dt_plume/dx) * D_plume_east(i,j) * du_deta_northeast(i,j)
                   endif

                   if (D_plume(i,j-1) >= H_cavity(i,j-1)) then
                      A_plume(0,-1,i,j) = A_plume(0,-1,i,j) + (dt_plume/dx) * D_plume_east(i,j) * du_deta_southwest(i,j)
                   endif

                   if (D_plume(i+1,j-1) >= H_cavity(i+1,j-1)) then
                      A_plume(1,-1,i,j) = A_plume(1,-1,i,j) + (dt_plume/dx) * D_plume_east(i,j) * du_deta_southeast(i,j)
                   endif

                else

                   if (D_plume(i,j) >= H_cavity(i,j)) then
                      A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dx) * D_plume_east(i,j) * c_east(i,j)
                   endif
                
                   if (D_plume(i+1,j) >= H_cavity(i+1,j)) then
                      A_plume(1,0,i,j) = A_plume(1,0,i,j) - (dt_plume/dx) * D_plume_east(i,j) * c_east(i,j)
                   endif
                   
                endif  ! apply_jacobian

                ! West edge
                ! These are the same as for the east edge, but with all i indices shifted to i-1
                ! Note: Minus sign applies because an increase of u(i-1,j) implies a decrease in the divergence.
                ! first the 2 terms in my original expression, without the dependence on U = plume_speed_east
                
                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'i-1, j, c_east, du_deta_west, du_deta_east:', &
                        i-1, j, c_east(i-1,j), du_deta_west(i-1,j), du_deta_east(i-1,j)
                endif
                
                if (apply_jacobian) then

                   if (D_plume(i-1,j) >= H_cavity(i-1,j)) then
!!                   if (plume_mask_cell(i-1,j) == 1) then
                      A_plume(-1,0,i,j) = A_plume(-1,0,i,j) - (dt_plume/dx) * D_plume_east(i-1,j) * du_deta_west(i-1,j)
                   endif

                   if (D_plume(i,j) >= H_cavity(i,j)) then
!!                   if (plume_mask_cell(i,j) == 1) then
                      A_plume(0,0,i,j) = A_plume(0,0,i,j) - (dt_plume/dx) * D_plume_east(i-1,j) * du_deta_east(i-1,j)
                   endif
             
                   if (D_plume(i-1,j+1) >= H_cavity(i-1,j+1)) then
                      A_plume(-1,1,i,j) = A_plume(-1,1,i,j) - (dt_plume/dx) * D_plume_east(i-1,j) * du_deta_northwest(i-1,j)
                   endif

                   if (D_plume(i,j+1) >= H_cavity(i,j+1)) then
                      A_plume(0,1,i,j) = A_plume(0,1,i,j) - (dt_plume/dx) * D_plume_east(i-1,j) * du_deta_northeast(i-1,j)
                   endif

                   if (D_plume(i-1,j-1) >= H_cavity(i-1,j-1)) then
                      A_plume(-1,-1,i,j) = A_plume(-1,-1,i,j) - (dt_plume/dx) * D_plume_east(i-1,j) * du_deta_southwest(i-1,j)
                   endif

                   if (D_plume(i,j-1) >= H_cavity(i,j-1)) then
                      A_plume(0,-1,i,j) = A_plume(0,-1,i,j) - (dt_plume/dx) * D_plume_east(i-1,j) * du_deta_southeast(i-1,j)
                   endif
                   
                else

                   if (D_plume(i,j) >= H_cavity(i,j)) then
                      A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dx) * D_plume_east(i-1,j) * c_east(i-1,j)
                   endif
             
                   if (D_plume(i-1,j) >= H_cavity(i-1,j)) then
                      A_plume(-1,0,i,j) = A_plume(-1,0,i,j) - (dt_plume/dx) * D_plume_east(i-1,j) * c_east(i-1,j)
                   endif
                   
                endif  ! apply_jacobian

                ! North edge

                ! first the 2 terms in my original expression, without the dependence on U = plume_speed_east

                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'i, j, c_north, dv_deta_south, dv_deta_north:', &
                        i, j, c_north(i,j), dv_deta_south(i,j), dv_deta_north(i,j)
                endif

                if (apply_jacobian) then

                   if (D_plume(i,j) >= H_cavity(i,j)) then
!!                   if (plume_mask_cell(i,j) == 1) then
                      A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dy) * D_plume_north(i,j) * dv_deta_south(i,j)
                   endif
                   
                   if (D_plume(i,j+1) >= H_cavity(i,j+1)) then
!!                   if (plume_mask_cell(i,j+1) == 1) then
                      A_plume(0,1,i,j) = A_plume(0,1,i,j) + (dt_plume/dy) * D_plume_north(i,j) * dv_deta_north(i,j)
                   endif

                   if (D_plume(i-1,j) >= H_cavity(i-1,j)) then
                      A_plume(-1,0,i,j) = A_plume(-1,0,i,j) + (dt_plume/dy) * D_plume_east(i,j) * dv_deta_southwest(i,j)
                   endif

                   if (D_plume(i-1,j+1) >= H_cavity(i-1,j+1)) then
                      A_plume(-1,1,i,j) = A_plume(-1,1,i,j) + (dt_plume/dy) * D_plume_east(i,j) * dv_deta_northwest(i,j)
                   endif

                   if (D_plume(i+1,j) >= H_cavity(i+1,j)) then
                      A_plume(1,0,i,j) = A_plume(1,0,i,j) + (dt_plume/dy) * D_plume_east(i,j) * dv_deta_southeast(i,j)
                   endif

                   if (D_plume(i+1,j+1) >= H_cavity(i+1,j+1)) then
                      A_plume(1,1,i,j) = A_plume(1,1,i,j) + (dt_plume/dy) * D_plume_east(i,j) * dv_deta_northeast(i,j)
                   endif

                else

                   if (D_plume(i,j) >= H_cavity(i,j)) then
                      A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dy) * D_plume_north(i,j) * c_north(i,j)
                   endif
                   
                   if (D_plume(i,j+1) >= H_cavity(i,j+1)) then
                      A_plume(0,1,i,j) = A_plume(0,1,i,j) - (dt_plume/dy) * D_plume_north(i,j) * c_north(i,j)
                   endif

                endif  ! apply_jacobian

                ! South edge
                
                ! first the 2 terms in my original expression, without the dependence on U = plume_speed_east
                ! Note: Minus sign applies because an increase of v(i,j-1) implies a decrease in the divergence.

                if (i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'i, j-1, c_north, dv_deta_south, dv_deta_north:', &
                        i, j, c_north(i,j-1), dv_deta_south(i,j-1), dv_deta_north(i,j-1)
                endif

                if (apply_jacobian) then

                   if (D_plume(i,j) >= H_cavity(i,j)) then
!!                   if (plume_mask_cell(i,j) == 1) then
                      A_plume(0,0,i,j) = A_plume(0,0,i,j) - (dt_plume/dy) * D_plume_north(i,j-1) * dv_deta_north(i,j-1)
                   endif
                   
                   if (D_plume(i,j-1) >= H_cavity(i,j-1)) then
!!                   if (plume_mask_cell(i,j-1) == 1) then
                      A_plume(0,-1,i,j) = A_plume(0,-1,i,j) - (dt_plume/dy) * D_plume_north(i,j-1) * dv_deta_south(i,j-1)
                   endif

                   if (D_plume(i-1,j) >= H_cavity(i-1,j)) then
                      A_plume(-1,0,i,j) = A_plume(-1,0,i,j) + (dt_plume/dy) * D_plume_east(i,j) * dv_deta_northwest(i,j-1)
                   endif

                   if (D_plume(i-1,j-1) >= H_cavity(i-1,j-1)) then
                      A_plume(-1,-1,i,j) = A_plume(-1,-1,i,j) + (dt_plume/dy) * D_plume_east(i,j) * dv_deta_southwest(i,j-1)
                   endif

                   if (D_plume(i+1,j) >= H_cavity(i+1,j)) then
                      A_plume(1,0,i,j) = A_plume(1,0,i,j) + (dt_plume/dy) * D_plume_east(i,j) * dv_deta_northeast(i,j-1)
                   endif

                   if (D_plume(i+1,j-1) >= H_cavity(i+1,j-1)) then
                      A_plume(1,-1,i,j) = A_plume(1,-1,i,j) + (dt_plume/dy) * D_plume_east(i,j) * dv_deta_southeast(i,j-1)
                   endif

                else

                   if (D_plume(i,j) >= H_cavity(i,j)) then
                      A_plume(0,0,i,j) = A_plume(0,0,i,j) + (dt_plume/dy) * D_plume_north(i,j-1) * c_north(i,j-1)
                   endif
                   
                   if (D_plume(i,j-1) >= H_cavity(i,j-1)) then
                      A_plume(0,-1,i,j) = A_plume(0,-1,i,j) - (dt_plume/dy) * D_plume_north(i,j-1) * c_north(i,j-1)
                   endif
                
                endif  ! apply_jacobian

             endif   ! nonlinear method = Newton

             if (i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, ' '
                print*, 'After adding c terms:'
                print*, 'i, j, A_plume, rhs:', &
                     i, j, A_plume(0,0,i,j), A_plume(1,0,i,j), A_plume(-1,0,i,j), A_plume(0,1,i,j), A_plume(0,-1,i,j), &
                     rhs(cellID(i,j))
                print*, ' '
             endif

400          continue

             ! add fterms
             !WHL - Not sure that these help.  Look more closely for bugs.

             if (nonlinear_method == 'Newton') then

                go to 500  ! skip for now

                !WHL - went ahead and removed '_up' suffixes
                !TODO - Make sure these aren't missing factors of 0.5

                ! east edge

                if (D_plume(i,j) >= H_cavity(i,j)) then
                   A_plume(0,0,i,j) = A_plume(0,0,i,j) &
                        + (dt_plume/dx) * D_plume_east(i,j) * f_east(i,j) * divu_mask_north(i,j)  &
                        - (dt_plume/dx) * D_plume_east(i,j) * f_east(i,j) * divu_mask_north(i,j-1)
                endif

                if (D_plume(i+1,j) >= H_cavity(i+1,j)) then
                   A_plume(1,0,i,j) = A_plume(1,0,i,j) &
                        + (dt_plume/dx) * D_plume_east(i,j) * f_east(i,j) * divu_mask_north(i+1,j)  &
                        - (dt_plume/dx) * D_plume_east(i,j) * f_east(i,j) * divu_mask_north(i+1,j-1)
                endif
                
                if (D_plume(i,j+1) >= H_cavity(i,j+1)) then
                   A_plume(0,1,i,j) = A_plume(0,1,i,j)  &
                        - (dt_plume/dx) * D_plume_east(i,j) * f_east(i,j) * divu_mask_north(i,j)
                endif
                
                if (D_plume(i,j-1) >= H_cavity(i,j-1)) then
                   A_plume(0,-1,i,j) = A_plume(0,-1,i,j)  &
                        + (dt_plume/dx) * D_plume_east(i,j) * f_east(i,j) * divu_mask_north(i,j-1)
                endif
                
                if (D_plume(i+1,j+1) >= H_cavity(i+1,j+1)) then
                   A_plume(1,1,i,j) = A_plume(1,1,i,j)  &
                        - (dt_plume/dx) * D_plume_east(i,j) * f_east(i,j) * divu_mask_north(i+1,j)
                endif
                
                if (D_plume(i+1,j-1) >= H_cavity(i+1,j-1)) then
                   A_plume(1,-1,i,j) = A_plume(1,-1,i,j)  &
                        + (dt_plume/dx) * D_plume_east(i,j) * f_east(i,j) * divu_mask_north(i+1,j-1)
                endif
                
                ! west edge
                
                if (D_plume(i-1,j) >= H_cavity(i-1,j)) then
                   A_plume(-1,0,i,j) = A_plume(-1,0,i,j) &
                        - (dt_plume/dx) * D_plume_east(i-1,j) * f_east(i-1,j) * divu_mask_north(i-1,j)  &
                        + (dt_plume/dx) * D_plume_east(i-1,j) * f_east(i-1,j) * divu_mask_north(i-1,j-1)
                endif
                
                if (D_plume(i,j) >= H_cavity(i,j)) then
                   A_plume(0,0,i,j) = A_plume(0,0,i,j) &
                        - (dt_plume/dx) * D_plume_east(i-1,j) * f_east(i-1,j) * divu_mask_north(i,j)  &
                        + (dt_plume/dx) * D_plume_east(i-1,j) * f_east(i-1,j) * divu_mask_north(i,j-1)
                endif
                
                if (D_plume(i-1,j+1) >= H_cavity(i-1,j+1)) then
                   A_plume(-1,1,i,j) = A_plume(-1,1,i,j)  &
                        + (dt_plume/dx) * D_plume_east(i-1,j) * f_east(i-1,j) * divu_mask_north(i-1,j)
                endif
                
                if (D_plume(i-1,j-1) >= H_cavity(i-1,j-1)) then
                   A_plume(-1,-1,i,j) = A_plume(-1,-1,i,j)  &
                        - (dt_plume/dx) * D_plume_east(i-1,j) * f_east(i-1,j) * divu_mask_north(i-1,j-1)
                endif
                
                if (D_plume(i,j+1) >= H_cavity(i,j+1)) then
                   A_plume(0,1,i,j) = A_plume(0,1,i,j)  &
                        + (dt_plume/dx) * D_plume_east(i-1,j) * f_east(i-1,j) * divu_mask_north(i,j)
                endif
                
                if (D_plume(i,j-1) >= H_cavity(i,j-1)) then
                   A_plume(0,-1,i,j) = A_plume(0,-1,i,j)  &
                        - (dt_plume/dx) * D_plume_east(i-1,j) * f_east(i-1,j) * divu_mask_north(i,j-1)
                endif
                
                ! north edge
                
                if (D_plume(i,j+1) >= H_cavity(i,j+1)) then
                   A_plume(0,1,i,j) = A_plume(0,1,i,j) &
                        - (dt_plume/dy) * D_plume_north(i,j) * f_north(i,j) * divu_mask_east(i,j+1)  &
                        + (dt_plume/dy) * D_plume_north(i,j) * f_north(i,j) * divu_mask_east(i-1,j+1)
                endif
                
                if (D_plume(i,j) >= H_cavity(i,j)) then
                   A_plume(0,0,i,j) = A_plume(0,0,i,j) &
                        - (dt_plume/dy) * D_plume_north(i,j) * f_north(i,j) * divu_mask_east(i,j)  &
                        + (dt_plume/dy) * D_plume_north(i,j) * f_north(i,j) * divu_mask_east(i-1,j)
                endif
                
                if (D_plume(i+1,j+1) >= H_cavity(i+1,j+1)) then
                   A_plume(1,1,i,j) = A_plume(1,1,i,j)  &
                        + (dt_plume/dy) * D_plume_north(i,j) * f_north(i,j) * divu_mask_east(i,j+1)
                endif
                
                if (D_plume(i-1,j+1) >= H_cavity(i-1,j+1)) then
                   A_plume(-1,1,i,j) = A_plume(-1,1,i,j)  &
                        - (dt_plume/dy) * D_plume_north(i,j) * f_north(i,j) * divu_mask_east(i-1,j+1)
                endif
                
                if (D_plume(i+1,j) >= H_cavity(i+1,j)) then
                   A_plume(1,0,i,j) = A_plume(1,0,i,j)  &
                        + (dt_plume/dy) * D_plume_north(i,j) * f_north(i,j) * divu_mask_east(i,j)
                endif
                
                if (D_plume(i-1,j) >= H_cavity(i-1,j)) then
                   A_plume(-1,0,i,j) = A_plume(-1,0,i,j)  &
                        - (dt_plume/dy) * D_plume_north(i,j) * f_north(i,j) * divu_mask_east(i-1,j)
                endif
                
                ! south edge
                
                if (D_plume(i,j) >= H_cavity(i,j)) then
                   A_plume(0,0,i,j) = A_plume(0,0,i,j) &
                        + (dt_plume/dy) * D_plume_north(i,j-1) * f_north(i,j-1) * divu_mask_east(i,j)  &
                        - (dt_plume/dy) * D_plume_north(i,j-1) * f_north(i,j-1) * divu_mask_east(i-1,j)
                endif
                
                if (D_plume(i,j-1) >= H_cavity(i,j-1)) then
                   A_plume(0,-1,i,j) = A_plume(0,-1,i,j) &
                        + (dt_plume/dy) * D_plume_north(i,j-1) * f_north(i,j-1) * divu_mask_east(i,j-1)  &
                        - (dt_plume/dy) * D_plume_north(i,j-1) * f_north(i,j-1) * divu_mask_east(i-1,j-1)
                endif
                
                if (D_plume(i+1,j) >= H_cavity(i+1,j)) then
                   A_plume(1,0,i,j) = A_plume(1,0,i,j)  &
                        - (dt_plume/dy) * D_plume_north(i,j-1) * f_north(i,j-1) * divu_mask_east(i,j)
                endif
                
                if (D_plume(i-1,j) >= H_cavity(i-1,j)) then
                   A_plume(-1,0,i,j) = A_plume(-1,0,i,j)  &
                        + (dt_plume/dy) * D_plume_north(i,j-1) * f_north(i,j-1) * divu_mask_east(i-1,j)
                endif
                
                if (D_plume(i+1,j-1) >= H_cavity(i+1,j-1)) then
                   A_plume(1,-1,i,j) = A_plume(1,-1,i,j)  &
                        - (dt_plume/dy) * D_plume_north(i,j-1) * f_north(i,j-1) * divu_mask_east(i,j-1)
                endif
                
                if (D_plume(i-1,j-1) >= H_cavity(i-1,j-1)) then
                   A_plume(-1,-1,i,j) = A_plume(-1,-1,i,j)  &
                        + (dt_plume/dy) * D_plume_north(i,j-1) * f_north(i,j-1) * divu_mask_east(i-1,j-1)
                endif
                
500             continue
                
             endif  ! nonlinear method = Newton

          endif  ! plume_mask_cell
          
       enddo
    enddo

    !WHL - Put a halo update here when running in parallel

!    print*, 'min, max A:', minval(A_plume), maxval(A_plume)
!    print*, 'min, max rhs:', minval(rhs), maxval(rhs)
!    print*, 'SLAP format'

    ! place nonzero elements in SLAP matrix format
    count = 0

    do n = 1, matrix_order

       i = iCellIndex(n)
       j = jCellIndex(n)

       if (plume_mask_cell(i,j) == 1) then

          ! loop over neighbor cells that can contribute terms to this matrix row

          do jA = -1,1
             do iA = -1,1

                if (A_plume(iA,jA,i,j) /= 0.0d0) then
                   count = count + 1
                   matrix%row(count) = n
                   matrix%col(count) = cellID(i+iA,j+jA)
                   matrix%val(count) = A_plume(iA,jA,i,j)
                      
                   if (matrix%col(count) == 0) then
                      print*, 'Bad matrix column: i, j, iA, jA =', i, j, iA, jA
                      stop
                   endif
                   
                   if (j==jtest) then
!!                         print*, 'i, j, iA, jA, row, col, val, rhs:', &
!!                              i, j, iA, jA, matrix%row(count), matrix%col(count), matrix%val(count), rhs(cellID(i,j))
                   endif
                      
                endif

             enddo   ! iA
          enddo  ! jA
          
       endif  ! plume_mask_cell

    enddo  ! n

    ! Set other matrix parameters
    matrix%order = matrix_order
    matrix%nonzeros = count
    matrix%symmetric = .false.

    ! call the SLAP solver

    call sparse_easy_solve(matrix,  rhs,      answer,   &
                           err,     niters,   whichsparse)
       
    print*, 'Called sparse_easy_solve: niters, err =', niters, err

    ! Update D_plume, given delta_D_plume (in the answer vector) from the solver.
    do n = 1, matrix_order
       i = iCellIndex(n)
       j = jCellIndex(n)
       D_plume(i,j) = D_plume(i,j) + answer(n)

       if (i==itest .and. j==jtest) then
          print*, ' '
          print*, 'After solve, i, j =:', itest, jtest
          print*, ' '
          print*, 'i, j: H_cavity, D_plume_latest, D_plume:',  &
               H_cavity(i,j), D_plume_latest(i,j), D_plume(i,j)
       endif

    enddo


    eta_plume(:,:) = max(D_plume(:,:) - H_cavity(:,:), 0.0d0)

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             if (i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Before relaxation, i, j =:', itest, jtest
                print*, ' '
                print*, 'i, j: H_cavity, D_plume_latest, D_plume, eta_latest, eta:', &
                     H_cavity(i,j), D_plume_latest(i,j), D_plume(i,j), eta_plume_latest(i,j), eta_plume(i,j)
!                print*, 'deta(i-1,j+1) =', eta_plume(i-1,j+1) - eta_plume_latest(i-1,j+1)
!                print*, 'deta(i,j+1) =', eta_plume(i,j+1) - eta_plume_latest(i,j+1)
!                print*, 'deta(i+1,j+1) =', eta_plume(i+1,j+1) - eta_plume_latest(i+1,j+1)
!                print*, 'deta(i-1,j) =', eta_plume(i-1,j) - eta_plume_latest(i-1,j)
!                print*, 'deta(i,j) =', eta_plume(i,j) - eta_plume_latest(i,j)
!                print*, 'deta(i+1,j) =', eta_plume(i+1,j) - eta_plume_latest(i+1,j)
                print*, 'predicted du_east(i,j) =', &
                       du_deta_west(i,j)*(eta_plume(i,j) - eta_plume_latest(i,j)) &
                     + du_deta_east(i,j)*(eta_plume(i+1,j) - eta_plume_latest(i+1,j)) &
                     + du_deta_northwest(i,j)*(eta_plume(i,j+1) - eta_plume_latest(i,j+1)) &
                     + du_deta_northeast(i,j)*(eta_plume(i+1,j+1) - eta_plume_latest(i+1,j+1)) &
                     + du_deta_southwest(i,j)*(eta_plume(i,j-1) - eta_plume_latest(i,j-1)) &
                     + du_deta_southeast(i,j)*(eta_plume(i+1,j-1) - eta_plume_latest(i+1,j-1))
                print*, 'predicted du_east(i-1,j) =', &
                       du_deta_east(i-1,j)*(eta_plume(i,j) - eta_plume_latest(i,j)) &
                     + du_deta_west(i-1,j)*(eta_plume(i-1,j) - eta_plume_latest(i-1,j)) &
                     + du_deta_northwest(i-1,j)*(eta_plume(i-1,j+1) - eta_plume_latest(i-1,j+1)) &
                     + du_deta_northeast(i-1,j)*(eta_plume(i,j+1) - eta_plume_latest(i,j+1)) &
                     + du_deta_southwest(i-1,j)*(eta_plume(i-1,j-1) - eta_plume_latest(i-1,j+1)) &
                     + du_deta_southeast(i-1,j)*(eta_plume(i,j-1) - eta_plume_latest(i,j-1))
             endif
          endif
       enddo
    enddo

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (plume_mask_cell(i,j) == 1) then

             !WHL - new try

             ! Relax D, regardless of eta.
!!             D_plume(i,j) = D_plume_latest(i,j) + relax_D*(D_plume(i,j) - D_plume_latest(i,j))

             ! Prevent the plume from inflating to a large value in one time step,
             !  due to the absence of du/deta terms in the matrix.
!!             if (D_plume_latest(i,j) < H_cavity(i,j) .and. D_plume(i,j) > H_cavity(i,j)) then
!!                D_plume(i,j) = H_cavity(i,j)
!!             endif
         
             ! If the relaxed D_plume falls below H_cavity (after being above H_cavity on the
             ! previous iteration), then reset it to H_cavity.
!!             if (D_plume_latest(i,j) > H_cavity(i,j) .and. D_plume(i,j) < H_cavity(i,j)) then
!!                D_plume(i,j) = H_cavity(i,j)
!!             endif

  
             eta_plume(i,j) = max(D_plume(i,j) - H_cavity(i,j), 0.d0)
             if (i==itest .and. j==jtest) then
                print*, 'Unrelaxed: i, j, new D_plume, eta_plume =', i, j, D_plume(i,j), eta_plume(i,j)
                print*, 'Unrelaxed change in D =', D_plume_latest(i,j) - D_plume(i,j)
             endif

             ! Prevent the plume from inflating to a large value in one time step,
             !  due to the absence of du/deta terms in the matrix.
             ! This is done by limiting D_plume to H_cavity on this iteration.
             !  The free surface can inflate in the next iteration, when du/deta terms are present.

!             if (D_plume_latest(i,j) < H_cavity(i,j) .and. D_plume(i,j) > H_cavity(i,j)) then
!                D_plume(i,j) = H_cavity(i,j)
!             endif

             if (D_plume_latest(i,j) < H_cavity(i,j)) then
                if (D_plume(i,j) > H_cavity(i,j)) then
                   D_plume(i,j) = H_cavity(i,j)
!!!                else
!!!                   D_plume(i,j) = D_plume_latest(i,j) + relax_D*(D_plume(i,j) - D_plume_latest(i,j))
                endif
             endif

             ! Do a simple relaxation to prevent eta_plume from oscillating back and forth
             ! If relax_eta = 1, there is no adjustment of D_plume.
             ! If relax_eta < 1, we take a weighted average of the previous and current values.
!!             if (D_plume_latest(i,j) > H_cavity(i,j) .and. D_plume(i,j) > H_cavity(i,j)) then
!!             if (D_plume_latest(i,j) > H_cavity(i,j)) then

             if (D_plume_latest(i,j) >= H_cavity(i,j)) then
                D_plume(i,j) = D_plume_latest(i,j) + relax_eta*(D_plume(i,j) - D_plume_latest(i,j))
             endif
             
             ! If the relaxed D_plume falls below H_cavity (after being above H_cavity on the
             ! previous iteration), then reset it to H_cavity.
             ! D_plume can then fall below H_cavity in the next iteration.
!!             if (D_plume_latest(i,j) > H_cavity(i,j) .and. D_plume(i,j) < H_cavity(i,j)) then
!!                D_plume(i,j) = H_cavity(i,j)
!!             endif

             ! This will prevent D_plume from going below H_cavity if relax_eta < 1.
             if (D_plume_latest(i,j) > H_cavity(i,j) .and. &
                  D_plume(i,j) < D_plume_latest(i,j) + relax_eta*(H_cavity(i,j) - D_plume_latest(i,j)) ) then
                D_plume(i,j) = D_plume_latest(i,j) + relax_eta*(H_cavity(i,j) - D_plume_latest(i,j))
             endif

             ! recompute eta_plume
             eta_plume(i,j) = max(D_plume(i,j) - H_cavity(i,j), 0.0d0)
             if (i==itest .and. j==jtest) then
                print*, 'Relaxed: i, j, new D_plume, eta_plume =', i, j, D_plume(i,j), eta_plume(i,j)
             endif

             !WHL - Play with different values of threshold
             if (D_plume_latest(i,j) > H_cavity(i,j) .and. eta_plume(i,j) < 1.0d-8) then
                D_plume(i,j) = H_cavity(i,j)
                eta_plume(i,j) = 0.0d0
             endif

          endif
       enddo
    enddo

    i = itest
    j = jtest
    print*, ' '
    print*, 'H_cavity, relaxed D, eta =', H_cavity(i,j), D_plume(i,j), eta_plume(i,j)
    print*, ' '
    print*, 'After relaxation:'
    print*, ' '
    print*, 'i, j: H_cavity, D_plume_latest, D_plume, eta_latest, eta:', &
         H_cavity(i,j), D_plume_latest(i,j), D_plume(i,j), eta_plume_latest(i,j), eta_plume(i,j)
    if (apply_jacobian) then
       print*, 'predicted du_east(i,j) =', &
              du_deta_west(i,j)*(eta_plume(i,j) - eta_plume_latest(i,j)) &
            + du_deta_east(i,j)*(eta_plume(i+1,j) - eta_plume_latest(i+1,j)) &
            + du_deta_northwest(i,j)*(eta_plume(i,j+1) - eta_plume_latest(i,j+1)) &
            + du_deta_northeast(i,j)*(eta_plume(i+1,j+1) - eta_plume_latest(i+1,j+1)) &
            + du_deta_southwest(i,j)*(eta_plume(i,j-1) - eta_plume_latest(i,j-1)) &
            + du_deta_southeast(i,j)*(eta_plume(i+1,j-1) - eta_plume_latest(i+1,j-1))
       print*, 'predicted du_east(i-1,j) =', &
              du_deta_east(i-1,j)*(eta_plume(i,j) - eta_plume_latest(i,j)) &
            + du_deta_west(i-1,j)*(eta_plume(i-1,j) - eta_plume_latest(i-1,j)) &
            + du_deta_northwest(i-1,j)*(eta_plume(i-1,j+1) - eta_plume_latest(i-1,j+1)) &
            + du_deta_northeast(i-1,j)*(eta_plume(i,j+1) - eta_plume_latest(i,j+1)) &
            + du_deta_southwest(i-1,j)*(eta_plume(i-1,j-1) - eta_plume_latest(i-1,j+1)) &
            + du_deta_southeast(i-1,j)*(eta_plume(i,j-1) - eta_plume_latest(i,j-1))
       print*, 'predicted dv_north(i,j) =', &
              dv_deta_south(i,j)*(eta_plume(i,j) - eta_plume_latest(i,j)) &
            + dv_deta_north(i,j)*(eta_plume(i,j+1) - eta_plume_latest(i,j+1)) &
            + dv_deta_northwest(i,j)*(eta_plume(i-1,j+1) - eta_plume_latest(i-1,j+1)) &
            + dv_deta_northeast(i,j)*(eta_plume(i+1,j+1) - eta_plume_latest(i+1,j+1)) &
            + dv_deta_southwest(i,j)*(eta_plume(i-1,j) - eta_plume_latest(i-1,j)) &
            + dv_deta_southeast(i,j)*(eta_plume(i+1,j) - eta_plume_latest(i+1,j))
       print*, 'predicted dv_south(i,j) =', &
              dv_deta_north(i,j-1)*(eta_plume(i,j) - eta_plume_latest(i,j)) &
            + dv_deta_south(i,j-1)*(eta_plume(i,j-1) - eta_plume_latest(i,j-1)) &
            + dv_deta_northwest(i,j-1)*(eta_plume(i-1,j) - eta_plume_latest(i-1,j)) &
            + dv_deta_northeast(i,j-1)*(eta_plume(i+1,j) - eta_plume_latest(i+1,j)) &
            + dv_deta_southwest(i,j-1)*(eta_plume(i-1,j-1) - eta_plume_latest(i-1,j-1)) &
            + dv_deta_southeast(i,j-1)*(eta_plume(i+1,j-1) - eta_plume_latest(i+1,j-1))
    else
       print*, 'predicted du_east(i,j) =', &
              c_east(i,j)*(eta_plume(i,j) - eta_plume_latest(i,j)) &
            - c_east(i,j)*(eta_plume(i+1,j) - eta_plume_latest(i+1,j))
       print*, 'predicted du_east(i-1,j) =', &
              c_east(i-1,j)*(eta_plume(i-1,j) - eta_plume_latest(i-1,j)) &
            - c_east(i-1,j)*(eta_plume(i,j) - eta_plume_latest(i,j))
    endif

  end subroutine compute_plume_thickness_old

!****************************************************

  subroutine compute_plume_melt_rate(&
       nx,         ny,      &
       gammaT,              &
       gammaS,              &
       plume_mask_cell,     &
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
         plume_mask_cell        ! = 1 for cells where scalar plume variables are computed

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
!       print*, 'Trial cubic solution =', solution
!       print*, 'True solution =', (10.d0 + sqrt(108.d0))**(1.d0/3.d0) - (-10.d0 + sqrt(108.d0))**(1.d0/3.d0) + 5.d0


    ! Loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          
          if (plume_mask_cell(i,j) == 1 .and. entrainment(i,j) > 0.0d0) then
             
             ! Interpolate the plume speed to the cell center, and compute the friction velocity ustar.
             
             u_plume = (u_plume_east(i,j) + u_plume_east(i-1,j)) / 2.0d0
             v_plume = (v_plume_north(i,j) + v_plume_north(i,j-1)) / 2.0d0
             plume_speed = sqrt(u_plume**2 + v_plume**2 + u_tidal**2)
             ustar_plume(i,j) = sqrt(c_drag) * plume_speed

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
                print*, 'pressure (Pa) =', pressure(i,j)
                print*, 'T_factor (m/s/deg), S_factor (m/s)=', T_factor, S_factor
                print*, 'entrainment (m/s) =', entrainment(i,j)
                print*, 'm1 (m/s/psu) =', m1
                print*, 'm2 (m/s) =', m2
                print*, 'denom =', denom
                print*, 'a, b, c, d =', ma, mb, mc, md
                print*, 'residual of cubic solve =', ma*bmlt_float(i,j)**3 + mb*bmlt_float(i,j)**2 + mc*bmlt_float(i,j) + md
             endif
             
             ! Given the melt rate, compute Sb and Tb
!               S_basal(i,j) = (S_factor * entrainment(i,j) * S_ambient(i,j)) /  &
!                               ( (bmlt_float(i,j) + S_factor) * (bmlt_float(i,j) + entrainment(i,j)) )
             S_basal(i,j) = (bmlt_float(i,j) - m2) / m1
             T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)

             ! Given m, compute S and T for plume
             T_plume(i,j) = T_ambient(i,j) - (lhci/(spec_heat_water*entrainment(i,j))) * bmlt_float(i,j)
             S_plume(i,j) = S_ambient(i,j) * entrainment(i,j) / (bmlt_float(i,j) + entrainment(i,j))

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

          else    ! plume_mask_cell = 0
             
             bmlt_float(i,j) = 0.0d0
             
             S_plume(i,j) = S_ambient(i,j)
             T_plume(i,j) = T_ambient(i,j)
             
             S_basal(i,j) = S_ambient(i,j)
             T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
             
          endif   ! plume_mask_cell and entrainment > 0
          
       enddo   ! i
    enddo   ! j

  end subroutine compute_plume_melt_rate

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
