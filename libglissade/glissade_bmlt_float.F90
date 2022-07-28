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
  use glimmer_paramets, only: unphys_val
  use glimmer_log
  use glide_types
  use cism_parallel, only: this_rank, main_task, nhalo, &
       parallel_type, parallel_halo, parallel_globalindex, parallel_boundary_value, &
       parallel_reduce_sum, parallel_reduce_min, parallel_reduce_max

  implicit none
  
  private
  public :: verbose_bmlt_float, glissade_basal_melting_float, &
       glissade_bmlt_float_thermal_forcing_init, glissade_bmlt_float_thermal_forcing

!!    logical :: verbose_bmlt_float = .false.
    logical :: verbose_bmlt_float = .true.

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

    ! loop limits for debug diagnostics
    integer :: kmin_diag = 1
    integer :: kmax_diag = 1

  contains

!****************************************************

  subroutine glissade_basal_melting_float(whichbmlt_float,              &
                                          parallel,                     &
                                          ewn,         nsn,             &
                                          dew,         dns,             &
                                          itest,       jtest,    rtest, &
                                          x1,                           &
                                          thck,        lsrf,            &
                                          topg,        eus,             &
                                          basal_melt,  ocean_data)

    use glissade_masks, only: glissade_get_masks
    use glimmer_paramets, only: tim0, thk0

    ! Compute the rate of basal melting for floating ice by one of several methods.

    !-----------------------------------------------------------------
    ! Input/output arguments
    !-----------------------------------------------------------------

    integer, intent(in) :: whichbmlt_float            ! method for computing melt rate of floating ice

    type(parallel_type), intent(in) :: &
         parallel                ! info for parallel communication

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

    type(glide_ocean_data), intent(inout) :: &
         ocean_data              ! derived type with fields and parameters related to ocean data

    !-----------------------------------------------------------------
    ! Note: The basal_melt derived type includes the 2D output field bmlt_float,
    !        along with a number of prescribed parameters for MISMIP+:
    !
    !       MISMIP+ Ice1
    !       - bmlt_float_omega     ! lapse rate for basal melting (m/s/m = s-1), default = 0.2/yr
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
    ! Local variables and pointers set to components of basal_melt and plume derived type
    !----------------------------------------------------------------      

    real(dp), dimension(:,:), pointer :: &
         bmlt_float             ! basal melt rate for floating ice (m/s) (> 0 for melt, < 0 for freeze-on)

    real(dp), dimension(:,:), pointer :: &
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

    integer  :: i, j
    real(dp) :: h_cavity        ! depth of ice cavity beneath floating ice (m)
    real(dp) :: z_draft         ! draft of floating ice (m below sea level)

    real(dp) :: frz_ramp_factor    ! multiplying factor for linear ramp at depths with basal freezing
    real(dp) :: melt_ramp_factor   ! multiplying factor for linear ramp at depths with basal melting

    !-----------------------------------------------------------------
    ! Compute the basal melt rate for floating ice
    !-----------------------------------------------------------------

    if (main_task .and. verbose_bmlt_float) print*, 'Computing bmlt_float, whichbmlt_float =', whichbmlt_float

    ! Set bmlt_float pointer and initialize
    bmlt_float  => basal_melt%bmlt_float
    bmlt_float(:,:) = 0.0d0

    ! Compute masks:
    ! - ice_mask = 1 where thck > 0
    ! - floating_mask = 1 where thck > 0 and ice is floating;
    ! - ocean_mask = 1 where topg is below sea level and ice is absent
    !Note: The '0.0d0' argument is thklim. Here, any ice with thck > 0 gets ice_mask = 1.
    !TODO: Modify glissade_get_masks so that 'parallel' is not needed.

    call glissade_get_masks(ewn,           nsn,            &
                            parallel,                      &
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
!                if (j == jtest .and. verbose_bmlt_float) then
!                   print*, 'cavity, tanh, thck, draft, melt rate (m/yr):', i, j, h_cavity, &
!                         tanh(h_cavity/basal_melt%bmlt_float_h0), thck(i,j), z_draft, bmlt_float(i,j)*scyr
!                endif

             endif   ! ice is present and floating

          enddo
       enddo

       !WHL - debug
       if (verbose_bmlt_float .and. this_rank == rtest) then
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

                ! Note: bmlt_float will be reduced further in shallow cavities if basal_melt%bmlt_cavity_h0 > 0.
                !       This reduction is done in subroutine glissade_bmlt_float_solve.
             endif   ! ice is present and floating

          enddo
       enddo

       !debug
       if (verbose_bmlt_float .and. this_rank == rtest) then
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

       ! TODO: Develop a new plume model.  I removed the old one, leaving just some utility subroutines.
       ! This option is not supported; the code aborts at startup if the user selects it.

       ! Compute melt rates using a plume model, given vertical profiles of T and S in the ambient ocean
       !
       ! See this paper for details:
       ! X. S. Asay-Davis et al. (2016), Experimental design for three interrelated 
       !    marine ice sheet and ocean model intercomparison projects: 
       !    MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
       !    Geosci. Model Devel., 9, 2471-2497, doi: 10.5194/gmd-9-2471-2016.

       ! Assign local pointers and variables to components of the plume derived type

!       T0 = plume%T0
!       Tbot = plume%Tbot
!       S0 = plume%S0
!       Sbot = plume%Sbot
!       zbed_deep = plume%zbed_deep
!       gammaT = plume%gammaT
!       gammaS = plume%gammaS

       ! the following fields are used or computed by the plume model
!       T_ambient     => plume%T_ambient
!       S_ambient     => plume%S_ambient

       ! Given the ice draft in each floating grid cell, compute the ambient ocean T and S
       !  using the prescribed MISOMIP profile.

!       where (floating_mask == 1) 

          ! MISOMIP+ profiles, Eqs. 21 and 22 
!          T_ambient = T0 + (Tbot - T0) * (lsrf / zbed_deep)
!          S_ambient = S0 + (Sbot - S0) * (lsrf / zbed_deep)

!       elsewhere
!          T_ambient = T0
!          S_ambient = S0
!       endwhere

       ! Note: The plume model expects floating_mask, T_ambient and S_ambient to be correct in halo cells.
       !       This is likely the case already, but do halo updates just in case.
!       call parallel_halo(floating_mask, parallel)
!       call parallel_halo(T_ambient, parallel)
!       call parallel_halo(S_ambient, parallel)

    endif   ! whichbmlt_float

  end subroutine glissade_basal_melting_float

!****************************************************

  subroutine glissade_bmlt_float_thermal_forcing_init(model, ocean_data)

    use glimmer_paramets, only: thk0, len0, tim0, unphys_val
    use glissade_masks, only : glissade_get_masks

    ! Initialization for basal melting based on ocean thermal forcing

    type(glide_global_type), intent(inout) :: model   !> derived type holding ice-sheet info

    type(glide_ocean_data), intent(inout) ::  &
         ocean_data                                   !> derived type holding ocean input data
    
    integer, dimension(model%general%ewn, model%general%nsn) ::  &
         ice_mask,                   &  ! = 1 if ice is present (thck > 0)
         ocean_mask,                 &  ! = 1 if topg < 0 and ice is absent
         land_mask                      ! = 1 if topg >= 0

    real(dp), dimension(:), allocatable :: &
         deltaT_basin_ismip6            ! prescribed deltaT_basin values for each of 18 basins

    integer :: itest, jtest, rtest      ! coordinates of diagnostic point
    integer :: i, j, k, nb
    integer :: ewn, nsn

    integer :: basin_number_min         ! global minval of the basin_number field

    logical :: simple_init = .false.

    type(parallel_type) :: parallel   ! info for parallel communication

    parallel = model%parallel

    ewn = model%general%ewn
    nsn = model%general%nsn

    ! Set debug diagnostics
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    allocate(deltaT_basin_ismip6(ocean_data%nbasin))

    if (verbose_bmlt_float .and. main_task) then
       print*, 'In glissade_bmlt_float_thermal_forcing_init'
       print*, 'bmlt_float_thermal_forcing_param =', model%options%bmlt_float_thermal_forcing_param
    endif

    !WHL - debug - some simple initializations for testing
    ! In config file, set nbasin = 4 (for 4 cores on Mac) and nzocn = 10

    if (simple_init) then

       ! Assign basin numbers based on this_rank (0 to 3 on a Mac)
       ocean_data%basin_number(:,:) = this_rank

       ! Set deltaT_basin in a similar way based on this_rank
       ! Will have more melting with larger rank
       ocean_data%deltaT_ocn(:,:) = 0.50d0 * this_rank

       ! Use Xylar's median value (m/yr) for gamma0
       ocean_data%gamma0 = 15000.d0

       ! Let the transient thermal forcing be steady in time, increasing from surface to bed
       do k = 1, ocean_data%nzocn
          ocean_data%zocn(k) = -100.0d0 * k   ! ocean level every 100 m
          ocean_data%thermal_forcing(k,:,:) = -ocean_data%zocn(k) / 500.0d0  ! 2 K/km
       enddo

       return

    endif  ! simple_init

    if (model%options%is_restart == RESTART_FALSE) then

       if (model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL .or.  &
           model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or. &
           model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

          ! Initialize deltaT_ocn based on deltaT_basin_ismip6, if needed for the ISMIP6 option
          ! For other options, deltaT_ocn(:,:) = 0 initially or has already been read in

          if (model%options%which_ho_bmlt_basin == HO_BMLT_BASIN_ISMIP6) then

             if (main_task) then
                print*, 'Assign deltaT_basin from ismip6'
             endif

             ! Note: For now, these values are hardwired for the standard 16 ISMIP6 basins
             if (ocean_data%nbasin /= 16) then
                call write_log('Error, ISMIP6 deltaT_basin values are set for exactly 16 Antarctic basins', GM_FATAL)
             endif

             ! Set values computed by Nico Jourdain to match observed basin-scale mean melt.
             ! See Jourdain et al. (2019) and Lipscomb et al. (2021).
             ! Note: Uncommented values are for the MeanAnt calibration; commented values are for PIGL.
             !       It would be possible to make MeanAnt v. PIGL a config option,
             !        but for now a new compile is needed to use the PIGL numbers.
             if (model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL) then

                ! MeanAnt
                deltaT_basin_ismip6 = (/ 0.68,  0.15,  0.62,  0.87,  0.36,  0.05, -0.11,  0.51,  &
                                         1.28, -0.13, -0.95, -0.13, -0.17, -0.05,  0.12, -0.34 /)

                ! PIGL
!                deltaT_basin_ismip6 = (/-0.04, -0.24,  0.06, -0.13, -0.17, -0.56, -0.27, -0.34, &
!                                        -0.14, -1.17, -2.01, -0.74, -0.38, -0.27, -0.11, -1.04 /)

             elseif (model%options%bmlt_float_thermal_forcing_param== BMLT_FLOAT_TF_ISMIP6_NONLOCAL) then

                ! MeanAnt
                deltaT_basin_ismip6 = (/ 0.57,  0.13,  0.51,  0.70,  0.27,  0.08, -0.12,  0.43,  &
                                         1.07, -0.01, -0.66, -0.06, -0.12, -0.06,  0.10, -0.16 /)

                ! PIGL
!                deltaT_basin_ismip6 = (/-0.19, -0.22, -0.10, -0.39, -0.30, -0.39, -0.28, -0.39,  &
!                                        -0.43, -0.70, -1.43, -0.37, -0.27, -0.27, -0.12, -0.46 /)

             elseif (model%options%bmlt_float_thermal_forcing_param== BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

                if (main_task) print*, '   Assign nonlocal-slope values'

                ! MeanAnt
                deltaT_basin_ismip6 = (/ 0.36, -0.03,  0.45,  0.05,  0.02, -0.22, -0.01,  0.37,  &
                                         0.64, -0.03, -0.58, -0.10, -0.11, -0.01,  0.14, -0.15 /)

                ! PIGL
!                deltaT_basin_ismip6 = (/ 0.03, -0.18,  0.13, -0.31, -0.21, -0.37, -0.14, -0.05,  &
!                                        -0.03, -0.42, -1.02, -0.27, -0.19, -0.15,  0.00, -0.32 /)

             endif

             ! Assign the numbers above to each grid cell, given its basin number
             do j = 1, nsn
                do i = 1, ewn
                   nb = ocean_data%basin_number(i,j)
                   if (nb >= 1) then
                      ocean_data%deltaT_ocn(i,j) = deltaT_basin_ismip6(nb)
                   endif
                enddo
             enddo

          endif   ! ho_bmlt_basin_ismip6

          !WHL - In earlier code, nonzero values of gamma0 could either be set in the config file,
          !       read from the input file, or assigned here based on the ISMIP6 parameterization.
          !      This led to errors because with multiple ways of setting gamma0, it was unclear
          !       which value would actually be used.
          !      Now, nonzero values of gamma0 must be set in the config file.

          if (verbose_bmlt_float .and. this_rank==rtest) then
             print*, ' '
             print*, 'Initialize ISMIP6 sub-shelf melting'
             print*, ' '
             print*, 'k, zocn:'
             do k = 1, ocean_data%nzocn
                print*, k, ocean_data%zocn(k)
             enddo
             print*, ' '
             print*, 'gamma0 =', ocean_data%gamma0
             print*, ' '
             print*, 'basin_number, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(i10)',advance='no') ocean_data%basin_number(i,j)
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'deltaT_ocn'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.4)',advance='no') ocean_data%deltaT_ocn(i,j)
                enddo
                write(6,*) ' '
             enddo
             do k = kmin_diag, kmax_diag
                print*, ' '
                print*, 'thermal_forcing, k =', k
                do j = jtest+3, jtest-3, -1
                   write(6,'(i6)',advance='no') j
                   do i = itest-3, itest+3
                      write(6,'(f10.3)',advance='no') ocean_data%thermal_forcing(k,i,j)
                   enddo
                   write(6,*) ' '
                enddo
             enddo
          endif  ! verbose_bmlt_float

          ! Fill halos (might not be needed)
          ! TODO: Remove these halo updates?
          call parallel_halo(ocean_data%basin_number, parallel)
          call parallel_halo(ocean_data%deltaT_ocn, parallel)
          call parallel_halo(ocean_data%thermal_forcing, parallel)

          ! Make sure every cell is assigned a basin number >= 1.
          ! If not, then extrapolate the current basin numbers to fill the grid.
          ! Note: Could remove this code if guaranteed that the basin number in the input file
          !       has valid values everywhere.
          basin_number_min = minval(model%ocean_data%basin_number)
          basin_number_min = parallel_reduce_min(basin_number_min)

          if (basin_number_min < 1) then

             if (verbose_bmlt_float .and. main_task) then
                print*, 'Extrapolate basin numbers'
             endif

             call basin_number_extrapolate(&
                  ewn,             nsn,     &
                  parallel,                 &
                  model%ocean_data%nbasin,  &
                  model%ocean_data%basin_number)
          endif

       endif   ! ISMIP6 thermal forcing option
 
    endif  ! restart_false

  end subroutine glissade_bmlt_float_thermal_forcing_init

!****************************************************

  subroutine glissade_bmlt_float_thermal_forcing(&
       bmlt_float_thermal_forcing_param, &
       ocean_data_extrapolate,    &
       parallel,                  &
       nx,        ny,             &
       dew,       dns,            &
       itest,     jtest,   rtest, &
       ice_mask,                  &
       ocean_mask,                &
       marine_connection_mask,    &
       f_ground_cell,             &
       thck,                      &
       lsrf,                      &
       topg,                      &
       ocean_data,                &
       bmlt_float,                &
       tf_anomaly_in,             &
       tf_anomaly_basin_in)

    use glimmer_paramets, only: thk0, unphys_val
    use glissade_grid_operators, only: glissade_slope_angle
    use glissade_utils, only: glissade_basin_average

    ! Compute a 2D field of sub-ice-shelf melting given a 3D thermal forcing field
    !  and the current lower ice surface, using either a local or nonlocal melt parameterization.

    integer, intent(in) :: &
         bmlt_float_thermal_forcing_param, & !> melting parameterization used to derive melt rate from thermal forcing;
                                             !> current options are quadratic and ISMIP6 local, nonlocal and nonlocal_slope
         ocean_data_extrapolate              !> = 1 if TF is to be extrapolated to sub-shelf cavities, else = 0

    type(parallel_type), intent(in) :: &
         parallel                            !> info for parallel communication

    integer, intent(in) :: &
         nx, ny                              !> number of grid cells in each direction

    real(dp), intent(in) :: &
         dew, dns                            !> grid cell size (m)

    integer, intent(in) :: itest, jtest, rtest  !> coordinates of diagnostic point

    integer, dimension(nx,ny), intent(in) :: &
         ice_mask,               & !> = 1 where ice is present (H > 0) else = 0
         marine_connection_mask    !> = 1 for cells with a marine connection to the ocean, else = 0
                                   !> Note: marine_connection_mask includes paths through grounded marine-based cells

    integer, dimension(nx,ny), intent(inout) :: &
         ocean_mask                !> = 1 for ice-free ocean, else = 0;
                                   !> can be set to 0 below where there is no valid ocean data

    real(dp), dimension(nx,ny), intent(in) ::  &
         f_ground_cell,          & !> fraction of grounded ice in each cell
         thck,                   & !> ice thickness (m)
         lsrf,                   & !> ice lower surface elevation (m), negative below sea level
         topg                      !> bed topography (m), negative below sea level

    type(glide_ocean_data), intent(inout) :: &
         ocean_data         !> derived type with fields and parameters related to ocean thermal forcing;
                            !> includes the following fields:
                            !> thermal_forcing = input 3D thermal forcing (deg C) from ocean data;
                            !>    can be modified by extrapolation
                            !> thermal_forcing_lsrf = 2D TF applied at lower ice surface
                            !> data_mask = mask that indicates where ocean data is valid (and need not be extrapolated)
                            !> nzocn = number of ocean levels
                            !> zocn = ocean levels (m)
                            !> nbasin = number of ocean basins
                            !> basin_number = integer assigned to each basin
                            !> gamma0 = basal melt rate coefficient for ISMIP6 melt parameterization
                            !> deltaT_ocn = ocean temperature corrections for ISMIP6 melt parameterization

    real(dp), dimension(:,:), intent(out) :: &
         bmlt_float         !> basal melt rate for floating ice (m/s)

    real(dp), intent(in), optional ::  &
         tf_anomaly_in      !> uniform thermal forcing anomaly (deg C), applied everywhere

    integer, intent(in), optional :: &
         tf_anomaly_basin_in  !> basin where anomaly is applied; for default value of 0, apply to all basins

    ! local variables

    integer :: i, j, k, nb
    integer :: iglobal, jglobal

    character(len=256) :: message

    ! Note: thermal_forcing_mask is similar but not identical to floating_mask.
    !       * floating_mask = 1 where ice is present, and thck satisfies a flotation condition
    !       * thermal_forcing_mask = 1 where ice is present (thck > 0) and f_ground_cell < 1, with lakes excluded

    integer, dimension(nx,ny) ::  &
         thermal_forcing_mask,          & ! = 1 where thermal forcing and bmlt_float can be nonzero, else = 0
         new_mask                         ! temporary mask

    real(dp), dimension(ocean_data%nzocn,nx,ny) ::  &
         thermal_forcing_in               ! TF passed to subroutine interpolate_thermal_forcing_to_lsrf;
                                          ! optionally corrected for nonzero tf_anomaly
    real(dp), dimension(nx,ny) ::  &
         theta_slope,                   & ! sub-shelf slope angle (radians)
         f_float                          ! weighting function for computing basin averages, in range [0,1]

    ! Note: Ocean basins are indexed from 1 to nbasin (previously indexed from 0 to nbasin-1)
    real(dp), dimension(ocean_data%nbasin) :: &
         thermal_forcing_basin,        &  ! basin average thermal forcing (K) at current time
         deltaT_basin_avg                 ! basin average value of deltaT_ocn

    real(dp) :: &
         tf_anomaly                       ! local version of tf_anomaly_in

    integer ::  &
         tf_anomaly_basin                 ! local version of tf_anomaly_basin_in

    ! Note: This range ought to cover all regions where ice is present, but could be modified if desired.
    real(dp), parameter ::  &
         thermal_forcing_max = 20.d0,  &  ! max allowed value of thermal forcing (K)
         thermal_forcing_min = -5.d0      ! min allowed value of thermal forcing (K)

    !TODO - Make H0_float a config parameter?
    real(dp), parameter ::  &
         H0_float = 50.d0                 ! thickness scale (m) for floating ice; used to reduce weights when H < H0_float

    if (verbose_bmlt_float .and. main_task) then
       print*, ' '
       print*, 'In subroutine glissade_bmlt_float_thermal_forcing'
       print*, '   bmlt_float_thermal_forcing_param =', bmlt_float_thermal_forcing_param
       print*, '   ocean_data_extrapolate =', ocean_data_extrapolate
       print*, '   nbasin =', ocean_data%nbasin
    endif

    if (present(tf_anomaly_in)) then
       tf_anomaly = tf_anomaly_in
    else
       tf_anomaly = 0.0d0
    endif

    if (present(tf_anomaly_basin_in)) then
       tf_anomaly_basin = tf_anomaly_basin_in
    else
       tf_anomaly_basin = 0
    endif

    ! initialize the output
    bmlt_float = 0.0d0

    ! Make sure thermal_forcing is up to date in halo cells.
    call parallel_halo(ocean_data%thermal_forcing, parallel)

    ! Insert an unphysical value at the global boundary.
    ! This is done to handle the case that global_bc = no_ice,
    !  which puts zeroes in global boundary cells.
    ! We do not want these zeroes to be interpreted as realistic thermal_forcing values.
    call parallel_boundary_value(ocean_data%thermal_forcing, unphys_val, parallel)

    ! Set thermal_forcing_mask
    ! This mask identifies cells where we could have basal melting and need valid TF data.
    where (ice_mask == 1 .and. f_ground_cell < 1.0d0 .and. marine_connection_mask == 1)
       thermal_forcing_mask = 1
    elsewhere
       thermal_forcing_mask = 0
    endwhere

    ! Extend the mask to include neighboring ice-free ocean cells.
    ! These cells could be ice-filled (and subject to melting) after transport.

    new_mask = thermal_forcing_mask

    do j = 2, ny-1
       do i = 2, nx-1
          if (thermal_forcing_mask(i,j) == 0 .and. ocean_mask(i,j) == 1) then
             if (thermal_forcing_mask(i-1,j) == 1 .or. thermal_forcing_mask(i+1,j) == 1 .or. &
                 thermal_forcing_mask(i,j-1) == 1 .or. thermal_forcing_mask(i,j+1) == 1) then
                new_mask(i,j) = 1
             endif
          endif
       enddo
    enddo

    thermal_forcing_mask = new_mask
    call parallel_halo(thermal_forcing_mask, parallel)

    !-----------------------------------------------
    ! Optionally, extrapolate the ocean forcing data to ice shelf cavities.
    ! Note: For coupled modeling (OCEAN_DATA_GLAD), thermal forcing is received once per mass balance time step.
    !       Typically, it is computed only on the ocean grid and needs to be extrapolated to cavities.
    !       During the first ice sheet dynamics time step within a mass balance time step,
    !        several tens of iterations typically are needed to extrapolate open-ocean values
    !        into large shelf cavities.
    !       On subsequent steps, only a few iterations are needed to extrapolate into newly floating cells.
    !-----------------------------------------------

    if (ocean_data_extrapolate == OCEAN_DATA_EXTRAPOLATE_TRUE) then

       if (verbose_bmlt_float .and. this_rank == rtest) then
          print*, ' '
          print*, 'Extrapolating TF data beneath ice shelves: rank, i, j =', rtest, itest, jtest
          print*, ' '
          print*, 'thermal_forcing_mask:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') thermal_forcing_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'marine_connection_mask:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') marine_connection_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'ocean mask:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') ocean_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'lsrf:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') lsrf(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'topg:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') topg(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'TF before extrapolating:'
          do k = kmin_diag, kmax_diag
             print*, ' '
             print*, 'kocn =', k
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') ocean_data%thermal_forcing(k,i,j)
                enddo
                write(6,*) ' '
             enddo
          enddo
       endif  ! verbose_bmlt_float

       ! Extrapolate the 3D thermal forcing field to sub-shelf cavities.
       ! Note: For now, unfilled cells retain values of unphys_val.
       ! TODO: Replace unphys_val with an integer mask?
 
       call glissade_thermal_forcing_extrapolate(&
            nx,        ny,                     &
            parallel,                          &
            itest,     jtest,     rtest,       &
            ocean_data%nzocn,                  &
            ocean_data%zocn,                   &  ! m
            lsrf,                              &  ! m
            topg,                              &  ! m
            thermal_forcing_mask,              &
            marine_connection_mask,            &
            unphys_val,                        &  ! identifies unfilled cells on input
            unphys_val,                        &  ! default value given to unfilled cells on output
            ocean_data%thermal_forcing)

       if (verbose_bmlt_float .and. this_rank == rtest) then
          print*, ' '
          print*, 'TF after extrapolating, rank, i, j =', rtest, itest, jtest
          do k = kmin_diag, kmax_diag
             print*, ' '
             print*, 'kocn =', k
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') ocean_data%thermal_forcing(k,i,j)
                enddo
                write(6,*) ' '
             enddo
          enddo
       endif

    else    ! ocean data are already given everywhere; do not extrapolate

       if (verbose_bmlt_float .and. this_rank == rtest) then
          print*, ' '
          print*, 'TF to interpolate to lsrf, rank, i, j =', rtest, itest, jtest
          do k = kmin_diag, kmax_diag
             print*, ' '
             print*, 'kocn =', k
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') ocean_data%thermal_forcing(k,i,j)
                enddo
                write(6,*) ' '
             enddo
          enddo
       endif

    endif   ! ocean_data_extrapolate

    !-----------------------------------------------
    ! Interpolate the thermal forcing to the lower ice surface.
    !-----------------------------------------------

    ! Optionally, add a uniform anomaly (= 0 by default) to the thermal forcing.
    thermal_forcing_in = ocean_data%thermal_forcing
    if (tf_anomaly /= 0.0d0) then
       if (tf_anomaly_basin >= 1 .and. tf_anomaly_basin <= ocean_data%nbasin) then  ! apply to one basin
          do j = 1, ny
             do i = 1, nx
                if (ocean_data%basin_number(i,j) == tf_anomaly_basin) then
                   thermal_forcing_in(:,i,j) = thermal_forcing_in(:,i,j) + tf_anomaly
                endif
             enddo
          enddo
       else   ! apply to all basins
          thermal_forcing_in = thermal_forcing_in + tf_anomaly
       endif
    endif

    call interpolate_thermal_forcing_to_lsrf(&
         nx,                ny,              &
         ocean_data%nzocn,                   &
         ocean_data%zocn,                    &
         thermal_forcing_mask,               &
         lsrf,                               &
         thermal_forcing_in,                 &
         ocean_data%thermal_forcing_lsrf)

    ! Bug check: Make sure there are no extreme values of thermal forcing.
    !            This could happen, for example, if the input thermal forcing has special values
    !             that are not overwritten with realistic values via extrapolation.

    do j = 1, ny
       do i = 1, nx
          if (ocean_data%thermal_forcing_lsrf(i,j) > thermal_forcing_max) then
             call parallel_globalindex(i, j, iglobal, jglobal, parallel)
             write(message,*) &
                  'Ocean thermal forcing error: extreme TF at rank, i, j, iglobal, jglobal, lsrf, TF =', &
                  this_rank, i, j, iglobal, jglobal, lsrf(i,j), ocean_data%thermal_forcing_lsrf(i,j)
             call write_log(message, GM_FATAL)
          elseif (ocean_data%thermal_forcing_lsrf(i,j) < thermal_forcing_min) then
             call parallel_globalindex(i, j, iglobal, jglobal, parallel)
             write(message,*) &
                  'Ocean thermal forcing error: extreme TF at rank, i, j, iglobal, jglobal, lsrf, TF =', &
                  this_rank, i, j, iglobal, jglobal, lsrf(i,j), ocean_data%thermal_forcing_lsrf(i,j)
             call write_log(message, GM_FATAL)
          endif
       enddo
    enddo

    if (verbose_bmlt_float .and. this_rank==rtest) then
       print*, ' '
       print*, 'basin number =', ocean_data%basin_number(itest,jtest)
       print*, ' '
       print*, 'lsrf (m)'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') lsrf(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thermal_forcing_lsrf (deg C)'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') ocean_data%thermal_forcing_lsrf(i,j)
          enddo
          write(6,*) ' '
       enddo
       if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL .or.  &
           bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or.  &
           bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then
          print*, ' '
          print*, 'deltaT_ocn (deg C)'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') ocean_data%deltaT_ocn(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'corrected TF (deg C)'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') &
                     ocean_data%thermal_forcing_lsrf(i,j) + ocean_data%deltaT_ocn(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif
    endif

    ! For ISMIP6 nonlocal parameterizations, compute the average thermal forcing for the basin.

    if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or.  &
        bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

       ! nonlocal parameterization
       ! Melt rate is a quadratic function of the local thermal forcing and basin-average thermal forcing

       ! Compute a weighting function that is proportional to the floating fraction of ice-filled cells,
       !  and also tapers linearly to zero for thin floating ice.
       ! This function is used to ensure smooth changes in the basin averages as cells
       !  transition between grounded and floating, or between ice-free and thick.

       f_float = 1.0d0 - f_ground_cell

       if (H0_float > 0.0d0) then
          where (thck > 0.0d0)
             f_float = f_float * min(thck/H0_float, 1.0d0)
          elsewhere
             f_float = 0.0d0
          endwhere
       endif

       if (verbose_bmlt_float .and. this_rank == rtest) then
          print*, ' '
          print*, 'thermal_forcing_mask:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') thermal_forcing_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'f_float'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') f_float(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thermal_forcing_mask * f_float:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thermal_forcing_mask(i,j) * f_float(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       ! Compute the average thermal forcing for each basin.
       ! The average is taken over grid cells with thermal_forcing_mask = 1,
       !  with reduced weights for partly grounded cells and thin floating cells.

       call glissade_basin_average(&
            nx,        ny,                   &
            ocean_data%nbasin,               &
            ocean_data%basin_number,         &
            thermal_forcing_mask * f_float,  &
            ocean_data%thermal_forcing_lsrf, &
            thermal_forcing_basin,           &
            itest, jtest, rtest)

       ! For diagnostics, compute the average value of deltaT_ocn in each basin.
       ! Note: Each cell in the basin should have this average value.

       call glissade_basin_average(&
            nx,        ny,                   &
            ocean_data%nbasin,               &
            ocean_data%basin_number,         &
            thermal_forcing_mask * f_float,  &
            ocean_data%deltaT_ocn,           &
            deltaT_basin_avg)

       if (verbose_bmlt_float .and. this_rank==rtest) then
          print*, ' '
          print*, 'thermal_forcing_basin:'
          do nb = 1, ocean_data%nbasin
             print*, nb, thermal_forcing_basin(nb)
          enddo
          print*, ' '
          print*, 'deltaT_basin_avg:'
          do nb = 1, ocean_data%nbasin
             print*, nb, deltaT_basin_avg(nb)
          enddo
       endif

    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL) then

       thermal_forcing_basin = 0.0d0
       deltaT_basin_avg = 0.0d0

    endif

    !-----------------------------------------------
    ! Compute the basal melt rate for each grid cell.
    ! Note: The output bmlt_float has units of m/yr.
    !-----------------------------------------------

    if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_QUADRATIC) then

       if (verbose_bmlt_float .and. main_task) then
          print*, 'Compute basal melt rate from quadratic parameterization'
          print*, 'rank, i, j:', rtest, itest, jtest
       endif

       call quadratic_bmlt_float(&
            nx,                ny,           &
            ocean_data%thermal_forcing_lsrf, &
            thermal_forcing_mask,            &
            bmlt_float)

    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL .or. &
            bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or. &
            bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

       if (verbose_bmlt_float .and. this_rank == rtest) then
          print*, 'Compute basal melt rate from ISMIP6 thermal forcing'
          print*, 'rank, i, j, basin number:', rtest, itest, jtest, ocean_data%basin_number(itest,jtest)
       endif

       ! Compute the basal melt rate based on an ISMIP6 thermal forcing parameterization.
       ! Note: bmlt_float is nonzero only for cells with thermal_forcing_mask = 1.

       call ismip6_bmlt_float(&
            bmlt_float_thermal_forcing_param, &
            nx,                ny,            &
            itest,   jtest,    rtest,         &
            ocean_data%nbasin,                &
            ocean_data%basin_number,          &
            ocean_data%gamma0,                &
            ocean_data%thermal_forcing_lsrf,  &
            ocean_data%deltaT_ocn,            &
            thermal_forcing_basin,            &
            deltaT_basin_avg,                 &
            thermal_forcing_mask,             &
            bmlt_float)

       if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

          ! Compute the angle between the lower ice shelf surface and the horizontal.
          ! This option can be used to concentrate basal melting near the grounding line,
          !  where slopes are typically larger, and to reduce melting near the calving front
          !  where slopes are small.

          call glissade_slope_angle(nx,       ny,     &
                                    dew,      dns,    &  ! m
                                    lsrf,             &  ! m
                                    theta_slope,      &  ! radians
                                    slope_mask_in = ice_mask)

          call parallel_halo(theta_slope, parallel)

          if (verbose_bmlt_float .and. this_rank==rtest) then
             print*, ' '
             print*, 'sin(theta_slope)'
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.5)',advance='no') sin(theta_slope(i,j))
                   enddo
                write(6,*) ' '
             enddo
          endif

          ! Make the melt rate proportional to sin(theta_slope)
          bmlt_float = bmlt_float * sin(theta_slope)

       endif

    endif   ! bmlt_float_thermal_forcing_param

    if (verbose_bmlt_float .and. this_rank==rtest) then
       print*, ' '
       print*, 'bmlt_float (m/yr), before thin ice adjustment'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') bmlt_float(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Reduce the melt rate in cells with thin floating ice,
    !  to reflect that these cells are only partly ice-filled.
    ! Note: This code gives bmlt_float = 0 in ice-free ocean cells,
    !       giving ice a chance to accumulate.
    if (H0_float > 0.0d0) then
       where (f_ground_cell < 1.0d0)
          bmlt_float = bmlt_float * min(thck/H0_float, 1.0d0)
       elsewhere
          bmlt_float = 0.0d0
       endwhere
    endif

    if (verbose_bmlt_float .and. this_rank==rtest) then
       print*, ' '
       print*, 'bmlt_float (m/yr), end of glissade_bmlt_float_thermal_forcing:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') bmlt_float(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Convert from m/yr to m/s for output.
    bmlt_float(:,:) = bmlt_float(:,:) / scyr

  end subroutine glissade_bmlt_float_thermal_forcing

!****************************************************

  subroutine glissade_thermal_forcing_extrapolate(&
       nx,              ny,                    &
       parallel,                               &
       itest,           jtest,        rtest,   &
       nzocn,           zocn,                  &
       lsrf,            topg,                  &
       thermal_forcing_mask,                   &
       marine_connection_mask,                 &
       unphys_val,      default_val,           &
       thermal_forcing)

    use glimmer_physcon, only : dtocnfrz_dz  ! value from Beckmann & Goosse (2003), eq. 2

    ! Extrapolate ocean thermal forcing from an ocean domain that typically excludes
    !  ice shelf cavities, to an expanded domain that includes cavities.
    ! The algorithm works as follows:
    ! First, for each grid cell, determine the levels where we would like to have filled values.
    ! - In ice-free ocean cells, this includes the sea surface down to the bed topography.
    ! - In sub-shelf cavities, this includes the levels ktop:kbot that extend from just above lsrf
    !   to just below topg.
    ! Then loop through each grid cell.  When we get to a cell with levels that are unfilled
    !  and are targeted for filling, we take the average value from filled neighbor cells at the same level.
    !  In this way we extend the domain of filled cells by one set of neighbors per iteration.
    !  We then extrapolate vertically to fill levels in the range ktop:kbot that are not yet filled.
    ! The iteration ends when there are no more cells and levels to fill.
    ! Note: The input thermal_forcing should be set to unphys_val for cells and levels without valid data.
    ! Note: Cells are filled only if they have a connection to the ocean through floating cells.
    !       Interior lakes are not filled.
    ! The algorithm is similar to that of subroutine glissade_scalar_extrapolate but more complex;
    !  there are multiple levels, and data can be extrapolated from a given level to levels above or below.

    integer, intent(in) ::  &
         nx, ny,               & ! grid dimensions
         itest, jtest, rtest,  & ! coordinates of diagnostic point
         nzocn                   ! number of ocean levels

    type(parallel_type), intent(in) :: &
         parallel                ! info for parallel communication

    real(dp), dimension(nzocn), intent(in) :: &
         zocn                    ! ocean levels (m, negative below sea level)

    !TODO - Pass in eus as well as topg?
    real(dp), dimension(nx,ny), intent(in) ::  &
         lsrf,                 & ! lower ice surface elevation (m)
         topg                    ! bed elevation (m)

    integer, dimension(nx,ny), intent(in) ::  &
         thermal_forcing_mask, & ! = 1 where thermal forcing and bmlt_float are potentially nonzero, else = 0
         marine_connection_mask  ! = 1 for cells with marine connection to the ocean, else = 0
                                 ! Note: marine_connection_mask includes paths through grounded marine-based cells

   real(dp), intent(in) :: &
         unphys_val,           & ! unphysical value given to cells/levels not yet filled
         default_val             ! default value given to unfilled cells on output;
                                 ! might be different from unphys_val

    ! The thermal forcing field should be indexed with k = 1 at the top, and k = nzocn at the bottom

    real(dp), dimension(nzocn,nx,ny), intent(inout) :: &
         thermal_forcing         ! 3D ocean thermal forcing to be extrapolated

    ! local variables

    real(dp), dimension(nzocn,nx,ny) :: &
         phi                     ! temporary copy of thermal_forcing

    integer, dimension(nzocn,nx,ny) :: &
         mask                    ! = 1 for filled cells/levels, = 0 for unfilled cells/levels

    integer, dimension(nx,ny) ::  &
         ktop,                 & ! top ocean layer in a cell to be filled, if possible
         kbot                    ! bottom ocean layer in a cell to be filled, if possible

    integer :: &
         sum_mask                ! sum of mask over neighbor cells at a given level

    real(dp) ::  &
         sum_phi                 ! sum of mask*thermal_forcing over neighbor cells at a given level
       
    integer :: &
         max_iter,             & ! max(nx,ny) * max(ewtasks, nxtasks)
         local_count,          & ! local counter for filled values
         global_count,         & ! global counter for filled values
         global_count_save       ! global counter for filled values from previous iteration

    integer :: i, j, k, iter
    integer :: iglobal, jglobal
    integer :: kw, ke, ks, kn           ! ocean level in neighbor cells
    real(dp) :: phiw, phie, phin, phis  ! field value in neighbor cells
    character(len=128) :: message
    logical :: filled                   ! true if filled with a valid value

    logical, parameter :: verbose_extrapolate = .false.  ! set to T to follow progress of each iteration

    ! Note: If thermal forcing is close to but not quite equal to unphys_val (e.g., because of roundoff error),
    !       it is interpreted as equal to unphys_val.
    real(dp), parameter :: tf_roundoff_threshold = 1.0d0  ! roundoff error threshold for thermal_forcing (deg K)

    ! Count the number of filled levels/cells with valid values in the input thermal_forcing.

    local_count = 0
    do j = 1+nhalo, ny-nhalo
       do i = 1+nhalo,  nx-nhalo
          do k = 1, nzocn
             if (abs(thermal_forcing(k,i,j) - unphys_val) < tf_roundoff_threshold) then
                thermal_forcing(k,i,j) = unphys_val
             else  ! real physical value
                local_count = local_count + 1
             endif
          enddo
       enddo
    enddo

    global_count_save = parallel_reduce_sum(local_count)

    ! For each marine-connected cell, compute the top and bottom layers where we need ocean data
    ! (either in the original input field, or extrapolated).
    ! Note: We define ktop and kbot >= 1 for all marine-connected cells, including 
    !       grounded marine-based cells with lsrf = topg
    !       (which might have an ocean connection via a thin water layer near the bed).

    ktop(:,:) = 0
    kbot(:,:) = 0

    do j = 1+nhalo, ny-nhalo
       do i = 1+nhalo,  nx-nhalo
          if (marine_connection_mask(i,j) == 1) then   ! there is a marine path to the ocean

             ! ktop is the layer just above lsrf
             if (lsrf(i,j) <= zocn(nzocn)) then ! lsrf lies at or below the bottom ocean layer
                ktop(i,j) = nzocn
             else
                do k = 1, nzocn
                   if (lsrf(i,j) > zocn(k)) then  ! this ocean layer lies below the lower ice surface
                      ktop(i,j) = max(1,k-1)      ! want ocean data from the layer above
                      exit
                   endif
                enddo
             endif

             ! kbot is the layer just below topg
             if (topg(i,j) >= zocn(1)) then     ! topg lies at or above the top ocean layer
                kbot(i,j) = 1
             else
                do k = nzocn, 1, -1
                   if (topg(i,j) < zocn(k)) then  ! this ocean layer lies above the bed topography
                      kbot(i,j) = min(nzocn,k+1)  ! want ocean data from the layer below
                      exit
                   endif
                enddo
             endif

          endif
       enddo   ! i
    enddo   ! j

    !TODO - Are these halo updates needed?
    call parallel_halo(ktop, parallel)
    call parallel_halo(kbot, parallel)

    if (verbose_bmlt_float .and. this_rank == rtest) then
       print*, ' ' 
       print*, 'ktop: itest, jtest, rtest =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') ktop(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' ' 
       print*, 'kbot: itest, jtest, rtest =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') kbot(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif
    
    ! Compute the maximum number of iterations.
    ! In the worst case, the initial field is filled only at one corner of the global domain and
    !  must be extrapolated to the opposite corner, one cell at a time.
    ! Halo updates after each iteration communicate filled values to neighboring tasks.
    ! Typically for Antarctica, the number of iterations required is ~100
    !  on an 8 km grid (to reach the farthest corners of the Ross Ice Shelf), with
    !  counts doubling for each halving of grid size.

    max_iter = max(parallel%ewtasks, parallel%nstasks) * max(nx-2*nhalo, ny-2*nhalo)

    ! Extrapolate the data horizontally.

    do iter = 1, max_iter

       ! Create a mask, = 1 for filled cells/levels and = 0 for unfilled cells/levels

       where (thermal_forcing == unphys_val)
          mask = 0
       elsewhere
          mask = 1
       endwhere

       if (verbose_extrapolate .and. this_rank == rtest) then
          print*, ' ' 
          print*, 'Iteration =', iter
          do k = kmin_diag, kmax_diag
             print*, ' '
             print*, 'thermal_forcing: k, zocn =', k, zocn(k)
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.2)',advance='no') thermal_forcing(k,i,j)
                enddo
                write(6,*) ' '
             enddo
          enddo
       endif

       ! Create a temporary copy of the TF field.
       phi(:,:,:) = thermal_forcing(:,:,:)
            
       ! Loop through all locally owned cells, filling levels ktop:kbot in unfilled cells
       !  that have one or more filled neighbors at the corresponding levels.
       ! In the end, all cells connected to the ocean via a path through cell edges should
       !  have levels filled between lsrf and topg.  Disconnected inland lakes will not be filled.
       ! Can think of the extrapolation as a crude version of circulation in a C-grid ocean model.

       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo,  nx-nhalo
             if (marine_connection_mask(i,j) == 1) then   ! there is a marine path to the ocean

                do k = ktop(i,j), kbot(i,j)
                   if (mask(k,i,j) == 0) then   ! not yet filled

                      filled = .false.

                      ! Set thermal_forcing(k,i,j) to the mean value in filled neighbors at the same level

                      if (k == ktop(i,j)) then

                         ! Need extra logic to allow spreading of values downward to greater depth,
                         !  in case kbot in a neighbor cell lies above ktop in this cell
                         ! Note: Thermal forcing is corrected for the elevation difference,
                         !       assuming linear dependence of Tf on zocn.

                         if (kbot(i-1,j) > 0 .and. kbot(i-1,j) < k) then  ! kbot in west neighbor lies above k in this cell
                            kw = kbot(i-1,j)
                            phiw = phi(kw,i-1,j) - dtocnfrz_dz * (zocn(kw) - zocn(k))
                         else
                            kw = k
                            phiw = phi(k,i-1,j)
                         endif

                         if (kbot(i+1,j) > 0 .and. kbot(i+1,j) < k) then  ! kbot in east neighbor lies above k in this cell
                            ke = kbot(i+1,j)
                            phie = phi(ke,i+1,j) - dtocnfrz_dz * (zocn(ke) - zocn(k))
                         else
                            ke = k
                            phie = phi(k,i+1,j)
                         endif

                         if (kbot(i,j-1) > 0 .and. kbot(i,j-1) < k) then  ! kbot in south neighbor lies above k in this cell
                            ks = kbot(i,j-1)
                            phis = phi(ks,i,j-1) - dtocnfrz_dz * (zocn(ks) - zocn(k))
                         else
                            ks = k
                            phis = phi(k,i,j-1)
                         endif

                         if (kbot(i,j+1) > 0 .and. kbot(i,j+1) < k) then  ! kbot in north neighbor lies above k in this cell
                            kn = kbot(i,j+1)
                            phin = phi(kn,i,j+1) - dtocnfrz_dz * (zocn(kn) - zocn(k))
                         else
                            kn = k
                            phin = phi(k,i,j+1)
                         endif

                         sum_mask = mask(kw,i-1,j) + mask(ke,i+1,j) + mask(ks,i,j-1) + mask(kn,i,j+1)

                         if (sum_mask > 0) then
                            sum_phi = mask(kw,i-1,j)*phiw + mask(ke,i+1,j)*phie &
                                    + mask(ks,i,j-1)*phis + mask(kn,i,j+1)*phin
                            thermal_forcing(k,i,j) = sum_phi / real(sum_mask, dp)
                            filled = .true.
                         endif

                      endif   ! k = ktop and mask = 0

                      ! Note: In the case k = ktop = kbot, we allow for both upward and downward spreading.
                      !       If there is no downward spreading to ktop (handled by the logic above,
                      !        then we look for upward spreading to kbot (handled by the logic below).

                      if (k == kbot(i,j) .and. .not.filled) then

                         ! need extra logic to allow spreading of values upward to shallower depth,
                         ! in case ktop in a neighbor cell lies below kbot in this cell

                         if (ktop(i-1,j) > k) then  ! ktop in west neighbor lies below k in this cell
                            kw = ktop(i-1,j)
                            phiw = phi(kw,i-1,j) - dtocnfrz_dz * (zocn(kw) - zocn(k))
                         else
                            kw = k
                            phiw = phi(k,i-1,j)
                         endif

                         if (ktop(i+1,j) > k) then  ! ktop in east neighbor lies below k in this cell
                            ke = ktop(i+1,j)
                            phie = phi(ke,i+1,j) - dtocnfrz_dz * (zocn(ke) - zocn(k))
                         else
                            ke = k
                            phie = phi(k,i+1,j)
                         endif

                         if (ktop(i,j-1) > k) then  ! ktop in south neighbor lies below k in this cell
                            ks = ktop(i,j-1)
                            phis = phi(ks,i,j-1) - dtocnfrz_dz * (zocn(ks) - zocn(k))
                         else
                            ks = k
                            phis = phi(k,i,j-1)
                         endif

                         if (ktop(i,j+1) > k) then  ! ktop in north neighbor lies below k in this cell
                            kn = ktop(i,j+1)
                            phin = phi(kn,i,j+1) - dtocnfrz_dz * (zocn(kn) - zocn(k))
                         else
                            kn = k
                            phin = phi(k,i,j+1)
                         endif

                         sum_mask = mask(kw,i-1,j) + mask(ke,i+1,j) + mask(ks,i,j-1) + mask(kn,i,j+1)

                         if (sum_mask > 0) then
                            sum_phi = mask(kw,i-1,j)*phiw + mask(ke,i+1,j)*phie &
                                    + mask(ks,i,j-1)*phis + mask(kn,i,j+1)*phin
                            thermal_forcing(k,i,j) = sum_phi / real(sum_mask, dp)
                            filled = .true.
                         endif

                      endif   ! k = kbot and mask = 0

                      if (k /= ktop(i,j) .and. k /= kbot(i,j) .and. .not.filled) then

                         ! simpler case; look only at neighbor levels with the same k value

                         sum_mask = mask(k,i-1,j) + mask(k,i+1,j) + mask(k,i,j-1) + mask(k,i,j+1)

                         if (sum_mask > 0) then
                            sum_phi = mask(k,i-1,j)*phi(k,i-1,j) + mask(k,i+1,j)*phi(k,i+1,j)   &
                                    + mask(k,i,j-1)*phi(k,i,j-1) + mask(k,i,j+1)*phi(k,i,j+1)
                            thermal_forcing(k,i,j) = sum_phi / real(sum_mask, dp)
                            filled = .true.
                         endif

                      endif

                   endif   ! mask(k,i,j) = 0

                enddo   ! k
             endif   ! ktop >=1, kbot >= 1
          enddo   ! i
       enddo  ! j

       ! Extend TF downward in columns that contain unfilled cells below filled cells.

       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo,  nx-nhalo
             if (marine_connection_mask(i,j) == 1) then   ! there is a marine path to the ocean
                do k = ktop(i,j)+1, kbot(i,j)
                   if (thermal_forcing(k,i,j) == unphys_val .and. thermal_forcing(k-1,i,j) /= unphys_val) then
                      thermal_forcing(k,i,j) = thermal_forcing(k-1,i,j) - dtocnfrz_dz * (zocn(k) - zocn(k-1))
                   endif
                enddo   ! k
             endif   ! ktop >=1, kbot >= 1
          enddo   ! i
       enddo  ! j

       ! Extend TF upward in columns that contain unfilled cells above filled cells.
       ! This is less common than having unfilled cells below filled cells, but it happens occasionally.

       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo,  nx-nhalo
             if (marine_connection_mask(i,j) == 1) then   ! there is a marine path to the ocean
                do k = kbot(i,j)-1, ktop(i,j), -1
                   if (thermal_forcing(k,i,j) == unphys_val .and. thermal_forcing(k+1,i,j) /= unphys_val) then
                      thermal_forcing(k,i,j) = thermal_forcing(k+1,i,j) - dtocnfrz_dz * (zocn(k) - zocn(k+1))
                   endif
                enddo   ! k
             endif   ! ktop >=1, kbot >= 1
          enddo   ! i
       enddo  ! j

       call parallel_halo(thermal_forcing, parallel)

       ! Insert an unphysical value at the global boundary.
       ! This is done to handle the case that global_bc = no_ice,
       !  which puts zeroes in global boundary cells.
       ! We do not want these zeroes to be interpreted as realistic thermal_forcing values.
       call parallel_boundary_value(thermal_forcing, unphys_val, parallel)

       ! Every several iterations, count the total number of filled cells/levels in the global domain.
       ! If this number has not increased since the previous iteration, then exit the loop.
       ! Note: Typically there are several ice sheet dynamics time steps per mass balance time step.
       !       The thermal forcing field from the ocean model (for OCEAN_DATA_GLAD) is updated
       !        at the start of the mass balance time step.  Just after this update, the extrapolation
       !        typically requires 100 or more iterations to converge.  But subsequent updates within this
       !        mass balance time step may require < 5 iterations.  For this reason, we check frequently
       !        for convergence in the first few steps, and then less frequently as the number of iterations grows.

       ! Check after iterations 1, 2, and 5, then after 10, 20, etc.
       if (iter == 1 .or. iter == 2 .or. iter == 5 .or. mod(iter, 10) == 0) then
          local_count = 0
          do j = 1+nhalo, ny-nhalo
             do i = 1+nhalo,  nx-nhalo
                if (marine_connection_mask(i,j) == 1) then   ! there is a marine path to the ocean
                   do k = ktop(i,j), kbot(i,j)
                      if (thermal_forcing(k,i,j) /= unphys_val) local_count = local_count + 1
                   enddo
                endif
             enddo
          enddo

          global_count = parallel_reduce_sum(local_count)

          if (global_count == global_count_save) then
             if (verbose_bmlt_float .and. this_rank == rtest) &
                  print*, 'Extrapolation converged: iter, global_count =', iter, global_count
             exit
          else
             if (verbose_bmlt_float .and. this_rank == rtest) &
                  print*, 'Extrapolation convergence check: iter, global_count =', iter, global_count
             global_count_save = global_count
          endif
          
       endif   ! time for a convergence check

       if (iter == max_iter) then
          print*, 'iter = max_iter:', max_iter
          call write_log('Ocean extrapolation error; number of filled cells has not plateaued', GM_FATAL)
       endif

    enddo   ! max_iter

    ! Bug check:
    ! Make sure all levels from ktop to kbot are filled in cells with thermal_forcing_mask = 1.

    do j = 1, ny
       do i = 1, nx
          if (thermal_forcing_mask(i,j) == 1) then
             do k = ktop(i,j), kbot(i,j)
                if (thermal_forcing(k,i,j) == unphys_val) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   print*, 'i, j, ktop, kbot =', i, j, ktop(i,j), kbot(i,j)
                   write(message,*) 'Ocean data extrapolation error: unphys value in level k, i, j:', &
                        k, iglobal, jglobal, thermal_forcing(k,i,j)
                   call write_log(message, GM_FATAL)
                endif
             enddo   ! k
          endif   ! thermal_forcing_mask
       enddo   ! i
    enddo   ! j

    ! Set the thermal forcing to a default value in the remaining unfilled cell/levels.
    ! For example, we might set thermal_forcing = 0 to get cleaner diagnostics.

    if (default_val /= unphys_val) then
       where (thermal_forcing == unphys_val)
          thermal_forcing = default_val
       endwhere
    endif

  end subroutine glissade_thermal_forcing_extrapolate

!****************************************************

  subroutine interpolate_thermal_forcing_to_lsrf(&
       nx,          ny,          &
       nzocn,       zocn,        &
       thermal_forcing_mask,     &
       lsrf,                     &
       thermal_forcing,          &
       thermal_forcing_lsrf)

    ! Interpolate the ocean thermal forcing field to the lower ice surface.

    integer, intent(in) :: &
         nx, ny                    !> number of grid cells in each dimension

    integer, intent(in) :: &
         nzocn                     !> number of ocean levels

    real(dp), dimension(nzocn), intent(in) :: &
         zocn                      !> ocean levels (m) where forcing is provided, negative below sea level

    integer, dimension(nx,ny), intent(in) :: &
         thermal_forcing_mask      !> = 1 if ice is present and floating, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         lsrf                      !> ice lower surface elevation (m), negative below sea level

    real(dp), dimension(nzocn,nx,ny), intent(in) :: &
         thermal_forcing           !> thermal forcing field at ocean levels

    real(dp), dimension(nx,ny), intent(out) :: &
         thermal_forcing_lsrf      !> thermal forcing at the lower ice surface

    ! local veriables

    integer :: i, j, k
    integer :: iglobal, jglobal
    real(dp) :: dtf, dzocn, dzice  ! terms used in linear interpolation

    ! Compute the thermal forcing at the lower ice surface.
    ! Above the top ocean level, use the TF value at the top level.
    ! Below the bottom ocean level, use the TF value at the bottom level.
    ! Use linear interpolation in between.

    do j = 1, ny
       do i = 1, nx
          if (thermal_forcing_mask(i,j) == 1) then
             if (lsrf(i,j) >= zocn(1)) then
                thermal_forcing_lsrf(i,j) = thermal_forcing(1,i,j)
             elseif (lsrf(i,j) < zocn(nzocn)) then
                thermal_forcing_lsrf(i,j) = thermal_forcing(nzocn,i,j)
             else
                do k = 1, nzocn-1
                   if (lsrf(i,j) < zocn(k) .and. lsrf(i,j) >= zocn(k+1)) then
                      dtf = thermal_forcing(k+1,i,j) - thermal_forcing(k,i,j)
                      dzocn = zocn(k+1) - zocn(k)
                      dzice = lsrf(i,j) - zocn(k)
                      thermal_forcing_lsrf(i,j) = thermal_forcing(k,i,j) + (dzice/dzocn) * dtf
                      exit
                   endif
                enddo
             endif
          else  ! not a floating cell connected to the ocean
             thermal_forcing_lsrf(i,j) = 0.0d0
          endif

       enddo
    enddo

  end subroutine interpolate_thermal_forcing_to_lsrf

!****************************************************

  subroutine ismip6_bmlt_float(&
       bmlt_float_thermal_forcing_param, &
       nx,         ny,            &
       itest,   jtest,    rtest,  &
       nbasin,                    &
       basin_number,              &
       gamma0,                    &
       thermal_forcing_lsrf,      &
       deltaT_ocn,                &
       thermal_forcing_basin,     &
       deltaT_basin_avg,          &
       thermal_forcing_mask,      &
       bmlt_float)

    ! Compute the basal melt rate as a quadratic function of thermal forcing, using either
    !  a local or nonlocal parameterization as specified for ISMIP6.
    ! Note: The input thermal forcing fields are nonnegative, but the deltaT correction can be negative
    !       giving negative effective TF, which is nonphysical with a quadratic paramterization.
    !       In this case, we set effective TF = 0.

    integer, intent(in) :: &
         bmlt_float_thermal_forcing_param  !> kind of melting parameterization, local or nonlocal

    integer, intent(in) :: &
         nx, ny                   !> number of grid cells in each dimension

    integer, intent(in) :: &
         itest, jtest, rtest      !> coordinates of diagnostic point

    integer, intent(in) :: &
         nbasin                   !> number of basins

    integer, dimension(nx,ny), intent(in) :: &
         basin_number             !> integer ID for each basin

    real(dp), intent(in) :: &
         gamma0                   !> basal melt rate coefficient (m/yr)

    real(dp), dimension(nx,ny), intent(in) :: &
         thermal_forcing_lsrf,  & !> thermal forcing (K) at lower ice surface
         deltaT_ocn               !> thermal forcing correction factor (deg C)

    real(dp), dimension(nbasin), intent(in) :: &
         thermal_forcing_basin, & !> thermal forcing averaged over each basin (K)
         deltaT_basin_avg         !> thermal forcing correction factor for each basin (deg C)

    integer, dimension(nx,ny), intent(in) :: &
         thermal_forcing_mask     !> = 1 where TF-driven bmlt_float can be > 0

    real(dp), dimension(nx,ny), intent(out) :: &
         bmlt_float               !> basal melt rate (m/yr) at lower ice surface

    ! local variables

    integer :: i, j, nb

    real(dp) :: coeff         ! constant coefficient = [(rhow*cp)/(rhoi*Lf)]^2, with units deg^(-2)

    ! ISMIP6 prescribed parameters
    real(dp), parameter ::  &
         rhoi_ismip6 = 918.0d0,    & ! ice density (kg/m^3)
         rhosw_ismip6 = 1028.0d0,  & ! seawater density (kg/m^3)
         Lf_ismip6 = 3.34d5,              & ! latent heat of fusion (J/kg)
         cpw_ismip6 = 3974.d0               ! specific heat of seawater (J/kg/K)

    real(dp) :: &
         eff_thermal_forcing,      & ! effective local thermal forcing, after deltaT correction
         eff_thermal_forcing_basin   ! effective basin thermal forcing, after deltaT correction

    ! initialize
    bmlt_float(:,:) = 0.0d0

    coeff = ( (rhosw_ismip6*cpw_ismip6)/(rhoi_ismip6*Lf_ismip6) )**2

    if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL) then

       ! local parameterization
       ! melt rate is a quadratic function of local thermal forcing

       do j = 1, ny
          do i = 1, nx
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing = max(0.0d0, thermal_forcing_lsrf(i,j) + deltaT_ocn(i,j))
                bmlt_float(i,j) = coeff * gamma0 * eff_thermal_forcing**2
             endif
          enddo
       enddo

    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or. &
            bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

       ! nonlocal parameterization
       ! melt rate is a quadratic function of local thermal forcing and basin-average thermal forcing
       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (thermal_forcing_mask(i,j) == 1) then
                ! Note: Can have bmlt_float < 0 where thermal_forcing_lsrf + deltaT_ocn < 0
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn(i,j)
                eff_thermal_forcing_basin = max(0.0d0, thermal_forcing_basin(nb) + deltaT_basin_avg(nb))
                bmlt_float(i,j) = coeff * gamma0 * eff_thermal_forcing * eff_thermal_forcing_basin

                !WHL - debug
                if (verbose_bmlt_float .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'In ismip6_bmlt_float, r, i, j, nb =', rtest, itest, jtest, nb
                   print*, 'gamma0, coeff =', gamma0, coeff
                   print*, 'thermal_forcing_lsrf =', thermal_forcing_lsrf(i,j)
                   print*, 'deltaT_ocn =', deltaT_ocn(i,j)
                   print*, 'thermal_forcing_basin =', thermal_forcing_basin(nb)
                   print*, 'deltaT_basin_avg =', deltaT_basin_avg(nb)
                   print*, 'eff_TF, eff_TF_basin =', eff_thermal_forcing, eff_thermal_forcing_basin
                   print*, 'bmlt_float =', bmlt_float(i,j)
                endif

             endif
          enddo
       enddo

    endif  ! local or nonlocal

  end subroutine ismip6_bmlt_float

!****************************************************

  subroutine quadratic_bmlt_float(&
       nx,              ny,   &
       thermal_forcing_lsrf,  &
       thermal_forcing_mask,  &
       bmlt_float)

    ! GL: 04-29-2019
    ! Compute the basal melt rate as a quadratic function of thermal_forcing,
    !  similarly to Pollard & DeConto (2012) and Martin et al. (2011).

    !TODO - Add specific heat of seawater to glimmer_physcon?
    use glimmer_physcon, only : rhoi, rhoo, lhci

    integer, intent(in) :: &
         nx, ny                   !> number of grid cells in each dimension

    real(dp), dimension(nx,ny), intent(in) :: &
         thermal_forcing_lsrf     !> thermal forcing (deg C) at lower ice surface

    integer, dimension(nx,ny), intent(in), optional :: &
         thermal_forcing_mask     !> = 1 where thermal forcing and bmlt_float can be nonzero

    real(dp), dimension(nx,ny), intent(out) :: &
         bmlt_float               !> basal melt rate (m/yr) at lower ice surface

    ! local variables

    integer :: i, j

    ! prescribed parameters

    real(dp), parameter ::    &
         cpo     = 3974.d0,   & ! specific heat of seawater (J/kg/K)
         gamma_t = 1.0d-4,    & ! thermal exchange velocity of ocean water (m/s)
         Fm      = 5.0d-3       ! dimensionless parameter for tuning purposes

    real(dp) :: &
         thermal_forcing        ! applied thermal forcing, constrained to be >= 0

    ! initialize
    bmlt_float(:,:) = 0.0d0

    ! Compute basal melt rates
    ! Note: gamma_t has units of m/s, but bmlt_float is given in m/yr, so multiply by scyr..
    ! Note: The applied thermal forcing must be non-negative.  TF < 0 is inappropriate for a quadratic function.
    do j = 1, ny
       do i = 1, nx
          if (thermal_forcing_mask(i,j) == 1) then
             thermal_forcing = max(thermal_forcing_lsrf(i,j), 0.0d0)
             bmlt_float(i,j) = rhoo * cpo * (gamma_t*scyr) * Fm * thermal_forcing**2 / (lhci * rhoi)
          endif
       enddo
    enddo

  end subroutine quadratic_bmlt_float

!****************************************************

  subroutine basin_number_extrapolate(&
           nx,         ny,  &
           parallel,        &
           nbasin,          &
           basin_number)

    integer, intent(in) :: &
         nx, ny                     !> number of grid cells in each direction

    type(parallel_type), intent(in) :: &
         parallel                   !> info for parallel communication

    integer, intent(in) :: &
         nbasin                     !> number of basins

    integer, dimension(nx,ny), intent(inout) :: &
         basin_number               !> integer basin ID for each cell


    integer :: basin_number_min     ! global minval for basin_number

    ! local variables

    integer :: local_count          ! number of cells with valid values
    integer :: global_count         ! number of global cells with valid values
    integer :: global_count_save    ! number of global cells with valid values in previous iteration
    integer :: iter                 ! iteration counter
    integer :: max_iter             ! max number of iterations given the domain size

    integer :: i, j

    integer, dimension(nx,ny) :: &
         valid_mask                 ! mask of cells valid basin numbers

    real(dp), dimension(nx,ny) :: &
         basin_number_new           ! work array for basin number

    logical, parameter ::  &
         verbose_basin_number = .true.

    ! Count the number of cells with valid basin numbers

    valid_mask = 0
    local_count = 0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1,  nx-nhalo
          if (basin_number(i,j) >= 1 .and. basin_number(i,j) <= nbasin) then
             valid_mask(i,j) = 1
             local_count = local_count + 1
          endif
       enddo
    enddo

    global_count_save = parallel_reduce_sum(local_count)

    call parallel_halo(valid_mask, parallel)
    call parallel_halo(basin_number, parallel)

    ! Compute the maximum number of iterations.
    ! In the worst case, the initial field is filled only at one corner of the global domain and
    !  must be extrapolated to the opposite corner, one cell at a time.
    ! Halo updates after each iteration communicate valid values to neighboring tasks.

    max_iter = max(parallel%ewtasks, parallel%nstasks) * max(nx-2*nhalo, ny-2*nhalo)

    if (verbose_basin_number .and. main_task) then
       print*, 'Extrapolating basin numbers to cells with invalid values'
       print*, 'Initial count of valid values:', global_count_save
       print*, 'max_iter =', max_iter
    endif

    ! Extrapolate the data horizontally

    do iter = 1, max_iter

       ! Extrapolate valid values by one cell in each direction
       basin_number_new = basin_number

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1,  nx-nhalo
             if (basin_number(i,j) < 1 .or. basin_number(i,j) > nbasin) then  ! invalid value
                if (valid_mask(i+1,j) == 1) then
                   basin_number_new(i,j) = basin_number(i+1,j)
                elseif (valid_mask(i-1,j) == 1) then
                   basin_number_new(i,j) = basin_number(i-1,j)
                elseif (valid_mask(i,j+1) == 1) then
                   basin_number_new(i,j) = basin_number(i,j+1)
                elseif (valid_mask(i,j-1) == 1) then
                   basin_number_new(i,j) = basin_number(i,j-1)
                endif
             endif
          enddo
       enddo

       basin_number = basin_number_new
       call parallel_halo(basin_number, parallel)

       ! Count the number of valid values and recompute the mask

       valid_mask = 0
       local_count = 0

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1,  nx-nhalo
             if (basin_number(i,j) >= 1 .and. basin_number(i,j) <= nbasin) then
                valid_mask(i,j) = 1
                local_count = local_count + 1
             endif
          enddo
       enddo

       global_count = parallel_reduce_sum(local_count)
       call parallel_halo(valid_mask, parallel)

       if (verbose_basin_number .and. main_task) then
          print*, iter, 'Basin number count =', global_count
       endif

       if (global_count == global_count_save) then
          if (verbose_basin_number .and. main_task) then
             print*, 'Exiting basin_number_extrapolate, iter =', iter
          endif
          exit
       else
          global_count_save = global_count
       endif

       if (iter == max_iter) then
          if (main_task) print*, 'Error: Exiting basin_number_extrapolate, max_iter =', max_iter
          call write_log('Error: Exiting basin_number_extrapolate without converging', GM_FATAL)
       endif

    enddo   ! iter

  end subroutine basin_number_extrapolate

!****************************************************

  ! Note: Old plume subroutines were here, but most have been removed.
  !       I kept some utility subroutines that might be useful later.

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
    
    character(len=128) :: message

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
    
  end subroutine compute_plume_velocity

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

