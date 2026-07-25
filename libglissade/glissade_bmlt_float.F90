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
  use glimmer_physcon, only: rhoi, rhoo, rhow, grav, lhci, cpw, scyr
  use glimmer_paramets, only: iulog, unphys_val
  use glimmer_log
  use glide_types
  use glimmer_utils, only: point_diag
  use cism_parallel, only: this_rank, main_task, nhalo, &
       parallel_type, parallel_halo, parallel_globalindex, parallel_boundary_value, &
       parallel_reduce_sum, parallel_reduce_min, parallel_reduce_max, &
       parallel_global_sum, parallel_global_sum_patch, parallel_is_zero

  implicit none
  
  private
  public :: verbose_bmlt_float, glissade_basal_melting_float, &
       glissade_bmlt_float_init, glissade_bmlt_float_solve

    ! Max and min allowed values for thermal forcing
    real(dp), parameter ::  &
         thermal_forcing_max = 20.d0,  &  ! max allowed value of thermal forcing (K)
         thermal_forcing_min = -5.d0      ! min allowed value of thermal forcing (K)

    integer, parameter :: &
         cavity_buffer = 0     ! distance from ice edge (measured in number of grid cells) over which ocean TF values
                               ! are discarded before starting TF extrapolation; must be <= nhalo

    ! loop limits for debug diagnostics
    integer :: kmin_diag = 1
    integer :: kmax_diag = 2

    logical, parameter :: verbose_bmlt_float = .true.

  contains

!****************************************************

  !TODO - Break up this subroutine into shorter subroutines called from glissade_bmlt_float_solve
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

    ! Compute the rate of basal melting for floating ice by one of several methods.
    ! For bmlt_float based on thermal forcing, see subroutine glissade_bmlt_float_thermal_forcing.
    ! TODO: Merge with subroutine glissade_bmlt_float_solve.

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

    if (main_task .and. verbose_bmlt_float) write(iulog,*) 'Computing bmlt_float, whichbmlt_float =', whichbmlt_float

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
                   bmlt_float(i,j) = basal_melt%bmlt_float_const/scyr
                endif

                !WHL - debug
                if (j==jtest .and. this_rank==rtest) then
!!                if (i==itest .and. j==jtest .and. this_rank==rtest) then
!!                   write(iulog,*) 'rank, i, j, bmlt_float:', this_rank, i, j, bmlt_float(i,j)*scyr
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
                ! Note: bmlt_float_omega has units of 1/yr; dividing by scyr converts to 1/s so bmlt_float has units of m/s
                bmlt_float(i,j) = (basal_melt%bmlt_float_omega/scyr) * tanh(h_cavity/basal_melt%bmlt_float_h0) &
                                * max(basal_melt%bmlt_float_z0 - z_draft, 0.0d0)
                !debug
!                if (j == jtest .and. verbose_bmlt_float) then
!                   write(iulog,*) 'cavity, tanh, thck, draft, melt rate (m/yr):', i, j, h_cavity, &
!                         tanh(h_cavity/basal_melt%bmlt_float_h0), thck(i,j), z_draft, bmlt_float(i,j)*scyr
!                endif

             endif   ! ice is present and floating

          enddo
       enddo

       if (verbose_bmlt_float) then
          call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
          call point_diag(bmlt_float*scyr, 'bmlt_float (m/yr)', itest, jtest, rtest, 7, 7)
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

       ! convert from m/yr to m/s
       bmlt_float = bmlt_float/scyr

       if (verbose_bmlt_float) then
          call point_diag(topg, 'topg (m)', itest, jtest, rtest, 7, 7)
          call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
          call point_diag(min(lsrf-eus, 0.0d0), 'z_draft (m)', itest, jtest, rtest, 7, 7)
          call point_diag(max(lsrf-topg, 0.0d0), 'h_cavity (m)', itest, jtest, rtest, 7, 7)
          call point_diag(bmlt_float*scyr, 'bmlt_float (m/yr)', itest, jtest, rtest, 7, 7)
       endif

!    elseif (whichbmlt_float == BMLT_FLOAT_PLUME) then

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

  !TODO - Maybe change to glissade_thermal_forcing_init, since it is called only when
  !       whichbmlt_float = BMLT_FLOAT_THERMAL_FORCING

  subroutine glissade_bmlt_float_init(model, ocean_data)

    use glimmer_paramets, only: unphys_val
    use glissade_masks, only : glissade_get_masks
    use glissade_grounding_line, only: glissade_grounded_fraction
    use glissade_utils, only: glissade_basin_average

    ! Initialization for basal melting
    ! Currently called only if model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING)

    type(glide_global_type), intent(inout) :: model   !> derived type holding ice-sheet info

    type(glide_ocean_data), intent(inout) ::  &
         ocean_data                                   !> derived type holding ocean input data
    
    integer, dimension(model%general%ewn, model%general%nsn) ::  &
         ice_mask,                    & ! = 1 if ice is present (thck > 0)
         floating_mask,               & ! = 1 where ice is present and floating, else = 0
         ocean_mask,                  & ! = 1 if topg is below sea level and ice is absent, else = 0
         land_mask                      ! = 1 if topg is at or above sea level, else = 0

    real(dp), dimension(model%ocean_data%nbasin) :: &
         thermal_forcing_basin,      &  ! basin average thermal forcing (K)
         dthck_dt_basin                 ! basin average of dthck_dt_obs (m/s)

    real(dp) :: tot_bmlt_float_target   ! target for total global bmlt_float (Gt/yr)
    real(dp), dimension(model%ocean_data%nbasin) :: bmlt_float_target_basin   ! basin targets (Gt/yr)
    integer :: itest, jtest, rtest      ! coordinates of diagnostic point
    integer :: i, j, k, nb
    integer :: ewn, nsn
    real(dp) :: dew, dns
    integer :: imin, imax, jmin, jmax
    integer :: buffer

    integer :: basin_number_min         ! global minval of the basin_number field

    type(parallel_type) :: parallel     ! info for parallel communication

    parallel = model%parallel

    ewn = model%general%ewn
    nsn = model%general%nsn

    dew = model%numerics%dew
    dns = model%numerics%dns

    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    if (verbose_bmlt_float .and. main_task) then
       write(iulog,*) 'In glissade_bmlt_float_thermal_forcing_init'
       write(iulog,*) 'bmlt_float_thermal_forcing_param =', model%options%bmlt_float_thermal_forcing_param
    endif

    ! update some masks
    call glissade_get_masks(&
         ewn,                 nsn,                   &
         parallel,                                   &
         model%geometry%thck, model%geometry%topg,   &
         model%climate%eus,   0.0d0,                 &  ! thklim = 0
         ice_mask,                                   &
         floating_mask = floating_mask,              &
         land_mask = land_mask,                      &
         ocean_mask = ocean_mask)

    ! update the grounded fraction, f_ground_cell
    ! used for some but not all of the options below; compute here in case it's needed
    call glissade_grounded_fraction(&
         ewn,                nsn,          &
         parallel,                      &
         itest, jtest, rtest,           &  ! diagnostic only
         model%geometry%thck,           &
         model%geometry%topg,           &
         model%climate%eus,             &
         ice_mask,                      &
         floating_mask,                 &
         land_mask,                     &
         model%options%which_ho_ground, &
         model%options%which_ho_flotation_function, &
         model%options%which_ho_fground_no_glp,     &
         model%geometry%f_flotation,    &
         model%geometry%f_ground,       &
         model%geometry%f_ground_cell)

    ! Note: Most of the following code is ignored on a standard restart but applied with a hybrid restart.

    if (model%options%is_restart == NO_RESTART .or. model%options%is_restart == HYBRID_RESTART) then

       !----------------------------------------------------------------------
       ! Optionally, prepare the ocean thermal forcing for extrapolation to ice shelf cavities.
       ! For standalone CISM runs, this is done at initialization (but not at restart) in 2 steps:
       ! (1) Set thermal_forcing = unphys_val in cavities. I.e., remove any values that are already present
       !     and would otherwise be interpreted as valid values.
       ! (2) Apply the extrapolation algorithm. Thermal forcing values beyond the calving front
       !     are extrapolated into the cavity, as far as the grounding line.
       ! We do (1) now. Then during the first timestep, we do (2). On subsequent timesteps we repeat (2),
       !  but starting from existing values in the cavity, so as not to repeat the computationally
       !  intensive extrapolation from scratch.
       ! This assumes that the far-field ocean values do not change during the run.
       ! If applying ocean anomalies in the far field and extrapolating the new thermal forcing into cavities,
       !  we would need some new logic to go back to step (1).
       !
       ! For coupled ESM runs with ocean_data_domain = OCEAN_DATA_GLAD), thermal forcing is received
       !  once per mass balance time step. In this case, we would want to apply step (1) before (2)
       !  each time CISM gets new ocean data. This will take some additional logic.
       !----------------------------------------------------------------------

       if (model%options%ocean_data_extrapolate == OCEAN_DATA_EXTRAPOLATE_TRUE) then

          ! Set TF = unphys_val everywhere except open ocean.
          ! Optionally, use a buffer to set TF = unphys_val for ocean cells
          !        within 1 or 2 cells of the ice edge.
          !       As the code is written now, the buffer cannot be larger than nhalo.

          buffer = min(cavity_buffer, nhalo)

          do j = 1+nhalo, nsn-nhalo
             do i = 1+nhalo, ewn-nhalo
                imin = i-buffer; imax = i+buffer
                jmin = j-buffer; jmax = j+buffer
                if (any(ice_mask(imin:imax,jmin:jmax) == 1) .or. any(ocean_mask(imin:imax,jmin:jmax) == 0)) then
                   ocean_data%thermal_forcing(:,i,j) = unphys_val
                end if
             end do
          end do

          call parallel_halo(ocean_data%thermal_forcing, parallel)

          if (verbose_bmlt_float) then
             if (this_rank ==rtest) write(iulog,*) 'Set TF = unphys_val in cavities, buffer =', buffer
             do k = kmin_diag, kmax_diag
                if (this_rank == rtest) write(iulog,*) 'k =', k
                call point_diag(ocean_data%thermal_forcing(k,:,:), 'TF before extrapolating', itest, jtest, rtest, 7, 7)
             enddo
          endif

       endif  ! ocean_data_extrapolate

       ! initialization for the ISMIP6 basal melt options

       if (model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL .or.  &
           model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or. &
           model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

          if (verbose_bmlt_float) then
             if (this_rank==rtest) then
                write(iulog,*) ' '
                write(iulog,*) 'Initialize ISMIP sub-shelf melting'
                write(iulog,*) ' '
                write(iulog,*) 'k, zocn:'
                do k = 1, ocean_data%nzocn
                   write(iulog,*) k, ocean_data%zocn(k)
                enddo
                write(iulog,*) 'gamma0 =', ocean_data%gamma0
             endif
             call point_diag(ocean_data%basin_number(:,:), 'basin_number', itest, jtest, rtest, 7, 7)
             call point_diag(ocean_data%deltaT_ocn(:,:), 'deltaT_ocn', itest, jtest, rtest, 7, 7)
             do k = kmin_diag, kmax_diag
                if (this_rank == rtest) write(iulog,*) 'k =', k
                call point_diag(ocean_data%thermal_forcing(k,:,:), 'thermal_forcing', itest, jtest, rtest, 7, 7)
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
                write(iulog,*) 'Extrapolate basin numbers'
             endif

             call basin_number_extrapolate(&
                  ewn,             nsn,     &
                  parallel,                 &
                  model%ocean_data%nbasin,  &
                  model%ocean_data%basin_number)

          endif

          ! If bmb_float was read in, then set a melt-rate target and write some diagnostics.
          ! Note: Both bmb_float and bmlt_float_target are defined as positive for melting

          if (.not.parallel_is_zero(model%ocean_data%bmb_float)) then
             ! convert kg/m2/yr w.e. to m/s of ice melt
             model%basal_melt%bmlt_float_target = model%ocean_data%bmb_float / (rhoi*scyr)
             if (verbose_bmlt_float) then
                call point_diag(model%ocean_data%bmb_float, 'bmb_float (kg/m2/yr)', itest, jtest, rtest, 7, 7)
                call point_diag(model%basal_melt%bmlt_float_target*scyr, 'bmlt_float_target (m/yr ice)', itest, jtest, rtest, 7, 7)
                tot_bmlt_float_target = parallel_global_sum(model%basal_melt%bmlt_float_target, parallel)
                bmlt_float_target_basin = parallel_global_sum_patch(model%basal_melt%bmlt_float_target, &
                     model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
                ! convert m/s to Gt/yr
                tot_bmlt_float_target = tot_bmlt_float_target * dew*dns*rhoi*scyr/1.d12
                bmlt_float_target_basin = bmlt_float_target_basin * dew*dns*rhoi*scyr/1.d12
                if (main_task) then
                   write(iulog,*) ' '
                   write(iulog,*) 'tot_bmlt_float_target (Gt/yr) =', tot_bmlt_float_target
                   write(iulog,*) 'Basin targets:'
                   do nb = 1, model%ocean_data%nbasin
                      write(iulog,*) nb, bmlt_float_target_basin(nb)
                   enddo
                endif
             endif
          endif

          ! Optionally, compute basin-scale values of deltaT_ocn to match a basin-average melt rate target.

          if (model%options%which_ho_deltaT_ocn == HO_DELTAT_OCN_CALIBRATE_BASIN) then

             if (parallel_is_zero(model%basal_melt%bmlt_float_target)) then
                call write_log('This value of which_ho_deltaT_ocn requires an input bmlt_float_target', GM_FATAL)
             endif

             ! Call the driver subroutine to calibrate deltaT_ocn and compute the resulting bmlt_float.
             ! With calibrate = T, the driver subroutine will call compute_bmlt_float_thermal_forcing
             !  with which_ho_deltaT_ocn == HO_DELTAT_OCN_CALIBRATE_BASIN.
             ! This is the signal to calibrate deltaT_ocn in the process of computing bmlt_float.
             ! Note: f_ground_cell must be up to date; that's done above

             call glissade_bmlt_float_solve(model, calibrate_in = .true.)

             if (verbose_bmlt_float) then
                call point_diag(model%ocean_data%deltaT_ocn, 'calibrated deltaT_ocn', itest, jtest, rtest, 7, 7)
                call point_diag(model%basal_melt%bmlt_float*scyr, 'bmlt_float after calibration', itest, jtest, rtest, 7, 7)
             endif

          endif   ! calibrate basins

          if (model%options%bmlt_float_init) then

             ! Compute bmlt_float now; then it can be written to output files at t = 0
             call glissade_bmlt_float_solve(model, calibrate_in = .false.)

             if (verbose_bmlt_float) then
                call point_diag(model%ocean_data%deltaT_ocn, 'deltaT_ocn at initialization', itest, jtest, rtest, 7, 7)
                call point_diag(model%basal_melt%bmlt_float*scyr, 'bmlt_float at initialization', itest, jtest, rtest, 7, 7)
             endif

             if (verbose_bmlt_float) then

                ! Compute the average deltaT_ocn in each basin
                call glissade_basin_average(&
                     ewn,             nsn,                        &
                     parallel,                                    &
                     model%ocean_data%nbasin,                     &
                     model%ocean_data%basin_number,               &
                     model%basal_melt%thermal_forcing_mask*1.0d0, &
                     model%ocean_data%deltaT_ocn,                 &
                     model%ocean_data%deltaT_ocn_basin)

                ! Compute the average thermal forcing in each basin (including deltaT_ocn)
                call glissade_basin_average(&
                     ewn,             nsn,                        &
                     parallel,                                    &
                     model%ocean_data%nbasin,                     &
                     model%ocean_data%basin_number,               &
                     model%basal_melt%thermal_forcing_mask*1.0d0, &
                     model%ocean_data%thermal_forcing_lsrf + model%ocean_data%deltaT_ocn, &
                     thermal_forcing_basin)

                ! Compute the average basal melt rate in each basin
                call glissade_basin_average(&
                     ewn,             nsn,                        &
                     parallel,                                    &
                     model%ocean_data%nbasin,                     &
                     model%ocean_data%basin_number,               &
                     model%basal_melt%thermal_forcing_mask*1.0d0, &
                     model%basal_melt%bmlt_float,                 &
                     model%scalars%bmlt_float_basin)

                if (this_rank == rtest) then
                   write(iulog,*) ' '
                   write(iulog,*) 'basin #, deltaT_ocn, TF (w/dT_ocn), bmlt_float (m/yr):'
                   do nb = 1, model%ocean_data%nbasin
                      write(iulog,'(i6,3f14.6)') nb, model%ocean_data%deltaT_ocn_basin(nb), &
                           thermal_forcing_basin(nb), model%scalars%bmlt_float_basin(nb)*scyr
                   enddo
                endif

             endif   ! verbose

          endif   ! compute initial bmlt_float

       endif   ! ISMIP6 melt scheme

       !TODO - Superseded by HO_DELTAT_OCN_CALIBRATE_BASIN?
       ! Optionally, compute the basin average of dthck_dt_obs, the observed rate of thickening/thinning.
       ! When inverting for deltaT_ocn, we can correct acab by applying (-dthck_dt_obs_basin).
       ! This induces a basal melt rate that will drive thinning when the correction is removed.
       ! On restart, dthck_dt_obs_basin is read from the restart file.
       !TODO: Is dthck_dt_obs needed in the restart file after dthck_dt_obs_basin is computed?

       if (model%options%enable_acab_dthck_dt_correction) then

          if (verbose_bmlt_float) then
             call point_diag(model%ocean_data%basin_number, 'basin_number', itest, jtest, rtest, 7, 7)
             call point_diag(floating_mask, 'floating_mask', itest, jtest, rtest, 7, 7)
             call point_diag(model%geometry%dthck_dt_obs*scyr, 'dthck_dt_obs (m/yr)', itest, jtest, rtest, 7, 7)
          endif

          call glissade_basin_average(&
               model%general%ewn, model%general%nsn,   &
               parallel,                               &
               model%ocean_data%nbasin,                &
               model%ocean_data%basin_number,          &
               floating_mask * 1.0d0,                  &   ! real mask
               model%geometry%dthck_dt_obs,            &   ! m/s
               dthck_dt_basin)                             ! m/s

          if (verbose_bmlt_float .and. main_task) then
             write(iulog,*) ' '
             write(iulog,*) 'nb, dthck_dt_basin (m/yr)'
             do nb = 1, model%ocean_data%nbasin
                write(iulog,*) nb, dthck_dt_basin(nb)*scyr
             enddo
          endif

          ! Make sure the basin average <= 0
          dthck_dt_basin(:) = min(dthck_dt_basin(:), 0.0d0)

          ! Assign the basin average to a 2D array
          do j = 1, model%general%nsn
             do i = 1, model%general%ewn
                nb = model%ocean_data%basin_number(i,j)
                model%geometry%dthck_dt_obs_basin(i,j) = dthck_dt_basin(nb)
             enddo
          enddo

       endif   ! enable_acab_dthck_dt_correction

    endif  ! no restart or hybrid restart

    if (model%options%which_ho_deltat_ocn == HO_DELTAT_OCN_DTHCK_DT) then

       ! Set deltaT_ocn based on dthck_dt_obs.
       ! This is done within the subroutine used to compute bmlt_float from thermal forcing.
       ! But instead of computing bmlt_float from TF, we find the value of deltaT_ocn
       !  that will increase TF as needed to match negative values of dthck_dt_obs.
       ! Note: This subroutine would usually be called during the initial diagnostic solve
       !       of the restart following a spin-up, without taking any prognostic timesteps.

       if (parallel_is_zero(model%geometry%dthck_dt_obs)) then
          call write_log('This which_ho_deltaT_ocn option requires an input dthck_dt_obs', GM_FATAL)
       endif

       call compute_bmlt_float_thermal_forcing(&
            model%options%bmlt_float_thermal_forcing_param, &
            model%options%ocean_data_extrapolate,     &
            parallel,                                 &
            ewn,       nsn,                           &
            dew,       dns,                           &   ! m
            itest,     jtest,   rtest,                &
            ice_mask,                                 &
            floating_mask,                            &
            ocean_mask,                               &
            model%geometry%marine_connection_mask,    &
            model%geometry%f_ground_cell,             &
            model%geometry%thck,                      &   ! m
            model%geometry%lsrf,                      &   ! m
            model%geometry%topg,                      &   ! m
            model%ocean_data,                         &
            model%basal_melt%thermal_forcing_mask,    &
            model%basal_melt%bmlt_float,              &   ! m/s
            which_ho_deltaT_ocn = model%options%which_ho_deltaT_ocn,  &
            dthck_dt_obs = model%geometry%dthck_dt_obs)   ! m/s

       if (verbose_bmlt_float) then
          call point_diag(model%geometry%dthck_dt_obs*scyr, 'dthck_dt_obs', itest, jtest, rtest, 7, 7)
          call point_diag(model%basal_melt%bmlt_float*scyr, 'bmlt_float (m/yr)', itest, jtest, rtest, 7, 7)
          call point_diag(model%ocean_data%deltaT_ocn, 'deltaT_ocn', itest, jtest, rtest, 7, 7)
       endif

    endif   ! HO_DELTAT_OCN_DTHCK_DT

  end subroutine glissade_bmlt_float_init

!****************************************************

  subroutine glissade_bmlt_float_solve(model, calibrate_in)

    ! Solve for basal melting beneath floating ice.

    use glimmer_paramets, only: eps08, eps11
    use glimmer_physcon, only: scyr
    use glissade_plume, only: glissade_plume_driver
    use glissade_mass_balance, only: glissade_add_2d_anomaly
    use glissade_masks, only: glissade_get_masks
    use glissade_utils, only: glissade_basin_average
    use cism_parallel, only:  parallel_reduce_max, parallel_is_zero, &
         parallel_global_sum, parallel_global_sum_patch

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    logical, intent(in), optional :: &
         calibrate_in             ! if true, then calibrate melt rates; false by default

    ! Local variables

    logical :: calibrate          ! local version of calibrate_in

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

    real(dp), dimension(model%ocean_data%nbasin) :: &
         tf_basin_sum,          & ! sum of thermal_forcing over all layers in each basin
         thermal_forcing_basin, & ! basin average thermal forcing (K)
         total_bmlt_float_basin   ! sum of bmlt_float in each basin (Gt/yr)

    real(dp) :: time_from_start   ! time (yr) since the start of applying the anomaly
    real(dp) :: anomaly_fraction  ! fraction of full anomaly to apply
    real(dp) :: tf_anomaly        ! uniform thermal forcing anomaly (deg C), applied everywhere
    integer  :: tf_anomaly_basin  ! basin number where anomaly is applied;
                                  ! for default value of 0, apply to all basins
    real(dp) :: total_bmlt_float  ! global sum of bmlt_float (Gt/yr)

    real(dp) :: local_maxval, global_maxval   ! max values of a given variable
    integer :: i, j, k, nb
    integer :: ewn, nsn
    real(dp) :: dew, dns
    integer :: itest, jtest, rtest
    real(dp) :: factor

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

    if (present(calibrate_in)) then
       calibrate = calibrate_in
    else
       calibrate = .false.
    endif

    ! ------------------------------------------------------------------------
    ! Compute the basal melt rate beneath floating ice.
    ! Note: model%basal_melt is a derived type with various fields and parameters
    ! ------------------------------------------------------------------------

    !WHL - Call simple options from this subroutine instead of first calling glissade_basal_melting_float?

    if (main_task .and. verbose_bmlt_float) write(iulog,*) 'In glissade_bmlt_float_solve'

    ! Compute masks:
    ! Note: The '0.0d0' argument is thklim. Any ice with thck > 0 gets ice_mask = 1.

    call glissade_get_masks(ewn,                 nsn,                   &
                            parallel,                                   &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   0.0d0,                 &  ! thklim = 0
                            ice_mask,                                   &
                            floating_mask = floating_mask,              &
                            ocean_mask = ocean_mask,                    &
                            land_mask = land_mask)

    !TODO - Update f_ground_cell here to be on the safe side?

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

    elseif (model%options%whichbmlt_float == BMLT_FLOAT_PLUME) then

       if (this_rank == rtest .and. verbose_bmlt_float) then
          write(iulog,*) ' '
          write(iulog,*) 'Compute bmlt_float from a sub-shelf plume model'
       endif

       call glissade_plume_driver(model, model%ocean_data, model%plume)

    elseif (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then

       if (this_rank == rtest .and. verbose_bmlt_float) then
          write(iulog,*) ' '
          write(iulog,*) 'Compute bmlt_float from ocean thermal forcing'
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

       if (parallel_is_zero(model%ocean_data%thermal_forcing)) then
          call write_log('thermal forcing = 0 everywhere, GM_WARNING')
       endif

       ! Repeat for individual basins. If a basin has TF = 0 everywhere, then assume it has no valid values.
       ! Set deltaT_ocn = 0 in the basin so that we will compute bmlt_float = 0 everywhere in the basin.

       ! Sum the thermal forcing in each basin
       tf_basin_sum = 0.0d0   ! TF sum per basin
       do k = 1, model%ocean_data%nzocn
          tf_basin_sum = tf_basin_sum +  parallel_global_sum_patch(model%ocean_data%thermal_forcing(k,:,:), &
               model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       enddo

       do nb = 1, model%ocean_data%nbasin
!          if (verbose_bmlt_float .and. main_task) then
!             write(iulog,*) 'nb, TF basin sum =', nb, tf_basin_sum(nb)
!          endif
!!          if (abs(tf_basin_sum(nb)) < eps11) then
          if (abs(tf_basin_sum(nb)) < 30000.d0) then   ! temporary value to exclude basins 9, 11 and 15
             where (model%ocean_data%basin_number == nb)
                model%ocean_data%deltaT_ocn = 0.0d0
             endwhere
             model%ocean_data%deltaT_ocn_basin(nb) = 0.0d0
          endif
       enddo

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
             write(iulog,*) 'time_from_start (yr):', time_from_start
             write(iulog,*) 'ocean_data%thermal forcing anomaly  (deg):', model%ocean_data%thermal_forcing_anomaly
             write(iulog,*) 'timescale (yr):', model%ocean_data%thermal_forcing_anomaly_timescale
             write(iulog,*) 'fraction:', anomaly_fraction
             write(iulog,*) 'current TF anomaly (deg):', tf_anomaly
             if (model%ocean_data%thermal_forcing_anomaly_basin /= 0) then
                write(iulog,*) 'anomaly applied to basin number', model%ocean_data%thermal_forcing_anomaly_basin
             endif
          endif
       else
          tf_anomaly = 0.0d0
          tf_anomaly_basin = 0
       endif

       if (calibrate) then   ! one-time calibration of deltaT_ocn based on a melt rate target

          call compute_bmlt_float_thermal_forcing(&
               model%options%bmlt_float_thermal_forcing_param, &
               model%options%ocean_data_extrapolate,  &
               parallel,                              &
               ewn,                nsn,               &
               dew,                dns,               & ! m
               itest,     jtest,   rtest,             &
               ice_mask,                              &
               floating_mask,                         &
               ocean_mask,                            &
               model%geometry%marine_connection_mask, &
               model%geometry%f_ground_cell,          &
               model%geometry%thck,                   & ! m
               model%geometry%lsrf,                   & ! m
               model%geometry%topg,                   & ! m
               model%ocean_data,                      &
               model%basal_melt%thermal_forcing_mask, &
               model%basal_melt%bmlt_float,           & ! m/s
               which_ho_deltaT_ocn = model%options%which_ho_deltaT_ocn,  &  ! K
               bmlt_float_target = model%basal_melt%bmlt_float_target)      ! m/s

       else   ! usual runtime computation of bmlt_float

          call compute_bmlt_float_thermal_forcing(&
               model%options%bmlt_float_thermal_forcing_param, &
               model%options%ocean_data_extrapolate,  &
               parallel,                              &
               ewn,                nsn,               &
               dew,                dns,               & ! m
               itest,     jtest,   rtest,             &
               ice_mask,                              &
               floating_mask,                         &
               ocean_mask,                            &
               model%geometry%marine_connection_mask, &
               model%geometry%f_ground_cell,          &
               model%geometry%thck,                   & ! m
               model%geometry%lsrf,                   & ! m
               model%geometry%topg,                   & ! m
               model%ocean_data,                      &
               model%basal_melt%thermal_forcing_mask, &
               model%basal_melt%bmlt_float,           & ! m/s
               tf_anomaly_in = tf_anomaly,            & ! K
               tf_anomaly_basin_in = tf_anomaly_basin)

       endif  ! calibrate

    else  ! other options include BMLT_FLOAT_CONSTANT, BMLT_FLOAT_MISMIP, BMLT_FLOAT_DEPTH
          !TODO - Call separate subroutines for each of these options?

       call glissade_basal_melting_float(model%options%whichbmlt_float,                         &
                                         parallel,                                              &
                                         ewn,                        nsn,                       &
                                         dew,                        dns,                       &
                                         itest,                      jtest,                     &
                                         rtest,                                                 &
                                         model%general%x1,                                      & ! m
                                         model%geometry%thck,                                   & ! m
                                         model%geometry%lsrf,                                   & ! m
                                         model%geometry%topg,                                   & ! m
                                         model%climate%eus,                                     & ! m
                                         model%basal_melt,                                      & ! bmlt_float in m/s
                                         model%ocean_data)

    endif  ! whichbmlt_float


    ! If desired, add a bmlt_anomaly field.
    ! This is done for the initMIP Greenland and Antarctic experimennts.

    if (model%options%enable_bmlt_anomaly) then

       ! Add the bmlt_float anomaly where ice is present and floating
       call glissade_add_2d_anomaly(&
            model%basal_melt%bmlt_float,              &   !
            model%basal_melt%bmlt_float_anomaly,      &   !
            model%basal_melt%bmlt_anomaly_tstart,     &   ! yr
            model%basal_melt%bmlt_anomaly_timescale,  &   ! yr
            model%numerics%time)                          ! yr

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

    if (verbose_bmlt_float) then
       call point_diag(model%basal_melt%bmlt_float*scyr, 'bmlt_float after f_ground adjustment (m/yr)', &
            itest, jtest, rtest, 7, 7)
    endif

    ! Reduce basal melting in shallow cavities if bmlt_cavity_h0 > 0.
    ! The tanh function follows Asay-Davis et al. (2016), Eqs. 14 and 17.
    ! Note: model%basal_melt%bmlt_cavity_h0 has units of m.
    ! Note: For BMLT_FLOAT_MISMIP, this reduction is done in subroutine glissade_basal_melting_float
    !       based on model%basal_melt%bmlt_float_h0 and should not be repeated here.

    if (model%basal_melt%bmlt_cavity_h0 > 0.0d0 .and.  &
        model%options%whichbmlt_float /= BMLT_FLOAT_MISMIP) then

       ! TODO: Make sure lsrf is up to date. Add eus term.

       h_cavity = max(model%geometry%lsrf - model%geometry%topg, 0.0d0)  ! cavity thickness (m)

       if (verbose_bmlt_float) then
          if (this_rank == rtest) then
             write(iulog,*) 'Reduce bmlt_float in shallow cavities, bmlt_cavity_h0 (m) =', &
                  model%basal_melt%bmlt_cavity_h0
          endif
          call point_diag(h_cavity, 'h_cavity (m)', itest, jtest, rtest, 7, 7)
       endif

       where (h_cavity > 0.0d0)
          model%basal_melt%bmlt_float = model%basal_melt%bmlt_float * &
               tanh(h_cavity/model%basal_melt%bmlt_cavity_h0)
          ! WHL - Uncomment the following (and comment the line above) to replace the tanh function with a linear ramp.
!          model%basal_melt%bmlt_float = model%basal_melt%bmlt_float * &
!               min(h_cavity/model%basal_melt%bmlt_cavity_h0, 1.0d0)
       elsewhere
          model%basal_melt%bmlt_float = 0.0d0
       endwhere

    endif   ! bmlt_cavity_h0 > 0

    ! local diagnostics

    if (verbose_bmlt_float) then
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'After glissade_bmlt_float_solve, which_ho_ground_bmlt =', model%options%which_ho_ground_bmlt
       endif
       if (model%options%which_ho_ground == HO_GROUND_GLP_DELUXE) then
          call point_diag(1.0d0 - model%geometry%f_ground_cell, '1 - f_ground_cell', itest, jtest, rtest, 7, 7)
       else
          call point_diag(floating_mask, 'floating_mask', itest, jtest, rtest, 7, 7)
       endif
       call point_diag(model%basal_melt%bmlt_float*scyr, 'Final bmlt_float (m/yr)', &
            itest, jtest, rtest, 7, 7)
    endif  ! verbose_bmlt_float

    ! basin-scale diagnostics

    if (verbose_bmlt_float .and. model%ocean_data%nbasin > 1) then

       !TODO - Adjust mask to exclude ocean cells? Or does not matter?
       ! Compute the average deltaT_ocn in each basin
       call glissade_basin_average(&
            ewn,             nsn,                        &
            parallel,                                    &
            model%ocean_data%nbasin,                     &
            model%ocean_data%basin_number,               &
            model%basal_melt%thermal_forcing_mask*1.0d0, &
            model%ocean_data%deltaT_ocn,                 &
            model%ocean_data%deltaT_ocn_basin)

       !TODO - Adjust mask to exclude ocean cells? Or does not matter?
       ! Compute the average thermal forcing in each basin (including deltaT_ocn)
       call glissade_basin_average(&
            ewn,             nsn,                        &
            parallel,                                    &
            model%ocean_data%nbasin,                     &
            model%ocean_data%basin_number,               &
            model%basal_melt%thermal_forcing_mask*1.0d0, &
            model%ocean_data%thermal_forcing_lsrf + model%ocean_data%deltaT_ocn, &
            thermal_forcing_basin)

       ! Compute the average basal melt rate in each basin
       call glissade_basin_average(&
            ewn,             nsn,                        &
            parallel,                                    &
            model%ocean_data%nbasin,                     &
            model%ocean_data%basin_number,               &
            model%basal_melt%thermal_forcing_mask*1.0d0, &
            model%basal_melt%bmlt_float,                 &
            model%scalars%bmlt_float_basin)

       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'basin #, deltaT_ocn, TF (w/dT_ocn), bmlt_float (m/yr):'
          do nb = 1, model%ocean_data%nbasin
             write(iulog,'(i6,3f14.6)') nb, model%ocean_data%deltaT_ocn_basin(nb), &
                  thermal_forcing_basin(nb), model%scalars%bmlt_float_basin(nb)*scyr
          enddo
       endif

       total_bmlt_float = parallel_global_sum(model%basal_melt%bmlt_float*rhoi*dew*dns, parallel)
       total_bmlt_float_basin = parallel_global_sum_patch(model%basal_melt%bmlt_float*rhoi*dew*dns, &
                 model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       factor = (scyr/1.0d12)  ! convert kg/s to Gt/yr
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'total_bmlt_float (Gt/yr) =', total_bmlt_float*factor
          write(iulog,*) ' '
          write(iulog,*) 'basin #, basin sum:'
          do nb = 1, model%ocean_data%nbasin
             write(iulog,'(i6,2f12.6)') nb, total_bmlt_float_basin(nb)*factor
          enddo
       endif

    endif   ! verbose and nbasin > 1

  end subroutine glissade_bmlt_float_solve

!=======================================================================

  subroutine compute_bmlt_float_thermal_forcing(&
       bmlt_float_thermal_forcing_param, &
       ocean_data_extrapolate,    &
       parallel,                  &
       nx,        ny,             &
       dew,       dns,            &
       itest,     jtest,   rtest, &
       ice_mask,                  &
       floating_mask,             &
       ocean_mask,                &
       marine_connection_mask,    &
       f_ground_cell,             &
       thck,                      &
       lsrf,                      &
       topg,                      &
       ocean_data,                &
       thermal_forcing_mask,      &
       bmlt_float,                &
       tf_anomaly_in,             &
       tf_anomaly_basin_in,       &
       which_ho_deltaT_ocn,       &
       dthck_dt_obs,              &
       bmlt_float_target)

    use glimmer_paramets, only: unphys_val
    use glissade_grid_operators, only: glissade_slope_angle
    use glissade_utils, only: glissade_basin_average, glissade_interpolate_3d_ocean_field_to_lsrf

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
         floating_mask,          & !> = 1 where ice is present and floating, else = 0
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

    ! Note: thermal_forcing_mask is similar but not identical to floating_mask.
    !       * floating_mask = 1 where ice is present, and thck satisfies a flotation condition
    !       * thermal_forcing_mask = 1 where ice is present (thck > 0) and f_ground_cell < 1, with lakes excluded
    integer, dimension(nx,ny), intent(out) ::  &
         thermal_forcing_mask   !> = 1 where thermal forcing and bmlt_float can be nonzero, else = 0

    real(dp), dimension(nx,ny), intent(out) :: &
         bmlt_float             !> basal melt rate for floating ice (m/s)

    real(dp), intent(in), optional ::  &
         tf_anomaly_in      !> uniform thermal forcing anomaly (deg C), applied everywhere

    integer, intent(in), optional :: &
         tf_anomaly_basin_in    !> basin where anomaly is applied; for default value of 0, apply to all basins

    integer, intent(in), optional :: &
         which_ho_deltaT_ocn    !> option to compute deltaT_ocn; relevant here if = HO_DELTAT_OCN_DTHCK_DT

    real(dp), dimension(nx,ny), intent(in), optional :: &
         dthck_dt_obs,        & !> observed rate of thickness change (m/s)
         bmlt_float_target      !> target basal melt rate (m/s), used to calibrate deltaT_ocn

    ! local variables

    integer :: i, j, k, nb
    integer :: iglobal, jglobal
    integer :: iter

    character(len=256) :: message

    integer, dimension(nx,ny) ::  &
         new_mask                         ! temporary mask

    real(dp), dimension(ocean_data%nzocn,nx,ny) ::  &
         thermal_forcing_in               ! TF passed to subroutine interpolate_thermal_forcing_to_lsrf;
                                          ! optionally corrected for nonzero tf_anomaly
    real(dp), dimension(nx,ny) ::  &
         deltaT_ocn_init,               & ! initial value of deltaT_ocn
         theta_slope,                   & ! sub-shelf slope angle (radians)
         f_float                          ! weighting function for computing basin averages, in range [0,1]

    real(dp), dimension(ocean_data%nbasin) :: &
         thermal_forcing_basin,        &  ! basin average thermal forcing (K) at current time
         thermal_forcing_basin_old,    &  ! old value of thermal_forcing_basin
         deltaT_basin_avg                 ! basin average value of deltaT_ocn

    real(dp) :: &
         tf_anomaly                       ! local version of tf_anomaly_in

    integer ::  &
         tf_anomaly_basin                 ! local version of tf_anomaly_basin_in

    real(dp), parameter ::  &
!!         H0_float = 50.d0                 ! thickness scale (m) for floating ice; used to reduce weights when H < H0_float
         H0_float = 0.0d0                 ! thickness scale (m) for floating ice; used to reduce weights when H < H0_float

    integer, parameter :: &
         cavity_buffer = 0     ! distance from ice edge (measured in number of grid cells) over which ocean TF values
                               ! are discarded before starting TF extrapolation; must be <= nhalo
    real(dp) :: buffer
    integer :: imin, imax, jmin, jmax
    logical, save :: first_call = .true.

    if (verbose_bmlt_float .and. main_task) then
       write(iulog,*) ' '
       write(iulog,*) 'In subroutine glissade_bmlt_float_thermal_forcing'
       write(iulog,*) '   bmlt_float_thermal_forcing_param =', bmlt_float_thermal_forcing_param
       write(iulog,*) '   ocean_data_extrapolate =', ocean_data_extrapolate
       write(iulog,*) '   nbasin =', ocean_data%nbasin
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

    ! Make sure the key input fields are up to date in halos
    call parallel_halo(ice_mask, parallel)
    call parallel_halo(ocean_mask, parallel)
    call parallel_halo(marine_connection_mask, parallel)
    call parallel_halo(f_ground_cell, parallel)
    call parallel_halo(ocean_data%thermal_forcing, parallel)

    !TODO - Remove this line?
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

       if (verbose_bmlt_float) then
          if (this_rank == rtest) write(iulog,*) 'Extrapolating TF data beneath ice shelves'
          call point_diag(thermal_forcing_mask, 'thermal_forcing_mask', itest, jtest, rtest, 7, 7)
          call point_diag(marine_connection_mask, 'marine_connection_mask', itest, jtest, rtest, 7, 7)
          call point_diag(ocean_mask, 'ocean_mask', itest, jtest, rtest, 7, 7)
          call point_diag(lsrf, 'lsrf (m)', itest, jtest, rtest, 7, 7)
          call point_diag(topg, 'topg (m)', itest, jtest, rtest, 7, 7)
          do k = kmin_diag, kmax_diag
             if (this_rank == rtest) write(iulog,*) 'k =', k
             call point_diag(ocean_data%thermal_forcing(k,:,:), 'TF before extrapolating', itest, jtest, rtest, 7, 7)
          enddo
       endif

       ! Extrapolate the 3D thermal forcing field to sub-shelf cavities.
       ! Note: For now, unfilled cells retain values of unphys_val on output.

       call thermal_forcing_extrapolate(&
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

       if (verbose_bmlt_float) then
          do k = kmin_diag, kmax_diag
             if (this_rank == rtest) write(iulog,*) 'k =', k
             call point_diag(ocean_data%thermal_forcing(k,:,:), 'TF after extrapolating', &
                  itest, jtest, rtest, 7, 7)
          enddo
       endif

    else    ! ocean data are already given everywhere; do not extrapolate

       if (verbose_bmlt_float) then
          do k = kmin_diag, kmax_diag
             if (this_rank == rtest) write(iulog,*) 'k =', k
             call point_diag(ocean_data%thermal_forcing(k,:,:), 'TF to interpolate to lsrf', &
                  itest, jtest, rtest, 7, 7)
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

    call glissade_interpolate_3d_ocean_field_to_lsrf(&
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

    if (verbose_bmlt_float) then
       if (this_rank == rtest) write(iulog,*) 'basin number =', ocean_data%basin_number(itest,jtest)
       call point_diag(lsrf, 'lsrf (m)', itest, jtest, rtest, 7, 7)
       call point_diag(ocean_data%thermal_forcing_lsrf, 'thermal_forcing_lsrf (degC)', itest, jtest, rtest, 7, 7)
       if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL .or.  &
           bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or.  &
           bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then
          call point_diag(ocean_data%deltaT_ocn, 'deltaT_ocn (degC)', itest, jtest, rtest, 7, 7)
          call point_diag(ocean_data%thermal_forcing_lsrf + ocean_data%deltaT_ocn, 'corrected TF (degC)', &
               itest, jtest, rtest, 7, 7)
       endif
    endif

    ! For the ISMIP6 parameterizations, compute the average thermal forcing for the basin.
    ! Note: For the local scheme, the basin-scale thermal forcing is not used,
    !       but is computed for diagnostics.

    if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL .or.  &
        bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or.  &
        bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

       !TODO - Remove the H0_float logic, which predates the subgrid CF parameterization.
       !       For now it's simply commented out.
!       ! Compute a weighting function that is proportional to the floating fraction of ice-filled cells,
!       !  and also tapers linearly to zero for thin floating ice.
!       ! This function is used to ensure smooth changes in the basin averages as cells
!       !  transition between grounded and floating, or between ice-free and thick.
!       f_float = 1.0d0 - f_ground_cell
!       if (H0_float > 0.0d0) then
!          where (thck > 0.0d0)
!             f_float = f_float * min(thck/H0_float, 1.0d0)
!          elsewhere
!             f_float = 0.0d0
!          endwhere
!       endif

       ! Compute a weighting function that is proportional to the floating fraction of ice-filled cells
       !  and is zero for ice-free cells.
       !TODO - Does it matter if this isn't smooth? Should this be weighted by effective_areafrac?
       where (thermal_forcing_mask == 1 .and. ice_mask == 1)
          f_float = 1.0d0 - f_ground_cell
       elsewhere
          f_float = 0.0d0
       endwhere

       if (verbose_bmlt_float) then
          call point_diag(thermal_forcing_mask, 'thermal_forcing_mask', itest, jtest, rtest, 7, 7)
          call point_diag(f_float, 'f_float', itest, jtest, rtest, 7, 7)
       endif

       ! Compute the average thermal forcing for each basin.
       ! The average is taken over grid cells with thermal_forcing_mask = 1,
       !  with reduced weights for partly grounded cells and thin floating cells.
       ! Note: The basin average includes deltaT_ocn corrections.

       call glissade_basin_average(&
            nx,        ny,                   &
            parallel,                        &
            ocean_data%nbasin,               &
            ocean_data%basin_number,         &
            thermal_forcing_mask * f_float,  &   !TODO - f_float only?
            ocean_data%thermal_forcing_lsrf + ocean_data%deltaT_ocn, &
            thermal_forcing_basin)

       ! For diagnostics, compute the average value of deltaT_ocn in each basin.

       call glissade_basin_average(&
            nx,        ny,                   &
            parallel,                        &
            ocean_data%nbasin,               &
            ocean_data%basin_number,         &
            thermal_forcing_mask * f_float,  &   ! TODO - f_float only?
            ocean_data%deltaT_ocn,           &
            deltaT_basin_avg)

       if (verbose_bmlt_float .and. this_rank==rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'thermal_forcing_basin (including deltaT_ocn corrections):'
          do nb = 1, ocean_data%nbasin
             write(iulog,*) nb, thermal_forcing_basin(nb)
          enddo
          write(iulog,*) ' '
          write(iulog,*) 'deltaT_basin_avg:'
          do nb = 1, ocean_data%nbasin
             write(iulog,*) nb, deltaT_basin_avg(nb)
          enddo
       endif

       ! Compute the angle between the lower ice shelf surface and the horizontal.
       ! This option can be used to concentrate basal melting near the grounding line,
       !  where slopes are typically larger, and to reduce melting near the calving front
       !  where slopes are small.
       ! Note: The slope is currently used only for the nonlocal-slope scheme.

       call glissade_slope_angle(&
            nx,           ny,    &
            dew,          dns,   &  ! m
            itest, jtest, rtest, &
            lsrf,                &  ! m
            theta_slope,         &  ! radians
            slope_mask_in = ice_mask)

       call parallel_halo(theta_slope, parallel)

       if (verbose_bmlt_float) then
          call point_diag(sin(theta_slope), 'sin(theta_slope)', itest, jtest, rtest, 7, 7, '(f10.5)')
       endif

    endif   ! ISMIP6 melt schemes

    !-----------------------------------------------
    ! Optionally, compute deltaT_ocn to fit dthck_dt_obs
    ! Typically, this would be called only once, during the first diagnostic solve
    ! following a spin-up.
    ! The following call of ismip6_bmlt_float checks that the computation works.
    !-----------------------------------------------

    if (present(which_ho_deltaT_ocn)) then

       if (which_ho_deltaT_ocn == HO_DELTAT_OCN_CALIBRATE_BASIN) then

          ! compute deltaT_ocn_basin to match bmlt_float_target in each basin
          ! Note: After this call, the call below to ismip6_bmlt_float is redundant,
          !       but the logic is simpler if we call it anyway.

          call calibrate_deltaT_ocn_basin(&
               bmlt_float_thermal_forcing_param,     &
               nx,         ny,                       &
               dew,        dns,                      &
               itest,   jtest,    rtest,             &
               parallel,                             &
               ocean_data%nbasin,                    &
               ocean_data%basin_number,              &
               ocean_data%gamma0,                    &
               ocean_data%thermal_forcing_basin_min, &
               ocean_data%thermal_forcing_lsrf,      &
               theta_slope,                          &
               thermal_forcing_mask,                 &
               f_float,                              &
               bmlt_float_target,                    &
               thermal_forcing_basin,                &
               bmlt_float,                           &
               ocean_data%deltaT_ocn,                &
               ocean_data%deltaT_ocn_basin)

       elseif (which_ho_deltaT_ocn == HO_DELTAT_OCN_DTHCK_DT) then  ! compute deltaT_ocn to fit dthck_dt_obs

          if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL) then

             deltaT_ocn_init = ocean_data%deltaT_ocn

             call set_deltaT_ocn_for_dthck_dt(&
                  bmlt_float_thermal_forcing_param, &
                  nx,         ny,                   &
                  itest,   jtest,    rtest,         &
                  ocean_data%nbasin,                &
                  ocean_data%basin_number,          &
                  ocean_data%gamma0,                &
                  ocean_data%thermal_forcing_lsrf,  &
                  theta_slope,                      &
                  thermal_forcing_basin,            &  ! K
                  thermal_forcing_mask,             &  ! K
                  dthck_dt_obs,                     &  ! m/s
                  deltaT_ocn_init,                  &  ! K
                  ocean_data%deltaT_ocn)               ! K

          elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or. &
                  bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

             ! Iterate, recomputing the basin average thermal forcing after each iteration.
             ! Could modify to check for convergence after each iteration, but in practice,
             !  ten iterations are enough to converge to ~4 decimal places.

             deltaT_ocn_init = ocean_data%deltaT_ocn

             do iter = 1, 10

                call set_deltaT_ocn_for_dthck_dt(&
                     bmlt_float_thermal_forcing_param, &
                     nx,         ny,                   &
                     itest,   jtest,    rtest,         &
                     ocean_data%nbasin,                &
                     ocean_data%basin_number,          &
                     ocean_data%gamma0,                &
                     ocean_data%thermal_forcing_lsrf,  &
                     theta_slope,                      &
                     thermal_forcing_basin,            &  ! K
                     thermal_forcing_mask,             &  ! K
                     dthck_dt_obs,                     &  ! m/s
                     deltaT_ocn_init,                  &  ! K
                     ocean_data%deltaT_ocn)               ! K

                thermal_forcing_basin_old = thermal_forcing_basin

                call glissade_basin_average(&
                     nx,        ny,                   &
                     parallel,                        &
                     ocean_data%nbasin,               &
                     ocean_data%basin_number,         &
                     thermal_forcing_mask * f_float,  &
                     ocean_data%thermal_forcing_lsrf + ocean_data%deltaT_ocn, &
                     thermal_forcing_basin)

                ! To reduce oscillations, go halfway from the old value to the value just computed
                thermal_forcing_basin = &
                     0.5d0 * (thermal_forcing_basin_old + thermal_forcing_basin)

             enddo   ! iteration

          endif   ! bmlt_float parameterization

       endif  ! which_ho_deltaT_ocn

    endif  ! present(which_ho_deltaT_ocn)

    !-----------------------------------------------
    ! Compute the basal melt rate for each grid cell.
    !-----------------------------------------------

    if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_QUADRATIC) then

       if (verbose_bmlt_float .and. main_task) then
          write(iulog,*) 'Compute basal melt rate from quadratic parameterization'
          write(iulog,*) 'rank, i, j:', rtest, itest, jtest
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
          write(iulog,*) 'Compute basal melt rate from ISMIP6 thermal forcing'
          write(iulog,*) 'rank, i, j, basin number:', rtest, itest, jtest, ocean_data%basin_number(itest,jtest)
       endif

       ! Compute the basal melt rate based on an ISMIP6 thermal forcing parameterization.
       ! Note: bmlt_float is nonzero only for cells with thermal_forcing_mask = 1.
       call ismip6_bmlt_float(&
            bmlt_float_thermal_forcing_param,     &
            nx,                ny,                &
            itest,   jtest,    rtest,             &
            ocean_data%nbasin,                    &
            ocean_data%basin_number,              &
            ocean_data%gamma0,                    &
            ocean_data%thermal_forcing_lsrf,      &
            ocean_data%deltaT_ocn,                &
            ocean_data%thermal_forcing_basin_min, &
            theta_slope,                          &
            thermal_forcing_basin,                &
            thermal_forcing_mask,                 &
            bmlt_float)

    endif   ! bmlt_float_thermal_forcing_param

    if (verbose_bmlt_float) then
       call point_diag(bmlt_float*scyr, 'bmlt_float (m/yr)', itest, jtest, rtest, 7, 7)
    endif

    ! Reduce the melt rate in cells with thin floating ice,
    !  to reflect that these cells are only partly ice-filled.
    ! Note: This code gives bmlt_float = 0 in ice-free ocean cells,
    !       giving ice a chance to accumulate.
    !TODO - Is this needed?
    if (H0_float > 0.0d0) then
       where (f_ground_cell < 1.0d0)
          bmlt_float = bmlt_float * min(thck/H0_float, 1.0d0)
       elsewhere
          bmlt_float = 0.0d0
       endwhere
    endif

  end subroutine compute_bmlt_float_thermal_forcing

!****************************************************

  subroutine thermal_forcing_extrapolate(&
       nx,              ny,                    &
       parallel,                               &
       itest,           jtest,        rtest,   &
       nzocn,           zocn,                  &
       lsrf,            topg,                  &
       thermal_forcing_mask,                   &
       marine_connection_mask,                 &
       unphys_val,      default_val,           &
       thermal_forcing)

    use glimmer_physcon, only: dtocnfrz_dz  ! value from Beckmann & Goosse (2003), eq. 2

    ! Extrapolate ocean thermal forcing from an ocean domain that typically excludes
    !  ice shelf cavities, to an expanded domain that includes cavities.
    ! The algorithm works as follows:
    ! First, for each grid cell, determine the levels where we would like to have filled values.
    ! - In ice-free ocean cells, this includes the sea surface down to the bed topography.
    ! - In sub-shelf cavities, this includes the levels ktop:kbot that extend from just above lsrf
    !   to just below topg.
    ! Then iteratively fill the unfilled values.
    ! (1) Horizontal iteration: Loop through each level in each grid cell. For levels that are unfilled
    !      and are targeted for filling, take the average value from filled neighbor cells at the same level.
    !     Continue this horizontal averaging and filling until no  more cells can be filled.
    ! (2) In columns with at least one filled value, but other values that are not filled,
    !     fill all the unfilled values in the range ktop:kbot.
    ! Repeat (1) and (2) until no more cells and levels can be filled.
    ! (3) Final iteration: In each level that is still unfilled, average the values from adjacent cells
    !     at all levels.
    ! Note: The input thermal_forcing should be set to unphys_val for cells and levels without valid data.
    ! Note: Cells are filled only if they have a connection to the ocean through marine-based cells.
    !       Interior lakes are not filled.

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

    ! The thermal forcing field should be indexed with k = 1 at the top and k = nzocn at the bottom
    real(dp), dimension(nzocn,nx,ny), intent(inout) :: &
         thermal_forcing         ! 3D ocean thermal forcing to be extrapolated

    ! local variables

    integer, dimension(nzocn,nx,ny) :: &
         marine_connection_mask_3d, & ! = 1 where marine_conection_mask = 1 and k in range ktop:kbot, else = 0
         filled_mask,               & ! = 1 for filled cells/levels, = 0 for unfilled cells/levels
         count_smoothing            ! counts the number of times each cell or level has been smoothed

    integer, dimension(nx,ny) ::  &
         ktop,                 & ! top ocean layer in a cell to be filled, if possible
         kbot                    ! bottom ocean layer in a cell to be filled, if possible

    integer :: &
         sum_mask                ! sum of mask over neighbor cells at a given level

    real(dp) ::  &
         sum_thermal_forcing     ! sum of mask*thermal_forcing over neighbor cells at a given level

    real(dp), dimension(-1:1,-1:1) ::  &
         thermal_forcing_kmin, & ! TF in a 3x3 group of cells in the layer above
         thermal_forcing_kmax    ! TF in a 3x3 group of cells in the layer below

    integer :: &
         max_iter,             & ! max(nx,ny) * max(ewtasks, nxtasks)
         global_count,         & ! global counter for filled values
         global_count_save,    & ! global counter for filled values from previous iteration
         global_count_horiz,   & ! global counter for filled values during horizontal iteration
         global_count_horiz_save ! global counter for filled values from previous horizontal iteration

    integer :: i, j, k, ii, jj, kk
    integer :: iglobal, jglobal
    integer :: iter, iter_horiz
    integer :: imin, imax, jmin, jmax, kmin, kmax

    character(len=128) :: message

    ! Note: If thermal forcing is close to but not quite equal to unphys_val (e.g., because of roundoff error),
    !       it is interpreted as equal to unphys_val.
    real(dp), parameter :: tf_roundoff_threshold = 1.0d0  ! roundoff error threshold for thermal_forcing (deg K)

    ! Note: A stencil size of 27 (3 vertical levels, each 3 x 3) seems to give the best results with the fewest iterations
!    integer, parameter :: stencil_size = 9
!    integer, parameter :: stencil_size = 25
    integer, parameter :: stencil_size = 27

    integer, parameter :: max_count_smoothing = 2    ! how many times to apply the Laplacian smoother at a given cell and level
                                                     ! 1 => initial fill only
                                                     ! 2 => initial fill and one additional smoothing

    integer, parameter :: &
         max_iter_finish = 10    ! max iterations for the short finishing stage

    logical, parameter :: verbose_extrapolate = .false.  ! set to T to follow progress of each iteration

    ! For each marine-connected cell, compute the top and bottom layers where we need ocean data
    ! (either in the original input field, or extrapolated).
    ! Note: We define ktop and kbot >= 1 for all marine-connected cells, including 
    !       grounded marine-based cells with lsrf = topg
    !       (which might have an ocean connection via a thin water layer near the bed).

    ktop = 0
    kbot = 0
    marine_connection_mask_3d = 0

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

             marine_connection_mask_3d(ktop(i,j):kbot(i,j),i,j) = 1

          endif
       enddo   ! i
    enddo   ! j

    call parallel_halo(ktop, parallel)
    call parallel_halo(kbot, parallel)
    call parallel_halo(marine_connection_mask_3d, parallel)

    if (verbose_extrapolate) then
       call point_diag(ktop, 'ktop', itest, jtest, rtest, 7, 7)
       call point_diag(kbot, 'kbot', itest, jtest, rtest, 7, 7)
    endif

    !Michele: I have realised that the remapped ocean is taking values where it shouldn't - points that are ice
    !shelves in CISM and should not have ocean values - this could be the reason for weird patterns of remapped
    !ocean fields. So the idea could be, before starting to counting unphys_val, to set to unphys val all points
    !that are ice covered.
    !
    !WHL: For standalone CISM, I think we should do this only at startup.
    !     For coupled CISM, we should do this only at the start of a coupling time step,
    !       when there is new ocean data to extrapolate.
    !     Otherwise, we will have to do many iterations on every time step, instead of just filling a few values
    !       in newly floating cells.
    !     For this reason, I moved the code below to the higher level, before this subroutine is called.

!       do j = 1+nhalo, ny-nhalo
!          do i = 1+nhalo,  nx-nhalo
!             if (ice_mask(i,j) == 1 .or. ocean_mask(i,j) == 0) then
!                thermal_forcing(:,i,j) = unphys_val
!             end if
!          end do
!       end do

    ! Check that the stencil size is supported
    if (.not. (stencil_size==9 .or. stencil_size==25 .or. stencil_size==27)) then
       call write_log('Error: Choose a supported stencil size for ocean data smoothing (9, 25, 27)', GM_FATAL)
    endif

    ! Set TF = unphys_val outside the vertical range (ktop:kbot)
    ! This is to prevent values below the sea floor from being extrapolated into the cavity.
    do j = 1+nhalo, ny-nhalo
       do i = 1+nhalo,  nx-nhalo
          do k = 1, nzocn
             if (k < ktop(i,j) .or. k > kbot(i,j)) thermal_forcing(k,i,j) = unphys_val
          enddo
       enddo
    enddo

    ! Correct thermal_forcing for roundoff-level diffs
    where (abs(thermal_forcing - unphys_val) < tf_roundoff_threshold)
       thermal_forcing = unphys_val
    endwhere

    ! Create a mask, = 1 for filled cells/levels and = 0 for unfilled cells/levels
    where (thermal_forcing == unphys_val)
       filled_mask = 0
    elsewhere
       filled_mask = 1
    endwhere

    ! Count the number of filled cells and levels
    global_count_save = parallel_global_sum(filled_mask, parallel, marine_connection_mask_3d)
    global_count_horiz_save = global_count_save

    ! Compute the maximum number of iterations.
    ! In the worst case, the initial field is filled only at one corner of the global domain and
    !  must be extrapolated to the opposite corner, one cell at a time.
    ! Halo updates after each iteration communicate filled values to neighboring tasks.
    ! Typically for Antarctica, the number of iterations required is ~100
    !  on an 8 km grid (to reach the farthest corners of the Ross Ice Shelf), with
    !  counts doubling for each halving of grid size.

    max_iter = max(parallel%ewtasks, parallel%nstasks) * max(nx-2*nhalo, ny-2*nhalo)
    if (verbose_extrapolate .and. this_rank == rtest) then
       write(iulog,*) 'Max number of iterations =', max_iter
       write(iulog,*) 'Stencil size =', stencil_size
    endif

    ! Extrapolate the data into ice-shelf cavities
    ! Loop through all locally owned cells, filling levels ktop:kbot in unfilled cells
    !  that have one or more filled neighbors at the corresponding levels.
    ! Repeat the process until the number of filled cells and levels has converged.
    ! In the end, all cells connected to the ocean should have all levels filled between lsrf and topg.
    ! Disconnected inland lakes will not be filled.
    ! Can think of the extrapolation as a crude version of ocean circulation.

    ! Note: Typically, there are several ice sheet dynamics time steps per mass balance time step.
    !       The thermal forcing field from the ocean model (for OCEAN_DATA_GLAD) is updated
    !        at the start of the mass balance time step.  Just after this update, the extrapolation
    !        requires many iterations to converge, but subsequent updates within this
    !        mass balance time step may require < 5 iterations.

    ! Initialize the smoothing counter. Any values already filled will receive no further smoothing.
    where (filled_mask == 0)
       count_smoothing = 0
    elsewhere
       count_smoothing = max_count_smoothing
    endwhere

    ! Note: With a 1-layer stencil, many outer iterations may be required, since the horizontal fill gets stuck repeatedly.
    !       With a 3-layer stencil, only one horizontal fill is needed, so we do not need multiple outer iterations.
    !       TODO: Use a 3-layer stencil only?
    ! Note: In principle we could call glissade_laplacian_smoother.
    !       However, it's more efficient to do this inline because of the position of the k index and the detailed logic.

    do iter = 1, max_iter  ! outer loop

       if (verbose_extrapolate .and. this_rank == rtest) write(iulog,*) 'Iteration =', iter

       do iter_horiz = 1, max_iter   ! inner loop

          do j = 1+nhalo, ny-nhalo
             do i = 1+nhalo,  nx-nhalo
                if (marine_connection_mask(i,j) == 1) then   ! there is a marine path to the ocean
                   do k = ktop(i,j), kbot(i,j)
                      if (filled_mask(k,i,j) == 0 .or. count_smoothing(k,i,j) < max_count_smoothing) then

                         if (stencil_size == 9) then  ! 9-point horizontal stencil, 1 vertical level
                            ! Apply 9-point Laplacian smoother
                            imin = i-1; imax = i+1
                            jmin = j-1; jmax = j+1
                            if (any(filled_mask(k,imin:imax,jmin:jmax) == 1)) then

                               sum_mask = &
                                    1*(filled_mask(k,i-1,j-1) + filled_mask(k,i-1,j+1)  + &
                                       filled_mask(k,i+1,j-1) + filled_mask(k,i+1,j+1)) + &
                                    2*(filled_mask(k,i-1,j)   + filled_mask(k,i,j-1)    + &
                                       filled_mask(k,i+1,j)   + filled_mask(k,i,j+1))   + &
                                    4*(filled_mask(k,i,j))

                               if (sum_mask > 0) then
                                  sum_thermal_forcing = &
                                       1.d0*(filled_mask(k,i-1,j-1)*thermal_forcing(k,i-1,j-1)  + &
                                             filled_mask(k,i-1,j+1)*thermal_forcing(k,i-1,j+1)  + &
                                             filled_mask(k,i+1,j-1)*thermal_forcing(k,i+1,j-1)  + &
                                             filled_mask(k,i+1,j+1)*thermal_forcing(k,i+1,j+1)) + &
                                       2.d0*(filled_mask(k,i-1,j)*thermal_forcing(k,i-1,j)  + &
                                             filled_mask(k,i,j-1)*thermal_forcing(k,i,j-1)  + &
                                             filled_mask(k,i+1,j)*thermal_forcing(k,i+1,j)  + &
                                             filled_mask(k,i,j+1)*thermal_forcing(k,i,j+1)) + &
                                       4.d0*(filled_mask(k,i,j)*thermal_forcing(k,i,j))
                                  thermal_forcing(k,i,j) = sum_thermal_forcing / real(sum_mask,dp)
                                  count_smoothing(k,i,j) = count_smoothing(k,i,j) + 1
                               endif   ! sum_mask

                            endif   ! any are filled

                         elseif (stencil_size == 25) then  ! 5 x 5, single level
                            ! Apply 25-point Laplacian smoother
                            ! Generated the stencil by starting with a 1D 5 point stencil, row1 = (1 4 6 4 1),
                            !  and building it into a 2x2 matrix, with row5 = row1, row2 = row4 = 4*row1, and row3 = 6*row1
                            ! Note: This requires two rows of halo cells.
                            imin = i-2; imax = i+2
                            jmin = j-2; jmax = j+2
                            if (any(filled_mask(k,imin:imax,jmin:jmax) == 1)) then

                               sum_mask = &
                                    1*(filled_mask(k,i-2,j-2) + filled_mask(k,i-2,j+2)  + &
                                       filled_mask(k,i+2,j-2) + filled_mask(k,i+2,j+2)) + &
                                    4*(filled_mask(k,i-1,j-2) + filled_mask(k,i-2,j-1)  + &
                                       filled_mask(k,i-1,j+2) + filled_mask(k,i-2,j+1)  + &
                                       filled_mask(k,i+1,j-2) + filled_mask(k,i+2,j-1)  + &
                                       filled_mask(k,i+1,j+2) + filled_mask(k,i+2,j+1)) + &
                                    6*(filled_mask(k,i-2,j)   + filled_mask(k,i+2,j)    + &
                                       filled_mask(k,i,j-2)   + filled_mask(k,i,j+2))   + &
                                   16*(filled_mask(k,i-1,j-1) + filled_mask(k,i-1,j+1)  + &
                                       filled_mask(k,i+1,j-1) + filled_mask(k,i+1,j+1)) + &
                                   24*(filled_mask(k,i-1,j)   + filled_mask(k,i+1,j)    + &
                                       filled_mask(k,i,j-1)   + filled_mask(k,i,j+1))   + &
                                   36*(filled_mask(k,i,j))

                               if (sum_mask > 0) then
                                  sum_thermal_forcing = &
                                       1.d0*(filled_mask(k,i-2,j-2)*thermal_forcing(k,i-2,j-2)  + &
                                             filled_mask(k,i-2,j+2)*thermal_forcing(k,i-2,j+2)  + &
                                             filled_mask(k,i+2,j-2)*thermal_forcing(k,i+2,j-2)  + &
                                             filled_mask(k,i+2,j+2)*thermal_forcing(k,i+2,j+2)) + &
                                       4.d0*(filled_mask(k,i-1,j-2)*thermal_forcing(k,i-1,j-2)  + &
                                             filled_mask(k,i-2,j-1)*thermal_forcing(k,i-2,j-1)  + &
                                             filled_mask(k,i-1,j+2)*thermal_forcing(k,i-1,j+2)  + &
                                             filled_mask(k,i-2,j+1)*thermal_forcing(k,i-2,j+1)  + &
                                             filled_mask(k,i+1,j-2)*thermal_forcing(k,i+1,j-2)  + &
                                             filled_mask(k,i+2,j-1)*thermal_forcing(k,i+2,j-1)  + &
                                             filled_mask(k,i+1,j+2)*thermal_forcing(k,i+1,j+2)  + &
                                             filled_mask(k,i+2,j+1)*thermal_forcing(k,i+2,j+1)) + &
                                       6.d0*(filled_mask(k,i-2,j)*thermal_forcing(k,i-2,j)      + &
                                             filled_mask(k,i+2,j)*thermal_forcing(k,i+2,j)      + &
                                             filled_mask(k,i,j-2)*thermal_forcing(k,i,j-2)      + &
                                             filled_mask(k,i,j+2)*thermal_forcing(k,i,j+2))     + &
                                      16.d0*(filled_mask(k,i-1,j-1)*thermal_forcing(k,i-1,j-1)  + &
                                             filled_mask(k,i-1,j+1)*thermal_forcing(k,i-1,j+1)  + &
                                             filled_mask(k,i+1,j-1)*thermal_forcing(k,i+1,j-1)  + &
                                             filled_mask(k,i+1,j+1)*thermal_forcing(k,i+1,j+1)) + &
                                      24.d0*(filled_mask(k,i-1,j)*thermal_forcing(k,i-1,j)      + &
                                             filled_mask(k,i+1,j)*thermal_forcing(k,i+1,j)      + &
                                             filled_mask(k,i,j-1)*thermal_forcing(k,i,j-1)      + &
                                             filled_mask(k,i,j+1)*thermal_forcing(k,i,j+1))     + &
                                      36.d0*(filled_mask(k,i,j)*thermal_forcing(k,i,j))
                                  thermal_forcing(k,i,j) = sum_thermal_forcing / real(sum_mask,dp)
                                  count_smoothing(k,i,j) = count_smoothing(k,i,j) + 1
                               endif   ! sum_mask

                            endif   ! any are filled

                         elseif (stencil_size == 27) then  ! 9-point horizontal stencil, 3 vertical levels
                            ! Apply 27-point Laplacian smoother
                            ! Layers above and below use the 9-point stencil above
                            ! Central layer uses the same stencil, multiplied by 2
                            imin = i-1; imax = i+1
                            jmin = j-1; jmax = j+1
                            kmin = max(k-1,ktop(i,j))  ! layer above, if potentially filled
                            kmax = min(k+1,kbot(i,j))  ! layer below, if potentially filled
                            if (any(filled_mask(kmin:kmax,imin:imax,jmin:jmax) == 1)) then

                               sum_mask = &
                                    2*(filled_mask(k,i-1,j-1) + filled_mask(k,i-1,j+1)  + &
                                       filled_mask(k,i+1,j-1) + filled_mask(k,i+1,j+1)) + &
                                    4*(filled_mask(k,i-1,j)   + filled_mask(k,i,j-1)    + &
                                       filled_mask(k,i+1,j)   + filled_mask(k,i,j+1))   + &
                                    8*(filled_mask(k,i,j))
                               if (kmin < k) then
                                  sum_mask = sum_mask + &
                                       1*(filled_mask(kmin,i-1,j-1) + filled_mask(kmin,i-1,j+1)  + &
                                          filled_mask(kmin,i+1,j-1) + filled_mask(kmin,i+1,j+1)) + &
                                       2*(filled_mask(kmin,i-1,j)   + filled_mask(kmin,i,j-1)    + &
                                          filled_mask(kmin,i+1,j)   + filled_mask(kmin,i,j+1))   + &
                                       4*(filled_mask(kmin,i,j))
                               endif
                               if (kmax > k) then
                                  sum_mask = sum_mask + &
                                       1*(filled_mask(kmax,i-1,j-1) + filled_mask(kmax,i-1,j+1)  + &
                                          filled_mask(kmax,i+1,j-1) + filled_mask(kmax,i+1,j+1)) + &
                                       2*(filled_mask(kmax,i-1,j)   + filled_mask(kmax,i,j-1)    + &
                                          filled_mask(kmax,i+1,j)   + filled_mask(kmax,i,j+1))   + &
                                       4*(filled_mask(kmax,i,j))
                               endif

                               if (sum_mask > 0) then
                                  sum_thermal_forcing = &
                                       2.d0*(filled_mask(k,i-1,j-1)*thermal_forcing(k,i-1,j-1)  + &
                                             filled_mask(k,i-1,j+1)*thermal_forcing(k,i-1,j+1)  + &
                                             filled_mask(k,i+1,j-1)*thermal_forcing(k,i+1,j-1)  + &
                                             filled_mask(k,i+1,j+1)*thermal_forcing(k,i+1,j+1)) + &
                                       4.d0*(filled_mask(k,i-1,j)*thermal_forcing(k,i-1,j)  + &
                                             filled_mask(k,i,j-1)*thermal_forcing(k,i,j-1)  + &
                                             filled_mask(k,i+1,j)*thermal_forcing(k,i+1,j)  + &
                                             filled_mask(k,i,j+1)*thermal_forcing(k,i,j+1)) + &
                                       8.d0*(filled_mask(k,i,j)*thermal_forcing(k,i,j))
                                  if (kmin < k) then
                                     ! Compute TF of water in the layer above, if that water were moved down to this layer
                                     thermal_forcing_kmin(:,:) = thermal_forcing(kmin,i-1:i+1,j-1:j+1) &
                                          - dtocnfrz_dz * (zocn(k) - zocn(kmin))
                                     sum_thermal_forcing = sum_thermal_forcing + &
                                          1.d0*(filled_mask(kmin,i-1,j-1)*thermal_forcing_kmin(-1,-1)  + &
                                                filled_mask(kmin,i-1,j+1)*thermal_forcing_kmin(-1, 1)  + &
                                                filled_mask(kmin,i+1,j-1)*thermal_forcing_kmin( 1,-1)  + &
                                                filled_mask(kmin,i+1,j+1)*thermal_forcing_kmin( 1, 1)) + &
                                          2.d0*(filled_mask(kmin,i-1,j)*thermal_forcing_kmin(-1, 0)  + &
                                                filled_mask(kmin,i,j-1)*thermal_forcing_kmin( 0,-1)  + &
                                                filled_mask(kmin,i+1,j)*thermal_forcing_kmin( 1, 0)  + &
                                                filled_mask(kmin,i,j+1)*thermal_forcing_kmin( 0, 1)) + &
                                          4.d0*(filled_mask(kmin,i,j)*thermal_forcing_kmin(0,0))
                                  endif
                                  if (kmax > k) then
                                     ! Compute TF of water in the layer below, if that water were moved up to this layer
                                     thermal_forcing_kmax(:,:) = thermal_forcing(kmax,i-1:i+1,j-1:j+1) &
                                          - dtocnfrz_dz * (zocn(k) - zocn(kmax))
                                     sum_thermal_forcing = sum_thermal_forcing + &
                                          1.d0*(filled_mask(kmax,i-1,j-1)*thermal_forcing_kmax(-1,-1)  + &
                                                filled_mask(kmax,i-1,j+1)*thermal_forcing_kmax(-1, 1)  + &
                                                filled_mask(kmax,i+1,j-1)*thermal_forcing_kmax( 1,-1)  + &
                                                filled_mask(kmax,i+1,j+1)*thermal_forcing_kmax( 1, 1)) + &
                                          2.d0*(filled_mask(kmax,i-1,j)*thermal_forcing_kmax(-1, 0)  + &
                                                filled_mask(kmax,i,j-1)*thermal_forcing_kmax( 0,-1)  + &
                                                filled_mask(kmax,i+1,j)*thermal_forcing_kmax( 1, 0)  + &
                                                filled_mask(kmax,i,j+1)*thermal_forcing_kmax( 0, 1)) + &
                                          4.d0*(filled_mask(kmax,i,j)*thermal_forcing_kmax(0,0))
                                  endif
                                  thermal_forcing(k,i,j) = sum_thermal_forcing / real(sum_mask,dp)
                                  count_smoothing(k,i,j) = count_smoothing(k,i,j) + 1
                               endif   ! sum_mask

                            endif   ! any are filled

                         endif   ! stencil_size
                      endif   ! not filled or fully smoothed
                   enddo   ! k
                endif   ! marine_connection mask
             enddo   ! i
          enddo  ! j

          call parallel_halo(thermal_forcing, parallel)

          ! update the mask
          where (thermal_forcing == unphys_val)
             filled_mask = 0
          elsewhere
             filled_mask = 1
          endwhere

          ! Every few iterations, check for convergence of the horizontal iteration loop.
          ! Check after iterations 1, 2, and 5, then after 10, 20, etc.
          if (iter_horiz == 1 .or. iter_horiz == 2 .or. iter_horiz == 5 .or. mod(iter_horiz, 10) == 0) then

             if (verbose_extrapolate) then
                do k = kmin_diag, kmax_diag
                   if (this_rank == rtest) write(iulog,*) 'k, zocn =', k, zocn(k)
                   call point_diag(filled_mask(k,:,:), 'filled_mask', itest, jtest, rtest, 7, 7)
                   call point_diag(thermal_forcing(k,:,:), 'thermal_forcing', itest, jtest, rtest, 7, 7)
                enddo
             endif

             global_count_horiz = parallel_global_sum(filled_mask, parallel, marine_connection_mask_3d)

             if (global_count_horiz == global_count_horiz_save) then
                if (verbose_extrapolate) then
                   if (this_rank == rtest) then
                      write(iulog,*) 'Horizontal extrapolation converged: iter, global_count =', &
                        iter_horiz, global_count_horiz
                   endif
                   write(message,*) 'Horizontal extrapolation converged: iter, global_count =', &
                        iter_horiz, global_count_horiz
                   call write_log(message)
                endif
                exit
             else   ! not converged
                if (verbose_extrapolate) then
                   if (this_rank == rtest) then
                      write(iulog,*) 'Horizontal extrapolation convergence check: iter, global_count =', &
                           iter_horiz, global_count_horiz
                   endif
                   write(message,*) 'Horizontal extrapolation convergence check: iter, global_count =', &
                        iter_horiz, global_count_horiz
                   call write_log(message)
                endif
                global_count_horiz_save = global_count_horiz
             endif

          endif   ! time for a convergence check

          if (iter_horiz == max_iter) then
             if (this_rank == rtest) write(iulog,*) 'iter_horiz = max_iter:', max_iter
             write(message,*) 'Ocean TF extrapolation error, too many iterations: max_iter =', max_iter
             call write_log(message, GM_FATAL)
          endif

       enddo ! do iter_horiz

       ! Note: If using a 3-level Laplacian stencil, the following downward and upward fills are unnecessary,
       !       and the outer iteration loop is not needed.

       ! Extend TF downward in columns that contain unfilled cells below filled cells.
       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo,  nx-nhalo
             if (marine_connection_mask(i,j) == 1) then
                do k = ktop(i,j)+1, kbot(i,j)
                   if (thermal_forcing(k,i,j) == unphys_val .and. thermal_forcing(k-1,i,j) /= unphys_val) then
                      thermal_forcing(k,i,j) = thermal_forcing(k-1,i,j) - dtocnfrz_dz * (zocn(k) - zocn(k-1))
                   endif
                enddo   ! k
             endif   ! marine_connection_mask
          enddo   ! i
       enddo  ! j

       ! Extend TF upward in columns that contain unfilled cells above filled cells.

       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo,  nx-nhalo
             if (marine_connection_mask(i,j) == 1) then
                do k = kbot(i,j)-1, ktop(i,j), -1
                   if (thermal_forcing(k,i,j) == unphys_val .and. thermal_forcing(k+1,i,j) /= unphys_val) then
                      thermal_forcing(k,i,j) = thermal_forcing(k+1,i,j) - dtocnfrz_dz * (zocn(k) - zocn(k+1))
                   endif
                enddo   ! k
             endif   ! marine_connection_mask
          enddo   ! i
       enddo  ! j

       call parallel_halo(thermal_forcing, parallel)

       if (verbose_extrapolate) then
          if (this_rank == rtest) write(iulog,*) 'After vertical fill:'
          do k = kmin_diag, kmax_diag
             if (this_rank == rtest) write(iulog,*) 'k, zocn =', k, zocn(k)
             call point_diag(filled_mask(k,:,:), 'filled_mask', itest, jtest, rtest, 7, 7)
             call point_diag(thermal_forcing(k,:,:), 'thermal_forcing', itest, jtest, rtest, 7, 7)
          enddo
       endif

       ! update the mask
       where (thermal_forcing == unphys_val)
          filled_mask = 0
       elsewhere
          filled_mask = 1
       endwhere

       ! Every few iterations, check for convergence of the outer iteration loop.
       ! Check after iterations 1, 2, and 5, then after 10, 15, etc.
       if (iter == 1 .or. iter == 2 .or. mod(iter, 5) == 0) then

          global_count = parallel_global_sum(filled_mask, parallel, marine_connection_mask_3d)

          if (global_count == global_count_save) then
             if (verbose_extrapolate) then
                if (this_rank == rtest) &
                  write(iulog,*) 'Extrapolation converged: iter, global_count =', iter, global_count
                write(message,*) 'Extrapolation converged: iter, global_count =', iter, global_count
                call write_log(message)
             endif
             exit
          else
             if (verbose_extrapolate) then
                if (this_rank == rtest) &
                  write(iulog,*) 'Extrapolation convergence check: iter, global_count =', iter, global_count
                write(message,*) 'Extrapolation convergence check: iter, global_count =', iter, global_count
                call write_log(message)
             endif
             global_count_save = global_count
          endif
          
       endif   ! time for a convergence check

       if (iter == max_iter) then
          if (this_rank == rtest) write(iulog,*) 'iter = max_iter:', max_iter
          write(message,*) 'Ocean TF extrapolation error, too many iterations: max_iter =', max_iter
          call write_log(message, GM_FATAL)
       endif

    enddo   ! max_iter

    ! Make sure all levels from ktop to kbot are filled in cells with thermal_forcing_mask = 1.
    ! If not yet filled, then check for a value in a nearby cell.
    ! If there is no such value, then there should be in a subsequent iteration.

    if (verbose_extrapolate .and. this_rank == rtest) write(iulog,*) 'Final extrapolation to outlier cells'

    do iter = 1, max_iter_finish

       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo, nx-nhalo
             if (thermal_forcing_mask(i,j) == 1) then
                do k = ktop(i,j), kbot(i,j)
                   if (thermal_forcing(k,i,j) == unphys_val) then
                      call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                      if (verbose_extrapolate) then
                         write(iulog,*) 'Ocean data extrapolation issue: unphys value in level k, i, j:', &
                              k, iglobal, jglobal, thermal_forcing(k,i,j)
                      endif
                      !Note: no risk of i or j out of bounds, since the loop above is limited to locally owned cells
                      ! Average valid values from all levels of adjacent columns
                      ! Include a depth adjustment for values from different levels
                      sum_mask = 0
                      sum_thermal_forcing = 0.0d0
                      do jj = j-1, j+1
                         do ii = i-1, i+1
                            do kk = 1,nzocn
                               if (filled_mask(kk,ii,jj) == 1) then
                                  sum_mask = sum_mask + 1
                                  sum_thermal_forcing = sum_thermal_forcing &
                                       + thermal_forcing(kk,ii,jj) - dtocnfrz_dz*(zocn(k) - zocn(kk))
                               endif
                            enddo
                         enddo
                      enddo
                      if(sum_mask > 0) then
                         thermal_forcing(k,i,j) = sum_thermal_forcing / real(sum_mask,dp)
                         if (verbose_extrapolate) write(iulog,*) '   Assigned value is ', thermal_forcing(k,i,j)
                      endif
                   endif   ! unphys_val
                enddo   ! k
             endif   ! thermal_forcing_mask
          enddo   ! i
       enddo   ! j

       call parallel_halo(thermal_forcing, parallel)

       where (thermal_forcing == unphys_val)
          filled_mask = 0
       elsewhere
          filled_mask = 1
       endwhere

       ! check for convergence

       global_count = parallel_global_sum(filled_mask, parallel, marine_connection_mask_3d)

       if (global_count == global_count_save) then
          if (verbose_extrapolate) then
             if (this_rank == rtest) &
                  write(iulog,*) 'Final extrapolation converged: iter, global_count =', iter, global_count
             write(message,*) 'Final extrapolation converged: iter, global_count =', iter, global_count
             call write_log(message)
          endif
          exit
       else
          if (verbose_extrapolate) then
             if (this_rank == rtest) &
                  write(iulog,*) 'Final extrapolation convergence check: iter, global_count =', iter, global_count
             write(message,*) 'Final extrapolation convergence check: iter, global_count =', iter, global_count
             call write_log(message)
          endif
          global_count_save = global_count
       endif

    enddo   ! max_iter_finish

    ! Bug check - Make sure there are valid values wherever thermal_forcing_mask = 1

    do j = 1+nhalo, ny-nhalo
       do i = 1+nhalo, nx-nhalo
          if (thermal_forcing_mask(i,j) == 1) then
             do k = ktop(i,j), kbot(i,j)
                if (thermal_forcing(k,i,j) == unphys_val) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   write(iulog,*) 'TF extrapolation error, no neighbours with valid values, k, i, j =',  &
                        k, iglobal, jglobal
                   write(message,*) 'TF extrapolation error, no neighbours with valid values, k, i, j =',  &
                        k, iglobal, jglobal
                   call write_log(message, GM_FATAL)
                endif
             enddo
          endif
       enddo
    enddo

    ! Set the thermal forcing to a default value in the remaining unfilled cell/levels.
    ! For example, we might set thermal_forcing = 0 to get cleaner diagnostics.

    if (default_val /= unphys_val) then
       where (thermal_forcing == unphys_val)
          thermal_forcing = default_val
       endwhere
    endif

  end subroutine thermal_forcing_extrapolate

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
       thermal_forcing_basin_min, &
       theta_slope,               &
       thermal_forcing_basin,     &
       thermal_forcing_mask,      &
       bmlt_float)

    ! Compute the basal melt rate as a quadratic function of thermal forcing, using either
    ! a local or nonlocal parameterization as specified for ISMIP6.
    ! Note: The input thermal forcing fields are nonnegative, but the deltaT correction can be negative
    !        giving negative effective TF.
    !       For the local scheme, negative effective TF is not allowed.
    !       For the nonlocal schemes, negative effective TF locally can be multiplied by non-negative
    !        basin-scale TF, giving basal freezing.

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
         theta_slope              !> sub-shelf slope angle (radians)

    real(dp), dimension(nx,ny), intent(in) :: &
         deltaT_ocn               !> thermal forcing correction factor (deg C)

    ! Note: If thermal_forcing_basin_min > 0, there will always be some thermal forcing in the basin
    !        even if the basin is too cold to melt ice. This will promote basal freezing in cells with local TF < 0.
    !       In the ISMIP6 dataset, the Ronne cavity has TF ~ 0.5 C near the GL; Filchner has TF ~0.3 C.
    real(dp), intent(in) :: &
         thermal_forcing_basin_min  !> min basin-scale TF; can be applied to nonlocal and nonlocal-slope schemes

    real(dp), dimension(nbasin), intent(in) :: &
         thermal_forcing_basin    !> thermal forcing averaged over each basin (deg C)

    integer, dimension(nx,ny), intent(in) :: &
         thermal_forcing_mask     !> = 1 where TF-driven bmlt_float can be > 0

    real(dp), dimension(nx,ny), intent(out) :: &
         bmlt_float               !> basal melt rate (m/s) at lower ice surface

    ! local variables

    integer :: i, j, nb

    real(dp) :: coeff         ! constant coefficient = [(rhoo*cpw)/(rhoi*lhci)]^2, with units deg^(-2)

    real(dp) :: &
         eff_thermal_forcing,      & ! effective local thermal forcing, after deltaT correction
         eff_thermal_forcing_basin   ! effective basin thermal forcing, after deltaT correction


    ! initialize
    bmlt_float(:,:) = 0.0d0

    ! Note: gamma0 has units of m/yr; divide by scyr so bmlt_float has units of m/s
    coeff = gamma0/scyr * ( (rhoo*cpw)/(rhoi*lhci) )**2

    if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL) then

       ! local parameterization
       ! melt rate is a quadratic function of local thermal forcing

       do j = 1, ny
          do i = 1, nx
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing = max(0.0d0, thermal_forcing_lsrf(i,j) + deltaT_ocn(i,j))
                bmlt_float(i,j) = coeff * eff_thermal_forcing**2
             endif
          enddo
       enddo

    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL) then

       ! nonlocal parameterization
       ! melt rate is a quadratic function of local thermal forcing and basin-average thermal forcing

       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn(i,j)
                eff_thermal_forcing_basin = max(thermal_forcing_basin_min, thermal_forcing_basin(nb))
                bmlt_float(i,j) = coeff * eff_thermal_forcing * eff_thermal_forcing_basin
             endif
          enddo
       enddo

    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

       ! same as nonlocal, but with larger gamma0, and multiplied by sin(theta_slope)

       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn(i,j)
                eff_thermal_forcing_basin = max(thermal_forcing_basin_min, thermal_forcing_basin(nb))
                bmlt_float(i,j) = coeff * sin(theta_slope(i,j)) * eff_thermal_forcing * eff_thermal_forcing_basin
             endif
          enddo
       enddo

    endif  ! local or nonlocal

  end subroutine ismip6_bmlt_float

!****************************************************

  subroutine calibrate_deltaT_ocn_basin(&
       bmlt_float_thermal_forcing_param, &
       nx,         ny,            &
       dew,        dns,           &
       itest,   jtest,    rtest,  &
       parallel,                  &
       nbasin,                    &
       basin_number,              &
       gamma0,                    &
       thermal_forcing_basin_min, &
       thermal_forcing_lsrf,      &
       theta_slope,               &
       thermal_forcing_mask,      &
       f_float,                   &
       bmlt_float_target,         &
       thermal_forcing_basin,     &
       bmlt_float,                &
       deltaT_ocn,                &
       deltaT_ocn_basin)

    ! This subroutine adjusts deltaT_basin to minimize the difference between the model and target melt rates
    !  in each basin. Typically the target melt rate is based on observations.
    ! The melt rate is averaged over cells that (1) are floating in CISM and (2) have valid melt rate targets.
    ! Negative melt rates (i.e., freeze-on) are allowed.
    ! Note: This subroutine is called only at initialization.
    !       If the run continues after initialization, it will use the deltaT_ocn values computed here.

    use glissade_utils, only: glissade_basin_average

    integer, intent(in) :: &
         bmlt_float_thermal_forcing_param  !> kind of melting parameterization, local or nonlocal

    integer, intent(in) :: &
         nx, ny                   !> number of grid cells in each dimension

    real(dp), intent(in) :: &
         dew, dns                 !> grid cell size (m)

    integer, intent(in) :: &
         itest, jtest, rtest      !> coordinates of diagnostic point

    type(parallel_type), intent(in) :: &
         parallel                 !> info for parallel communication

    integer, intent(in) :: &
         nbasin                   !> number of basins

    integer, dimension(nx,ny), intent(in) :: &
         basin_number             !> integer ID for each basin

    real(dp), intent(in) :: &
         gamma0                   !> basal melt rate coefficient (m/yr)

    real(dp), intent(in) :: &
         thermal_forcing_basin_min  !> min basin-scale TF; can be applied to nonlocal and nonlocal-slope schemes
                                    !> default = 0.; this means the basin-scale thermal forcing can vanish in cold basins

    integer, dimension(nx,ny), intent(in) :: &
         thermal_forcing_mask     !> = 1 where TF-driven bmlt_float can be > 0


    real(dp), dimension(nx,ny), intent(in) :: &
         f_float,               & ! weighting function for computing basin averages, in range [0,1];
                                  ! = (1 - f_ground_cell) for ice-filled cells, else = 0
         thermal_forcing_lsrf,  & !> thermal forcing (K) at lower ice surface
         theta_slope,           & !> sub-shelf slope angle (radians)
         bmlt_float_target        !> target value of basal melt rate (m/s)

    real(dp), dimension(nx,ny), intent(out) :: &
         bmlt_float,            & !> basal melt rate (m/s)
         deltaT_ocn               !> 2D thermal forcing correction (deg K); uniform within each basin

    real(dp), dimension(nbasin), intent(out) :: &
         thermal_forcing_basin, & ! basin-average thermal forcing (K)
         deltaT_ocn_basin         !> basin-scale thermal forcing correction (deg K);
                                  !> minimizes difference between computed and target bmlt_float

    ! local variables

    integer :: i, j, nb, iter

    real(dp), dimension(nx,ny) :: &
         bmlt_float_mask               ! = 1 where bmlt_float and bmlt_float_target have valid values

    real(dp), dimension(nbasin) :: &
         bmlt_float_target_basin,    & ! target for basin-average melt rate (m/s)
         bmlt_float_basin,           & ! basin-average melt rate (m/s)
         dT_basin_increment            ! incremental change in deltaT_ocn_basin

    real(dp) :: &
         coeff,                      & ! coefficient in basal melt formulas
         eff_thermal_forcing_basin,  & ! basin-avearge effective thermal forcing (K)
         numer, denom,               & !
         diff, max_diff,             & ! difference between model and target melt rates (m/yr)
         factor,                     & ! unit conversion factor
         total_bmlt_float,           & ! global sum of bmlt_float (kg/s)
         total_bmlt_float_target       ! global sum of bmlt_float_target (kg/s)

    real(dp), parameter :: &
         bmlt_float_threshold = 1.e-4  ! threshold for melt rate convergence (m/yr)

    integer, parameter :: max_iter = 25

    if (verbose_bmlt_float .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'In calibrate deltaT_ocn_basin, ISMIP6 TF param =', &
            bmlt_float_thermal_forcing_param
    endif

    ! initialize

    coeff = (gamma0/scyr) * ( (rhoo*cpw)/(rhoi*lhci) )**2
    deltaT_ocn = 0.0d0
    deltaT_ocn_basin = 0.0d0
    bmlt_float = 0.0d0

    ! Compute a real-valued mask of cells with valid values of bmlt_float and bmlt_float_target:
    ! bmlt_float_mask = f_float (the floating fraction in cells where ice is present)
    !  provided bmlt_float_target /= 0, else bmlt_float_mask = 0.
    ! The goal of the calculation below is to compute deltaT_ocn such that
    !  the basin-average values of bmlt_float and bmlt_float_target are equal,
    !  where the average is taken over this mask.

    where (bmlt_float_target /= 0.0d0)
       bmlt_float_mask = f_float
    elsewhere
       bmlt_float_mask = 0.0d0
    endwhere

    ! Compute the target melt rate for each basin
    ! Note: The mask is over the shared region.

    call glissade_basin_average(&
         nx,             ny,        &
         parallel,                  &
         nbasin,                    &
         basin_number,              &
         bmlt_float_mask,           &
         bmlt_float_target,         &
         bmlt_float_target_basin)

    ! iterative loop

    do iter = 1, max_iter

       ! Compute the average thermal forcing for each basin, given the latest deltaT_ocn.

       call glissade_basin_average(&
            nx,             ny,                 &
            parallel,                           &
            nbasin,                             &
            basin_number,                       &
            bmlt_float_mask,                    &
            thermal_forcing_lsrf + deltaT_ocn,  &
            thermal_forcing_basin)

       ! Update bmlt_float, depending on the melt parameterization

       call ismip6_bmlt_float(&
            bmlt_float_thermal_forcing_param, &
            nx,         ny,            &
            itest,   jtest,    rtest,  &
            nbasin,                    &
            basin_number,              &
            gamma0,                    &
            thermal_forcing_lsrf,      &
            deltaT_ocn,                &
            thermal_forcing_basin_min, &
            theta_slope,               &
            thermal_forcing_basin,     &
            thermal_forcing_mask,      &
            bmlt_float)

       ! Update the average melt rate in each basin

       call glissade_basin_average(&
            nx,             ny,        &
            parallel,                  &
            nbasin,                    &
            basin_number,              &
            bmlt_float_mask,           &
            bmlt_float,                &
            bmlt_float_basin)

       ! Adjust deltaT_ocn_basin.
       ! The goal is to obtain a basin-average melt rate that matches the basin-average target rate.

       do nb = 1, nbasin

          ! Note: The factor of 2 in the denom comes from differentiating the quadratic term
          numer = bmlt_float_target_basin(nb) - bmlt_float_basin(nb)
          denom = 2.0d0 * coeff * thermal_forcing_basin(nb)

          ! Some trial and error showed that multiplying (numer/denom) by 0.9
          ! leads to faster convergence than taking the full increment;
          ! converges in ~7 iterations with a threshold of 1.e-4 m/yr
          if (denom > 0.0d0) then
             dT_basin_increment(nb) = 0.9d0 * numer/denom
          else
             dT_basin_increment(nb) = 0.0d0
          endif
          deltaT_ocn_basin(nb) = deltaT_ocn_basin(nb) + dT_basin_increment(nb)

       enddo   ! nbasin

       ! Copy the basin values to the 2D deltaT_ocn array
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             nb = basin_number(i,j)
             deltaT_ocn(i,j) = deltaT_ocn_basin(nb)
          enddo
       enddo

       ! Check for convergence
       ! Note: Convert m/s to m/yr when computing differences
       max_diff = 0.0d0
       do nb = 1, nbasin
          diff = abs(bmlt_float_target_basin(nb) - bmlt_float_basin(nb))*scyr
          max_diff = max(diff, max_diff)
       enddo
       if (verbose_bmlt_float .and. this_rank == rtest) then
          write(iulog,*) 'iter, max diff b/w model and target bmlt_float (m/yr):', iter, max_diff
       endif

       if (max_diff < bmlt_float_threshold) then
          if (verbose_bmlt_float .and. this_rank == rtest) write(iulog,*) 'Iteration converged'
          exit
       elseif (iter == max_iter) then
          call write_log('Exceeded max # of iterations for deltaT_ocn calibration', GM_FATAL)
       endif

    enddo   ! iter

    call parallel_halo(deltaT_ocn, parallel)

    ! Optional global and basin-scale diagnostics

    if (verbose_bmlt_float) then
       ! Note: bmlt_float_basin and bmlt_float_target_basin should agree to within bmlt_float_threshold
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'basin #, TF_basin, dT_ocn_basin, bmlt_float, bmlt_float_tgt:'
          do nb = 1, nbasin
             write(iulog,'(i6,4f12.6)') nb, thermal_forcing_basin(nb), deltaT_ocn_basin(nb), &
                  bmlt_float_basin(nb)*scyr, bmlt_float_target_basin(nb)*scyr
          enddo
       endif

       total_bmlt_float = parallel_global_sum(bmlt_float*f_float*rhoi*dew*dns, parallel)
       total_bmlt_float_target = parallel_global_sum(bmlt_float_target*f_float*rhoi*dew*dns, parallel)
       bmlt_float_basin = parallel_global_sum_patch(bmlt_float*f_float*rhoi*dew*dns, &
                 nbasin, basin_number, parallel)
       bmlt_float_target_basin = parallel_global_sum_patch(bmlt_float_target*f_float*rhoi*dew*dns, &
                 nbasin, basin_number, parallel)
       factor = (scyr/1.0d12)  ! convert kg/s to Gt/yr
       ! Note: bmlt_float_basin and bmlt_float_target_basin will generally not agree.
       !       This is because the bmlt_float sums include regions where bmlt_float_target = 0.
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'total_bmlt_float (Gt/yr) =', total_bmlt_float*factor
          write(iulog,*) 'total_bmlt_float_target (Gt/yr) =', total_bmlt_float_target*factor
          write(iulog,*) ' '
          write(iulog,*) 'basin #, basin sum, target sum:'
          do nb = 1, nbasin
             write(iulog,'(i6,2f12.6)') nb, bmlt_float_basin(nb)*factor, bmlt_float_target_basin(nb)*factor
          enddo
       endif
    endif

  end subroutine calibrate_deltaT_ocn_basin

!****************************************************

  subroutine set_deltaT_ocn_for_dthck_dt(&
       bmlt_float_thermal_forcing_param, &
       nx,         ny,            &
       itest,   jtest,    rtest,  &
       nbasin,                    &
       basin_number,              &
       gamma0,                    &
       thermal_forcing_lsrf,      &
       theta_slope,               &
       thermal_forcing_basin,     &
       thermal_forcing_mask,      &
       dthck_dt_target,           &
       deltaT_ocn_init,           &
       deltaT_ocn_new)

    ! This subroutine adjusts deltaT_ocn to match a target basal melt rate = -dH/dt.
    ! Typically the target rate comes from observations. Where the target melt rate < 0 (i.e., dH/dt > 0),
    !  no correction is computed.
    ! It is assumed that dH/dt = 0 for deltaT_ocn_init.
    !
    ! The adjustment is made for all cells that can have melting driven by thermal forcing (thermal_forcing_mask = 1).
    ! This includes cells that are currently grounded, but might be floating in a forward run.

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

    real(dp), dimension(nbasin), intent(in) :: &
         thermal_forcing_basin    !> thermal forcing averaged over each basin (deg C)

    integer, dimension(nx,ny), intent(in) :: &
         thermal_forcing_mask     !> = 1 for cells with melting driven by thermal forcing, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         thermal_forcing_lsrf,  & !> thermal forcing (K) at lower ice surface
         theta_slope,           & !> sub-shelf slope angle (radians)
         dthck_dt_target,       & !> target value of dthck_dt (m/s)
         deltaT_ocn_init          !> initial thermal forcing correction factor (deg C)

    real(dp), dimension(nx,ny), intent(out) :: &
         deltaT_ocn_new           !> new thermal forcing correction factor (deg C)

    ! local variables

    integer :: i, j, nb

    real(dp) :: coeff         ! constant coefficient = [(rhow*cp)/(rhoi*Lf)]^2, with units deg^(-2)

    real(dp), dimension(nx,ny) :: &
         bmlt_float_init,            & ! initial melt rate (m/s)  before adding dTocn
         bmlt_float_new,             & ! new melt rate (m/s) after adding dTocn
         dbmlt_float,                & ! additional melting needed (m/s)
         dTocn                         ! ocean warming term (deg C), added to deltaT_ocn_init

    real(dp) :: &
         eff_thermal_forcing,        & ! effective local thermal forcing (deg C), before adding dTocn
         eff_thermal_forcing_basin     ! effective basin thermal forcing (deg C), before adding dTocn

    ! initialize

    ! Note: gamma0 has units of m/yr; dividing by scyr gives melt rates in m/s
    coeff = (gamma0/scyr) * ( (rhoo*cpw)/(rhoi*lhci) )**2
    dTocn = 0.0d0
    bmlt_float_init = 0.0d0
    bmlt_float_new = 0.0d0

    ! Ignore cells where dthck_dt_target > 0
    dbmlt_float = max(-1.d0*dthck_dt_target, 0.0d0)

    if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL) then

       ! local parameterization
       ! melt rate is a quadratic function of local thermal forcing:
       !    m = C * max(0, F0)^2,
       ! where C = coeff and F0 = initial thermal forcing
       ! If F0 > 0, the increase in melt rate due to ocean warming dT is given by
       !    dm = 2C * F0 * dT + C * dT^2,
       ! a quadratic equation that we solve for dT:
       !    dT = -F0 + sqrt(F0^2 + dm/C).
       ! If F0 < 0, then we solve
       !    dm = C * (F0 + dT)^2, giving
       !    dT = sqrt(dm/C) - F0.

       do j = 1, ny
          do i = 1, nx
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn_init(i,j)
                if (eff_thermal_forcing > 0.0d0) then  ! increment positive TF
                   bmlt_float_init(i,j) = coeff * eff_thermal_forcing**2
                   dTocn(i,j) = -eff_thermal_forcing   &
                        + sqrt(eff_thermal_forcing**2 + dbmlt_float(i,j)/coeff)
                else  ! add enough warming to change TF from negative to positive
                   dTocn(i,j) = sqrt(dbmlt_float(i,j)/coeff) - eff_thermal_forcing
                endif
             endif

             if (verbose_bmlt_float .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'In ismip6_set_deltaT_ocn_basin, r, i, j =', rtest, itest, jtest
                write(iulog,*) 'dthck_dt_target (m/yr)=', dthck_dt_target(i,j)*scyr
                write(iulog,*) 'thermal_forcing_lsrf =', thermal_forcing_lsrf(i,j)
                write(iulog,*) 'deltaT_ocn_init =', deltaT_ocn_init(i,j)
                write(iulog,*) 'dTocn adjustment =', dTocn(i,j)
                write(iulog,*) 'deltaT_ocn_new =', deltaT_ocn_init(i,j) + dTocn(i,j)
             endif

          enddo
       enddo

    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL) then

       ! nonlocal parameterization
       ! melt rate is a quadratic function of local thermal forcing and basin-average thermal forcing
       ! m = C * Fb * F0,
       !    where C = coeff * gamma0, Fb = basin-average thermal forcing, and F0 = local thermal forcing
       ! If Fb > 0 and F0 > 0, the increase in melt rate due to ocean warming dT is given by
       !    dm = C * Fb * dT, implying dT = C * Fb / dm
       ! If F0 < 0, then we have
       !    dm = C * Fb * (F0 + dT), implying dT = dm/(C*Fb) - F0
       ! Since Fb is a function of F0 throughout the basin, the subroutine should be called iteratively.

       dTocn(i,j) = 0.0d0

       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing_basin = max(0.0d0, thermal_forcing_basin(nb))
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn_init(i,j)
                if (eff_thermal_forcing > 0.0d0) then  ! increment positive TF
                   bmlt_float_init(i,j) = coeff * eff_thermal_forcing_basin * eff_thermal_forcing
                   dTocn(i,j) = dbmlt_float(i,j) / (coeff * eff_thermal_forcing_basin)
                else  ! add enough warming to change TF from negative to positive
                   dTocn(i,j) = dbmlt_float(i,j) / (coeff * eff_thermal_forcing_basin) - eff_thermal_forcing
                endif
             endif

             if (verbose_bmlt_float .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'In ismip6_set_deltaT_ocn_for_dthck_dt, r, i, j, nb =', rtest, itest, jtest, nb
                write(iulog,*) 'thermal_forcing_lsrf =', thermal_forcing_lsrf(i,j)
                write(iulog,*) 'thermal_forcing_basin =', thermal_forcing_basin(nb)
                write(iulog,*) 'dthck_dt_target (m/yr)=', dthck_dt_target(i,j)*scyr
                write(iulog,*) 'deltaT_ocn_init =', deltaT_ocn_init(i,j)
                write(iulog,*) 'dTocn adjustment =', dTocn(i,j)
                write(iulog,*) 'deltaT_ocn_new =', deltaT_ocn_init(i,j) + dTocn(i,j)
             endif

          enddo
       enddo

    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then

       ! same as nonlocal, but with larger gamma0, and multiplied by sin(theta_slope)

       dTocn(i,j) = 0.0d0

       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing_basin = max(0.0d0, thermal_forcing_basin(nb))
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn_init(i,j)
                if (eff_thermal_forcing > 0.0d0) then  ! increment positive TF
                   bmlt_float_init(i,j) = &
                        coeff * sin(theta_slope(i,j)) * eff_thermal_forcing_basin * eff_thermal_forcing
                   dTocn(i,j) = dbmlt_float(i,j) / &
                        (coeff * sin(theta_slope(i,j)) * eff_thermal_forcing_basin)
                else  ! add enough warming to change TF from negative to positive
                   dTocn(i,j) = dbmlt_float(i,j) / &
                        (coeff * sin(theta_slope(i,j)) * eff_thermal_forcing_basin) - eff_thermal_forcing
                endif
             endif
          enddo
       enddo

    endif  ! bmlt_float_thermal_forcing_param

    ! Adjust deltaT_ocn
    deltaT_ocn_new = deltaT_ocn_init + dTocn

    ! For diagnostics, compute the new melt rate

    if (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL) then
       do j = 1, ny
          do i = 1, nx
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn_new(i,j)
                if (eff_thermal_forcing > 0.0d0) then  ! increment positive TF
                   bmlt_float_new(i,j) = coeff * eff_thermal_forcing**2
                endif
             endif
          enddo
       enddo
    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL) then
       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing_basin = max(0.0d0, thermal_forcing_basin(nb))
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn_new(i,j)
                if (eff_thermal_forcing > 0.0d0) then  ! increment positive TF
                   bmlt_float_new(i,j) = coeff * eff_thermal_forcing_basin * eff_thermal_forcing
                endif
             endif
          enddo
       enddo
    elseif (bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then
       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (thermal_forcing_mask(i,j) == 1) then
                eff_thermal_forcing_basin = max(0.0d0, thermal_forcing_basin(nb))
                eff_thermal_forcing = thermal_forcing_lsrf(i,j) + deltaT_ocn_new(i,j)
                if (eff_thermal_forcing > 0.0d0) then  ! increment positive TF
                   bmlt_float_new(i,j) = &
                        coeff * sin(theta_slope(i,j)) * eff_thermal_forcing_basin * eff_thermal_forcing
                endif
             endif
          enddo
       enddo
    endif

    if (verbose_bmlt_float) then
       call point_diag(thermal_forcing_lsrf, 'thermal_forcing_lsrf (degC)', itest, jtest, rtest, 7, 7)
       call point_diag(deltaT_ocn_init, 'Initial deltaT_ocn', itest, jtest, rtest, 7, 7)
       call point_diag(thermal_forcing_lsrf + deltaT_ocn_init, 'Initial effective TF', itest, jtest, rtest, 7, 7)
       call point_diag(bmlt_float_init*scyr, 'Initial melt rate (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(dTocn, 'deltaT_ocn adjustment', itest, jtest, rtest, 7, 7)
       call point_diag(deltaT_ocn_new, 'New deltaT_ocn', itest, jtest, rtest, 7, 7)
       call point_diag(thermal_forcing_lsrf + deltaT_ocn_new, 'New effective TF', itest, jtest, rtest, 7, 7)
       call point_diag(bmlt_float_new*scyr, 'New melt rate (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag((bmlt_float_new - bmlt_float_init)*scyr, 'Melt difference (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(dthck_dt_target*scyr, 'dthck_dt_target (m/yr)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine set_deltaT_ocn_for_dthck_dt

!****************************************************

  !TODO - Deprecate this subroutine?
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
         bmlt_float               !> basal melt rate (m/s) at lower ice surface

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

    ! Compute basal melt rates (m/s)
    ! Note: The applied thermal forcing must be non-negative.  TF < 0 is inappropriate for a quadratic function.
    do j = 1, ny
       do i = 1, nx
          if (thermal_forcing_mask(i,j) == 1) then
             thermal_forcing = max(thermal_forcing_lsrf(i,j), 0.0d0)
             bmlt_float(i,j) = rhoo * cpo * gamma_t * Fm * thermal_forcing**2 / (lhci * rhoi)
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
    !TODO - Replace local_count with calls to parallel_global_sum?
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

    logical, parameter :: verbose_basin_number = .false.

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
       write(iulog,*) 'Extrapolating basin numbers to cells with invalid values'
       write(iulog,*) 'Initial count of valid values:', global_count_save
       write(iulog,*) 'max_iter =', max_iter
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
!!          write(iulog,*) iter, 'Basin number count =', global_count
       endif

       if (global_count == global_count_save) then
          if (verbose_basin_number .and. main_task) then
             write(iulog,*) 'Exiting basin_number_extrapolate, iter =', iter
          endif
          exit
       else
          global_count_save = global_count
       endif

       if (iter == max_iter) then
          if (main_task) write(iulog,*) 'Error: Exiting basin_number_extrapolate, max_iter =', max_iter
          call write_log('Error: Exiting basin_number_extrapolate without converging', GM_FATAL)
       endif

    enddo   ! iter

  end subroutine basin_number_extrapolate

!****************************************************

end module glissade_bmlt_float

!****************************************************

