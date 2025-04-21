!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_inversion.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2017
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

module glissade_inversion

  use glimmer_physcon, only: scyr, grav
!!  use glimmer_paramets, only: eps08, tim0, len0, vel0, thk0
  use glimmer_paramets, only: eps08, vel0, thk0
  use glimmer_log
  use glide_types
  use glide_thck, only: glide_calclsrf
  use glide_diagnostics, only: point_diag
  use cism_parallel, only: this_rank, main_task, nhalo, &
       parallel_type, parallel_halo, staggered_parallel_halo, &
       parallel_reduce_min, parallel_reduce_max

  implicit none

  private
  public :: verbose_inversion, glissade_inversion_init, glissade_inversion_solve

  !-----------------------------------------------------------------------------
  ! Subroutines to invert for basal fields (including basal friction beneath
  ! grounded ice and basal melting beneath floating ice) by relaxing toward
  ! a target ice thickness field.
  !-----------------------------------------------------------------------------

!    logical, parameter :: verbose_inversion = .false.
    logical, parameter :: verbose_inversion = .true.

!***********************************************************************

contains

!***********************************************************************

  subroutine glissade_inversion_init(model)

    ! Initialize inversion for fields of basal friction and basal melting
    ! Should be called after usrf and thck have been input and (possibly) modified by initial calving

    use glissade_masks, only: glissade_get_masks
    use glissade_grid_operators, only: glissade_stagger
    use glissade_basal_traction, only: set_coulomb_c_elevation
    use glissade_utils, only: glissade_usrf_to_thck, glissade_thck_to_usrf

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    integer :: i, j
    integer :: nb                   ! basin number
    integer :: itest, jtest, rtest  ! local diagnostic point

    real(dp) :: var_maxval          ! max value of a given real variable; = 0.0 if not yet read in
    real(dp) :: var1_maxval, var2_maxval  ! like var_maxval
    integer :: var_maxval_int       ! max value of a given integer variable; = 0 if not yet read in

    character(len=100) :: message

    integer, dimension(model%general%ewn, model%general%nsn) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         ocean_mask,           & ! = 1 where topg is below sea level and ice is absent
         land_mask               ! = 1 where topg is at or above sea level

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         f_flotation,          & ! flotation function (m)
         thck_obs                ! observed ice thickness, derived from usrf_obs and topg

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         coulomb_c_icegrid      ! initial coulomb_c at cell centers based on masks

    real(dp) :: h_obs, h_flotation, h_buff   ! thck_obs, flotation thickness, and thck_flotation_buffer scaled to m
    real(dp) :: dh                           ! h_obs - h_flotation
    real(dp) :: dh_decimal                   ! decimal part remaining after subtracting the truncation of dh

    integer :: ewn, nsn

    type(parallel_type) :: parallel    ! info for global communication

    parallel = model%parallel

    ewn = model%general%ewn
    nsn = model%general%nsn

    ! Set local diagnostic point
    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    !----------------------------------------------------------------------
    ! If inverting for Cp or Cc, then set the target elevation, usrf_obs.
    !----------------------------------------------------------------------

    if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or.  &
        model%options%which_ho_coulomb_c  == HO_COULOMB_C_INVERSION  .or.  &
        model%options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION) then

       ! We are likely trying to match usrf_obs, so check whether it has been read in already.
       ! If not, set it to the initial usrf field.

       var_maxval = maxval(model%geometry%usrf_obs)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing
       else
          ! initialize to the initial upper surface elevation
          model%geometry%usrf_obs(:,:) = model%geometry%usrf(:,:)
       endif

       ! Given usrf_obs and topg, compute thck_obs.

       call glissade_usrf_to_thck(&
            model%geometry%usrf_obs,  &
            model%geometry%topg,      &
            model%climate%eus,        &
            thck_obs)

       if (model%options%is_restart == NO_RESTART) then

          ! At the start of the run, adjust thck_obs so that the observational target is not too close to thck_flotation.
          ! The reason for this is that if we restore H to values very close to thck_flotation,
          !  it can be easier for cells to flip between grounded and floating in the forward run.
          ! Note: The f_ground algorithm can have problems if we have 4 adjacent cells with a checkerboard pattern
          !        of equal and opposite values of dh = thck_obs - thck_flotation.
          !       To prevent this, do not make the new dh exactly equal to thck_flotation_buffer,
          !        but rather increment it up or down in integer steps until we exceed the buffer.
          !       E.g., if the buffer is 1 m, then a cavity of 0.357 m is increased to 1.357 m.

          do j = 1, nsn
             do i = 1, ewn
                if (thck_obs(i,j) > 0.0d0 .and.  &
                     model%geometry%topg(i,j) - model%climate%eus < 0.0d0) then   ! marine-based ice
                   ! convert to meters (can skip the conversion when code scaling is removed)
                   h_obs = thck_obs(i,j) * thk0
                   h_flotation = -(rhoo/rhoi) * (model%geometry%topg(i,j) - model%climate%eus) * thk0
                   h_buff = model%inversion%thck_flotation_buffer * thk0
                   dh = h_obs - h_flotation
                   if (abs(dh) < h_buff) then
                      if (dh > 0.0d0) then
                         dh_decimal = dh - floor(dh)
                         h_obs = h_flotation + h_buff + dh_decimal
                      else   ! dh <= 0.0
                         dh_decimal = ceiling(dh) - dh
                         h_obs = h_flotation - h_buff - dh_decimal
                      endif
                      thck_obs(i,j) = h_obs / thk0
                   endif
                endif
                thck_obs(i,j) = max(thck_obs(i,j), 0.0d0)
             enddo
          enddo

          ! Where thck_obs < inversion_thck_threshold, set it to zero.
          ! One reason to do this is to avoid restoring ice to small values at the calving front.
          ! Probably not necessary if doing basin-scale inversion for floating ice instead of inversion in each grid cell.
          ! The factor of eps08 is included as a buffer against rounding errors.
          ! For instance, if the threshold is 1.0 m, we don't want to remove ice that has H = 1.0 - 1.e-13 due to rounding.
          do j = nhalo+1, nsn-nhalo
             do i = nhalo+1, ewn-nhalo
                if (verbose_inversion .and. thck_obs(i,j) > 0.0d0 .and. &
                     thck_obs(i,j) + eps08 < model%inversion%thck_threshold) then
                   !WHL - debug
!!                   write(6,*) 'thck_obs < threshold, rank, i, j, thck:', this_rank, i, j, thck_obs(i,j)*thk0
                endif
             enddo
          enddo

          model%inversion%thck_threshold = max(model%inversion%thck_threshold, model%numerics%thklim)
          where (thck_obs + eps08 <= model%inversion%thck_threshold)
             thck_obs = 0.0d0
          endwhere

          ! Set thck to be consistent with thck_obs
          model%geometry%thck = thck_obs

          ! Reset usrf_obs to be consistent with thck_obs.
          ! (usrf itself will be recomputed later in glissade_initialise)
          call glissade_thck_to_usrf(&
               thck_obs,  &
               model%geometry%topg,      &
               model%climate%eus,        &
               model%geometry%usrf_obs)

       endif   ! not a restart

       call parallel_halo(model%geometry%usrf_obs, parallel)
       call parallel_halo(thck_obs, parallel)

    endif  ! inversion for Cp, Cc or deltaT_ocn

    !----------------------------------------------------------------------
    ! If inverting for the flow enhancement factor E, then set the target
    ! surface ice speed, velo_sfc_obs.
    !----------------------------------------------------------------------

    if (model%options%which_ho_flow_enhancement_factor == HO_FLOW_ENHANCEMENT_FACTOR_INVERSION) then

       ! Make sure that either (1) usfc_obs and vsfc_obs were read in
       ! or (2) velo_sfc_obs was read in
       var_maxval = maxval(model%velocity%velo_sfc_obs)
       var_maxval = parallel_reduce_max(var_maxval)

       var1_maxval = maxval(model%velocity%usfc_obs)
       var1_maxval = parallel_reduce_max(var1_maxval)

       var2_maxval = maxval(model%velocity%vsfc_obs)
       var2_maxval = parallel_reduce_max(var2_maxval)

       if (var_maxval > 0.0d0) then
          ! velo_sfc_obs was read in; do nothing

       elseif (var1_maxval > 0.0d0 .and. var2_maxval > 0.0d0) then
          ! usfc_obs and vsfc_obs were read in; compute velo_sfc_obs
          model%velocity%velo_sfc_obs(:,:) = &
               sqrt(model%velocity%usfc_obs(:,:)**2 + model%velocity%vsfc_obs(:,:)**2)

       else ! abort
          call write_log('Error, must read in velo_sfc_obs, usfc_obs and/or vsfc_obs to invert for E', GM_FATAL)
       endif

    endif   ! flow enhancement factor inversion

    if (verbose_inversion) then
       call point_diag(model%temper%flow_enhancement_factor, 'init_inversion: flow_enhancement E =', itest, jtest, rtest, 7, 7)
    endif

    ! Set masks that are used below
    ! Modify glissade_get_masks so that 'parallel' is not needed
    call glissade_get_masks(ewn,                 nsn,                   &
                            parallel,                                   &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   model%numerics%thklim, &
                            ice_mask,                                   &
                            floating_mask = floating_mask,              &
                            ocean_mask = ocean_mask,                    &
                            land_mask = land_mask)

    !----------------------------------------------------------------------
    ! computations specific to powerlaw_c (Cp) and coulomb_c (Cc) inversion
    !----------------------------------------------------------------------

    if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION) then

       ! initialize powerlaw_c, if not already read in
       var_maxval = maxval(model%basal_physics%powerlaw_c)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; powerlaw_c has been read in already (e.g., when restarting)
       else
          ! initialize to a uniform value (which can be set in the config file)
          model%basal_physics%powerlaw_c(:,:) = model%basal_physics%powerlaw_c_const
       endif  ! var_maxval > 0

       if (verbose_inversion) then
          call point_diag(model%basal_physics%powerlaw_c, 'init_inversion for powerlaw_c', itest, jtest, rtest, 7, 7)
       endif

    elseif (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION) then

       ! initialize coulomb_c, if not already read in
       var_maxval = maxval(model%basal_physics%coulomb_c)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; coulomb_c has been read in already (e.g., when restarting)
       else

          if (model%options%which_ho_coulomb_c_relax == HO_COULOMB_C_RELAX_NONE .or.  &
              model%options%which_ho_coulomb_c_relax == HO_COULOMB_C_RELAX_CONSTANT) then

             !Note: The initial value used to depend on whether cells were floating or grounded.
             !      Now, we initialize coulomb_c to the same value everywhere.

             model%basal_physics%coulomb_c = model%basal_physics%coulomb_c_const
             model%basal_physics%coulomb_c_relax = model%basal_physics%coulomb_c_const

          elseif (model%options%which_ho_coulomb_c_relax == HO_COULOMB_C_RELAX_ELEVATION) then

             ! Set coulomb_c_relax based on bed elevation, and set coulomb_c = coulomb_c_relax.
             ! Note: If the bed topography is fixed, coulomb_c_relax could be set once and for all.
             !       If isostasy is on, coulomb_c_relax needs to be reset as the bed evolves.

             call set_coulomb_c_elevation(&
                  ewn,                nsn,                   &
                  model%geometry%topg*thk0,                  &  ! m
                  model%climate%eus*thk0,                    &  ! m
                  model%basal_physics%coulomb_c_relax_min,   &
                  model%basal_physics%coulomb_c_relax_max,   &
                  model%basal_physics%coulomb_c_bedmin,      &  ! m
                  model%basal_physics%coulomb_c_bedmax,      &  ! m
                  model%basal_physics%coulomb_c)

             model%basal_physics%coulomb_c_relax = model%basal_physics%coulomb_c

          endif

       endif  ! var_maxval > 0

       if (verbose_inversion) then
          call point_diag(model%basal_physics%coulomb_c, 'init_inversion for coulomb_c', itest, jtest, rtest, 7, 7)
       endif

    endif   ! Cp or Cc inversion

    !----------------------------------------------------------------------
    ! computations specific to flow_enhancement_factor inversion
    !----------------------------------------------------------------------

    if (model%options%which_ho_flow_enhancement_factor == HO_FLOW_ENHANCEMENT_FACTOR_INVERSION) then

       ! initialize flow_enhancement_factor, if not already read in
       var_maxval = maxval(model%temper%flow_enhancement_factor)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; flow_enhancement_factor has been read in already (e.g., when restarting)
       else

          ! initialize to the default values for grounded and floating ice
          ! For ice-free ocean, flow_enhancement_factor = 0 for now, but will change if ice-covered

          where (floating_mask == 1)
             model%temper%flow_enhancement_factor = model%paramets%flow_enhancement_factor_float
          elsewhere (ice_mask == 1 .or. land_mask == 1)  ! grounded ice or land
             model%temper%flow_enhancement_factor = model%paramets%flow_enhancement_factor_ground
          endwhere

       endif   ! var_maxval > 0

    endif   ! flow_enhancement_factor inversion

    !----------------------------------------------------------------------
    ! computations specific to basin-scale ocean temperature inversion
    !----------------------------------------------------------------------

    if (model%options%which_ho_deltaT_basin == HO_DELTAT_BASIN_INVERSION) then

       if (model%options%is_restart == NO_RESTART) then

          ! Set floating_thck_target for floating ice and lightly grounded ice.
          ! Here, "lightly grounded" means that the magnitude of f_flotation = (-topg - eus) - (rhoi/rhoo)*thck
          !  is less than a prescribed threshold.  (Recall f_flotation < 0 for grounded ice.)
          ! The inversion will nudge the ice thickness toward this target in a basin-average sense.
          ! Positive volume biases will be corrected with ocean warming or ice softening,
          !  and negative biases with ocean cooling or ice stiffening.

          do j = 1, nsn
             do i = 1, ewn
                f_flotation(i,j) = (-(model%geometry%topg(i,j) - model%climate%eus)  &
                              - (rhoi/rhoo)*model%geometry%thck(i,j)) * thk0    ! f_flotation < 0 for grounded ice
                if (model%geometry%thck(i,j) > 0.0d0 .and. &
                    model%geometry%marine_connection_mask(i,j) == 1 .and. &
                    f_flotation(i,j) > -model%inversion%basin_flotation_threshold) then
                   model%inversion%floating_thck_target(i,j) = model%geometry%thck(i,j)
                else
                   model%inversion%floating_thck_target(i,j) = 0.0d0
                endif
             enddo
          enddo

          if (verbose_inversion) then
             call point_diag(model%inversion%floating_thck_target*thk0, &
                  'After init_inversion, floating_thck_target', itest, jtest, rtest, 7, 7)
             call point_diag(model%geometry%thck*thk0, 'thck', itest, jtest, rtest, 7, 7)
             call point_diag(f_flotation, 'f_flotation (m)', itest, jtest, rtest, 7, 7)
          endif   ! verbose

       endif   ! not a restart

       call parallel_halo(model%inversion%floating_thck_target, parallel)

    endif  ! which_ho_deltaT_basin_inversion

    if (verbose_inversion) then
       call point_diag(model%geometry%usrf_obs*thk0, &
            'After init_inversion, usrf_obs  (m)', itest, jtest, rtest, 7, 7)
       call point_diag(thck_obs*thk0, 'thck_obs  (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_inversion_init

!***********************************************************************

  subroutine glissade_inversion_solve(model)

    ! Invert for one or more of several quantities, including powerlaw_c and coulomb_c
    ! (beneath grounded ice), deltaT_ocn (beneath floating ice), and flow_enhancement_factor.

    use glissade_masks, only: glissade_get_masks
    use glissade_bmlt_float, only: glissade_bmlt_float_thermal_forcing
    use glissade_grounding_line, only: glissade_grounded_fraction
    use glissade_utils, only: glissade_usrf_to_thck, glissade_basin_average
    use glissade_grid_operators, only: glissade_stagger
    use glissade_basal_traction, only: set_coulomb_c_elevation

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ice_mask,           & ! = 1 where thck > thklim, else = 0
         floating_mask,      & ! = 1 where ice is present and floating, else = 0
         ocean_mask,         & ! = 1 where topg is below sea level and ice is absent
         land_mask             ! = 1 where topg is at or above sea level

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck_obs,           & ! observed thickness target (m)
         f_ground_cell_obs,  & ! f_ground_cell as a function of thck_obs (instead of current thck)
         f_ground_obs,       & ! f_ground as a function of thck_obs (instead of current thck)
         f_flotation_obs       ! f_flotation_obs as a function of thck_obs (instead of current thck)

    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) ::   &
         stag_thck,          & ! ice thickness on staggered grid
         stag_dthck_dt,      & ! dthck_dt on staggered grid
         stag_thck_obs         ! thck_obs on staggered grid

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         deltaT_ocn_relax                      ! relax deltaT_ocn toward this value

    real(dp), dimension(model%ocean_data%nbasin) :: &
         deltaT_ocn_basin_avg                  ! basin average of deltaT_ocn (degC)

    logical :: invert_coulomb_c, invert_powerlaw_c

    type(parallel_type) :: parallel  ! info for parallel communication
    integer :: ewn, nsn
    integer :: itest, jtest, rtest   ! local diagnostic point
    integer :: i, j, nb

    real(dp) :: phaseout_factor      ! multiplying factor for inversion timescales
    real(dp) :: babc_timescale, deltaT_ocn_timescale, flow_enhancement_timescale  ! copies of model parameters

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

    call glissade_get_masks(ewn,                 nsn,                   &
                            parallel,                                   &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   model%numerics%thklim, &
                            ice_mask,                                   &
                            floating_mask = floating_mask,              &
                            ocean_mask = ocean_mask,                    &
                            land_mask = land_mask)

    ! Make copies of the inversion timescales
    babc_timescale = model%inversion%babc_timescale
    deltaT_ocn_timescale = model%inversion%deltaT_ocn_timescale
    flow_enhancement_timescale = model%inversion%flow_enhancement_timescale

    ! Optionally, increase the timescales so that the inversion is phased out over time.
    ! Note: The terms proportional to dH/dt do not depend on this timescale, so the
    !       inversion will keep acting to prevent the thickness from changing further.

    if (model%inversion%phaseout_timescale > 0.0d0) then

       ! Increase the timescales
       phaseout_factor = exp(model%numerics%time/model%inversion%phaseout_timescale)
       babc_timescale = babc_timescale * phaseout_factor
       deltaT_ocn_timescale = deltaT_ocn_timescale * phaseout_factor
       flow_enhancement_timescale = flow_enhancement_timescale * phaseout_factor
       if (verbose_inversion .and. this_rank == rtest) then
          print*, 'Inversion phaseout factor, babc_timescale (yr)=', phaseout_factor, babc_timescale/scyr
       endif
    endif

    ! If inverting for Cp = powerlaw_c or Cc = coulomb_c, then update it here.
    ! If running with glaciers, inversion for powerlaw_c is done elsewhere,
    !  in subroutine glissade_glacier_update.

    if ( model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or. &
         model%options%which_ho_coulomb_c  == HO_COULOMB_C_INVERSION) then

       ! Given the surface elevation target, compute the thickness target.
       ! (This can change in time if the bed topography is dynamic.)

       call glissade_usrf_to_thck(&
            model%geometry%usrf_obs,  &
            model%geometry%topg,      &
            model%climate%eus,        &
            thck_obs)

       ! Interpolate the thickness fields to the staggered grid

       call glissade_stagger(ewn,         nsn,              &
                             thck_obs,    stag_thck_obs)

       call glissade_stagger(ewn,                  nsn,             &
                             model%geometry%thck,  stag_thck)

       call glissade_stagger(ewn,                      nsn,             &
                             model%geometry%dthck_dt,  stag_dthck_dt)

       call staggered_parallel_halo(stag_thck_obs, parallel)
       call staggered_parallel_halo(stag_thck, parallel)
       call staggered_parallel_halo(stag_dthck_dt, parallel)

       ! Invert for powerlaw_c or coulomb_c
       ! The logic is the same for each; only the max and min values and the in/out field are different.

       if ( model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION) then

          if (verbose_inversion) then
             call point_diag(model%basal_physics%powerlaw_c, 'powerlaw_c', itest, jtest, rtest, 7, 7)
          endif

          !TODO - Add an option to set based on elevation, like coulomb_c?
          model%basal_physics%powerlaw_c_relax = model%basal_physics%powerlaw_c_const

          call invert_basal_friction(&
!!               model%numerics%dt*tim0,                   &  ! s
               model%numerics%dt,                        &  ! s
               ewn,               nsn,                   &
!!               model%numerics%dew*len0,                  &  ! m
!!               model%numerics%dns*len0,                  &  ! m
               model%numerics%dew,                       &  ! m
               model%numerics%dns,                       &  ! m
               itest,    jtest,   rtest,                 &
               model%inversion%babc_thck_scale,          &  ! m
               babc_timescale,                           &  ! s
               model%inversion%babc_length_scale,        &  ! m
               model%inversion%babc_relax_factor,        &
               model%basal_physics%powerlaw_c_max,       &
               model%basal_physics%powerlaw_c_min,       &
               model%geometry%f_ground,                  &
               stag_thck*thk0,                           &  ! m
               stag_thck_obs*thk0,                       &  ! m
               stag_dthck_dt,                            &  ! m/s
               model%basal_physics%powerlaw_c_relax,     &
               model%basal_physics%powerlaw_c)

          ! halo update
          call staggered_parallel_halo(model%basal_physics%powerlaw_c, parallel)

          if (verbose_inversion) then
             call point_diag(model%basal_physics%powerlaw_c, 'powerlaw_c', itest, jtest, rtest, 7, 7)
          endif

       elseif ( model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION) then

          ! Set the relaxation target, coulomb_c_relax

          if (model%options%which_ho_coulomb_c_relax == HO_COULOMB_C_RELAX_CONSTANT) then

             model%basal_physics%coulomb_c_relax = model%basal_physics%coulomb_c_const

          elseif (model%options%which_ho_coulomb_c_relax == HO_COULOMB_C_RELAX_ELEVATION) then
             ! set coulomb_c_relax based on bed elevation
             ! Note: Could be called once at initialization, if the bed topography is fixed

             call set_coulomb_c_elevation(&
                  ewn,                nsn,                   &
                  model%geometry%topg*thk0,                  &  ! m
                  model%climate%eus*thk0,                    &  ! m
                  model%basal_physics%coulomb_c_relax_min,   &
                  model%basal_physics%coulomb_c_relax_max,   &
                  model%basal_physics%coulomb_c_bedmin,      &  ! m
                  model%basal_physics%coulomb_c_bedmax,      &  ! m
                  model%basal_physics%coulomb_c_relax)

          else

             model%basal_physics%coulomb_c_relax = 0.0d0  ! no relaxation

          endif

          if (verbose_inversion) then
             call point_diag(model%basal_physics%coulomb_c, 'Old coulomb_c', itest, jtest, rtest, 7, 7, '(f10.5)')
             call point_diag(model%basal_physics%effecpress_stag, 'effecpress_stag', &
                  itest, jtest, rtest, 7, 7, '(f10.1)')
             call point_diag(rhoi*grav*stag_thck*thk0, 'overburden', itest, jtest, rtest, 7, 7, '(f10.1)')
             call point_diag((model%geometry%thck - thck_obs)*thk0, 'thck - thck_obs (m)', &
                  itest, jtest, rtest, 7, 7)
             call point_diag(model%basal_physics%coulomb_c_relax, 'coulomb_c_relax', itest, jtest, rtest, 7, 7, '(f10.5)')
          endif

          call invert_basal_friction(&
!!               model%numerics%dt*tim0,                   &  ! s
               model%numerics%dt,                        &  ! s
               ewn,               nsn,                   &
!!               model%numerics%dew*len0,                  &  ! m
!!               model%numerics%dns*len0,                  &  ! m
               model%numerics%dew,                       &  ! m
               model%numerics%dns,                       &  ! m
               itest,    jtest,   rtest,                 &
               model%inversion%babc_thck_scale,          &  ! m
               babc_timescale,                           &  ! s
               model%inversion%babc_length_scale,        &  ! m
               model%inversion%babc_relax_factor,        &
               model%basal_physics%coulomb_c_max,        &
               model%basal_physics%coulomb_c_min,        &
               model%geometry%f_ground,                  &
               stag_thck*thk0,                           &  ! m
               stag_thck_obs*thk0,                       &  ! m
               stag_dthck_dt,                            &  ! m/s
               model%basal_physics%coulomb_c_relax,      &
               model%basal_physics%coulomb_c)

          ! halo update
          call staggered_parallel_halo(model%basal_physics%coulomb_c, parallel)

          if (verbose_inversion .and. this_rank == rtest) then
             call point_diag(model%basal_physics%coulomb_c, 'New coulomb_c', itest, jtest, rtest, 7, 7, '(f10.5)')
          endif   ! verbose_inversion

       endif  ! invert for powerlaw_c or coulomb_c

    else   ! do not invert for powerlaw_c or coulomb_c; just print optional diagnostics

       if (verbose_inversion) then
          call point_diag(model%geometry%f_ground, 'f_ground at vertices', itest, jtest, rtest, 7, 7, '(f10.4)')
          call point_diag(model%basal_physics%powerlaw_c, 'powerlaw_c', itest, jtest, rtest, 7, 7, '(f10.2)')
          call point_diag(model%basal_physics%coulomb_c, 'coulomb_c', itest, jtest, rtest, 7, 7, '(f10.4)')
       endif

    endif   ! invert for powerlaw_c or coulomb_c

    ! Replace zeroes (if any) with small nonzero values to avoid divzeroes.
    ! Note: The current algorithm initializes Cc to a nonzero value everywhere and never sets Cp = 0;
    !       this code is just to be on the safe side.

    if (model%options%which_ho_powerlaw_c /= HO_POWERLAW_C_CONSTANT) then
       where (model%basal_physics%powerlaw_c == 0.0d0)
          model%basal_physics%powerlaw_c = model%basal_physics%powerlaw_c_min
       endwhere
    endif

    if (model%options%which_ho_coulomb_c /= HO_COULOMB_C_CONSTANT) then
       where (model%basal_physics%coulomb_c == 0.0d0)
          model%basal_physics%coulomb_c = model%basal_physics%coulomb_c_min
       endwhere
    endif


    ! If inverting for deltaT_ocn at the basin level, then update it here

    if ( model%options%which_ho_deltaT_basin == HO_DELTAT_BASIN_INVERSION) then

       call glissade_inversion_deltaT_basin(&
!!            model%numerics%dt * tim0,                  &  ! s
            model%numerics%dt,                         &  ! s
            ewn, nsn,                                  &
!!            model%numerics%dew * len0,                 &  ! m
!!            model%numerics%dns * len0,                 &  ! m
            model%numerics%dew,                        &  ! m
            model%numerics%dns,                        &  ! m
            itest, jtest, rtest,                       &
            model%ocean_data%nbasin,                   &
            model%ocean_data%basin_number,             &
            model%geometry%thck*thk0,                  &  ! m
            model%geometry%dthck_dt,                   &  ! m/s
            model%inversion%floating_thck_target*thk0, &  ! m
            model%inversion%deltaT_ocn_thck_scale,     &  ! m
            deltaT_ocn_timescale,                      &  ! s
            model%inversion%deltaT_ocn_temp_scale,     &  ! degC
            model%inversion%deltaT_basin_relax,        &  ! degC
            model%inversion%basin_mass_correction,     &
            model%inversion%basin_number_mass_correction, &
            model%ocean_data%deltaT_ocn)

    endif   ! which_ho_deltaT_basin

    ! If inverting for deltaT_ocn based on observed ice thickness, then update it here.

    if ( model%options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION) then

       ! Given the surface elevation target, compute the thickness target.
       ! This can change in time if the bed topography is dynamic.

       call glissade_usrf_to_thck(&
            model%geometry%usrf_obs,  &
            model%geometry%topg,      &
            model%climate%eus,        &
            thck_obs)

       ! Set the value toward which we relax deltaT_ocn during inversion.
       ! If running with ocean basins, we relax toward the basin-average value,
       !  computed over all floating cells in the basin.

       deltaT_ocn_relax(:,:) = 0.0d0

       if (model%ocean_data%nbasin > 1) then
          call glissade_basin_average(&
               model%general%ewn, model%general%nsn,  &
               model%ocean_data%nbasin,               &
               model%ocean_data%basin_number,         &
               floating_mask * 1.0d0,                 &   ! real mask
               model%ocean_data%deltaT_ocn,           &
               deltaT_ocn_basin_avg)
       endif

       do j = nhalo+1, nsn-nhalo
          do i = nhalo+1, ewn-nhalo
             nb = model%ocean_data%basin_number(i,j)
             deltaT_ocn_relax(i,j) = deltaT_ocn_basin_avg(nb)
          enddo
       enddo

       ! Given the thickness target, invert for deltaT_ocn

       call glissade_inversion_deltaT_ocn(&
!!            model%numerics%dt * tim0,              &  ! s
            model%numerics%dt,                     &  ! s
            ewn,           nsn,                    &
!!            model%numerics%dew*len0,               &   ! m
!!            model%numerics%dns*len0,               &   ! m
            model%numerics%dew,                    &   ! m
            model%numerics%dns,                    &   ! m
            itest, jtest,  rtest,                  &
            model%inversion%deltaT_ocn_thck_scale, &  ! m
            deltaT_ocn_timescale,                  &  ! s
            model%inversion%deltaT_ocn_temp_scale, &  ! degC
            model%inversion%deltaT_ocn_length_scale,& ! m
            deltaT_ocn_relax,                      &  ! degC
            model%geometry%f_ground_cell,          &
            model%geometry%thck * thk0,            &  ! m
            thck_obs * thk0,                       &  ! m
            model%geometry%dthck_dt,               &  ! m/s
            model%ocean_data%deltaT_ocn)              ! degC

       call parallel_halo(model%ocean_data%deltaT_ocn, parallel)

    endif   ! which_ho_deltaT_ocn

    ! If setting deltaT_ocn based on observed dthck_dt, then do so here.
    ! TODO - Deprecate this option?
    if (model%options%which_ho_deltat_ocn == HO_DELTAT_OCN_DTHCK_DT) then

       ! Set deltaT_ocn based on dthck_dt_obs.
       ! This is done within the subroutine used to compute bmlt_float from thermal forcing.
       ! But instead of computing bmlt_float from TF, we find the value of deltaT_ocn
       !  that will increase TF as needed to match negative values of dthck_dt_obs.
       ! Note: This subroutine would usually be called during the initial diagnostic solve
       !       of the restart following a spin-up, without taking any prognostic timesteps.

       call glissade_bmlt_float_thermal_forcing(&
            model%options%bmlt_float_thermal_forcing_param, &
            model%options%ocean_data_extrapolate,     &
            parallel,                                 &
            ewn,       nsn,                           &
!!            model%numerics%dew*len0,                  &   ! m
!!            model%numerics%dns*len0,                  &   ! m
            model%numerics%dew,                       &   ! m
            model%numerics%dns,                       &   ! m
            itest,     jtest,   rtest,                &
            ice_mask,                                 &
            ocean_mask,                               &
            model%geometry%marine_connection_mask,    &
            model%geometry%f_ground_cell,             &
            model%geometry%thck*thk0,                 &   ! m
            model%geometry%lsrf*thk0,                 &   ! m
            model%geometry%topg*thk0,                 &   ! m
            model%ocean_data,                         &
            model%basal_melt%bmlt_float,              &
            which_ho_deltaT_ocn = model%options%which_ho_deltaT_ocn,  &
            dthck_dt_obs = model%geometry%dthck_dt_obs)   ! m/yr

    endif   ! which_ho_deltaT_ocn

    !WHL - debug
    ! For testing subgrid CF schemes: Do not invert for deltaT_ocn where calving_mask = 1,
    ! because this will prevent CF advance. In these cells, set deltaT_ocn = 0.
    if (model%options%which_ho_calving_front == HO_CALVING_FRONT_SUBGRID .and. &
        model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then
       where (model%calving%calving_mask == 1) model%ocean_data%deltaT_ocn = 0.0d0
    endif

    ! If inverting for flow_enhancement_factor, then update it here

    if ( model%options%which_ho_flow_enhancement_factor == HO_FLOW_ENHANCEMENT_FACTOR_INVERSION) then

       ! Given the surface elevation target, compute the thickness target.
       ! This can change in time if the bed topography is dynamic.

       call glissade_usrf_to_thck(&
            model%geometry%usrf_obs,  &
            model%geometry%topg,      &
            model%climate%eus,        &
            thck_obs)

       ! Compute f_ground_cell based on thck_obs instead of thck.
       ! This is done so that the relaxation target is based on whether the target ice
       !  (not the current ice) is grounded or floating.
       ! Note: f_flotation_obs and f_ground_obs are not used, but they
       !       are required output arguments for the subroutine.
       ! Note: This call is not needed if the target for grounded and floating ice
       !       have the same target for E.

       call glissade_grounded_fraction(&
            ewn,          nsn,             &
            parallel,                      &
            itest, jtest, rtest,           &  ! diagnostic only
            thck_obs*thk0,                 &
            model%geometry%topg*thk0,      &
            model%climate%eus*thk0,        &
            ice_mask,                      &
            floating_mask,                 &
            land_mask,                     &
            model%options%which_ho_ground, &
            model%options%which_ho_flotation_function, &
            model%options%which_ho_fground_no_glp,     &
            f_flotation_obs,               &
            f_ground_obs,                  &
            f_ground_cell_obs)

       call glissade_inversion_flow_enhancement_factor(&
!!            model%numerics%dt * tim0,                         &
            model%numerics%dt,                                &  ! s
            ewn, nsn,                                         &
!!            model%numerics%dew*len0,                          &  ! m
!!            model%numerics%dns*len0,                          &  ! m
            model%numerics%dew,                               &  ! m
            model%numerics%dns,                               &  ! m
            itest, jtest, rtest,                              &
            model%velocity%velo_sfc * vel0,                   &  ! m/s
            model%velocity%velo_sfc_obs * vel0,               &  ! m/s
            model%geometry%dthck_dt,                          &  ! m/s
            ice_mask,                                         &
            model%geometry%f_ground_cell,                     &
            f_ground_cell_obs,                                &
            model%paramets%flow_enhancement_factor_ground,    &
            model%paramets%flow_enhancement_factor_float,     &
            model%inversion%flow_enhancement_velo_scale,      &  ! m/s
            flow_enhancement_timescale,                       &  ! s
            model%inversion%flow_enhancement_thck_scale,      &  ! m
            model%inversion%flow_enhancement_length_scale,    &  ! m
            model%inversion%flow_enhancement_relax_factor,    &
            model%temper%flow_enhancement_factor)

       call parallel_halo(model%temper%flow_enhancement_factor, parallel)

    endif   ! which_ho_flow_enhancement_factor

  end subroutine glissade_inversion_solve

!***********************************************************************

  subroutine invert_basal_friction(&
       dt,                        &
       nx,            ny,         &
       dx,            dy,         &
       itest, jtest,  rtest,      &
       babc_thck_scale,           &
       babc_timescale,            &
       babc_length_scale,         &
       babc_relax_factor,         &
       friction_c_max,            &
       friction_c_min,            &
       f_ground,                  &
       stag_thck,                 &
       stag_thck_obs,             &
       stag_dthck_dt,             &
       friction_c_relax,          &
       friction_c)

    use glissade_grid_operators, only: glissade_laplacian_stagvar

    ! Compute a spatially varying basal friction field defined at cell vertices.
    ! Here, the field has the generic name 'friction_c', which could be either powerlaw_c or coulomb_c.
    ! The method is similar to that of Pollard & DeConto (TC, 2012) and is applied to all grounded ice.
    ! Adjustments are based on a thickness target:
    !    Where stag_thck > stag_thck_obs, friction_c is reduced to increase sliding.
    !    Where stag_thck < stag_thck_obs, friction_c is increased to reduce sliding.
    ! The resulting friction_c is constrained to lie within a prescribed range, [friction_c_min, friction_c_max].
    ! Note: For grounded ice with fixed bed topography, inversion based on thck is equivalent to inversion based on usrf.
    !       But for ice that is partly floating, it seems better to invert based on thck, because thck errors
    !        for shelves are much larger than usrf errors, and we do not want to underweight the errors.
    ! Note, March 2025: Rewriting the math in terms of d(logC)/dt instead of dC/dt.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    real(dp), intent(in) :: &
         dx, dy                  ! grid cell length (m)

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         babc_thck_scale,      & ! inversion inversion scale (m)
         babc_timescale,       & ! inversion timescale (s); must be > 0
         babc_length_scale,    & ! diffusive length scale (m) for inversion
         babc_relax_factor,    & ! controls strength of relaxation to default values
         friction_c_max,       & ! upper bound for friction_c (units correspond to powerlaw_c or coulomb_c)
         friction_c_min          ! lower bound for friction_c

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
         f_ground,             & ! grounded fraction at vertices, 0 to 1
         stag_thck,            & ! ice thickness at vertices (m)
         stag_thck_obs,        & ! observed ice thickness at vertices (m)
         stag_dthck_dt,        & ! rate of change of ice thickness at vertices (m/s)
         friction_c_relax        ! basal friction field to which we (optionally) relax

    real(dp), dimension(nx-1,ny-1), intent(inout) ::  &
         friction_c              ! basal friction field to be adjusted (powerlaw_c or coulomb_c)

   ! local variables

    real(dp), dimension(nx-1,ny-1) ::  &
         stag_dthck,           & ! stag_thck - stag_thck_obs
         logC,                 & ! log_10(friction_c)
         dlogC,                & ! change in log_10(friction_c)
         logC_relax,           & ! log_10(friction_c_relax)
         del2_logC               ! Laplacian term, d2(logC)/dx2 + d2(logC)/dy2

    integer, dimension(nx-1,ny-1) :: &
         del2_mask               ! mask for Laplacian smoothing

    real(dp) ::  &
         thck_target,          & ! local target for ice thickness (m)
         term_thck,            & ! tendency term based on thickness error
         term_dHdt,            & ! tendency term based on dH/dt
         term_laplacian,       & ! tendency term based on Laplacian smoothing
         term_relax              ! tendency term based on relaxation to default value

    integer :: i, j

    real(dp), parameter :: logmin = -99.d0   ! arbitrary negative value;
                                             ! values of log(c) below logmin are considered non-physical

    !TODO - Make sure logmin is not entering any calculations

    ! Compute the log (base 10) of the current friction_c field.
    ! We work with log(C) instead of C itself, because the physical effects of changing C
    ! by an amount dC are much greater at low C than at high C.
    where (friction_c > 0.0d0)
       logC = log10(friction_c)
    elsewhere
       logC = logmin
    endwhere

    ! initialize
    dlogC(:,:) = 0.0d0
    where (friction_c_relax > 0.0d0)
       logc_relax = log10(friction_c_relax)
    elsewhere
       logc_relax = logmin
    endwhere

    ! Compute the difference between the current and target thickness
    stag_dthck(:,:) = stag_thck(:,:) - stag_thck_obs(:,:)

    ! Compute the Laplacian of logc.
    ! Do this for all cells, including floating and ice-free cells,
    ! to get a smooth transition at the ice margin.
    del2_mask = 1

    call glissade_laplacian_stagvar(&
         nx,           ny,          &
         dx,           dy,          &
         logC,         del2_logC,   &
         del2_mask)

    ! Loop over vertices
    ! Note: f_ground should be computed before transport. Thus, if a vertex is grounded
    !       before transport and is fully floating afterward, friction_c is computed here.

    ! Compute the rate of change of log(C).
    ! For a thickness target H_obs, the rate is given by
    !     d(logC)/dt = (H - H_obs)/(H0*tau0) + (dH/dt)*2/H0 + del2(C)*L0^2/tau0 - r*(log(C) - log(Cr))/tau0,
    ! where tau0 = babc_timescale, H0 = babc_thck_scale, r = babc_relax_factor,
    !  C_r is a relaxation target, and L0 is a diffusive length scale.
    
    ! Without the relaxation and smoothing terms, this equation is similar to that of a damped harmonic oscillator:
    !     m * d2x/dt2 = -k*x - c*dx/dt
    ! where m is the mass, k is a spring constant, and c is a damping term.
    ! A harmonic oscillator is critically damped when c = 2*sqrt(m*k).
    !  In this case the system reaches equilibrium as quickly as possible without oscillating.
    ! Assuming unit mass (m = 1) and critical damping with k = 1/(tau^2), we obtain
    !   d2x/dt2 = -1/tau * (x/tau - 2*dx/dt)
    ! If we identify (H - H_obs)/(H0*tau) with x/tau; (2/H0)*dH/dt with 2*dx/dt; and (1/C)*dC/dt with d2x/dt2,
    !  we obtain the equation similiar to the one solved here.
    
    do j = 1, ny-1
       do i = 1, nx-1

          ! initialize terms
          term_thck = 0.0d0
          term_dHdt = 0.0d0
          term_relax = 0.0d0
          term_laplacian = 0.0d0

          if (f_ground(i,j) > 0.0d0) then  ! ice is at least partly grounded

             ! Compute tendency terms based on the thickness target
             if (babc_thck_scale > 0.0d0) then
                term_thck = -stag_dthck(i,j) / (babc_thck_scale*babc_timescale)
                term_dHdt = -stag_dthck_dt(i,j) * 2.0d0 / babc_thck_scale
             endif

          endif  ! f_ground > 0

          ! Note: If the ice is fully floating, there is no reason to compute term_thck and term_dHdt,
          !       since the ice thickness is unrelated to friction.
          !       However, we still compute Laplacian and relax terms, to give a smooth transition
          !        between grounded and floating cells. This makes the inversion more robust under
          !        grounding-line migration.

          ! Add a Laplacian smoothing term, which will discourage large curvature and make the field more linear.
          ! del2(logC) < 0 for peaks, which are lowered by this term, and del2(logC) > 0 for valleys, which are raised.
          term_laplacian = del2_logC(i,j) * babc_length_scale**2 / babc_timescale

          ! Add a term to relax C toward a target value, friction_c_relax
          if (logC(i,j) > logmin) then
             term_relax = -babc_relax_factor * (logC(i,j) - logC_relax(i,j)) / babc_timescale
          else
             term_relax = 0.0d0
          endif

          ! The remaining logic applied whether or not the cell is grounded

          ! Sum the terms
          dlogC(i,j) = (term_thck + term_dHdt + term_laplacian + term_relax) * dt

          ! Limit to prevent a large change in one step
          if (abs(dlogC(i,j)) > 0.1d0 * dt/scyr) then
             if (dlogC(i,j) > 0.0d0) then
                dlogC(i,j) =  0.1d0 * dt/scyr
             else
                dlogC(i,j) = -0.1d0 * dt/scyr
             endif
          endif

          ! Update log(C)
          logC(i,j) = logC(i,j) + dlogC(i,j)

          ! Convert log(C) back to C
          if (logC(i,j) > logmin) then
             friction_c(i,j) = 10.d0**(logC(i,j))
          else
             friction_c(i,j) = 0.0d0
          endif

          ! Limit to a physically reasonable range
          friction_c(i,j) = min(friction_c(i,j), friction_c_max)
          friction_c(i,j) = max(friction_c(i,j), friction_c_min)

          !WHL - debug
          if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
             print*, ' '
             print*, 'Increment friction_c: rank, i, j =', rtest, itest, jtest
             print*, 'dx, dy, length_scale (m)=', dx, dy, babc_length_scale
             print*, 'thck (m), thck_obs, dthck, dthck_dt (m/yr):', &
                  stag_thck(i,j), stag_thck_obs(i,j), stag_dthck(i,j), stag_dthck_dt(i,j)*scyr
             print*, 'dH term, dH/dt term, laplacian term, relax term, sum =', &
                  term_thck*dt, term_dHdt*dt, term_laplacian*dt, term_relax*dt, &
                  (term_thck + term_dHdt + term_laplacian + term_relax)*dt
             print*, 'dlogC, new friction_c =', dlogc(i,j), friction_c(i,j)
          endif

       enddo  ! i
    enddo  ! j

    ! optional diagnostics
    if (verbose_inversion) then
       call point_diag(stag_dthck, 'stag_thck - stag_thck_obs', itest, jtest, rtest, 7, 7)
       call point_diag(stag_dthck_dt*scyr, 'stag_dthck_dt (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(f_ground, 'f_ground', itest, jtest, rtest, 7, 7)
       call point_diag(del2_logc, 'del2(logC)', itest, jtest, rtest, 7, 7, '(e12.3)')
       call point_diag(logC, 'logC', itest, jtest, rtest, 7, 7)
       call point_diag(dlogc, 'dlogC', itest, jtest, rtest, 7, 7, '(e12.3)')
    endif

  end subroutine invert_basal_friction

!***********************************************************************

  subroutine glissade_inversion_deltaT_basin(&
       dt,                          &
       nx,            ny,           &
       dx,            dy,           &
       itest, jtest,  rtest,        &
       nbasin,                      &
       basin_number,                &
       thck,                        &
       dthck_dt,                    &
       floating_thck_target,        &
       deltaT_ocn_thck_scale,       &
       deltaT_ocn_timescale,        &
       deltaT_ocn_temp_scale,       &
       deltaT_basin_relax,          &
       basin_mass_correction,       &
       basin_number_mass_correction,&
       deltaT_ocn)

    use glissade_utils, only: glissade_basin_average

    ! For the case that bmlt_float is computed based on thermal_forcing,
    !  adjust deltaT_ocn at the basis scale, where deltaT_ocn is a bias corrrection
    !  or tuning parameter for the thermal forcing parameterization.
    ! In each basin, we compute the mean thickness of floating or lightly grounded ice
    !  and compare to a target thickness (usually based on observations).
    ! In basins where this ice is too thick, we increase deltaT_ocn uniformly across the basin.
    !  and where this ice is too thin, we decrease deltaT_ocn.
    ! Note: Other possible targets include the total floating area or grounded area.
    !       One reason not to use the total floating area is that the deltaT_ocn
    !        correction can become entangled with the calving scheme.
    !       One reason not to use the total grounded area is that the relative change
    !        in grounded area associated with GL advance or retreat will be very small
    !        in some basins compared to the total grounded area; also, we don't want
    !        growth of ice on beds above sea level to influence the correction.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    real(dp), intent(in) :: &
         dx, dy                  ! grid cell size in each direction (m)

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    integer, intent(in) :: &
         nbasin                  ! number of basins

    integer, dimension(nx,ny), intent(in) :: &
         basin_number            ! basin ID for each grid cell

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,             &     ! ice thickness (m)
         dthck_dt,         &     ! dH/dt (m/s)
         floating_thck_target    ! target thickness for floating ice (m);
                                 ! includes lightly grounded ice, based on basin_flotation_threshold
    real(dp), intent(in) :: &
         deltaT_ocn_thck_scale,& ! inversion thickness scale (m); must be > 0
         deltaT_ocn_timescale, & ! inversion timescale (s); must be > 0
         deltaT_ocn_temp_scale,& ! inversion temperature scale (degC)
         deltaT_basin_relax,   & ! value toward which we relax each basin (degC)
         basin_mass_correction   ! optional mass correction (Gt) for a selected basin

    integer, intent(in) :: &
         basin_number_mass_correction ! integer ID for the basin receiving the correction

    real(dp), dimension(nx,ny), intent(inout) ::  &
         deltaT_ocn              ! deltaT correction to thermal forcing in each basin (deg C)

    ! local variables

    real(dp), dimension(nbasin) :: &
         floating_thck_target_basin,      &   ! floating mean thickness target in each basin (m^3)
         floating_thck_basin,             &   ! current mean ice thickness in each basin (m)
         floating_dthck_dt_basin,         &   ! rate of change of basin mean ice thickness (m/s)
         dT_basin_dt,                     &   ! rate of change of deltaT_basin in each basin (degC/s)
         deltaT_basin                         ! basin-level averages of deltaT_ocn (degC)

    integer :: i, j
    integer :: nb      ! basin number
    real(dp) :: dthck  ! difference between modeled and observed mean thickness
    real(dp) :: term_thck, term_dHdt, term_relax  ! terms in prognostic equation for deltaT_basin
    real(dp), dimension(nx,ny) :: mask        ! temporary mask, = 1 everywhere

    ! Note: In some basins, the floating ice volume may be too small no matter how much we lower deltaT_ocn,
    !        since the basal melt rate drops to zero and can go no lower.
    !       To prevent large negative values, the deltaT_ocn correction is capped at a moderate negative value.
    !       A positive cap might not be needed but is included to be on the safe side.

    ! TODO: Make these config parameters?
    real(dp), parameter :: &
         deltaT_basin_maxval = 2.0d0,       & ! max magnitude of basin temperature correction (deg C)
         dT_basin_dt_maxval = 1.0d0/scyr      ! max allowed magnitude of d(deltaT_basin)/dt (deg/yr converted to deg/s)

    ! For each basin, compute the current and target mean ice thickness for the target region.
    ! Also compute the current rate of mean thickness change.

    call get_basin_targets(&
         nx,          ny,                     &
         dx,          dy,                     &
         nbasin,      basin_number,           &
         thck,        dthck_dt,               &
         floating_thck_target,                &
         basin_number_mass_correction,        &
         basin_mass_correction,               &
         floating_thck_target_basin,          &
         floating_thck_basin,                 &
         floating_dthck_dt_basin)

    ! Compute the current deltaT in each basin.
    ! Since deltaT_basin(nb) = deltaT_ocn for all cells in a given basin,
    !  it suffices to compute the average value of deltaT_ocn.

    mask = 1  ! do not mask out any points

    call glissade_basin_average(&
         nx,          ny,                   &
         nbasin,      basin_number,         &
         mask,                              &
         deltaT_ocn,  deltaT_basin)

    ! header for optional diagnostics
    if (verbose_inversion .and. this_rank == rtest) then
       print*, ' '
       print*, 'basin, term_thck*dt, term_dHdt*dt, term_relx*dt, new deltaT_basin:'
    endif

    ! Warm the basin where the ice is too thick, and cool where the ice is too thin.
    ! The calculation is basically the same as for local deltaT_ocn inversion,
    !  but with basin-average thickness replacing local thickness.

    do nb = 1, nbasin

       ! Compute d/dt(T_basin)
       dthck = floating_thck_basin(nb) - floating_thck_target_basin(nb)
       term_thck = (dthck/deltaT_ocn_thck_scale) * (deltaT_ocn_temp_scale/deltaT_ocn_timescale)
       term_dHdt = deltaT_ocn_temp_scale * floating_dthck_dt_basin(nb) * 2.0d0 / deltaT_ocn_thck_scale
       term_relax = -(deltaT_basin(nb) - deltaT_basin_relax) / deltaT_ocn_timescale
       dT_basin_dt(nb) = term_thck + term_dHdt + term_relax

       ! Limit dT_basin/dt to a prescribed range
       ! This prevents rapid changes in basins with small volume targets.
       dT_basin_dt(nb) = min(dT_basin_dt(nb),  dT_basin_dt_maxval)
       dT_basin_dt(nb) = max(dT_basin_dt(nb), -dT_basin_dt_maxval)

       ! Update deltaT_basin and limit to a prescribed range
       deltaT_basin(nb) = deltaT_basin(nb) + dT_basin_dt(nb) * dt
       deltaT_basin(nb) = min(deltaT_basin(nb),  deltaT_basin_maxval)
       deltaT_basin(nb) = max(deltaT_basin(nb), -deltaT_basin_maxval)

       ! deltaT_basin diagnostics
       if (verbose_inversion .and. this_rank == rtest) then
          write(6,'(i6,4f14.7)') nb, term_thck*dt, term_dHdt*dt, term_relax*dt, deltaT_basin(nb)
       endif

    enddo

    ! Increment deltaT_ocn in each grid cell.
    ! Note: deltaT_ocn is a 2D field, but here its value is uniform in each basin.
    do j = 1, ny
       do i = 1, nx
          nb = basin_number(i,j)
          if (nb >= 1 .and. nb <= nbasin) then
             deltaT_ocn(i,j) = deltaT_basin(nb)
          endif
       enddo
    enddo

  end subroutine glissade_inversion_deltaT_basin

!***********************************************************************

  subroutine glissade_inversion_deltaT_ocn(&
       dt,                       &
       nx,            ny,        &
       dx,            dy,        &
       itest, jtest,  rtest,     &
       deltaT_ocn_thck_scale,    &
       deltaT_ocn_timescale,     &
       deltaT_ocn_temp_scale,    &
       deltaT_ocn_length_scale,  &
       deltaT_ocn_relax,         &
       f_ground_cell,            &
       thck,                     &
       thck_obs,                 &
       dthck_dt,                 &
       deltaT_ocn)

    ! Compute spatially varying temperature correction factors at cell centers.
    ! Adjustments are made in floating grid cells, typically based on a thickness target:
    !    Where thck > thck_obs, deltaT_ocn is increased to increase basal melting.
    !    Where thck < thck_obs, deltaT_ocn is reduced to reduce basal melting.
    ! Note: deltaT_ocn is constrained to lie within a prescribed range,
    ! [deltaT_ocn_min, deltaT_ocn_max].

    use glissade_grid_operators, only: glissade_laplacian

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    real(dp), intent(in) :: &
         dx, dy                  ! grid cell length (m)

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         deltaT_ocn_thck_scale,&   ! inversion thickness scale (m); must be > 0
         deltaT_ocn_timescale, &   ! inversion timescale (s); must be > 0
         deltaT_ocn_temp_scale,&   ! inversion temperature scale (degC)
         deltaT_ocn_length_scale   ! diffusive length scale (m) for inversion

    real(dp), dimension(nx,ny), intent(in) ::  &
         deltaT_ocn_relax,     &   ! deltaT_ocn field toward which we relax
         f_ground_cell,        &   ! grounded fraction at cell centers, 0 to 1
         thck,                 &   ! ice thickness (m)
         thck_obs,             &   ! observed ice thickness (m)
         dthck_dt                  ! rate of change of ice thickness (m/s)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         deltaT_ocn              ! temperature correction factor (degC)

    ! local variables

    real(dp), dimension(nx,ny) ::  &
         dthck,                & ! thck - thck_obs
         del2_deltaT_ocn         ! Laplacian term, d2(dT_ocn)/dx2 + d2(dT_ocn)/dy2

    integer, dimension(nx,ny) :: &
         del2_mask               ! mask for Laplacian smoothing

    real(dp) ::  &
         term_thck,            & ! tendency term based on thickness target
         term_dHdt,            & ! tendency term based on dH/dt
         term_laplacian,       & ! tendency term based on Laplacian smoothing
         term_relax,           & ! term that relaxes deltaT_ocn toward base value
         term_sum                ! sum of the terms above

    integer :: i, j

    real(dp), parameter :: &
         deltaT_ocn_maxval = 4.0d0   ! max magnitude of the local temperature correction (deg C)

    ! Check for positive scales

    if (deltaT_ocn_thck_scale <= 0.0d0) then
       call write_log('Error, deltaT_ocn_thck_scale must be > 0', GM_FATAL)
    endif

    if (deltaT_ocn_timescale <= 0.0d0) then
       call write_log('Error, deltaT_ocn timescale must be > 0', GM_FATAL)
    endif

    ! Compute difference between current and target thickness.
    dthck(:,:) = thck(:,:) - thck_obs(:,:)

    ! Compute the Laplacian of deltaT_ocn.
    ! Note: Grounded cells are included in the mask. In these cells, we relax deltaT_ocean
    !       toward the basin average, with Laplacian smoothing to give a smooth transition
    !       between grounded and floating cells.
    where (thck > 0.0d0)
       del2_mask = 1
    elsewhere
       del2_mask = 0
    endwhere

    call glissade_laplacian(&
         nx,           ny,              &
         dx,           dy,              &
         deltaT_ocn,   del2_deltaT_ocn, &
         del2_mask)

    ! Loop over cells where f_ground_cell < 1
    ! Note: f_ground_cell should be computed before transport, so that if a cell is at least
    !       partly floating before transport and fully grounded afterward, deltaT_ocn is computed.

    do j = 1, ny
       do i = 1, nx

          ! Compute the rate of change of deltaT_ocn.
          ! For a thickness target H_obs, the rate is given by
          !     dTc/dt = -T0 * [(H - H_obs)/(H0]^n / tau0 + (dH/dt)*2/H0 + del2(dTc)*L0^2/tau0 + (T_r - T)/tau0]
          ! where Tc = deltaT_ocn, tau0 = deltaT_ocn_timescale, H0 = deltaT_ocn_thck_scale,
          !  T0 = deltaT_ocn_temp_scale, L0 = deltaT_ocn_length_scale, T_r is a relaxation target,
          !  and n is an exponent.
          ! T0 should be similar in magnitude to the max deltaT_ocn we will accept when dthck ~ H0.
          ! T0 plays a role similar to relax_factor in the inversions for Cc, Cp and E;
          !  it controls the size of the dH and dH/dt terms compared to the relaxation term.
          ! Increasing T0 makes the relaxation relatively weaker.

          ! initialize terms
          term_thck = 0.0d0
          term_dHdt = 0.0d0
          term_relax = 0.0d0
          term_laplacian = 0.0d0
          
          if (thck(i,j) > 0.0d0 .and. f_ground_cell(i,j) < 1.0d0) then  ! ice is present and at least partly floating

             term_thck = (dthck(i,j)/deltaT_ocn_thck_scale) * (deltaT_ocn_temp_scale/deltaT_ocn_timescale)
             term_dHdt = deltaT_ocn_temp_scale * dthck_dt(i,j) * 2.0d0 / deltaT_ocn_thck_scale

          endif
          
          ! Note: If the ice is fully grounded, there is no reason to compute term_thck and term_dHdt,
          !       since the ice thickness is unrelated to ocean temperature.
          !       However, we still compute Laplacian and relax terms, to give a smooth transition
          !        between grounded and floating cells. This makes the inversion more robust under
          !        grounding-line migration.

          ! Compute a Laplacian smoothing term, which will make the field more linear.
          ! del2(dT_ocn) < 0 for peaks, which are lowered by this term, and del2(dT_ocn) > 0 for valleys, which are raised.
          term_laplacian = del2_deltaT_ocn(i,j) * deltaT_ocn_length_scale**2 / deltaT_ocn_timescale

          ! Compute a term to relax C toward a target value, deltaT_ocn_relax.
          ! This could be either deltaT_ocn_basin_avg, or otherwise a default value of 0.
          term_relax = (deltaT_ocn_relax(i,j) - deltaT_ocn(i,j)) / deltaT_ocn_timescale

          ! Update deltatT_ocn
          deltaT_ocn(i,j) = deltaT_ocn(i,j) + (term_thck + term_dHdt + term_laplacian + term_relax) * dt

          ! Limit to a physically reasonable range
          deltaT_ocn(i,j) = min(deltaT_ocn(i,j),  deltaT_ocn_maxval)
          deltaT_ocn(i,j) = max(deltaT_ocn(i,j), -deltaT_ocn_maxval)

          if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
             print*, ' '
             print*, 'Increment deltaT_ocn: rank, i, j =', rtest, itest, jtest
             print*, 'thck scale (m), temp scale (degC), timescale (yr):', &
                  deltaT_ocn_thck_scale, deltaT_ocn_temp_scale, deltaT_ocn_timescale/scyr
             print*, 'thck, thck_obs, err thck (m), dthck_dt (m/yr):', &
                  thck(i,j), thck_obs(i,j), dthck(i,j), dthck_dt(i,j)*scyr
             print*, 'term_thck, term_dHdt, term_laplacian, term_relax:', &
                  term_thck*dt, term_dHdt*dt, term_laplacian*dt, term_relax*dt
             print*, 'term_sum, new dT_ocn:', &
                  (term_thck + term_dHdt + term_laplacian + term_relax)*dt, deltaT_ocn(i,j)
          endif

       enddo  ! i
    enddo  ! j

    ! optional diagnostics
    if (verbose_inversion) then
       call point_diag(f_ground_cell, 'f_ground_cell', itest, jtest, rtest, 7, 7)
       call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(dthck, 'err thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(dthck_dt*scyr, 'dthck/dt (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(deltaT_ocn, 'New deltaT_ocn (deg)', itest, jtest, rtest, 7, 7, '(f10.5)')
    endif

  end subroutine glissade_inversion_deltaT_ocn

!***********************************************************************

  subroutine glissade_inversion_flow_enhancement_factor(&
       dt,                                 &
       nx,            ny,                  &
       dx,            dy,                  &
       itest, jtest,  rtest,               &
       velo_sfc,                           &
       velo_sfc_obs,                       &
       dthck_dt,                           &
       ice_mask,                           &
       f_ground_cell,                      &
       f_ground_cell_obs,                  &
       flow_enhancement_factor_ground,     &
       flow_enhancement_factor_float,      &
       flow_enhancement_velo_scale,        &
       flow_enhancement_timescale,         &
       flow_enhancement_thck_scale,        &
       flow_enhancement_length_scale,      &
       flow_enhancement_relax_factor,      &
       flow_enhancement_factor)

    use glissade_grid_operators, only: glissade_unstagger, glissade_laplacian

    ! Compute a spatially varying field of flow enhancement factors at cell centers.
    ! This is an empirical factor, often denoted as E, that multiplies the
    !  temperature-dependent flow factor A in the equation for effective viscosity.
    ! Larger E corresponds to softer ice and faster flow.
    ! The CISM default for grounded ice is 1.0.  Higher values are typical in SIA models,
    !  and lower values are often needed to match observed speeds in ice shelves.
    ! This subroutine adjusts E based on a surface speed target:
    !    Where velo_sfc > velo_sfc_obs, E is decreased to slow the ice.
    !    Where velo_sfc < velo_sfc_obs, E is increased to speed up the ice.
    ! E is constrained to lie within a prescribed range.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    real(dp), intent(in) :: &
         dx, dy                  ! grid cell length (m)

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
         velo_sfc,             & ! surface ice speed (m/s)
         velo_sfc_obs            ! observed surface ice speed (m/s)

    real(dp), dimension(nx,ny), intent(in) ::  &
         dthck_dt,             & ! rate of change of thickness (m/s)
         f_ground_cell,        & ! grounded fraction at cell centers, based on current thck
         f_ground_cell_obs       ! grounded fraction at cell centers, based on thck_obs

    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask                ! = 1 where ice is present, else = 0

    real(dp), intent(in) :: &
         flow_enhancement_factor_ground,      & ! default flow_enhancement_factor for grounded ice
         flow_enhancement_factor_float,       & ! default flow_enhancement_factor for floating ice
         flow_enhancement_velo_scale,         & ! velocity scale for adjusting flow_enhancement_factor (m/s)
         flow_enhancement_thck_scale,         & ! thickness scale for adjusting flow_enhancement_factor (m)
         flow_enhancement_length_scale,       & ! diffusive length scale for adjusting flow_enhancement_factor (m)
         flow_enhancement_timescale,          & ! timescale for adjusting flow_enhancement_factor (s)
         flow_enhancement_relax_factor          ! controls strength of relaxation (unitless)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         flow_enhancement_factor          ! flow enhancement factor (unitless)

    ! local variables

    integer :: i, j

    real(dp), dimension(nx,ny) ::  &
         velo_sfc_cell,        & ! velo_sfc interpolated to cell centers
         velo_sfc_obs_cell,    & ! velo_sfc_obs interpolated to cell centers
         dvelo,                & ! velo_sfc_cell - velo_sfc_obs_cell
         logE,                 & ! log_10(flow_enhancement_factor)
         dlogE,                & ! change in log C
         del2_logE,            & ! Laplacian term, d2(logE)/dx2 + d2(logE)/dy2
         E_relax,              & ! value toward which E is relaxed
         logE_relax              ! log10(E_relax)

    integer, dimension(nx,ny) :: &
         del2_mask               ! mask for Laplacian smoothing

    real(dp) ::  &
         term_velo,            & ! tendency term based on speed target
         term_dHdt,            & ! tendency term based on dthck/dt
         term_laplacian,       & ! tendency term based on Laplacian smoothing
         term_relax,           & ! term that relaxes E toward a default value
         dflow_factor_dt         ! sum of these three terms

    ! TODO: Make these config parameters?
    real(dp), parameter :: &
         flow_enhancement_factor_max = 10.0d0, & ! max allowed value of flow_enhancement_factor (unitless)
         flow_enhancement_factor_min = 0.10d0    ! min allowed value of flow_enhancement_factor (unitless)

    real(dp), parameter :: logmin = -99.d0   ! arbitrary negative value;
                                             ! values of log(E) below logmin are considered non-physical

    ! Make sure E has a nonzero value in all ice-covered cells.
    ! This is needed for cells that have filled with ice since the previous call.
    ! Also, set E to a default value of flow_enhancement_factor_ground in ice-free cells.
    ! A value E = 0 is problematic if the cell later becomes ice-filled (e.g., after turning off inversion)
    ! because it leads to flwa = 0, giving NaNs in the velocity solver.

    where (flow_enhancement_factor == 0.0d0)
       where (ice_mask == 1)
          flow_enhancement_factor = f_ground_cell  * flow_enhancement_factor_ground  &
                         + (1.0d0 - f_ground_cell) * flow_enhancement_factor_float
       elsewhere (f_ground_cell_obs > 0.0d0)
          flow_enhancement_factor = flow_enhancement_factor_ground
       elsewhere   ! f_ground_cell_obs = 0; likely floating or ice-free ocean
          flow_enhancement_factor = flow_enhancement_factor_float
       endwhere
    endwhere

    ! Compute the log (base 10) of the current flow_enhancement_factor field E.
    ! We work with log(E) instead of E itself, because the physical effects of changing E
    ! by an amount dE are much greater at low E than at high E.
    where (flow_enhancement_factor > 0.0d0)
       logE = log10(flow_enhancement_factor)
    elsewhere
       logE = logmin
    endwhere

    dlogE(:,:) = 0.0d0

    ! Initialize the relaxation target
    ! This is the value we would ideally choose, depending on whether the ice is grounded or floating.
    ! Note: Relax toward the floating value for observed ice-free ocean as well as observed floating ice.
    !TODO - Make sure f_ground_cell_obs > 0 for ice-free land.
    E_relax(:,:) = f_ground_cell_obs  * flow_enhancement_factor_ground  &
        + (1.0d0 - f_ground_cell_obs) * flow_enhancement_factor_float

    where (E_relax > 0.0d0)
       logE_relax = log10(E_relax)
    elsewhere
       logE_relax = logmin
    endwhere

    ! Interpolate velo_sfc and velo_sfc_obs to cell centers
    call glissade_unstagger(&
         nx,           ny,         &
         velo_sfc,     velo_sfc_cell)

    call glissade_unstagger(&
         nx,           ny,         &
         velo_sfc_obs, velo_sfc_obs_cell)

    ! Compute difference between current and target speed
    ! Note: For ice-covered cells with ice-free targets, velo_sfc_cell > 0 and velo_sfc_obs_cell = 0.
    !       In that case, velo_sfc_cell - velo_sfc_obs_cell > 0.
    !       However, we set dvelo = 0 to prevent E from decreasing, which could slow and thicken the ice.
    where (velo_sfc_obs_cell > 0.0d0)
       dvelo = velo_sfc_cell - velo_sfc_obs_cell
    elsewhere
       dvelo = 0.0d0
    endwhere

    ! Compute the Laplacian of logE.
    ! Ignore values in ice-free cells.
    where (ice_mask == 1)
       del2_mask = 1
    elsewhere
       del2_mask = 0
    endwhere

    call glissade_laplacian(&
         nx,           ny,          &
         dx,           dy,          &
         logE,         del2_logE,   &
         del2_mask)

    ! Loop over cells where ice is present
    do j = 1, ny
       do i = 1, nx

          ! Compute the rate of change of log(E).
          ! In general, positive speed errors drive a decrease in E (so the ice becomes more viscous).
          ! For a speed target v_obs, the rate is given by
          !     dlogE/dt =  (v - v_obs)/(V0*tau0) + (dH/dt)*2/H0 + del2(logE)*L0^2/tau0 - r*ln(E/E_r)/tau0
          ! where tau = flow_enhancement_timescale, V0 = flow_enhancement_velo_scale,
          !  H0 = flow_enhancement_thck_scale, dH/dt = thickness tendency,
          !  r = flow_enhancement_relax_factor, L0 = flow_enhancement_length_scale, and E_r is a relaxation target.
          ! It may seem more natural for the second term in brackets to be (dv/dt)*2/V0.
          !  However, this leads to oscillations, since the velocity is very sensitive to E.
          ! Using dH/dt instead allows a balance between term_velo and the other terms as the flow
          !  approaches a steady state.

          ! initialize terms
          term_velo = 0.0d0
          term_dHdt = 0.0d0
          term_relax = 0.0d0
          term_laplacian = 0.0d0

          if (ice_mask(i,j) == 1) then

             if (flow_enhancement_thck_scale > 0.0d0) then
                if (flow_enhancement_velo_scale > 0.0d0) then
                   term_velo = -dvelo(i,j) / (flow_enhancement_velo_scale * flow_enhancement_timescale)
                endif
                term_dHdt = dthck_dt(i,j) * 2.0d0 / flow_enhancement_thck_scale
             endif

          endif

          ! Compute a Laplacian smoothing term, which will discourage large curvature and make the field more linear.
          ! del2(logE) < 0 for peaks, which are lowered by this term, and del2(logE) > 0 for valleys, which are raised.
          if (flow_enhancement_timescale > 0.0d0) then
             term_laplacian = del2_logE(i,j) * flow_enhancement_length_scale**2 / flow_enhancement_timescale
          endif

          ! Compute a relaxation term that nudges flow_enhancement_factor toward a target value.
          ! With flow_enhancement_relax_factor = 1, we have term_relax = -1 (or 1) when the factor is within
          !  a factor of e (or 1/e) of its target.
          if (flow_enhancement_timescale > 0.0d0 .and. logE(i,j) > logmin) then
             term_relax = flow_enhancement_relax_factor * (logE_relax(i,j) - logE(i,j)) / flow_enhancement_timescale
          endif

          ! Sum the terms
          dlogE(i,j) = (term_velo + term_dHdt + term_laplacian + term_relax) * dt

          ! Limit to prevent a large change in one step
          if (abs(dlogE(i,j)) > 0.1d0 * dt/scyr) then
             if (dlogE(i,j) > 0.0d0) then
                dlogE(i,j) =  0.1d0 * dt/scyr
             else
                dlogE(i,j) = -0.1d0 * dt/scyr
             endif
          endif

          ! Update log(E)
          logE(i,j) = logE(i,j) + dlogE(i,j)

          ! Convert log(E) back to E
          if (logE(i,j) > logmin) then
             flow_enhancement_factor(i,j) = 10.d0**(logE(i,j))
          else
             flow_enhancement_factor(i,j) = 0.0d0
          endif

          ! Limit to a physically reasonable range
          flow_enhancement_factor(i,j) = min(flow_enhancement_factor(i,j), flow_enhancement_factor_max)
          flow_enhancement_factor(i,j) = max(flow_enhancement_factor(i,j), flow_enhancement_factor_min)

          if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
             print*, ' '
             print*, 'Increment flow_enhancement_factor: rank, i, j =', rtest, itest, jtest
             print*, 'dx, dy, length_scale =', dx, dy, flow_enhancement_length_scale
             print*, 'velo scale (m/yr), timescale (yr):', &
                  flow_enhancement_velo_scale*scyr, flow_enhancement_timescale/scyr
             print*, 'velo_sfc_cell, velo_sfc_obs, dvelo, dH_dt (m/yr):', &
                  velo_sfc_cell(i,j)*scyr, velo_sfc_obs_cell(i,j)*scyr, dvelo(i,j)*scyr, dthck_dt(i,j)*scyr
             print*, 'init flow enhancement factor =', flow_enhancement_factor(i,j)
             print*, 'dvelo term, dthck/dt term, laplacian term, relax term, sum =', &
                  term_velo*dt, term_dHdt*dt, term_laplacian*dt, term_relax*dt, &
                  (term_velo + term_dHdt + term_laplacian + term_relax)*dt
             print*, 'dlogE, new E =', dlogE(i,j), flow_enhancement_factor(i,j)
          endif

       enddo  ! i
    enddo  ! j

    ! optional diagnostics
    if (verbose_inversion) then
       call point_diag(f_ground_cell, 'Invert for E, f_ground_cell', itest, jtest, rtest, 7, 7)
       call point_diag(f_ground_cell_obs, 'f_ground_cell_obs', itest, jtest, rtest, 7, 7)
       call point_diag(velo_sfc_obs_cell*scyr, 'velo_sfc_obs (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(velo_sfc_cell*scyr, 'velo_sfc (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(dvelo*scyr, 'dvelo (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(dthck_dt*scyr, 'dthck_dt (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(del2_logE, 'del2(logE)', itest, jtest, rtest, 7, 7, '(e12.3)')
       call point_diag(logE, 'logE', itest, jtest, rtest, 7, 7)
       call point_diag(dlogE, 'dlogE', itest, jtest, rtest, 7, 7, '(e12.3)')
       call point_diag(flow_enhancement_factor, 'New flow_enhancement_factor E', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_inversion_flow_enhancement_factor

!***********************************************************************

  subroutine get_basin_targets(&
         nx,          ny,                     &
         dx,          dy,                     &
         nbasin,      basin_number,           &
         thck,        dthck_dt,               &
         floating_thck_target,                &
         basin_number_mass_correction,        &
         basin_mass_correction,               &
         floating_thck_target_basin,          &
         floating_thck_basin,                 &
         floating_dthck_dt_basin)

    ! For each basin, compute the current ice area and volume and the target ice area and volume
    ! of cells included in the floating thickness target.
    ! Derive the current and target mean ice thickness for each basin, along with the
    ! current rate of change.
    ! Optionally, the volume target in a single basin can be adjusted relative to observations.

    use glissade_utils, only: glissade_basin_sum

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    real(dp), intent(in) :: &
         dx, dy                  ! grid cell size in each direction (m)

    integer, intent(in) :: &
         nbasin                  ! number of basins

    integer, dimension(nx,ny), intent(in) :: &
         basin_number            ! basin ID for each grid cell

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,             &     ! ice thickness (m)
         dthck_dt,         &     ! dH/dt (m/s)
         floating_thck_target    ! target thickness for floating ice (m)

    real(dp), intent(in) :: &
         basin_mass_correction   ! optional mass correction (Gt) for a selected basin

    integer, intent(in) :: &
         basin_number_mass_correction ! integer ID for the basin receiving the correction

    real(dp), dimension(nbasin), intent(out) :: &
         floating_thck_target_basin,      & ! floating mean thickness target in each basin (m^3)
         floating_thck_basin,             & ! current mean ice thickness in each basin (m)
         floating_dthck_dt_basin            ! rate of change of basin mean ice thickness (m/s)


    ! local variables

    real(dp), dimension(nbasin) :: &
         floating_area_target_basin,      & ! floating ice area target in each basin (m^3)
         floating_volume_target_basin,    & ! floating ice volume target in each basin (m^3)
         floating_volume_basin,           & ! current floating ice volume in each basin (m^3)
         floating_dvolume_dt_basin          ! rate of change of basin volume (m^3/s)

    real(dp), dimension(nx,ny) ::  &
         floating_target_rmask,            &! real mask, = 1.0 where floating_thck_target > 0, else = 0.0
         cell_area                          ! area of grid cells (m^2)

    integer :: nb

    cell_area(:,:) = dx*dy

    ! Compute a mask for cells with a nonzero floating ice target

    where (floating_thck_target > 0.0d0)
       floating_target_rmask = 1.0d0
    elsewhere
       floating_target_rmask = 0.0d0
    endwhere

    ! For each basin, compute the area of the cells with floating_target_rmask = 1.

    call glissade_basin_sum(nx,         ny,                &
         nbasin,     basin_number,      &
         floating_target_rmask,         &
         cell_area,                     &
         floating_area_target_basin)

    ! For each basin, compute the target total ice volume in cells with floating_target_rmask = 1.
    ! Note: We could compute floating_volume_target_basin just once and write it to restart,
    !       but it is easy enough to recompute here.

    call glissade_basin_sum(&
         nx,         ny,                &
         nbasin,     basin_number,      &
         floating_target_rmask,         &
         floating_thck_target*dx*dy,    &
         floating_volume_target_basin)

    ! For each basin, compute the current total ice volume in cells with floating_target_rmask = 1.

    call glissade_basin_sum(&
         nx,         ny,                &
         nbasin,     basin_number,      &
         floating_target_rmask,         &
         thck*dx*dy,                    &
         floating_volume_basin)

    ! For each basin, compute the rate of change of the current volume in cells with floating_target_rmask = 1.

    call glissade_basin_sum(&
         nx,         ny,                &
         nbasin,     basin_number,      &
         floating_target_rmask,         &
         dthck_dt*dx*dy,                &
         floating_dvolume_dt_basin)

    ! Optionally, apply a correction to the ice volume target in a selected basin.
    ! Note: This option could in principle be applied to multiple basins, but currently is supported for one basin only.
    !       In practice, this basin is likely to be the Amundsen Sea Embayment (ISMIP6 basin #9).

    if (abs(basin_mass_correction) > 0.0d0 .and. basin_number_mass_correction > 0) then

       nb = basin_number_mass_correction
       floating_volume_target_basin(nb) = floating_volume_target_basin(nb) + &
            basin_mass_correction * (1.0d12/rhoi)   ! Gt converted to m^3
       if (verbose_inversion .and. main_task) then
          print*, ' '
          print*, 'Basin with mass correction:', basin_number_mass_correction
          print*, 'mass correction (Gt)     =', basin_mass_correction
          print*, 'volume correction (km^3) =', basin_mass_correction * (1.0d3/rhoi)
          print*, 'New volume target (km^3) =', floating_volume_target_basin(nb) / 1.0d9
       endif
    endif   ! basin_mass correction

    ! For each basin, compute the current and target mean ice thickness, and the rate of change of mean ice thickness.
    where (floating_area_target_basin > 0.0d0)
       floating_thck_target_basin = floating_volume_target_basin / floating_area_target_basin
       floating_thck_basin = floating_volume_basin / floating_area_target_basin
       floating_dthck_dt_basin = floating_dvolume_dt_basin / floating_area_target_basin
    elsewhere
       floating_thck_target_basin = 0.0d0
       floating_thck_basin = 0.0d0
       floating_dthck_dt_basin = 0.0d0
    endwhere

    if (verbose_inversion .and. main_task) then
       print*, ' '
       print*, 'basin, area target (km^2), vol target (km^3), mean H target (m):'
       do nb = 1, nbasin
          write(6,'(i6,3f12.3)') nb, floating_area_target_basin(nb)/1.d6, &
               floating_volume_target_basin(nb)/1.d9, floating_thck_target_basin(nb)
       enddo
       print*, ' '
       print*, 'basin, mean thickness (m), thickness diff (m), dthck_dt (m/yr):'
       do nb = 1, nbasin
          write(6,'(i6,3f12.3)') nb, floating_thck_basin(nb), &
               (floating_thck_basin(nb) - floating_thck_target_basin(nb)), &
               floating_dthck_dt_basin(nb)*scyr
       enddo
    endif

  end subroutine get_basin_targets

!=======================================================================

end module glissade_inversion

!=======================================================================
