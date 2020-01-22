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

  use glimmer_physcon, only: scyr
  use glimmer_paramets, only: thk0, len0
  use glimmer_log
  use glide_types
  use glide_thck, only: glide_calclsrf
  use parallel

  implicit none

  private
  public :: verbose_inversion, glissade_init_inversion, &
            glissade_inversion_bmlt_float,      &
            glissade_inversion_basal_friction,  &
            glissade_inversion_bmlt_basin

  !-----------------------------------------------------------------------------
  ! Subroutines to invert for basal fields (including basal friction beneath
  ! grounded ice and basal melting beneath floating ice) by relaxing toward
  ! a target ice thickness field.
  !-----------------------------------------------------------------------------

!!    logical, parameter :: verbose_inversion = .false.
    logical, parameter :: verbose_inversion = .true.

    real(dp), parameter :: eps08 = 1.0d-08  ! small number

!***********************************************************************

contains

!***********************************************************************

  subroutine glissade_init_inversion(model)

    ! Initialize inversion for fields of basal friction and basal melting
    ! Should be called after usrf and thck have been input and (possibly) modified by initial calving

    use glissade_masks, only: glissade_get_masks
    use glissade_calving, only: glissade_marine_connection_mask  !TODO - Move to mask module
    use glissade_bmlt_float, only: basin_sum
    use glissade_grounding_line, only: glissade_grounded_fraction

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    integer :: i, j
    integer :: nb                   ! basin number
    integer :: itest, jtest, rtest  ! local diagnostic point

    real(dp) :: var_maxval          ! max value of a given real variable; = 0.0 if not yet read in
    integer :: var_maxval_int       ! max value of a given integer variable; = 0 if not yet read in

    character(len=100) :: message

    integer, dimension(model%general%ewn, model%general%nsn) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         ocean_mask,           & ! = 1 where topg is below sea level and ice is absent
         land_mask,            & ! = 1 where topg is at or above sea level
         marine_connection_mask  ! = 1 for cells with a marine connection to the ocean

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck_obs                ! observed ice thickness, derived from usrf_obs and topg

    real(dp) :: h_cavity                     ! cavity thickness (m); < 0 for grounded ice
    real(dp) :: h_obs, h_flotation, h_buff   ! thck_obs, flotation thickness, and thck_flotation_buffer scaled to m
    real(dp) :: dh                           ! h_obs - f_flotation
    real(dp) :: dh_decimal                   ! decimal part remaining after subtracting the truncation of dh

    integer :: ewn, nsn

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
    ! If inverting for Cp or bmlt_float, then set the target elevation, usrf_obs.
    !----------------------------------------------------------------------

    if (model%options%which_ho_cp_inversion == HO_CP_INVERSION_COMPUTE .or.  &
        model%options%which_ho_bmlt_inversion == HO_BMLT_INVERSION_COMPUTE) then

       ! We are inverting for usrf_obs, so check whether it has been read in already.
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
       call usrf_to_thck(model%geometry%usrf_obs,  &
                         model%geometry%topg,      &
                         model%climate%eus,        &
                         thck_obs)

       ! Optionally, adjust the initial thickness and then reset usrf_obs.

       if (model%options%is_restart == RESTART_FALSE) then

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

          model%inversion%thck_threshold = max(model%inversion%thck_threshold, model%numerics%thklim)
          where (thck_obs <= model%inversion%thck_threshold)
             thck_obs = 0.0d0
          endwhere

          ! Set thck to be consistent with thck_obs
          model%geometry%thck = thck_obs

          ! Reset usrf_obs to be consistent with thck_obs.
          ! (usrf itself will be recomputed later in glissade_initialise)
          call thck_to_usrf(thck_obs,  &
                            model%geometry%topg,      &
                            model%climate%eus,        &
                            model%geometry%usrf_obs)

       endif   ! not a restart

       call parallel_halo(model%geometry%usrf_obs)
       call parallel_halo(thck_obs)

    endif  ! which_ho_cp_inversion or which_ho_bmlt_inversion

    ! Set masks that are used below
    call glissade_get_masks(ewn,                 nsn,                   &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   model%numerics%thklim, &
                            ice_mask,                                   &
                            floating_mask = floating_mask,              &
                            ocean_mask = ocean_mask,                    &
                            land_mask = land_mask)

    !----------------------------------------------------------------------
    ! computations specific to bmlt_float inversion
    !----------------------------------------------------------------------

    if (model%options%which_ho_bmlt_inversion == HO_BMLT_INVERSION_COMPUTE) then

       if (model%options%is_restart == RESTART_FALSE) then

          ! At the start of the run, compute a mask that determines where bmlt_float_inversion can be applied.
          ! This mask includes only cells with a path to the ocean through marine-based cells.

          !TODO - Would need to call this subroutine repeatedly if topography is changing at runtime.
          call glissade_marine_connection_mask(ewn,          nsn,            &
                                               itest, jtest, rtest,          &
                                               ocean_mask,   land_mask,      &
                                               model%inversion%bmlt_float_inversion_mask)

       endif  ! not a restart

       ! Note: There is no initialization of bmlt_float_save.
       ! If restarting, it should have been read in already.
       ! If not restarting, it will have been set to zero, which is an appropriate initial value.

       call parallel_halo(model%inversion%bmlt_float_inversion_mask)
       call parallel_halo(model%inversion%bmlt_float_save)

    endif   ! which_ho_bmlt_inversion

    !----------------------------------------------------------------------
    ! computations specific to powerlaw_c (= Cp) inversion
    !----------------------------------------------------------------------

    if (model%options%which_ho_cp_inversion == HO_CP_INVERSION_COMPUTE) then

       ! initialize powerlaw_inversion, if not already read in
       var_maxval = maxval(model%inversion%powerlaw_c_inversion)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; powerlaw_c_inversion has been read in already (e.g., when restarting)
       else
          ! initialize to a uniform value (which can be set in the config file)
          model%inversion%powerlaw_c_inversion(:,:) = model%basal_physics%powerlaw_c
       endif  ! var_maxval > 0

       if (verbose_inversion .and. this_rank == rtest) then
          print*, ' '
          print*, 'glissade_init_inversion: powerlaw_c_inversion:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.1)',advance='no') model%inversion%powerlaw_c_inversion(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    endif   ! which_ho_cp_inversion

    !----------------------------------------------------------------------
    ! computations specific to inversion of deltaT_basin
    !----------------------------------------------------------------------

    if (model%options%which_ho_bmlt_basin_inversion == HO_BMLT_BASIN_INVERSION_COMPUTE) then

       if (model%options%is_restart == RESTART_FALSE) then

          ! Compute a mask that determines where thermal forcing can be applied.
          ! This mask includes only cells with a path to the ocean through marine-based cells.
          ! Cells without such a connection are not given a thickness target.

          call glissade_marine_connection_mask(ewn,          nsn,            &
                                               itest, jtest, rtest,          &
                                               ocean_mask,   land_mask,      &
                                               marine_connection_mask)

          ! Set floating_thck_target to the thickness of lightly floating and lightly grounded ice.
          ! Here, "lightly" means that the absolute value of f_flotation = (-topg - eus) - (rhoi/rhoo)*thck
          !  is less than a prescribed threshold.
          ! Thus we include both ice that is floating but might ground (leading to
          !  a positive volume bias that will be corrected with ocean warming) and ice
          !  that is grounded but might float (leading to a negative volume bias
          !  that will be corrected with ocean cooling).

          do j = 1, nsn
             do i = 1, ewn
                h_cavity = (-(model%geometry%topg(i,j) - model%climate%eus)  &
                            - (rhoi/rhoo)*model%geometry%thck(i,j)) * thk0    ! h_cavity < 0 for grounded ice
                if (model%geometry%thck(i,j) > 0.0d0 .and. marine_connection_mask(i,j) == 1 .and. &
                     abs(h_cavity) < model%inversion%bmlt_basin_cavity_threshold) then
                   model%inversion%floating_thck_target(i,j) = model%geometry%thck(i,j)
                else
                   model%inversion%floating_thck_target(i,j) = 0.0d0
                endif
             enddo
          enddo

          if (verbose_inversion .and. this_rank == rtest) then
             print*, 'bmlt_basin_cavity_threshold =', model%inversion%bmlt_basin_cavity_threshold
             print*, ' '
             print*, 'After init_inversion, floating_thck_target (m):'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%inversion%floating_thck_target(i,j)*thk0
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'thck (m):'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'h_cavity (m):'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') &
                        (-(model%geometry%topg(i,j) - model%climate%eus)  &
                        - (rhoi/rhoo)*model%geometry%thck(i,j)) * thk0
                enddo
                write(6,*) ' '
             enddo
          endif   ! verbose

       endif   ! not a restart

       call parallel_halo(model%inversion%floating_thck_target)

    endif  ! which_ho_bmlt_basin_inversion

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'After init_inversion, usrf_obs (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%usrf_obs(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'After init_inversion, thck_obs (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_obs(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_init_inversion

!***********************************************************************

  subroutine glissade_inversion_bmlt_float(model,               &
                                           thck_new_unscaled,   &
                                           ice_mask,            &
                                           floating_mask)

    use glimmer_paramets, only: tim0, thk0
    use glimmer_physcon, only: scyr

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    real(dp), dimension(model%general%ewn, model%general%nsn), intent(in) ::   &
         thck_new_unscaled       ! ice thickness expected after mass balance, without applying bmlt_float_inversion (m)

    !Note: These masks are not part of the model derived type, and they are computed before transport
    !      based on the old ice thickness, so they cannot be computed here.
    !TODO - Make these masks part of the model derived type, so they do not need to be passed in?

    integer, dimension(model%general%ewn, model%general%nsn), intent(in) ::   &
         ice_mask,             & ! = 1 if thck > 0, else = 0
         floating_mask           ! = 1 where ice is present and floating, else = 0

    ! --- Local variables ---

    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
         topg_unscaled,        & ! bedrock topography (m)
         bmlt_float_new,       & ! newly computed value of bmlt_float, per unit grid cell area (m/s)
         dthck_dt_inversion,   & ! newly computed value of dthck_dt (m/s)
         bmlt_weight,          & ! weighting factor that reduces bmlt_float in partly grounded cells and shallow cavities
         thck_obs,             & ! observed ice thickness, derived from usrf_obs and topg
         thck_projected          ! projected thickness after appyling bmlt_float_save * bmlt_weight

    real(dp) ::  &
         nudging_factor,       & ! factor in range [0,1], used for inversion of bmlt_float
         weaning_time            ! time since the start of weaning (numerics%time - inversion%wean_tstart)

    integer :: i, j
    integer :: ewn, nsn
    integer :: itest, jtest, rtest

    logical :: first_time = .true.

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

    ! Compute a weighting factor that reduces the applied basal melting in partly or fully grounded cells.

    if (model%options%which_ho_ground == HO_GROUND_GLP_DELUXE .and.  &
        model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_FLOATING_FRAC) then

       ! Note: weight_float_cell = 1 - f_ground_cell in the limit H_cavity >> h0_cavity (or h0_cavity = 0)
       !       Otherwise, weight_float_cell is reduced for shallow cavities,
       !        to represent the difficulty for warm ocean water to enter the cavity.

       bmlt_weight = model%geometry%weight_float_cell

    else

       ! Compute a weighting factor proportional to the floating cell fraction.
       ! Depending on which_ho_ground and which_ho_ground_bmlt, basal melting may or may not be allowed
       !  in partly grounded cells.

       call get_float_fraction_factor(&
            model%options%which_ho_ground,       &
            model%options%which_ho_ground_bmlt,  &
            ice_mask,                            &
            floating_mask,                       &
            model%geometry%f_ground_cell,        &
            bmlt_weight)

    endif

    if (model%options%which_ho_bmlt_inversion == HO_BMLT_INVERSION_COMPUTE) then

       ! Invert for bmlt_float, adjusting the melt rate to relax toward the observed thickness.
       ! Note: Other kinds of sub-shelf basal melting are handled in subroutine glissade_bmlt_float_solve.
       !       Inversion is done here, after transport, when there is an updated ice thickness.
       !       Then bmlt_float_inversion is added to the previously computed bmlt.
       ! Note: Typically, whichbmlt_float = 0 when doing a model spin-up with inversion.
       !       However, we might want to add an anomaly to fields already computed by inversion.

       ! Compute the time scale for nudging.
       ! The idea (for now) is that nudging is associated with a timescale.
       ! With strong nudging we have a short timescale, ~ 1 yr, given by inversion%bmlt_timescale.
       ! With short nudging we have a long timescale, 100+ yr.
       ! Here we compute a nudging factor between 0 and 1.
       ! Then we divide bmlt_timescale by nudging_factor to get the timescale used during this timestep.
       ! Notes:
       ! * model%numerics%time = time in years since start of run
       ! * nudging_factor = 1 from the start of the run until inversion%wean_bmlt_float_tstart.
       ! * Then nudging_factor falls off until we reach nudging_factor_min, which is a floor for nudging.
       ! * If t > wean_bmlt_float_tend, we stop nudging entirely.

       !TODO - Do away with nudging_factor, and work directly with bmlt_timescale.

       if (model%inversion%wean_bmlt_float_tend > 0.0d0) then
          nudging_factor = 1.0d0  ! full nudging at start of run
       else
          nudging_factor = 0.0d0  ! no nudging if wean_bmlt_float_tend = 0
       endif

       if (model%inversion%wean_bmlt_float_tend > 0.0d0 .and. &
            model%numerics%time >= model%inversion%wean_bmlt_float_tstart) then
          if (model%numerics%time < model%inversion%wean_bmlt_float_tend) then
             weaning_time = model%numerics%time - model%inversion%wean_bmlt_float_tstart
             ! exponentially weighted nudging commented out.
!!                nudging_factor = exp(-weaning_time / model%inversion%wean_bmlt_float_timescale)
             ! Let nudging_factor fall off as 1/weaning_time.  As a result, bmlt_timescale will increase
             ! in proportion to weaning_time: by 1 yr for every 10 model years.
             ! The increase in bmlt_timescale is faster with this scaling than with exponential scaling.
             !TODO - Make the hardwired constant of 10 a config parameter.
             nudging_factor = 10.d0 / weaning_time  ! Make 10 = bmlt_timescale multiplier?
             nudging_factor = min(nudging_factor, 1.0d0)
             ! Optionally, do not allow the nudging factor (if > 0) to fall below a prescribed minimum value.
             ! This allows us to exclude nudging that is so small as to have virtually no effect.
             nudging_factor = max(nudging_factor, model%inversion%nudging_factor_min)
          else
             nudging_factor = 0.0d0
          endif
       endif

       if (verbose_inversion .and. this_rank == rtest) then
          print*, ' '
          print*, 'tstep_count, time, bmlt_float nudging_factor =', &
               model%numerics%tstep_count, model%numerics%time, nudging_factor
       endif

       ! Compute the new thickness, assuming application of bmlt_float_save * bmlt_weight.
       ! The correction to bmlt_float_save is based on the difference between this projected thickness
       !  and the observational target.

       thck_projected = thck_new_unscaled  &
                      - (model%inversion%bmlt_float_save * bmlt_weight * model%numerics%dt*tim0)

       ! thickness tendency dH/dt from one step to the next (m/s)
       ! Note: model%inversion%thck_save is in the restart file as needed for exact restart.

       if (verbose_inversion .and. this_rank == rtest) then
          print*, 'time, tstart:', model%numerics%time, model%numerics%tstart
       endif

!!       if (first_time .and. model%options%is_restart == RESTART_FALSE) then
       if (first_time) then

          first_time = .false.

          if (model%options%is_restart == RESTART_TRUE) then
             dthck_dt_inversion = (thck_projected - model%inversion%thck_save) / (model%numerics%dt * tim0)

             if (verbose_inversion .and. this_rank == rtest) then
                i = itest
                j = jtest
                print*, 'Compute dH/dt for inversion; restart = T'
                print*, 'rank, i, j, H_proj, H_proj_save, dH/dt (m/yr):', &
                     this_rank, i, j, thck_projected(i,j), model%inversion%thck_save(i,j), dthck_dt_inversion(i,j)*scyr
                print*, 'bmlt_weight:', bmlt_weight(i,j)
             endif
          else
             dthck_dt_inversion = 0.0d0  ! default to 0 on first step of the run
          endif

       else

          dthck_dt_inversion = (thck_projected - model%inversion%thck_save) / (model%numerics%dt * tim0)

          if (verbose_inversion .and. this_rank == rtest) then
             i = itest
             j = jtest
             print*, 'Compute dH/dt for inversion'
             print*, 'rank, i, j, H_proj, H_proj_save, dH/dt (m/yr):', &
                  this_rank, i, j, thck_projected(i,j), model%inversion%thck_save(i,j), dthck_dt_inversion(i,j)*scyr
             print*, 'bmlt_weight:', bmlt_weight(i,j)
          endif

       endif   ! first time

       ! Given the surface elevation target, compute the thickness target.
       ! (This can change in time if the bed topography is dynamic.)

       call usrf_to_thck(model%geometry%usrf_obs,  &
                         model%geometry%topg,      &
                         model%climate%eus,        &
                         thck_obs)

       ! Compute the new value of bmlt_float

       if (nudging_factor > eps08) then

          call invert_bmlt_float(model%numerics%dt * tim0,                      &    ! s
                                 ewn,               nsn,                        &
                                 itest,   jtest,    rtest,                      &
                                 thck_projected,                                &    ! m
                                 thck_obs*thk0,                                 &    ! m
                                 model%geometry%topg*thk0,                      &    ! m
                                 model%climate%eus*thk0,                        &    ! m
                                 ice_mask,                                      &
                                 dthck_dt_inversion,                            &    ! m/s
                                 model%inversion%bmlt_float_inversion_mask,     &
                                 model%inversion%thck_flotation_buffer*thk0,    &    ! m
                                 model%inversion%bmlt_timescale/nudging_factor, &    ! s
                                 model%inversion%bmlt_float_save,               &    ! m/s
                                 bmlt_weight,                                   &    ! [0,1]
                                 bmlt_float_new)                                     ! m/s

          ! Limit bmlt_float_new to physically reasonable values.
          ! Typically, bmlt_max_melt is greater in magnitude than bmlt_max_freeze.
          ! Note: These parameters have been scaled to have units of m/s.
          !       They are ignored if equal to zero.

          if (model%inversion%bmlt_max_melt*scyr > eps08) then

             bmlt_float_new = min (bmlt_float_new, model%inversion%bmlt_max_melt)

            !WHL - This formula will give a smoother transition to the max rate.
!            bmlt_float_new = min (bmlt_float_new, &
!                 model%inversion%bmlt_max_melt * (1.0d0 - exp(-bmlt_float_new/model%inversion%bmlt_max_melt)))

          endif  ! bmlt_max_melt

          if (model%inversion%bmlt_max_freeze*scyr > eps08) then

             bmlt_float_new = max (bmlt_float_new, -model%inversion%bmlt_max_freeze)

             !WHL - This formula will give a smoother transition to the max rate.
!            bmlt_float_new = max (bmlt_float_new &
!                 -model%inversion%bmlt_max_freeze * (1.0d0 - exp(bmlt_float_new/model%inversion%bmlt_max_freeze)))

          endif  ! bmlt_max_freeze

          ! save the value just computed
          ! This value represents a melting potential based on ocean conditions at the lower ice surface.
          ! The applied melting (per unit grid cell area) is reduced where cavities are shallow and/or
          !  ice is partly grounded.

          model%inversion%bmlt_float_save = bmlt_float_new

          model%inversion%bmlt_float_inversion = bmlt_float_new * bmlt_weight

          model%inversion%thck_save = thck_new_unscaled  &
                                   - (model%inversion%bmlt_float_inversion * model%numerics%dt*tim0)

          if (verbose_inversion .and. this_rank == rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'Inverting for bmlt_float: rank, i, j =', rtest, i, j
             print*, 'thck_projected (m), thck_obs (m):', &
                  thck_projected(i,j), thck_obs(i,j)*thk0
             print*, 'bmlt_float (per floating area), bmlt_float (per cell area), nudging factor:', &
                  model%inversion%bmlt_float_save(i,j)*scyr, &
                  model%inversion%bmlt_float_inversion(i,j)*scyr, nudging_factor
             print*, ' '
             print*, 'Inversion, f_ground_cell:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.4)',advance='no') model%geometry%f_ground_cell(i,j)
                enddo
             write(6,*) ' '
             enddo
             print*, ' '
             print*, 'dH_dt_inversion (m/yr)'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.4)',advance='no') dthck_dt_inversion(i,j)*scyr
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'bmlt_weight:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.4)',advance='no') bmlt_weight(i,j)
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'saved bmlt_float before weighting:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.4)',advance='no') model%inversion%bmlt_float_save(i,j)*scyr
                enddo
                write(6,*) ' '
             enddo
          endif ! verbose_inversion

       else  ! no nudging

          ! Note: To hold the inverted bmlt_float fixed, we would typically set which_ho_inversion = HO_INVERSION_APPLY.
          !       Alternatively, if running with which_ho_inversion = HO_INVERSION_COMPUTE,
          !        nudging is turned off when time > wean_bmlt_float_tend.

          if (verbose_inversion .and. main_task) print*, 'Apply saved value of bmlt_float inversion'

          model%inversion%bmlt_float_inversion = model%inversion%bmlt_float_save * bmlt_weight

       endif   ! nudging is turned on

    elseif (model%options%which_ho_bmlt_inversion == HO_BMLT_INVERSION_APPLY) then

       if (verbose_inversion .and. main_task) print*, 'Apply saved value of bmlt_float inversion'

       model%inversion%bmlt_float_inversion = model%inversion%bmlt_float_save * bmlt_weight

    endif   ! which_ho_inversion

    if (verbose_inversion .and. this_rank == rtest) then
       print*, ' '
       print*, 'new bmlt_float_inversion (m/yr)'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') model%inversion%bmlt_float_inversion(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_inversion_bmlt_float

!***********************************************************************

  subroutine get_float_fraction_factor(&
       which_ho_ground,                  &
       which_ho_ground_bmlt,             &
       ice_mask,                         &
       floating_mask,                    &
       f_ground_cell,                    &
       float_fraction_factor)

    ! Based on the grounding-line options and the floating mask or fraction field,
    !  compute a weighting factor for reducing bmlt_float in cells containing the GL.

    integer, intent(in) :: &
         which_ho_ground,      & ! option to use GLP for vertices and/or cells
         which_ho_ground_bmlt    ! determines which cells can have nonzero bmlt_float

    integer, dimension(:,:), intent(in) :: &
         ice_mask,             & ! = 1 where ice is present (thk > 0), else = 0
         floating_mask           ! = 1 where ice is present and floating, else = 0

    real(dp), dimension(:,:), intent(in) ::  &
         f_ground_cell           ! grounded fraction of grid cell, in range [0,1]

    real(dp), dimension(:,:), intent(out) ::  &
         float_fraction_factor   ! fraction of the cell where basal melting/freezing is allowed

    ! initialize
    float_fraction_factor = 0.0d0

    ! Compute float_fraction_factor based on GL options

    if (which_ho_ground == HO_GROUND_GLP_DELUXE) then

       if (which_ho_ground_bmlt == HO_GROUND_BMLT_FLOATING_FRAC) then

          where (ice_mask == 1 .and. f_ground_cell <  1.0d0 - eps08)
             float_fraction_factor = 1.0d0 - f_ground_cell
          endwhere

       elseif (which_ho_ground_bmlt == HO_GROUND_BMLT_ZERO_GROUNDED) then

          where (ice_mask == 1 .and. f_ground_cell < eps08)
             float_fraction_factor = 1.0d0
          endwhere

       elseif (which_ho_ground_bmlt == HO_GROUND_BMLT_NO_GLP) then

          where (floating_mask == 1)
             float_fraction_factor = 1.0d0
          endwhere

       endif   ! which_ho_ground_bmlt

    else  ! other HO_GROUND_GLP options

       where (floating_mask == 1)
          float_fraction_factor = 1.0d0
       endwhere

    endif  ! which_ho_ground

  end subroutine get_float_fraction_factor

!***********************************************************************

  subroutine invert_bmlt_float(dt,                           &
                               nx,            ny,            &
                               itest, jtest,  rtest,         &
                               thck_projected,               &
                               thck_obs,                     &
                               topg,                         &
                               eus,                          &
                               ice_mask,                     &
                               dthck_dt,                     &
                               bmlt_float_inversion_mask,    &
                               thck_flotation_buffer,        &
                               bmlt_timescale,               &
                               bmlt_float_save,              &
                               bmlt_weight,                  &
                               bmlt_float_new)

    ! Compute spatially varying bmlt_float by inversion.
    ! Apply a melt/freezing rate that will restore the ice in floating grid cells
    !  (and grounding-line adjacent grid cells) to the target surface elevation.
    ! Note: bmlt_float_inversion is defined as positive for melting, negative for freezing.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    ! Note: thck and usrf should be the expected values after applying the mass balance
    !       (although the mass balance may not yet have been applied)
    real(dp), dimension(nx,ny), intent(in) ::  &
         thck_projected,       & ! ice thickness (m) expected after mass balance, before adjusting bmlt_float_inversion
         thck_obs,             & ! observed ice thickness target (m)
         topg                    ! bedrock topography (m) (diagnostic only)

    real(dp), intent(in) :: &
         eus                     ! eustatic sea level (m) (diagnostic only)

   ! Note: When this subroutine is called, ice_mask = 1 where thck > 0, not thck > thklim.
    integer, dimension(nx,ny), intent(in) ::  &
         bmlt_float_inversion_mask, & ! = 1 for cells where bmlt_float is potentially computed and applied, else = 0;
                                      ! computed at startup based on bed topography; must have a marine connection to ocean
         ice_mask                     ! = 1 where ice is present, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         dthck_dt                ! rate of change of ice thickness (m/s) in previous timestep

    real(dp), intent(in) ::  &
         thck_flotation_buffer,& ! buffer thickness (m) to prevent thck very close to thck_flotation
         bmlt_timescale          ! timescale (s) for relaxing toward observations by changing bmlt_float

    real(dp), dimension(nx,ny), intent(in) ::  &
         bmlt_float_save,      & ! previous value of bmlt_float, before weighting by f_ground_cell (m/s)
         bmlt_weight             ! weighting factor for cells that are partly grounded or have shallow cavities, in range [0,1]

    real(dp), dimension(nx,ny), intent(out) ::  &
         bmlt_float_new          ! new value of bmlt_float (m/s), based on relaxation to observed thickness

    ! local variables

    integer, dimension(nx,ny) ::  &
         bmlt_float_mask         ! = 1 for cells where bmlt_float is computed and applied, else = 0
                                 ! start with bmlt_float_inversion mask; then remove grounded cells

    real(dp), dimension(nx,ny):: &
         thck_flotation,       & ! thickness at which ice becomes afloat (m)
         thck_cavity,          & ! thickness of ocean cavity beneath floating ice (m)
         thck_target,          & ! thickness target (m); = thck_obs unless thck_obs > thck_flotation
         dthck,                & ! thck - thck_target
         dbmlt_float             ! change in bmlt_float (m/s)

    integer :: i, j, ii, jj, iglobal, jglobal

    character(len=100) :: message

    real(dp) :: &
         term1,           & ! adjustment term for bmlt_float, proportional to thck - thck_target
         term2              ! adjustment term for bmlt_float, proportional to dthck/dt

    real(dp), parameter :: max_dbmlt_factor = 1.0d6  ! max multiplier for dbmlt_float, allowing for bmlt_weight ~ 0

    ! For floating cells, adjust the basal melt rate (or freezing rate, if bmlt < 0)
    !  so as to restore the upper surface to a target based on observations.

    ! Compute the flotation thickness
    where (topg - eus < 0.0d0)
       thck_flotation = -(rhoo/rhoi) *(topg - eus)
    elsewhere
       thck_flotation = 0.0d0
    endwhere

    ! Compute the ocean cavity thickness beneath floating ice (diagnostic only)
    thck_cavity = -(topg - eus) - (rhoi/rhoo)*thck_projected

    ! initialize
    thck_target(:,:) = 0.0d0
    dthck(:,:) = 0.0d0
    dbmlt_float(:,:) = 0.0d0
    bmlt_float_new(:,:) = 0.0d0

    ! Note: bmlt_float_inversion_mask is based on the initial geometry.
    !       Where this mask = 0, we never invert for bmlt_float.
    !       Where this mask = 1, we invert for bmlt_float in cells that satisfy the floating criterion.

    bmlt_float_mask = bmlt_float_inversion_mask

    ! Eliminate ice-free cells and fully grounded cells (based on bmlt_weight)
    where (ice_mask == 0 .or. bmlt_weight < tiny(0.0d0))
       bmlt_float_mask = 0
    endwhere

    ! For cells with bmlt_float_mask = 1, compute bmlt_float_inversion that will restore the thickness
    !  to the observed target.
    ! TODO: If not using the basal melting GLP, then restoring all the way to the grounded target will lead to oscillations,
    !       and we may need to use a buffer to prevent over-restoring.  For now, focus on runs with a basal melting GLP.


    ! loop over cells
    do j = 1, ny
       do i = 1, nx

          if (bmlt_float_mask(i,j) == 1) then  ! at least partly floating

             if (thck_obs(i,j) > thck_flotation(i,j)) then  ! grounded target

                thck_target(i,j) = thck_obs(i,j)

             else  ! floating target

                !TODO - Assuming we have assigned the buffer correctly at the beginning, can we just set thck_target = thck_obs?
                thck_target(i,j) = thck_obs(i,j)
                thck_target(i,j) = min(thck_target(i,j), thck_flotation(i,j) - thck_flotation_buffer)

             endif

             thck_target(i,j) = max(thck_target(i,j), 0.0d0)

             ! compute the difference between the projected thickness and the target thickness
             dthck(i,j) = thck_projected(i,j) - thck_target(i,j)

             ! Compute the rate of change of the melt rate.
             ! This rate of change is equal to the sum of two terms:
             !     dmb/dt = (H_new - H_target)/tau^2 + (2/tau) * dH/dt
             ! where mb is the basal melt rate, and tau = bmlt_timescale.
             ! This equation is similar to that of a damped harmonic oscillator:
             !     m * d2x/dt2 = -k*x - c*dx/dt
             ! A harmonic oscillator is critically damped when c = 2*sqrt(m*k).
             !  In this case the system is damped as strongly as possible without oscillating.
             ! Assuming unit mass (m = 1) and critical damping with k = 1/(tau^2), we obtain
             !   d2x/dt2 = -x/tau^2 - (2*dx/dt)/tau
             ! If we identify (H_new - H_target) with x; dH/dt with dx/dt; and d2x/dt2 with dmb/dt,
             !  we obtain the equation solved here.

             if (bmlt_timescale > dt) then
                term1 = dthck(i,j) / (bmlt_timescale)**2
                term2 = dthck_dt(i,j) * 2.0d0/bmlt_timescale
             else
                term1 = dthck(i,j) / (dt**2)
                term2 = 0.0d0  ! nonzero dH/dt term leads to oscillations when bmlt_timescale = dt
             endif

             dbmlt_float(i,j) = (term1 + term2) * dt

             ! Reduce the magnitude of dbmlt_float in cells with bmlt_weight < 1
             !  (partly grounded and/or shallow cavity)
             if (bmlt_weight(i,j) > 0.0d0) then  ! should have bmlt_weight > 0 where bmlt_float_mask = 1
                dbmlt_float(i,j) = dbmlt_float(i,j) / bmlt_weight(i,j)
             endif

             ! Increment bmlt_float
             bmlt_float_new(i,j) = bmlt_float_save(i,j) + dbmlt_float(i,j)

             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for bmlt_float_inversion: rank, i, j =', rtest, itest, jtest
                print*, 'bmlt_timescale (yr):', bmlt_timescale/scyr
                print*, 'bmlt_float_save, bmlt_weight:', &
                     bmlt_float_save(i,j)*scyr, bmlt_weight(i,j)
                print*, 'projected H, Hobs, dH:', thck_projected(i,j), thck_target(i,j), dthck(i,j)
                print*, 'dH/dt (m/yr):', dthck_dt(i,j)*scyr
                print*, 'dthck term, dthck/dt term, sum (m/yr):', &
                     term1*dt*scyr, term2*dt*scyr, dbmlt_float(i,j)*scyr
                print*, 'bmlt_float_new (m/yr):', bmlt_float_new(i,j)*scyr
             endif

          endif   ! bmlt_float_mask = 1

       enddo   ! i
    enddo   ! j

    call parallel_halo(bmlt_float_mask) ! diagnostic only
    call parallel_halo(thck_target)     ! diagnostic only
    call parallel_halo(bmlt_float_new)

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'bmlt_float mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') bmlt_float_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck_flotation (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') thck_flotation(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck_cavity (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') thck_cavity(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'H_target (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') thck_target(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'dthck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') dthck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'projected H, current bmlt_float (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') thck_projected(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'dH_dt_inversion (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') dthck_dt(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'term1 * dt (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') dthck(i,j)/(bmlt_timescale)**2 * dt * scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'term2 * dt (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') dthck_dt(i,j) * (2.0d0/bmlt_timescale) * dt * scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'dbmlt_float (m/yr), before weighting:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') bmlt_weight(i,j)*dbmlt_float(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'dbmlt_float (m/yr), after weighting:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') dbmlt_float(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine invert_bmlt_float

!***********************************************************************

  subroutine glissade_inversion_basal_friction(model)

    use glimmer_paramets, only: tim0, thk0
    use glimmer_physcon, only: scyr
    use glissade_grid_operators, only: glissade_stagger

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! --- Local variables ---

    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
         thck_obs                ! observed ice thickness, derived from usrf_obs and topg

    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) ::   &
         stag_thck,            & ! ice thickness on staggered grid
         stag_thck_obs,        & ! thck_obs on staggered grid
         stag_dthck_dt           ! dthck_dt on staggered grid

    integer :: i, j
    integer :: ewn, nsn
    integer :: itest, jtest, rtest

    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_unscaled

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

    if (model%options%which_ho_cp_inversion == HO_CP_INVERSION_COMPUTE) then

       ! Compute the new value of powerlaw_c_inversion.

       ! Given the surface elevation target, compute the thickness target.
       ! (This can change in time if the bed topography is dynamic.)
       call usrf_to_thck(model%geometry%usrf_obs,  &
                         model%geometry%topg,      &
                         model%climate%eus,        &
                         thck_obs)

       ! Interpolate thck_obs to the staggered grid
       ! Note: For this and the following fields, the interpolation will use values in all four neighbor cells,
       !       including ice-free cells.
       call glissade_stagger(ewn,         nsn,              &
                             thck_obs,    stag_thck_obs)

       ! Interpolate thck to the staggered grid
       call glissade_stagger(ewn,                  nsn,             &
                             model%geometry%thck,  stag_thck)

       ! Interpolate dthck_dt to the staggered grid
       call glissade_stagger(ewn,                      nsn,             &
                             model%geometry%dthck_dt,  stag_dthck_dt)

       ! Invert for powerlaw_c_inversion
       call invert_basal_friction(model%numerics%dt*tim0,                 &  ! s
                                  ewn,               nsn,                 &
                                  itest,    jtest,   rtest,               &
                                  model%inversion%babc_timescale,         &  ! s
                                  model%inversion%babc_thck_scale,        &  ! m
                                  model%inversion%powerlaw_c_max,         &
                                  model%inversion%powerlaw_c_min,         &
                                  model%geometry%f_ground,                &
                                  stag_thck*thk0,                         &  ! m
                                  stag_thck_obs*thk0,                     &  ! m
                                  stag_dthck_dt,                          &  ! m/s
                                  model%inversion%powerlaw_c_inversion)

    else

       ! do not adjust powerlaw_c_inversion; just print optional diagnostics
       if (verbose_inversion .and. this_rank == rtest) then
          print*, ' '
          print*, 'f_ground at vertices:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') model%geometry%f_ground(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'powerlaw_c_inversion:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') model%inversion%powerlaw_c_inversion(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    endif   ! which_ho_inversion

    ! Replace zeroes (if any) with small nonzero values to avoid divzeroes.
    ! Note: The current algorithm initializes Cp to a nonzero value everywhere and never sets Cp = 0;
    !       this check is just to be on the safe side.

    where (model%inversion%powerlaw_c_inversion == 0.0d0)
       model%inversion%powerlaw_c_inversion = model%inversion%powerlaw_c_min
    endwhere

  end subroutine glissade_inversion_basal_friction

!***********************************************************************

  subroutine invert_basal_friction(dt,                       &
                                   nx,            ny,        &
                                   itest, jtest,  rtest,     &
                                   babc_timescale,           &
                                   babc_thck_scale,          &
                                   powerlaw_c_max,           &
                                   powerlaw_c_min,           &
                                   f_ground,                 &
                                   stag_thck,                &
                                   stag_thck_obs,            &
                                   stag_dthck_dt,            &
                                   powerlaw_c_inversion)

    ! Compute a spatially varying basal friction field, powerlaw_c_inversion, defined at cell vertices.
    ! The method is similar to that of Pollard & DeConto (TC, 2012), and is applied to all grounded ice.
    ! Where stag_thck > stag_thck_obs, powerlaw_c is reduced to increase sliding.
    ! Where stag_thck < stag_thck_obs, powerlaw_c is increased to reduce sliding.
    ! Note: powerlaw_c is constrained to lie within a prescribed range.
    ! Note: For grounded ice with fixed topography, inversion based on thck is equivalent to inversion based on usrf.
    !       But for ice that is partly floating, it seems better to invert based on thck, because thck errors
    !        errors are greater in magnitude than errors in usrf, and we do not want to underweight the errors.
    !       With dynamic topography, we would either invert based on usrf, or else adjust thck_obs to match usrf_obs.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         babc_timescale,       & ! inversion timescale (s); must be > 0
         babc_thck_scale,      & ! thickness inversion scale (m); must be > 0
         powerlaw_c_max,       & ! upper bound for powerlaw_c, Pa (m/yr)^(-1/3)
         powerlaw_c_min          ! lower bound for powerlaw_c, Pa (m/yr)^(-1/3)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
         f_ground,             & ! grounded fraction at vertices, 0 to 1
         stag_thck,            & ! ice thickness at vertices (m)
         stag_thck_obs,        & ! observed ice thickness at vertices (m)
         stag_dthck_dt           ! rate of change of ice thickness at vertices (m/s)

    real(dp), dimension(nx-1,ny-1), intent(inout) ::  &
         powerlaw_c_inversion    ! powerlaw_c_inversion field to be adjusted

    ! local variables

    real(dp), dimension(nx-1,ny-1) ::  &
         stag_dthck,           & ! stag_thck - stag_thck_obs
         dpowerlaw_c             ! change in powerlaw_c

    real(dp) :: term1, term2
    integer :: i, j

    ! Initialize
    dpowerlaw_c(:,:) = 0.0d0

    ! Compute difference between current and target thickness
    stag_dthck(:,:) = stag_thck(:,:) - stag_thck_obs(:,:)

    ! optional diagnostics
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Old powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') powerlaw_c_inversion(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'stag_thck - stag_thck_obs:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') stag_dthck(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'stag_dthck_dt (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') stag_dthck_dt(i,j)*scyr
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'f_ground'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') f_ground(i,j)
          enddo
          print*, ' '
       enddo
    endif

    ! Loop over vertices where f_ground > 0
    ! Note: f_ground should be computed before transport, so that if a vertex is grounded
    !       before transport and fully floating afterward, powerlaw_c_inversion is computed here.

    do j = 1, ny-1
       do i = 1, nx-1

          if (f_ground(i,j) > 0.0d0) then  ! ice is at least partly grounded

             ! Compute the rate of change of powerlaw_c, based on stag_dthck and stag_dthck_dt.
             ! This rate of change is proportional to the sum of two terms:
             !     dCp/dt = -Cp * (1/tau) * (H - H_obs)/H0 + (2*tau/H0) * dH/dt
             ! where tau = babc_timescale and H0 = babc_thck_scale.
             ! This equation is similar to that of a damped harmonic oscillator:
             !     m * d2x/dt2 = -k*x - c*dx/dt
             ! where m is the mass, k is a spring constant, and c is a damping term.
             ! A harmonic oscillator is critically damped when c = 2*sqrt(m*k).
             !  In this case the system reaches equilibrium as quickly as possible without oscillating.
             ! Assuming unit mass (m = 1) and critical damping with k = 1/(tau^2), we obtain
             !   d2x/dt2 = -1/tau * (x/tau - 2*dx/dt)
             ! If we identify (H - H_obs)/(H0*tau) with x/tau; (2/H0)*dH/dt with 2*dx/dt; and (1/Cp)*dCp/dt with d2x/dt2,
             !  we obtain the equation solved here.

             term1 = -stag_dthck(i,j) / (babc_thck_scale * babc_timescale)
             term2 = -stag_dthck_dt(i,j) * 2.0d0 / babc_thck_scale

             dpowerlaw_c(i,j) = powerlaw_c_inversion(i,j) * (term1 + term2) * dt

             ! Limit to prevent a large relative change in one step
             if (abs(dpowerlaw_c(i,j)) > 0.05d0 * powerlaw_c_inversion(i,j)) then
                if (dpowerlaw_c(i,j) > 0.0d0) then
                   dpowerlaw_c(i,j) =  0.05d0 * powerlaw_c_inversion(i,j)
                else
                   dpowerlaw_c(i,j) = -0.05d0 * powerlaw_c_inversion(i,j)
                endif
             endif

             ! Update powerlaw_c
             powerlaw_c_inversion(i,j) = powerlaw_c_inversion(i,j) + dpowerlaw_c(i,j)

             ! Limit to a physically reasonable range
             powerlaw_c_inversion(i,j) = min(powerlaw_c_inversion(i,j), powerlaw_c_max)
             powerlaw_c_inversion(i,j) = max(powerlaw_c_inversion(i,j), powerlaw_c_min)

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                 print*, ' '
                print*, 'Invert for powerlaw_c: rank, i, j =', rtest, itest, jtest
                print*, 'thck, thck_obs, dthck, dthck_dt:', &
                     stag_thck(i,j), stag_thck_obs(i,j), stag_dthck(i,j), stag_dthck_dt(i,j)*scyr
                print*, 'dthck term, dthck_dt term, sum =', term1*dt, term2*dt, (term1 + term2)*dt
                print*, 'dpowerlaw_c, newpowerlaw_c =', dpowerlaw_c(i,j), powerlaw_c_inversion(i,j)
             endif

          else   ! f_ground = 0

             ! do nothing; keep the old value

          endif  ! f_ground > 0

       enddo  ! i
    enddo  ! j

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'New powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') powerlaw_c_inversion(i,j)
          enddo
          print*, ' '
       enddo
    endif   ! verbose_inversion

  end subroutine invert_basal_friction

!***********************************************************************

  subroutine glissade_inversion_bmlt_basin(dt,                          &
                                           nx,            ny,           &
                                           dx,            dy,           &
                                           itest, jtest,  rtest,        &
                                           nbasin,                      &
                                           basin_number,                &
                                           thck,                        &
                                           dthck_dt,                    &
                                           floating_thck_target,        &
                                           dbmlt_dtemp_scale,           &
                                           bmlt_basin_timescale,        &
                                           deltaT_basin)

    use glissade_bmlt_float, only: basin_sum

    ! For the case that bmlt_float is computed based on thermal_forcing,
    !  adjust deltaT_basin, which can be thought of as a bias corrrection
    !  or tuning parameter for the thermal forcing parameterization.
    ! In each  basin, we compute the area of marine-grounded ice and compare
    !  to the target area (usually based on observations).
    ! Where there is too much marine-grounded ice, we increase deltaT_basin,
    !  and where there is too little, we decrease deltaT_basin.
    ! Note: Other possible targets include the total floating area or grounded area.
    !       One reason not to use the total floating area is that the deltaT_basin
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
         floating_thck_target    ! target thickness for floating ice (m)

    real(dp), intent(in) :: &
         dbmlt_dtemp_scale,    & ! scale for rate of change of bmlt w/temperature, (m/s)/degC
         bmlt_basin_timescale    ! timescale for adjusting deltaT_basin (s)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         deltaT_basin            ! deltaT correction to thermal forcing in each basin (deg C)

    ! local variables

    real(dp), dimension(nx,ny) ::  &
         floating_target_rmask, &! real mask, = 1.0 where floating_thck_target > 0, else = 0.0
         cell_area               ! area of grid cells (m^2)

    real(dp), dimension(nbasin) :: &
         floating_area_target_basin,      &   ! floating ice area target in each basin (m^3)
         floating_volume_target_basin,    &   ! floating ice volume target in each basin (m^3)
         floating_thck_target_basin,      &   ! floating mean thickness target in each basin (m^3)
         floating_volume_basin,           &   ! current floating ice volume in each basin (m^3)
         floating_thck_basin,             &   ! current mean ice thickness in each basin (m)
         floating_dvolume_dt_basin,       &   ! rate of change of basin volume (m^3/s)
         floating_dthck_dt_basin,         &   ! rate of change of basin mean ice thickness (m/s)
         dTbasin_dt,                      &   ! rate of change of deltaT_basin (degC/s)
         basin_max, basin_min,            &   ! min and max of deltaT_basin in each basin
                                              ! (all cells in the basin should have the same value of deltaT_basin)
         deltaT_basin_nb                      ! same as deltaT_basin, but with dimension nbasin

    integer :: i, j
    integer :: nb      ! basin number
    real(dp) :: term1, term2

    ! Note: In some basins, the floating ice volume may be too small no matter how much we lower deltaT_basin,
    !        since the basal melt rate drops to zero and can go no lower.
    !       To prevent large negative values, the deltaT_basin correction is capped at a moderate negative value.
    !       A positive cap might not be needed but is included to be on the safe side.

    real(dp), parameter :: &
         deltaT_basin_maxval = 2.0d0,  &   ! max allowed magnitude of deltaT_basin (deg C)
         dTbasin_dt_maxval = 1.0d0/scyr     ! max allowed magnitude of d(deltaT_basin)/dt (deg/yr converted to deg/s)

    cell_area(:,:) = dx*dy

    ! Compute a mask for cells with a nonzero floating ice target

    where (floating_thck_target > 0.0d0)
       floating_target_rmask = 1.0d0
    elsewhere
       floating_target_rmask = 0.0d0
    endwhere

    ! For each basin, compute the area of the cells with floating_target_rmask = 1.

    call basin_sum(nx,         ny,                &
                   nbasin,     basin_number,      &
                   floating_target_rmask,         &
                   cell_area,                     &
                   floating_area_target_basin)

    ! For each basin, compute the target total ice volume in cells with floating_target_rmask = 1.
    ! Note: We could compute floating_volume_target_basin just once and write it to restart,
    !       but it is easy enough to recompute here.

    call basin_sum(nx,         ny,                &
                   nbasin,     basin_number,      &
                   floating_target_rmask,         &
                   floating_thck_target*dx*dy,    &
                   floating_volume_target_basin)

    ! For each basin, compute the current total ice volume in cells with floating_target_rmask = 1.

    call basin_sum(nx,         ny,                &
                   nbasin,     basin_number,      &
                   floating_target_rmask,         &
                   thck*dx*dy,                    &
                   floating_volume_basin)

    ! For each basin, compute the rate of change of the current volume in cells with floating_target_rmask = 1.

    call basin_sum(nx,         ny,                &
                   nbasin,     basin_number,      &
                   floating_target_rmask,         &
                   dthck_dt*dx*dy,                &
                   floating_dvolume_dt_basin)

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

    ! Compute the rate of change of deltaT_basin for each basin.
    ! Warm the basin where diff_ratio > 0 (too much ice) and cool where diff_ratio < 0 (too little ice).
    ! Note: deltaT_basin is a 2D field, but its value is uniform in each basin.

    do nb = 1, nbasin
       term1 = (1.0d0/dbmlt_dtemp_scale) * &
            (floating_thck_basin(nb) - floating_thck_target_basin(nb)) / (bmlt_basin_timescale**2)
       term2 = (1.0d0/dbmlt_dtemp_scale) * 2.0d0 * floating_dthck_dt_basin(nb) / bmlt_basin_timescale
       dTbasin_dt(nb) = term1 + term2
    enddo

    ! Limit the dTbasin)/dt to a prescribed range
    ! This prevents rapid changes in basins with small volume targets, where diff_ratio_basin can be large.
    where (dTbasin_dt > dTbasin_dt_maxval)
       dTbasin_dt = dTbasin_dt_maxval
    elsewhere (dTbasin_dt < -dTbasin_dt_maxval)
       dTbasin_dt = -dTbasin_dt_maxval
    endwhere

    ! Increment deltaT_basin
    do j = 1, ny
       do i = 1, nx
          nb = basin_number(i,j)
          if (nb >= 1 .and. nb <= nbasin) then
             deltaT_basin(i,j) = deltaT_basin(i,j) + dTbasin_dt(nb) * dt
          endif
       enddo
    enddo

    ! Limit deltaT_basin to a prescribed range
    where (deltaT_basin > deltaT_basin_maxval)
       deltaT_basin =  deltaT_basin_maxval
    elsewhere (deltaT_basin < -deltaT_basin_maxval)
       deltaT_basin = -deltaT_basin_maxval
    endwhere

    ! deltaT_basin diagnostics for each basin

    if (verbose_inversion) then

       !Note: Some variables are 2D fields rather than basin-only fields.
       !      The logic below extracts the basin values from the 2D fields.
       !      TODO: Write a subroutine to do this?

       basin_min(:) = 0.0d0
       basin_max(:) = 0.0d0

       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (nb >= 1 .and. nb <= nbasin) then
                basin_min(nb) = min(basin_min(nb), deltaT_basin(i,j))
                basin_max(nb) = max(basin_max(nb), deltaT_basin(i,j))
             endif
          enddo
       enddo

       do nb = 1, nbasin
          basin_min(nb) = parallel_reduce_min(basin_min(nb))
          basin_max(nb) = parallel_reduce_max(basin_max(nb))
       enddo

       deltaT_basin_nb = 0.0d0
       where (basin_min < 0.0d0)
          deltaT_basin_nb = basin_min
       elsewhere (basin_max > 0.0d0)
          deltaT_basin_nb = basin_max
       endwhere

       if (this_rank == rtest) then
          print*, 'bmlt_basin_timescale (yr) =', bmlt_basin_timescale/scyr
          print*, 'dbmlt_dtemp_scale (m/yr/degC) =', dbmlt_dtemp_scale
          print*, ' '
          print*, 'basin number, area target (km^2), volume target (km^3), mean thickness target (m):'
          do nb = 1, nbasin
             write(6,'(i6,3f12.3)') nb, floating_area_target_basin(nb)/1.d6, &
                  floating_volume_target_basin(nb)/1.d9, floating_thck_target_basin(nb)
          enddo
          print*, ' '
          print*, 'basin number, mean thickness (m), thickness diff (m), dthck_dt (m/yr):'
          do nb = 1, nbasin
             write(6,'(i6,3f12.3)') nb, floating_thck_basin(nb), &
                  (floating_thck_basin(nb) - floating_thck_target_basin(nb)), &
                  floating_dthck_dt_basin(nb)*scyr
          enddo
          print*, ' '
          print*, 'basin number, term1*dt, term2*dt, dTbasin, new deltaT_basin:'
          do nb = 1, nbasin
             write(6,'(i6,4f12.6)') nb, &
                  dt/dbmlt_dtemp_scale * (floating_thck_basin(nb) - floating_thck_target_basin(nb)) / (bmlt_basin_timescale**2), &
                  dt/dbmlt_dtemp_scale * 2.0d0 * floating_dthck_dt_basin(nb) / bmlt_basin_timescale, &
                  dt*dTbasin_dt(nb), deltaT_basin_nb(nb)
          enddo
       endif

    endif   ! verbose_inversion

  end subroutine glissade_inversion_bmlt_basin

!***********************************************************************

  !TODO - Move the two following subroutines to a utility module?

  subroutine usrf_to_thck(usrf, topg, eus, thck)

    ! Given the bed topography and upper ice surface elevation, compute the ice thickness.
    ! The ice is assumed to satisfy a flotation condition.
    ! That is, if topg - eus < 0 (marine-based ice), and if the upper surface is too close
    !  to sea level to ground the ice, then the ice thickness is chosen to satisfy
    !  rhoi*H = -rhoo*(topg-eus).
    ! Note: usrf, topg, eus and thck must all have the same units (often but not necessarily meters).

    use glimmer_physcon, only : rhoo, rhoi

    real(dp), dimension(:,:), intent(in) :: &
         usrf,           & ! ice upper surface elevation
         topg              ! elevation of bedrock topography

    real(dp), intent(in) :: &
         eus               ! eustatic sea level

    real(dp), dimension(:,:), intent(out) :: &
         thck              ! ice thickness

    ! initialize
    thck(:,:) = 0.0d0

    where (usrf > (topg - eus))   ! ice is present, thck > 0
       where (topg - eus < 0.0d0)   ! marine-based ice
          where ((topg - eus) * (1.0d0 - rhoo/rhoi) > usrf)  ! ice is floating
             thck = usrf / (1.0d0 - rhoi/rhoo)
          elsewhere   ! ice is grounded
             thck = usrf - (topg - eus)
          endwhere
       elsewhere   ! land-based ice
          thck = usrf - (topg - eus)
       endwhere
    endwhere

  end subroutine usrf_to_thck

!***********************************************************************

  subroutine thck_to_usrf(thck, topg, eus, usrf)

    ! Given the bed topography and ice thickness, compute the upper surface elevation.
    ! The ice is assumed to satisfy a flotation condition.
    ! That is, if topg - eus < 0 (marine-based ice), and if the ice is too thin to be grounded,
    !  then the upper surface is chosen to satisfy rhoi*H = rhoo*(H - usrf),
    !  or equivalently usrf = (1 - rhoi/rhoo)*H.
    ! Note: usrf, topg, eus and thck must all have the same units (often but not necessarily meters).

    use glimmer_physcon, only : rhoo, rhoi

    real(dp), dimension(:,:), intent(in) :: &
         thck,           & ! ice thickness
         topg              ! elevation of bedrock topography

    real(dp), intent(in) :: &
         eus               ! eustatic sea level

    real(dp), dimension(:,:), intent(out) :: &
         usrf              ! ice upper surface elevation

    where ((topg - eus) < -(rhoi/rhoo)*thck)
       usrf = (1.0d0 - rhoi/rhoo)*thck   ! ice is floating
    elsewhere   ! ice is grounded
       usrf = (topg - eus) + thck
    endwhere

  end subroutine thck_to_usrf

!=======================================================================

end module glissade_inversion

!=======================================================================
