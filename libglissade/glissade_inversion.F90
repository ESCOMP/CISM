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
  use glimmer_paramets, only: thk0
  use glimmer_log
  use glide_types
  use glide_thck, only: glide_calclsrf
  use parallel

  implicit none

  private
  public :: verbose_inversion, glissade_init_inversion, &
            glissade_inversion_bmlt_float, glissade_inversion_basal_traction

  ! All subroutines in this module are public

  !-----------------------------------------------------------------------------
  ! Subroutines to invert for basal fields (including basal traction beneath
  ! grounded ice and basal melting beneath floating ice) by relaxing toward
  ! a target ice thickness field.
  !-----------------------------------------------------------------------------

    logical, parameter :: verbose_inversion = .false.
!!    logical, parameter :: verbose_inversion = .true.

    real(dp), parameter :: eps08 = 1.0d-08  ! small number

!***********************************************************************

contains

!***********************************************************************

  subroutine glissade_init_inversion(model)

    ! Initialize inversion for fields of basal traction and basal melting
    ! Should be called after usrf and thck have been input and (possibly) modified by initial calving

    use glissade_masks, only: glissade_get_masks
    use glissade_calving, only: glissade_ocean_connection_mask  !TODO - Move to mask module
    use glissade_grounding_line, only: glissade_grounded_fraction
    use parallel

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    integer :: i, j

    integer :: itest, jtest, rtest  ! local diagnostic point

    real(dp) :: var_maxval          ! max value of a given real variable; = 0.0 if not yet read in
    integer :: var_maxval_int       ! max value of a given integer variable; = 0 if not yet read in

    character(len=100) :: message

    integer, dimension(model%general%ewn, model%general%nsn) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         ocean_mask,           & ! = 1 where topg is below sea level and ice is absent
         land_mask,            & ! = 1 where topg is at or above sea level
         lake_mask               ! = 1 for inland lakes

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck_flotation          ! flotation thickness

    real(dp) :: h_obs, h_flotation, h_buff   ! thck_obs, thck_flotation, and thck_flotation_buffer scaled to m
    real(dp) :: dh                           ! h_obs - f_flotation
    real(dp) :: dh_decimal                   ! decimal part remaining after subtracting the truncation of dh

    integer :: ewn, nsn


    real(dp) :: powerlaw_c_init_high         ! Cp for high topography
    real(dp) :: powerlaw_c_init_low          ! Cp for low topography

    !TODO - Make some of these parameters user-configurable?

    !WHL - The following values give good results for a 4-km simulation in Dec. 2018.
    !      There is a narrow trough flowing into the Shackleton shelf that overthickens and
    !      can cause a crash around year 400 (with dt = 0.20 yr) unless the bed is
    !      initialized to have faster sliding in the trough and slower sliding outside the trough.
    !TODO - Revisit these values; see if we can find a single initial value of powerlaw_c
    !       that works well for all ice-covered cells.

    real(dp) :: powerlaw_c_init_high_frac = 0.80d0   ! fraction of Cp_max for cells at high topography
    real(dp) :: powerlaw_c_init_low_frac  = 0.20d0   ! fraction of Cp_max for cells at low topography
    real(dp) :: powerlaw_c_init_topg_high = 500.d0   ! high topography value (m)
    real(dp) :: powerlaw_c_init_topg_low = -500.d0   ! low topography value (m)

    real(dp), parameter :: &
         f_ground_threshold = tiny(0.0d0)  ! used to define ocean connection mask; path to ocean must pass
                                           ! through cells with f_ground less than threshold
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

    ! Get masks
    call glissade_get_masks(ewn,                 nsn,                   &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   model%numerics%thklim, &
                            ice_mask,                                   &
                            land_mask = land_mask)

    !TODO - Check that the 'is_restart' logic is correct.

    if (model%options%which_ho_inversion == HO_INVERSION_COMPUTE) then

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
       ! We are not inverting for thck_obs, but it can be a useful diagnostic.

       thck_flotation(:,:) = 0.0d0
       model%geometry%thck_obs(:,:) = 0.0d0

       !TODO - Write a subroutine to compute thck, given usrf, topg, and eus?
       where (model%geometry%thck > 0.0d0)  ! ice is present
          where (model%geometry%topg - model%climate%eus < 0.0d0)  ! marine-based ice
             thck_flotation = -(rhoo/rhoi) * (model%geometry%topg - model%climate%eus)
             where (model%geometry%topg + thck_flotation < model%geometry%usrf_obs)  ! grounded
                model%geometry%thck_obs = model%geometry%usrf_obs - (model%geometry%topg - model%climate%eus)
             elsewhere  ! floating
                model%geometry%thck_obs = model%geometry%usrf_obs / (1.0d0 - rhoi/rhoo)
             endwhere
          elsewhere   ! land-based ice
             thck_flotation = 0.0d0
             model%geometry%thck_obs = model%geometry%usrf_obs - (model%geometry%topg - model%climate%eus)
          endwhere
       endwhere

       if (model%options%is_restart == RESTART_FALSE) then

          ! At the start of the run, adjust thck_obs so that the observational target is not too close to thck_flotation.
          ! The reason for this is that if we restore H to values very close to thck_flotation,
          !  it can be easier for cells to flip between grounded and floating in the forward run.
          ! Note: The f_ground algorithm can have problems if we have 4 adjacent cells with a checkerboard pattern
          !        of equal and opposite values of dh = thck_obs - thck_flotation.
          !       To prevent this, do not make the new dh exactly equal to thck_flotation_buffer,
          !        but rather increment it up or down in integer steps until we exceed the buffer.
          !       E.g., if the buffer is 1 m, then a cavity of 0.357 m is increased to 1.357 m.

          do j = 1, model%general%nsn
             do i = 1, model%general%ewn
                if (model%geometry%thck_obs(i,j) > 0.0d0 .and.  &
                     model%geometry%topg(i,j) - model%climate%eus < 0.0d0) then
                   ! convert to m (can skip the conversion when code scaling is removed)
                   h_obs = model%geometry%thck_obs(i,j) * thk0
                   h_flotation = thck_flotation(i,j) * thk0
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
                      model%geometry%thck_obs(i,j) = h_obs / thk0
                   endif
                endif
                model%geometry%thck_obs(i,j) = max(model%geometry%thck_obs(i,j), 0.0d0)
             enddo
          enddo

          ! Where thck_obs < inversion_thck_threshold, set it to zero.
          ! One reason to do this is to avoid restoring ice to small values at the calving front.

          model%inversion%thck_threshold = max(model%inversion%thck_threshold, model%numerics%thklim)
          where (model%geometry%thck_obs <= model%inversion%thck_threshold)
             model%geometry%thck_obs = 0.0d0
          endwhere

          ! Set thck to be consistent with thck_obs
          model%geometry%thck = model%geometry%thck_obs

          ! Reset usrf_obs to be consistent with thck_obs.
          ! (usrf itself will be reset later in glissade_initialise)
          !TODO - Use glide_calclsrf/usrf instead
          where (model%geometry%topg - model%climate%eus < (-rhoi/rhoo) * model%geometry%thck_obs)
             model%geometry%usrf_obs = (1.0d0 - rhoi/rhoo) * model%geometry%thck_obs  ! floating
          elsewhere
             model%geometry%usrf_obs = model%geometry%topg +  model%geometry%thck_obs ! grounded
          endwhere

       endif   ! not a restart

       call parallel_halo(model%geometry%thck_obs)
       call parallel_halo(model%geometry%usrf_obs)

       if (model%options%is_restart == RESTART_FALSE) then

          ! At the start of the run, compute a mask that determines where bmlt_float_inversion can be applied.
          ! This mask includes only cells with a direct path to the ocean through cells without grounded ice or ice-free land.
          ! During the simulation, additional cells may ground, in which case bmlt_float_inversion will not
          !  be applied in those cells either.

          ! Find inland lakes: cells that are floating but have no connection through other floating cells to the ocean.

          ! Recompute masks
          call glissade_get_masks(ewn,                 nsn,                   &
                                  model%geometry%thck, model%geometry%topg,   &
                                  model%climate%eus,   model%numerics%thklim, &
                                  ice_mask,                                   &
                                  floating_mask = floating_mask,              &
                                  ocean_mask = ocean_mask,                    &
                                  land_mask = land_mask)

          if (model%options%which_ho_ground == HO_GROUND_GLP_DELUXE) then

             ! use f_ground_cell to find ocean connections, instead of floating_mask as computed above

             call glissade_grounded_fraction(ewn,      nsn,                 &
                                             itest, jtest, rtest,           &  ! diagnostic only
                                             model%geometry%thck*thk0,      &
                                             model%geometry%topg*thk0,      &
                                             model%climate%eus*thk0,        &
                                             ice_mask,                      &
                                             floating_mask,                 &
                                             land_mask,                     &
                                             model%options%which_ho_ground, &
                                             model%options%which_ho_flotation_function, &
                                             model%options%which_ho_fground_no_glp,     &
                                             model%geometry%f_flotation,    &
                                             model%geometry%f_ground,       &
                                             model%geometry%f_ground_cell)

             ! Use f_ground_cell to adjust floating_mask as input to glissade_ocean_connection_mask.
             ! Cells with f_ground_cell greater than a prescribed threshold do not serve as paths to the ocean.

             where (model%geometry%f_ground_cell > f_ground_threshold)
                floating_mask = 0
             endwhere

          endif  ! which_ho_ground

          call glissade_ocean_connection_mask(ewn,          nsn,            &
                                              itest, jtest, rtest,          &
                                              ice_mask,     floating_mask,  &
                                              ocean_mask,   land_mask,      &
                                              model%inversion%bmlt_float_inversion_mask)

       endif  ! not a restart

       call parallel_halo(model%inversion%bmlt_float_inversion_mask)

       !TODO - Remove this code if not inverting for topography
       ! Check whether topg_obs has been read in already.
       ! If not, then set topg_obs to the initial topography.
!       var_maxval = maxval(model%geometry%topg_obs)
!       var_maxval = parallel_reduce_max(var_maxval)
!       if (var_maxval > 0.0d0) then
          ! do nothing
!       else
          ! initialize to the input topography
!          model%geometry%topg_obs(:,:) = model%geometry%topg(:,:)
!       endif
!       call parallel_halo(model%geometry%topg_obs)

       ! Initialize powerlaw_c_save, if not already read in.
       ! This is the value saved after each time step, which optionally can be included
       !  in a weighted average during the following time step.
       var_maxval = maxval(model%inversion%powerlaw_c_save)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; powerlaw_c_save has been read in already (e.g., when restarting)
       else
!!          ! setting to a large value so that basal flow starts slow and gradually speeds up as needed
!!          model%inversion%powerlaw_c_save(:,:) = model%inversion%powerlaw_c_max

          powerlaw_c_init_high = model%inversion%powerlaw_c_max * powerlaw_c_init_high_frac
          powerlaw_c_init_low  = model%inversion%powerlaw_c_max * powerlaw_c_init_low_frac

          where (model%geometry%thck_obs > 0.0d0)
             ! uncomment to set to a large value everywhere
!!             model%inversion%powerlaw_c_save = model%inversion%powerlaw_c_max * 0.5d0

             !TODO - Try setting to a lower value and see if flow is still stable?
!!             model%inversion%powerlaw_c_save = model%inversion%powerlaw_c_land

             ! make it a function of topography
             where (model%geometry%topg * thk0 > powerlaw_c_init_topg_high)
                model%inversion%powerlaw_c_save = powerlaw_c_init_high
             elsewhere (model%geometry%topg * thk0 < powerlaw_c_init_topg_low)
                model%inversion%powerlaw_c_save = powerlaw_c_init_low
             elsewhere (model%geometry%topg * thk0 >= powerlaw_c_init_topg_low .and.  &
                        model%geometry%topg * thk0 <= powerlaw_c_init_topg_high)
                model%inversion%powerlaw_c_save = powerlaw_c_init_low  &
                     + (powerlaw_c_init_high - powerlaw_c_init_low)  &
                     * (model%geometry%topg*thk0 - powerlaw_c_init_topg_low) / &
                       (powerlaw_c_init_topg_high - powerlaw_c_init_topg_low)
             endwhere

          elsewhere (land_mask == 1)
             model%inversion%powerlaw_c_save = model%inversion%powerlaw_c_land
          elsewhere   ! ice-free ocean
             model%inversion%powerlaw_c_save = model%inversion%powerlaw_c_marine
          endwhere

       endif  ! var_maxval > 0

       call parallel_halo(model%inversion%powerlaw_c_save)

       ! Note: There is no initialization of bmlt_float_save.
       ! If restarting, it should have been read in already.
       ! If not restarting, it will have been set to zero, which is an appropriate initial value.

       call parallel_halo(model%inversion%bmlt_float_save)

       !WHL - temporary setting - Copy bmlt_float_save to bmlt_float_inversion at restart.
       !                          Then no longer use bmlt_float_save.
       model%inversion%bmlt_float_inversion = model%inversion%bmlt_float_save

    endif  ! which_ho_inversion

    !WHL - debug
    if (this_rank == rtest) then
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
             write(6,'(f10.3)',advance='no') model%geometry%thck_obs(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_init_inversion

!***********************************************************************

  subroutine glissade_inversion_bmlt_float(model,               &
                                           nudging_factor,      &
                                           thck_new_unscaled,   &
                                           ice_mask,            &
                                           floating_mask,       &
                                           land_mask)

    use parallel

    use glimmer_paramets, only: tim0, thk0
    use glimmer_physcon, only: scyr

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    real(dp), intent(in) :: &
         nudging_factor          ! factor in range [0,1]
                                 ! = 1 to nudge fully toward the new value, = 0 to keep the current value

    real(dp), dimension(model%general%ewn, model%general%nsn), intent(in) ::   &
         thck_new_unscaled       ! ice thickness expected after mass balance, without inversion (m)

    !Note: These masks are not part of the model derived type, and they are computed before transport
    !      based on the old ice thickness, so they cannot be computed here.
    !TODO - Make these masks part of the model derived type, so they do not need to be passed in?

    integer, dimension(model%general%ewn, model%general%nsn), intent(in) ::   &
         ice_mask,             & ! = 1 if thck > 0, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask               ! = 1 if topg is at or above sea level, else = 0

    ! --- Local variables ---

    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
         thck_unscaled,        & ! ice thickness (m)
         topg_unscaled,        & ! bedrock topography (m)
         lsrf_new_unscaled,    & ! expected new lower surface elevation (m)
         usrf_new_unscaled,    & ! expected new upper surface elevation (m)
         bmlt_float_new,       & ! newly computed value of bmlt_float (m/s)
         correction              ! correction term in range [-1,1]

    !WHL - Make this a config parameter?
    real(dp) :: correction_time_scale = 2.0d0    ! time scale (yr) for correction term proportional to dH/dt

    integer :: i, j
    integer :: ewn, nsn
    integer :: itest, jtest, rtest

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

    topg_unscaled = model%geometry%topg * thk0

    ! Calculate the expected new lower and upper ice surface
    ! Note: usrf_new_unscaled is used in inversion calculations, but model%geometry%usrf is not updated
    !       until after the mass balance calculation. 
    call glide_calclsrf(thck_new_unscaled, topg_unscaled, model%climate%eus*thk0, lsrf_new_unscaled)
    usrf_new_unscaled = max(0.d0, thck_new_unscaled + lsrf_new_unscaled)

    !TODO - Pass this in, so not computed in two different places?
    if (model%options%which_ho_inversion == HO_INVERSION_COMPUTE) then

       ! Invert for bmlt_float_inversion, adjusting the melt rate to relax toward the observed thickness.
       ! Note: basal_melt%bmlt_float_inversion is passed out with units of m/s

       ! Note: Other kinds of sub-shelf basal melting are handled in subroutine glissade_bmlt_float_solve.
       !       Inversion is done here, after transport, when there is an updated ice thickness.
       !       Then bmlt_float_inversion is added to the previously computed bmlt.
       ! Note: Typically, whichbmlt_float = 0 when doing a model spin-up with inversion.
       !       However, we might want to add an anomaly to fields already computed by inversion.
       ! Note: If the basal melt GLP is turned on, it sets bmlt_float = 0 in partly floating cells.
       !       However, it does not limit bmlt_float_inversion, which is applied to all floating cells,
       !       including partly floating cells (in order to match observed thicknesses at the grounding line).

       ! Compute the new value of bmlt_float_inversion.
       ! Optionally, the new value is then combined with the saved value in a weighted average.

       call invert_bmlt_float(model%numerics%dt * tim0,               &    ! s
                              ewn,               nsn,                 &
                              itest,   jtest,    rtest,               &
                              model%inversion,                        &
                              thck_new_unscaled,                      &    ! m
                              model%geometry%usrf_obs*thk0,           &    ! m
                              topg_unscaled,                          &    ! m
                              model%climate%eus*thk0,                 &    ! m
                              model%options%which_ho_ground,          &
                              model%options%which_ho_ground_bmlt,     &
                              model%geometry%f_ground_cell,           &
                              ice_mask,                               &
                              floating_mask,                          &
                              bmlt_float_new)

       ! Depending on the value of which_ho_ground_bmlt, reduce or zero out bmlt_float_new where the ice is grounded.
       ! TODO - Move this code inside invert_bmlt_float?
       !        Then it comes out already weighted by (1 - f_ground_cell).

       if (model%options%which_ho_ground == HO_GROUND_GLP_DELUXE) then

          ! reduce bmlt_float_inversion based on f_ground_cell

          if (model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_FLOATING_FRAC) then

             where (ice_mask == 1 .and. model%geometry%f_ground_cell > 0.0d0)
                bmlt_float_new = bmlt_float_new * (1.0d0 - model%geometry%f_ground_cell)
             endwhere

          elseif (model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_ZERO_GROUNDED) then

             where (ice_mask == 1 .and. model%geometry%f_ground_cell >= eps08)
               bmlt_float_new = 0.0d0
             endwhere

          endif   ! which_ho_ground_bmlt

       else

          ! reduce bmlt_float_inversion based on floating_mask

          where (ice_mask == 1 .and. floating_mask == 0)
             bmlt_float_new = 0.0d0
          endwhere

       endif  ! which_ho_ground


       !WHL - debug
       if (verbose_inversion .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'Inverting for bmlt_float: rank, i, j =', rtest, i, j
          print*, 'usrf (m), usrf_obs (m), new bmlt_float (m/yr):', usrf_new_unscaled(i,j), &
               model%geometry%usrf_obs(i,j)*thk0, bmlt_float_new(i,j)*scyr
          print*, 'bmlt_float old value, new value, nudging factor:', &
               model%inversion%bmlt_float_inversion(i,j)*scyr, bmlt_float_new(i,j)*scyr, nudging_factor
          print*, ' '
          print*, 'new bmlt_float inversion, limited by f_ground_cell (m/yr):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') bmlt_float_new(i,j)*scyr
             enddo
             write(6,*) ' '
          enddo
       endif

       model%inversion%bmlt_float_inversion(:,:) = nudging_factor  * bmlt_float_new(:,:) &
                                        + (1.0d0 - nudging_factor) * model%inversion%bmlt_float_inversion(:,:)


       ! Adjust bmlt_float_inversion if dH/dt is substantial.
       ! This term damps cycling due to perturbations in m lagging perturbations in H.
       !TODO - Combine with nudging_factor above?

       correction(:,:) = 0.0d0

       where (bmlt_float_new /= 0.0d0)
          correction = model%geometry%dthck_dt * (model%numerics%dt * tim0/scyr) / correction_time_scale
          model%inversion%bmlt_float_inversion = model%inversion%bmlt_float_inversion + correction
       endwhere

       if (verbose_inversion .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, 'bmlt_float nudged value =', model%inversion%bmlt_float_inversion(i,j)*scyr
       endif

       !TODO - Remove bmlt_float_save when no longer needed; but do a copy for now.
       model%inversion%bmlt_float_save = model%inversion%bmlt_float_inversion

       if (verbose_inversion .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'f_ground_cell:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.6)',advance='no') model%geometry%f_ground_cell(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'dH/dt (m/yr):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') model%geometry%dthck_dt(i,j)*scyr
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'correction term:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') correction(i,j)*scyr
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'nudged bmlt_float inversion (m/yr):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') model%inversion%bmlt_float_inversion(i,j)*scyr
!                write(6,'(e10.3)',advance='no') model%inversion%bmlt_float_inversion(i,j)*scyr ! to view small values
             enddo
             write(6,*) ' '
          enddo
       endif   ! verbose_inversion

    endif   ! which_ho_inversion

  end subroutine glissade_inversion_bmlt_float

!***********************************************************************

  subroutine glissade_inversion_basal_traction(model,          &
                                               nudging_factor, &
                                               ice_mask,       &
                                               floating_mask,  &
                                               land_mask)

    use parallel

    use glimmer_paramets, only: tim0, thk0
    use glimmer_physcon, only: scyr
    use glissade_grid_operators, only: glissade_stagger, glissade_stagger_real_mask

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    real(dp), intent(in) :: &
         nudging_factor          ! factor in range [0,1]
                                 ! = 1 to nudge fully toward the new value, = 0 to keep the current value

    !TODO - Compute these locally?
    integer, dimension(model%general%ewn, model%general%nsn), intent(in) ::   &
       ice_mask,             & ! = 1 if thck > 0, else = 0
       floating_mask,        & ! = 1 where ice is present and floating, else = 0
       land_mask               ! = 1 if topg is at or above sea level, else = 0

    ! --- Local variables ---

    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
       powerlaw_c_new          ! newly computed value of powerlaw_c, Pa (m/yr)^(-1/3)

    real(dp) :: alpha          ! shorthand for inversion%babc_time_smoothing, in range [0,1]

    integer :: i, j
    integer :: ewn, nsn
    integer :: itest, jtest, rtest

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

    !TODO - Check exact restart.  
    if (model%options%which_ho_inversion == HO_INVERSION_COMPUTE) then

       ! Update the upper ice surface
       call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
       model%geometry%usrf = max(00.d0, model%geometry%lsrf + model%geometry%thck)

       ! Optionally, compute an exponential moving average of usrf.
       ! The smaller the factor, the more rapidly earlier values are discounted.
       ! The default is alpha = 0, in which case earlier values are ignored.
       !TODO - Is the moving average of usrf still needed?

       alpha = model%inversion%babc_time_smoothing
       alpha = min(alpha, 1.0d0 - 1.0d0/real(model%numerics%tstep_count,dp))  ! decrease smoother for first few time steps
       alpha = min(1.0d0, max(alpha,0.0d0))  ! limit to [0,1]
       if (alpha < 1.0d0) then
          ! take moving averages of usrf with contributions from previous values
          model%inversion%usrf_inversion(:,:) = (1.d0 - alpha) * model%geometry%usrf(:,:)*thk0  &
                                                      + alpha  * model%inversion%usrf_inversion(:,:)
       else
          ! simply copy the latest values
          model%inversion%usrf_inversion(:,:) = model%geometry%usrf(:,:)
       endif   ! alpha < 1

       !WHL - debug
       if (verbose_inversion .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'Computed moving averages: rank, i, j, alpha =', rtest, i, j, alpha
          print*, ' '
          print*, 'current usrf (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') model%geometry%usrf(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'moving average usrf:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') model%inversion%usrf_inversion(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'dH/dt from previous time step (m/yr):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') model%geometry%dthck_dt(i,j) * scyr
             enddo
             write(6,*) ' '
          enddo
       endif

       ! Determine the basal friction field, powerlaw_c_inversion, if desired.
       ! Notes: (1) For inversion purposes, ice_mask = 1 where thck > 0.0 (not where thck > thklim).
       !        (2) usrf_unscaled is the expected new value after applying the mass balance.
       !        (3) These masks are computed before horizontal transport. So for instance, if a cell
       !            is grounded before transport and floating afterward, it is treated as grounded.

       call invert_basal_traction(model%numerics%dt*tim0,                 &  ! s
                                  ewn,               nsn,                 &
                                  itest,    jtest,   rtest,               &
                                  model%inversion,                        &
                                  ice_mask,                               &
                                  land_mask,                              &
                                  floating_mask,                          & 
                                  model%options%which_ho_ground,          &
                                  model%geometry%f_ground_cell,           &
                                  model%inversion%usrf_inversion,         &  ! m
                                  model%geometry%usrf_obs*thk0,           &  ! m
                                  model%geometry%dthck_dt,                &  ! m/s
                                  powerlaw_c_new)

       if (verbose_inversion .and. this_rank == rtest) then
          print*, ' '
          print*, 'tstep_count, time, basal friction nudging_factor =', &
            model%numerics%tstep_count, model%numerics%time, nudging_factor
       endif

       !Note: The saved value is adjusted only where the new value is nonzero
       !      (i.e., where the ice is at least partly grounded).
       !      This means that zero values when a cell is ice-free or floating are not included in the saved average.
       !TODO - Replace powerlaw_c_save with powerlaw_c_inversion.

       where (powerlaw_c_new > 0.0d0)
          model%inversion%powerlaw_c_save(:,:) = nudging_factor  * powerlaw_c_new(:,:) &
                                      + (1.0d0 - nudging_factor) * model%inversion%powerlaw_c_save(:,:)
       endwhere

       if (verbose_inversion .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, 'powerlaw_c saved value, new value:', model%inversion%powerlaw_c_save(i,j), powerlaw_c_new(i,j)
          print*, 'powerlaw_c nudged value =', model%inversion%powerlaw_c_save(i,j)
       endif

       ! Set powerlaw_c_inversion = powerlaw_c_save where the ice is at least partly grounded.
       ! Elsewhere, set powerlaw_c_inversion = 0.

       where (land_mask == 1 .or. model%geometry%f_ground_cell > 0.0d0) 
          model%inversion%powerlaw_c_inversion = model%inversion%powerlaw_c_save
       elsewhere
          model%inversion%powerlaw_c_inversion = 0.0d0
       endwhere

       ! Interpolate powerlaw_c_inversion to the staggered grid.

       if (model%options%which_ho_ground == HO_GROUND_GLP_DELUXE) then

          ! For the staggered averaging, give each cell a weight of f_ground_cell

          call glissade_stagger_real_mask(ewn,            nsn,                         &
                                          model%inversion%powerlaw_c_inversion,        &
                                          model%inversion%stag_powerlaw_c_inversion,   &
                                          model%geometry%f_ground_cell)
       else

          ! For the staggered averaging, give each cell a weight of 1 if ice-covered, else 0
          ! stagger_margin_in = 1: Interpolate using only the values in ice-covered cells
          !TODO - Weigh based on floating_mask?

          call glissade_stagger(ewn,             nsn,                        &
                                model%inversion%powerlaw_c_inversion,       &
                                model%inversion%stag_powerlaw_c_inversion,  &
                                ice_mask,                                   &
                                stagger_margin_in = 1)

       endif   ! which_ho_ground

       ! Replace zeroes with default values to avoid divzeroes
       where (model%inversion%stag_powerlaw_c_inversion == 0.0d0)
          model%inversion%stag_powerlaw_c_inversion = model%inversion%powerlaw_c_min
       endwhere

       if (verbose_inversion .and. this_rank == rtest) then
          print*, ' '
          print*, 'stag_powerlaw_c:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') model%inversion%stag_powerlaw_c_inversion(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    endif   ! which_ho_inversion

  end subroutine glissade_inversion_basal_traction

!***********************************************************************

  !TODO - Remove this subroutine?
  subroutine invert_basal_topography(dt,                       &
                                     nx,            ny,        &
                                     itest, jtest,  rtest,     &
                                     ice_mask,                 &
                                     grounding_line_mask,      &
                                     usrf,                     &
                                     usrf_obs,                 &
                                     topg,                     &
                                     topg_obs,                 &
                                     eus)

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    integer, dimension(nx,ny), intent(in) :: &
         ice_mask,             & ! = 1 where ice is present (thk > 0), else = 0
         grounding_line_mask     ! = 1 if a cell is adjacent to the grounding line, else = 0

    ! Note: usrf should be the expected new value of usrf after applying the mass balance
    !       (although the mass balance may not yet have been applied)
    real(dp), dimension(nx,ny), intent(in) ::  &
         usrf,                 & ! upper surface elvation (m)
         usrf_obs,             & ! observed upper surface elvation (m)
         topg_obs                ! observed basal topography (m)

    real(dp), intent(in) :: &
         eus                     ! eustatic sea level (m)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         topg                    ! basal topography (m)

    ! local variables

    !TODO - Make these config parameters?
    real(dp), parameter :: &
!!         topg_inversion_timescale = 1000.d0*scyr,  & ! timescale for topg inversion, yr converted to s
         topg_inversion_timescale = 100.d0*scyr,  & ! timescale for topg inversion, yr converted to s
         topg_maxcorr = 100.d0                       ! max allowed correction in topg, compared to obs (m)
 
    real(dp) :: dtopg_dt

    integer :: i, j

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Before topg adjustment, topg:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') topg(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    !TODO - Apply to grounded cells only?
    ! Adjust basal topography in cells adjacent to the grounding line.
    ! Raise the topg where usrf < usrf_obs, and lower topg where usrf > usrf_obs
    ! Note: The grounding_line mask is computed before horizontal transport.
    !       It includes grounded cells adjacent to at least one floating or ice-free ocean cell,
    !        and floating cells adjacent to at least one grounded cell. 
    do j = 1, ny
       do i = 1, nx
          if (ice_mask(i,j) == 1 .and. grounding_line_mask(i,j) == 1) then
             dtopg_dt = -(usrf(i,j) - usrf_obs(i,j)) / topg_inversion_timescale
             topg(i,j) = topg(i,j) + dtopg_dt*dt
             topg(i,j) = min(topg(i,j), topg_obs(i,j) + topg_maxcorr)
             topg(i,j) = max(topg(i,j), topg_obs(i,j) - topg_maxcorr)
          endif
       enddo
    enddo

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'After topg adjustment, topg:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') topg(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine invert_basal_topography

!***********************************************************************

  subroutine invert_basal_traction(dt,                       &
                                   nx,            ny,        &
                                   itest, jtest,  rtest,     &
                                   inversion,                &
                                   ice_mask,                 &
                                   land_mask,                &
                                   floating_mask,            &
                                   which_ho_ground,          &
                                   f_ground_cell,            &
                                   usrf,                     &
                                   usrf_obs,                 &
                                   dthck_dt,                 &
                                   powerlaw_c_new)

    ! Compute a spatially varying basal traction field, powerlaw_c_inversion.
    ! The method is similar to that of Pollard & DeConto (TC, 2012), and is applied to all grounded ice.
    ! Where usrf > usrf_obs, powerlaw_c is reduced to increase sliding.
    ! Where usrf < usrf_obs, powerlaw_c is increased to reduce sliding.
    ! Note: powerlaw_c is constrained to lie within a prescribed range.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_inversion), intent(inout) :: &
         inversion               ! inversion object

    integer, dimension(nx,ny), intent(in) :: &
         ice_mask,             & ! = 1 where ice is present (thk > 0), else = 0
         land_mask,            & ! = 1 if topg > eus, else = 0
         floating_mask           ! = 1 where ice is present and floating, else = 0

    integer, intent(in) :: &
         which_ho_ground         ! ! option to use GLP for vertices and/or cells

    real(dp), dimension(nx,ny), intent(in) ::  &
         f_ground_cell,        & ! grounded fraction of grid cell, 0 to 1
         usrf,                 & ! ice upper surface elevation (m)
         usrf_obs,             & ! observed upper surface elevation (m)
         dthck_dt                ! rate of change of ice thickness (m/s)

    real(dp), dimension(nx,ny), intent(out) ::  &
         powerlaw_c_new          ! new value of powerlaw_c_inverison

    ! local variables

    integer, dimension(nx,ny) :: &
         powerlaw_c_inversion_mask  ! = 1 where we invert for powerlaw_c, else = 0

    real(dp), dimension(nx,ny) ::  &
         dusrf,                & ! usrf - usrf_obs on ice grid
         temp_powerlaw_c,      & ! temporary value of powerlaw_c_inversion (before smoothing)
         dpowerlaw_c             ! change in powerlaw_c

    real(dp) :: term1, term2
    real(dp) :: factor
    real(dp) :: dpowerlaw_c_smooth
    real(dp) :: sum_powerlaw_c

    integer :: i, j, ii, jj
    integer :: count

    ! parameters in inversion derived type:
    ! * powerlaw_c_max        = upper bound for powerlaw_c, Pa (m/yr)^(-1/3)
    ! * powerlaw_c_min        = lower bound for powerlaw_c, Pa (m/yr)^(-1/3)
    ! * babc_timescale        = inversion timescale (s); must be > 0
    ! * babc_thck_scale       = thickness inversion scale (m); must be > 0
    ! * babc_dthck_dt_scale   = dthck_dt inversion scale (m/s); must be > 0
    ! * babc_space_smoothing  = factor for spatial smoothing of powerlaw_c_inversion; larger => more smoothing
    !
    ! Note on babc_space_smoothing: A smoothing factor of 1/8 gives a 4-1-1-1-1 smoother.
    !       This is numerically well behaved, but may oversmooth in bowl-shaped regions;
    !        a smaller value may be better as H converges toward H_obs.

    dpowerlaw_c(:,:) = 0.0d0

    ! Compute difference between current and target upper surface elevation
    dusrf(:,:) = usrf(:,:) - usrf_obs(:,:)

    ! Compute a mask of cells where we invert for powerlaw_c.
    ! The mask includes land-based cells, as well as marine-based cells that are
    !  fully or partly grounded based on f_ground_cell.
    ! Note: f_ground_cell should be computed before transport, so that if a cell is grounded
    !       before transport and fully floating afterward, powerlaw_c_inversion is computed here
    !       rather than being set to zero.

    if (which_ho_ground == HO_GROUND_GLP_DELUXE) then

       ! Determine cells for inversion using f_ground_cell
       where (land_mask == 1 .or. f_ground_cell > 0.0d0)
          powerlaw_c_inversion_mask = 1
       elsewhere
          powerlaw_c_inversion_mask = 0
       endwhere

    else

       ! Determine cells for inversion using floating_mask
       where (land_mask == 1 .or. (ice_mask == 1 .and. floating_mask == 0))
          powerlaw_c_inversion_mask = 1
       elsewhere
          powerlaw_c_inversion_mask = 0
       endwhere

    endif   ! which_ho_ground

    call parallel_halo(powerlaw_c_inversion_mask)

    ! Check for newly grounded cells that have powerlaw_c = 0 (from when they were ice-free or floating).
    ! Give these cells a sensible default value (either land or marine).
    !TODO - This code may now be unnecessary, since powerlaw_c_save is initialized in all cells, including ice-free cells.
    do j = 1, ny
       do i = 1, nx
          if (powerlaw_c_inversion_mask(i,j) == 1) then  ! ice is land-based and/or grounded

             if (inversion%powerlaw_c_save(i,j) == 0.0d0) then
                ! set to a sensible default
                ! If on land, set to a typical land value
                ! If grounded marine ice, set to a smaller value
                if (land_mask(i,j) == 1) then
                   inversion%powerlaw_c_save(i,j) = inversion%powerlaw_c_land
                else
                   inversion%powerlaw_c_save(i,j) = inversion%powerlaw_c_marine
                endif
             endif  ! powerlaw_c_save = 0

          endif  ! powerlaw_c_inversion_mask = 1
       enddo  ! i
    enddo  ! j

    call parallel_halo(inversion%powerlaw_c_save)

    ! Loop over cells
    ! Note: powerlaw_c_inversion is computed at cell centers where usrf and thck are located.
    !       Later, it is interpolated to vertices where beta and basal velocity are located.

    do j = 1, ny
       do i = 1, nx
          if (powerlaw_c_inversion_mask(i,j) == 1) then  ! ice is land-based and/or grounded

             ! Invert for powerlaw_c based on dthck and dthck_dt
             term1 = -dusrf(i,j) / inversion%babc_thck_scale
             term2 = -dthck_dt(i,j) / inversion%babc_dthck_dt_scale

             !WHL - debug - Trying to turn off a potential unstable feedback:
             ! (1) dH/dt < 0, so Cp increases
             ! (2) Increased Cp results in dH/dt > 0, so Cp decreases
             ! (3) Amplify and repeat until the model crashes
             !TODO - Check whether this cycle occurs with a simple power law (as opposed to Schoof law)
             term2 = min(term2,  1.0d0)
             term2 = max(term2, -1.0d0)

             dpowerlaw_c(i,j) = (dt/inversion%babc_timescale) &
                  * inversion%powerlaw_c_save(i,j) * (term1 + term2)

             ! Limit to prevent huge change in one step
             if (abs(dpowerlaw_c(i,j)) > 0.05 * inversion%powerlaw_c_save(i,j)) then
                if (dpowerlaw_c(i,j) > 0.0d0) then
                   dpowerlaw_c(i,j) =  0.05d0 * inversion%powerlaw_c_save(i,j)
                else
                   dpowerlaw_c(i,j) = -0.05d0 * inversion%powerlaw_c_save(i,j)
                endif
             endif

             powerlaw_c_new(i,j) = inversion%powerlaw_c_save(i,j) + dpowerlaw_c(i,j)

             ! Limit to a physically reasonable range
             powerlaw_c_new(i,j) = min(powerlaw_c_new(i,j), inversion%powerlaw_c_max)
             powerlaw_c_new(i,j) = max(powerlaw_c_new(i,j), inversion%powerlaw_c_min)

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for powerlaw_c: rank, i, j =', rtest, itest, jtest
                print*, 'usrf, usrf_obs, dusrf, dthck_dt:', usrf(i,j), usrf_obs(i,j), dusrf(i,j), dthck_dt(i,j)*scyr
                print*, '-dusrf/usrf_scale, -dthck_dt/dthck_dt_scale, sum =', &
                     term1, term2, & 
                     term1 + term2
                print*, 'dpowerlaw_c, newpowerlaw_c =', dpowerlaw_c(i,j), powerlaw_c_new(i,j)
             endif

          else  ! powerlaw_c_inversion_mask = 0

             ! set powerlaw_c = 0
             ! Note: Zero values are ignored when interpolating powerlaw_c to vertices,
             !       and in forward runs where powerlaw_c is prescribed from a previous inversion.
             ! Warning: If a cell is grounded some of the time and floating the rest of the time,
             !           the time-averaging routine will accumulate zero values as if they are real.
             !          Time-average fields should be used with caution.
             !TODO - Ignore zero values only if no ice is present?  I.e., incorporate zeroes for floating ice.

             powerlaw_c_new(i,j) = 0.0d0

          endif  ! powerlaw_c_inversion_mask
       enddo  ! i
    enddo  ! j

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Old powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') inversion%powerlaw_c_save(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'usrf - usrf_obs:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') dusrf(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'dthck_dt (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') dthck_dt(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'powerlaw_c inversion mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') powerlaw_c_inversion_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Before smoothing, powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') powerlaw_c_new(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif
 
    if (inversion%babc_space_smoothing > 0.0d0) then

       ! Save the value just computed
       temp_powerlaw_c(:,:) = powerlaw_c_new(:,:)

       ! Apply Laplacian smoothing.
       ! Since powerlaw_c lives at cell centers but is interpolated to vertices, smoothing can damp checkerboard noise.
       !TODO - Write an operator for Laplacian smoothing?
       do j = 2, ny-1
          do i = 2, nx-1
             if (powerlaw_c_inversion_mask(i,j) == 1) then  ! ice is at least partly grounded

                dpowerlaw_c_smooth = -4.0d0 * inversion%babc_space_smoothing * temp_powerlaw_c(i,j)
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if ((ii == i .or. jj == j) .and. (ii /= i .or. jj /= j)) then  ! edge neighbor
                         if (powerlaw_c_inversion_mask(ii,jj) == 1) then  ! neighbor is at least partly grounded
                            dpowerlaw_c_smooth = dpowerlaw_c_smooth &
                                 + inversion%babc_space_smoothing*temp_powerlaw_c(ii,jj)
                         else
                            dpowerlaw_c_smooth = dpowerlaw_c_smooth &
                                 + inversion%babc_space_smoothing*temp_powerlaw_c(i,j)
                         endif
                      endif
                   enddo
                enddo

                ! Note: If smoothing is too strong, it can reverse the sign of the change in powerlaw_c.
                !       The logic below ensures that if powerlaw_c is increasing, the smoothing can reduce
                !        the change to zero, but not cause powerlaw_c to decrease relative to old_powerlaw_c
                !        (and similarly if powerlaw_c is decreasing).
                !TODO - Replace old_powerlaw_c with powerlaw_c_save?

                if (dpowerlaw_c(i,j) > 0.0d0) then
                   if (temp_powerlaw_c(i,j) + dpowerlaw_c_smooth > inversion%powerlaw_c_save(i,j)) then
                      powerlaw_c_new(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth
                   else
                      ! allow the smoothing to hold Cp at its old value, but not reduce Cp
                      powerlaw_c_new(i,j) = inversion%powerlaw_c_save(i,j)
                   endif
                elseif (dpowerlaw_c(i,j) < 0.0d0) then
                   if (temp_powerlaw_c(i,j) + dpowerlaw_c_smooth < inversion%powerlaw_c_save(i,j)) then
                      powerlaw_c_new(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth
                   else
                      ! allow the smoothing to hold Cp at its old value, but not increase Cp
                      powerlaw_c_new(i,j) = inversion%powerlaw_c_save(i,j)
                   endif
                endif  ! dpowerlaw_c > 0

             endif  ! powerlaw_c_inversion_mask = 1
          enddo   ! i
       enddo   ! j

    endif  ! smoothing factor > 0

    call parallel_halo(powerlaw_c_new)

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       if (inversion%babc_space_smoothing > 0.0d0) then
          print*, ' '
          print*, 'After smoothing, powerlaw_c:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') powerlaw_c_new(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif
    endif

  end subroutine invert_basal_traction

!***********************************************************************

  subroutine invert_bmlt_float(dt,                           &
                               nx,            ny,            &
                               itest, jtest,  rtest,         &
                               inversion,                    &
                               thck,                         &
                               usrf_obs,                     &
                               topg,                         &
                               eus,                          &
                               which_ho_ground,              &
                               which_ho_ground_bmlt,         &
                               f_ground_cell,                &
                               ice_mask,                     &
                               floating_mask,                &
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

    type(glide_inversion), intent(inout) :: &
         inversion               ! inversion object

    ! Note: thck and usrf should be the expected values after applying the mass balance
    !       (although the mass balance may not yet have been applied)
    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                 & ! ice thickness (m)
         usrf_obs,             & ! observed upper surface elevation (m)
         topg,                 & ! bedrock topography (m)
         f_ground_cell           ! fractional area of each cell which is grounded, for purposes of bmlt_float
                                 ! = 0 or 1 without a GLP; can have intermediate values with a GLP
    real(dp), intent(in) :: &
         eus                     ! eustatic sea level (m)

    integer, intent(in) :: &
         which_ho_ground,      & ! option to use GLP for vertices and/or cells
         which_ho_ground_bmlt    ! determines which cells can have nonzero bmlt_float

   ! Note: When this subroutine is called, ice_mask = 1 where thck > 0, not thck > thklim.
    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask           ! = 1 where ice is present and floating, else = 0

    real(dp), dimension(nx,ny), intent(out) ::  &
         bmlt_float_new          ! new value of bmlt_float_inversion (m/s), based on current ice geometry

    ! local variables

    integer, dimension(nx,ny) ::  &
         bmlt_float_mask         ! = 1 for cells where bmlt_float is computed and applied, else = 0
                                 ! start with inversion%bmlt_float_inversion mask; then remove grounded cells

    real(dp), dimension(nx,ny):: &
         thck_flotation,       & ! thickness at which ice becomes afloat (m)
         thck_cavity,          & ! thickness of ocean cavity beneath floating ice (m)
         thck_target             ! thickness target (m); = thck_obs unless thck_obs > thck_flotation

    real(dp) :: thck_obs         ! observed thickness
                                 ! TODO: Pass this in?

    integer :: i, j, ii, jj, iglobal, jglobal

    character(len=100) :: message

    real(dp) :: &
         dthck          ! thck - thck_target

    ! For floating cells, adjust the basal melt rate (or freezing rate, if bmlt < 0)
    !  so as to restore the upper surface to a target based on observations.

    ! Compute the flotation thickness
    where (topg - eus < 0.0d0)
       thck_flotation = -(rhoo/rhoi) *(topg - eus)
    elsewhere
       thck_flotation = 0.0d0
    endwhere

    ! Compute the ocean cavity thickness beneath floating ice (diagnostic only)
    where (floating_mask == 1)
       thck_cavity = -(topg - eus) - (rhoi/rhoo)*thck
    elsewhere
       thck_cavity = 0.0d0
    endwhere

    ! initialize
    thck_target(:,:) = 0.0d0
    bmlt_float_new(:,:) = 0.0d0

    ! Note: inversion%bmlt_float_inversion_mask is based on the initial geometry
    !       Where this mask = 0, we never invert for bmlt_float.
    !       Where this mask = 1, we invert for bmlt_float provided the cell satisfies the floating criterion.
    bmlt_float_mask(:,:) = inversion%bmlt_float_inversion_mask

    ! Identify cells that can have nonzero bmlt_float_inversion, based on which_ho_ground and which_ho_ground_bmlt

    if (which_ho_ground == HO_GROUND_GLP_DELUXE) then

       ! Determine bmlt_float_mask using f_ground_cell

       if (which_ho_ground_bmlt == HO_GROUND_BMLT_FLOATING_FRAC) then

          ! No basal melt in fully grounded cells
          !TODO - Allow inversion where ice_mask = 0 but ice is present in obs? 
          where (ice_mask == 0 .or. f_ground_cell >= 1.0d0 - tiny(0.0d0))
             bmlt_float_mask = 0
          endwhere

       elseif (which_ho_ground_bmlt == HO_GROUND_BMLT_ZERO_GROUNDED) then

          ! No basal melt in cells that are even partly grounded
          where (ice_mask == 0 .or. f_ground_cell > tiny(0.0d0))
             bmlt_float_mask = 0
          endwhere

       endif  ! which_ho_ground_bmlt

    else

       ! No basal melt in grounded cells
       where (floating_mask == 0)
          bmlt_float_mask = 0
       endwhere

    endif   ! which_ho_ground


    ! For cells with bmlt_float_mask = 1, compute bmlt_float_inversion that will restore the thickness
    !  to the observed target.
    ! TODO: If not using the basal melting GLP, then restoring all the way to the grounded target will lead to oscillations,
    !       and we may need to use a buffer to prevent over-restoring.  For now, focus on runs with a basal melting GLP.

    ! loop over cells
    do j = 1, ny
       do i = 1, nx

          if (bmlt_float_mask(i,j) == 1) then  ! at least partly floating

             if (usrf_obs(i,j) - (topg(i,j) - eus) > thck_flotation(i,j)) then  ! grounded target

                !TODO - Pass in thck_obs
                thck_obs = usrf_obs(i,j) - (topg(i,j) - eus)

                ! target is the observed thickness (assuming topg is fixed)
                thck_target(i,j) = thck_obs

                ! Do not allow the target to be too strongly grounded
                thck_target(i,j) = min(thck_target(i,j), &
                                       thck_flotation(i,j) + inversion%bmlt_max_thck_above_flotation*thk0)

             else  ! floating target
                thck_target(i,j) = usrf_obs(i,j) * rhoo/(rhoo - rhoi)  ! target is the observed thickness
                thck_target(i,j) = min(thck_target(i,j), thck_flotation(i,j) - inversion%thck_flotation_buffer*thk0)
             endif

             thck_target(i,j) = max(thck_target(i,j), 0.0d0)

             dthck = thck(i,j) - thck_target(i,j)

             !WHL - debug - Some cells are thinning dynamically, and it requires a very large freezing rate
             !               (bmlt_float < 0) to restore them to the target thickness.
             !              For these cells, limit the freezing rate to a prescribed maximum, thus allowing
             !               the cells to remain thinner than the target.
             ! Note: inversion%bmlt_freeze_max has been scaled to have units of m/s

             if (inversion%bmlt_freeze_max*scyr > eps08) then
                if (dthck < 0.0d0) then
                   ! Reduces to -bmlt_float_freeze_max when dthck/dt is large, and to dthck/dt when dthck/dt is small.
                   ! Minus sign because bmlt_float is positive for melting
                   bmlt_float_new(i,j) = -inversion%bmlt_freeze_max * &
                        (1.0d0 - exp((dthck/dt)/inversion%bmlt_freeze_max))
                else
                   bmlt_float_new(i,j) = dthck / dt
                endif
             else
                bmlt_float_new(i,j) = dthck / dt
             endif

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for bmlt_float_inversion: rank, i, j =', rtest, itest, jtest
                print*, 'topg - eus, usrf_obs:', topg(i,j) - eus, usrf_obs(i,j)
                print*, 'thck, buffer, thck_target:', thck(i,j), inversion%thck_flotation_buffer, thck_target(i,j)
                print*, 'dthck/dt, bmlt_float_new:', (dthck/dt)*scyr, bmlt_float_new(i,j)*scyr
             endif

          endif   ! bmlt_float_mask = 1

       enddo   ! i
    enddo   ! j

    call parallel_halo(bmlt_float_mask) ! diagnostic only
    call parallel_halo(thck_target)  ! diagnostic only
    call parallel_halo(bmlt_float_new)
 
    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'bmlt_float inversion mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') bmlt_float_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'f_ground_cell:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') f_ground_cell(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'floating_mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') floating_mask(i,j)
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
       print*, 'expected thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck_target (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') thck_target(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'new (unlimited) bmlt_float inversion (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') bmlt_float_new(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine invert_bmlt_float

!=======================================================================

end module glissade_inversion

!=======================================================================
