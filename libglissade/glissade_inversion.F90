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
  use cism_parallel, only: this_rank, main_task, nhalo, &
       parallel_type, parallel_halo, staggered_parallel_halo, &
       parallel_reduce_min, parallel_reduce_max

  implicit none

  private
  public :: verbose_inversion, glissade_init_inversion, glissade_inversion_basal_friction, &
            glissade_inversion_bmlt_basin, glissade_inversion_flow_factor_basin

  !-----------------------------------------------------------------------------
  ! Subroutines to invert for basal fields (including basal friction beneath
  ! grounded ice and basal melting beneath floating ice) by relaxing toward
  ! a target ice thickness field.
  !-----------------------------------------------------------------------------

!!    logical, parameter :: verbose_inversion = .false.
    logical, parameter :: verbose_inversion = .true.

!***********************************************************************

contains

!***********************************************************************

  subroutine glissade_init_inversion(model)

    ! Initialize inversion for fields of basal friction and basal melting
    ! Should be called after usrf and thck have been input and (possibly) modified by initial calving

    use glissade_masks, only: glissade_get_masks
    use glissade_grid_operators, only: glissade_stagger

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
         land_mask               ! = 1 where topg is at or above sea level

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck_obs                ! observed ice thickness, derived from usrf_obs and topg

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         coulomb_c_icegrid      ! initial coulomb_c at cell centers based on masks

    real(dp) :: f_flotation                  ! flotation function (m); < 0 for grounded ice, > 0 for floating ice
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
    ! If inverting for Cp or Cc, then set the target elevation, usrf_obs,
    !  and the target surface ice speed, velo_sfc_obs.
    ! Note: Must read in usfc_obs and vsfc_obs to set velo_sfc_obs correctly.
    !       Typically, the inversion is based only on the surface elevation, usrf_obs,
    !        but compute velo_sfc_obs regardless since it is a useful diagnostic.
    !----------------------------------------------------------------------

    if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or.  &
        model%options%which_ho_coulomb_c  == HO_COULOMB_C_INVERSION) then

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
          ! Probably not necessary if doing basin-scale inversion for floating ice instead of inversion in each grid cell.
          !WHL - debug
          do j = nhalo+1, nsn-nhalo
             do i = nhalo+1, ewn-nhalo
                if (verbose_inversion .and. thck_obs(i,j) > 0.0d0 .and. &
                     thck_obs(i,j) < model%inversion%thck_threshold) then
                   print*, 'thck_obs < threshold, rank, i, j, thck:', this_rank, i, j, thck_obs(i,j)*thk0
                endif
             enddo
          enddo

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

       call parallel_halo(model%geometry%usrf_obs, parallel)
       call parallel_halo(thck_obs, parallel)

       ! Set the surface speed target, velo_sfc_obs
       if (model%options%is_restart == RESTART_FALSE) then
          model%velocity%velo_sfc_obs(:,:) = &
               sqrt(model%velocity%usfc_obs(:,:)**2 + model%velocity%vsfc_obs(:,:)**2)
       endif

       ! If inverting based on a velocity target, check that nonzero values were read in
       if (model%inversion%babc_velo_scale > 0.0d0) then
          var_maxval = maxval(model%velocity%velo_sfc_obs)
          var_maxval = parallel_reduce_max(var_maxval)
          if (var_maxval == 0.0d0) then
             call write_log &
                  ('Error: velo_sfc_obs = 0 everywhere, when babc_velo_scale > 0', GM_FATAL)
          endif
       endif

    endif  ! inversion for Cp or Cc

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

       ! initialize powerlaw_c_inversion, if not already read in
       var_maxval = maxval(model%basal_physics%powerlaw_c)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; powerlaw_c_inversion has been read in already (e.g., when restarting)
       else
          ! initialize to a uniform value (which can be set in the config file)
          model%basal_physics%powerlaw_c(:,:) = model%basal_physics%powerlaw_c_const
       endif  ! var_maxval > 0

       if (verbose_inversion .and. this_rank == rtest) then
          print*, ' '
          print*, 'glissade_init_inversion: powerlaw_c_inversion:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.1)',advance='no') model%basal_physics%powerlaw_c(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    elseif (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION) then

       ! initialize coulomb_c_inversion, if not already read in
       var_maxval = maxval(model%basal_physics%coulomb_c)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; coulomb_c has been read in already (e.g., when restarting)
       else
          ! Set a low initial value for cells that are floating or ice-free ocean
          ! Set a higher value for cells that are ground ice and/or land-covered
          where (ocean_mask == 1 .or. floating_mask == 1)
             coulomb_c_icegrid = model%basal_physics%coulomb_c_min
          elsewhere
             coulomb_c_icegrid = model%basal_physics%coulomb_c_const
          endwhere

          ! Interpolate to the staggered grid
          call glissade_stagger(ewn,         nsn,       &
                                coulomb_c_icegrid,      &
                                model%basal_physics%coulomb_c)
       endif  ! var_maxval > 0

       if (verbose_inversion .and. this_rank == rtest) then
          print*, ' '
          print*, 'glissade_init_inversion: coulomb_c:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%basal_physics%coulomb_c(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    endif   ! Cp or Cc inversion

    !----------------------------------------------------------------------
    ! computations specific to inversion of deltaT_basin or flow_factor_basin
    !----------------------------------------------------------------------

    if (model%options%which_ho_bmlt_basin == HO_BMLT_BASIN_INVERSION .or. &
        model%options%which_ho_flow_factor_basin == HO_FLOW_FACTOR_BASIN_INVERSION) then

       if (model%options%is_restart == RESTART_FALSE) then

          ! Set floating_thck_target for floating ice and lightly grounded ice.
          ! Here, "lightly grounded" means that the magnitude of f_flotation = (-topg - eus) - (rhoi/rhoo)*thck
          !  is less than a prescribed threshold.  (Recall f_flotation < 0 for grounded ice.)
          ! The inversion will nudge the ice thickness toward this target in a basin-average sense.
          ! Positive volume biases will be corrected with ocean warming or ice softening,
          !  and negative biases with ocean cooling or ice stiffening.

          do j = 1, nsn
             do i = 1, ewn
                f_flotation = (-(model%geometry%topg(i,j) - model%climate%eus)  &
                              - (rhoi/rhoo)*model%geometry%thck(i,j)) * thk0    ! f_flotation < 0 for grounded ice
                if (model%geometry%thck(i,j) > 0.0d0 .and. &
                    model%geometry%marine_connection_mask(i,j) == 1 .and. &
                    f_flotation > -model%inversion%basin_flotation_threshold) then
                   model%inversion%floating_thck_target(i,j) = model%geometry%thck(i,j)
                else
                   model%inversion%floating_thck_target(i,j) = 0.0d0
                endif
             enddo
          enddo

          ! If inverting for the flow factor, then initialize to a constant
          if (model%options%which_ho_flow_factor_basin == HO_FLOW_FACTOR_BASIN_INVERSION) then
             model%temper%flow_factor_basin(:,:) = model%paramets%flow_enhancement_factor_float
          endif

          if (verbose_inversion .and. this_rank == rtest) then
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
             print*, 'f_flotation (m):'
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

       call parallel_halo(model%inversion%floating_thck_target, parallel)

    endif  ! which_ho_bmlt_basin_inversion or which_ho_flow_factor_basin_inversion

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

  subroutine glissade_inversion_basal_friction(model)

    use glimmer_paramets, only: tim0, thk0, vel0
    use glimmer_physcon, only: scyr, grav
    use glissade_grid_operators, only: glissade_stagger, glissade_stagger_real_mask

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! --- Local variables ---

    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
         thck_obs                ! observed ice thickness, derived from usrf_obs and topg

    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) ::   &
         stag_thck,            & ! ice thickness on staggered grid
         stag_dthck_dt,        & ! dthck_dt on staggered grid
         stag_thck_obs,        & ! thck_obs on staggered grid
         velo_sfc                ! surface ice speed

    !WHL - debug
    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) ::   &
         stag_topg,            &
         stag_thck_flotation

    integer :: i, j
    integer :: ewn, nsn
    integer :: itest, jtest, rtest

    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_unscaled

    logical :: &
         f_ground_weight = .true.   ! if true, then weigh ice thickness by f_ground_cell for staggered interpolation
                                    ! Found that unweighted staggering can lead to low-frequency thickness oscillations
                                    !  in Antarctic runs, because of large dH/dt in floating cells

    logical :: invert_coulomb_c, invert_powerlaw_c

    type(parallel_type) :: parallel

    parallel = model%parallel

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

    ! Set logical variables

    invert_coulomb_c = .false.
    invert_powerlaw_c = .false.

    if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION) then
       invert_powerlaw_c = .true.
    elseif (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION) then
       invert_coulomb_c = .true.
    endif

    if (invert_powerlaw_c .or. invert_coulomb_c) then

       ! Compute the new value of powerlaw_c or coulomb_c at each vertex

       ! Given the surface elevation target, compute the thickness target.
       ! (This can change in time if the bed topography is dynamic.)
       call usrf_to_thck(model%geometry%usrf_obs,  &
                         model%geometry%topg,      &
                         model%climate%eus,        &
                         thck_obs)

       ! Interpolate the thickness fields to the staggered grid

       if (f_ground_weight) then  ! give a greater weight to cells that are fully grounded

          ! Interpolate thck_obs to the staggered grid
          call glissade_stagger_real_mask(&
               ewn,         nsn,               &
               thck_obs,    stag_thck_obs,     &
               model%geometry%f_ground_cell)

          ! Interpolate thck to the staggered grid
          call glissade_stagger_real_mask(&
               ewn,                  nsn,       &
               model%geometry%thck,  stag_thck, &
               model%geometry%f_ground_cell)

          ! Interpolate dthck_dt to the staggered grid
          call glissade_stagger_real_mask(&
               ewn,                      nsn,           &
               model%geometry%dthck_dt,  stag_dthck_dt, &
               model%geometry%f_ground_cell)

       else   ! equally weight the values in all four neighbor cells, including ice-free cells

          ! Interpolate thck_obs to the staggered grid
          call glissade_stagger(ewn,         nsn,              &
                                thck_obs,    stag_thck_obs)

          ! Interpolate thck to the staggered grid
          call glissade_stagger(ewn,                  nsn,             &
                                model%geometry%thck,  stag_thck)

          ! Interpolate dthck_dt to the staggered grid
          call glissade_stagger(ewn,                      nsn,             &
                                model%geometry%dthck_dt,  stag_dthck_dt)

       endif   ! f_ground_weight

       call staggered_parallel_halo(stag_thck_obs, parallel)
       call staggered_parallel_halo(stag_thck, parallel)
       call staggered_parallel_halo(stag_dthck_dt, parallel)

       ! Given the ice velocity, compute the surface speed

       velo_sfc(:,:) = sqrt(model%velocity%uvel(1,:,:)**2 + model%velocity%vvel(1,:,:)**2)
       call staggered_parallel_halo(velo_sfc, parallel)

       ! Invert for powerlaw_c or coulomb_c
       ! The logic is the same for each; only the max and min values and the in/out field are different.

       if (invert_powerlaw_c) then

          if (verbose_inversion .and. this_rank == rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'Old powerlaw_c:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.2)',advance='no') model%basal_physics%powerlaw_c(i,j)
                enddo
                print*, ' '
             enddo
          endif   ! verbose_inversion

          call invert_basal_friction(model%numerics%dt*tim0,                   &  ! s
                                     ewn,               nsn,                   &
                                     itest,    jtest,   rtest,                 &
                                     model%inversion%babc_timescale,           &  ! s
                                     model%inversion%babc_thck_scale,          &  ! m
                                     model%inversion%babc_velo_scale,          &  ! m/yr
                                     model%basal_physics%powerlaw_c_max,       &
                                     model%basal_physics%powerlaw_c_min,       &
                                     model%geometry%f_ground,                  &
                                     stag_thck*thk0,                           &  ! m
                                     stag_thck_obs*thk0,                       &  ! m
                                     stag_dthck_dt,                            &  ! m/s
                                     velo_sfc*(vel0*scyr),                     &  ! m/yr
                                     model%velocity%velo_sfc_obs*(vel0*scyr),  &  ! m/yr
                                     model%basal_physics%powerlaw_c)

          if (verbose_inversion .and. this_rank == rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'New powerlaw_c:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.2)',advance='no') model%basal_physics%powerlaw_c(i,j)
                enddo
                print*, ' '
             enddo
          endif   ! verbose_inversion

       elseif (invert_coulomb_c) then

          if (verbose_inversion .and. this_rank == rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'Old coulomb_c:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.5)',advance='no') model%basal_physics%coulomb_c(i,j)
                enddo
                print*, ' '
             enddo
             print*, ' '
             print*, 'effecpress/overburden:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   if (stag_thck(i,j) > 0.0d0) then
                      write(6,'(f10.4)',advance='no') &
                           model%basal_physics%effecpress_stag(i,j) / &
                           (rhoi * grav * stag_thck(i,j) * thk0)
                   else
                      write(6,'(f10.4)',advance='no') 0.0d0
                   endif
                enddo
                print*, ' '
             enddo

          endif   ! verbose_inversion

          ! Compute flotation thickness, given by H = (rhoo/rhoi)*|b|

          ! Interpolate topg to the staggered grid
          call glissade_stagger(ewn,                   nsn,             &
                                model%geometry%topg,   stag_topg)

          ! correct for eus (if nonzero) and convert to meters
          stag_topg = (stag_topg - model%climate%eus) * thk0
          stag_thck_flotation = (rhoo/rhoi) * max(-stag_topg, 0.0d0)

          call invert_basal_friction(model%numerics%dt*tim0,                   &  ! s
                                     ewn,               nsn,                   &
                                     itest,    jtest,   rtest,                 &
                                     model%inversion%babc_timescale,           &  ! s
                                     model%inversion%babc_thck_scale,          &  ! m
                                     model%inversion%babc_velo_scale,          &  ! m/yr
                                     model%basal_physics%coulomb_c_max,        &
                                     model%basal_physics%coulomb_c_min,        &
                                     model%geometry%f_ground,                  &
                                     stag_thck*thk0,                           &  ! m
                                     stag_thck_obs*thk0,                       &  ! m
                                     stag_dthck_dt,                            &  ! m/s
                                     velo_sfc*(vel0*scyr),                     &  ! m/yr
                                     model%velocity%velo_sfc_obs*(vel0*scyr),  &  ! m/yr
                                     model%basal_physics%coulomb_c,            &
                                     stag_thck_flotation = stag_thck_flotation, &
                                     p_ocean = model%basal_physics%p_ocean_penetration)

          if (verbose_inversion .and. this_rank == rtest) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'New coulomb_c:'
             do j = jtest+3, jtest-3, -1
                do i = itest-3, itest+3
                   write(6,'(f10.5)',advance='no') model%basal_physics%coulomb_c(i,j)
                enddo
                print*, ' '
             enddo
          endif   ! verbose_inversion

       endif  ! invert for powerlaw_c or coulomb_c

    else   ! do not invert for powerlaw_c or coulomb_c; just print optional diagnostics

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
          print*, 'powerlaw_c:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') model%basal_physics%powerlaw_c(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'coulomb_c:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.4)',advance='no') model%basal_physics%coulomb_c(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    endif   ! invert_powerlaw_c or invert_coulomb_c

    ! Replace zeroes (if any) with small nonzero values to avoid divzeroes.
    ! Note: The current algorithm initializes Cc to a nonzero value everywhere and never sets Cc = 0;
    !       this code is just to be on the safe side.

    if (model%options%which_ho_powerlaw_c /= HO_POWERLAW_C_CONSTANT) then
       where (model%basal_physics%powerlaw_c == 0.0d0)
          model%basal_physics%powerlaw_c = model%basal_physics%powerlaw_c_min
       endwhere
    endif

    if (model%options%which_ho_powerlaw_c /= HO_COULOMB_C_CONSTANT) then
       where (model%basal_physics%coulomb_c == 0.0d0)
          model%basal_physics%coulomb_c = model%basal_physics%coulomb_c_min
       endwhere
    endif

  end subroutine glissade_inversion_basal_friction

!***********************************************************************

  subroutine invert_basal_friction(dt,                       &
                                   nx,            ny,        &
                                   itest, jtest,  rtest,     &
                                   babc_timescale,           &
                                   babc_thck_scale,          &
                                   babc_velo_scale,          &
                                   friction_c_max,           &
                                   friction_c_min,           &
                                   f_ground,                 &
                                   stag_thck,                &
                                   stag_thck_obs,            &
                                   stag_dthck_dt,            &
                                   velo_sfc,                 &
                                   velo_sfc_obs,             &
                                   friction_c,               &
                                   stag_thck_flotation,      &
                                   p_ocean)

    ! Compute a spatially varying basal friction field defined at cell vertices.
    ! Here, the field has the generic name 'friction_c', which could be either powerlaw_c or coulomb_c.
    ! The method is similar to that of Pollard & DeConto (TC, 2012), and is applied to all grounded ice.
    ! Adjustments can be based on a thickness target, a surface ice target, or both:
    !    Where stag_thck > stag_thck_obs, friction_c is reduced to increase sliding.
    !    Where stag_thck < stag_thck_obs, friction_c is increased to reduce sliding.
    !    Where velo_sfc > velo_sfc_obs, friction_c is increased to reduce sliding.
    !    Where velo_sfc < velo_sfc_obs, friction_c is decreased to increase sliding.
    ! Note: friction_c is constrained to lie within a prescribed range, [friction_c_min, friction_c_max].
    ! Note: For grounded ice with fixed topography, inversion based on thck is equivalent to inversion based on usrf.
    !       But for ice that is partly floating, it seems better to invert based on thck, because thck errors
    !        are greater in magnitude than usrf errors, and we do not want to underweight the errors.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         babc_timescale,       & ! inversion timescale (s); must be > 0
         babc_thck_scale,      & ! thickness inversion scale (m)
         babc_velo_scale,      & ! velocity inversion scale (m/yr)
         friction_c_max,       & ! upper bound for friction_c (units correspond to powerlaw_c or coulomb_c)
         friction_c_min          ! lower bound for friction_c

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
         f_ground,             & ! grounded fraction at vertices, 0 to 1
         stag_thck,            & ! ice thickness at vertices (m)
         stag_thck_obs,        & ! observed ice thickness at vertices (m)
         stag_dthck_dt,        & ! rate of change of ice thickness at vertices (m/s)
         velo_sfc,             & ! ice surface speed at vertices (m/yr)
         velo_sfc_obs            ! observed ice surface speed at vertices (m/yr)

    real(dp), dimension(nx-1,ny-1), intent(inout) ::  &
         friction_c              ! basal friction field to be adjusted (powerlaw_c or coulomb_c)

    real(dp), dimension(nx-1,ny-1), intent(in), optional :: &
         stag_thck_flotation     ! flotation thickness (m) on staggered grid; used for term3_thck

    real(dp), intent(in), optional :: p_ocean

    ! local variables

    real(dp), dimension(nx-1,ny-1) ::  &
         stag_dthck,           & ! stag_thck - stag_thck_obs
         dvelo_sfc,            & ! velo_sfc - velo_sfc_obs
         dfriction_c             ! change in friction_c

    real(dp), dimension(nx-1, ny-1) :: &
         term1_thck, term2_thck, term3_thck, & ! tendency terms based on thickness target
         term1_velo                            ! tendency term based on surface speed target

    integer :: i, j

    ! Initialize
    dfriction_c(:,:) = 0.0d0

    ! Compute difference between current and target thickness and surface speed
    ! Note: Where the target cell is ice-free, stag_dthck will be > 0, to encourage thinning.
    !       Where the target speed = 0 (because of missing data, or because the target
    !        is ice-free), there is no nudging toward a target speed.

    stag_dthck(:,:) = stag_thck(:,:) - stag_thck_obs(:,:)

    where (velo_sfc_obs > 0.0d0)
       dvelo_sfc = velo_sfc - velo_sfc_obs
    elsewhere
       dvelo_sfc = 0.0d0
    endwhere

    ! optional diagnostics
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
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
             write(6,'(f10.3)',advance='no') stag_dthck_dt(i,j)*scyr
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'velo_sfc - velo_sfc_obs (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') dvelo_sfc(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'f_ground'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') f_ground(i,j)
          enddo
          print*, ' '
       enddo
    endif

    ! Initialize the tendency terms
    term1_thck(:,:) = 0.0d0
    term2_thck(:,:) = 0.0d0
    term3_thck(:,:) = 0.0d0
    term1_velo(:,:) = 0.0d0

    ! Loop over vertices where f_ground > 0
    ! Note: f_ground should be computed before transport, so that if a vertex is grounded
    !       before transport and fully floating afterward, friction_c is computed here.

    do j = 1, ny-1
       do i = 1, nx-1

          if (f_ground(i,j) > 0.0d0) then  ! ice is at least partly grounded

             ! Compute the rate of change of friction_c, based on stag_dthck and stag_dthck_dt,
             !  and/or dvelo_sfc.
             ! For a thickness target, the rate of change is proportional to the sum of two terms:
             !     dC/dt = -C * (1/tau) * [(H - H_obs)/H0 + (2*tau/H0) * dH/dt]
             ! where tau = babc_timescale and H0 = babc_thck_scale.
             ! This equation is similar to that of a damped harmonic oscillator:
             !     m * d2x/dt2 = -k*x - c*dx/dt
             ! where m is the mass, k is a spring constant, and c is a damping term.
             ! A harmonic oscillator is critically damped when c = 2*sqrt(m*k).
             !  In this case the system reaches equilibrium as quickly as possible without oscillating.
             ! Assuming unit mass (m = 1) and critical damping with k = 1/(tau^2), we obtain
             !   d2x/dt2 = -1/tau * (x/tau - 2*dx/dt)
             ! If we identify (H - H_obs)/(H0*tau) with x/tau; (2/H0)*dH/dt with 2*dx/dt; and (1/C)*dC/dt with d2x/dt2,
             !  we obtain the equation solved here.
             ! With a surface speed target (babc_velo_scale > 0), we add a term proportional to (u - u_obs)/u0.
             ! However, there is no tendency term associated with velocity changes du/dt.
             ! Note: babc_thck_scale and babc_velo_scale have default values of 0.
             !       Setting either or both to positive values in the config file will activate the inversion.

             ! Compute tendency terms based on the thickness target

             if (babc_thck_scale > 0.0d0) then
                term1_thck(i,j) = -stag_dthck(i,j) / (babc_thck_scale * babc_timescale)
                term2_thck(i,j) = -stag_dthck_dt(i,j) * 2.0d0 / babc_thck_scale
             endif

             ! Third tendency term added for coulomb_c inversion.
             ! The origin of the term is as follows:  Basal shear stress is proportional to N * C_c.
             ! We want to relax (N * C_c) such that H approaches a steady value without oscillating.
             !  The prognostic equation is:
             !      1/(N*C_c) * d(N*C_c)/dt = -(1/tau) * [(H - H_obs)/H0 + (2*tau/H0) * dH/dt]
             ! Using the product rule on the LHS gives a term of the form (1/N)(dN/dt).
             ! Move this term to the RHS and set dN/dt = dN/dh * dh/dt.
             ! With N = (rhoi*g*H) * (1 - Hf/H)^p, we can show dN/dh = N * [(1 - p)/H + p/(H - Hf)],
             ! giving the term below.
             ! The result is to increase the rate of change of C_c near the grounding line
             !  when ice is thinning and retreating.  Ideally, this should help with stability
             !  by reducing the chance of overshooting the GL.  In practice, it may only delay
             !  the onset of instability, if the thickness target is an unstable state.

             if (present(p_ocean)) then
                if (stag_thck(i,j) > 0.0d0 .and. &
                     stag_thck(i,j) > stag_thck_flotation(i,j)) then
                   term3_thck(i,j) = -stag_dthck_dt(i,j) *  &
                        ( (1.0d0 - p_ocean)/stag_thck(i,j) &
                        + p_ocean / (stag_thck(i,j) - stag_thck_flotation(i,j)) )
                endif
             endif

             ! Compute tendency terms based on the surface speed target
             ! Note: I tried adding a term2_velo in analogy to term2_thck (Dec. 2021),
             !       but it triggers oscillations in friction_c without improving accuracy.

             if (babc_velo_scale > 0.0d0) then
                term1_velo(i,j) = dvelo_sfc(i,j) / (babc_velo_scale * babc_timescale)
             endif

             if (verbose_inversion .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                print*, ' '
                print*, 'Increment C_c, rank, i, j =', this_rank, i, j
                print*, 'dt*term1_thck =', dt*term1_thck(i,j)
                print*, 'dt*term2_thck =', dt*term2_thck(i,j)
                if (present(p_ocean)) print*, 'dt*term3_thck =', dt*term3_thck(i,j)
                if (babc_velo_scale > 0.0d0) print*, 'dt*term1_velo =', dt*term3_thck(i,j)
             endif

             dfriction_c(i,j) = friction_c(i,j) * &
                  (term1_thck(i,j) + term2_thck(i,j) + term3_thck(i,j) + term1_velo(i,j)) * dt

             ! Limit to prevent a large relative change in one step
             if (abs(dfriction_c(i,j)) > 0.05d0 * friction_c(i,j)) then
                if (dfriction_c(i,j) > 0.0d0) then
                   dfriction_c(i,j) =  0.05d0 * friction_c(i,j)
                else
                   dfriction_c(i,j) = -0.05d0 * friction_c(i,j)
                endif
             endif

             ! Update friction_c
             friction_c(i,j) = friction_c(i,j) + dfriction_c(i,j)

             ! Limit to a physically reasonable range
             friction_c(i,j) = min(friction_c(i,j), friction_c_max)
             friction_c(i,j) = max(friction_c(i,j), friction_c_min)

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for friction_c: rank, i, j =', rtest, itest, jtest
                print*, 'thck, thck_obs, dthck, dthck_dt:', &
                     stag_thck(i,j), stag_thck_obs(i,j), stag_dthck(i,j), stag_dthck_dt(i,j)*scyr
                print*, 'velo_sfc, velo_sfc_obs, dvelo_sfc:', velo_sfc(i,j), velo_sfc_obs(i,j), dvelo_sfc(i,j)
                print*, 'dthck term, dthck_dt term, sum =', &
                     term1_thck(i,j)*dt, term2_thck(i,j)*dt, (term1_thck(i,j) + term2_thck(i,j))*dt
                if (present(p_ocean)) print*, 'dN/dH term:', term3_thck(i,j)*dt
                if (babc_velo_scale > 0.0d0) print*, 'dvelo term =', term1_velo(i,j)*dt
                print*, 'dfriction_c, new friction_c =', dfriction_c(i,j), friction_c(i,j)
             endif

          else   ! f_ground = 0

             ! do nothing; keep the old value

          endif  ! f_ground > 0

       enddo  ! i
    enddo  ! j

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
                                           basin_mass_correction,       &
                                           basin_number_mass_correction,&
                                           dbmlt_dtemp_scale,           &
                                           bmlt_basin_timescale,        &
                                           deltaT_basin)

    ! For the case that bmlt_float is computed based on thermal_forcing,
    !  adjust deltaT_basin, which can be thought of as a bias corrrection
    !  or tuning parameter for the thermal forcing parameterization.
    ! In each basin, we compute the mean thickness of floating or lightly grounded ice
    !  and compare to a target thickness (usually based on observations).
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
         basin_mass_correction   ! optional mass correction (Gt) for a selected basin

    integer, intent(in) :: &
         basin_number_mass_correction ! integer ID for the basin receiving the correction

    real(dp), intent(in) :: &
         dbmlt_dtemp_scale,    & ! scale for rate of change of bmlt w/temperature, (m/s)/degC
         bmlt_basin_timescale    ! timescale for adjusting deltaT_basin (s)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         deltaT_basin            ! deltaT correction to thermal forcing in each basin (deg C)

    ! local variables

    real(dp), dimension(nbasin) :: &
         floating_thck_target_basin,      &   ! floating mean thickness target in each basin (m^3)
         floating_thck_basin,             &   ! current mean ice thickness in each basin (m)
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

    ! TODO: Make these config parameters
    real(dp), parameter :: &
         deltaT_basin_maxval = 2.0d0,      & ! max allowed magnitude of deltaT_basin (deg C)
         dTbasin_dt_maxval = 1.0d0/scyr      ! max allowed magnitude of d(deltaT_basin)/dt (deg/yr converted to deg/s)

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

    ! Determine the rate of change of deltaT_basin for each basin.
    ! Warm the basin where the ice is too thick, and cool where the ice is too thin.
    ! Note: deltaT_basin is a 2D field, but its value is uniform in each basin.

    do nb = 1, nbasin
       term1 = (1.0d0/dbmlt_dtemp_scale) * &
            (floating_thck_basin(nb) - floating_thck_target_basin(nb)) / (bmlt_basin_timescale**2)
       term2 = (1.0d0/dbmlt_dtemp_scale) * 2.0d0 * floating_dthck_dt_basin(nb) / bmlt_basin_timescale
       dTbasin_dt(nb) = term1 + term2
    enddo

    ! Limit dTbasin/dt to a prescribed range
    ! This prevents rapid changes in basins with small volume targets.
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

       if (main_task) then
          print*, ' '
          print*, 'basin, term1*dt, term2*dt, dTbasin, new deltaT_basin:'
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

  subroutine glissade_inversion_flow_factor_basin(&
       dt,                           &
       nx,            ny,            &
       dx,            dy,            &
       itest, jtest,  rtest,         &
       nbasin,                       &
       basin_number,                 &
       thck,                         &
       dthck_dt,                     &
       floating_thck_target,         &
       basin_mass_correction,        &
       basin_number_mass_correction, &
       flow_factor_basin_thck_scale, &
       flow_factor_basin_timescale,  &
       flow_factor_basin)

    ! For the case that we are trying to match a thickness target for floating ice
    !  in each basin, adjust a basin-specific flow enhancement factor called flow_factor_basin.
    ! The logic is similar to that for deltaT_basin inversion above.
    ! In each basin, we compute the mean thickness of floating or lightly grounded ice
    !  and compare to a target thickness (usually based on observations).
    ! Where the ice is thicker than the target, we increase flow_factor_basin,
    !  which makes the ice less viscous and faster-flowing.
    ! Where the ice is thinner than the targer, we reduce flow_factor_basin.

    use glissade_utils, only: glissade_basin_sum

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
         basin_mass_correction            ! optional mass correction (Gt) for a selected basin

    integer, intent(in) :: &
         basin_number_mass_correction     ! integer ID for the basin receiving the correction

    real(dp), intent(in) :: &
         flow_factor_basin_thck_scale,  & ! thickness scale for adjusting flow_factor_basin (s)
         flow_factor_basin_timescale      ! timescale for adjusting flow_factor_basin (s)

    real(dp), dimension(nx,ny), intent(inout) ::  &
         flow_factor_basin                ! flow enhancement factor for floating ice in each basin (unitless)

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
         dflow_factor_basin_dt,           &   ! rate of change of flow_factor_basin (1/s)
         basin_max, basin_min,            &   ! min and max of flow_factor_basin in each basin
                                              ! (all cells in the basin should have the same value of flow_factor_basin)
         flow_factor_basin_nb                 ! same as flow_factor_basin, but with dimension nbasin

    integer :: i, j
    integer :: nb      ! basin number
    real(dp) :: term1, term2

    ! Note: Max and min values are somewhat arbitrary.
    ! TODO: Make these config parameters
    real(dp), parameter :: &
         flow_factor_basin_maxval = 5.0d0,  & ! max allowed magnitude of flow_factor_basin (unitless)
         flow_factor_basin_minval = 0.2d0,  & ! min allowed magnitude of flow_factor_basin (unitless)
         dflow_factor_basin_dt_maxval = 0.1d0/scyr  ! max allowed magnitude of d(flow_factor_basin)/dt (1/yr converted to 1/s)

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

    ! Diagnose the current flow_factor_basin for each basin.
    ! This assumes that all cells in a basin have the same value.

    basin_max(:) = 0.0d0
    do j = 1, ny
       do i = 1, nx
          nb = basin_number(i,j)
          if (nb >= 1 .and. nb <= nbasin) then
             basin_max(nb) = max(basin_max(nb), flow_factor_basin(i,j))
          endif
       enddo
    enddo

    do nb = 1, nbasin
       flow_factor_basin_nb(nb) = parallel_reduce_max(basin_max(nb))
    enddo

    ! Determine the rate of change of flow_factor_basin for each basin.
    ! Raise the factor (i.e., lower the viscosity) where the ice is too thick,
    !  and lower the factor (raise the viscosity) where the ice is too thin.
    ! The prognostic equation is similar to that for coulomb_c, another scalar.
    ! Note: flow_factor_basin is a 2D field, but its value is uniform in each basin.

    do nb = 1, nbasin
       term1 = (floating_thck_basin(nb) - floating_thck_target_basin(nb)) / &
            (flow_factor_basin_thck_scale * flow_factor_basin_timescale)
       term2 = 2.0d0 * floating_dthck_dt_basin(nb) / flow_factor_basin_thck_scale
       dflow_factor_basin_dt(nb) = flow_factor_basin_nb(nb) * (term1 + term2)
    enddo

    ! Limit flow_factor_basin/dt to a prescribed range
    ! This prevents rapid changes in basins with small volume targets.
    where (dflow_factor_basin_dt > dflow_factor_basin_dt_maxval)
       dflow_factor_basin_dt = dflow_factor_basin_dt_maxval
    elsewhere (dflow_factor_basin_dt < -dflow_factor_basin_dt_maxval)
       dflow_factor_basin_dt = -dflow_factor_basin_dt_maxval
    endwhere

    ! Increment flow_factor_basin
    do j = 1, ny
       do i = 1, nx
          nb = basin_number(i,j)
          if (nb >= 1 .and. nb <= nbasin) then
             flow_factor_basin(i,j) = flow_factor_basin(i,j) + dflow_factor_basin_dt(nb) * dt
          endif
       enddo
    enddo

    ! Limit flow_factor_basin to a prescribed range
    where (flow_factor_basin > flow_factor_basin_maxval)
       flow_factor_basin =  flow_factor_basin_maxval
    elsewhere (flow_factor_basin < flow_factor_basin_minval)
       flow_factor_basin = flow_factor_basin_minval
    endwhere

    ! flow_factor_basin diagnostics for each basin

    if (verbose_inversion) then

       !Note: Some variables are 2D fields rather than basin-only fields.
       !      The logic below extracts the basin values from the 2D fields.
       !      TODO: Write a subroutine to do this?

       basin_max(:) = 0.0d0

       do j = 1, ny
          do i = 1, nx
             nb = basin_number(i,j)
             if (nb >= 1 .and. nb <= nbasin) then
                basin_max(nb) = max(basin_max(nb), flow_factor_basin(i,j))
             endif
          enddo
       enddo

       do nb = 1, nbasin
          flow_factor_basin_nb(nb) = parallel_reduce_max(basin_max(nb))
       enddo

       if (verbose_inversion .and. this_rank == rtest) then
          print*, ' '
          print*, 'basin, term1*dt, term2*dt, new flow_factor_basin:'
          do nb = 1, nbasin
             write(6,'(i6,4f12.6)') nb, &
                  dt * (floating_thck_basin(nb) - floating_thck_target_basin(nb)) / &
                  (flow_factor_basin_thck_scale * flow_factor_basin_timescale),  &
                  dt * 2.0d0 * floating_dthck_dt_basin(nb) / flow_factor_basin_thck_scale, &
                  flow_factor_basin_nb(nb)
          enddo
       endif

    endif   ! verbose_inversion

  end subroutine glissade_inversion_flow_factor_basin

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
