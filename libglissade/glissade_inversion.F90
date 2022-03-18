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
  public :: verbose_inversion, glissade_init_inversion, &
            glissade_inversion_basal_friction,  &
            glissade_inversion_bmlt_basin

  !-----------------------------------------------------------------------------
  ! Subroutines to invert for basal fields (including basal friction beneath
  ! grounded ice and basal melting beneath floating ice) by relaxing toward
  ! a target ice thickness field.
  !-----------------------------------------------------------------------------

    logical, parameter :: verbose_inversion = .false.

!***********************************************************************

contains

!***********************************************************************

  subroutine glissade_init_inversion(model)

    ! Initialize inversion for fields of basal friction and basal melting
    ! Should be called after usrf and thck have been input and (possibly) modified by initial calving

    use glissade_masks, only: glissade_get_masks
    use glissade_utils, only: glissade_usrf_to_thck, glissade_thck_to_usrf

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
       call glissade_usrf_to_thck(&
            model%geometry%usrf_obs,  &
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
          call glissade_thck_to_usrf(&
               thck_obs,  &
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
          ! initialize to a uniform value (which can be set in the config file)
          model%basal_physics%coulomb_c(:,:) = model%basal_physics%coulomb_c_const
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
    ! computations specific to inversion of deltaT_basin
    !----------------------------------------------------------------------

    if (model%options%which_ho_bmlt_basin_inversion == HO_BMLT_BASIN_INVERSION_COMPUTE) then

       if (model%options%is_restart == RESTART_FALSE) then

          ! Set floating_thck_target for floating ice and lightly grounded ice.
          ! Here, "lightly grounded" means that the magnitude of f_flotation = (-topg - eus) - (rhoi/rhoo)*thck
          !  is less than a prescribed threshold.  (Recall f_flotation < 0 for grounded ice.)
          ! The inversion will nudge the ice thickness toward this target in a basin-average sense.
          ! Positive volume biases will be corrected with ocean warming, and negative biases with cooling.

          do j = 1, nsn
             do i = 1, ewn
                f_flotation = (-(model%geometry%topg(i,j) - model%climate%eus)  &
                              - (rhoi/rhoo)*model%geometry%thck(i,j)) * thk0    ! f_flotation < 0 for grounded ice
                if (model%geometry%thck(i,j) > 0.0d0 .and. &
                    model%geometry%marine_connection_mask(i,j) == 1 .and. &
                    f_flotation > -model%inversion%bmlt_basin_flotation_threshold) then
                   model%inversion%floating_thck_target(i,j) = model%geometry%thck(i,j)
                else
                   model%inversion%floating_thck_target(i,j) = 0.0d0
                endif
             enddo
          enddo

          if (verbose_inversion .and. this_rank == rtest) then
             print*, 'bmlt_basin_flotation_threshold =', model%inversion%bmlt_basin_flotation_threshold
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

  subroutine glissade_inversion_basal_friction(model)

    use glimmer_paramets, only: tim0, thk0, vel0
    use glimmer_physcon, only: scyr
    use glissade_grid_operators, only: glissade_stagger, glissade_stagger_real_mask
    use glissade_utils, only: glissade_usrf_to_thck

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
       call glissade_usrf_to_thck(&
            model%geometry%usrf_obs,  &
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
          endif   ! verbose_inversion

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
                                     model%basal_physics%coulomb_c)

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
                                   friction_c)

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

    ! local variables

    real(dp), dimension(nx-1,ny-1) ::  &
         stag_dthck,           & ! stag_thck - stag_thck_obs
         dvelo_sfc,            & ! velo_sfc - velo_sfc_obs
         dfriction_c             ! change in friction_c

    real(dp) :: term1_thck, term2_thck  ! tendency terms based on thickness target
    real(dp) :: term1_velo              ! tendency term based on surface speed target

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

    ! Loop over vertices where f_ground > 0
    ! Note: f_ground should be computed before transport, so that if a vertex is grounded
    !       before transport and fully floating afterward, friction_c is computed here.

    do j = 1, ny-1
       do i = 1, nx-1

          if (f_ground(i,j) > 0.0d0) then  ! ice is at least partly grounded

             ! Compute the rate of change of friction_c, based on stag_dthck and stag_dthck_dt,
             !  and/or dvelo_sfc.
             ! For a thickness target, the rate of change is proportional to the sum of two terms:
             !     dC/dt = -C * (1/tau) * (H - H_obs)/H0 + (2*tau/H0) * dH/dt
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
                term1_thck = -stag_dthck(i,j) / (babc_thck_scale * babc_timescale)
                term2_thck = -stag_dthck_dt(i,j) * 2.0d0 / babc_thck_scale
             else
                term1_thck = 0.0d0
                term2_thck = 0.0d0
             endif

             ! Compute tendency terms based on the surface speed target
             ! Note: I tried adding a term2_velo in analogy to term2_thck (Dec. 2021),
             !       but it triggers oscillations in friction_c without improving accuracy.

             if (babc_velo_scale > 0.0d0) then
                term1_velo = dvelo_sfc(i,j) / (babc_velo_scale * babc_timescale)
             else
                term1_velo = 0.0d0
             endif

             dfriction_c(i,j) = friction_c(i,j) * (term1_thck + term2_thck + term1_velo) * dt

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
                     term1_thck*dt, term2_thck*dt, (term1_thck + term2_thck)*dt
                print*, 'dvelo term =', term1_velo*dt
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
                                           dbmlt_dtemp_scale,           &
                                           bmlt_basin_timescale,        &
                                           deltaT_basin,                &
                                           bmlt_basin_mass_correction,  &
                                           bmlt_basin_number_mass_correction)

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

    real(dp), intent(in), optional :: &
         bmlt_basin_mass_correction   ! optional mass correction (Gt) for a selected basin

    integer, intent(in), optional :: &
         bmlt_basin_number_mass_correction ! integer ID for the basin receiving the correction

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

    ! Optionally, apply a correction to the ice volume target in a selected basin.
    ! Note: This option could in principle be applied to multiple basins, but currently is supported for one basin only.
    !       In practice, this basin is likely to be the Amundsen Sea Embayment (ISMIP6 basin #9).

    if (present(bmlt_basin_mass_correction) .and. present(bmlt_basin_number_mass_correction)) then
       if (abs(bmlt_basin_mass_correction) > 0.0d0 .and. bmlt_basin_number_mass_correction > 0) then

          nb = bmlt_basin_number_mass_correction
          floating_volume_target_basin(nb) = floating_volume_target_basin(nb) + &
               bmlt_basin_mass_correction * (1.0d12/rhoi)   ! Gt converted to m^3
          if (verbose_inversion .and. this_rank == rtest) then
             print*, ' '
             print*, 'Basin with mass correction:', bmlt_basin_number_mass_correction
             print*, 'mass correction (Gt) =', bmlt_basin_mass_correction
             print*, 'volume correction (km^3) =', bmlt_basin_mass_correction * (1.0d3/rhoi)
             print*, 'New volume target (km^3) =', floating_volume_target_basin(nb) / 1.0d9
          endif
       endif
    endif   ! present(bmlt_basin_mass_correction)

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

    if (verbose_inversion .and. this_rank == rtest) then
       if (present(bmlt_basin_mass_correction) .and. present(bmlt_basin_number_mass_correction)) then
          if (abs(bmlt_basin_mass_correction) > 0.0d0 .and. bmlt_basin_number_mass_correction > 0) then
             nb = bmlt_basin_number_mass_correction
             print*, 'New basin thickness target =', floating_thck_target_basin(nb)
          endif
       endif
    endif

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

       if (verbose_inversion .and. this_rank == rtest) then
          print*, 'bmlt_basin_timescale (yr) =', bmlt_basin_timescale/scyr
          print*, 'dbmlt_dtemp_scale (m/yr/degC) =', dbmlt_dtemp_scale
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
          print*, ' '
          print*, 'basin, term1*dt, term2*dt, dTbasin, new deltaT_basin:'
          do nb = 1, nbasin
             write(6,'(i6,4f12.6)') nb, &
                  dt/dbmlt_dtemp_scale * (floating_thck_basin(nb) - floating_thck_target_basin(nb)) / &
                  (bmlt_basin_timescale**2), &
                  dt/dbmlt_dtemp_scale * 2.0d0 * floating_dthck_dt_basin(nb) / bmlt_basin_timescale, &
                  dt*dTbasin_dt(nb), deltaT_basin_nb(nb)
          enddo
       endif

    endif   ! verbose_inversion

  end subroutine glissade_inversion_bmlt_basin

!=======================================================================

end module glissade_inversion

!=======================================================================
