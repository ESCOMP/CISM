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
  use parallel

  implicit none

  ! All subroutines in this module are public

  !-----------------------------------------------------------------------------
  ! Subroutines to invert for basal fields (including basal traction beneath
  ! grounded ice and basal melting beneath floating ice) by relaxing toward
  ! a target ice thickness field.
  !-----------------------------------------------------------------------------

    logical, parameter :: verbose_inversion = .false.

!***********************************************************************

contains

!***********************************************************************

  subroutine glissade_init_inversion(model)

    ! Initialize inversion for fields of basal traction and basal melting

    use glissade_masks, only: glissade_get_masks
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
         ocean_mask,           & ! = 1 where ice is absent and topg < eus, else = 0
         land_mask               ! = 1 where topg >= eus, else = 0

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         thck_flotation          ! flotation thickness (m)

    real(dp) :: dthck

    ! Set local diagnostic point
    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif


    if (model%options%which_ho_inversion == HO_INVERSION_COMPUTE) then

       ! Save the initial ice thickness, if it will be used as the observational target for inversion.
       ! Note: If calving is done at initialization, the target is the post-calving thickness.
       !       The inversion will not try to put ice where, e.g., initial icebergs are removed.

       ! Check whether thck_obs has been read in already.
       ! If not, then set thck_obs to the initial thickness (possibly modified by initial calving).
       var_maxval = maxval(model%geometry%thck_obs)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; thck_obs has been read in already (e.g., after restart)
       else

          ! initialize to the input thickness
          model%geometry%thck_obs(:,:) = model%geometry%thck(:,:)

          !TODO - Delete the following code if not needed?
          ! If adjusting topography, then the target thck_obs is not uniquely grounded or floating.

          ! Adjust thck_obs so that the observational target is not too close to thck_flotation.
          ! The reason for this is that if we restore H to values very close to thck_flotation,
          !  it is easy for cells to flip between grounded and floating in the forward run.

!          where (model%geometry%topg - model%climate%eus < 0.0d0)
!             thck_flotation = -(rhoo/rhoi) * (model%geometry%topg - model%climate%eus) * thk0   ! convert to m
!          elsewhere
!             thck_flotation = 0.0d0
!          endwhere

!          do j = 1, model%general%nsn
!             do i = 1, model%general%ewn
!                if (model%geometry%thck_obs(i,j) > 0.0d0) then
!                   dthck = model%geometry%thck_obs(i,j)*thk0 - thck_flotation(i,j)  ! difference (m)
!                   if (abs(dthck) < inversion%bmlt_thck_buffer) then
!                      if (dthck > 0.0d0) then
!                         model%geometry%thck_obs(i,j) = (thck_flotation(i,j) + inversion%bmlt_thck_buffer) / thk0
!                      else
!                         model%geometry%thck_obs(i,j) = (thck_flotation(i,j) - inversion%bmlt_thck_buffer) / thk0
!                      endif
!                   endif
!                endif
!                model%geometry%thck_obs(i,j) = max(model%geometry%thck_obs(i,j), 0.0d0)
!             enddo
!          enddo

          ! Another adjustment: Make sure the observational target thickness is not too close to thklim.
          ! If thck_obs is close to thklim, there is a greater chance that the modeled thickness will
          !  flicker on either side of thklim.

!          where (model%geometry%thck_obs > 0.0d0)
!             model%geometry%thck_obs = max(model%geometry%thck_obs, &
!                                           model%numerics%thklim + inversion%bmlt_thck_buffer/thk0)
!          endwhere

       endif  ! var_maxval

       ! Check whether usrf_obs has been read in already.
       ! If not, then set usrf_obs to the initial upper surface elevation (possibly modified by initial calving).
       var_maxval = maxval(model%geometry%usrf_obs)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing
       else
          ! initialize to the input upper surface elevation
          model%geometry%usrf_obs(:,:) = model%geometry%usrf(:,:)
       endif

       ! Check whether topg_obs has been read in already.
       ! If not, then set topg_obs to the initial topography.
       var_maxval = maxval(model%geometry%topg_obs)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing
       else
          ! initialize to the input topography
          model%geometry%topg_obs(:,:) = model%geometry%topg(:,:)
       endif

!       call parallel_halo(model%geometry%thck_obs)
       call parallel_halo(model%geometry%usrf_obs)
       call parallel_halo(model%geometry%topg_obs)

       ! Check whether powerlaw_c_inversion has been read in already.
       ! If not, then set to a constant value.
       var_maxval = maxval(model%inversion%powerlaw_c_inversion)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; powerlaw_c_inversion has been read in already (e.g., after restart)
       else
          ! setting to a large value so that basal flow starts slow and gradually speeds up as needed
          model%inversion%powerlaw_c_inversion(:,:) = model%inversion%powerlaw_c_max
!!          model%inversion%powerlaw_c_inversion(:,:) = model%inversion%powerlaw_c
       endif

       call parallel_halo(model%inversion%powerlaw_c_inversion)

    elseif (model%options%which_ho_inversion == HO_INVERSION_PRESCRIBE) then

       ! prescribing basal friction coefficient and basal melting from previous inversion

       ! Check that the required fields from the inversion are present: powerlaw_c_inversion and bmlt_float_inversion.

       ! Note: A good way to supply powerlaw_c_prescribed is to compute powerlaw_c_inversion
       !        over some period at the end of the inversion run, after the ice is spun up.
       !       After the inversion run, rename powerlaw_c_inversion_tavg as powerlaw_c_presribed and
       !        copy it to the input file for the prescribed run.
       !       And similarly for bmlt_float_inversion and bmlt_float_prescribed

       var_maxval = maxval(model%inversion%powerlaw_c_prescribed)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! powerlaw_c_prescribed has been read in as required
          write(message,*) 'powerlaw_c_prescribed has been read from input file'
          call write_log(trim(message))
       else
          write(message,*) 'ERROR: Must read powerlaw_c_prescribed from input file to use this inversion option'
          call write_log(trim(message), GM_FATAL)
       endif

       call parallel_halo(model%inversion%powerlaw_c_prescribed)

       var_maxval = maxval(abs(model%inversion%bmlt_float_prescribed))
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! bmlt_float_prescribed has been read in as required
          write(message,*) 'bmlt_float_prescribed has been read from input file'
          call write_log(trim(message))
       else
          write(message,*) 'ERROR: Must read bmlt_float_prescribed from input file to use this inversion option'
          call write_log(trim(message), GM_FATAL)
       endif

       call parallel_halo(model%inversion%bmlt_float_prescribed)

       ! If not a restart, then initialize powerlaw_c_inversion and bmlt_float_inversion to presribed values.
       ! If a restart run, both fields typically are read from the restart file.
       !  An exception would be if we are starting an inversion run in restart mode, using a restart file
       !  from the end of a spin-up with inversion. In this case the restart file would contain the fields
       !  powerlaw_c_prescribed and bmlt_float_prescribed, and we still need to initialize powerlaw_c_inversion
       !   and bmlt_float_inversion.
 
       ! Note: powerlaw_c_inversion is adjusted at runtime where either
       !       (1) Ice is grounded in the forward run but powerlaw_c was not computed in the inversion run, or
       !       (2) Ice is floating in the forward run

       var_maxval = maxval(abs(model%inversion%powerlaw_c_inversion))
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! powerlaw_c_inversion has been read from a restart file; nothing to do here
       else
          ! initialize powerlaw_c_inversion
          model%inversion%powerlaw_c_inversion(:,:) = model%inversion%powerlaw_c_prescribed(:,:)
       endif

       call parallel_halo(model%inversion%powerlaw_c_inversion)

       var_maxval = maxval(abs(model%inversion%bmlt_float_inversion))
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! bmlt_float_inversion has been read from a restart file; nothing to do here
       else
          ! initialize bmlt_float_inversion
          model%inversion%bmlt_float_inversion(:,:) = model%inversion%bmlt_float_prescribed(:,:)
       endif

       call parallel_halo(model%inversion%bmlt_float_inversion)

   endif  ! which_ho_inversion

  end subroutine glissade_init_inversion

!***********************************************************************

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
                                   floating_mask,            &
                                   land_mask,                &
                                   grounding_line_mask,      &
                                   usrf,                     &
                                   usrf_obs,                 &
                                   dthck_dt)

    ! Compute a spatially varying basal traction field, powerlaw_c_inversion.
    ! The method is similar to that of Pollard & DeConto (TC, 2012), and is applied to all grounded ice.
    ! Where thck > thck_obs, powerlaw_c is reduced to increase sliding.
    ! Where thck < thck_obs, powerlaw_c is increased to reduce sliding.
    ! Note: powerlaw_c is constrained to lie within a prescribed range.
    ! Note: The grounding line mask is computed before horizontal transport.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_inversion), intent(inout) :: &
         inversion               ! inversion object

    integer, dimension(nx,ny), intent(in) :: &
         ice_mask,             & ! = 1 where ice is present (thk > 0), else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask,            & ! = 1 if topg > eus, else = 0
         grounding_line_mask     ! = 1 if a cell is adjacent to the grounding line, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         usrf,                 & ! ice upper surface elevation (m)
         usrf_obs,             & ! observed upper surface elevation (m)
         dthck_dt                ! rate of change of ice thickness (m/s)

    ! local variables

    integer, dimension(nx,ny) :: &
         powerlaw_c_inversion_mask  ! = 1 where we invert for powerlaw_c, else = 0

    real(dp), dimension(nx,ny) ::  &
         dusrf,                & ! usrf - usrf_obs on ice grid
         old_powerlaw_c,       & ! old value of powerlaw_c_inversion (start of timestep)
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
    ! * babc_time_smoothing   = factor for exponential moving average of usrf_inversion and dthck_dt_inversion
    !                                     (used to adjust powerlaw_c_inversion; larger => more smoothing)
    !
    ! Note on babc_space_smoothing: A smoothing factor of 1/8 gives a 4-1-1-1-1 smoother.
    !       This is numerically well behaved, but may oversmooth in bowl-shaped regions;
    !        a smaller value may be better as H converges toward H_obs.

    dpowerlaw_c(:,:) = 0.0d0

    ! Compute difference between current and target upper surface elevation
    dusrf(:,:) = usrf(:,:) - usrf_obs(:,:)

    ! Compute a mask of cells where we invert for powerlaw_c.
    ! The mask includes land-based cells, as well as marine-based cells that are grounded 
    !  or are adjacent to the grounding line.
    ! (Floating cells are GL-adjacent if they have at least one grounded neighbor.)
    ! The mask should be computed before transport, so that (for instance) if a cell is grounded
    !  during transport and floating afterward, powerlaw_c_inversion is computed here
    !  rather than being set to zero.

    where ( land_mask == 1 .or. (ice_mask == 1 .and. floating_mask == 0) .or. grounding_line_mask == 1 )
       powerlaw_c_inversion_mask = 1
    elsewhere
       powerlaw_c_inversion_mask = 0
    endwhere

    call parallel_halo(powerlaw_c_inversion_mask)

    ! Check for newly grounded cells that have powerlaw_c = 0 (from when they were ice-free or floating).
    ! Give these cells a sensible default value (either land or marine).
    do j = 1, ny
       do i = 1, nx
          if (powerlaw_c_inversion_mask(i,j) == 1) then  ! ice is land-based, grounded or GL-adjacent

             if (inversion%powerlaw_c_inversion(i,j) == 0.0d0) then
                ! set to a sensible default
                ! If on land, set to a typical land value
                ! If grounded marine ice, set to a smaller value
                if (land_mask(i,j) == 1) then
                   inversion%powerlaw_c_inversion(i,j) = inversion%powerlaw_c_land
                else
                   inversion%powerlaw_c_inversion(i,j) = inversion%powerlaw_c_marine
                endif
             endif  ! powerlaw_c_inversion = 0

          endif  ! powerlaw_c_inversion_mask = 1
       enddo  ! i
    enddo  ! j

    call parallel_halo(inversion%powerlaw_c_inversion)

    ! Loop over cells
    ! Note: powerlaw_c_inversion is computed at cell centers where usrf and thck are located.
    !       Later, it is interpolated to vertices where beta and basal velocity are located.

    do j = 1, ny
       do i = 1, nx
          if (powerlaw_c_inversion_mask(i,j) == 1) then  ! ice is land-based, grounded or GL-adjacent

             ! Save the starting value
             old_powerlaw_c(i,j) = inversion%powerlaw_c_inversion(i,j)

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
                  * inversion%powerlaw_c_inversion(i,j) * (term1 + term2)

             ! Limit to prevent huge change in one step
             if (abs(dpowerlaw_c(i,j)) > 0.05 * inversion%powerlaw_c_inversion(i,j)) then
                if (dpowerlaw_c(i,j) > 0.0d0) then
                   dpowerlaw_c(i,j) =  0.05d0 * inversion%powerlaw_c_inversion(i,j)
                else
                   dpowerlaw_c(i,j) = -0.05d0 * inversion%powerlaw_c_inversion(i,j)
                endif
             endif

             inversion%powerlaw_c_inversion(i,j) = inversion%powerlaw_c_inversion(i,j) + dpowerlaw_c(i,j)

             ! Limit to a physically reasonable range
             inversion%powerlaw_c_inversion(i,j) = min(inversion%powerlaw_c_inversion(i,j), &
                                                       inversion%powerlaw_c_max)
             inversion%powerlaw_c_inversion(i,j) = max(inversion%powerlaw_c_inversion(i,j), &
                                                       inversion%powerlaw_c_min)

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for powerlaw_c and coulomb_c: rank, i, j =', rtest, itest, jtest
                print*, 'usrf, usrf_obs, dusrf, dthck_dt:', usrf(i,j), usrf_obs(i,j), dusrf(i,j), dthck_dt(i,j)*scyr
                print*, '-dusrf/usrf_scale, -dthck_dt/dthck_dt_scale, sum =', &
                     term1, term2, & 
                     term1 + term2
                print*, 'dpowerlaw_c, newpowerlaw_c =', dpowerlaw_c(i,j), inversion%powerlaw_c_inversion(i,j)
             endif

          else  ! powerlaw_c_inversion_mask = 0

             ! set powerlaw_c = 0
             ! Note: Zero values are ignored when interpolating powerlaw_c to vertices,
             !       and in forward runs where powerlaw_c is prescribed from a previous inversion.
             ! Warning: If a cell is grounded some of the time and floating the rest of the time,
             !           the time-averaging routine will accumulate zero values as if they are real.
             !          Time-average fields should be used with caution.

             inversion%powerlaw_c_inversion(i,j) = 0.0d0

          endif  ! powerlaw_c_inversion_mask
       enddo  ! i
    enddo  ! j

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Before smoothing, powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') inversion%powerlaw_c_inversion(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif
 
    if (inversion%babc_space_smoothing > 0.0d0) then

       ! Save the value just computed
       temp_powerlaw_c(:,:) = inversion%powerlaw_c_inversion(:,:)

       ! Apply Laplacian smoothing.
       ! Since powerlaw_c lives at cell centers but is interpolated to vertices, smoothing can damp checkerboard noise.
       !TODO - Write an operator for Laplacian smoothing?
       do j = 2, ny-1
          do i = 2, nx-1
             if (powerlaw_c_inversion_mask(i,j) == 1) then  ! ice is grounded or GL-adjacent

                dpowerlaw_c_smooth = -4.0d0 * inversion%babc_space_smoothing * temp_powerlaw_c(i,j)
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if ((ii == i .or. jj == j) .and. (ii /= i .or. jj /= j)) then  ! edge neighbor
                         if (powerlaw_c_inversion_mask(ii,jj) == 1) then  ! neighbor is grounded or GL-adjacent
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

                if (dpowerlaw_c(i,j) > 0.0d0) then
                   if (temp_powerlaw_c(i,j) + dpowerlaw_c_smooth > old_powerlaw_c(i,j)) then
                      inversion%powerlaw_c_inversion(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth
                   else
                      ! allow the smoothing to hold Cp at its old value, but not reduce Cp
                      inversion%powerlaw_c_inversion(i,j) = old_powerlaw_c(i,j)
                   endif
                elseif (dpowerlaw_c(i,j) < 0.0d0) then
                   if (temp_powerlaw_c(i,j) + dpowerlaw_c_smooth < old_powerlaw_c(i,j)) then
                      inversion%powerlaw_c_inversion(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth
                   else
                      ! allow the smoothing to hold Cp at its old value, but not increase Cp
                      inversion%powerlaw_c_inversion(i,j) = old_powerlaw_c(i,j)
                   endif
                endif  ! dpowerlaw_c > 0

             endif  ! powerlaw_c_inversion_mask = 1
          enddo   ! i
       enddo   ! j

    endif  ! smoothing factor > 0

    call parallel_halo(inversion%powerlaw_c_inversion)

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
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
       if (inversion%babc_space_smoothing > 0.0d0) then
          print*, ' '
          print*, 'After smoothing, powerlaw_c:'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') inversion%powerlaw_c_inversion(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif
    endif

  end subroutine invert_basal_traction

!***********************************************************************

  subroutine prescribe_basal_traction(nx,            ny,            &
                                      itest, jtest,  rtest,         &
                                      inversion,                    &
                                      ice_mask,                     &
                                      floating_mask,                &
                                      land_mask,                    &
                                      grounding_line_mask)

    ! Compute Cp = powerlaw_c when Cp is prescribed from a previous inversion run.
    ! - For cells where the ice is grounded and a prescribed Cp exists,
    !   we simply have Cp = Cp_prescribed.
    ! - For cells where the ice is grounded and the prescribed Cp = 0 (since the cell
    !   was floating or ice-free in the inversion run), we set Cp to a sensible default
    !   based on whether the cell is land-based or marine-based.
    ! - For cells where the ice is floating (whether or not a prescribed Cp exists),
    !   we set Cp = 0.

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_inversion), intent(inout) :: &
         inversion               ! inversion object

    integer, dimension(nx,ny), intent(in) :: &
         ice_mask,             & ! = 1 where ice is present (thck > 0), else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask,            & ! = 1 if topg > eus, else = 0
         grounding_line_mask     ! = 1 if a cell is adjacent to the grounding line, else = 0

    ! local variables

    integer, dimension(nx,ny) :: &
         powerlaw_c_inversion_mask  ! = 1 where we invert for powerlaw_c, else = 0

    integer :: i, j, ii, jj

    ! Compute a mask of cells where powerlaw_c is nonzero.
    ! The mask includes cells that are grounded and/or are adjacent to the grounding line.
    ! Floating cells are GL-adjacent if they have at least one grounded neighbor.

    where ( (ice_mask == 1 .and. floating_mask == 0) .or. grounding_line_mask == 1 )
       powerlaw_c_inversion_mask = 1
    elsewhere
       powerlaw_c_inversion_mask = 0
    endwhere

    call parallel_halo(powerlaw_c_inversion_mask)

    ! Assign values of powerlaw_c

    do j = 1, ny
       do i = 1, nx
          if (powerlaw_c_inversion_mask(i,j) == 1) then

             if (inversion%powerlaw_c_prescribed(i,j) > 0.0d0) then ! use the prescribed value

                inversion%powerlaw_c_inversion(i,j) = inversion%powerlaw_c_prescribed(i,j)

             else  ! assign a sensible default

                if (land_mask(i,j) == 1) then
                   inversion%powerlaw_c_inversion(i,j) = inversion%powerlaw_c_land
                else
                   inversion%powerlaw_c_inversion(i,j) = inversion%powerlaw_c_marine
                endif

             endif  ! powerlaw_c_prescribed > 0

          endif  ! powerlaw_c_inversion_mask
       enddo  ! i
    enddo  ! j

    call parallel_halo(inversion%powerlaw_c_inversion)

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'floating_mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') floating_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'powerlaw_c_prescribed:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') inversion%powerlaw_c_prescribed(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'powerlaw_c_inversion:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') inversion%powerlaw_c_inversion(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif  ! verbose

  end subroutine prescribe_basal_traction

  !***********************************************************************

  subroutine invert_bmlt_float(dt,                           &
                               nx,            ny,            &
                               itest, jtest,  rtest,         &
                               inversion,                    &
                               thck,                         &
                               usrf_obs,                     &
                               topg,                         &
                               eus,                          &
                               ice_mask,                     &
                               floating_mask,                &
                               land_mask)

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
         topg                    ! bedrock topography (m)

    real(dp), intent(in) :: &
         eus                     ! eustatic sea level (m)

   ! Note: When this subroutine is called, ice_mask = 1 where thck > 0, not thck > thklim.
    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask               ! = 1 where topg >= eus, else = 0

    ! local variables

    integer, dimension(nx,ny) ::  &
         bmlt_inversion_mask     ! = 1 for cells where bmlt_float is computed and applied, else = 0

    real(dp), dimension(nx,ny):: &
         thck_flotation,       & ! thickness at which ice becomes afloat (m)
         thck_cavity,          & ! thickness of ocean cavity beneath floating ice (m)
         thck_target             ! thickness target (m); = thck_obs unless thck_obs > thck_flotation

    integer :: i, j, ii, jj, iglobal, jglobal

    character(len=100) :: message

    real(dp) :: bmlt_factor          ! factor for reducing basal melting

    real(dp), parameter :: inversion_bmlt_timescale = 0.0d0*scyr  ! timescale for freezing in cavities (m/s)

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

    ! For floating and weakly grounded cells, compute a target thickness based on the target surface elevation.
    ! Note: Cells with a floating target are always restored to that target.
    !       Cells with a grounded target are restored to a thickness slightly greater than thck_flotation,
    !        provided they are currently floating or weakly grounded.
    !       These cells are not restored all the way to the grounded target.

    ! initialize
    bmlt_inversion_mask(:,:) = 0
    thck_target(:,:) = 0.0d0
    inversion%bmlt_float_inversion(:,:) = 0.0d0

    ! loop over cells
    do j = 1, ny
       do i = 1, nx

          if (land_mask(i,j) == 1) then

             ! do nothing; bmlt_float_inversion = 0

          elseif (usrf_obs(i,j) - (topg(i,j) - eus) > thck_flotation(i,j)) then  ! grounded target

             ! If the ice is now floating or very weakly grounded (thck < thck_flotation + thck_buffer),
             !  then compute bmlt_float < 0 to restore to thck_flotation + thck_buffer.

             if (thck(i,j) < thck_flotation(i,j) + inversion%bmlt_thck_buffer) then ! floating or very weakly grounded

                ! Restore to a thickness slightly greater than thck_flotation
                ! (generally not all the way to the observed thickness, since we would prefer
                !  for the basal traction inversion to achieve this by roughening the bed).
                bmlt_inversion_mask(i,j) = 1
                thck_target(i,j) = thck_flotation(i,j) + inversion%bmlt_thck_buffer
                
             else   ! strongly grounded

                ! do nothing; bmlt_float_inversion = 0

             endif

          elseif (usrf_obs(i,j) > 0.0d0) then  ! floating target

             ! Note: With usrf_obs > 0 requirement, we do not melt columns that are ice-free in observations.
             !       This allows the calving front to advance (if not using a no-advance calving mask).

             bmlt_inversion_mask(i,j) = 1
             thck_target(i,j) = usrf_obs(i,j) * rhoo/(rhoo - rhoi)

          endif

          if (bmlt_inversion_mask(i,j) == 1) then
             inversion%bmlt_float_inversion(i,j) = (thck(i,j) - thck_target(i,j)) / dt

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for bmlt_float_inversion: rank, i, j =', rtest, itest, jtest
                print*, 'topg - eus, usrf_obs:', topg(i,j) - eus, usrf_obs(i,j)
                print*, 'thck, thck_target, bmlt_float:', & 
                     thck(i,j), thck_target(i,j), inversion%bmlt_float_inversion(i,j)*dt
             endif

          endif   ! bmlt_inversion_mask = 1

       enddo   ! i
    enddo   ! j

    !TODO - Test the following code further, or delete it?  So far, I haven't found a timescale to work well.
    ! If a nonzero timescale is specified, then multiply bmlt_float_inversion by a factor
    !  proportional to dt/timescale, in the range (0,1].

    if (inversion_bmlt_timescale > 0.0d0) then
       bmlt_factor = min(dt/inversion_bmlt_timescale, 1.0d0)
       inversion%bmlt_float_inversion(:,:) = inversion%bmlt_float_inversion(:,:) * bmlt_factor
    endif

    call parallel_halo(bmlt_inversion_mask) ! diagnostic only
    call parallel_halo(thck_target)  ! diagnostic only

    call parallel_halo(inversion%bmlt_float_inversion)

    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'bmlt_inversion_mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') bmlt_inversion_mask(i,j)
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
             write(6,'(f10.3)',advance='no') thck_flotation(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck_cavity (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_cavity(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck_target (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_target(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bmlt_float_inversion (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') inversion%bmlt_float_inversion(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine invert_bmlt_float

!***********************************************************************

  subroutine prescribe_bmlt_float(dt,                           &
                                  nx,            ny,            &
                                  itest, jtest,  rtest,         &
                                  inversion,                    &
                                  thck,                         &
                                  topg,                         &
                                  eus,                          &
                                  ice_mask,                     &
                                  floating_mask,                &
                                  land_mask,                    &
                                  grounding_line_mask)

    ! Prescribe bmlt_float based on the value computed from inversion.
    ! Note: bmlt_float_inversion is defined as positive for melting, negative for freezing.
    ! This field is applied only in floating and grounding-line adjacent cells.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_inversion), intent(inout) :: &
         inversion               ! inversion object

    ! Note: thck should be the expected values after applying the mass balance
    !       (although the mass balance may not yet have been applied) 
    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                 & ! ice thickness (m)
         topg                    ! bedrock elevation (m)

    real(dp), intent(in) :: &
         eus                     ! eustatic sea level (m)

   ! Note: When this subroutine is called, ice_mask = 1 where thck > 0, not thck > thklim.
    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask,            & ! = 1 where topg >= eus, else = 0
         grounding_line_mask     ! = 1 if a cell is adjacent to the grounding line, else = 0

    ! local variables

    integer, dimension(nx,ny) ::  &
         bmlt_inversion_mask     ! = 1 for cells where bmlt_float is computed and applied, else = 0

    real(dp), dimension(nx,ny):: &
         thck_flotation,       & ! flotation thickness (m)
         thck_cavity,          & ! thickness (m) of ocean cavity (diagnostic only)
         thck_final              ! final thickness (m) if full melt rate is applied

    integer :: i, j

    ! Note: This subroutine should be called after other mass-balance terms have been applied,
    !  after horizontal transport, and preferably after calving.

    if (verbose_inversion .and. main_task) then
       print*, ' '
       print*, 'In prescribe_bmlt_float'
    endif

    ! Compute the flotation thickness
    where (topg - eus < 0.0d0)
       thck_flotation = -(rhoo/rhoi) * (topg - eus)
    elsewhere
       thck_flotation = 0.0d0
    endwhere

    ! Compute the ocean cavity thickness beneath floating ice (diagnostic only)
    where (floating_mask == 1)
       thck_cavity = -(topg - eus) - (rhoi/rhoo)*thck
    elsewhere
       thck_cavity = 0.0d0
    endwhere

    ! Compute a mask of floating cells and marine-based grounding-line cells, 
    !  where bmlt_float_inversion can potentially be applied.
    ! Where bmlt_inversion_mask = 1, apply bmlt_float_prescribed.
    ! Note: The land mask may not be needed, since land cells are excluded from the inversion.
    !       But this mask is included for generality, in case of dynamic topography.

    bmlt_inversion_mask(:,:) = 0
    inversion%bmlt_float_inversion(:,:) = 0.0d0
    thck_final(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx
          if (land_mask(i,j) == 1) then

             ! do nothing; bmlt_float_inversion = 0

          elseif (floating_mask(i,j) == 1 .or.  &
                  (ice_mask(i,j) == 1 .and. grounding_line_mask(i,j) == 1) ) then

             bmlt_inversion_mask(i,j) = 1
             inversion%bmlt_float_inversion(i,j) = inversion%bmlt_float_prescribed(i,j)

             ! Make sure the final thickness is non-negative.

             thck_final(i,j) = thck(i,j) - inversion%bmlt_float_inversion(i,j)*dt

             if (thck_final(i,j) < 0.0d0) then
                thck_final(i,j) = 0.0d0
                inversion%bmlt_float_inversion(i,j) = (thck(i,j) - thck_final(i,j)) / dt
             endif

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Prescribe bmlt_float_inversion: rank, i, j =', rtest, itest, jtest
                print*, 'thck, thck_final, bmlt_float_inversion*dt:', thck(i,j), thck_final(i,j), &
                     inversion%bmlt_float_inversion(i,j)*dt
             endif

          endif   ! masks
       enddo   ! i
    enddo   ! j

    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'bmlt_float_prescribed (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') inversion%bmlt_float_prescribed(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bmlt_inversion_mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') bmlt_inversion_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'grounding_line_mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') grounding_line_mask(i,j)
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
             write(6,'(f10.3)',advance='no') thck_flotation(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck_cavity (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_cavity(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thck_final:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_final(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bmlt_float_inversion (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') inversion%bmlt_float_inversion(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine prescribe_bmlt_float

!=======================================================================

end module glissade_inversion

!=======================================================================
