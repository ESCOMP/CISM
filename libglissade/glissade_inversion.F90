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

  !-----------------------------------------------------------------------------
  ! Subroutines to invert for basal fields (including basal traction beneath
  ! grounded ice and basal melting beneath floating ice) by relaxing toward
  ! a target ice thickness field.
  !-----------------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************

  subroutine glissade_init_inversion(model)

    use glissade_masks, only: glissade_get_masks

    ! Initialize inversion for fields of basal traction and basal melting

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    real(dp) :: var_maxval          ! max value of a given variable; = 0 if not yet read in

    character(len=100) :: message

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
          model%geometry%thck_obs(:,:) = model%geometry%thck(:,:)
       endif

       call parallel_halo(model%geometry%thck_obs)

       ! Check whether powerlaw_c_2d has been read in already.
       ! If not, then set to a constant value.
       var_maxval = maxval(model%basal_physics%powerlaw_c_2d)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; powerlaw_c_2d has been read in already (e.g., after restart)
       else
          model%basal_physics%powerlaw_c_2d(:,:) = model%basal_physics%powerlaw_c
       endif

       call parallel_halo(model%basal_physics%powerlaw_c_2d)

    elseif (model%options%which_ho_inversion == HO_INVERSION_PRESCRIBED) then

       ! prescribing basal friction coefficient and basal melting from previous inversion

       ! Check that the required fields from the inversion are present: powerlaw_c_2d and bmlt_float_inversion.

       ! Note: A good way to supply powerlaw_c_2d is to compute powerlaw_c_2d_tavg
       !        over some period at the end of the inversion run, after the ice is spun up.
       !       To output this field from the inversion run, uncomment 'average: 1' under
       !         powerlaw_c_2d in glide_vars.def, then configure and rebuild the code.
       !       After the inversion run, rename powerlaw_c_2d_tavg as powerlaw_c_2d and
       !        copy it to the input file for the prescribed run.
       !       And similarly for bmlt_float_inversion_tavg and bmlt_float_inversion

       var_maxval = maxval(model%basal_physics%powerlaw_c_2d)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! powerlaw_c_2d has been read in as required
          write(message,*) 'powerlaw_c_2d has been read from input file'
          call write_log(trim(message))
       else
          write(message,*) 'ERROR: Must read powerlaw_c_2d from input file to use this inversion option'
          call write_log(trim(message), GM_FATAL)
       endif

       call parallel_halo(model%basal_physics%powerlaw_c_2d)

       var_maxval = maxval(abs(model%basal_melt%bmlt_float_inversion))
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! bmlt_float_inversion has been read in as required
          write(message,*) 'bmlt_float_inversion has been read from input file'
          call write_log(trim(message))
       else
          write(message,*) 'ERROR: Must read bmlt_float_inversion from input file to use this inversion option'
          call write_log(trim(message), GM_FATAL)
       endif

       call parallel_halo(model%basal_melt%bmlt_float_inversion)

    endif  ! which_ho_inversion

  end subroutine glissade_init_inversion

!***********************************************************************

  !TODO - Add code to set powerlaw_c for prescribed case.
  !       Use prescribed values where available, and otherwise extrapolate from nearby values.
  !       In this way, we can avoid having very wrong values where the GL has advanced.

!***********************************************************************

  subroutine invert_basal_traction(dt,                           &
                                   nx,            ny,            &
                                   itest, jtest,  rtest,         &
                                   basal_physics,                &
                                   ice_mask,      floating_mask, &
                                   thck,          dthck_dt,      &
                                   thck_obs)

    ! Compute spatially varying fields, powerlaw_c_2d and coulomb_c_2d, by inversion.
    ! The method is similar to that of Pollard & DeConto (TC, 2012), and is applied to all grounded ice.
    ! Where thck > thck_obs, powerlaw_c and coulomb_c are reduced to increase sliding.
    ! Where thck < thck_obs, powerlaw_c and coulomb_c are increased to reduce sliding.
    ! Note: powerlaw_c is constrained to lie within a prescribed range.
    !       The ratio of powerlaw_c to coulomb_c is fixed (except that coulomb_c must be <= 1).

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_basal_physics), intent(inout) :: &
         basal_physics           ! basal physics object

    integer, dimension(nx,ny), intent(in) :: &
         ice_mask,             & ! = 1 where ice is present (thk > thklim), else = 0
         floating_mask           ! = 1 where ice is present and floating, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                 & ! ice thickness (m)
         dthck_dt,             & ! rate of change of ice thickness (m/s)
         thck_obs                ! observed thickness (m)

    ! local variables

    real(dp), dimension(nx,ny) ::  &
         dthck,                & ! thck - thck_obs on ice grid
         old_powerlaw_c,       & ! old value of powerlaw_c_2d (start of timestep)
         temp_powerlaw_c,      & ! temporary value of powerlaw_c_2d (before smoothing)
         dpowerlaw_c             ! change in powerlaw_c

    real(dp) :: term1, term2
    real(dp) :: factor
    real(dp) :: dpowerlaw_c_smooth

    integer :: i, j
    integer :: ii, jj

    ! inversion parameters in basal_physics derived type:
    ! * powerlaw_c_max                  = upper bound for powerlaw_c, Pa (m/yr)^(-1/3)
    ! * powerlaw_c_min                  = lower bound for powerlaw_c, Pa (m/yr)^(-1/3)
    ! * powerlaw_coulomb_ratio          = powerlaw_c/coulomb_c (same units as powerlaw_c)
    ! * inversion_babc_timescale        = inversion timescale (s); must be > 0
    ! * inversion_babc_thck_scale       = thickness inversion scale (m); must be > 0
    ! * inversion_babc_dthck_dt_scale   = dthck_dt inversion scale (m/s); must be > 0
    ! * inversion_babc_smoothing_factor = factor for smoothing powerlaw_c_2d; higher => more smoothing
    !
    ! Note on smoothing: A smoothing factor of 1/8 gives a 4-1-1-1-1 smoother.
    !       This is numerically well behaved, but may oversmooth in bowl-shaped regions;
    !        a smaller value may be better as H converges toward H_obs.

    logical, parameter :: verbose_inversion = .false.

    ! Save the starting value
    old_powerlaw_c(:,:) = basal_physics%powerlaw_c_2d(:,:)
    dpowerlaw_c(:,:) = 0.0d0

    ! Compute difference between current and target thickness
    dthck(:,:) = thck(:,:) - thck_obs(:,:)

    ! Loop over cells
    ! Note: powerlaw_c_2d and coulomb_c_2d are computed at cell centers where thck is located.
    !       Later, they are interpolated to vertices where beta and basal velocity are located.

    do j = 1, ny
       do i = 1, nx

          if (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) then  ! ice is present and grounded

             ! Invert for powerlaw_c_2d and coulomb_c based on dthck and dthck_dt
             term1 = -dthck(i,j) / basal_physics%inversion_babc_thck_scale
             term2 = -dthck_dt(i,j) / basal_physics%inversion_babc_dthck_dt_scale

             dpowerlaw_c(i,j) = (dt/basal_physics%inversion_babc_timescale) &
                  * basal_physics%powerlaw_c_2d(i,j) * (term1 + term2)

             ! Limit to prevent huge change in one step
             if (abs(dpowerlaw_c(i,j)) > 0.05 * basal_physics%powerlaw_c_2d(i,j)) then
                if (dpowerlaw_c(i,j) > 0.0d0) then
                   dpowerlaw_c(i,j) =  0.05d0 * basal_physics%powerlaw_c_2d(i,j)
                else
                   dpowerlaw_c(i,j) = -0.05d0 * basal_physics%powerlaw_c_2d(i,j)
                endif
             endif

             basal_physics%powerlaw_c_2d(i,j) = basal_physics%powerlaw_c_2d(i,j) + dpowerlaw_c(i,j)

             ! Limit to a physically reasonable range
             basal_physics%powerlaw_c_2d(i,j) = min(basal_physics%powerlaw_c_2d(i,j), basal_physics%powerlaw_c_max)
             basal_physics%powerlaw_c_2d(i,j) = max(basal_physics%powerlaw_c_2d(i,j), basal_physics%powerlaw_c_min)

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for powerlaw_c and coulomb_c: rank, i, j =', rtest, itest, jtest
                print*, 'thck, thck_obs, dthck, dthck_dt:', thck(i,j), thck_obs(i,j), dthck(i,j), dthck_dt(i,j)*scyr
                print*, '-dthck/thck_scale, -dthck_dt/dthck_dt_scale, sum =', &
                     -dthck(i,j)/basal_physics%inversion_babc_thck_scale, &
                     -dthck_dt(i,j)/basal_physics%inversion_babc_dthck_dt_scale, &
                     term1 + term2
                print*, 'dpowerlaw_c, newpowerlaw_c =', dpowerlaw_c(i,j), basal_physics%powerlaw_c_2d(i,j)
             endif

          else  ! ice_mask = 0 or floating_mask = 1

             ! set to default value
             basal_physics%powerlaw_c_2d(i,j) = basal_physics%powerlaw_c

          endif  ! ice_mask = 1 and floating_mask = 0

       enddo  ! i
    enddo  ! j

    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Before smoothing, powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') basal_physics%powerlaw_c_2d(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif
 
    ! Save the value just computed
    temp_powerlaw_c(:,:) = basal_physics%powerlaw_c_2d(:,:)

    ! Apply Laplacian smoothing to C_p.
    ! Since C_p is at cell centers but is interpolated to vertices, smoothing can damp checkerboard noise.
    !TODO - Write an operator for Laplacian smoothing?
    do j = 2, ny-1
       do i = 2, nx-1
          if (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) then  ! cell (i,j) is grounded

             dpowerlaw_c_smooth = -4.0d0 * basal_physics%inversion_babc_smoothing_factor * temp_powerlaw_c(i,j)
             do jj = j-1, j+1
                do ii = i-1, i+1
                   if ((ii == i .or. jj == j) .and. (ii /= i .or. jj /= j)) then  ! edge neighbor
                      if (ice_mask(ii,jj) == 1 .and. floating_mask(ii,jj) == 0) then   ! cell (ii,jj) is grounded
                         dpowerlaw_c_smooth = dpowerlaw_c_smooth &
                              + basal_physics%inversion_babc_smoothing_factor*temp_powerlaw_c(ii,jj)
                      else
                         dpowerlaw_c_smooth = dpowerlaw_c_smooth &
                              + basal_physics%inversion_babc_smoothing_factor*temp_powerlaw_c(i,j)
                      endif
                   endif
                enddo
             enddo

             ! Note: If smoothing is too strong, it can reverse the sign of the change in powerlaw_c.
             !       The logic below ensures that if powerlaw_c_2d is increasing, the smoothing can reduce
             !        the change to zero, but not cause powerlaw_c to decrease relative to old_powerlaw_c
             !        (and similarly if powerlaw_c_2d is decreasing).

             if (dpowerlaw_c(i,j) > 0.0d0) then
                if (temp_powerlaw_c(i,j) + dpowerlaw_c_smooth > old_powerlaw_c(i,j)) then
                   basal_physics%powerlaw_c_2d(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth
                else
                  ! allow the smoothing to hold Cp at its old value, but not reduce Cp
                   basal_physics%powerlaw_c_2d(i,j) = old_powerlaw_c(i,j)
                endif
             elseif (dpowerlaw_c(i,j) < 0.0d0) then
                if (temp_powerlaw_c(i,j) + dpowerlaw_c_smooth < old_powerlaw_c(i,j)) then
                   basal_physics%powerlaw_c_2d(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth
                else
                  ! allow the smoothing to hold Cp at its old value, but not increase Cp
                   basal_physics%powerlaw_c_2d(i,j) = old_powerlaw_c(i,j)
                endif
             endif  ! dpowerlaw_c > 0

             ! The next 5 lines are commented out. If used in place of the limiting above,
             !  this code not only prevents the sign of the change from reversing, but also
             !  prevents the smoothing from more than doubling the original change.
             ! It would take more testing to determine whether or not this is a good idea.

!             if (abs(dpowerlaw_c_smooth) > abs(dpowerlaw_c(i,j))) then
!                factor = abs(dpowerlaw_c(i,j)) / abs(dpowerlaw_c_smooth)
!                dpowerlaw_c_smooth = dpowerlaw_c_smooth * factor
!             endif
!             basal_physics%powerlaw_c_2d(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth

          endif  ! cell is grounded

          if (verbose_inversion .and. this_rank==rtest .and. i==itest .and. j==jtest) then
             print*, 'Smoothing correction, new powerlaw_c:', dpowerlaw_c_smooth, basal_physics%powerlaw_c_2d(i,j)
          endif

       enddo
    enddo

    call parallel_halo(basal_physics%powerlaw_c_2d)

    ! Set coulomb_c assuming a fixed ratio of powerlaw_c/coulomb_c
    basal_physics%coulomb_c_2d(:,:) = basal_physics%powerlaw_c_2d(:,:) / basal_physics%powerlaw_coulomb_ratio

    ! Limit coulomb_c to be <= 1, so that basal stress <= effective pressure N
    basal_physics%coulomb_c_2d(:,:) = min(basal_physics%coulomb_c_2d(:,:), 1.0d0)

    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then

       i = itest
       j = jtest
       print*, 'thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, 'thck - thck_obs:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') dthck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, 'dthck_dt (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') dthck_dt(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'After smoothing, powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') basal_physics%powerlaw_c_2d(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, 'coulomb_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.4)',advance='no') basal_physics%coulomb_c_2d(i,j)
          enddo
          write(6,*) ' '
       enddo

    endif

  end subroutine invert_basal_traction

  !***********************************************************************

  subroutine invert_bmlt_float(dt,                           &
                               nx,            ny,            &
                               itest, jtest,  rtest,         &
                               basal_melt,                   &
                               thck,                         &
                               thck_obs,                     &
                               ice_mask,                     &
                               floating_mask,                &
                               ocean_mask,                   &
                               land_mask)

    ! Compute spatially varying bmlt_float by inversion.
    ! Where thck > thck_obs, bmlt_float_inversion is increased.
    ! Where thck < thck_obs, bmlt_float_inversion is decreased.
    ! Note: bmlt_float_inversion is defined as positive for melting, negative for freezing.

    !TODO - Move this subroutine?
    use glissade_calving, only: glissade_find_lakes

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_basal_melt), intent(inout) :: &
         basal_melt              ! basal melt object

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                 & ! ice thickness (m)
         thck_obs                ! observed thickness (m)

    ! Note: When this subroutine is called, ice_mask = 1 where thck > 0, not thck > thklim.
    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         ocean_mask,           & ! = 1 where ice is absent and topg < eus, else = 0
         land_mask               ! = 1 where topg >= eus, else = 0

    ! local variables

    integer, dimension(nx,ny) ::  &
         lake_mask               ! = 1 for floating cells disconnected from the ocean

    real(dp), dimension(nx,ny) ::  &
         dthck,                & ! thck - thck_obs on ice grid
         old_bmlt_float,       & ! old value of bmlt_float_inversion (start of timestep)
         temp_bmlt_float,      & ! temporary value of bmlt_float_inversion (before smoothing)
         dbmlt_float             ! change in bmlt_float_inversion

    real(dp) :: term1, dbmlt_float_smooth

    integer :: i, j, ii, jj

    ! Where the observed ice is floating, adjust the basal melt rate (or freezing rate, if bmlt < 0)
    !  so as to relax the ice thickness toward the observed target.
    ! Note: This subroutine should be called after other mass-balance terms have been applied,
    !  and after horizontal transport.
    ! We compute the difference (H - bmlt_float_inversion*dt) - H_obs,
    !  which is the thickness error that would remain after applying the current bmlt_float_inversion.
    ! We then increase or decrease bmlt_float_inversion with a characteristic timescale,
    !  thereby reducing the thickness error.
    ! As the timescale approaches zero, the adjusted bmlt_float_inversion will approach the value
    !  needed to give H = H_obs.

    logical, parameter :: verbose_inversion = .false.

    if (verbose_inversion .and. main_task) then
       print*, ' '
       print*, 'In invert_bmlt_float'
    endif

    ! Identify lake cells: floating interior cells that will not be restored
    ! to the target thickness

    call glissade_find_lakes(nx,               ny,                &
                             itest,   jtest,   rtest,             &
                             ice_mask,         floating_mask,     &
                             ocean_mask,       lake_mask)

    ! Compute a mask of cells where bmlt_float_inversion will be computed
    !TODO - Make bmlt_inversion_mask a local field?

    basal_melt%bmlt_inversion_mask(:,:) = 0

    do j = 2, ny-1
       do i = 2, nx-1
          if (floating_mask(i,j) == 1) then
             ! check for land neighbors
             if (land_mask(i-1,j) == 1 .or. land_mask(i+1,j) == 1 .or. &
                 land_mask(i,j-1) == 1 .or. land_mask(i,j+1) == 1) then
                ! mask = 0; do not invert for bmlt_float
             elseif (lake_mask(i,j) == 1) then
                ! mask = 0; do not invert for bmlt_float
             else
                basal_melt%bmlt_inversion_mask(i,j) = 1
             endif
          endif
       enddo
    enddo

    call parallel_halo(basal_melt%bmlt_inversion_mask)

    ! Save the starting value of bmlt_float_inversion
    old_bmlt_float(:,:) = basal_melt%bmlt_float_inversion(:,:)
    dbmlt_float(:,:) = 0.0d0

    ! Compute difference between the current and target thickness
    dthck(:,:) = thck(:,:) - thck_obs(:,:)

    ! Loop over cells
    do j = 1, ny
       do i = 1, nx

          if (basal_melt%bmlt_inversion_mask(i,j) == 1) then

             if (basal_melt%inversion_bmlt_timescale > 0.0d0) then
                ! Adjust bmlt_float_inversion to reduce the thickness error
                dbmlt_float(i,j) = (dthck(i,j) - basal_melt%bmlt_float_inversion(i,j)*dt)  &
                                  / basal_melt%inversion_bmlt_timescale 
             else
                ! Set bmlt_float_inversion such that thck = thck_obs after inversion
                dbmlt_float(i,j) = dthck(i,j)/dt - basal_melt%bmlt_float_inversion(i,j)
             endif

             basal_melt%bmlt_float_inversion(i,j) = basal_melt%bmlt_float_inversion(i,j) + dbmlt_float(i,j)

             !WHL - I think this may not be needed
             ! Limit to a physically reasonable range
!             basal_melt%bmlt_float_inversion(i,j) = min(basal_melt%bmlt_float_inversion(i,j), basal_melt%bmlt_float_inversion_max)
!             basal_melt%bmlt_float_inversion(i,j) = max(basal_melt%bmlt_float_inversion(i,j), basal_melt%bmlt_float_inversion_min)

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for bmlt_float_inversion: rank, i, j =', rtest, itest, jtest
                print*, 'thck, thck_obs, dthck, bmlt*dt:', &
                     thck(i,j), thck_obs(i,j), dthck(i,j), basal_melt%bmlt_float_inversion(i,j)*dt
                print*, 'dbmlt_float, new bmlt_float (m/yr) =', dbmlt_float(i,j)*scyr, basal_melt%bmlt_float_inversion(i,j)*scyr
             endif

          else  ! bmlt_inversion_mask = 0

             basal_melt%bmlt_float_inversion(i,j) = 0.0d0

          endif  ! bmlt_inversion_mask

       enddo  ! i
    enddo  ! j

    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Before smoothing, bmlt_float (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') basal_melt%bmlt_float_inversion(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Save the value just computed
    temp_bmlt_float(:,:) = basal_melt%bmlt_float_inversion(:,:)

    ! Apply Laplacian smoothing to bmlt_float_inversion.
    !TODO - Write an operator for Laplacian smoothing?
    do j = 2, ny-1
       do i = 2, nx-1
          if (basal_melt%bmlt_inversion_mask(i,j) == 1) then

             dbmlt_float_smooth = -4.0d0 * basal_melt%inversion_bmlt_smoothing_factor * temp_bmlt_float(i,j)
             do jj = j-1, j+1
                do ii = i-1, i+1
                   if ((ii == i .or. jj == j) .and. (ii /= i .or. jj /= j)) then  ! edge neighbor
                      if (basal_melt%bmlt_inversion_mask(ii,jj) == 1) then   ! inverting for bmlt_float in cell (ii,jj)
                         dbmlt_float_smooth = dbmlt_float_smooth &
                              + basal_melt%inversion_bmlt_smoothing_factor*temp_bmlt_float(ii,jj)
                      else
                         dbmlt_float_smooth = dbmlt_float_smooth &
                              + basal_melt%inversion_bmlt_smoothing_factor*temp_bmlt_float(i,j)
                      endif
                   endif
                enddo
             enddo

             ! Note: If smoothing is too strong, it can reverse the sign of the change in bmlt_float.
             !       The logic below ensures that if bmlt_float is increasing, the smoothing can reduce
             !        the change to zero, but not cause bmlt_float to decrease relative to old_bmlt_float
             !        (and similarly if bmlt_float is decreasing).

             if (dbmlt_float(i,j) > 0.0d0) then
                if (temp_bmlt_float(i,j) + dbmlt_float_smooth > old_bmlt_float(i,j)) then
                   basal_melt%bmlt_float_inversion(i,j) = temp_bmlt_float(i,j) + dbmlt_float_smooth
                else
                  ! allow the smoothing to hold bmlt_float at its old value, but not reduce bmlt_float
                   basal_melt%bmlt_float_inversion(i,j) = old_bmlt_float(i,j)
                endif
             elseif (dbmlt_float(i,j) < 0.0d0) then
                if (temp_bmlt_float(i,j) + dbmlt_float_smooth < old_bmlt_float(i,j)) then
                   basal_melt%bmlt_float_inversion(i,j) = temp_bmlt_float(i,j) + dbmlt_float_smooth
                else
                   ! allow the smoothing to hold bmlt_float at its old value, but not increase bmlt_float
                   basal_melt%bmlt_float_inversion(i,j) = old_bmlt_float(i,j)
                endif
             endif  ! dbmlt_float > 0

          endif  ! bmlt_inversion_mask = 1

          if (verbose_inversion .and. this_rank==rtest .and. i==itest .and. j==jtest) then
             print*, 'Smoothing correction, new bmlt_float:', dbmlt_float_smooth*scyr, basal_melt%bmlt_float_inversion(i,j)*scyr
          endif

       enddo
    enddo

    call parallel_halo(basal_melt%bmlt_float_inversion)

    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, 'thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, 'thck - bmlt*dt - thck_obs:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') dthck(i,j) - basal_melt%bmlt_float_inversion(i,j)*dt
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'After smoothing, bmlt_float_inversion (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') basal_melt%bmlt_float_inversion(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif


  end subroutine invert_bmlt_float
  
!=======================================================================

end module glissade_inversion

!=======================================================================
