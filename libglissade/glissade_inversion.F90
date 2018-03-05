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

!!    logical, parameter :: verbose_inversion = .false.
    logical, parameter :: verbose_inversion = .true.

    !TODO - Make these config parameters?
    real(dp), parameter :: &
         powerlaw_c_land = 20000.d0, &
         powerlaw_c_marine = 1000.d0

    real(dp), parameter :: &
         bmlt_inversion_thck_over_flot = 1.0d0   ! Ice restored from floating to grounded is this much thicker than thck_flotation (m)

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
          model%geometry%thck_obs(:,:) = model%geometry%thck(:,:)
       endif

       call parallel_halo(model%geometry%thck_obs)

       ! Check whether powerlaw_c_inversion has been read in already.
       ! If not, then set to a constant value.
       var_maxval = maxval(model%basal_physics%powerlaw_c_inversion)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! do nothing; powerlaw_c_inversion has been read in already (e.g., after restart)
       else
          ! setting to a large value so that basal flow starts slow and gradually speeds up as needed
          model%basal_physics%powerlaw_c_inversion(:,:) = model%basal_physics%powerlaw_c_max
       endif

       call parallel_halo(model%basal_physics%powerlaw_c_inversion)

    elseif (model%options%which_ho_inversion == HO_INVERSION_PRESCRIBED) then

       ! prescribing basal friction coefficient and basal melting from previous inversion

       ! Check that the required fields from the inversion are present: powerlaw_c_inversion and bmlt_float_inversion.

       ! Note: A good way to supply powerlaw_c_prescribed is to compute powerlaw_c_inversion
       !        over some period at the end of the inversion run, after the ice is spun up.
       !       After the inversion run, rename powerlaw_c_inversion_tavg as powerlaw_c_presribed and
       !        copy it to the input file for the prescribed run.
       !       And similarly for bmlt_float_inversion and bmlt_float_prescribed

       var_maxval = maxval(model%basal_physics%powerlaw_c_prescribed)
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! powerlaw_c_prescribed has been read in as required
          write(message,*) 'powerlaw_c_prescribed has been read from input file'
          call write_log(trim(message))
       else
          write(message,*) 'ERROR: Must read powerlaw_c_prescribed from input file to use this inversion option'
          call write_log(trim(message), GM_FATAL)
       endif

       call parallel_halo(model%basal_physics%powerlaw_c_prescribed)

       var_maxval = maxval(abs(model%basal_melt%bmlt_float_prescribed))
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! bmlt_float_prescribed has been read in as required
          write(message,*) 'bmlt_float_prescribed has been read from input file'
          call write_log(trim(message))
       else
          write(message,*) 'ERROR: Must read bmlt_float_prescribed from input file to use this inversion option'
          call write_log(trim(message), GM_FATAL)
       endif

       call parallel_halo(model%basal_melt%bmlt_float_prescribed)

       ! If not a restart, then initialize powerlaw_c_inversion and bmlt_float_inversion to presribed values.
       ! If a restart run, both fields typically are read from the restart file.
       !  An exception would be if we are starting an inversion run in restart mode, using a restart file
       !  from the end of a spin-up with inversion. In this case the restart file would contain the fields
       !  powerlaw_c_prescribed and bmlt_float_prescribed, and we still need to initialize powerlaw_c_inversion
       !   and bmlt_float_inversion.
 
       ! Note: powerlaw_c_inversion is adjusted at runtime where either
       !       (1) Ice is grounded in the forward run but powerlaw_c was not computed in the inversion run, or
       !       (2) Ice is floating in the forward run

       var_maxval = maxval(abs(model%basal_physics%powerlaw_c_inversion))
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! powerlaw_c_inversion has been read from a restart file; nothing to do here
       else
          ! initialize powerlaw_c_inversion
          model%basal_physics%powerlaw_c_inversion(:,:) = model%basal_physics%powerlaw_c_prescribed(:,:)
       endif

       call parallel_halo(model%basal_physics%powerlaw_c_inversion)

       var_maxval = maxval(abs(model%basal_melt%bmlt_float_inversion))
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval > 0.0d0) then
          ! bmlt_float_inversion has been read from a restart file; nothing to do here
       else
          ! initialize bmlt_float_inversion
          model%basal_melt%bmlt_float_inversion(:,:) = model%basal_melt%bmlt_float_prescribed(:,:)
       endif

       call parallel_halo(model%basal_melt%bmlt_float_inversion)

   endif  ! which_ho_inversion

  end subroutine glissade_init_inversion

!***********************************************************************

  subroutine invert_basal_traction(dt,                       &
                                   nx,            ny,        &
                                   itest, jtest,  rtest,     &
                                   basal_physics,            &
                                   ice_mask,                 &
                                   floating_mask,            &
                                   land_mask,                &
                                   thck,                     &
                                   dthck_dt,                 &
                                   thck_obs)

    ! Compute a spatially varying basal traction field, powerlaw_c_inversion.
    ! The method is similar to that of Pollard & DeConto (TC, 2012), and is applied to all grounded ice.
    ! Where thck > thck_obs, powerlaw_c is reduced to increase sliding.
    ! Where thck < thck_obs, powerlaw_c is increased to reduce sliding.
    ! Note: powerlaw_c is constrained to lie within a prescribed range.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_basal_physics), intent(inout) :: &
         basal_physics           ! basal physics object

    integer, dimension(nx,ny), intent(in) :: &
         ice_mask,             & ! = 1 where ice is present (thk > 0), else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask               ! = 1 if topg > eus, else = 0

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                 & ! ice thickness (m)
         dthck_dt,             & ! rate of change of ice thickness (m/s)
         thck_obs                ! observed thickness (m)

    ! local variables

    real(dp), dimension(nx,ny) ::  &
         dthck,                & ! thck - thck_obs on ice grid
         old_powerlaw_c,       & ! old value of powerlaw_c_inversion (start of timestep)
         temp_powerlaw_c,      & ! temporary value of powerlaw_c_inversion (before smoothing)
         dpowerlaw_c             ! change in powerlaw_c

    real(dp) :: term1, term2
    real(dp) :: factor
    real(dp) :: dpowerlaw_c_smooth
    real(dp) :: sum_powerlaw_c

    integer :: i, j, ii, jj
    integer :: count

    ! inversion parameters in basal_physics derived type:
    ! * powerlaw_c_max                  = upper bound for powerlaw_c, Pa (m/yr)^(-1/3)
    ! * powerlaw_c_min                  = lower bound for powerlaw_c, Pa (m/yr)^(-1/3)
    ! * inversion_babc_timescale        = inversion timescale (s); must be > 0
    ! * inversion_babc_thck_scale       = thickness inversion scale (m); must be > 0
    ! * inversion_babc_dthck_dt_scale   = dthck_dt inversion scale (m/s); must be > 0
    ! * inversion_babc_smoothing_factor = factor for smoothing powerlaw_c_inversion; higher => more smoothing
    !
    ! Note on smoothing: A smoothing factor of 1/8 gives a 4-1-1-1-1 smoother.
    !       This is numerically well behaved, but may oversmooth in bowl-shaped regions;
    !        a smaller value may be better as H converges toward H_obs.

    dpowerlaw_c(:,:) = 0.0d0

    ! Compute difference between current and target thickness
    dthck(:,:) = thck(:,:) - thck_obs(:,:)

    ! Check for newly grounded cells that have powerlaw_c = 0 (from when they were ice-free or floating).
    ! Give these cells a sensible default value by extrapolating from neighbor cells.
    do j = 2, ny-1
       do i = 2, nx-1
          if (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) then  ! grounded ice
             if (basal_physics%powerlaw_c_inversion(i,j) == 0.0d0) then
                ! set to a sensible default
                ! If on land, set to a typical land value
                ! If grounded marine ice, set to a smaller value
                if (land_mask(i,j) == 1) then
                   basal_physics%powerlaw_c_inversion(i,j) = powerlaw_c_land
                else
                   basal_physics%powerlaw_c_inversion(i,j) = powerlaw_c_marine
                endif
             endif  ! powerlaw_c_inversion = 0
          endif  ! grounded ice
       enddo  ! i
    enddo  ! j

    call parallel_halo(basal_physics%powerlaw_c_inversion)

    ! Loop over cells
    ! Note: powerlaw_c_inversion is computed at cell centers where thck is located.
    !       Later, it is interpolated to vertices where beta and basal velocity are located.

    do j = 1, ny
       do i = 1, nx

          if (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) then  ! ice is present and grounded

             ! Save the starting value
             old_powerlaw_c(i,j) = basal_physics%powerlaw_c_inversion(i,j)

             ! Invert for powerlaw_c based on dthck and dthck_dt
             term1 = -dthck(i,j) / basal_physics%inversion_babc_thck_scale
             term2 = -dthck_dt(i,j) / basal_physics%inversion_babc_dthck_dt_scale

             dpowerlaw_c(i,j) = (dt/basal_physics%inversion_babc_timescale) &
                  * basal_physics%powerlaw_c_inversion(i,j) * (term1 + term2)

             ! Limit to prevent huge change in one step
             if (abs(dpowerlaw_c(i,j)) > 0.05 * basal_physics%powerlaw_c_inversion(i,j)) then
                if (dpowerlaw_c(i,j) > 0.0d0) then
                   dpowerlaw_c(i,j) =  0.05d0 * basal_physics%powerlaw_c_inversion(i,j)
                else
                   dpowerlaw_c(i,j) = -0.05d0 * basal_physics%powerlaw_c_inversion(i,j)
                endif
             endif

             basal_physics%powerlaw_c_inversion(i,j) = basal_physics%powerlaw_c_inversion(i,j) + dpowerlaw_c(i,j)

             ! Limit to a physically reasonable range
             basal_physics%powerlaw_c_inversion(i,j) = min(basal_physics%powerlaw_c_inversion(i,j), &
                                                           basal_physics%powerlaw_c_max)
             basal_physics%powerlaw_c_inversion(i,j) = max(basal_physics%powerlaw_c_inversion(i,j), &
                                                           basal_physics%powerlaw_c_min)

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for powerlaw_c and coulomb_c: rank, i, j =', rtest, itest, jtest
                print*, 'thck, thck_obs, dthck, dthck_dt:', thck(i,j), thck_obs(i,j), dthck(i,j), dthck_dt(i,j)*scyr
                print*, '-dthck/thck_scale, -dthck_dt/dthck_dt_scale, sum =', &
                     -dthck(i,j)/basal_physics%inversion_babc_thck_scale, &
                     -dthck_dt(i,j)/basal_physics%inversion_babc_dthck_dt_scale, &
                     term1 + term2
                print*, 'dpowerlaw_c, newpowerlaw_c =', dpowerlaw_c(i,j), basal_physics%powerlaw_c_inversion(i,j)
             endif

          else  ! ice_mask = 0 or floating_mask = 1

             ! set powerlaw_c = 0
             ! Note: Zero values are ignored when interpolating powerlaw_c to vertices,
             !       and in forward runs where powerlaw_c is prescribed from a previous inversion.
             ! Warning: If a cell is grounded some of the time and floating the rest of the time,
             !           the time-averaging routine will accumulate zero values as if they are real.
             !          Time-average fields should be used with caution.

             basal_physics%powerlaw_c_inversion(i,j) = 0.0d0

          endif  ! grounded ice

       enddo  ! i
    enddo  ! j

    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Before smoothing, powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') basal_physics%powerlaw_c_inversion(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif
 
    ! Save the value just computed
    temp_powerlaw_c(:,:) = basal_physics%powerlaw_c_inversion(:,:)

    ! Apply Laplacian smoothing.
    ! Since powerlaw_c lives at cell centers but is interpolated to vertices, smoothing can damp checkerboard noise.
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
             !       The logic below ensures that if powerlaw_c is increasing, the smoothing can reduce
             !        the change to zero, but not cause powerlaw_c to decrease relative to old_powerlaw_c
             !        (and similarly if powerlaw_c is decreasing).

             if (dpowerlaw_c(i,j) > 0.0d0) then
                if (temp_powerlaw_c(i,j) + dpowerlaw_c_smooth > old_powerlaw_c(i,j)) then
                   basal_physics%powerlaw_c_inversion(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth
                else
                  ! allow the smoothing to hold Cp at its old value, but not reduce Cp
                   basal_physics%powerlaw_c_inversion(i,j) = old_powerlaw_c(i,j)
                endif
             elseif (dpowerlaw_c(i,j) < 0.0d0) then
                if (temp_powerlaw_c(i,j) + dpowerlaw_c_smooth < old_powerlaw_c(i,j)) then
                   basal_physics%powerlaw_c_inversion(i,j) = temp_powerlaw_c(i,j) + dpowerlaw_c_smooth
                else
                  ! allow the smoothing to hold Cp at its old value, but not increase Cp
                   basal_physics%powerlaw_c_inversion(i,j) = old_powerlaw_c(i,j)
                endif
             endif  ! dpowerlaw_c > 0

          endif  ! cell is grounded
       enddo
    enddo

    call parallel_halo(basal_physics%powerlaw_c_inversion)

    ! Set coulomb_c to a constant
    !TODO - Switch from array to constant field in basal traction subroutine
    basal_physics%coulomb_c_inversion(:,:) = basal_physics%coulomb_c

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
             write(6,'(f10.2)',advance='no') basal_physics%powerlaw_c_inversion(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine invert_basal_traction

  !***********************************************************************

  subroutine prescribe_basal_traction(nx,            ny,            &
                                      itest, jtest,  rtest,         &
                                      basal_physics,                &
                                      ice_mask,                     &
                                      floating_mask,                &
                                      land_mask)

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

    type(glide_basal_physics), intent(inout) :: &
         basal_physics           ! basal physics object

    integer, dimension(nx,ny), intent(in) :: &
         ice_mask,             & ! = 1 where ice is present (thck > 0), else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask               ! = 1 where topg >= eus, else = 0

    ! local variables

    real(dp), dimension(nx,ny) :: &
         new_powerlaw_c          ! new powerlaw_c values extrapolated from existing values

    integer :: i, j, ii, jj

    integer :: count              ! counter
    real(dp) :: sum_powerlaw_c    ! sum of powerlaw_c in neighbor cells

    ! Zero out powerlaw_c where ice is not present (thck > 0) and grounded
    where (ice_mask == 0 .or. floating_mask == 1)
       basal_physics%powerlaw_c_inversion = 0.0d0
    endwhere

    ! Assign values of powerlaw_c in newly grounded cells

    do j = 2, ny-1
       do i = 2, nx-1
          if (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) then  ! grounded ice

             if (basal_physics%powerlaw_c_inversion(i,j) > 0.0d0) then

                ! nothing to do here; cell was already grounded

             elseif (basal_physics%powerlaw_c_prescribed(i,j) > 0.0d0) then ! use the prescribed value

                basal_physics%powerlaw_c_inversion(i,j) = basal_physics%powerlaw_c_prescribed(i,j)

             else  ! assign a sensible default

                if (land_mask(i,j) == 1) then
                   basal_physics%powerlaw_c_inversion(i,j) = powerlaw_c_land
                else
                   basal_physics%powerlaw_c_inversion(i,j) = powerlaw_c_marine
                endif

             endif  ! powerlaw_c > 0

          endif  ! grounded
       enddo  ! i
    enddo  ! j

    call parallel_halo(basal_physics%powerlaw_c_inversion)

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
             write(6,'(f10.2)',advance='no') basal_physics%powerlaw_c_prescribed(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'current powerlaw_c:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.2)',advance='no') basal_physics%powerlaw_c_inversion(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif  ! verbose

  end subroutine prescribe_basal_traction

!***********************************************************************

  subroutine invert_bmlt_float(dt,                           &
                               nx,            ny,            &
                               itest, jtest,  rtest,         &
                               basal_melt,                   &
                               thck,                         &
                               thck_obs,                     &
                               topg,                         &
                               acab,                         &
                               bmlt,                         &
                               ice_mask,                     &
                               floating_mask,                &
                               land_mask)

    ! Compute spatially varying bmlt_float by inversion.
    ! Apply a melt/freezing rate that will restore the ice in floating grid cells to the target thickness.
    ! Cells that are floating and should be grounded are thickened enough to be lightly grounded,
    !  but generally not all the way to the target thickness.
    ! Note: bmlt_float_inversion is defined as positive for melting, negative for freezing.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_basal_melt), intent(inout) :: &
         basal_melt              ! basal melt object

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                 & ! ice thickness (m)
         thck_obs,             & ! observed thickness (m)
         topg,                 & ! bedrock topography (m)
         acab,                 & ! surface mass balance (m/s), including runtime adjustments
         bmlt

   ! Note: When this subroutine is called, ice_mask = 1 where thck > 0, not thck > thklim.
    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask               ! = 1 where topg >= eus, else = 0

    ! local variables

    integer, dimension(nx,ny) ::  &
         bmlt_inversion_mask,  & ! = 1 for cells where bmlt_float is computed and applied, else = 0
         floating_mask_bmlt      ! = 1 for cells that are durably floating, else = 0
                                 ! (where "durably floating" is defined as floating both before and after transport)

    real(dp), dimension(nx,ny):: &
         thck_flotation,       & ! thickness at which ice becomes afloat (m)
         thck_cavity,          & ! thickness of ocean cavity beneath floating ice (m)
         thck_target             ! thickness target (m); = thck_obs unless thck_obs > thck_flotation

    integer :: i, j, ii, jj, iglobal, jglobal

    character(len=100) :: message

    real(dp) :: bmlt_factor          ! factor for reducing basal melting

    real(dp), parameter :: inversion_bmlt_timescale = 0.0d0*scyr  ! timescale for freezing in cavities (m/s)

    ! Where the observed ice is floating, adjust the basal melt rate (or freezing rate, if bmlt < 0)
    !  so as to relax the ice thickness toward a target thickness based on observations.
    ! Note: This subroutine should be called after other mass-balance terms have been applied,
    !  after horizontal transport, and preferably after calving.

    if (verbose_inversion .and. main_task) then
       print*, ' '
       print*, 'In invert_bmlt_float'
    endif

    ! Compute the flotation thickness
    where (topg < 0.0d0)
       thck_flotation = -(rhoo/rhoi)*topg
    elsewhere
       thck_flotation = 0.0d0
    endwhere

    ! Compute the ocean cavity thickness beneath floating ice
    where (floating_mask == 1)
       thck_cavity = -topg - (rhoi/rhoo)*thck
    elsewhere
       thck_cavity = 0.0d0
    endwhere

    ! Restore selected cells to a target thickness based on observations.
    ! The rules are:
    ! (1) Any ice-filled cell that is floating in observations is restored after transport
    !     to its target thickness.
    !     This includes cells that are either floating or grounded after transport.
    ! (2) Any ice-filled cell that is grounded in observations but floating (or very lightly grounded)
    !      after transport is restored to a target grounded thickness.  This grounded thickness, however,
    !      is not the observed thickness, but rather Hf + thck_over_flot.
    !     Here, thck_over_float is a small thickness (1 m by default), and "very lightly grounded"
    !      means that Hf < H < Hf + thck_over_flot.
    !     The reason for not restoring grounded cells all the way to the observed value is that
    !      we would like the basal sliding coefficients to adjust to help ground the ice more firmly,
    !      rather than rely on large negative basal melt rates.

    ! Compute a mask of floating cells where bmlt_float_inversion will be computed.
    ! The mask excludes land-based cells.
    ! It includes cells that are
    ! (1) floating in observations (even if they are grounded in the model, after transport), or
    ! (2) grounded in observations but floating in the model (after transport)

    ! initialize
    bmlt_inversion_mask(:,:) = 0
    thck_target(:,:) = 0.0d0

    ! loop over cells
    do j = 2, ny-1
       do i = 2, nx-1
          if (land_mask(i,j) == 1) then

             ! do nothing; bmlt_float_inversion = 0

          !TODO - Relax thck(i,j) > 0 requirement?
          elseif (thck(i,j) > 0.0d0 .and. thck_obs(i,j) > 0.0d0) then  ! ice-covered marine cell in obs

             if (thck_obs(i,j) < thck_flotation(i,j)) then

                ! floating in obs; restore to thck_obs

                bmlt_inversion_mask(i,j) = 1
                thck_target(i,j) = thck_obs(i,j)

                !TODO - Remove the error check when satisfied that things are working.
                if (thck(i,j) >= thck_flotation(i,j) .and. basal_melt%grounded_mask_start(i,j) == 1) then
                   ! This is not supposed to happen, so throw a fatal error.
                   call parallel_globalindex(i, j, iglobal, jglobal)

                   print*, 'Error, floating cell has grounded:, task, i, j, H_obs, H_f, H:',  &
                        this_rank, i, j, thck_obs(i,j), thck_flotation(i,j), thck(i,j)
                   print*, 'iglobal, jglobal =', iglobal, jglobal
                   write(message,*) 'ERROR in invert_bmlt_float: Cell that should be floating has grounded'
                   call write_log(message, GM_FATAL)
                endif

             elseif (thck_obs(i,j) >= thck_flotation(i,j) .and. &
                     thck(i,j) < thck_flotation(i,j) + bmlt_inversion_thck_over_flot) then

                ! grounded in obs but currently floating; reground but do not set thck = thck_obs.
                ! The reason for this is that we would prefer to adjust basal sliding parameters to bring
                !  the thickness closer to observations, instead of relying on an artificially large bmlt_float.
                bmlt_inversion_mask(i,j) = 1
                thck_target(i,j) = thck_flotation(i,j) + bmlt_inversion_thck_over_flot

             endif

          endif
       enddo   ! i
    enddo   ! j

    call parallel_halo(bmlt_inversion_mask)
    call parallel_halo(thck_target)


    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'thck_target:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_target(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
    endif

    ! Now compute bmlt_float_inversion based on the thickness target.
    ! Account for the surface mass balance, recalling that acab and bmlt have opposite sign conventions.
    ! (acab > 0 for positive SMB, whereas bmlt > 0 for negative BMB.)
    ! Typically, the background bmlt = 0 for floating cells when computing bmlt_float_inversion,
    !  but we might be computing bmlt_float_inversion for cells that were grounded at the start
    !  of the time step and thus have nonzero bmlt.

    basal_melt%bmlt_float_inversion(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx
          if (bmlt_inversion_mask(i,j) == 1) then

             basal_melt%bmlt_float_inversion(i,j) = (thck(i,j) - thck_target(i,j))/dt + acab(i,j) - bmlt(i,j)

             ! Adjust to account for the surface mass balance.

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Invert for bmlt_float_inversion: rank, i, j =', rtest, itest, jtest
                print*, 'thck, thck_obs, acab*dt, bmlt*dt:', &
                     thck(i,j), thck_obs(i,j), acab(i,j)*dt, &
                     basal_melt%bmlt_float_inversion(i,j)*dt
             endif

          endif  ! bmlt_inversion_mask = 1

       enddo   ! i
    enddo   ! j

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
    endif

    !TODO - Test this code further, or delete it?  So far, I haven't found a timescale to work well.
    ! If a nonzero timescale is specified, then multiply bmlt_float_inversion by a factor
    !  proportional to dt/timescale, in the range (0,1].

    if (inversion_bmlt_timescale > 0.0d0) then
       bmlt_factor = min(dt/inversion_bmlt_timescale, 1.0d0)
       basal_melt%bmlt_float_inversion(:,:) = basal_melt%bmlt_float_inversion(:,:) * bmlt_factor
    endif

    call parallel_halo(basal_melt%bmlt_float_inversion)

    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'thck_flotation (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_flotation(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Before bmlt inversion, floating_mask:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') floating_mask(i,j)
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
       print*, 'thck_obs (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_obs(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bmlt_float_inversion (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') basal_melt%bmlt_float_inversion(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine invert_bmlt_float

!***********************************************************************

  subroutine prescribe_bmlt_float(dt,                           &
                                  nx,            ny,            &
                                  itest, jtest,  rtest,         &
                                  basal_melt,                   &
                                  thck,                         &
                                  topg,                         &
                                  ice_mask,                     &
                                  floating_mask,                &
                                  land_mask)

    ! Prescribe bmlt_float based on the value computed from inversion.
    ! Note: bmlt_float_inversion is defined as positive for melting, negative for freezing.
    ! This field is applied only beneath floating ice.

    real(dp), intent(in) ::  dt  ! time step (s)

    integer, intent(in) :: &
         nx, ny                  ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(glide_basal_melt), intent(inout) :: &
         basal_melt              ! basal melt object

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                 & ! ice thickness (m)
         topg                    ! bedrock elevation (m)

   ! Note: When this subroutine is called, ice_mask = 1 where thck > 0, not thck > thklim.
    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,             & ! = 1 where ice is present, else = 0
         floating_mask,        & ! = 1 where ice is present and floating, else = 0
         land_mask               ! = 1 where topg >= eus, else = 0

    ! local variables

    integer, dimension(nx,ny) ::  &
         bmlt_inversion_mask     ! = 1 for cells where bmlt_float is computed and applied, else = 0

    real(dp), dimension(nx,ny):: &
         thck_flotation,       & ! flotation thickness (m)
         thck_final              ! final thickness (m) if full melt rate is applied

    integer :: i, j, ii, jj

    real(dp) :: dthck                ! thickness change (m)

    ! Note: This subroutine should be called after other mass-balance terms have been applied,
    !  after horizontal transport, and preferably after calving.

    if (verbose_inversion .and. main_task) then
       print*, ' '
       print*, 'In prescribe_bmlt_float'
    endif

    ! Make a mask to identify cells that were floating at the start of the time step (before transport)
    !  or are floating now (after transport). These are cells to which

    ! Compute the flotation thickness
    where (topg < 0.0d0)
       thck_flotation = -(rhoo/rhoi)*topg
    elsewhere
       thck_flotation = 0.0d0
    endwhere

    ! Compute a mask of floating cells where bmlt_float_inversion can potentially be computed.
    ! The rule is that bmlt_float_inversion can be applied to any grid cell that is afloat either
    !  before or after transport, with the exception of ice-free cells (with thck = 0) and land-based cells.
    ! Note: The land mask is probably not needed, since land cells are excluded from the inversion.
    !       But this mask is included for generality, in case of dynamic topography.

    bmlt_inversion_mask(:,:) = 0.0d0
    thck_final(:,:) = 0.0d0

    do j = 2, ny-1
       do i = 2, nx-1
          if (land_mask(i,j) == 1) then

             ! do nothing; bmlt_float_inversion = 0

          elseif ( ice_mask(i,j) == 1 .and. &
                   (basal_melt%floating_mask_start(i,j) == 1 .or. floating_mask(i,j) == 1) ) then

             bmlt_inversion_mask(i,j) = 1
             dthck = -basal_melt%bmlt_float_prescribed(i,j) * dt
             thck_final(i,j) = min(thck(i,j) + dthck, thck_flotation(i,j) + bmlt_inversion_thck_over_flot)

          endif
       enddo
    enddo

    call parallel_halo(bmlt_inversion_mask)
    call parallel_halo(thck_final)

    ! Now compute bmlt_float_inversion based on the final thickness.

    basal_melt%bmlt_float_inversion(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx
          if (bmlt_inversion_mask(i,j) == 1) then

             basal_melt%bmlt_float_inversion(i,j) = (thck(i,j) - thck_final(i,j)) / dt

             !WHL - debug
             if (verbose_inversion .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Prescribe bmlt_float_inversion: rank, i, j =', rtest, itest, jtest
                print*, 'thck, mask, bmlt*dt:', thck(i,j), bmlt_inversion_mask(i,j), &
                     basal_melt%bmlt_float_inversion(i,j)*dt
             endif

          endif  ! bmlt_inversion_mask = 1

       enddo   ! i
    enddo   ! j

    !WHL - debug
    if (verbose_inversion .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'thck_final:'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_final(i,j)
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
       print*, 'prescribed bmlt_float (m/yr):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') basal_melt%bmlt_float_inversion(i,j)*scyr
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine prescribe_bmlt_float

!=======================================================================

end module glissade_inversion

!=======================================================================
