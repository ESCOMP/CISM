!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   glissade_utils.F90 - part of the Community Ice Sheet Model (CISM)
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
! This module holds some utility subroutines for the Glissade dynamical core
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glissade_utils

  use glimmer_global, only: dp
  use glimmer_paramets, only: iulog, eps11, eps08
  use glimmer_log
  use glide_types
  use cism_parallel, only: this_rank, main_task

  implicit none

  private
  public :: glissade_adjust_thickness, glissade_smooth_usrf, &
       glissade_smooth_topography, glissade_adjust_topography, &
       glissade_basin_sum, glissade_basin_average, &
       glissade_usrf_to_thck, glissade_thck_to_usrf, &
       glissade_edge_fluxes, glissade_input_fluxes, &
       glissade_quadrant_sum, glissade_bounding_box, &
       glissade_rms_error, write_array_to_file, &
       glissade_cleanup_tiny_thickness, glissade_cleanup_icefree_cells

  interface write_array_to_file
     module procedure write_array_to_file_real8_2d
     module procedure write_array_to_file_real8_3d
  end interface

contains

!****************************************************************************

  subroutine glissade_adjust_thickness(model)

    ! Optionally, check for spurious surface depressions that could arise in the following case:
    ! (1) usrf, thck, and topg have all been read in.  (Recall that usrf is an optional input.)
    ! (2) There are interior lakes: regions disconnected from the ocean, where (usrf - thck) > topg.
    ! (3) The ice in these interior lake regions is too thick to float.
    ! In this case, the default behavior is to reset usrf = topg + thck, possibly leading to
    !  steep surface depressions and unstable flow.
    ! The alternative is to set thck = usrf - topg in grounded regions, maintaining the observed usrf.
    !
    !TODO: In this and the next two subroutines, we could pass in thck, topg, etc. instead of the model derived type.

    use cism_parallel, only: parallel_reduce_max

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    ! local variables

    real(dp) :: usrf_max
    real(dp) :: topg
    real(dp) :: thck_flot

    integer :: i, j
    integer :: nx, ny
    integer :: itest, jtest, rtest
    ! The following variables give the boundaries of a box going from itest-3 to itest+3,
    ! and jtest-3 to jtest+3, but limited to stay within the range of the local points
    ! owned by rdiag_local
    integer :: itest_m3, itest_p3, jtest_m3, jtest_p3

    logical, parameter :: verbose_adjust_thickness = .true.

    ! Copy some model variables to local variables

    nx = model%general%ewn
    ny = model%general%nsn

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
       itest_m3 = max(itest-3, 1)
       itest_p3 = min(itest+3, nx)
       jtest_m3 = max(jtest-3, 1)
       jtest_p3 = min(jtest+3, ny)
    endif

    ! Make sure ursf was read in with nonzero values.
    ! Otherwise, we cannot use usrf to adjust the ice thickness.

    usrf_max = maxval(model%geometry%usrf)
    usrf_max = parallel_reduce_max(usrf_max)

    if (usrf_max > tiny(0.0d0)) then

       do j = 1, ny
          do i = 1, nx
             topg = model%geometry%topg(i,j) - model%climate%eus  ! shorthand for relative bed topography
             if (model%geometry%usrf(i,j) - model%geometry%thck(i,j) > topg) then
                thck_flot = -(rhoo/rhoi) * topg
                if (model%geometry%thck(i,j) >= thck_flot) then  ! grounded
                   ! increase thck to remove the sub-ice cavity
                   model%geometry%thck(i,j) = model%geometry%usrf(i,j) - topg
                else   ! floating
                   ! do nothing; keep the existing thickness
                endif
             elseif (model%geometry%usrf(i,j) - model%geometry%thck(i,j) < topg) then
                ! reduce thck so that lsrf = topg
                model%geometry%thck(i,j) = model%geometry%usrf(i,j) - topg
             endif
          enddo
       enddo

    else   ! usrf_max < tiny

       call write_log('Error: Must read in usrf to use adjust_input_thickness option', GM_FATAL)

    endif   ! usrf_max > tiny

  end subroutine glissade_adjust_thickness

!****************************************************************************

  subroutine glissade_smooth_usrf(model, nsmooth)

    ! Use a Laplacian smoother to smooth the upper surface elevation,
    !  and compute a thickness consistent with this new elevation.
    ! This can be useful if the input thickness and topography are inconsistent,
    !  such that their sum has large gradients.

    use glide_thck, only: glide_calclsrf
    use glissade_masks, only: glissade_get_masks
    use glissade_grid_operators, only: glissade_laplacian_smoother
    use cism_parallel, only: parallel_halo

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    integer, intent(in), optional :: nsmooth     ! number of smoothing passes

    ! local variables

    real(dp), dimension(model%general%ewn, model%general%nsn) :: &
         topg,               &  ! bed topography (m)
         thck,               &  ! thickness (m)
         usrf,               &  ! surface elevation (m)
         usrf_smoothed          ! surface elevation after smoothing

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ice_mask,           &  ! = 1 if ice is present (thck > 0, else = 0
         floating_mask,      &  ! = 1 if ice is present (thck > 0) and floating, else = 0
         ocean_mask             ! = 1 if topg < 0 and ice is absent, else = 0

    integer :: n_smoothing_passes   ! local version of nsmooth
    integer :: i, j, n
    integer :: nx, ny
    integer :: itest, jtest, rtest

!    logical, parameter :: verbose_smooth_usrf = .false.
    logical, parameter :: verbose_smooth_usrf = .true.

    ! Initialize

    if (present(nsmooth)) then
       n_smoothing_passes = nsmooth
    else
       n_smoothing_passes = 1
    endif

    ! Copy some model variables to local variables

    nx = model%general%ewn
    ny = model%general%nsn

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ! compute the initial upper surface elevation
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    ! Save input fields
    topg = (model%geometry%topg - model%climate%eus)
    thck = model%geometry%thck
    usrf = model%geometry%usrf

    ! compute initial masks
    call glissade_get_masks(nx,                  ny,                    &
                            model%parallel,                             &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   0.0d0,                 &  ! thklim = 0
                            ice_mask,                                   &
                            floating_mask = floating_mask,              &
                            ocean_mask = ocean_mask)

    do n = 1, n_smoothing_passes

       call glissade_laplacian_smoother(nx,     ny,              &
                                        usrf,   usrf_smoothed,   &
                                        npoints_stencil = 9)

       ! Force usrf = topg on ice-free land
       where (topg > 0.0d0 .and. ice_mask == 0) usrf_smoothed = topg

       ! Force usrf = unsmoothed value for floating ice and ice-free ocean, to avoid advancing the calving front
       where (floating_mask == 1 .or. ocean_mask == 1)
          usrf_smoothed = usrf
       endwhere

       ! Force usrf >= topg
       usrf_smoothed = max(usrf_smoothed, topg)

       usrf = usrf_smoothed
       call parallel_halo(usrf, model%parallel)

    enddo

    ! Given the smoothed usrf, adjust the input thickness such that topg is unchanged.
    ! Do this only where ice is present.  Elsewhere, usrf = topg.

    where (usrf > topg)     ! ice is present
       where (topg < 0.0d0)    ! marine-based ice
          where (topg*(1.0d0 - rhoo/rhoi) > usrf)  ! ice is floating
             thck = usrf / (1.0d0 - rhoi/rhoo)
          elsewhere   ! ice is grounded
             thck = usrf - topg
          endwhere
       elsewhere   ! land-based ice
          thck = usrf - topg
       endwhere
    endwhere

    ! Copy the new thickness and usrf to the model derived type
    model%geometry%thck = thck
    model%geometry%usrf = usrf

  end subroutine glissade_smooth_usrf

!****************************************************************************

  subroutine glissade_smooth_topography(model)

    ! Use a Laplacian smoother to smooth the input bed topography
    !TODO - This smoothing needs some more testing.  In particular, it is unclear how best to treat
    !        the ice thickness in regions that transition from grounded to floating
    !        when the topography is smoothed. Is it better to preserve thickness, or to
    !        increase thickness to keep the ice grounded?

    use glide_thck, only: glide_calclsrf
    use glissade_masks, only: glissade_get_masks
    use glissade_grid_operators, only: glissade_laplacian_smoother

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    ! local variables

    integer, dimension(model%general%ewn, model%general%ewn) :: &
         ice_mask,           &  !
         floating_mask

    real(dp), dimension(model%general%ewn, model%general%nsn) :: &
         topg_smoothed          ! bed topography after smoothing

    integer :: i, j
    integer :: nx, ny
    integer :: itest, jtest, rtest
    ! The following variables give the boundaries of a box going from itest-3 to itest+3,
    ! and jtest-3 to jtest+3, but limited to stay within the range of the local points
    ! owned by rdiag_local
    integer :: itest_m3, itest_p3, jtest_m3, jtest_p3

    logical, parameter :: verbose_smooth_topg = .false.

    ! Copy some model variables to local variables

    nx = model%general%ewn
    ny = model%general%nsn

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
       itest_m3 = max(itest-3, 1)
       itest_p3 = min(itest+3, nx)
       jtest_m3 = max(jtest-3, 1)
       jtest_p3 = min(jtest+3, ny)
    endif

    ! compute the initial upper surface elevation (to be held fixed under smoothing of bed topography)
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    ! compute initial mask
    ! Modify glissade_get_masks so that 'parallel' is not needed
    call glissade_get_masks(nx,                  ny,                    &
                            model%parallel,                             &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   0.0d0,                 &  ! thklim = 0
                            ice_mask,                                   &
                            floating_mask = floating_mask)

    call glissade_laplacian_smoother(model%general%ewn, model%general%nsn,  &
                                     model%geometry%topg, topg_smoothed,    &
                                     npoints_stencil = 5)

    !WHL - debug - Try doing less smoothing than the smoother computes
    model%geometry%topg = 0.50d0 * (model%geometry%topg + topg_smoothed)

    ! Given the smoothed topography, adjust the input thickness such that usrf is unchanged.
    where (model%geometry%topg - model%climate%eus < 0.0d0)  ! marine-based ice
       where (ice_mask == 1 .and. floating_mask == 0)
          ! Ice was grounded before smoothing of topography; assume it is still grounded.
          ! This means that where topg has been lowered, we will thicken the ice.
          model%geometry%thck = model%geometry%usrf - model%geometry%topg
       elsewhere
          ! Ice was floating before smoothing of topography.
          ! It may now be grounded where topg has been raised, in which case we move lsrf up to meet the topography.
          model%geometry%lsrf = max(model%geometry%lsrf, model%geometry%topg)
          model%geometry%thck = model%geometry%usrf - model%geometry%lsrf
       endwhere
    elsewhere   ! land-based ice
       model%geometry%thck = model%geometry%usrf - model%geometry%topg
    endwhere

    !WHL - usrf for debugging only
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

  end subroutine glissade_smooth_topography

!****************************************************************************

  subroutine glissade_adjust_topography(model)

    ! Adjust the input bed topography in a specified region.
    ! For example, we may want to raise the topography close to the surface in a region
    !  where the ice is not sufficiently grounded, and the data are not well constrained.
    ! Note: So far, this subroutine has been used to raise eastern Thwaites topography.
    !       It has not been used to lower topography.

    use glide_thck, only: glide_calclsrf  ! TODO - Make this a glissade subroutine (e.g., in this module)

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    ! local variables

    real(dp), dimension(:,:), allocatable :: &
         topg           ! topography with units of m

    integer :: i, j
    integer :: nx, ny
    integer :: itest, jtest, rtest
    ! The following variables give the boundaries of a box going from itest-3 to itest+3,
    ! and jtest-3 to jtest+3, but limited to stay within the range of the local points
    ! owned by rdiag_local
    integer :: itest_m3, itest_p3, jtest_m3, jtest_p3

    real(dp) :: factor

    real(dp) :: &
         xmin, xmax, ymin, ymax  ! x and y boundaries of adjusted region

    ! Note: There are two ways to use these parameters.
    ! (1) If topg_max_adjust > topg_no_adjust, then we change high topography (topg > topg_max_adjust)
    !     by topg_delta and leave low topography (topg < topg_no_adjust) unchanged.
    ! (2) If topg_max_adjust < topg_no_adjust, then we change low topography (topg < topg_max_adjust)
    !     by topg_delta and leave high topography (topg > topg_no_adjust) unchanged.
    ! Between topg_no_adjust and topg_max_adjust, the adjustment is phased in linearly.

    real(dp) :: &
         topg_no_adjust, &    ! elevation (m) beyond which there is no adjustment
         topg_max_adjust, &   ! elevation (m) beyond which there is full adjustment (by topg_delta)
         topg_delta           ! max change in topography (m); can be either sign

    logical, parameter :: verbose_adjust_topg = .true.

    ! Copy some model variables to local variables

    nx = model%general%ewn
    ny = model%general%nsn

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
       itest_m3 = max(itest-3, 1)
       itest_p3 = min(itest+3, nx)
       jtest_m3 = max(jtest-3, 1)
       jtest_p3 = min(jtest+3, ny)
    endif

    xmin = model%paramets%adjust_topg_xmin
    xmax = model%paramets%adjust_topg_xmax
    ymin = model%paramets%adjust_topg_ymin
    ymax = model%paramets%adjust_topg_ymax
    topg_no_adjust = model%paramets%adjust_topg_no_adjust
    topg_max_adjust = model%paramets%adjust_topg_max_adjust
    topg_delta = model%paramets%adjust_topg_delta

    if (verbose_adjust_topg .and. this_rank == rtest) then
       i = itest
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'Adjust input topography, diag point: r, i ,j =', rtest, itest, jtest
       write(iulog,*) 'x1, y1 =', model%general%x1(i), model%general%y1(j)
       write(iulog,*) 'thck, topg =', model%geometry%thck(i,j), model%geometry%topg(i,j)
       write(iulog,*) 'xmin, xmax =', xmin, xmax
       write(iulog,*) 'ymin, ymax =', ymin, ymax
       write(iulog,*) 'topg_no_adjust, topg_max_adjust (m) =', topg_no_adjust, topg_max_adjust
       write(iulog,*) 'topg_delta =', topg_delta
    endif

    ! Compute the lower and upper ice surface before the adjustment
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    !TODO - Use model%geometry%topg - model%climate%eus?
    allocate(topg(model%general%ewn, model%general%nsn))
    topg = model%geometry%topg

    ! Apply the topographic correction.
    ! Case 1: topg_max_adjust > topg_no_adjust; change higher topography and leave lower topography unchanged
    ! Case 2: topg_max_adjust < topg_no_adjust; change lower topography and leave higher topography unchanged

    if (topg_max_adjust > topg_no_adjust) then

       ! Where topg > topg_max_adjust, apply the max correction, topg_delta.
       ! Where topg < topg_no_adjust, apply no correction.
       ! Where topg_no_adjust < topg < topg_max_adjust, phase in the correction linearly.

       do j = 1, ny
          do i = 1, nx
             if (model%general%x1(i) >= xmin .and. model%general%x1(i) <= xmax .and. &
                 model%general%y1(j) >= ymin .and. model%general%y1(j) <= ymax) then
                if (topg(i,j) > topg_no_adjust) then
                   factor = min((topg(i,j) - topg_no_adjust)/(topg_max_adjust - topg_no_adjust), 1.0d0)
                   topg(i,j) = topg(i,j) + factor * topg_delta
                endif
             endif
          enddo
       enddo

    elseif (topg_max_adjust < topg_no_adjust) then

       ! Where topg < topg_max_adjust, apply the max correction, topg_delta.
       ! Where topg > topg_no_adjust, apply no correction.
       ! Where topg_max_adjust < topg < topg_no_adjust, phase in the correction linearly.

       do j = 1, ny
          do i = 1, nx
             if (model%general%x1(i) >= xmin .and. model%general%x1(i) <= xmax .and. &
                 model%general%y1(j) >= ymin .and. model%general%y1(j) <= ymax) then
                if (topg(i,j) < topg_no_adjust) then
                   factor = min((topg_no_adjust - topg(i,j))/(topg_no_adjust - topg_max_adjust), 1.0d0)
                   topg(i,j) = topg(i,j) + factor * topg_delta
                endif
             endif
          enddo
       enddo

    endif

    model%geometry%topg = topg
    deallocate(topg)

    ! In some cells, the new lower surface (usrf - thck) may lie below the topography.
    ! In these cells, reduce the ice thickness such that lsrf = topg, preserving the input value of usrf.

    where (model%geometry%usrf - model%geometry%thck < model%geometry%topg)
       model%geometry%thck = model%geometry%usrf - model%geometry%topg
    endwhere

  end subroutine glissade_adjust_topography

!****************************************************

  !TODO - Calls to this subroutine could be replaced by inline calls to parallel_global_sum_patch
  subroutine glissade_basin_sum(&
       nx,           ny,            &
       parallel,                    &
       nbasin,       basin_number,  &
       rmask,                       &
       field_2d,                    &
       field_basin_sum)

    ! For a given 2D input field, compute the sum over a basin.
    ! The sum is taken over grid cells with mask = 1.
    ! All cells are weighted equally.

    use cism_parallel, only: parallel_global_sum_patch

    integer, intent(in) :: &
         nx, ny                    !> number of grid cells in each dimension

    type(parallel_type), intent(in) :: &
         parallel                  !> info for parallel communication

    integer, intent(in) :: &
         nbasin                    !> number of basins

    integer, dimension(nx,ny), intent(in) :: &
         basin_number              !> basin ID for each grid cell

    ! Note: For the next two fields, the dimension can be either (nx,ny) or (nx-1,ny-1)
    real(dp), dimension(:,:), intent(in) :: &
         rmask,                 &  !> real mask for weighting the input field
         field_2d                  !> input field to be averaged over basins

    real(dp), dimension(nbasin), intent(out) :: &
         field_basin_sum           !> basin-sum output field

    !TODO - Replace sumcell with sumarea, and pass in cell area.
    !       Current algorithm assumes all cells with mask = 1 have equal weight.

    field_basin_sum = parallel_global_sum_patch(rmask*field_2d, nbasin, basin_number, parallel)

  end subroutine glissade_basin_sum

!****************************************************

  subroutine glissade_basin_average(&
       nx,           ny,            &
       parallel,                    &
       nbasin,       basin_number,  &
       rmask,                       &
       field_2d,                    &
       field_basin_avg)

    ! For a given 2D input field, compute the average over a basin.
    ! The average is taken over grid cells with mask = 1.
    ! All cells are weighted equally.
    ! Note: This subroutine assumes an input field located at cell centers

    use cism_parallel, only: parallel_global_sum_patch

    integer, intent(in) :: &
         nx, ny                    !> number of grid cells in each dimension

    type(parallel_type), intent(in) :: &
         parallel                  !> info for parallel communication

    integer, intent(in) :: &
         nbasin                    !> number of basins

    integer, dimension(nx,ny), intent(in) :: &
         basin_number              !> basin ID for each grid cell

    ! Note: For the next two fields, the dimension can be either (nx,ny) or (nx-1,ny-1)
    real(dp), dimension(:,:), intent(in) :: &
         rmask,                  & !> real mask for weighting the value in each cell
         field_2d                  !> input field to be averaged over basins

    real(dp), dimension(nbasin), intent(out) :: &
         field_basin_avg           !> basin-average output field

    ! local variables

    integer :: nb

    !TODO - Replace sumcell with sumarea, and pass in cell area.
    !       Current algorithm assumes all cells with mask = 1 have equal weight.

    real(dp), dimension(nbasin) ::  &
         summask_global,         & ! sum of mask in each basin on full domain
         sumfield_global           ! sum of field over full domain

    summask_global = parallel_global_sum_patch(rmask, nbasin, basin_number, parallel)
    sumfield_global = parallel_global_sum_patch(rmask*field_2d, nbasin, basin_number, parallel)

    do nb = 1, nbasin
       if (summask_global(nb) > tiny(0.0d0)) then
          field_basin_avg(nb) = sumfield_global(nb)/summask_global(nb)
       else
          field_basin_avg(nb) = 0.0d0
       endif
    enddo

  end subroutine glissade_basin_average

!****************************************************************************

  subroutine glissade_rms_error(&
       nx,       ny,          &
       mask,     parallel,    &
       field,    field_ref,   &
       rmse)

    use cism_parallel, only: parallel_global_sum

    ! Compute the root-mean-square error of an input field relative to a reference field.
    ! Typically, the input field would be computed by CISM, with the reference field
    !  based on observations.

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny                  ! number of cells in x and y direction on input grid (global)

    integer, dimension(nx,ny), intent(in) :: &
         mask                    ! = 1 for the domain over which the rmse is computed

    type(parallel_type), intent(in) :: &
         parallel                ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) :: &
         field,                & ! 2D model field
         field_ref               ! reference field

    real(dp), intent(out) :: &
         rmse                    ! root-mean-square error

    ! local variables

    real(dp), dimension(nx,ny) :: &
         sq_diff              ! |field - field_ref|^2

    real(dp) :: &
         sum_sq_diff,       & ! global sum of sq_diff
         ncells               ! number of global cells with mask = 1

    ncells = parallel_global_sum(mask, parallel)

    sq_diff = (abs(field - field_ref))**2
    sum_sq_diff = parallel_global_sum(sq_diff, parallel, mask)

    if (ncells > 0.0d0) then
       rmse = sqrt(sum_sq_diff/ncells)
    else
       rmse = 0.0d0
    endif

  end subroutine glissade_rms_error

!***********************************************************************

  subroutine glissade_usrf_to_thck(usrf, topg, eus, thck)

    ! Given the bed topography and upper ice surface elevation, compute the ice thickness.
    ! The ice is assumed to satisfy a flotation condition.
    ! That is, if topg - eus < 0 (marine-based ice), and if the upper surface is too close
    !  to sea level to ground the ice, then the ice thickness is chosen to satisfy
    !  rhoi*H = -rhoo*(topg-eus).
    ! Note: usrf, topg, eus and thck must all have the same units (usually but not necessarily meters).

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

  end subroutine glissade_usrf_to_thck

!***********************************************************************

  subroutine glissade_thck_to_usrf(thck, topg, eus, usrf)

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

  end subroutine glissade_thck_to_usrf

!***********************************************************************

  subroutine glissade_edge_fluxes(&
        nx,        ny,        &
        dew,       dns,       &
        itest,     jtest,  rtest, &
        thck,                 &
        uvel,      vvel,      &
        flux_e,    flux_n)

    use cism_parallel, only: nhalo

    ! Compute ice volume fluxes across each cell edge

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                & ! number of cells in x and y direction on input grid (global)
         itest, jtest, rtest

    real(dp), intent(in) :: &
         dew, dns                 ! cell edge lengths in EW and NS directions (m)

    real(dp), dimension(nx,ny), intent(in) :: &
         thck                     ! ice thickness (m) at cell centers

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
         uvel, vvel               ! vertical mean velocity (m/s) at cell corners

    real(dp), dimension(nx,ny), intent(out) :: &
         flux_e, flux_n           ! ice volume fluxes (m^3/yr) at cell edges

    ! local variables

    integer :: i, j
    real(dp) :: thck_edge, u_edge, v_edge
    logical, parameter :: verbose_edge_fluxes = .false.

    ! loop over locally owned edges
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          ! east edge volume flux
          thck_edge = 0.5d0 * (thck(i,j) + thck(i+1,j))
          u_edge = 0.5d0 * (uvel(i,j-1) + uvel(i,j))
          flux_e(i,j) = thck_edge * u_edge * dns  ! m^3/yr

          ! north edge volume flux
          thck_edge = 0.5d0 * (thck(i,j) + thck(i,j+1))
          v_edge = 0.5d0 * (vvel(i-1,j) + vvel(i,j))
          flux_n(i,j) = thck_edge * v_edge * dew  ! m^3/yr

          if (verbose_edge_fluxes .and. this_rank == rtest .and. i==itest .and. j==jtest) then
             write(iulog,*) 'East  flux: rank, i, j, H, u, flx =', &
                  rtest, itest, jtest, thck_edge, u_edge, flux_e(i,j)
             write(iulog,*) 'North flux: rank, i, j, H, v, flx =', &
                  rtest, itest, jtest, thck_edge, v_edge, flux_n(i,j)
          endif

       enddo
    enddo

  end subroutine glissade_edge_fluxes

!***********************************************************************

  subroutine glissade_input_fluxes(&
        nx,      ny,            &
        dew,     dns,           &
        dt,                     &
        itest,   jtest,  rtest, &
        thck,                   &
        uvel,    vvel,          &
        flux_in,                &
        parallel)

    use glimmer_physcon, only: scyr
    use cism_parallel, only: nhalo, parallel_halo, staggered_parallel_halo

    ! Compute ice volume fluxes into a cell from each neighboring cell

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                & ! number of cells in x and y direction on input grid (global)
         itest, jtest, rtest

    real(dp), intent(in) :: &
         dew, dns,              & ! cell edge lengths in EW and NS directions (m)
         dt                       ! timestep (s)

    real(dp), dimension(nx,ny), intent(in) :: &
         thck                     ! ice thickness (m) at cell centers

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
         uvel, vvel               ! vertical mean velocity (m/s) at cell corners

    real(dp), dimension(-1:1,-1:1,nx,ny), intent(out) :: &
         flux_in                  ! ice volume fluxes (m^3/s) into cell from each neighbor cell

    type(parallel_type), intent(in) :: parallel   ! info for parallel communication

    ! local variables

    integer :: i, j, ii, jj

    real(dp) :: &
         u_sw, u_se, u_ne, u_nw,    & ! u velocity components at each vertex
         v_sw, v_se, v_ne, v_nw       ! u velocity components at each vertex

    real(dp) :: &
         area_w, area_s, area_e, area_n,   & ! area flux from each neighbor cell
         area_sw, area_se, area_ne, area_nw

    logical, parameter :: verbose_input_fluxes = .false.

    ! halo updates for thickness and velocity

    call parallel_halo(thck, parallel)
    call staggered_parallel_halo(uvel, parallel)
    call staggered_parallel_halo(vvel, parallel)

    ! initialize
    flux_in(:,:,:,:) = 0.0d0

    ! Estimate the ice volume flux into each cell from each neighbor.
    ! Note: flux_in(0,0,:,:) = 0 since there is no flux from a cell into itself.
    ! The loop includes one row of halo cells.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          ! Compute the upwind velocity components at each vertex
          u_sw = uvel(i-1,j-1)
          v_sw = vvel(i-1,j-1)
          u_se = uvel(i,j-1)
          v_se = vvel(i,j-1)
          u_ne = uvel(i,j)
          v_ne = vvel(i,j)
          u_nw = uvel(i-1,j)
          v_nw = vvel(i-1,j)

          ! Estimate the area fluxes (m^2/s) into this cells from each edge neighbor
          ! Note: The first line on the RHS accounts for the velocity component
          !        perpendicular to the edge; this is a rectangle area.
          !       The next two lines are corrections proportional to the velocity
          !        component parallel to the edge; these are triangle areas.
          area_w = 0.5d0*(max( u_nw,0.0d0) + max( u_sw,0.0d0))*dns   &
                 - 0.5d0* max( u_nw,0.0d0) * max( v_nw,0.0d0)*dt     &
                 - 0.5d0* max( u_sw,0.0d0) * max(-v_sw,0.0d0)*dt
          area_s = 0.5d0*(max( v_sw,0.0d0) + max( v_se,0.0d0))*dew   &
                 - 0.5d0* max(-u_sw,0.0d0) * max( v_sw,0.0d0)*dt     &
                 - 0.5d0* max( u_se,0.0d0) * max( v_se,0.0d0)*dt
          area_e = 0.5d0*(max(-u_se,0.0d0) + max(-u_ne,0.0d0))*dns   &
                 - 0.5d0* max(-u_se,0.0d0) * max(-v_se,0.0d0)*dt     &
                 - 0.5d0* max(-u_ne,0.0d0) * max( v_ne,0.0d0)*dt
          area_n = 0.5d0*(max(-v_ne,0.0d0) + max(-v_nw,0.0d0))*dew   &
                 - 0.5d0* max( u_ne,0.0d0) * max(-v_ne,0.0d0)*dt     &
                 - 0.5d0* max(-u_nw,0.0d0) * max(-v_nw,0.0d0)*dt

          ! Estimate the area fluxes (m^2/s) from each diagonal neighbor.
          ! These are rectangle areas.
          area_sw = max( u_sw,0.0d0)*max( v_sw,0.0d0)*dt
          area_se = max(-u_se,0.0d0)*max( v_se,0.0d0)*dt
          area_ne = max(-u_ne,0.0d0)*max(-v_ne,0.0d0)*dt
          area_nw = max( u_nw,0.0d0)*max(-v_nw,0.0d0)*dt

          ! Estimate the volume fluxes from each edge neighbor
          flux_in(-1, 0,i,j) = area_w * thck(i-1,j)
          flux_in( 0,-1,i,j) = area_s * thck(i,j-1)
          flux_in( 1, 0,i,j) = area_e * thck(i+1,j)
          flux_in( 0, 1,i,j) = area_n * thck(i,j+1)

          ! Estimate the volume fluxes from each diagonal neighbor
          flux_in(-1,-1,i,j) = area_sw * thck(i-1,j-1)
          flux_in( 1,-1,i,j) = area_se * thck(i+1,j-1)
          flux_in( 1, 1,i,j) = area_ne * thck(i+1,j+1)
          flux_in(-1, 1,i,j) = area_nw * thck(i-1,j+1)

          if (verbose_input_fluxes .and. this_rank == rtest .and. i==itest .and. j==jtest) then
             write(iulog,*) ' '
             write(iulog,*) 'upstream u (m/yr), this_rank, i, j:'
             write(iulog,'(3f18.12)') u_nw*scyr, u_ne*scyr
             write(iulog,'(3f18.12)') u_sw*scyr, u_se*scyr
             write(iulog,*) ' '
             write(iulog,*) 'upstream v (m/yr):'
             write(iulog,'(3f18.12)') v_nw*scyr, v_ne*scyr
             write(iulog,'(3f18.12)') v_sw*scyr, v_se*scyr
             write(iulog,*) ' '
             write(iulog,*) 'Input area fluxes (km^2/yr):'
             write(iulog,'(3f18.12)') area_nw*scyr/1.0d6, area_n*scyr/1.0d6, area_ne*scyr/1.0d6
             write(iulog,'(3f18.12)') area_w *scyr/1.0d6,       0.0d0,       area_e *scyr/1.0d6
             write(iulog,'(3f18.12)') area_sw*scyr/1.0d6, area_s*scyr/1.0d6, area_se*scyr/1.0d6
             write(iulog,*) 'Total =', &
                  (area_w + area_s + area_e + area_n + area_sw + area_se + area_ne + area_nw)*scyr/1.0d6
             write(iulog,*) ' '
             write(iulog,*) 'Estimated edge area fluxes:'
             area_w = 0.5d0*(max( u_nw,0.0d0) + max( u_sw,0.0d0))*dns
             area_s = 0.5d0*(max( v_sw,0.0d0) + max( v_se,0.0d0))*dew
             area_e = 0.5d0*(max(-u_se,0.0d0) + max(-u_ne,0.0d0))*dns
             area_n = 0.5d0*(max(-v_ne,0.0d0) + max(-v_nw,0.0d0))*dew
             write(iulog,*) 'area_w =', area_w*scyr/1.0e6
             write(iulog,*) 'area_s =', area_s*scyr/1.0e6
             write(iulog,*) 'area_e =', area_e*scyr/1.0e6
             write(iulog,*) 'area_n =', area_n*scyr/1.0e6
             write(iulog,*) 'Total =', (area_w + area_s + area_e + area_n)*scyr/1.0d6
             write(iulog,*) ' '
             write(iulog,*) 'Input ice volume fluxes (km^3/yr):'
             do jj = 1,-1,-1
                do ii = -1,1
                   write(iulog,'(f15.8)',advance='no') flux_in(ii,jj,i,j)*scyr/1.0d9
                enddo
                write(iulog,*) ' '
             enddo
          endif

       enddo   ! i
    enddo   ! j

    do jj = -1, 1
       do ii = -1, 1
          call parallel_halo(flux_in(ii,jj,:,:), parallel)
       enddo
    enddo

  end subroutine glissade_input_fluxes

!***********************************************************************

  subroutine glissade_quadrant_sum(&
       nx,        ny,            &
       parallel,                 &
       field,     quadrant_sum)

    ! Integrate a field over each of 4 quadrants.
    ! This can be useful in idealized experiments like CalvingMIP to check for
    !  violations of reflectional or rotational symmetry.
    ! Note: These sums are not independent of processor count
    ! TODO: Make them reproducible, using quadrant masks?

    use cism_parallel, only: nhalo, parallel_global_sum_patch, parallel_globalindex, gather_var

    ! Input/output arguments

    integer, intent(in) :: &
         nx, ny                   ! number of local cells in x and y direction on input grid

    type(parallel_type), intent(in) :: parallel   ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) :: &
         field                    ! 2D input field

    real(dp), dimension(4), intent(out) :: &
         quadrant_sum             ! global sum over each of 4 quadrants

    logical, parameter :: check_asymmetry = .true.

    ! Local variables

    integer :: i, j
    integer :: ig, jg             ! i and j indices on the global grid
    integer :: nxg, nyg           ! dimensions of global domain
    integer :: nx2, ny2           ! nx/2 and ny/2 (if nx and ny are even)
                                  ! (nx-1)/2 and (ny-1)/2 (if nx and ny are odd)

    integer, dimension(nx,ny) :: &
         quadrant_mask            ! mask assigning each cell to a quadrent (1, 2, 3 or 4)

    real(dp), dimension(:,:), allocatable :: field_global

    real(dp) :: meanval, diff
    real(dp), parameter :: symmetry_tol = 1.0d-5    ! tolerance level for asymmetry


    ! Compute a mask that assigns each cell to one of 4 quadrants.
    ! Note: If nx or ny is odd, the middle row or column is excluded from the quadrant sums.

    nxg = parallel%global_ewn
    nyg = parallel%global_nsn

    if (mod(nxg,2) == 0) then   ! global_ewn is even
       nx2 = nxg/2
    else   ! nx is odd
       nx2 = (nxg-1)/2
    endif

    if (mod(nyg,2) == 0) then   ! global_nsn is even
       ny2 = nyg/2
    else   ! ny is odd
       ny2 = (nyg-1)/2
    endif

    quadrant_mask(:,:) = 0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          call parallel_globalindex(i, j, ig, jg, parallel)
          if (ig > nx2) then
             if (jg > ny2) then   ! NE quadrant
                quadrant_mask(i,j) = 1
             else   ! jg <= ny2; SE quadrant
                quadrant_mask(i,j) = 4
             endif
          else   ! ig <= nx2
             if (jg > ny2) then   ! NW quadrant
                quadrant_mask(i,j) = 2
             else   ! jg <= ny2; SW quadrant
                quadrant_mask(i,j) = 3
             endif
          endif
       enddo
    enddo

    ! Compute the global sums
    ! Note: These sums are reproducible if reproducible_sums = .true.
    quadrant_sum(:) = parallel_global_sum_patch(field, 4, quadrant_mask, parallel)

    if (check_asymmetry) then

       call gather_var(field, field_global, parallel)

       if (main_task) then

          ! Identify asymmetries in reflection across the y-axis
          if (abs(quadrant_sum(1) + quadrant_sum(4) - quadrant_sum(2) - quadrant_sum(3)) > symmetry_tol) then
             do j = 1, nyg
                do i = 1, nx2
                   if (abs(field_global(i,j)) > eps11 .or. abs(field_global(nxg-i+1,j)) > eps11) then
                      meanval = 0.5d0 * (field_global(i,j) + field_global(nxg-i+1,j))
                      diff = abs(field_global(i,j) - field_global(nxg-i+1,j))
                      if (diff > meanval*symmetry_tol) then
                         write(iulog,*) 'Warning, y-reflection asymmetry: i, j, val(i,j), val(nxg-i+1,j), diff/mean:', &
                              i, j, field_global(i,j), field_global(nxg-i+1,j), diff/meanval
                      endif
                   endif
                enddo
             enddo
          endif

          ! Identify asymmetries in reflection across the x-axis
          if (abs(quadrant_sum(1) + quadrant_sum(2) - quadrant_sum(3) - quadrant_sum(4)) > symmetry_tol) then
             do j = 1, nyg
                do i = 1, nx2
                   if (abs(field_global(i,j)) > eps11 .or. abs(field_global(i,nyg-j+1)) > eps11) then
                      meanval = 0.5d0 * (field_global(i,j) + field_global(i,nyg-j+1))
                      diff = abs(field_global(i,j) - field_global(i,nyg-j+1))
                      if (diff > meanval*symmetry_tol) then
                         write(iulog,*) 'Warning, x-reflection asymmetry: i, j, val(i,j), val(i,nyg-j+1), diff/mean:', &
                              i, j, field_global(i,j), field_global(i,nyg-j+1), diff/meanval
                      endif
                   endif
                enddo
             enddo
          endif

          if (allocated(field_global)) deallocate(field_global)

       endif   ! main_task
    endif   ! check_asymmetry

  end subroutine glissade_quadrant_sum

!***********************************************************************

  subroutine glissade_bounding_box(&
       nx,              ny,           &
       dx,              dy,           &
       x1,              y1,           &
       x_point,         y_point,      &
       field1,                        &
       field1_at_point,               &
       field2,                        &
       field2_at_point)

    ! Input/output arguments

    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy
    real(dp), dimension(nx), intent(in) :: x1
    real(dp), dimension(ny), intent(in) :: y1
    real(dp), intent(in) :: x_point, y_point
    real(dp), dimension(nx,ny), intent(in) :: field1
    real(dp), intent(out) :: field1_at_point
    real(dp), dimension(nx,ny), intent(in), optional :: field2
    real(dp), intent(out), optional :: field2_at_point

    ! Local variables

    integer :: i, j, ipt, jpt

    real(dp), dimension(2,4) :: box_coords    ! x and y coordinates of 4 box corners
    real(dp), dimension(4) :: box_field
    integer, dimension(4) :: box_mask

    logical, parameter :: verbose_bounding_box = .false.
!!    logical, parameter :: verbose_bounding_box = .true.

    ipt = 0; jpt = 0
    do i = 1, nx-1
       if (x1(i) <= x_point .and. x1(i+1) > x_point) then
          ipt = i   ! i index for SW corner of box
       endif
    enddo
    do j = 1, ny-1
       if (y1(j) <= y_point .and. y1(j+1) > y_point) then
          jpt = j   ! j index for SW corner of box
       endif
    enddo
    if (ipt == 0 .or. jpt == 0) then
       write(iulog,*) 'glissade_bounding_box, bad location: rank, i, j =', this_rank, ipt, jpt
    endif
    if (verbose_bounding_box) then
       write(iulog,*) 'Point coordinates:', x_point, y_point
       write(iulog,*) '  Point is on rank', this_rank
       write(iulog,*) '  Point is bounded by i =', ipt, ipt+1
       write(iulog,*) '  Point is bounded by j =', jpt, jpt+1
    endif

    ! Copy coordinates into an array
    ! In arrays with 4 indices, the box corners are ordered (1) SW, (2) SE, (3) NE, (4) NW
    box_coords(1,1) = x1(ipt)       ! SW cell
    box_coords(2,1) = y1(jpt)
    box_coords(1,2) = x1(ipt+1)     ! SE cell
    box_coords(2,2) = y1(jpt)
    box_coords(1,3) = x1(ipt+1)     ! NE cell
    box_coords(2,3) = y1(jpt+1)
    box_coords(1,4) = x1(ipt)       ! NW cell
    box_coords(2,4) = y1(jpt+1)

    ! Copy field1 into an array
    ! Assume that any nonzero values are valid
    box_field(:) = 0.0d0
    box_mask(:) = 0
    if (field1(ipt,jpt) > eps11) then
       box_field(1) = field1(ipt,jpt)
       box_mask(1) = 1
    endif
    if (field1(ipt+1,jpt) > eps11) then
       box_field(2) = field1(ipt+1,jpt)
       box_mask(2) = 1
    endif
    if (field1(ipt+1,jpt+1) > eps11) then
       box_field(3) = field1(ipt+1,jpt+1)
       box_mask(3) = 1
    endif
    if (field1(ipt,jpt+1) > eps11) then
       box_field(4) = field1(ipt,jpt+1)
       box_mask(4) = 1
    endif

    ! Make sure at least one box corner has a nonzero value.
    ! If not, then extend the box to the south or north.
    ! A general solution to this problem would require careful logic,
    ! but the following logic works for the Caprona axes on a 5-km grid.

    if (sum(box_mask) == 0) then
       if (field1(ipt,jpt-1) > eps11 .or. field1(ipt+1,jpt-1) > eps11) then
          box_coords(2,1) = y1(jpt-1)  ! new SW cell
          box_coords(2,2) = y1(jpt-1)  ! new SE cell
          if (field1(ipt,jpt-1) > eps11) then
             box_field(1) = field1(ipt,jpt-1)
             box_mask(1) = 1
          endif
          if (field1(ipt+1,jpt-1) > eps11) then
             box_field(2) = field1(ipt+1,jpt-1)
             box_mask(2) = 1
          endif
       elseif (field1(ipt,jpt+1) > eps11 .or. field1(ipt+1,jpt+11) > eps11) then
          box_coords(2,1) = y1(jpt+1)  ! new NW cell
          box_coords(2,2) = y1(jpt+1)  ! new NE cell
          if (field1(ipt,jpt+1) > eps11) then
             box_field(1) = field1(ipt,jpt+1)
             box_mask(1) = 1
          endif
          if (field1(ipt+1,jpt+1) > eps11) then
             box_field(2) = field1(ipt+1,jpt+1)
             box_mask(2) = 1
          endif
       endif
    endif

    if (sum(box_mask) == 0) then
       call write_log('Warning, all corners of bounding box have field1 = 0', GM_WARNING)
    endif

    ! Compute field1 at the point inside the box
    call bounding_box_interpolate(&
         dx,              dy,       &
         x_point,         y_point,  &
         box_coords(:,:),           &
         box_field(:),              &
         box_mask(:),               &
         field1_at_point)

    ! Repeat for field2, if present

    if (present(field2) .and. present(field2_at_point)) then

       ! Copy field2 into an array
       ! Assume that any nonzero values are valid
       box_field(:) = 0.0d0
       box_mask(:) = 0
       if (field2(ipt,jpt) /= 0.0d0) then
          box_field(1) = field2(ipt,jpt)
          box_mask(1) = 1
       endif
       if (field2(ipt+1,jpt) /= 0.0d0) then
          box_field(2) = field2(ipt+1,jpt)
          box_mask(2) = 1
       endif
       if (field2(ipt+1,jpt+1) /= 0.0d0) then
          box_field(3) = field2(ipt+1,jpt+1)
          box_mask(3) = 1
       endif
       if (field2(ipt,jpt+1) /= 0.0d0) then
          box_field(4) = field2(ipt,jpt+1)
          box_mask(4) = 1
       endif

       if (sum(box_mask) == 0) then
          box_coords(2,1) = y1(jpt-1)  ! new SW cell
          box_coords(2,2) = y1(jpt-1)  ! new SE cell
          if (field2(ipt,jpt-1) > eps11) then
             box_field(1) = field2(ipt,jpt-1)
             box_mask(1) = 1
          endif
          if (field2(ipt+1,jpt-1) > eps11) then
             box_field(2) = field2(ipt+1,jpt-1)
             box_mask(2) = 1
          endif
       elseif (field2(ipt,jpt+1) > eps11 .or. field2(ipt+1,jpt+11) > eps11) then
          box_coords(2,1) = y1(jpt+1)  ! new NW cell
          box_coords(2,2) = y1(jpt+1)  ! new NE cell
          if (field2(ipt,jpt+1) > eps11) then
             box_field(1) = field2(ipt,jpt+1)
             box_mask(1) = 1
          endif
          if (field2(ipt+1,jpt+1) > eps11) then
             box_field(2) = field2(ipt+1,jpt+1)
             box_mask(2) = 1
          endif
       endif

       if (sum(box_mask) == 0) then
          call write_log('Warning, all corners of bounding box have field2 = 0', GM_WARNING)
       endif

       ! Compute field2 at the point inside the box
       call bounding_box_interpolate(&
            dx,              dy,       &
            x_point,         y_point,  &
            box_coords(:,:),           &
            box_field(:),              &
            box_mask(:),               &
            field2_at_point)

    endif   ! present(field2)

  end subroutine glissade_bounding_box

!***********************************************************************

  subroutine bounding_box_interpolate(&
       dx,           dy,      &
       x_point,      y_point, &
       corner_coords,         &
       corner_values,         &
       corner_mask,           &
       point_value)

    ! Given the values of a field at the four corners of a bounding box,
    !  make a linear approximation of the field value at a given point
    !  inside the box.
    ! This is cruder than a bilinear interpolation. It's intended to give
    !  an approximate answer, sometimes when cornerss are masked out
    !  (i.e., valid values are not available at all 4 corners.
    ! Note: This subroutine works for any distance units as long as units are consistent.

    ! Input/output arguments

    real(dp), intent(in) :: &
         dx, dy                  ! dimensions of the box

    real(dp), intent(in) :: &
         x_point, y_point        ! x and y coordinates of the point inside the box

    real(dp), dimension(2,4), intent(in) :: &
         corner_coords           ! x and y coordinates at each of 4 corners;
                                 ! ordering is SW, SE, NE, NW

    real(dp), dimension(4), intent(in) :: &
         corner_values           ! value of field at each corner; SW/SE/NE/NW ordering

    integer, dimension(4), intent(in) :: &
         corner_mask             ! = 1 for valid values, 0 for not valid

    real(dp), intent(out) :: &
         point_value             ! estimated field value at the selected point

    ! Local variables

    integer :: i, j
    real(dp) :: xp, yp                             ! coordinates of the point in the box
    real(dp) :: dxp, dyp                           ! coordinates of the point relative to a corner
    real(dp) :: x_sw, x_se, x_ne, x_nw             ! x coordinates for each corner
    real(dp) :: y_sw, y_se, y_ne, y_nw             ! y coordinates for each corner
    real(dp) :: f_sw, f_se, f_ne, f_nw             ! field values at each corner
    real(dp) :: f_e, f_w, f_n, f_s                 ! field values interpolated to edge midpoints
    real(dp) :: df_dx, df_dy                       ! field derivatives
    integer :: mask_sw, mask_se, mask_ne, mask_nw  ! mask values for each corner; = 1 for valid values, else 0
    integer :: mask_e, mask_w, mask_n, mask_s      ! mask values for each edge; = 1 for valid values, else 0

    logical, parameter :: verbose_bounding_box_interpolate = .false.

    ! Initialize
    ! These copies aren't strictly necessary, but the compass labels make things easier to visualize.

    xp = x_point
    yp = y_point

    f_sw = corner_values(1)
    f_se = corner_values(2)
    f_ne = corner_values(3)
    f_nw = corner_values(4)

    x_sw = corner_coords(1,1)
    y_sw = corner_coords(2,1)
    x_se = corner_coords(1,2)
    y_se = corner_coords(2,2)
    x_ne = corner_coords(1,3)
    y_ne = corner_coords(2,3)
    x_nw = corner_coords(1,4)
    y_nw = corner_coords(2,4)

    mask_sw = corner_mask(1)
    mask_se = corner_mask(2)
    mask_ne = corner_mask(3)
    mask_nw = corner_mask(4)

    ! assume edge values are valid unless both corner values are found to be masked out
    mask_e = 1
    mask_w = 1
    mask_s = 1
    mask_n = 1

    ! Interpolate field values to cell edges

    if (mask_se > 0 .and. mask_ne > 0) then
       f_e = 0.5d0 * (f_se + f_ne)
    elseif (mask_se > 0.0d0) then
       f_e = f_se
    elseif (mask_ne > 0.0d0) then
       f_e = f_ne
    else
       mask_e = 0
    endif

    if (mask_sw > 0 .and. mask_nw > 0) then
       f_w = 0.5d0 * (f_sw + f_nw)
    elseif (mask_sw > 0.0d0) then
       f_w = f_sw
    elseif (mask_nw > 0.0d0) then
       f_w = f_nw
    else
       mask_w = 0
    endif

    if (mask_nw > 0 .and. mask_ne > 0) then
       f_n = 0.5d0 * (f_nw + f_ne)
    elseif (mask_nw > 0.0d0) then
       f_n = f_nw
    elseif (mask_ne > 0.0d0) then
       f_n = f_ne
    else
       mask_n = 0
    endif

    if (mask_sw > 0 .and. mask_se > 0) then
       f_s = 0.5d0 * (f_sw + f_se)
    elseif (mask_sw > 0.0d0) then
       f_s = f_sw
    elseif (mask_se > 0.0d0) then
       f_s = f_se
    else
       mask_s = 0
    endif

    ! Estimate the derivatives
    ! Requires at least one valid value per edge to compute a derivative

    if (mask_e > 0 .and. mask_w > 0) then
       df_dx = (f_e - f_w)/dx
    else
       df_dx = 0
    endif

    if (mask_n > 0 .and. mask_s > 0) then
       df_dy = (f_n - f_s)/dy
    else
       df_dy = 0
    endif

    ! Estimate the value at the point inside the box.
    ! (Still computes a value if the corner is outside the box,
    !  but there's no guarantee the extrapolation will be accurate.))
    ! At least one corner should have a valid value.

    if (mask_sw > 0) then
       dxp = xp - x_sw
       dyp = yp - y_sw
       point_value = f_sw + df_dx*dxp + df_dy*dyp
    elseif (mask_se > 0) then
       dxp = xp - x_se
       dyp = yp - y_se
       point_value = f_se + df_dx*dxp + df_dy*dyp
    elseif (mask_ne > 0) then
       dxp = xp - x_ne
       dyp = yp - y_ne
       point_value = f_ne + df_dx*dxp + df_dy*dyp
    elseif (mask_nw > 0) then
       dxp = xp - x_nw
       dyp = yp - y_nw
       point_value = f_nw + df_dx*dxp + df_dy*dyp
    else
       point_value = 0.0d0
       call write_log('Warning, glissade_bounding_box: no valid values', GM_WARNING)
    endif

    if (verbose_bounding_box_interpolate) then
       write(6,*) 'In bounding_box_interpolate, rank =', this_rank
       write(6,*) 'point coordinates =', xp, yp
       write(6,*) 'df/dx, df/dy:', df_dx, df_dy
       write(6,*) 'dxp, dyp:', dxp, dyp
       write(6,*) 'point value =', point_value
    endif

  end subroutine bounding_box_interpolate

!***********************************************************************

  ! subroutines belonging to the write_array_to_file interface
  subroutine write_array_to_file_real8_2d(arr, fileunit, filename, parallel, write_binary)

    ! Copy the input array into a global array and write all values to an output file.
    ! This can be useful for debugging, if we want to find differences between two fields
    ! (e.g., in two different runs).
    ! This version writes out 64-bit character strings corresponding to the binary representation
    ! of each floating-point variable. This can be useful for BFB comparisons.
    ! Sometimes, two floating-point variables appear to have the same values in base 10,
    ! when the last few bits actually vary.
    !TODO - Allow either float or binary output

    use glimmer_utils, only: double_to_binary
    use cism_parallel, only: gather_var

    real(dp), dimension(:,:), intent(in) :: arr
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    type(parallel_type), intent(in) :: parallel
    logical, intent(in), optional :: write_binary

    integer :: i, j
    character(len=64) :: binary_str
    real(dp), dimension(:,:), allocatable :: arr_global
    logical :: binary_output

    if (present(write_binary)) then
       binary_output = write_binary
    else
       binary_output = .false.
    endif

    call gather_var(arr, arr_global, parallel)
    if (main_task) then
       open(unit=fileunit, file=trim(filename), status='replace', position='append')

       if (binary_output) then
          do j = 1, parallel%global_nsn
             do i = 1, parallel%global_ewn
                call double_to_binary(arr_global(i,j), binary_str)
                write (fileunit, '(2i6,a4,a64)') i, j, '    ', binary_str
             enddo
          enddo
       else
          do j = 1, parallel%global_nsn
             do i = 1, parallel%global_ewn
                write (fileunit, '(2i6,a4,f24.16)') i, j, '    ', arr_global(i,j)
             enddo
          enddo
       endif

       close(unit=fileunit)
       deallocate(arr_global)
    endif

  end subroutine write_array_to_file_real8_2d


  subroutine write_array_to_file_real8_3d(arr, fileunit, filename, parallel, write_binary, cycle_indices)

    ! Copy the input array into a global array and write all values to an output file.
    ! This can be useful for debugging, if we want to find differences between two fields
    ! (e.g., in two different runs).
    ! This version writes out 64-bit character strings corresponding to the binary representation
    ! of each floating-point variable. This can be useful for BFB comparisons.
    ! Sometimes, two floating-point variables appear to have the same values in base 10,
    ! when the last few bits actually vary.
    !TODO - Allow either float or binary output

    use glimmer_utils, only: double_to_binary
    use cism_parallel, only: gather_var

    real(dp), dimension(:,:,:), intent(in) :: arr    ! first two indices are i and j
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    type(parallel_type), intent(in) :: parallel
    logical, intent(in), optional :: write_binary
    logical, intent(in), optional :: cycle_indices   ! if true, then index 3->1, 1->2, 2->3

    integer :: i, j, k, kmax
    character(len=64) :: binary_str
    real(dp), dimension(:,:,:), allocatable :: arr_global
    real(dp), dimension(:,:,:), allocatable :: arr_cycle
    logical :: binary_output
    logical :: cycle_ind

    if (present(write_binary)) then
       binary_output = write_binary
    else
       binary_output = .false.
    endif

    if (present(cycle_indices)) then
       cycle_ind = cycle_indices
    else
       cycle_ind = .false.
    endif

    if (cycle_ind) then
       allocate(arr_cycle(size(arr,3), size(arr,1), size(arr,2)))
       kmax = size(arr,3)
       do j = 1, size(arr,2)
          do i = 1, size(arr,1)
             do k = 1, kmax
                arr_cycle(k,i,j) = arr(i,j,k)
             enddo
          enddo
       enddo
       call gather_var(arr_cycle, arr_global, parallel)
       deallocate(arr_cycle)
    else
       kmax = size(arr,1)
       call gather_var(arr, arr_global, parallel)
    endif

    if (main_task) then
       open(unit=fileunit, file=trim(filename), status='unknown')

       if (binary_output) then
          do j = 1, parallel%global_nsn
             do i = 1, parallel%global_ewn
                do k = 1, kmax
                   call double_to_binary(arr_global(k,i,j), binary_str)
                   write (fileunit, '(3i6,a4,a64)') i, j, k, '    ', binary_str
                enddo
             enddo
          enddo
       else
          do j = 1, parallel%global_nsn
             do i = 1, parallel%global_ewn
                do k = 1, kmax
                   write (fileunit, '(3i6,a4,f24.16)') i, j, k, '    ', arr_global(k,i,j)
                enddo
             enddo
          enddo
       endif

       close(unit=fileunit)
       deallocate(arr_global)
    endif

  end subroutine write_array_to_file_real8_3d

!=======================================================================

  subroutine glissade_cleanup_tiny_thickness(model, tiny_thck)

    ! Remove ice from cells with very small thicknesses.
    ! Add to the calving flux for now, but later put in the cleanup category

    use cism_parallel, only: parallel_halo

    type(glide_global_type), intent(inout) :: model   ! model instance
    real(dp), intent(in) :: tiny_thck    ! zero out where thck < tiny_thck

    integer :: nx, ny
    integer :: i, j

    type(parallel_type) :: parallel   ! info for parallel communication

    nx = model%general%ewn
    ny = model%general%nsn

    parallel = model%parallel

    ! Make sure the ice thickness is updated in halo cells
    call parallel_halo(model%geometry%thck, parallel)

    where (model%geometry%thck > 0.0d0 .and. model%geometry%thck < tiny_thck)
       model%calving%calving_thck = model%calving%calving_thck + model%geometry%thck
       model%geometry%thck = 0.0d0
    endwhere

  end subroutine glissade_cleanup_tiny_thickness

!=======================================================================

  subroutine glissade_cleanup_icefree_cells(model)

    ! Clean up prognostic variables in ice-free cells.
    ! This means seting most tracers to zero (or min(artm,0) for the case of temperature).

    use cism_parallel, only: parallel_halo

    type(glide_global_type), intent(inout) :: model   ! model instance

    integer :: nx, ny
    integer :: i, j

    type(parallel_type) :: parallel   ! info for parallel communication

    nx = model%general%ewn
    ny = model%general%nsn

    parallel = model%parallel

    ! Make sure the ice thickness is updated in halo cells
    call parallel_halo(model%geometry%thck, parallel)

    ! Set prognostic variables in ice-free columns to default values (usually zero).
    do j = 1, ny
       do i = 1, nx

          if (model%geometry%thck_old(i,j) > 0.0d0 .and. model%geometry%thck(i,j) == 0.0d0) then

             ! basal water
             model%basal_hydro%bwat(i,j) = 0.0d0

             ! thermal variables
             if (model%options%whichtemp == TEMP_INIT_ZERO) then
                model%temper%temp(:,i,j) = 0.0d0
             else
                model%temper%temp(:,i,j) = min(model%climate%artm(i,j), 0.0d0)
             endif

             if (model%options%whichtemp == TEMP_ENTHALPY) then
                model%temper%waterfrac(:,i,j) = 0.0d0
             endif

             ! other tracers
             ! Note: Tracers should be added here as they are added to the model

             if (model%options%whichcalving == CALVING_DAMAGE) then
                model%calving%damage(:,i,j) = 0.0d0
             endif

             if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
                model%geometry%ice_age(:,i,j) = 0.0d0
             endif

          endif    ! thck = 0

       enddo
    enddo

  end subroutine glissade_cleanup_icefree_cells

!****************************************************************************

!TODO - Other utility subroutines to add here?
!       E.g., calclsrf; subroutines to zero out tracers

!****************************************************************************

end module glissade_utils

!****************************************************************************
