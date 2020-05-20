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
  use glimmer_log
  use glide_types

  implicit none

  private
  public :: glissade_adjust_thickness, glissade_smooth_topography, glissade_adjust_topography
  public :: glissade_stdev, verbose_stdev

  logical, parameter :: verbose_stdev = .true.

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

    use glimmer_paramets, only: thk0
    use parallel, only: this_rank, parallel_reduce_max

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

    logical, parameter :: verbose_adjust_thickness = .false.

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

    ! Make sure ursf was read in with nonzero values.
    ! Otherwise, we cannot use usrf to adjust the ice thickness.

    usrf_max = maxval(model%geometry%usrf)
    usrf_max = parallel_reduce_max(usrf_max)

    if (usrf_max > tiny(0.0d0)) then

       if (verbose_adjust_thickness .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'adjust thck: itest, jtest, rank =', itest, jtest, rtest
          print*, ' '
          print*, 'Before thck adjustment, usrf (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%usrf(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'Before thck adjustment, thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'Before thck adjustment, topg (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%topg(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'Before thck adjustment, cavity thickness (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') (model%geometry%usrf(i,j) - model%geometry%thck(i,j)  &
                     - model%geometry%topg(i,j)) * thk0
             enddo
             write(6,*) ' '
          enddo
       endif   ! verbose

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

       if (verbose_adjust_thickness .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, ' '
          print*, 'adjust thck: itest, jtest, rank =', itest, jtest, rtest
          print*, ' '
          print*, 'After thck adjustment, usrf (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%usrf(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'After thck adjustment, thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'After thck adjustment, topg (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') model%geometry%topg(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'After thck adjustment, cavity thickness (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') (model%geometry%usrf(i,j) - model%geometry%thck(i,j)  &
                     - model%geometry%topg(i,j)) * thk0
             enddo
             write(6,*) ' '
          enddo
       endif   ! verbose

    else   ! usrf_max < tiny

       call write_log('Error: Must read in usrf to use adjust_input_thickness option', GM_FATAL)

    endif   ! usrf_max > tiny

  end subroutine glissade_adjust_thickness

!****************************************************************************

  subroutine glissade_smooth_topography(model)

    ! Use a Laplacian smoother to smooth the input bed topography
    !TODO - This smoothing needs some more testing.  In particular, it is unclear how best to treat
    !        the ice thickness in regions that transition from grounded to floating
    !        when the topography is smoothed. Is it better to preserve thickness, or to
    !        increase thickness to keep the ice grounded?

    use glimmer_paramets, only: thk0
    use parallel, only: this_rank
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
    endif

    ! compute the initial upper surface elevation (to be held fixed under smoothing of bed topography)
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    ! compute initial mask
    call glissade_get_masks(nx,                  ny,                    &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   0.0d0,                 &  ! thklim = 0
                            ice_mask,                                   &
                            floating_mask = floating_mask)

    if (verbose_smooth_topg .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'itest, jtest, rank =', itest, jtest, rtest
       print*, ' '
       print*, 'Before Laplacian smoother, topg (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%topg(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Before Laplacian smoother, usrf (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%usrf(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Before Laplacian smoother, thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

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

    if (verbose_smooth_topg .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'itest, jtest, rank =', itest, jtest, rtest
       print*, ' '
       print*, 'After Laplacian smoother, topg (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%topg(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'After Laplacian smoother, usrf (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%usrf(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'After Laplacian smoother, thck (m):'
       do j = jtest+3, jtest-3, -1
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_smooth_topography

!****************************************************************************

  subroutine glissade_adjust_topography(model)

    ! Adjust the input bed topography in a specified region.
    ! For example, we may want to raise the topography close to the surface in a region
    !  where the ice is not sufficiently grounded, and the data are not well constrained.
    ! Note: So far, this subroutine has been used to raise eastern Thwaites topography.
    !       It has not been used to lower topography.

    use glimmer_paramets, only: thk0
    use parallel, only: this_rank
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

    real(dp) :: factor
    real(dp) :: xmin, xmax, ymin, ymax
    real(dp) :: topg_lo, topg_hi, topg_delta

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
    endif

    xmin = model%paramets%adjust_topg_xmin
    xmax = model%paramets%adjust_topg_xmax
    ymin = model%paramets%adjust_topg_ymin
    ymax = model%paramets%adjust_topg_ymax
    topg_lo = model%paramets%adjust_topg_lo
    topg_hi = model%paramets%adjust_topg_hi
    topg_delta = model%paramets%adjust_topg_delta

    if (verbose_adjust_topg .and. this_rank == rtest) then
       i = itest
       j = jtest
       print*, ' '
       print*, 'Adjust input topography, diag point: r, i ,j =', rtest, itest, jtest
       print*, 'x1, y1 =', model%general%x1(i), model%general%y1(j)
       print*, 'thck, topg =', model%geometry%thck(i,j)*thk0, model%geometry%topg(i,j)*thk0
       print*, 'xmin, xmax =', xmin, xmax
       print*, 'ymin, ymax =', ymin, ymax
       print*, 'topg_lo, topg_hi (m) =', topg_lo, topg_hi
       print*, 'topg_delta =', topg_delta
    endif

    ! Compute the lower and upper ice surface before the adjustment
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    if (verbose_adjust_topg .and. this_rank == rtest) then
       print*, ' '
       print*, 'Input usrf (m):'
       print*, ' '
       do j = jtest+3, jtest-3, -1
          write(6,'(a10)',advance='no') '          '
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%usrf(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Input thck (m):'
       print*, ' '
       do j = jtest+3, jtest-3, -1
          write(6,'(a10)',advance='no') '          '
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Input lsrf (m):'
       print*, ' '
       do j = jtest+3, jtest-3, -1
          write(6,'(a10)',advance='no') '          '
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%lsrf(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

    if (verbose_adjust_topg .and. this_rank == rtest) then
       print*, ' '
       print*, 'Input topography (m):'
       print*, ' '
       write(6,'(a10)',advance='no') '  y1 \ x1 '
       do i = itest-3, itest+3
          write(6,'(f10.0)',advance='no') model%general%x1(i)
       enddo
       print*, ' '
       do j = jtest+3, jtest-3, -1
          write(6,'(f10.0)',advance='no') model%general%y1(j)
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%topg(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

    !TODO - Use model%geometry%topg - model%climate%eus?
    allocate(topg(model%general%ewn, model%general%nsn))
    topg = model%geometry%topg * thk0

    ! Apply the topographic correction.
    ! Where topg > topg_hi, apply the max correction, topg_delta.
    ! Where topg < topg_lo, apply no correction.
    ! Where topg_lo < topg < topg_hi, phase in the correction linearly.

    do j = 1, ny
       do i = 1, nx
          if (model%general%x1(i) >= xmin .and. model%general%x1(i) <= xmax .and. &
              model%general%y1(j) >= ymin .and. model%general%y1(j) <= ymax) then
             if (topg(i,j) > topg_lo) then
                factor = min((topg(i,j) - topg_lo)/(topg_hi - topg_lo), 1.0d0)
                topg(i,j) = topg(i,j) + factor * topg_delta
             endif
          endif
       enddo
    enddo

    model%geometry%topg = topg / thk0
    deallocate(topg)

    if (verbose_adjust_topg .and. this_rank == rtest) then
       print*, ' '
       print*, 'New topography (m):'
       print*, ' '
       write(6,'(a10)',advance='no') '  y1 \ x1 '
       do i = itest-3, itest+3
          write(6,'(f10.0)',advance='no') model%general%x1(i)
       enddo
       print*, ' '
       do j = jtest+3, jtest-3, -1
          write(6,'(f10.0)',advance='no') model%general%y1(j)
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%topg(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

    ! In some cells, the new lower surface (usrf - thck) may lie below the topography.
    ! In these cells, reduce the ice thickness such that lsrf = topg, preserving the input value of usrf.

    where (model%geometry%usrf - model%geometry%thck < model%geometry%topg)
       model%geometry%thck = model%geometry%usrf - model%geometry%topg
    endwhere

    if (verbose_adjust_topg .and. this_rank == rtest) then
       print*, ' '
       print*, 'Corrected thck (m):'
       print*, ' '
       do j = jtest+3, jtest-3, -1
          write(6,'(a10)',advance='no') '          '
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') model%geometry%thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_adjust_topography

!****************************************************************************

  subroutine glissade_stdev(&
        nx,        ny,             &
        grid_ratio,                &
        idiag,     jdiag,          &
        phi,       stdev)

    use parallel, only: main_task

    integer, intent(in) :: nx, ny                    ! number of cells in x and y direction on input grid (global)
    integer, intent(in) :: grid_ratio                ! nx/nx_coarse = ny/ny_coarse = ratio of fine to coarse grid
    integer, intent(in) :: idiag, jdiag              ! global coordinates of diagnostic cell
    real(dp), dimension(nx,ny), intent(in) :: phi    ! input field
    real(dp), dimension(nx,ny), intent(out) :: stdev ! standard deviation of input field, on the coarse grid

    integer :: i, j, ic, jc, ilo, ihi, jlo, jhi

    integer :: nx_coarse, ny_coarse      ! dimensions of coarse grid where we compute st.dev

    real(dp) :: &
         sumx,   &
         sumx2,  &
         xav,    &
         x2av,   &
         stdev_coarse

    ! Check that nx_coarse and ny_coarse divide evenly into nx and ny with the same quotient

    if (verbose_stdev .and. main_task) then
       print*, 'In verbose_stdev'
       print*, 'Fine nx, ny =', nx, ny
       print*, 'grid_ratio =', grid_ratio
    endif

    if (mod(nx,grid_ratio) == 0) then
       nx_coarse = nx/grid_ratio
    else
       if (main_task) print*, 'stdev error, nx/grid_ratio is not an integer'
       stop
    endif

    if (mod(ny,grid_ratio) == 0) then
       ny_coarse = ny/grid_ratio
       if (main_task) print*, 'nx_coarse, ny_coarse =', nx_coarse, ny_coarse
    else
       if (main_task) print*, 'stdev error, ny/grid_ratio is not an integer'
       stop
    endif

    stdev(:,:) = 0.0d0

    ! Loop over the coarse grid
    do jc = 1, ny_coarse
       do ic = 1, nx_coarse

          ! Compute the standard deviation of phi within one cell on the coarse grid

          ilo = (ic-1)*grid_ratio + 1
          ihi = ic * grid_ratio
          jlo = (jc-1)*grid_ratio + 1
          jhi = jc * grid_ratio

          sumx = 0.0d0
          sumx2 = 0.0d0
          do j = jlo, jhi
             do i = ilo, ihi
                sumx = sumx + phi(i,j)
                sumx2 = sumx2 + phi(i,j)**2
             enddo
          enddo
          xav = sumx / grid_ratio**2
          x2av = sumx2 / grid_ratio**2
          stdev_coarse = sqrt(x2av - xav*xav)

          ! Write this value to each cell on the fine grid.
          ! This is done because CISM doesn't have an easy way to work with two grids.
          !TODO: Read in two grids?
          stdev(ilo:ihi,jlo:jhi) = stdev_coarse

       enddo
    enddo

    if (verbose_stdev .and. main_task) then
       print*, ' '
       print*, 'Input field, idiag, jdiag:', idiag, jdiag
       print*, ' '
       do j = jdiag+3, jdiag-3, -1
          write(6,'(a10)',advance='no') '          '
          do i = idiag-3, idiag+3
             write(6,'(f10.3)',advance='no') phi(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, 'stdev:'
       print*, ' '
       do j = jdiag+3, jdiag-3, -1
          write(6,'(a10)',advance='no') '          '
          do i = idiag-3, idiag+3
             write(6,'(f10.3)',advance='no') stdev(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_stdev


!TODO - Other utility subroutines to add here?
!       E.g., tridiag; calclsrf; subroutines to zero out tracers

!****************************************************************************

end module glissade_utils

!****************************************************************************
