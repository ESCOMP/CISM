!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_calving.F90 - part of the Community Ice Sheet Model (CISM)  
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
!!#ifdef HAVE_CONFIG_H
!!#include "config.inc"
!!#endif

module glissade_calving

  use glide_types
  use glimmer_global, only: dp
  use glimmer_log
  use parallel

  use glimmer_paramets, only: thk0

  implicit none

  private
  public :: glissade_calving_mask_init, glissade_thck_calving_threshold_init, &
            glissade_calve_ice, glissade_cull_calving_front, &
            glissade_remove_icebergs, glissade_limit_cliffs
  public :: verbose_calving

  logical, parameter :: verbose_calving = .false.

contains

!-------------------------------------------------------------------------------

  subroutine glissade_calving_mask_init(dx,                dy,               &
                                        thck,              topg,             &
                                        eus,               thklim,           &
                                        calving_front_x,   calving_front_y,  &
                                        calving_mask)

    ! Compute an integer calving mask if needed for the CALVING_GRID_MASK option

    use glissade_masks, only: glissade_get_masks

    ! Input/output arguments

    real(dp), intent(in) :: dx, dy                 !> cell dimensions in x and y directions (m)
    real(dp), dimension(:,:), intent(in) :: thck   !> ice thickness (m)
    real(dp), dimension(:,:), intent(in) :: topg   !> present bedrock topography (m)
    real(dp), intent(in) :: eus                    !> eustatic sea level (m)
    real(dp), intent(in) :: thklim                 !> minimum thickness for dynamically active grounded ice (m)
    real(dp), intent(in) :: calving_front_x        !> for CALVING_GRID_MASK option, calve ice wherever abs(x) > calving_front_x (m)
    real(dp), intent(in) :: calving_front_y        !> for CALVING_GRID_MASK option, calve ice wherever abs(y) > calving_front_y (m)

    integer, dimension(:,:), intent(inout) :: calving_mask   !> output mask: calve floating ice wherever the mask = 1

    ! Local variables

    real(dp) :: xcell, ycell     ! global cell center coordinates (m)
    integer :: nx, ny            ! horizontal grid dimensions
    integer :: i, j              ! local cell indices
    integer :: iglobal, jglobal  ! global cell indices

    integer, dimension(:,:), allocatable :: &
         ice_mask,             & ! = 1 where ice is present
         ocean_mask              ! = 1 for ice-free ocean

    real(dp) :: mask_maxval      ! maxval of calving_mask

    nx = size(calving_mask,1)
    ny = size(calving_mask,2)

    mask_maxval = maxval(calving_mask)
    mask_maxval = parallel_reduce_max(mask_maxval)

    ! Compute the calving mask, if not read in at initialization
 
    if (mask_maxval > 0) then

       ! calving_mask was read from the input file; do not need to compute a mask here

       if (verbose_calving .and. main_task) print*, 'Calving_mask was read from the input file'

    elseif (calving_front_x > 0.0d0 .or. calving_front_y > 0.0d0) then

       if (main_task) print*, 'Computing calving_mask based on calving_front_x/y'

       ! initialize
       calving_mask(:,:) = 0   ! no calving by default

       if (calving_front_x > 0.0d0) then

          ! set calving_mask = 1 where abs(x) > calving_front_x

          do j = 1, ny
             do i = 1, nx

                ! find global i and j indices
                call parallel_globalindex(i, j, iglobal, jglobal)

                ! find cell center x coordinate
                xcell = (dble(iglobal) - 0.5d0) * dx

                ! set calving mask = 1 based on cell coordinates relative to the calving front
                ! Note: Using absolute value to support symmetry with respect to x = 0
                if (abs(xcell) > calving_front_x) then
                   calving_mask(i,j) = 1
                endif

             enddo   ! i
          enddo   ! j

       endif   ! calving_front_x > 0

       if (calving_front_y > 0.0d0) then

          ! set calving_mask = 1 where abs(y) > calving_front_y

          do j = 1, ny
             do i = 1, nx

                ! find global i and j indices
                call parallel_globalindex(i, j, iglobal, jglobal)

                ! find cell center y coordinate
                ycell = (dble(jglobal) - 0.5d0) * dy

                ! set calving mask = 1 based on cell coordinates relative to the calving front
                if (abs(ycell) > calving_front_y) then
                   calving_mask(i,j) = 1
                endif

             enddo   ! i
          enddo   ! j

       endif   ! calving_front_y > 0

    else  ! compute the calving mask based on the initial ice extent
 
       if (main_task) then
          print*, 'Computing calving_mask based on initial ice extent'
       endif

       ! initialize
       calving_mask(:,:) = 0  ! no calving by default

       ! Get an ocean mask
       allocate(ice_mask(nx,ny))
       allocate(ocean_mask(nx,ny))

       call glissade_get_masks(nx,            ny,             &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               ocean_mask = ocean_mask)

       ! Set calving_mask = 1 for ice-free ocean cells.
       ! Any ice entering these cells during the run will calve.
       do j = 1, ny
          do i = 1, nx
             if (ocean_mask(i,j) == 1) then
                calving_mask(i,j) = 1
             endif
          enddo
       enddo

       deallocate(ice_mask)
       deallocate(ocean_mask)

    endif  ! mask_maxval > 0

    call parallel_halo(calving_mask)

  end subroutine glissade_calving_mask_init

!-------------------------------------------------------------------------------

  subroutine glissade_thck_calving_threshold_init(&
       nx,      ny,               &
       itest,   jtest,   rtest,   &
       which_ho_calving_front,    &
       thck,    topg,             &
       eus,     thklim,           &
       marine_connection_mask,    &
       calving_minthck,           &
       thck_calving_threshold)

    ! Given the input ice thickness, identify calving-front cells and compute the CF thickness.
    ! This subroutine is called at initialization when whichcalving = CALVING_THCK_THRESHOLD.
    !
    ! There are two options:
    ! (1) Set thck_calving_threshold to a uniform value of calving_minthck.
    ! (2) Set thck_calving_threshold to thck_calving_front, as computed from the input ice thickness.
    ! Note: This subroutine uses SI units.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask
    use glissade_grid_operators, only: glissade_scalar_extrapolate

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer, intent(in) :: nx, ny                              !> horizontal grid dimensions
    integer, intent(in) :: itest, jtest, rtest                 !> coordinates of diagnostic point
    integer, intent(in) :: which_ho_calving_front              !> = 1 for subgrid calving-front scheme, else = 0

    real(dp), dimension(nx,ny), intent(in) :: thck             !> ice thickness (m)
    real(dp), dimension(nx,ny), intent(in) :: topg             !> present bedrock topography (m)
    real(dp), intent(in)                   :: eus              !> eustatic sea level (m)
    real(dp), intent(in)                   :: thklim           !> minimum thickness for dynamically active grounded ice (m)
    integer,  dimension(nx,ny), intent(in) :: marine_connection_mask  !> = 1 for cells with a marine connection to the ocean

    ! Note: calving_minthck = 0 by default.
    !       If calving_minthck > 0 is specified in the config file, then we set thck_calving_threshold = calving_minthck.
    !       Otherwise, we set thck_calving_threshold to the input ice thickness (thck)

    real(dp), intent(in)                   :: calving_minthck  !> uniform thickness threshold for calving ice (m)
    real(dp), dimension(nx,ny), intent(out) ::  &
         thck_calving_threshold                 !> ice in the calving domain will calve when thinner than this value (m)

    ! Local variables

    integer :: i, j, n

    integer, dimension(nx,ny) :: &
         ice_mask,                  & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,             & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,                & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,                 & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask           ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny) :: &
         thck_calving_front,             & ! effective ice thickness at the calving front
         thck_calving_threshold_smoothed   ! smoothed version of thck_calving_threshold

    real(dp), parameter :: eps08 = 1.0d-8   ! small number

    !TODO - Make the calving thresholds config parameters?
    real(dp), parameter :: &
         calving_threshold_min = 100.d0,  &! min allowed value (m) of thck_calving_threshold
         calving_threshold_max = 500.d0    ! max allowed value (m) of thck_calving_threshold

    if (calving_minthck > eps08) then

       ! Set thck_calving_threshold to a uniform value
       thck_calving_threshold(:,:) = calving_minthck

       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'Set thck_calving_threshold to calving_minthck (m):', calving_minthck
       endif

    else

       ! Set thck_calving_threshold based on the input ice thickness and calving-front thickness

       ! Get masks

       call glissade_get_masks(nx,            ny,             &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               land_mask = land_mask)

       ! Note: If using a subgrid CF scheme, thck_calving_front is the effective ice thickness at the calving front,
       !        equal to the mean thickness of marine interior neighbors.
       !       Without a subgrid CF scheme, thck_calving_front is set to thck in calving_front cells.

       call glissade_calving_front_mask(nx,            ny,                 &
                                        which_ho_calving_front,            &
                                        thck,          topg,               &
                                        eus,                               &
                                        ice_mask,      floating_mask,      &
                                        ocean_mask,    land_mask,          &
                                        calving_front_mask,                &
                                        thck_calving_front)

       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'Calving front masks:'
          print*, ' '
          print*, 'thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'floating_mask, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') floating_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'calving_front_mask, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') calving_front_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thck_calving_front (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck_calving_front(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       if (which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then

          ! Extrapolate and smooth thck_calving_front, starting with CF cells where it was just computed

          if (verbose_calving .and. this_rank == rtest) then
             print*, ' '
             print*, 'Set thck_calving_threshold by extrapolating thck_calving_front to marine-connected cells'
          endif

          ! Limit thck_calving_front to lie within a prescribed range.
          where (calving_front_mask == 1)
             thck_calving_front = max(thck_calving_front, calving_threshold_min)
             thck_calving_front = min(thck_calving_front, calving_threshold_max)
          endwhere

          if (verbose_calving .and. this_rank == rtest) then
             print*, ' '
             print*, 'After corrections:'
             print*, ' '
             print*, 'thck_calving_front (m), itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') thck_calving_front(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

          ! Extrapolate the CF thickness from cells with calving_front_mask = 1
          !  to all cells with marine_connection_mask = 1.
          ! Apply a Laplacian smoother during the extrapolation.

          call glissade_scalar_extrapolate(nx,    ny,                 &
                                           calving_front_mask,        &
                                           thck_calving_front,        &
                                           marine_connection_mask,    &
                                           thck_calving_threshold,    &
                                           npoints_stencil = 9,       &
                                           apply_smoother = .true.,   &
                                           itest = itest, jtest = jtest, rtest = rtest)

          call parallel_halo(thck_calving_threshold)

          if (verbose_calving .and. this_rank == rtest) then
             print*, ' '
             print*, 'Extrapolated thck_calving_front to interior marine-based cells'
             print*, ' '
             print*, 'thck_calving_threshold (m), itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') thck_calving_threshold(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

       elseif (which_ho_calving_front == HO_CALVING_FRONT_NO_SUBGRID) then

          ! Set thck_calving_threshold = thck in marine-connected cells.
          ! As a result, any floating ice at the calving front, once thinned by mass-balance
          !  or ice-dynamics changes, will experience increased calving.

          where (marine_connection_mask == 1)
             thck_calving_threshold = thck
          elsewhere
             thck_calving_threshold = 0.0d0
          endwhere

          if (verbose_calving .and. this_rank == rtest) then
             print*, ' '
             print*, 'Set thck_calving_threshold = thck in marine-connected cells'
             print*, ' '
             print*, 'thck_calving_threshold (m), itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') thck_calving_threshold(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

       endif   ! which_ho_calving_front

    endif  ! calving_minthck > eps08

  end subroutine glissade_thck_calving_threshold_init

!-------------------------------------------------------------------------------

  subroutine glissade_calve_ice(nx,               ny,    &
                                which_calving,           &
                                calving_domain,          &
                                which_ho_calving_front,  &
                                calving,                 &  ! calving derived type
                                itest,   jtest,   rtest, &
                                dt,                      &  ! s
                                dx,               dy,    &  ! m
                                sigma,                   &
                                thklim,                  &  ! m
                                thck,             relx,  &  ! m
                                topg,             eus)      ! m

    ! Calve ice according to one of several methods.
    ! Note: This subroutine uses SI units.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask

    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer, intent(in) :: nx, ny                  !> horizontal grid dimensions

    integer, intent(in) :: which_calving           !> option for calving law
    integer, intent(in) :: calving_domain          !> option for where calving can occur
                                                   !> = 0 if calving occurs at the ocean edge only
                                                   !> = 1 if calving occurs everywhere the calving criterion is met
                                                   !> = 2 if calving occurs where criterion is met and there is a connected path
                                                   !>     to the ocean through other cells where the criterion is met
    integer, intent(in) :: which_ho_calving_front  !> = 1 for subgrid calving-front scheme, else = 0

    type(glide_calving), intent(inout) :: calving !> calving object

!    Note: The calving object includes the following fields and parameters used in this subroutine:
!    real(dp), intent(in)                     :: marine_limit        !> lower limit on topography elevation at marine edge before ice calves
                                                                     !> Note: marine_limit (shared by Glide) has scaled model units
!    real(dp), intent(in)                     :: calving_fraction    !> fraction of ice lost at marine edge when calving; 
                                                                     !> used with CALVING_FLOAT_FRACTION
!    real(dp), intent(in)                     :: timescale           !> timescale (s) for calving; calving_thck = thck * max(dt/timescale, 1)
                                                                     !> if timescale = 0, then calving_thck = thck
!    real(dp), intent(in)                     :: thck_calving_threshold  !> calve ice in the calving domain if thck < thck_calving_threshold (m);
                                                                         !> used with CALVING_THCK_THRESHOLD, EIGENCALVING and CALVING_DAMAGE
!    real(dp), intent(in)                     :: eigencalving_constant   !> eigencalving constant; m/s (lateral calving rate) per Pa (tensile stress)
!    real(dp), intent(in)                     :: eigen2_weight       !> weight given to tau_eigen2 relative to tau_eigen1 in tau_eff (unitless)
!    real(dp), dimension(:,:), intent(in)     :: tau_eigen1          !> first eigenvalue of 2D horizontal stress tensor (Pa)
!    real(dp), dimension(:,:), intent(in)     :: tau_eigen2          !> second eigenvalue of 2D horizontal stress tensor (Pa)
!    real(dp), dimension(:,:), intent(inout)  :: tau_eff             !> effective stress (Pa) for calving; derived from tau_eigen1/2
!    real(dp), dimension(:,:,:), intent(inout):: damage              !> 3D scalar damage parameter
!    real(dp), intent(in)                     :: damage_threshold    !> threshold value where ice is sufficiently damaged to calve
!    real(dp), intent(in)                     :: damage_constant     !> rate of change of damage (1/s) per unit stress (Pa)
!    real(dp) :: intent(in)                   :: lateral_rate_max    !> max lateral calving rate (m/s) for damaged ice
!    real(dp), dimension(:,:), intent(inout)  :: lateral_rate        !> lateral calving rate (m/s) at the calving front
                                                                     !> used with EIGENCALVING and CALVING_DAMAGE
!    integer,  dimension(:,:), intent(in)     :: calving_mask        !> integer mask: calve ice where calving_mask = 1
!    real(dp), dimension(:,:), intent(out)    :: calving_thck        !> thickness lost due to calving in each grid cell (m)

    integer, intent(in) :: itest, jtest, rtest                   !> coordinates of diagnostic point
    real(dp), intent(in)                      :: dt                !> model timestep (s)
    real(dp), intent(in)                      :: dx, dy            !> grid cell size in x and y directions (m)
    real(dp), dimension(:), intent(in)        :: sigma             !> vertical sigma coordinate
    real(dp), intent(in)                      :: thklim            !> minimum thickness for dynamically active grounded ice (m)
    real(dp), dimension(nx,ny), intent(inout) :: thck              !> ice thickness (m)
    real(dp), dimension(nx,ny), intent(in)    :: relx              !> relaxed bedrock topography (m)
    real(dp), dimension(nx,ny), intent(in)    :: topg              !> present bedrock topography (m)
    real(dp), intent(in)                      :: eus               !> eustatic sea level (m)

    ! local variables

    integer :: nz          ! number of vertical levels
                           ! Note: number of ice layers = nz-1
    integer :: i, j, k, n
    integer :: ii, jj

    real(dp), dimension(nx,ny) ::  &
         thck_calving_front,     & ! effective ice thickness at the calving front
         tau1, tau2,             & ! tau_eigen1 and tau_eigen2 (Pa), modified for calving
         damage_column             ! 2D vertically integrated scalar damage parameter

    ! basic masks
    integer, dimension(nx,ny)  ::  &
         ice_mask,               & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,          & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,             & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,              & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    ! Note: Calving occurs in a cell if and only if (1) the calving law permits calving, 
    !       and (2) the cell is in the calving domain, as specified by the calving_domain option.
    !       The calving domain by default is limited to the ocean edge (CALVING_DOMAIN_OCEAN_EDGE), 
    !       but can be extended to include all ice-covered cells (CALVING_DOMAIN_EVERYWHERE).

    !TODO - Make these integer masks like the ones above?
    logical, dimension(nx,ny) ::  &
         calving_law_mask,    & ! = T where the calving law permits calving, else = F
         calving_domain_mask    ! = T in the domain where calving is allowed to occur (e.g., at ocean edge), else = F

    real(dp) :: &
         float_fraction_calve, & ! = calving_fraction for which_calving = CALVING_FLOAT_FRACTION
                                 ! = 1.0 for which_calving = CALVING_FLOAT_ZERO
         thinning_rate,        & ! vertical thinning rate (m/s)
         calving_frac,         & ! fraction of potential calving that is actually applied
         upstream_lateral_rate,& ! lateral calving rate (m/s) applied to upstream cell
         frac_lateral,         & ! lateral_rate / lateral_rate_max 
         areafrac,             & ! fractional ice-covered area in a calving_front cell
         dthck,                & ! thickness change (m)
         d_damage_dt,          & ! rate of change of damage scalar (1/s)
         factor                  ! factor in quadratic formula

    character(len=100) :: message
   
    ! initialize

    nz = size(sigma)

    if (which_calving == CALVING_NONE) then   ! do nothing
       if (verbose_calving .and. main_task) print*, 'No calving'
       return
    endif

    !WHL - debug
    if (verbose_calving .and. main_task) then
       print*, ' '
       print*, 'In glissade_calve_ice, which_calving =', which_calving
       print*, 'calving_domain =', calving_domain
    endif

    ! Set the thickness fraction to be removed in each calving cell
    ! Note: The CALVING_FLOAT_FRACTION option has been superseded by the calving%timescale variable,
    !       but is included here for consistency with Glide.
    ! TODO: Remove CALVING_FLOAT_FRACTION option?

    if (which_calving == CALVING_FLOAT_FRACTION) then

       !WHL - Changed definition of calving fraction; now it is the fraction lost
       !      rather than the fraction remaining
       float_fraction_calve = calving%calving_fraction
       
    else  ! other calving options

       if (calving%timescale == 0.0d0) then  ! calve the entire column for eligible columns (this is the default)
          float_fraction_calve = 1.0d0
       else  ! calve a fraction of the column based on the calving time scale
          float_fraction_calve = min(dt/calving%timescale, 1.0d0)
       endif
       
    endif
       
    ! Do the calving based on the value of which_calving

    if (which_calving == EIGENCALVING .or. which_calving == CALVING_DAMAGE) then

       ! These two methods have several features in common:
       ! (1) The eigenvalues of the 2D horizontal stress tensor are key fields controlling the calving rate.
       ! (2) A lateral calving rate is computed in calving-front cells, then converted to a thinning rate.
       ! (3) The thinning rate is applied to CF cells and, if sufficiently large, to adjacent interior cells.
       !
       ! The main difference is that for eigencalving, the lateral calving rate is based on current stresses
       !  at the calving front, whereas for damage-based calving, the lateral calving rate is based on damage,
       !  which accumulates in floating cells due to stresses and then is advected downstream to the calving front.
       !
       ! At some point, we may want to prognose damage in a way that depends on other factors such as mass balance.

       ! Get masks
       ! Need a calving_front_mask; calving/thinning is applied only to cells at the calving front.
       ! Here, thck_calving_front is the effective thickness at the calving front, equal to
       !  the minimum thickness of a marine-based neighbor that is not on the calving front.
       ! Note: Cells with calving_front_mask = 1 are dynamically inactive unless thck >= thck_calving_front.
       !       For calving purposes, all calving_front cells are treated identically, whether or not
       !        dynamically active. Inactive cells receive eigenvalues by extrapolation from active cells.
       ! We pass in which_ho_calving_front = HO_CALVING_FRONT_SUBGRID so the subroutine will compute calving_front_mask
       !  and thck_calving_front, which are needed for calving options with a thickness threshold.
       ! The actual value for the run might be HO_CALVING_FRONT_NO_SUBGRID, with thck_calving_front needed
       !  only to compute the thck_calving_threshold field.

       call glissade_get_masks(nx,            ny,             &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               land_mask = land_mask)


       call glissade_calving_front_mask(nx,            ny,              &
                                        which_ho_calving_front,         &
                                        thck,          topg,            &
                                        eus,                            &
                                        ice_mask,      floating_mask,   &
                                        ocean_mask,    land_mask,       &
                                        calving_front_mask,             &
                                        thck_calving_front)

       !WHL - Debug
       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'floating_mask, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') floating_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'calving_front_mask, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') calving_front_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thck_calving_front (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck_calving_front(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif   ! verbose_calving

       ! For each floating cell, compute an effective stress based on eigenvalues of the stress tensor.
       ! Ignore negative eigenvalues corresponding to compressive stresses

       tau1 = max(calving%tau_eigen1, 0.0d0)
       tau2 = max(calving%tau_eigen2, 0.0d0)

       ! Ignore values on grounded ice
       where (floating_mask == 0)
          tau1 = 0.0d0
          tau2 = 0.0d0
       endwhere

       call parallel_halo(tau1)
       call parallel_halo(tau2)

       ! Compute the effective stress.
       ! Note: By setting eigen2_weight > 1, we can give greater weight to the second principle stress.
       !       This may be useful in calving unbuttressed shelves that are spreading in both directions.

       calving%tau_eff(:,:) = sqrt(tau1(:,:)**2 + (calving%eigen2_weight * tau2(:,:))**2)
       !WHL -  Uncomment the next line if assuming dependence on tau2 alone
!       calving%tau_eff(:,:) = calving%eigen2_weight * tau2(:,:)

       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'tau1 (Pa), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') tau1(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'tau2 (Pa), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') tau2(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'tau_eff (Pa), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') calving%tau_eff(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       ! Use the effective stress either to directly compute a lateral calving rate (for eigencalving),
       ! or to accumulate damage which is then used to derive a lateral calving rate (for damage-based calving).

       calving%lateral_rate(:,:) = 0.0d0

       if (which_calving == EIGENCALVING) then

          ! Compute the lateral calving rate (m/s) from the effective tensile stress in calving_front cells

          do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask(i,j) == 1) then
                   calving%lateral_rate(i,j) = calving%eigencalving_constant * calving%tau_eff(i,j)
                endif
             enddo   ! i
          enddo   ! j

       elseif (which_calving == CALVING_DAMAGE) then

          ! Prognose changes in damage.
          ! For now, this is done using a simple scheme based on the effective tensile stress, calving%tau_eff
          ! The damage is subsequently advected downstream.
          ! Note: The damage is formally a 3D field, which makes it easier to advect, even though
          !       (in the current scheme) the damage source term is uniform in each column.

          do j = 2, ny-1
             do i = 2, nx-1
                if (floating_mask(i,j) == 1) then
                   d_damage_dt = calving%damage_constant * calving%tau_eff(i,j)  ! damage_constant has units of s^{-1}/(Pa)
                   calving%damage(:,i,j) = calving%damage(:,i,j) + d_damage_dt * dt
                   calving%damage(:,i,j) = min(calving%damage(:,i,j), 1.0d0)
                   calving%damage(:,i,j) = max(calving%damage(:,i,j), 0.0d0)
                else  ! set damage to zero for grounded ice
                   calving%damage(:,i,j) = 0.0d0
                endif
             enddo
          enddo

          ! Compute the vertically integrated damage in each column.
          damage_column(:,:) = 0.0d0

          do j = 1, ny
             do i = 1, nx
                do k = 1, nz-1
                   damage_column(i,j) = damage_column(i,j) + calving%damage(k,i,j) * (sigma(k+1) - sigma(k))
                enddo
             enddo
          enddo

          ! Convert damage in CF cells to a lateral calving rate (m/s).
          ! Note: Although eigenprod = 0 in inactive calving-front cells, these cells can have significant damage
          !       advected from upstream, so in general we should not have to interpolate damage from upstream.
          !TODO - Verify this.
          do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask(i,j) == 1) then
                   frac_lateral = (damage_column(i,j) - calving%damage_threshold) / (1.0d0 - calving%damage_threshold)
                   frac_lateral = max(0.0d0, min(1.0d0, frac_lateral))
                   calving%lateral_rate(i,j) = calving%lateral_rate_max * frac_lateral  ! m/s
                endif
             enddo
          enddo

          if (verbose_calving .and. this_rank==rtest) then
             print*, ' '
             print*, 'damage increment, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.6)',advance='no') calving%damage_constant * calving%tau_eff(i,j) * dt
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'new damage, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.6)',advance='no') damage_column(i,j)
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
          endif
          
       endif   ! EIGENCALVING or CALVING_DAMAGE

       ! The following operations are shared by eigencalving and damage-based calving.

       call parallel_halo(calving%lateral_rate)

       ! Convert the lateral calving rate to a vertical thinning rate, conserving volume.
       ! Note: The calved volume is proportional to the effective shelf-edge thickness (thck_calving_front),
       !        not the nominal ice thickness (thck).
       !TODO: Change variable names? E.g., thinning_rate is really a volume loss rate.

       do j = 2, ny-1
          do i = 2, nx-1
             if (calving%lateral_rate(i,j) >  0.0d0) then

                thinning_rate = calving%lateral_rate(i,j) * thck_calving_front(i,j) / sqrt(dx*dy)  ! m/s
                dthck = thinning_rate * dt  ! m

                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   print*, ' '
                   print*, 'Calving: r, i, j =', rtest, itest, jtest
                   print*, 'dx (m), dt (yr) =', sqrt(dx*dy), dt/scyr
                   print*, 'lateral calving rate (m/yr) =', calving%lateral_rate(i,j)*scyr
                   print*, 'dthck (m) =', dthck
                endif

                ! Compute the new ice thickness
                if (dthck > thck(i,j)) then
                   calving%calving_thck(i,j) = calving%calving_thck(i,j) + thck(i,j)
                   thck(i,j) = 0.0d0
                else   ! dthck <= thck
                   thck(i,j) = thck(i,j) - dthck
                   calving%calving_thck(i,j) = calving%calving_thck(i,j) + dthck
                endif

             endif   ! calving%lateral_rate > 0
          enddo   ! i
       enddo   ! j

       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'Finished eigencalving or damage-based calving, task =', this_rank
          print*, ' '
          print*, 'lateral calving rate (m/yr), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') calving%lateral_rate(i,j) * scyr
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'calving_thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') calving%calving_thck(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'new thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    endif  ! eigencalving or damage-based calving


    if (which_calving == CALVING_THCK_THRESHOLD .or. which_calving == EIGENCALVING    &
                                                .or. which_calving == CALVING_DAMAGE) then

       ! Note: Eigencalving or damage-based calving, if done above, is followed by thickness-based calving.
       !       This helps get rid of thin ice near the CF where stress eigenvalues might be small.

       ! Get masks
       ! For eigencalving, masks were computed above, but should be recomputed before doing more calving

       call glissade_get_masks(nx,            ny,             &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               land_mask = land_mask)

       call glissade_calving_front_mask(&
                               nx,            ny,                 &
                               which_ho_calving_front,            &
                               thck,          topg,               &
                               eus,                               &
                               ice_mask,      floating_mask,      &
                               ocean_mask,    land_mask,          &
                               calving_front_mask,                &
                               thck_calving_front)

       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'Thickness-based calving:'
          print*, ' '
          print*, 'calving_front_mask, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') calving_front_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thck_calving_front (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck_calving_front(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thck_calving_threshold (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') calving%thck_calving_threshold(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       ! Apply thinning in calving-front cells whose effective thickness (H_e = thck_calving_front) is less than
       !  a prescribed minimum value (Hc_min = thck_calving_threshold).
       !
       ! The effective thinning rate is given by
       !
       !    dH_e/dt = -(Hc_min - H_e) / tau_c  where Hc_min > H_e
       !    dH_e/dt = 0 elsewhere
       !
       ! where tau_c = calving%timescale.
       !
       ! The thinning rate applied to the mean cell thickness (thck) is given by
       !
       !    dH/dt = min(H/H_e, 1) * dH_e/dt
       !
       ! Thus, any ice with H_e < Hc_min is removed on a time scale given by tau_c.

       if (calving%timescale <= 0.0d0) then
          write(message,*) 'Must set calving timescale to a positive nonzero value for this calving option'
          call write_log(message, GM_FATAL)
       endif

       do j = 2, ny-1
          do i = 2, nx-1
             if (calving_front_mask(i,j) == 1 .and. &
                  thck_calving_front(i,j) > 0.0d0 .and. thck_calving_front(i,j) <= calving%thck_calving_threshold(i,j)) then
                
!!                if (verbose_calving .and. thck(i,j) > 0.0d0) &
!!                     print*, 'Calve thin floating ice: task, i, j, thck =', this_rank, i, j, thck(i,j)

                ! calving%timescale has units of s
                thinning_rate = (calving%thck_calving_threshold(i,j) - thck_calving_front(i,j)) / calving%timescale
                !WHL - Do not weight by areafrac
!!                areafrac = min(thck(i,j)/thck_calving_front(i,j), 1.0d0)
!!                dthck = areafrac*thinning_rate * dt
                dthck = thinning_rate * dt

                !WHL - debug
                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   print*, ' '
                   print*, 'Thinning: r, i, j =', rtest, itest, jtest
                   print*, 'thck:', thck(i,j)
                   print*, 'thck_calving_front (m) =', thck_calving_front(i,j)
                   print*, 'thck_calving_threshold (m) =', calving%thck_calving_threshold(i,j)
                   print*, 'thinning rate (m/yr) =', thinning_rate * scyr
                   print*, 'dthck (m) =', dthck
                endif

                if (dthck > thck(i,j)) then
                   calving%calving_thck(i,j) = calving%calving_thck(i,j) + thck(i,j)
                   thck(i,j) = 0.0d0
                else
                   thck(i,j) = thck(i,j) - dthck
                   calving%calving_thck(i,j) = calving%calving_thck(i,j) + dthck
                endif

             endif   ! thck_calving_front < thck_calving_threshold in calving_front cell
          enddo   ! i
       enddo   ! j

       !WHL - debug
       if (verbose_calving .and. this_rank == rtest) then

          print*, ' '
          print*, 'Did thickness-based calving, task =', this_rank 
         print*, ' '
          print*, 'new thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif  ! verbose

    else   ! other calving options

       ! Get masks.
       ! Use thickness limit of 0.0 instead of thklim so as to remove ice from any cell
       !  that meets the calving criteria, not just dynamically active ice.

       call glissade_get_masks(nx,            ny,             &
                               thck,          topg,           &
                               eus,           0.0d0,          &   ! thklim = 0.0
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask)

       ! set the calving-law mask
       ! Note: Cells that meet the calving-law criteria will be calved provided they also lie in the calving domain,
       !       as determined below.

       select case (which_calving)

       case(CALVING_FLOAT_ZERO, CALVING_FLOAT_FRACTION)     ! calve ice that is floating

          do j = 1, ny
             do i = 1, nx
                if (floating_mask(i,j) == 1) then
                   calving_law_mask(i,j) = .true.
                else
                   calving_law_mask(i,j) = .false.
                endif
             enddo
          enddo

          !NOTE: The Glide version of CALVING_FLOAT_ZERO calves all floating ice.
          !      Glissade calves floating ice only in the calving domain, which is CALVING_DOMAIN_OCEAN_EDGE by default.
          !      Must set calving_domain = CALVING_DOMAIN_EVERYWHERE to match the Glide behavior.
          !TODO: Change the default to calving_domain_everywhere?

       case(CALVING_RELX_THRESHOLD)   ! set thickness to zero if relaxed bedrock is below a given level

          !WHL - The Glide version of CALVING_RELX_THRESHOLD calves ice wherever the relaxed bedrock criterion is met.
          !      Must set calving_domain = CALVING_DOMAIN_EVERYWHERE to match the Glide behavior.
          ! Note: calving%marine_limit (a holdover from Glide) has scaled model units
          where (relx <= calving%marine_limit*thk0 + eus)   ! convert marine_limit from scaled units to m
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_TOPG_THRESHOLD)   ! set thickness to zero if present bedrock is below a given level

          where (topg < calving%marine_limit*thk0 + eus)    ! convert marine_limit from scaled units to m
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_HUYBRECHTS)    ! Huybrechts grounding line scheme for Greenland initialization

          if (eus > -80.d0) then
             where (relx <= 2.d0*eus)
                calving_law_mask = .true.
             elsewhere
                calving_law_mask = .false.
             end where
          elseif (eus <= -80.d0) then
             where (relx <= (2.d0*eus - 0.25d0*(eus + 80.d0)**2.d0))
                calving_law_mask = .true.
             elsewhere
                calving_law_mask = .false.
             end where
          end if

       end select

       ! halo update (may not be necessary if thck, damage, etc. are correct in halos, but including to be safe)
       call parallel_halo(calving_law_mask)

       ! set the calving domain mask

       if (calving_domain == CALVING_DOMAIN_OCEAN_EDGE) then  ! calving domain includes floating cells at margin only
                                                              !WHL - Could modify to include grounded marine cells at margin
          do j = 2, ny-1
             do i = 2, nx-1

                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   print*, 'task, i, j, ice_mask, floating_mask:',  &
                        this_rank, i, j, ice_mask(i,j), floating_mask(i,j)
                endif

                if ( floating_mask(i,j) == 1 .and.   &
                     (ocean_mask(i-1,j)==1 .or. ocean_mask(i+1,j)==1 .or. ocean_mask(i,j-1)==1 .or. ocean_mask(i,j+1)==1) ) then
                   calving_domain_mask(i,j) = .true.
                else
                   calving_domain_mask(i,j) = .false.
                endif
             enddo
          enddo

          ! halo update (since the loop above misses some halo cells)
          call parallel_halo(calving_domain_mask)

          if (verbose_calving .and. this_rank==rtest) then
             print*, ' '
             print*, 'calving_domain_mask, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(L10)',advance='no') calving_domain_mask(i,j)
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
          endif

       elseif (calving_domain == CALVING_DOMAIN_EVERYWHERE) then  ! calving domain includes all cells

          calving_domain_mask(:,:) = .true.

       endif   ! calving_domain

       ! Calve ice where calving_law_mask = T and calving_domain_mask = T
       do j = 1, ny
          do i = 1, nx
             if (calving_law_mask(i,j) .and. calving_domain_mask(i,j)) then

                if (verbose_calving .and. this_rank==rtest .and. thck(i,j) > 0.0d0) then
!!                   print*, 'Calve ice: task, i, j, calving_thck =', this_rank, i, j, float_fraction_calve * thck(i,j)
                endif

                calving%calving_thck(i,j) = calving%calving_thck(i,j) + float_fraction_calve * thck(i,j)
                thck(i,j) = thck(i,j) - float_fraction_calve * thck(i,j)
            endif
          enddo
       enddo

    endif   ! which_calving

    if (verbose_calving .and. this_rank==rtest) then
       print*, ' '
       print*, 'calving_thck, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') calving%calving_thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'After calving, new thck:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
    endif

  end subroutine glissade_calve_ice

!---------------------------------------------------------------------------

  subroutine glissade_cull_calving_front(&
       nx,           ny,          &
       itest, jtest, rtest,       &
       thck,         topg,        &
       eus,          thklim,      &
       which_ho_calving_front,    &
       ncull_calving_front,       &
       calving_thck)

    ! Optionally, remove all cells currently at the calving front are removed.
    ! This subroutine would typically be called at initialization, if at all.
    ! Culling can removed long, skinny floating peninsulas that can be dynamically unstable.
    !  Without this step, peninsulas up to two cells thick (with calving-front cells on each side)
    !  will typically be removed as icebergs (because there is no path back to grounded ice through active cells).
    ! With one round of culling, peninsulas up to four cells thick will be removed (two outer layers
    !  during the preliminary step, followed by two inner layers on the remove_iceberg step).
    ! If necessary, culling can be repeated to remove peninsulas with a thickness of 6 layers, 8 layers, etc.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask

    integer, intent(in) :: nx, ny                       !> horizontal grid dimensions
    integer, intent(in) :: itest, jtest, rtest          !> coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(inout) :: thck     !> ice thickness
    real(dp), dimension(nx,ny), intent(in)    :: topg     !> present bedrock topography
    real(dp), intent(in)    :: eus                      !> eustatic sea level
    real(dp), intent(in)    :: thklim                   !> minimum thickness for dynamically active grounded ice
    integer, intent(in)     :: which_ho_calving_front   !> = 1 for subgrid calving-front scheme, else = 0
    integer, intent(in) :: &
         ncull_calving_front           !> number of times to cull calving_front cells at initialization

    real(dp), dimension(nx,ny), intent(inout) :: calving_thck   !> thickness lost due to calving in each grid cell;
                                                              !> on output, includes ice that is culled here
    ! local variables

    integer :: i, j, n

    integer,  dimension(nx,ny) ::  &
         ice_mask,           & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,      & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,         & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,          & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask    ! = 1 where ice is floating and borders the ocean, else = 0

    real(dp),  dimension(nx,ny) ::  &
         thck_calving_front    ! effective ice thickness at the calving front

    do n = 1, ncull_calving_front

       ! calculate masks
       ! Note: Passing in thklim = 0.0 does not work because it erroneously counts thin floating cells as active.
       !       Then the algorithm can fail to identify floating regions that should be removed
       !       (since they are separated from any active cells).

       call glissade_get_masks(nx,            ny,             &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               land_mask = land_mask)

       call glissade_calving_front_mask(nx,            ny,                 &
                                        which_ho_calving_front,            &
                                        thck,          topg,               &
                                        eus,                               &
                                        ice_mask,      floating_mask,      &
                                        ocean_mask,    land_mask,          &
                                        calving_front_mask,                &
                                        thck_calving_front)

       if (main_task) then
          call write_log ('cull_calving_front: Removing ice from calving_front cells')
          print*, 'cull_calving_front: Removing ice from calving_front cells'
       endif

       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'calving_front_mask for culling, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') calving_front_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'cull_calving_front: Before removing CF cells, n =', n
          print*, ' '
          print*, 'thck:'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       do j = 1, ny
          do i = 1, nx
             if (calving_front_mask(i,j) == 1) then
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
             endif
          enddo
       enddo

       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'cull_calving_front: After removing CF cells, n =', n
          print*, ' '
          print*, 'thck:'
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    enddo  ! ncull_calving_front

  end subroutine glissade_cull_calving_front

!---------------------------------------------------------------------------

  subroutine glissade_remove_icebergs(&
       nx,           ny,            &
       itest, jtest, rtest,         &
       thck,                        &
       f_ground_cell,               &
       ice_mask,                    &
       land_mask,                   &
       active_ice_mask,             &
       calving_thck)

    ! Remove any icebergs. 
        
    ! The algorithm is as follows:
    ! (1) Mark all cells with ice (either active or inactive) with the initial color.
    !     Mark other cells with the boundary color.
    ! (2) Seed the fill by giving active grounded cells the fill color.
    ! (3) Recursively fill all cells that are connected to filled cells by a path
    !     that passes through active cells only.
    ! (4) Repeat the recursion as necessary to spread the fill to adjacent processors.
    ! (5) Once the fill is done, any floating cells that still have the initial color
    !     are considered to be icebergs and are removed.
    !
    ! Notes:
    ! (1) Grounded cells must have f_ground_cell > f_ground_threshold to seed the fill.
    ! (2) The recursive fill applies to edge neighbors, not corner neighbors.
    !     The path back to grounded ice must go through edges, not corners.
    ! (3) Inactive cells can be filled (if adjacent to active cells), but
    !     do not further spread the fill.
    ! (4) Should have thklim > 0.  With a limit of 0.0, very thin floating cells
    !     can be wrongly counted as active, and icebergs can be missed.
    ! (5) Land-based cells that still have the initial color are not marked as icebergs.

    use glissade_masks, only: glissade_fill_with_buffer, initial_color, fill_color, boundary_color

    integer :: nx, ny                                   !> horizontal grid dimensions
    integer, intent(in) :: itest, jtest, rtest          !> coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(inout) :: thck            !> ice thickness
    real(dp), dimension(nx,ny), intent(in)    :: f_ground_cell   !> grounded fraction in each grid cell
    integer,  dimension(nx,ny), intent(in)    :: ice_mask        !> = 1 where ice is present (thck > thklim), else = 0
    integer,  dimension(nx,ny), intent(in)    :: land_mask       !> = 1 where topg - eus >= 0, else = 0
    integer,  dimension(nx,ny), intent(in)    :: active_ice_mask !> = 1 for dynamically active cells
    real(dp), dimension(nx,ny), intent(inout) :: calving_thck    !> thickness lost due to calving in each grid cell;
                                                                 !> on output, includes ice in icebergs
    ! local variables

    integer :: i, j, iter

    integer :: &
         max_iter,             & ! max(ewtasks, nstasks)
         local_count,          & ! local counter for filled values
         global_count,         & ! global counter for filled values
         global_count_save       ! globalcounter for filled values from previous iteration

    integer,  dimension(nx,ny) ::  &
         color                 ! integer 'color' for identifying icebergs

    real(dp),  dimension(nx,ny) ::  &
         thck_calving_front    ! effective ice thickness at the calving front

    real(dp), parameter :: &   ! threshold for counting cells as grounded
         f_ground_threshold = 0.10d0

    if (verbose_calving .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_remove_icebergs'
       print*, ' '
       print*, 'thck, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'active_ice_mask:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') active_ice_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'f_ground_cell:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') f_ground_cell(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif
    
    ! Initialize iceberg removal
    ! Note: Any cell with ice, active or inactive, receives the initial color.
    !       Inactive cells can later receive the fill color (if adjacent to active cells)
    !        but cannot further spread the fill color.
    !       This protects inactive calving-front cells from removal, as desired.

    do j = 1, ny
       do i = 1, nx
          if (thck(i,j) > 0.0d0) then
             color(i,j) = initial_color
          else
             color(i,j) = boundary_color
          endif
       enddo
    enddo

    ! Loop through cells, identifying active cells with grounded ice.
    ! Fill each grounded cell and then recursively fill active neighbor cells, whether grounded or not.
    ! We may have to do this several times to incorporate connections between neighboring processors.

    max_iter = max(ewtasks,nstasks)
    global_count_save = 0

    do iter = 1, max_iter

       if (iter == 1) then   ! identify grounded cells that can seed the fill

          do j = 1, ny
             do i = 1, nx

                ! Active cells that are firmly grounded can seed the fill.
                ! It is somewhat arbitrary how to define "firmly grounded", but here we require
                !  that f_ground_cell exceeds a threshold value defined above.
                ! Note: If running without a GLP, then f_ground_cell is binary, either 0 or 1.

                if (active_ice_mask(i,j) == 1 .and. f_ground_cell(i,j) > f_ground_threshold) then  ! grounded ice

                   if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then

                      ! assign the fill color to this cell, and recursively fill neighbor cells
                      call glissade_fill_with_buffer(nx,    ny,    &
                                                     i,     j,     &
                                                     color, active_ice_mask)

                   endif

                endif
             enddo
          enddo

       else  ! count > 1

          ! Check for halo cells that were just filled on neighbor processors
          ! Note: In order for a halo cell to seed the fill on this processor, it must not only have the fill color,
          !       but also must be an active cell.

          call parallel_halo(color)

          ! west halo layer
          i = nhalo
          do j = 1, ny
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(nx,    ny,    &
                                               i+1,   j,     &
                                               color, active_ice_mask)
             endif
          enddo

          ! east halo layers
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(nx,    ny,    &
                                               i-1,   j,     &
                                               color, active_ice_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(nx,    ny,    &
                                               i,     j+1,   &
                                               color, active_ice_mask)
             endif
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(nx,    ny,    &
                                               i,     j-1,   &
                                               color, active_ice_mask)
             endif
          enddo

       endif  ! count = 1

       local_count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) local_count = local_count + 1
          enddo
       enddo

       !WHL - If running a large problem, may want to reduce the frequency of this global sum
       global_count = parallel_reduce_sum(local_count)

       if (global_count == global_count_save) then
          if (verbose_calving .and. main_task) &
               print*, 'Fill converged: iter, global_count =', iter, global_count
          exit
       else
          if (verbose_calving .and. main_task) &
               print*, 'Convergence check: iter, global_count =', iter, global_count
          global_count_save = global_count
       endif

    enddo  ! count

    ! Icebergs are cells that still have the initial color and are not on land.
    ! Remove ice in these cells, adding it to the calving field.
    ! Note: There is an exception for cells that are
    !       (1) adjacent to at least one ice-covered cell (sharing an edge), and
    !       (2) connected diagonally to an active cell with the fill color.
    !       Such cells are considered part of the inactive calving front and are
    !        allowed to continue filling instead of calving.

    do j = 2, ny-1
       do i = 2, nx-1
          if (color(i,j) == initial_color .and. land_mask(i,j) == 0) then
             if (  ( color(i-1,j+1)==fill_color .and. active_ice_mask(i-1,j+1)==1 .and. &
                       (ice_mask(i-1,j)==1 .or. ice_mask(i,j+1)==1) ) &
              .or. ( color(i+1,j+1)==fill_color .and. active_ice_mask(i+1,j+1)==1 .and. &
                       (ice_mask(i+1,j)==1 .or. ice_mask(i,j+1)==1) ) &
              .or. ( color(i-1,j-1)==fill_color .and. active_ice_mask(i-1,j-1)==1 .and. &
                       (ice_mask(i-1,j)==1 .or. ice_mask(i,j-1)==1) ) &
              .or. ( color(i+1,j-1)==fill_color .and. active_ice_mask(i+1,j-1)==1 .and. &
                       (ice_mask(i+1,j)==1 .or. ice_mask(i,j-1)==1) ) ) then
                ! do nothing; this cell is part of the inactive calving front
             else  ! not part of the inactive calving front; calve as an iceberg
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
                !TODO - Also handle tracers?  E.g., set damage(:,i,j) = 0.d0?
             endif  ! diagonally connected or not
          endif
       enddo
    enddo

    if (verbose_calving .and. this_rank == rtest) then
       print*, ' '
       print*, 'Done in glissade_remove_icebergs'
       print*, ' '
       print*, 'thck, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_remove_icebergs

!---------------------------------------------------------------------------

  subroutine glissade_limit_cliffs(&
       nx,             ny,              &
       itest,  jtest,  rtest,           &
       dt,                              &
       which_ho_calving_front,          &
       taumax_cliff,   cliff_timescale, &
       thck,           topg,            &
       eus,            thklim,          &
       calving_thck)

    ! Impose a thickness limit on marine ice cliffs.
    ! These are defined as grounded marine-based cells adjacent to inactive calving_front cells or ice-free ocean.
    ! Ice removed from cliffs is added to the calving flux.

    use glissade_masks

    integer, intent(in)  :: nx, ny                      !> horizontal grid dimensions
    integer, intent(in)  :: itest, jtest, rtest         !> coordinates of diagnostic point
    real(dp), intent(in) :: dt                          !> model timestep (s)

    integer, intent(in)     :: which_ho_calving_front   !> = 1 for subgrid calving-front scheme, else = 0
    real(dp), intent(in)    :: taumax_cliff             !> yield stress (Pa) for marine-based ice cliffs
    real(dp), intent(in)    :: cliff_timescale          !> timescale (s) for limiting marine cliff thickness
    real(dp), dimension(nx,ny), intent(inout) :: thck   !> ice thickness (m)
    real(dp), dimension(nx,ny), intent(in)    :: topg   !> present bedrock topography (m)
    real(dp), intent(in)    :: eus                      !> eustatic sea level (m)
    real(dp), intent(in)    :: thklim                   !> minimum thickness for dynamically active grounded ice (m)

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_thck             !> thickness (m) lost due to calving; on output, includes ice calved at marine cliffs

    ! local variables

    integer :: i, j

    integer, dimension(nx,ny) ::  &
         ice_mask,           & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,      & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,         & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,          & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask, & ! = 1 where ice is floating and borders the ocean, else = 0
         active_ice_mask,    & ! = 1 for dynamically active cells
         marine_cliff_mask     ! = 1 where ice is grounded and marine-based and borders at least one
                               !     ocean or inactive calving_front cell, else = 0

    real(dp), dimension(nx,ny) ::  &
         thck_calving_front    ! effective ice thickness at the calving front

    real(dp) :: &
         thinning_rate,        & ! vertical thinning rate (m/s)
         dthck,                & ! thickness change (m)
         thckmax_cliff,        & ! max stable ice thickness in marine_cliff cells
         factor

    ! Update masks, including the marine_cliff mask.
    ! Note: We do not use calving_front_mask or thck_calving_front directly.
    !       But to identify marine cliffs, we use active_ice_mask, which depends on whether there is a subgrid calving front.

    call glissade_get_masks(nx,            ny,             &
                            thck,          topg,           &
                            eus,           thklim,         &
                            ice_mask,                      &
                            floating_mask = floating_mask, &
                            ocean_mask = ocean_mask,       &
                            land_mask = land_mask)

    call glissade_calving_front_mask(nx,            ny,           &
                                     which_ho_calving_front,      &
                                     thck,          topg,               &
                                     eus,                               &
                                     ice_mask,      floating_mask,      &
                                     ocean_mask,    land_mask,          &
                                     calving_front_mask,                &
                                     thck_calving_front,                &
                                     active_ice_mask = active_ice_mask)

    call glissade_marine_cliff_mask(nx,            ny,                &
                                    ice_mask,      floating_mask,     &
                                    land_mask,     active_ice_mask,   &
                                    marine_cliff_mask)

    if (verbose_calving .and. this_rank==rtest) then
       print*, ' '
       print*, 'In glissade_limit_cliffs'
       print*, ' '
       print*, 'marine_cliff_mask, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') marine_cliff_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'thckmax_cliff, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             factor = taumax_cliff / (rhoi*grav)   ! units are Pa for taumax, m for factor
             thckmax_cliff = factor + sqrt(factor**2 + (rhoo/rhoi)*(topg(i,j))**2)  ! m
             write(6,'(f10.3)',advance='no') thckmax_cliff
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
    endif

    do j = 2, ny-1
       do i = 1, nx-1
          if (marine_cliff_mask(i,j) == 1) then

             ! Compute the max stable ice thickness in the cliff cell.
             ! This is eq. 2.10 in Bassis & Walker (2012)
             factor = taumax_cliff / (rhoi*grav)   ! units are Pa for taumax, m for factor
             thckmax_cliff = factor + sqrt(factor**2 + (rhoo/rhoi)*(topg(i,j))**2)  ! m

             !WHL - debug
             if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, ' '
                print*, 'Cliff thinning: r, i, j =', rtest, itest, jtest
                print*, 'thck, thckmax_cliff (m) =', thck(i,j), thckmax_cliff
             endif

             ! If thicker than the max stable thickness, then remove some ice and add it to the calving field
             ! Note: By default, cliff_timescale = 0, which means thck is reset to thckmax_cliff each timestep.
             !       Might want to try other values when looking at marine ice cliff instability.
             if (thck(i,j) > thckmax_cliff) then

                if (cliff_timescale > 0.0d0) then
                   thinning_rate = (thck(i,j) - thckmax_cliff) / cliff_timescale
                   dthck = min(thck(i,j) - thckmax_cliff, thinning_rate*dt)
                else
                   dthck = thck(i,j) - thckmax_cliff
                endif

                !WHL - debug
                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   print*, 'thinning rate (m/yr) =', thinning_rate * scyr
                   print*, 'dthck (m) =', dthck
                endif

                thck(i,j) = thck(i,j) - dthck
                calving_thck(i,j) = calving_thck(i,j) + dthck

             endif  ! thck > thckmax_cliff

          endif  ! marine_cliff cell
       enddo   ! i
    enddo   ! j

  end subroutine glissade_limit_cliffs

!---------------------------------------------------------------------------

end module glissade_calving

!---------------------------------------------------------------------------
