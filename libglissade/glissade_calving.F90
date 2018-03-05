!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_calving.F90 - part of the Community Ice Sheet Model (CISM)  
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
!!#ifdef HAVE_CONFIG_H
!!#include "config.inc"
!!#endif

module glissade_calving

  use glide_types
  use glimmer_global, only: dp
  use glimmer_log
  use parallel

  use glimmer_paramets, only: thk0, tim0

  implicit none

  ! colors for fill subroutine
  integer, parameter :: initial_color = 0   ! initial color, represented by integer
  integer, parameter :: fill_color = 1      ! fill color, represented by integer
  integer, parameter :: boundary_color = -1 ! boundary color, represented by integer

  !WHL - debug
    logical, parameter :: verbose_calving = .false.
!!    logical, parameter :: verbose_calving = .true.
  

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

  subroutine glissade_calve_ice(which_calving,           &
                                calving_domain,          &
                                which_ho_calving_front,  &
                                remove_icebergs,         &
                                limit_marine_cliffs,     &
                                itest,   jtest,   rtest, &
                                thck,             relx,  &
                                topg,             eus,   &
                                thklim,                  &
                                marine_limit,            &
                                calving_fraction,        &
                                calving_timescale,       &
                                dt,                      &
                                dx,               dy,    &
                                eigenprod,               &
                                eigencalving_constant,   &
                                calving_minthck,         &
                                taumax_cliff,            &
                                cliff_timescale,         &
                                calving_mask,            &
                                damage,                  &
                                damage_threshold,        &
                                damage_column,           &
                                sigma,                   &
                                calving_thck,            &
                                cull_calving_front_in,   &
                                ncull_calving_front_in)

    ! Calve ice according to one of several methods

    use glissade_masks

    implicit none

    !TODO - Convert input/output arguments to SI units
    !TODO - Shorten the argument list by passing in a calving derived type
    !---------------------------------------------------------------------
    ! Subroutine arguments
    ! Currently, thck, relx, topg, eus, marine_limit, calving_minthck and calving_thck are scaled by thk0
    ! TODO - Remove thickness scaling from this subroutine.  Would need to multiply several input arguments by thk0. 
    !---------------------------------------------------------------------

    integer, intent(in) :: which_calving          !> option for calving law
    integer, intent(in) :: calving_domain         !> option for where calving can occur
                                                  !> = 0 if calving occurs at the ocean edge only
                                                  !> = 1 if calving occurs everywhere the calving criterion is met
                                                  !> = 2 if calving occurs where criterion is met and there is a connected path
                                                  !>     to the ocean through other cells where the criterion is met
    integer, intent(in) :: which_ho_calving_front !> = 1 for subgrid calving-front scheme, else = 0
    logical, intent(in) :: remove_icebergs        !> if true, then remove icebergs after calving
    logical, intent(in) :: limit_marine_cliffs    !> if true, then limit the thickness of marine-based ice cliffs
    integer, intent(in) :: itest, jtest, rtest    !> coordinates of diagnostic point

    real(dp), dimension(:,:), intent(inout) :: thck              !> ice thickness
    real(dp), dimension(:,:), intent(in)    :: relx              !> relaxed bedrock topography
    real(dp), dimension(:,:), intent(in)    :: topg              !> present bedrock topography
    real(dp), intent(in)                    :: eus               !> eustatic sea level
    real(dp), intent(in)                    :: thklim            !> minimum thickness for dynamically active grounded ice
    real(dp), intent(in)                    :: marine_limit      !> lower limit on topography elevation at marine edge before ice calves
    real(dp), intent(in)                    :: calving_fraction  !> fraction of ice lost at marine edge when calving; 
                                                                 !> used with which_calving = CALVING_FLOAT_FRACTION
    real(dp), intent(in)                    :: calving_timescale !> time scale for calving; calving_thck = thck * max(dt/calving_timescale, 1)
                                                                 !> if calving_timescale = 0, then calving_thck = thck
    real(dp), intent(in)                    :: dt                !> model timestep (used with calving_timescale)
    real(dp), intent(in)                    :: dx, dy            !> grid cell size in x and y directions (m)
    real(dp), dimension(:,:), intent(inout) :: eigenprod         !> determinant (product of eigenvalues) of horizontal strain rate tensor (s^{-2})
    real(dp), intent(in)                    :: eigencalving_constant   !> dimensionless constant used for eigencalving
    real(dp), intent(in)                    :: calving_minthck   !> min thickness of floating ice before it calves;
                                                                 !> used with CALVING_THCK_THRESHOLD and EIGENCALVING
    real(dp), intent(in)                    :: taumax_cliff      !> yield stress (Pa) for marine-based ice cliffs
                                                                 !> used with limit_marine_cliffs option
    real(dp), intent(in)                    :: cliff_timescale   !> time scale for limiting marine cliff thickness
    integer,  dimension(:,:), intent(in)    :: calving_mask      !> integer mask: calve ice where calving_mask = 1
    real(dp), dimension(:,:,:), intent(in)  :: damage            !> 3D scalar damage parameter
!    real(dp), dimension(:,:,:), intent(inout)  :: damage         !> 3D scalar damage parameter
                                                                  !WHL - 'inout' if damage is updated below
    real(dp), dimension(:,:), intent(out)   :: damage_column     !> 2D vertically integrated scalar damage parameter
    real(dp), intent(in)                    :: damage_threshold  !> threshold value where ice is sufficiently damaged to calve
    real(dp), dimension(:), intent(in)      :: sigma             !> vertical sigma coordinate
    real(dp), dimension(:,:), intent(out)   :: calving_thck      !> thickness lost due to calving in each grid cell

    logical, intent(in), optional :: &
         cull_calving_front_in         !> if true, then cull calving_front cells to improve model stability
                                       !> generally should be called only at model initialization

    integer, intent(in), optional :: &
         ncull_calving_front_in        !> number of times to cull calving_front cells at initialization

    ! local variables

    logical :: cull_calving_front     ! local version of cull_calving_front_in
    integer :: ncull_calving_front    ! local version of ncull_calving_front_in

    integer :: nx, ny      ! horizontal grid dimensions
    integer :: nz          ! number of vertical levels
                           ! Note: number of ice layers = nz-1
    integer :: i, j, k, n
    integer :: ii, jj
    integer :: count, maxcount_fill  ! loop counters

    real(dp), dimension(:,:), allocatable ::  &
         thck_calving_front,     & ! effective ice thickness at the calving front
         thck_init                 ! value of thck before calving

    real(dp), dimension(:,:), allocatable ::  &
         calving_thck_init         ! debug diagnostic only

    integer,  dimension(:,:), allocatable ::  &
         color                ! integer 'color' for filling the calving domain (with CALVING_DOMAIN_OCEAN_CONNECT)

    ! basic masks
    integer, dimension(:,:), allocatable   ::  &
         ice_mask,               & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,          & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,             & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,              & ! = 1 where topg is at or above sea level, else = 0
         active_ice_mask,        & ! = 1 for cells that are dynamically active, else = 0 
         calving_front_mask,     & ! = 1 where ice is floating and borders at least one ocean cell, else = 0
         marine_cliff_mask         ! = 1 where ice is grounded and marine-based and borders at least
                                   !  one ocean or calving_front cell, else = 0

    ! Note: Calving occurs in a cell if and only if (1) the calving law permits calving, 
    !       and (2) the cell is in the calving domain, as specified by the calving_domain option.
    !       The calving domain by default is limited to the ocean edge (CALVING_DOMAIN_OCEAN_EDGE), 
    !       but can be extended to include all ice-covered cells (CALVING_DOMAIN_EVERYWHERE), or
    !       cells connected to the ocean through other cells that meet the calving criterion
    !       (CALVING_DOMAIN_OCEAN_CONNECT).
    ! TODO: Change the default to calving_domain_everywhere?

    !TODO - Make these integer masks like the ones above?
    logical, dimension(:,:), allocatable   ::  &
         calving_law_mask,    & ! = T where the calving law permits calving, else = F
         calving_domain_mask    ! = T in the domain where calving is allowed to occur (e.g., at ocean edge), else = F

    real(dp) :: &
         float_fraction_calve  ! = calving_fraction for which_calving = CALVING_FLOAT_FRACTION
                               ! = 1.0 for which_calving = CALVING_FLOAT_ZERO

    real(dp) :: &
         calving_rate,       & ! lateral calving rate (m/yr)
         thinning_rate,      & ! vertical thinning rate (m/yr, converted to scaled model units)
         calving_frac,       & ! fraction of potential calving that is actually applied
         areafrac,           & ! fractional ice-covered area in a calving_front cell
         dthck,              & ! thickness change (model units)
         thckmax_cliff,      & ! max stable ice thickness in marine_cliff cells
         factor                ! factor in quadratic formula

    real(dp), parameter :: &
         thinning_limit = 0.99d0  ! When ice not originally on the calving front is allowed to thin,
                                  ! the resulting thickness must be at least thinning_limit*thck_init

    character(len=100) :: message
   
    !WHL - debug
    integer :: sum_fill_local, sum_fill_global  ! number of filled cells

    integer :: iplot1, iplot2

    !default
!    iplot1 = nx-20
!    iplot2 = nx-1

    ! initialize

    if (present(cull_calving_front_in)) then
       cull_calving_front = cull_calving_front_in
    else
       cull_calving_front = .false.
    endif

    if (present(ncull_calving_front_in)) then
       ncull_calving_front = ncull_calving_front_in
    else
       ncull_calving_front = 1  ! only used if cull_calving_front = T
    endif

    calving_thck(:,:) = 0.d0

    nx = size(thck,1)
    ny = size(thck,2)
    nz = size(sigma)


    if (which_calving == CALVING_NONE) then

       ! remove icebergs (if desired) and return

       if (remove_icebergs) then

          if (verbose_calving .and. main_task) then
             print*, 'Remove icebergs'
          endif

          ! Set the thickness thresholds for cells to be counted as active.
          ! Floating regions will be identified as icebergs unless they are connected to grounded ice
          !  along a path consisting of active cells.

          call glissade_remove_icebergs(&
               itest,   jtest,   rtest, &
               thck,                    &
               topg,          eus,      &
               thklim,                  &
               which_ho_calving_front,  &
               calving_thck,            &
               cull_calving_front,      &
               ncull_calving_front)

       endif

       return

    endif

    ! allocate masks
    ! Not all of these are needed for all calving options, but it is simplest just to allocate them all
    allocate (calving_law_mask(nx,ny))
    allocate (calving_domain_mask(nx,ny))
    allocate (ice_mask(nx,ny))
    allocate (floating_mask(nx,ny))
    allocate (ocean_mask(nx,ny))
    allocate (land_mask(nx,ny))
    allocate (active_ice_mask(nx,ny))
    allocate (calving_front_mask(nx,ny))
    allocate (thck_calving_front(nx,ny))
    allocate (thck_init(nx,ny))
    allocate (marine_cliff_mask(nx,ny))

    !WHL - debug
    allocate(calving_thck_init(nx,ny))
    calving_thck_init(:,:) = thck(:,:)

    !WHL - debug
    if (verbose_calving .and. main_task) then
       print*, ' '
       print*, 'In glissade_calve_ice'
       print*, 'which_calving =', which_calving
       print*, 'calving_domain =', calving_domain
!       print*, 'i, relx, topg, thck, usfc:'
!       do i = iplot1, iplot2
!          print*, i, relx(i,j)*thk0, topg(i,j)*thk0, thck(i,j)*thk0, (topg(i,j) + thck(i,j))*thk0
!       enddo
    endif

    ! Set the thickness fraction to be removed in each calving cell
    ! Note: The CALVING_FLOAT_FRACTION option has been superseded by the calving_timescale variable,
    !       but is included here for consistency with Glide.
    ! TODO: Remove CALVING_FLOAT_FRACTION option?

    if (which_calving == CALVING_FLOAT_FRACTION) then

       !WHL - Changed definition of calving fraction; now it is the fraction lost
       !      rather than the fraction remaining
       float_fraction_calve = calving_fraction
       
    else  ! other calving options

       if (calving_timescale == 0.0d0) then  ! calve the entire column for eligible columns (this is the default)
          float_fraction_calve = 1.0d0
       else  ! calve a fraction of the column based on the calving time scale
          float_fraction_calve = min(dt/calving_timescale, 1.0d0)
       endif
       
    endif
       
    ! Do the calving based on the value of which_calving

    if (which_calving == EIGENCALVING) then

       ! Save the initial thickness, which is used below to identify upstream interior cells.
       thck_init(:,:) = thck(:,:)

       ! Get masks
       ! Need a calving_front_mask; calving/thinning is applied only to cells at the calving front.
       ! Here, thck_calving_front is the effective thickness at the calving front, equal to
       !  the minimum thickness of a marine-based neighbor that is not on the calving front.
       ! Note: Cells with calving_front_mask = 1 are dynamically inactive unless thck >= thck_calving_front.
       !       For calving purposes, all calving_front cells are treated identically, whether or not
       !        dynamically active. Inactive cells receive an eigenproduct by extrapolation
       !        from active cells.


       call glissade_get_masks(nx,            ny,             &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               calving_front_mask = calving_front_mask, &
                               thck_calving_front = thck_calving_front)

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
                write(6,'(f10.3)',advance='no') thck_calving_front(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          
       endif  ! verbose

       ! Loop over cells, computing calving for cells on the calving front.
       ! The calving rate is proportional to the eigenproduct in the cell, following Levermann et al. (2012).
       ! If the cell has zero eigenproduct (likely as a result of being dynamically inactive in the last velocity solve),
       !  then use the eigenproduct (if positive) from the upstream cell that determines thck_calving_front for this cell.

       ! Zero out the eigenproduct in grid cells that are not floating.
       ! The idea is that spreading in grounded cells is not an accurate approximation
       !  of spreading in adjacent floating cells.

       do j = 1, ny
          do i = 1, nx
             if (floating_mask(i,j) == 0) then
                eigenprod(i,j) = 0.0d0
             endif
          enddo
       enddo

       ! Extend eigenprod from active floating cells to adjacent (usually inactive) calving_front cells.
       do j = 2, ny-1
          do i = 2, nx-1
             if (calving_front_mask(i,j) == 1) then

                ! if eigenprod = 0 here, then look for positive eigenprod upstream
                if (eigenprod(i,j) == 0) then

                   ! loop over edge neighbors and choose the maximum spreading rate
                   do jj = j-1, j+1
                      do ii = i-1, i+1
                         if (ii == i .or. jj == j) then  ! edge neighbors only
                            if (thck_init(ii,jj) > 0.0d0 .and. calving_front_mask(ii,jj) == 0) then
                               eigenprod(i,j) = max(eigenprod(i,j), eigenprod(ii,jj))
                            endif
                         endif
                      enddo   ! ii
                   enddo   !! jj

                   if (eigenprod(i,j) == 0) then

                      ! If no edge neighors have positive spreading, then loop over corner neighbors.
                      ! This will put spreading values in calving_front cells at the corners of fronts.
                      do jj = j-1, j+1
                         do ii = i-1, i+1
                            if (abs(ii-i) == 1 .and. abs(jj-j) == 1) then ! corner neighbors
                               if (thck_init(ii,jj) > 0.0d0 .and. calving_front_mask(ii,jj) == 0) then
                                  eigenprod(i,j) = max(eigenprod(i,j), eigenprod(ii,jj))
                               endif
                            endif
                         enddo   ! ii
                      enddo   ! jj

                   endif   ! eigenprod = 0 after looping over edge neighbors

                endif  ! eigenprod = 0

                ! Do eigencalving in this calving_front cell

                if (eigenprod(i,j) > 0.0d0 .and. thck_calving_front(i,j) > 0.0d0) then

                   ! compute the calving rate as in Levermann et al. (2012)
                   ! Units: The eigencalving_constant K has units of m*yr.
                   !        The determinant has units of s^(-2); convert to yr^(-2).
                   !        Thus the calving rate has units of m/yr in the lateral direction.
                   !        Both dx and dy have units of m. Usually dx = dy, but take sqrt(dx*dy) for generality.
                   calving_rate = eigencalving_constant * eigenprod(i,j)*scyr*scyr

                   ! Convert the lateral calving rate to a vertical thinning rate, conserving volume.
                   ! Note: The calved volume is proportional to the effective shelf-edge thickness (thck_calving_front),
                   !       not the nominal ice thickness (thck).
                   thinning_rate = calving_rate * thck_calving_front(i,j)*thk0 / sqrt(dx*dy)  ! m/yr
                   dthck = thinning_rate * (tim0/scyr)/thk0  * dt  ! convert to model units

                   if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                      print*, ' '
                      print*, 'Eigencalving: r, i, j =', rtest, itest, jtest
                      print*, 'dx (m), dt (yr) =', sqrt(dx*dy), dt*tim0/scyr
                      print*, 'spreading rate (yr-2) =', eigenprod(i,j)*scyr*scyr
                      print*, 'calving rate (m/yr) =', calving_rate
                      print*, 'dthck (m) =', thinning_rate * dt*tim0/scyr
                   endif

                   ! Compute the new ice thickness
                   ! If the column calves completely, then determine how much potential calving still remains.

                   if (dthck > thck(i,j)) then
                      calving_frac = thck(i,j)/dthck
                      calving_rate = calving_rate * (1.0d0 - calving_frac)
                      calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                      thck(i,j) = 0.0d0

                      ! Apply some calving to the upstream cell with thck_init = thck_calving_front(i,j).
                      ! However, do not allow the upstream cell to end up thinner than thinning_limit*thck_init,
                      !  where thinning_limit is a scalar slightly less than 1.
                      ! The thinned upstream cell will usually be inactive during the next velocity solve.

                      do jj = j-1, j+1
                         do ii = i-1, i+1
                            if (thck_init(ii,jj) > 0.0d0 .and. thck_init(ii,jj) == thck_calving_front(i,j)) then
                               thinning_rate = calving_rate * thck_calving_front(i,j)*thk0 / sqrt(dx*dy)  ! m/yr
                               dthck = thinning_rate * (tim0/scyr)/thk0  * dt  ! convert to model units
                               dthck = min(dthck, thck(ii,jj) - thinning_limit*thck_init(ii,jj))
                               thck(ii,jj) = thck(ii,jj) - dthck
                               calving_thck(ii,jj) = calving_thck(ii,jj) + dthck
                            endif
                         enddo   ! ii
                      enddo   ! jj

                   else
                      thck(i,j) = thck(i,j) - dthck
                      calving_thck(i,j) = calving_thck(i,j) + dthck
                   endif

                endif   ! eigenprod > 0 and thck_calving_front > 0

             endif   ! calving_front_mask = 1
          enddo   ! i
       enddo   ! j

       if (verbose_calving .and. this_rank == rtest) then
          print*, ' '
          print*, 'eigenprod (yr-2), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(e10.2)',advance='no') eigenprod(i,j) * scyr*scyr
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'Finished eigencalving, task =', this_rank
          print*, ' '
          print*, 'calving_thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') calving_thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'new thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
       endif

    endif  ! eigencalving

    ! Note: Both the CALVING_THCK_THRESHOLD and EIGENCALVING options incorporate thickness-based calving.
    !       Eigencalving can be considered a special case in which thickness-based calving is preceded by
    !        thinning in cells with a positive strain rate determinant.

    if (which_calving == CALVING_THCK_THRESHOLD .or. which_calving == EIGENCALVING) then

       !WHL - debug
       calving_thck_init(:,:) = calving_thck(:,:)

       ! Save the initial thickness, which is used below to identify upstream interior cells.
       thck_init(:,:) = thck(:,:)

       ! Get masks
       ! For eigencalving, masks were computed above, but should be recomputed before doing more calving

       call glissade_get_masks(nx,            ny,             &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               calving_front_mask = calving_front_mask, &
                               thck_calving_front = thck_calving_front)

       !WHL - debug
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
                write(6,'(f10.3)',advance='no') thck_calving_front(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo

       endif

       ! Apply thinning in calving-front cells whose effective thickness H_e (thck_calving_front) is less than
       !  a prescribed minimum value Hc_min (calving_minthck).
       !
       ! The effective thinning rate is given by
       !
       !    dH_e/dt = -(Hc_min - H_e) / tau_c  where Hc_min > H_e
       !    dH_e/dt = 0 elsewhere
       !
       ! where tau_c = calving_timescale.
       !
       ! The thinning rate applied to the mean cell thickness (thck) is given by
       !
       !    dH/dt = max(H/H_e, 1) * dH_e/dt
       !
       ! Thus, any ice with H_e < Hc_min is removed on a time scale given by tau_c.

       if (calving_timescale <= 0.0d0) then
          write(message,*) 'Must set calving_timescale to a positive nonzero value for this calving option'
          call write_log(message, GM_FATAL)
       endif

       do j = 2, ny-1
          do i = 2, nx-1
             if (calving_front_mask(i,j) == 1 .and. &
                  thck_calving_front(i,j) > 0.0d0 .and. thck_calving_front(i,j) <= calving_minthck) then
                
!!                if (verbose_calving .and. thck(i,j) > 0.0d0) &
!!                     print*, 'Calve thin floating ice: task, i, j, thck =', this_rank, i, j, thck(i,j)*thk0

                ! Note: calving_minthck, thck_calving_front and calving_timescale have scaled model units
                thinning_rate = (calving_minthck - thck_calving_front(i,j)) / calving_timescale
                areafrac = min(thck(i,j)/thck_calving_front(i,j), 1.0d0)
                dthck = areafrac*thinning_rate*dt

                !WHL - debug
                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   print*, ' '
                   print*, 'Thinning: r, i, j =', rtest, itest, jtest
                   print*, 'thck:', thck(i,j)*thk0
                   print*, 'thck_calving_front (m) =', thck_calving_front(i,j)*thk0
                   print*, 'calving_minthck (m) =', calving_minthck*thk0
                   print*, 'areafrac =', areafrac
                   print*, 'thinning rate (m/yr) =', thinning_rate * thk0*scyr/tim0
                   print*, 'dt (yr) ', dt * tim0/scyr
                   print*, 'dthck (m) =', dthck * thk0
                endif

                if (dthck > thck(i,j)) then
                   calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                   thck(i,j) = 0.0d0

                   ! Apply a little bit of thinning to the upstream cell with thck_init = thck_calving_front(i,j)
                   !  (if the upstream cell is floating).
                   ! However, do not allow the upstream cell to end up thinner than thinning_limit*thck_init,
                   !  where thinning_limit is a scalar slightly less than 1.
                   ! The thinned upstream cell will usually be inactive during the next velocity solve.

                   do jj = j-1, j+1
                      do ii = i-1, i+1
                         if (thck_init(ii,jj) > 0.0d0 .and. thck_init(ii,jj) == thck_calving_front(i,j) &
                              .and. floating_mask(ii,jj) == 1) then
                            thinning_rate = (calving_minthck - thck_calving_front(i,j)) / calving_timescale
                            dthck = max(thinning_rate*dt, 0.0d0)
                            dthck = min(dthck, thck(ii,jj) - thinning_limit*thck_init(ii,jj))
                            thck(ii,jj) = thck(ii,jj) - dthck
                            calving_thck(ii,jj) = calving_thck(ii,jj) + dthck
                         endif
                      enddo   ! ii
                   enddo   ! jj

                else
                   thck(i,j) = thck(i,j) - dthck
                   calving_thck(i,j) = calving_thck(i,j) + dthck
                endif

             endif   ! thck_calving_front < calving_minthck in calving_front cell
          enddo   ! i
       enddo   ! j

       !WHL - debug
       if (verbose_calving .and. this_rank == rtest) then

          print*, ' '
          print*, 'Did thickness-based calving, task =', this_rank
          print*, ' '
          print*, 'Thickness-based calving_thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
!                write(6,'(f10.3)',advance='no') calving_thck(i,j)*thk0
                write(6,'(f10.3)',advance='no') (calving_thck(i,j) - calving_thck_init(i,j))*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'new thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo

       endif  ! verbose

       ! Clean up the calving front following eigencalving and/or thickness-based calving.
       ! Note: A cell at an advancing calving front can have thck > thck_calving_front.
       !  Such a cell will be active during the next velocity calculation. 
       !  In order to prevent unrealistic thicknesses at the calving front 
       !   (i.e., calving_front cells thicker than adjacent interior cells), we reset 
       !   the thickness to thck_calving_front (just enough to make the cell active) 
       !   and add the excess ice to the calving_thck field.

       calving_thck_init(:,:) = calving_thck(:,:)

       do j = 2, ny-1
          do i = 2, nx-1
             if (calving_front_mask(i,j) == 1 .and. &
                  thck_calving_front(i,j) > 0.0d0 .and. thck(i,j) > thck_calving_front(i,j)) then
                dthck = thck(i,j) - thck_calving_front(i,j)
                calving_thck(i,j) = calving_thck(i,j) + dthck
                thck(i,j) = thck_calving_front(i,j)
             endif
          enddo   ! j
       enddo   ! i

       if (verbose_calving .and. this_rank == rtest) then

          print*, ' '
          print*, 'Limited ice thickness at the calving front, task =', this_rank

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

          print*, ' '
          print*, 'thck_calving_front (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck_calving_front(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'calving_thck (m) from limiting, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') (calving_thck(i,j) - calving_thck_init(i,j))*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'new thck (m), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          
       endif

    elseif (which_calving == CALVING_GRID_MASK) then

       ! calve ice where the input calving mask = 1

       if (verbose_calving .and. this_rank==rtest) then
          print*, ' '
          print*, 'Limit advance of calving front'
          print*, ' '
          print*, 'starting thck, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
       endif

       if (verbose_calving .and. this_rank==rtest) then
          print*, ' '
          print*, 'calving_mask, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') calving_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
       endif
       
       do j = 1, ny
          do i = 1, nx
             if (thck(i,j) > 0.0d0 .and. calving_mask(i,j) == 1) then
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
                !TODO - Reset temperature and other tracers?
             endif
          enddo
       enddo

       if (verbose_calving .and. this_rank==rtest) then
          print*, ' '
          print*, 'new thck, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
       endif

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

          where (relx <= marine_limit + eus)
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_TOPG_THRESHOLD)   ! set thickness to zero if present bedrock is below a given level

          where (topg < marine_limit + eus)
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_HUYBRECHTS)    ! Huybrechts grounding line scheme for Greenland initialization

       !WHL - Previously, this code assumed that eus and relx have units of meters.
       !      Changed to be consistent with dimensionless thickness units.
!       if (eus > -80.d0) then
!          where (relx <= 2.d0*eus)
!             calving_thck = thck
!             thck = 0.0d0
!          end where
!       elseif (eus <= -80.d0) then
!          where (relx <= (2.d0*eus - 0.25d0*(eus + 80.d0)**2.d0))
!             calving_thck = thck
!             thck = 0.0d0
!          end where
!       end if
          if (eus*thk0 > -80.d0) then
             where (relx*thk0 <= 2.d0*eus*thk0)
                calving_law_mask = .true.
             elsewhere
                calving_law_mask = .false.
             end where
          elseif (eus*thk0 <= -80.d0) then
             where (relx*thk0 <= (2.d0*eus*thk0 - 0.25d0*(eus*thk0 + 80.d0)**2.d0))
                calving_law_mask = .true.
             elsewhere
                calving_law_mask = .false.
             end where
          end if

       case(CALVING_DAMAGE)   ! remove ice that is sufficiently damaged
                              !WHL - This is a rough initial implementation

          !WHL - debug - test damage field for MISMIP
!          print*, 'nx, ny =', nx, ny
!          print*, 'Prescribe damage:'
!          damage(:,:,:) = 0.d0
!          do j = nhalo+1, ny-nhalo
!             do i = nhalo+1, nx-nhalo
!                if (j == 3) then
!                   damage(:,i,j) = damage_threshold + 0.1d0
!                endif
!             enddo
!          enddo

          !TODO - Define ice_mask = 1 for dynamically active ice only (thck > thklim)?

          ! Diagnose the vertically integrated damage in each column,
          ! assuming the 3D damage field has been prognosed external to this subroutine.
          !WHL - For now, simply compute the thickness-weighted mean damage
          damage_column(:,:) = 0.0d0
          do j = 1, ny
             do i = 1, nx
                if (ice_mask(i,j) == 1) then
                   do k = 1, nz-1
                      damage_column(i,j) = damage_column(i,j) + damage(k,i,j) * (sigma(k+1) - sigma(k))
                   enddo
                endif
             enddo
          enddo
          
          ! set calving-law mask based on the vertically integrated damage
          where (damage_column > damage_threshold)  !WHL - could use '>=' instead of '>' if preferred
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

          !WHL - debug - print values of calving_law_mask
!          if (main_task) then
!             print*, 'i, j, damage, calving_law_mask:'
!             j = 3
!             do i = 1, nx
!                print*, i, j, damage_column(i,j), calving_law_mask(i,j)
!             enddo
!          endif

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

       elseif (calving_domain == CALVING_DOMAIN_OCEAN_CONNECT) then

          calving_domain_mask(:,:) = .false.

          ! initialize
          ! Assign the initial color to cells that meet the calving-law criteria and thus could calve,
          !  but only if they are connected to the ocean through other cells that meet the criteria.
          ! Assign the boundary color to cells that do not meet the calving-law criteria.

          allocate (color(nx,ny))
          do j = 1, ny
             do i = 1, nx
                if (calving_law_mask(i,j)) then
                   color(i,j) = initial_color
                else
                   color(i,j) = boundary_color
                endif
             enddo
          enddo

          ! Loop through cells, identifying cells that lie on the ocean margin.
          ! Fill each such cell with calving_law_mask = T, and recursively fill neighbor cells with calving_law_mask = T.
          ! We may have to do this several times to incorporate connections between neighboring processors.
          ! Alternatively, we could use a smaller value of maxcount_fill and allow calving to occur over multiple time steps,
          !  but this could make the calving evolution more dependent on ntasks than desired.

!          maxcount_fill = 1  ! setting maxcount_fill = 1 will ignore connections between neighboring processors.
          maxcount_fill = max(ewtasks,nstasks)

          do count = 1, maxcount_fill

             if (count == 1) then   ! identify margin cells that can seed the fill

                do j = 1, ny
                   do i = 1, nx
                      if ( floating_mask(i,j) == 1 .and. calving_law_mask(i,j) .and.  &
                           (ocean_mask(i-1,j)==1 .or. ocean_mask(i+1,j)==1 .or. &
                           ocean_mask(i,j-1)==1 .or. ocean_mask(i,j+1)==1) ) then
                         if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then
                            ! assign the fill color to this cell, and recursively fill neighbor cells
                            call glissade_fill(nx,    ny,       &
                                               i,     j,        &
                                               color, ice_mask)
                         endif
                      endif
                   enddo
                enddo

             else  ! count > 1; check for halo cells that were just filled on neighbor processors

                call parallel_halo(color)

                ! west halo layer
                i = nhalo
                do j = 1, ny
                   if (color(i,j) == fill_color) call glissade_fill(nx,    ny,    &
                                                                    i+1,   j,     &
                                                                    color, ice_mask)
                enddo

                ! east halo layer
                i = nx - nhalo + 1
                do j = 1, ny
                   if (color(i,j) == fill_color) call glissade_fill(nx,    ny,    &
                                                                    i-1,   j,     &
                                                                    color, ice_mask)
                enddo

                ! south halo layer
                j = nhalo
                do i = nhalo+1, nx-nhalo  ! already checked halo corners above
                   if (color(i,j) == fill_color) call glissade_fill(nx,    ny,    &
                                                                    i,     j+1,   &
                                                                    color, ice_mask)
                enddo

                ! north halo layer
                j = ny-nhalo+1
                do i = nhalo+1, nx-nhalo  ! already checked halo corners above
                   if (color(i,j) == fill_color) call glissade_fill(nx,    ny,    &
                                                                    i,     j-1,   &
                                                                    color, ice_mask)
                enddo

             endif  ! count = 1

             sum_fill_local = 0
             do j = nhalo+1, ny-nhalo
                do i = nhalo+1, nx-nhalo
                   if (color(i,j) == fill_color) sum_fill_local = sum_fill_local + 1
                enddo
             enddo

             !WHL - If running a large problem, may want to reduce the frequency of this global sum
             sum_fill_global = parallel_reduce_sum(sum_fill_local)

             if (verbose_calving) then
!!                print*, 'this_rank, sum_fill_local, sum_fill_global:', this_rank, sum_fill_local, sum_fill_global 
             endif

          enddo  ! count

          ! At this point, all cells with calving_law_mask = T should have the fill color if they
          !  are connected to the margin through other cells with calving_law_mask = T.
          !  These cells are now assigned to the calving domain.

          do j = 1, ny
             do i = 1, nx
                if (color(i,j) == fill_color) then
                   calving_domain_mask(i,j) = .true.
                else
                   calving_domain_mask(i,j) = .false.
                endif
             enddo
          enddo

          deallocate(color)

          ! Note: For this option, all cells in the calving domain have calving_law = T, so the logic below
          !  (calving_law_mask = T .and. calving_domain_mask = T) is redundant, but it does no harm.

       endif   ! calving_domain

       ! Calve ice where calving_law_mask = T and calving_domain_mask = T
       do j = 1, ny
          do i = 1, nx
             if (calving_law_mask(i,j) .and. calving_domain_mask(i,j)) then

                if (verbose_calving .and. this_rank==rtest .and. thck(i,j) > 0.0d0) then
!!                   print*, 'Calve ice: task, i, j, calving_thck =', this_rank, i, j, float_fraction_calve * thck(i,j)*thk0
                endif

                calving_thck(i,j) = calving_thck(i,j) + float_fraction_calve * thck(i,j)
                thck(i,j) = thck(i,j) - float_fraction_calve * thck(i,j)
            endif
          enddo
       enddo

       !WHL - debug
       if (verbose_calving .and. this_rank==rtest) then
          print*, ' '
          print*, 'After calving, new thck: itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
       endif

    endif   ! which_calving

    ! Optionally, impose a thickness limit on marine ice cliffs.
    ! These are defined as grounded marine-based cells adjacent to inactive calving_front cells or ice-free ocean.

    if (limit_marine_cliffs) then

       ! Update masks, including the marine_cliff mask

       call glissade_get_masks(nx,            ny,                 &
                               thck,          topg,               &
                               eus,           thklim,             &
                               ice_mask,                          &
                               floating_mask = floating_mask,     &
                               ocean_mask = ocean_mask,           &
                               land_mask = land_mask,             &
                               active_ice_mask = active_ice_mask, &
                               which_ho_calving_front = which_ho_calving_front, &
                               calving_front_mask = calving_front_mask, &
                               thck_calving_front = thck_calving_front, &
                               marine_cliff_mask = marine_cliff_mask)

       if (verbose_calving .and. this_rank==rtest) then
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
       endif

       if (verbose_calving .and. this_rank==rtest) then
          print*, ' '
          print*, 'thckmax_cliff, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                factor = taumax_cliff / (rhoi*grav)   ! units are Pa for taumax, m for factor
                thckmax_cliff = factor + sqrt(factor**2 + (rhoo/rhoi)*(topg(i,j)*thk0)**2)  ! m 
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
                thckmax_cliff = factor + sqrt(factor**2 + (rhoo/rhoi)*(topg(i,j)*thk0)**2)  ! m 
                thckmax_cliff = thckmax_cliff / thk0   ! convert to model units

                !WHL - debug
                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   print*, ' '
                   print*, 'Cliff thinning: r, i, j =', rtest, itest, jtest
                   print*, 'thck, thckmax_cliff (m) =', thck(i,j)*thk0, thckmax_cliff*thk0
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
!!                      print*, ' '
!!                      print*, 'r, i, j, thck, thckmax_cliff:', &
!!                           this_rank, i, j, thck(i,j)*thk0, thckmax_cliff*thk0 
!!                      print*, 'thinning rate (model units) =', thinning_rate
                      print*, 'thinning rate (m/yr) =', thinning_rate * thk0*scyr/tim0
                      print*, 'dthck (m) =', dthck * thk0
                   endif

                   thck(i,j) = thck(i,j) - dthck
                   calving_thck(i,j) = calving_thck(i,j) + dthck

                endif  ! thck > thckmax_cliff

             endif  ! marine_cliff cell
          enddo   ! i
       enddo   ! j

    endif   ! limit_marine_cliffs

    ! Remove any icebergs.
    ! Typically these will be removed by the calving scheme above, but if not, 
    !  then they need to be removed before calling the velocity solver,
    !  which will have problems in regions without any grounded ice.

    if (remove_icebergs) then

       if (verbose_calving .and. main_task) then
          print*, 'Remove icebergs'
       endif

       ! Set the thickness thresholds for cells to be counted as active.
       ! Floating regions will be identified as icebergs unless they are connected to grounded ice
       !  along a path consisting only of active cells.
       ! Note: A limit of 0.0 does not work because it erroneously counts very thin floating cells as active.
       !       Then the algorithm can fail to identify floating regions that should be removed.

       call glissade_remove_icebergs(&
            itest,   jtest,   rtest, &
            thck,                    &
            topg,          eus,      &
            thklim,                  &
            which_ho_calving_front,  &
            calving_thck,            &
            cull_calving_front,      &
            ncull_calving_front)

    endif

    if (verbose_calving .and. this_rank == rtest) then

       print*, ' '
       print*, 'Final calving_thck (m), itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') calving_thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'Final thck (m), itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo

    endif  ! verbose_calving

    ! cleanup
    deallocate (calving_law_mask)
    deallocate (calving_domain_mask)
    deallocate (ice_mask)
    deallocate (floating_mask)
    deallocate (ocean_mask)
    deallocate (land_mask)
    deallocate (active_ice_mask)
    deallocate (calving_front_mask)
    deallocate (thck_calving_front)
    deallocate (thck_init)
    deallocate (marine_cliff_mask)

    if (allocated(thck_init)) deallocate(thck_init)

  end subroutine glissade_calve_ice

!---------------------------------------------------------------------------

  subroutine glissade_remove_icebergs(&
       itest,   jtest,   rtest,     &
       thck,                        &
       topg,          eus,          &
       thklim,                      &
       which_ho_calving_front,      &
       calving_thck,                &
       cull_calving_front_in,       &
       ncull_calving_front_in)

    ! Remove any icebergs. 
        
    ! The algorithm is as follows:
    ! (1) Mark all cells with ice (either active or inactive) with the initial color.
    !     Mark other cells with the boundary color.
    ! (2) Seed the fill by giving all active grounded cells the fill color.
    ! (3) Recursively fill all cells that are connected to filled cells by a path
    !     that passes through active cells only.
    ! (4) Repeat the recursion as necessary to spread the fill to adjacent processors.
    ! (5) Once the fill is done, any cells that still have the initial color and
    !     are not on land are considered to be icebergs. They are removed.
    !
    ! Notes:
    ! (1) The recursive fill applies to edge neighbors, not corner neighbors.
    !     The path back to grounded ice must go through edges, not corners.
    ! (2) Inactive cells can be filled (if adjacent to active cells), but
    !     do not further spread the fill.
    ! (3) Grounded cells that still have the initial color are not removed.
    !     They are considered harmless.

    use glissade_masks
    use glimmer_paramets, only: thk0

    integer, intent(in) :: itest, jtest, rtest          !> coordinates of diagnostic point

    real(dp), dimension(:,:), intent(inout) :: thck     !> ice thickness
    real(dp), dimension(:,:), intent(in)    :: topg     !> present bedrock topography
    real(dp), intent(in)    :: eus                      !> eustatic sea level
    real(dp), intent(in)    :: thklim                   !> minimum thickness for dynamically active grounded ice
    integer, intent(in)     :: which_ho_calving_front   !> = 1 for subgrid calving-front scheme, else = 0
    real(dp), dimension(:,:), intent(inout) :: calving_thck   !> thickness lost due to calving in each grid cell;
                                                              !> on output, includes ice in icebergs

    logical, intent(in), optional :: &
         cull_calving_front_in         !> if true, remove peninsulas by first removing a layer of calving_front cells

    integer, intent(in), optional :: &
         ncull_calving_front_in        !> number of times to cull calving_front cells at initialization

    ! local variables

    logical :: cull_calving_front    ! local version of cull_calving_front_in
    integer :: ncull_calving_front   ! local version of ncull_calving_front_in

    integer :: nx, ny                ! horizontal grid dimensions

    integer :: i, j, n
    integer :: count, maxcount_fill  ! loop counters

    integer,  dimension(:,:), allocatable   ::  &
         ice_mask,           & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,      & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,         & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,          & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask, & ! = 1 where ice is floating and borders the ocean, else = 0
         active_ice_mask,    & ! = 1 for dynamically active cells
         color                 ! integer 'color' for identifying icebergs

    real(dp),  dimension(:,:), allocatable   ::  &
         thck_calving_front    ! effective ice thickness at the calving front

    !WHL - debug
    real(dp) :: sum_fill_local, sum_fill_global

    if (present(cull_calving_front_in)) then
       cull_calving_front = cull_calving_front_in
    else
       cull_calving_front = .false.
    endif

    if (present(ncull_calving_front_in)) then
       ncull_calving_front = ncull_calving_front_in
    else
       ncull_calving_front = 1
    endif

    nx = size(thck,1)
    ny = size(thck,2)

    allocate (ice_mask(nx,ny))
    allocate (floating_mask(nx,ny))
    allocate (ocean_mask(nx,ny))
    allocate (land_mask(nx,ny))
    allocate (calving_front_mask(nx,ny))
    allocate (thck_calving_front(nx,ny))
    allocate (active_ice_mask(nx,ny))
    allocate (color(nx,ny))

    ! calculate masks
    ! Note: Passing in thklim = 0.0 does not work because it erroneously counts thin floating cells as active.
    !       Then the algorithm can fail to identify floating regions that should be removed
    !       (since they are separated from any active cells).

    call glissade_get_masks(nx,            ny,                  &
                            thck,          topg,                &
                            eus,           thklim,              &
                            ice_mask,                           &
                            floating_mask = floating_mask,      &
                            ocean_mask = ocean_mask,            &
                            land_mask = land_mask,              &
                            active_ice_mask = active_ice_mask,  &
                            which_ho_calving_front = which_ho_calving_front, &
                            calving_front_mask = calving_front_mask,         &
                            thck_calving_front = thck_calving_front)

    !WHL - debug
    if (verbose_calving .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_remove_icebergs'
       print*, ' '
       print*, 'thck, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)*thk0
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
       print*, 'thck_calving_front, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck_calving_front(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'active_ice_mask, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') active_ice_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif
    
    !WHL - debug
    ! Optionally, do a preliminary step where all cells currently at the calving front are removed.
    ! Then recompute the masks.
    ! The results is that long, skinning floating peninsulas that can be dynamically unstable are more likely
    !  to be removed. Without this step, peninsulas that are two cells thick (with calving-front cells on each side)
    !  will typically be removed as icebergs (because there is no path back to grounded ice through active cells).
    !  With this step, peninsulas up to four cells thick will be removed (two outer layers during the preliminary step,
    !  followed by two inner layers on the remove_iceberg step).
    ! If necessary, this step could be repeated to remove peninsulas with a thickness of 6 layers, 8 layers, etc.

    if (cull_calving_front) then

       do n = 1, ncull_calving_front

          ! remove the calving_front cells just identified
          if (main_task) then
             call write_log ('cull_calving_front: Removing ice from calving_front cells')
             print*, 'cull_calving_front: Removing ice from calving_front cells'
          endif

          do j = 1, ny
             do i = 1, nx
                if (calving_front_mask(i,j) == 1) then
                   calving_thck(i,j) = calving_thck(i,j) + thck(i,j) 
                   thck(i,j) = 0.0d0
                endif
             enddo
          enddo

          ! update the masks
          ! Note: Some floating cells that were previously active interior cells are now calving_front cells.
          !       These will not be removed if they are adjacent to active cells with a path to grounded ice,
          !        but will be removed if they form peninsulas one or two cells thick.

          call glissade_get_masks(nx,            ny,                  &
                                  thck,          topg,                &
                                  eus,           thklim,              &
                                  ice_mask,                           &
                                  floating_mask = floating_mask,      &
                                  ocean_mask = ocean_mask,            &
                                  land_mask = land_mask,              &
                                  active_ice_mask = active_ice_mask,  &
                                  which_ho_calving_front = which_ho_calving_front, &
                                  calving_front_mask = calving_front_mask,         &
                                  thck_calving_front = thck_calving_front)

          !WHL - debug
          if (verbose_calving .and. this_rank == rtest) then
             print*, ' '
             print*, 'cull_calving_front: After removing CF cells, n =', n
             print*, ' '
             print*, 'thck, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') thck(i,j)*thk0
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
             print*, 'thck_calving_front, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.3)',advance='no') thck_calving_front(i,j)*thk0
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'active_ice_mask, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(i10)',advance='no') active_ice_mask(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif

       enddo  ! ncull_calving_front

    endif  ! cull_calving_front

    ! initialize
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

    maxcount_fill = max(ewtasks,nstasks)

    if (verbose_calving .and. main_task) then
       print*, 'maxcount_fill =', maxcount_fill
    endif

    do count = 1, maxcount_fill

       if (count == 1) then   ! identify grounded cells that can seed the fill

          do j = 1, ny
             do i = 1, nx
                if (active_ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) then  ! grounded ice
                   if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then
                      ! assign the fill color to this cell, and recursively fill neighbor cells
                      call glissade_fill(nx,    ny,    &
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
                call glissade_fill(nx,    ny,    &
                                   i+1,   j,     &
                                   color, active_ice_mask)
             endif
          enddo

          ! east halo layers
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i-1,   j,     &
                                   color, active_ice_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i,     j+1,   &
                                   color, active_ice_mask)
             endif
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i,     j-1,   &
                                   color, active_ice_mask)
             endif
          enddo

       endif  ! count = 1

       sum_fill_local = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) sum_fill_local = sum_fill_local + 1
          enddo
       enddo

       !WHL - If running a large problem, may want to reduce the frequency of this global sum
       sum_fill_global = parallel_reduce_sum(sum_fill_local)

       if (verbose_calving .and. main_task) then
!!          print*, 'this_rank, sum_fill_local, sum_fill_global:', this_rank, sum_fill_local, sum_fill_global
       endif

    enddo  ! count

    ! Icebergs are cells that still have the initial color and are not on land.
    ! Remove ice in these cells, adding it to the calving field.
    ! Note: Inactive land-based cells are not considered to be icebergs.
    ! Note: Another exception is that we do not remove cells that are
    !       (1) adjacent to at least one floating cell (sharing an edge), and
    !       (2) connected diagonally to active cells with the fill color.
    !       Such cells are considered part of the inactive calving front and are
    !        allowed to continue filling instead of calving.


    do j = 2, ny-1
       do i = 2, nx-1
          if (color(i,j) == initial_color .and. land_mask(i,j) == 0) then
             if (  ( color(i-1,j+1)==fill_color .and. active_ice_mask(i-1,j+1)==1 .and. &
                       (floating_mask(i-1,j)==1 .or. floating_mask(i,j+1)==1) ) &
              .or. ( color(i+1,j+1)==fill_color .and. active_ice_mask(i+1,j+1)==1 .and. &
                       (floating_mask(i+1,j)==1 .or. floating_mask(i,j+1)==1) ) &
              .or. ( color(i-1,j-1)==fill_color .and. active_ice_mask(i-1,j-1)==1 .and. &
                       (floating_mask(i-1,j)==1 .or. floating_mask(i,j-1)==1) ) &
              .or. ( color(i+1,j-1)==fill_color .and. active_ice_mask(i+1,j-1)==1 .and. &
                       (floating_mask(i+1,j)==1 .or. floating_mask(i,j-1)==1) ) ) then
                ! do nothing; this cell is part of the inactive calving front
             else  ! not part of the inactive calving front; calve as an iceberg
                !WHL - debug
!!                if (verbose_calving .and. thck(i,j) > 0.0 .and. this_rank == rtest) then
!!                   print*, 'Remove iceberg: task, i, j, thck =', this_rank, i, j, thck(i,j)*thk0
!!                endif
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
                !TODO - Also handle tracers?  E.g., set damage(:,i,j) = 0.d0?
             endif

          endif
       enddo
    enddo

    ! cleanup
    deallocate (ice_mask)
    deallocate (floating_mask)
    deallocate (ocean_mask)
    deallocate (calving_front_mask)
    deallocate (thck_calving_front)
    deallocate (active_ice_mask)
    deallocate (color)

    if (verbose_calving .and. this_rank == rtest) then
       print*, ' '
       print*, 'Done in glissade_remove_icebergs'
       print*, ' '
       print*, 'thck, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f10.3)',advance='no') thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_remove_icebergs

!****************************************************************************

  subroutine glissade_find_lakes(nx,           ny,             &
                                 itest, jtest, rtest,          &
                                 ice_mask,     floating_mask,  &
                                 ocean_mask,   lake_mask)

  ! Identify cells with basal lakes: i.e., cells that are floating but have
  ! no connection through other floating cells to the ocean.

  !TODO - Move this subroutine elsewhere (e.g., glissade_basal_traction)

    integer, intent(in) :: nx, ny                  !> horizontal grid dimensions

    integer, intent(in) :: itest, jtest, rtest     !> coordinates of diagnostic point

    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,               & !> = 1 where ice is present (thck > thklim), else = 0
         floating_mask,          & !> = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask                !> = 1 where topg is below sea level and ice is absent, else = 0

    integer, dimension(nx,ny), intent(out) ::  &
         lake_mask                 !> = 1 for floating cells disconnected from the ocean, else = 0

    ! local variables

    integer, dimension(nx,ny) ::  &
         color                     ! integer 'color' for identifying icebergs

    integer :: i, j
    integer :: count, maxcount_fill  ! loop counters

    logical, parameter :: verbose_lakes = .false.

    !WHL - debug
    real(dp) :: sum_fill_local, sum_fill_global
    integer :: ig, jg

    if (verbose_lakes .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_find_lakes, itest, jtest, rank =', itest, jtest, rtest
       print*, ' '
       print*, 'ice_mask'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') ice_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'floating_mask'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') floating_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! initialize
    ! Floating cells receive the initial color;
    !  grounded cells and ice-free cells receive the boundary color.

    do j = 1, ny
       do i = 1, nx
          if (floating_mask(i,j) == 1) then
             color(i,j) = initial_color
          else    ! grounded or ice-free
             color(i,j) = boundary_color
          endif
       enddo
    enddo

    ! Loop through cells, identifying floating cells that border the ocean.
    ! Fill each such floating cell, and then recursively fill floating neighbor cells.
    ! We may have to do this several times to incorporate connections between neighboring processors.

    maxcount_fill = max(ewtasks,nstasks)

    if (verbose_lakes .and. main_task) print*, 'maxcount_fill =', maxcount_fill

    do count = 1, maxcount_fill

       if (count == 1) then   ! identify floating cells adjacent to ocean cells, which can seed the fill

          do j = 2, ny-1
             do i = 2, nx-1
                if (floating_mask(i,j) == 1) then
                   if (ocean_mask(i-1,j) == 1 .or. ocean_mask(i+1,j) == 1 .or.   &
                       ocean_mask(i,j-1) == 1 .or. ocean_mask(i,j+1) == 1) then

                      if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then

                         ! assign the fill color to this cell, and recursively fill floating neighbor cells
                         call glissade_fill(nx,    ny,    &
                                            i,     j,     &
                                            color, floating_mask)
                      endif
                   endif  ! adjacent to ocean
                endif  ! floating
             enddo  ! i
          enddo  ! j

       else  ! count > 1

          ! Check for halo cells that were just filled on neighbor processors
          ! Note: In order for a halo cell to seed the fill on this processor, it must not only have the fill color,
          !       but also must be an active cell.
          !TODO - Only need to check whether color = fill_color?

          call parallel_halo(color)

          ! west halo layer
          i = nhalo
          do j = 1, ny
             if (color(i,j) == fill_color .and. floating_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i+1,   j,     &
                                   color, floating_mask)
             endif
          enddo

          ! east halo layers
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color .and. floating_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i-1,   j,     &
                                   color, floating_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. floating_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i,     j+1,   &
                                   color, floating_mask)
             endif
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. floating_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i,     j-1,   &
                                   color, floating_mask)
             endif
          enddo

       endif  ! count = 1

       sum_fill_local = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) sum_fill_local = sum_fill_local + 1
          enddo
       enddo

       !WHL - If running a large problem, may want to reduce the frequency of this global sum
       sum_fill_global = parallel_reduce_sum(sum_fill_local)

       if (verbose_lakes .and. main_task) then
          print*, 'this_rank, sum_fill_local, sum_fill_global:', this_rank, sum_fill_local, sum_fill_global
       endif

    enddo  ! count

    ! Identify lake cells: floating cells that still have the initial color.

    lake_mask(:,:) = 0

    do j = 1, ny
       do i = 1, nx
          if (color(i,j) == initial_color .and. floating_mask(i,j) == 1) then
             lake_mask(i,j) = 1

             if (verbose_lakes .and. this_rank == rtest) then
                call parallel_globalindex(i, j, ig, jg)
                print*, 'Lake cell: task, i, j, ig, jg =', this_rank, i, j, ig, jg
             endif

          endif
       enddo
    enddo

    call parallel_halo(lake_mask)

    if (verbose_lakes .and. this_rank == rtest) then
       print*, ' '
       print*, 'color, rank =', this_rank
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') color(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'floating_mask, rank =', this_rank
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') floating_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'lake_mask, rank =', this_rank
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') lake_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_find_lakes

!****************************************************************************

  recursive subroutine glissade_fill(nx,  ny,         &
                                     i,   j,          &
                                     color,           &
                                     active_ice_mask)

    ! Given a domain with an initial color, a boundary color and a fill color,
    ! assign the fill color to all cells that either (1) are prescribed to have
    ! the fill color or (2) are connected to cells with the fill color.

    integer, intent(in) :: nx, ny                       !> domain size
    integer, intent(in) :: i, j                         !> horizontal indices of current cell

    integer, dimension(nx,ny), intent(inout) :: &
         color                                          ! color (initial, fill or boundary)

    integer, dimension(nx,ny), intent(in) :: &
         active_ice_mask                                ! true for dynamically active ice

    if (color(i,j) /= fill_color .and. color(i,j) /= boundary_color) then

       ! assign the fill color to this cell
       color(i,j) = fill_color

       ! If a cell contains inactive ice, then fill this cell but do not call
       !  glissade_fill recursively
       if (active_ice_mask(i,j) == 0) then  ! this cell is inactive
          return    ! skip the recursion
       endif

       ! recursively call this subroutine for each neighbor to see if it should be filled
       !TODO - May want to rewrite this to avoid recursion, which can crash the code when
       !       the recursion stack is very large on fine grids.
       if (i > 1)  call glissade_fill(nx,    ny,  &
                                      i-1,   j,   &
                                      color, active_ice_mask)
       if (i < nx) call glissade_fill(nx,    ny,  &
                                      i+1,   j,   &
                                      color, active_ice_mask)
       if (j > 1)  call glissade_fill(nx,    ny, &
                                      i,     j-1, &
                                      color, active_ice_mask)
       if (j < ny) call glissade_fill(nx,    ny,  &
                                      i,     j+1, &
                                      color, active_ice_mask)

    endif   ! not fill color or boundary color

  end subroutine glissade_fill

!---------------------------------------------------------------------------

end module glissade_calving

!---------------------------------------------------------------------------
