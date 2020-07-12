!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_masks.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains routines for computing various masks used by the Glissade 
! velocity solver.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_masks

    use glimmer_global, only: dp
    use glimmer_log
    use glimmer_physcon, only: rhoi, rhoo
    use glide_types
    use parallel

    implicit none

    private
    public :: glissade_get_masks, glissade_calving_front_mask, &
              glissade_marine_cliff_mask, glissade_ice_sheet_mask,  &
              glissade_marine_connection_mask, glissade_lake_mask,  &
              glissade_extend_mask, glissade_fill, glissade_fill_with_buffer
    public :: initial_color, fill_color, boundary_color

    ! colors for fill subroutines
    integer, parameter :: initial_color = 0   ! initial color, represented by integer
    integer, parameter :: fill_color = 1      ! fill color, represented by integer
    integer, parameter :: boundary_color = -1 ! boundary color, represented by integer

  contains

!****************************************************************************

  subroutine glissade_get_masks(nx,          ny,          &
                                thck,        topg,        &
                                eus,         thklim,      &
                                ice_mask,                 &
                                floating_mask,            &
                                ocean_mask,               &
                                land_mask,                &
                                grounding_line_mask,      &
                                active_ice_mask)

    !----------------------------------------------------------------
    ! Compute various masks for the Glissade dycore.
    !
    ! There are different ways to handle masks.
    ! The approach in Glide was to define an integer cell_mask array,
    !  in which each bit encodes information about whether a cell is ice-covered
    !  or ice-free, floating or grounded, etc.
    ! The approach here is to compute a separate 2D integer array for
    !  each kind of mask. This uses more memory but also is more transparent.
    !
    ! The basic masks used for Glissade dynamic calculations are as follows:
    !
    ! (1) ice_mask = 1 where ice is present (thck > thklim), else = 0
    ! (2) floating_mask = 1 if ice is present (thck > thklim) and floating, else = 0
    ! (3) ocean_mask = 1 if the topography is below sea level (topg < eus) and thk <= thklim, else = 0
    ! (4) land_mask = 1 if the topography is at or above sea level (topg >= eus), else = 0
    ! (5) grounding_line_mask = 1 if a cell is adjacent to the grounding line, else = 0
    ! (6) active_ice_mask = 1 for dynamically active cells, else = 0
    !     With the subgrid calving front scheme, cells that lie on the calving front and have
    !     thck < thck_calving_front are inactive. Otherwise, all cells with ice_mask = 1 are active.
    !
    ! where thck = ice thickness
    !       thklim = threshold thickness for ice to be dynamically active
    !       topg = bed topography
    !       eus = eustatic sea level (= 0 by default)
    !
    ! Notes:
    ! (1) thck, thklim, topg and eus can either have units of length (e.g., meters)
    !     or be dimensionless, as long as they are defined consistently.
    ! (2) ice_mask is always computed; the other masks are optional.
    ! (3) Thermal calculations may have a different threshold, thklim_temp
    !     (where generally thklim_temp < thklim).
    !     This mask can be computed by replacing thklim with thklim_temp in the subroutine call.
    ! (4) For some calculations it may be useful to call this subroutine with thklim = 0
    !     so as to identify all cells with nonzero ice thickness, not just dynamically active cells.
    !----------------------------------------------------------------
    
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
         nx,  ny                ! number of grid cells in each direction

    ! Default dimensions are meters, but this subroutine will work for
    ! any units as long as thck, topg, eus and thklim have the same units.

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                & ! ice thickness (m)
         topg                   ! elevation of topography (m)

    real(dp), intent(in) :: &
         eus,                 & ! eustatic sea level (m), = 0. by default
         thklim                 ! minimum ice thickness for active cells (m)

    integer, dimension(nx,ny), intent(out) ::  &
         ice_mask               ! = 1 if thck > thklim, else = 0  

    integer, dimension(nx,ny), intent(out), optional ::  &
         floating_mask,       & ! = 1 if thck > thklim and ice is floating, else = 0
         ocean_mask,          & ! = 1 if topg is below sea level and thk <= thklim, else = 0
         land_mask,           & ! = 1 if topg is at or above sea level, else = 0
         grounding_line_mask, & ! = 1 if a cell is adjacent to the grounding line, else = 0
         active_ice_mask        ! = 1 if dynamically active, else = 0

    !----------------------------------------------------------------
    ! Local arguments
    !----------------------------------------------------------------

    integer :: i, j, ii, jj

    integer, dimension(nx,ny) ::  &
         grounded_mask          ! = 1 if ice is present and grounded, else = 0

    !----------------------------------------------------------------
    ! Compute masks in cells
    !----------------------------------------------------------------

    ice_mask(:,:) = 0

    do j = 1, ny
       do i = 1, nx

          if (thck(i,j) > thklim) then
             ice_mask(i,j) = 1
          else
             ice_mask(i,j) = 0
          endif

          if (present(ocean_mask)) then
             if (topg(i,j) < eus .and. ice_mask(i,j) == 0) then
                ocean_mask(i,j) = 1
             else
                ocean_mask(i,j) = 0
             endif
          endif

          if (present(floating_mask)) then
             if (topg(i,j) - eus < (-rhoi/rhoo)*thck(i,j) .and. ice_mask(i,j) == 1) then
                floating_mask(i,j) = 1
             else
                floating_mask(i,j) = 0
             endif
          endif

          if (present(land_mask)) then
             if (topg(i,j) >= eus) then
                land_mask(i,j) = 1
             else
                land_mask(i,j) = 0
             endif
          endif

          ! Note: active_ice_mask will be overwritten if the subgrid calving front scheme is used
          if (present(active_ice_mask)) then
             active_ice_mask(i,j) = ice_mask(i,j)
          endif

       enddo  ! i
    enddo  ! j

    ! halo updates
    ! Note: These are not strictly needed because the above loops include halo cells.
    !       However, they are included in case the user calls this subroutine without
    !        first updating thck in halo cells.

    call parallel_halo(ice_mask)
    if (present(floating_mask)) call parallel_halo(floating_mask)
    if (present(active_ice_mask)) call parallel_halo(active_ice_mask)

    ! Identify grounded cells; this mask is used in some calculations below
    if (present(floating_mask)) then
       where (ice_mask == 1 .and. floating_mask == 0)
          grounded_mask = 1
       elsewhere
          grounded_mask = 0
       endwhere
    endif

    ! Optionally, compute grounding line mask using grounded_mask, floating_mask and ocean_mask

    if (present(grounding_line_mask)) then

       if (.not.present(floating_mask) .or. .not.present(ocean_mask)) then
          call write_log('Need floating_mask and ocean_mask to compute grounding_line_mask', GM_FATAL)
       endif

       grounding_line_mask(:,:) = 0

       do j = 2, ny-1
          do i = 2, nx-1

             if (grounded_mask(i,j) == 1) then
                ! check whether one or more neighbors is a floating or ocean cell
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (floating_mask(ii,jj) == 1 .or. ocean_mask(ii,jj) == 1) then
                         grounding_line_mask(i,j) = 1
                      endif
                   enddo
                enddo
             elseif (floating_mask(i,j) == 1) then
                ! check whether one or more neighbors is a grounded cell
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (grounded_mask(ii,jj) == 1) then
                         grounding_line_mask(i,j) = 1
                      endif
                   enddo
                enddo
             endif   ! grounded_mask or floating_mask

          enddo   ! i
       enddo   ! j

       call parallel_halo(grounding_line_mask)

    endif   ! present(grounding_line_mask)

    ! Note: Halo calls are not included for the ocean and land masks.
    !       Halo values will still be correct, provided that topg is correct in halo cells.
    !       The reason not to include these calls is that for outflow global BCs,
    !        we may not want to set ocean_mask and land_mask = 0 in halo cells (as would be
    !        done automatically for outflow BCs). Instead, we want to compute ocean_mask
    !        and land_mask in halo cells based on topg (which for outflow BCs is extrapolated
    !        to halo cells from adjacent physical cells). 
    !       In particular, setting ocean_mask = 1 in the global halo ensures that calving_front
    !        cells are treated correctly just inside the global halo.

!!    if (present(ocean_mask)) call parallel_halo(ocean_mask)
!!    if (present(land_mask)) call parallel_halo(land_mask)

  end subroutine glissade_get_masks

!****************************************************************************

  subroutine glissade_calving_front_mask(&
       nx,                     ny,                   &
       which_ho_calving_front,                       &
       thck,                   topg,                 &
       eus,                                          &
       ice_mask,               floating_mask,        &
       ocean_mask,             land_mask,            &
       calving_front_mask,     thck_calving_front,   &
       active_ice_mask,        marine_interior_mask, &
       effective_areafrac,     marine_cliff_mask)

    ! Compute a calving_front mask, effective calving_front thickness, and related fields.
    ! Note: With the subgrid calving front scheme, cells that lie on the calving front and have
    !       thck < thck_calving_front are inactive. Otherwise, all cells with ice_mask = 1 are active.

    integer, intent(in) ::   &
         nx,  ny,              &  ! number of grid cells in each direction
         which_ho_calving_front   ! subgrid calving front option

    ! Default dimensions are meters, but this subroutine will work for any units
    !  as long as thck, topg, and eus have the same units.

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                  & ! ice thickness (m)
         topg                     ! elevation of topography (m)

    real(dp), intent(in) :: &
         eus                      ! eustatic sea level (m), = 0. by default

    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,              & ! = 1 if thck > thklim, else = 0
         floating_mask,         & ! = 1 if thck > thklim and ice is floating, else = 0
         ocean_mask,            & ! = 1 if topg is below sea level and thk <= thklim, else = 0
         land_mask                ! = 1 if topg is at or above sea level, else = 0

    integer, dimension(nx,ny), intent(out) ::  &
         calving_front_mask       ! = 1 if ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(out) :: &
         thck_calving_front       ! effective ice thickness at the calving front

    integer, dimension(nx,ny), intent(out), optional ::  &
         active_ice_mask,       & ! = 1 if dynamically active, else = 0
         marine_interior_mask,  & ! = 1 if ice is marine-based and borders no ocean cells, else = 0
         marine_cliff_mask        ! = 1 if ice is grounded and marine-based and borders at least one ocean
                                  !     or inactive calving_front cell, else = 0

    real(dp), dimension(nx,ny), intent(out), optional :: &
         effective_areafrac       ! effective ice-covered fraction, in range [0,1]
                                  ! 0 < f < 1 for partial calving-front cells

    !----------------------------------------------------------------
    ! Local arguments
    !----------------------------------------------------------------

    integer :: i, j, ii, jj

    integer, dimension(nx,ny) :: &
         interior_marine_mask   ! same as marine_interior mask; used internally

    real(dp), dimension(nx,ny) :: &
         thck_flotation         ! flotation thickness

    integer :: sum_cell         ! temporary sums
    real(dp) :: sum_thck

    ! Compute a calving front mask and effective calving front thickness.
    ! Optionally, compute some related fields.

    if (which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then

       calving_front_mask(:,:) = 0
       interior_marine_mask(:,:) = 0

       ! Identify calving front cells (floating cells that border ice-free ocean)
       ! and marine-based interior cells (marine-based cells not at the calving front).
       do j = 2, ny-1
          do i = 2, nx-1
             if (floating_mask(i,j) == 1) then
                if (ocean_mask(i-1,j) == 1 .or. ocean_mask(i+1,j) == 1 .or. &
                    ocean_mask(i,j-1) == 1 .or. ocean_mask(i,j+1) == 1) then
                   calving_front_mask(i,j) = 1
                else
                   interior_marine_mask(i,j) = 1
                endif
             elseif (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0 .and. topg(i,j) < eus) then  ! grounded marine-based ice
                interior_marine_mask(i,j) = 1
             endif
          enddo
       enddo

       call parallel_halo(calving_front_mask)
       call parallel_halo(interior_marine_mask)

       ! Compute thck_calving_front, an effective thickness for calving-front cells.
       ! It is set to the mean thickness in adjacent marine interior cells.
       ! Note: For CF cells without any marine interior neighbors, we return thck_calving_front = 0.
       ! TODO: Make sure this doesn't lead to numerical problems.

       thck_calving_front(:,:) = 0.0d0

       do j = 2, ny-1
          do i = 2, nx-1
             if (calving_front_mask(i,j) == 1) then

                sum_cell = 0
                sum_thck = 0.0d0

                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (interior_marine_mask(ii,jj) == 1) then
                         sum_cell = sum_cell + 1
                         sum_thck = sum_thck + thck(ii,jj)
                      endif
                   enddo
                enddo
                if (sum_cell > 0) then
                   thck_calving_front(i,j) = sum_thck/sum_cell
                endif

             endif
          enddo
       enddo

       call parallel_halo(thck_calving_front)

       ! Limit thck_calving_front so as not to exceed the flotation thickness
       where (thck_calving_front > 0.0d0)
          thck_flotation = -(rhoo/rhoi) * (topg - eus)
          thck_calving_front = min(thck_calving_front, thck_flotation)
       endwhere

       ! Optionally, copy interior_marine_mask to marine_interior_mask for output.
       ! The reason to have two copies of the same mask is to allow thck_calving_front to be computed,
       !   whether or not marine_interior_mask is present.

       if (present(marine_interior_mask)) then
          marine_interior_mask = interior_marine_mask
       endif

       ! Optionally, use the ratio thck/thck_calving_front to compute effective_areafrac.
       ! TODO - Think about whether we should have effective_areafrac = 1 for ice-free land.

       if (present(effective_areafrac)) then

          do j = 1, ny
             do i = 1, nx
                if (calving_front_mask(i,j) == 1 .and. thck_calving_front(i,j) > 0.0d0) then
                   effective_areafrac(i,j) = thck(i,j) / thck_calving_front(i,j)
                   effective_areafrac(i,j) = min(effective_areafrac(i,j), 1.0d0)
                elseif (ocean_mask(i,j) == 1) then
                   effective_areafrac(i,j) = 0.0d0
                else  ! non-CF ice-covered cells and/or land cells
                   effective_areafrac(i,j) = 1.0d0
                endif
             enddo
          enddo

       endif   ! present(effective_areafrac)

       ! Optionally, update the active_ice_mask so that CF cells with thck < thck_calving_front are inactive,
       ! but those with thck >= thck_calving_front are active.

       if (present(active_ice_mask)) then

          ! initialize
          active_ice_mask(:,:) = 0

          ! Mark ice-filled cells as active.
          ! Calving-front cells, however, are inactive, unless they have thck >= thck_calving front.

          do j = 2, ny-1
             do i = 2, nx-1
                if (ice_mask(i,j) == 1) then
                   if (calving_front_mask(i,j) == 0) then
                      active_ice_mask(i,j) = 1
                   elseif (calving_front_mask(i,j) == 1) then
                      !WHL - If two adjacent cells are being restored to the same thickness, there is a
                      !       chance of flickering here, with the CF cell alternately being restored to
                      !       slightly greater or slightly less than the thickness of its interior neighbor.
                      !      For this reason, let the cell be active if thck is very close to thck_calving front,
                      !       but slightly less.
                      if (thck_calving_front(i,j) > 0.0d0 .and. &
                          thck(i,j) >= 0.999d0*thck_calving_front(i,j)) then
                         active_ice_mask(i,j) = 1
                      endif
                   endif   ! calving_front_mask
                endif  ! ice_mask
             enddo
          enddo

          call parallel_halo(active_ice_mask)

       endif   ! present(active_ice_mask)

    else   ! no subgrid calving front scheme

       calving_front_mask(:,:) = 0
       interior_marine_mask(:,:) = 0
       thck_calving_front(:,:) = 0.0d0

       ! Identify calving front cells (floating cells that border ice-free ocean)
       ! and marine-based interior cells (marine-based cells not at the calving front).

       do j = 2, ny-1
          do i = 2, nx-1
             if (floating_mask(i,j) == 1) then
                if (ocean_mask(i-1,j) == 1 .or. ocean_mask(i+1,j) == 1 .or. &
                    ocean_mask(i,j-1) == 1 .or. ocean_mask(i,j+1) == 1) then
                   calving_front_mask(i,j) = 1
                   thck_calving_front(i,j) = thck(i,j)
                else
                   interior_marine_mask(i,j) = 1
                endif
             elseif (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0 .and. topg(i,j) < eus) then  ! grounded marine-based ice
                interior_marine_mask(i,j) = 1
             endif
          enddo
       enddo

       call parallel_halo(calving_front_mask)
       call parallel_halo(thck_calving_front)
       call parallel_halo(interior_marine_mask)

       ! Optionally, copy interior_marine_mask to marine_interior_mask for output.
       if (present(marine_interior_mask)) then
          marine_interior_mask = interior_marine_mask
       endif

       if (present(effective_areafrac)) then
          where (ice_mask == 1 .or. land_mask == 1)
             effective_areafrac = 1.0d0
          elsewhere
             effective_areafrac = 0.0d0
          endwhere
       endif

       if (present(active_ice_mask)) then
          active_ice_mask(:,:) = ice_mask(:,:)
       endif

    endif  ! which_ho_calving_front

  end subroutine glissade_calving_front_mask

!****************************************************************************

  subroutine glissade_marine_cliff_mask(&
       nx,                     ny,                   &
       ice_mask,               floating_mask,        &
       land_mask,              active_ice_mask,      &
       marine_cliff_mask)

    ! Compute a mask to identify marine cliff cells.
    ! These are defined as cells with grounded marine ice, adjacent to ice-free ocean cells
    !  and/or inactive calving front cells.

    integer, intent(in) ::   &
         nx,  ny                  ! number of grid cells in each direction

    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,              & ! = 1 if thck > thklim, else = 0
         floating_mask,         & ! = 1 if thck > thklim and ice is floating, else = 0
         land_mask,             & ! = 1 if topg is at or above sea level, else = 0
         active_ice_mask          ! = 1 if dynamically active, else = 0

    integer, dimension(nx,ny), intent(out), optional ::  &
         marine_cliff_mask        ! = 1 if ice is grounded and marine-based and borders at least one ocean
                                  !     or inactive calving_front cell, else = 0

    !----------------------------------------------------------------
    ! Local arguments
    !----------------------------------------------------------------

    integer :: i, j


    marine_cliff_mask(:,:) = 0

    do j = 2, ny-1
       do i = 2, nx-1
          if (ice_mask(i,j) == 1 .and. land_mask(i,j) == 0 .and. floating_mask(i,j) == 0) then ! grounded marine-based ice
             if ( (land_mask(i-1,j) == 0 .and. active_ice_mask(i-1,j) == 0) .or. &  ! adjacent to inactive CF or ocean
                  (land_mask(i+1,j) == 0 .and. active_ice_mask(i+1,j) == 0) .or. &
                  (land_mask(i,j-1) == 0 .and. active_ice_mask(i,j-1) == 0) .or. &
                  (land_mask(i,j+1) == 0 .and. active_ice_mask(i,j+1) == 0) ) then
                marine_cliff_mask(i,j) = 1
             endif   ! marine cliff cell
          endif  ! grounded marine-based ice
       enddo  ! i
    enddo   ! j

    call parallel_halo(marine_cliff_mask)

  end subroutine glissade_marine_cliff_mask

!****************************************************************************

  subroutine glissade_ice_sheet_mask(nx,            ny,     &
                                     itest, jtest,  rtest,  &
                                     ice_mask,      thck,   &
                                     ice_sheet_mask,        &
                                     ice_cap_mask)

    ! Define masks that identify the ice sheet as distinct from ice caps.
    ! An ice cap is defined as a patch of ice separate from the main ice sheet.

    ! The algorithm is as follows:
    ! (1) Mark all cells with ice (ice_mask = 1) with the initial color.
    !     Mark other cells with the boundary color.
    ! (2) Seed the fill by giving the fill color to some cells that are definitely
    !     part of the ice sheet (based on thck > minthck_ice_sheet).
    ! (3) Recursively fill all cells that are connected to filled cells by a path
    !     that passes through ice-covered cells only.
    ! (4) Repeat the recursion as necessary to spread the fill to adjacent processors.
    ! (5) Once the fill is done, any cells that still have the initial color and
    !     are on land are considered to be ice caps.
    !
    ! Note: The recursive fill applies to edge neighbors, not corner neighbors.
    !        The path back to ice sheet cells must go through edges, not corners.
    !       The ice sheet seeding criterion can be changed by adjusting minthck_ice_sheet.

    integer, intent(in) :: nx, ny                  !> horizontal grid dimensions

    integer, intent(in) :: itest, jtest, rtest     !> coordinates of diagnostic point

    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask               !> = 1 if ice is present (thck > thklim)

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck                   !> ice thickness (m)

    integer, dimension(nx,ny), intent(out) ::  &
         ice_sheet_mask         !> = 1 for ice sheet cells

    integer, dimension(nx,ny), intent(out) ::  &
         ice_cap_mask           !> = 1 for ice cap cells, separately from the main ice sheet

    real(dp), parameter :: &
         minthck_ice_sheet = 2000.d0  !> thickness threshold (m) for initializing ice sheet cells

    ! local variables

    integer :: i, j, iter

    integer :: &
         max_iter,             & ! max(ewtasks, nstasks)
         local_count,          & ! local counter for filled values
         global_count,         & ! global counter for filled values
         global_count_save       ! globalcounter for filled values from previous iteration

    integer, dimension(nx,ny) ::  &
         color                  !> color variable for the fill

    logical, parameter :: verbose_ice_sheet_mask = .false.

    ! initialize
    ! Note: Ice-covered cells receive the initial color, and ice-free cells receive the boundary color.

    do j = 1, ny
       do i = 1, nx
          if (ice_mask(i,j) == 1) then
             color(i,j) = initial_color
          else
             color(i,j) = boundary_color
          endif
       enddo
    enddo

    ! Loop through cells, identifying cells that are definitely part of the ice sheet
    !  based on a threshold ice thickness. 
    ! Fill these cells and then recursively fill ice-covered neighbors.
    ! We may have to do this several times to incorporate connections between neighboring processors.

    max_iter = max(ewtasks,nstasks)
    global_count_save = 0

    do iter = 1, max_iter

       if (iter == 1) then   ! identify ice sheet cells that can seed the fill

          do j = 1, ny
             do i = 1, nx
                if (color(i,j) == initial_color .and. thck(i,j) >= minthck_ice_sheet) then
                   ! assign the fill color to this cell, and recursively fill ice-covered neighbors
                   call glissade_fill(nx,    ny,    &
                                      i,     j,     &
                                      color, ice_mask)
                endif
             enddo
          enddo

       else  ! iter > 1

          ! Check for halo cells that were just filled on neighbor processors
          ! Note: In order for a halo cell to seed the fill on this processor, it must already have the fill color.

          call parallel_halo(color)

          ! west halo layer
          i = nhalo
          do j = 1, ny
             if (color(i,j) == fill_color) then
                call glissade_fill(nx,    ny,    &
                                   i+1,   j,     &
                                   color, ice_mask)
             endif
          enddo

          ! east halo layers
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color) then
                call glissade_fill(nx,    ny,    &
                                   i-1,   j,     &
                                   color, ice_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color) then
                call glissade_fill(nx,    ny,    &
                                   i,     j+1,   &
                                   color, ice_mask)
             endif
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color) then
                call glissade_fill(nx,    ny,    &
                                   i,     j-1,   &
                                   color, ice_mask)
             endif
          enddo

       endif  ! iter = 1

       ! Count the number of filled cells.  If converged, then exit the loop.

       local_count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) local_count = local_count + 1
          enddo
       enddo

       global_count = parallel_reduce_sum(local_count)

       if (global_count == global_count_save) then
          if (verbose_ice_sheet_mask .and. main_task) &
               print*, 'Fill converged: iter, global_count =', iter, global_count
          exit
       else
          if (verbose_ice_sheet_mask .and. main_task) &
               print*, 'Convergence check: iter, global_count =', iter, global_count
          global_count_save = global_count
       endif

    enddo  ! max_iter

    ! Any cells with the fill color are considered to be part of the land-based ice sheet.
    ! Any cells with the initial color are deemed to be ice caps.
 
    ice_sheet_mask(:,:) = 0.0d0
    ice_cap_mask(:,:) = 0.0d0

    call parallel_halo(color)

    do j = 1, ny
       do i = 1, nx
          if (color(i,j) == initial_color) then
             ice_cap_mask(i,j) = 1
          elseif (color(i,j) == fill_color) then
             ice_sheet_mask(i,j) = 1
          endif
       enddo
    enddo

  end subroutine glissade_ice_sheet_mask

!****************************************************************************

  subroutine glissade_marine_connection_mask(nx,           ny,             &
                                             itest, jtest, rtest,          &
                                             thck,          topg,          &
                                             eus,           thklim,        &
                                             marine_connection_mask)

  ! Identify cells that have a marine path to the ocean.
  ! The path can include grounded marine-based ice.
  ! The ocean is identified as ice-free cells with bed elevation below sea level,
  !  or optionally with bed elevation below some threshold value.

    ! subroutine arguments

    integer, intent(in) :: nx, ny                  !> horizontal grid dimensions

    integer, intent(in) :: itest, jtest, rtest     !> coordinates of diagnostic point

    ! Note: Input thck and topg need to be correct in halo cells
    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                   & !> ice thickness (m)
         topg                      !> elevation of topography (m)

    real(dp), intent(in) :: &
         eus,                    & !> eustatic sea level (m), = 0. by default
         thklim                    !> minimum ice thickness for active cells (m)

    integer, dimension(nx,ny), intent(out) ::  &
         marine_connection_mask    !> = 1 for ocean cells and cells connected to the ocean through marine-based ice

    ! local variables

    integer, dimension(nx,ny) ::  &
         ocean_mask,          &  !> = 1 where topg - eus is below sea level and ice is absent, else = 0
         ocean_mask_temp,     &  !> temporary version of ocean_mask
         marine_mask,         &  !> marine-based cells; topg - eus < 0
         color                   ! integer 'color' mask to mark filled cells

    integer :: i, j, iter

    integer :: &
         max_iter,             & ! max(ewtasks, nstasks)
         local_count,          & ! local counter for filled values
         global_count,         & ! global counter for filled values
         global_count_save       ! globalcounter for filled values from previous iteration

    !Note: could make this a config parameter if different values are desired for different grids
    real(dp), parameter :: &
         ocean_topg_threshold = -500.d0   !> ocean threshold elevation (m) to seed the fill; negative below sea level

    logical, parameter :: verbose_marine_connection = .true.

    ! Compute ocean_mask, which is used to seed the fill.
    ! If ocean_topg_threshold was passed in, then ocean_mask includes only cells
    !  with topg - eus < ocean_topg_threshold.

    where (thck <= thklim .and. topg - eus < ocean_topg_threshold)
       ocean_mask = 1
    elsewhere
       ocean_mask = 0
    endwhere

    ! Occasionally, e.g. in Greenland fjords, a cell could be marked as ice-free ocean even though it is landlocked.
    ! To make this less likely, set ocean_mask = 0 for any cells with non-ocean neighbors.
    ! This logic could be iterated if needed.

    ocean_mask_temp = ocean_mask
    do j = 2, ny-1
       do i = 2, nx-1
          if (ocean_mask_temp(i-1,j) == 0 .or. ocean_mask_temp(i+1,j) == 0 .or.  &
              ocean_mask_temp(i,j-1) == 0 .or. ocean_mask_temp(i,j+1) == 0) then
             ocean_mask(i,j) = 0
          endif
       enddo
    enddo

    call parallel_halo(ocean_mask)

    ! initialize
    ! Compute a marine mask; = 1 for all cells with topg - eus < 0.
    ! Marine-based cells receive the initial color; land cells receive the boundary color.
    ! Seed the fill with ocean cells.

    where (ocean_mask == 1)
       marine_mask = 1
       color = fill_color
    elsewhere (topg - eus < 0.0d0)  ! marine-based cells
       marine_mask = 1
       color = initial_color
    elsewhere   ! land cells
       marine_mask = 0
       color = boundary_color
    endwhere

    if (verbose_marine_connection .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_marine_connection_mask, itest, jtest, rank =', itest, jtest, rtest
       print*, ' '
       print*, 'marine_mask'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') marine_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Loop through cells, identifying marine-based cells that border the ocean.
    ! Fill each such cell, and then recursively fill marine-based neighbor cells.
    ! We may have to do this several times to incorporate connections between neighboring processors.
    ! The result is a mask that all marine-based cells connected to the ocean are filled.

    max_iter = max(ewtasks,nstasks)
    global_count_save = 0

    do iter = 1, max_iter

       if (iter == 1) then   ! identify marine-based cells adjacent to ocean cells, which can seed the fill

          do j = 2, ny-1
             do i = 2, nx-1
                if (marine_mask(i,j) == 1) then
                   if (ocean_mask(i-1,j) == 1 .or. ocean_mask(i+1,j) == 1 .or.   &
                       ocean_mask(i,j-1) == 1 .or. ocean_mask(i,j+1) == 1) then

                      if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then

                         ! assign the fill color to this cell, and recursively fill marine-based neighbor cells
                         call glissade_fill(nx,    ny,    &
                                            i,     j,     &
                                            color, marine_mask)
                      endif
                   endif  ! adjacent to ocean
                endif  ! marine-based
             enddo  ! i
          enddo  ! j

       else  ! iter > 1

          ! Check for halo cells that were just filled on neighbor processors
          ! Note: In order for a halo cell to seed the fill on this processor, it must not only have the fill color,
          !       but also must have marine_mask = 1.

          call parallel_halo(color)

          ! west halo layer
          i = nhalo
          do j = 1, ny
             if (color(i,j) == fill_color .and. marine_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i+1,   j,     &
                                   color, marine_mask)
             endif
          enddo

          ! east halo layers
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color .and. marine_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i-1,   j,     &
                                   color, marine_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. marine_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i,     j+1,   &
                                   color, marine_mask)
             endif
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. marine_mask(i,j) == 1) then
                call glissade_fill(nx,    ny,    &
                                   i,     j-1,   &
                                   color, marine_mask)
             endif
          enddo

       endif  ! iter = 1

       ! Count the number of filled cells.  If converged, then exit the loop.

       local_count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) local_count = local_count + 1
          enddo
       enddo

       global_count = parallel_reduce_sum(local_count)

       if (global_count == global_count_save) then
          if (verbose_marine_connection .and. main_task) &
               print*, 'Fill converged: iter, global_count =', iter, global_count
          exit
       else
          if (verbose_marine_connection .and. main_task) &
               print*, 'Convergence check: iter, global_count =', iter, global_count
          global_count_save = global_count
       endif

    enddo  ! max_iter

    call parallel_halo(color)

    ! Set the marine connection mask.  This includes:
    ! (1) cells that are already ocean
    ! (2) cells with the fill color, meaning they are marine-based cells connected to the ocean

    marine_connection_mask(:,:) = 0

    do j = 1, ny
       do i = 1, nx
          if (ocean_mask(i,j) == 1 .or. color(i,j) == fill_color) then
             marine_connection_mask(i,j) = 1
          endif
       enddo
    enddo

    call parallel_halo(marine_connection_mask)

    if (verbose_marine_connection .and. this_rank == rtest) then
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
       print*, 'marine_connection_mask, rank =', this_rank
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') marine_connection_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_marine_connection_mask

!****************************************************************************

  subroutine glissade_lake_mask(nx,           ny,             &
                                itest, jtest, rtest,          &
                                floating_mask,                &
                                ocean_mask,                   &
                                lake_mask,                    &
                                ocean_connection_mask)

  ! TODO - Rewrite this subroutine with marine_connection_mask as an input.
  !        Lake cells are floating cells without a marine connection.

  ! Identify interior lake cells: cells that are floating but are not connected
  !  to the ocean along a path through other floating cells.
  ! Optionally, identify cells with an ocean connection: either ice-free ocean,
  !  or floating and connected through other floating cells to the ocean.
  ! Note: The path to the ocean must pass through edge neighbors, not corner neighbors.

    integer, intent(in) :: nx, ny                  !> horizontal grid dimensions

    integer, intent(in) :: itest, jtest, rtest     !> coordinates of diagnostic point

    ! Note: The calling subroutine decides how to define floating_mask.
    !       E.g., it could be defined based on f_ground_cell with which_ho_ground = 2.
    ! Each input mask needs to be correct in halo cells

    integer, dimension(nx,ny), intent(in) ::  &
         floating_mask,          & !> = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask                !> = 1 where topg - eus is below sea level and ice is absent, else = 0

    integer, dimension(nx,ny), intent(out) ::  &
         lake_mask                 !> = 1 for floating cells disconnected from the ocean, else = 0

    integer, dimension(nx,ny), intent(out), optional ::  &
         ocean_connection_mask     !> = 1 for ocean cells, and cells connected to the ocean through floating ice

    ! local variables

    integer, dimension(nx,ny) ::  &
         color,                  & ! integer 'color' mask to mark filled cells
         border_mask               ! = 1 for grounded marine ice adjacent to ocean-connected cells

    integer :: i, j, iter

    integer :: &
         max_iter,             & ! max(ewtasks, nstasks)
         local_count,          & ! local counter for filled values
         global_count,         & ! global counter for filled values
         global_count_save       ! globalcounter for filled values from previous iteration

    logical, parameter :: verbose_lake = .false.

    integer :: ig, jg

    if (verbose_lake .and. this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_lake_mask, itest, jtest, rank =', itest, jtest, rtest
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
    ! The result is a mask that identifies (with the fill color) all floating cells connected to the ocean.

    max_iter = max(ewtasks,nstasks)
    global_count_save = 0

    do iter = 1, max_iter

       if (iter == 1) then   ! identify floating cells adjacent to ocean cells, which can seed the fill

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
          !       but also must have floating_mask = 1.

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

       endif  ! iter = 1

       ! Count the number of filled cells.  If converged, then exit the loop.

       local_count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) local_count = local_count + 1
          enddo
       enddo

       global_count = parallel_reduce_sum(local_count)

       if (global_count == global_count_save) then
          if (verbose_lake .and. main_task) &
               print*, 'Fill converged: iter, global_count =', iter, global_count
          exit
       else
          if (verbose_lake .and. main_task) &
               print*, 'Convergence check: iter, global_count =', iter, global_count
          global_count_save = global_count
       endif

    enddo  ! max_iter

    call parallel_halo(color)

    ! Identify lake cells: floating cells that still have the initial color.

    lake_mask(:,:) = 0

    do j = 1, ny
       do i = 1, nx
          if (color(i,j) == initial_color .and. floating_mask(i,j) == 1) then
             lake_mask(i,j) = 1

             if (verbose_lake .and. this_rank == rtest) then
                call parallel_globalindex(i, j, ig, jg)
                print*, 'Lake cell: task, i, j, ig, jg =', this_rank, i, j, ig, jg
             endif

          endif
       enddo
    enddo

    call parallel_halo(lake_mask)

    if (verbose_lake .and. this_rank == rtest) then
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

    if (present(ocean_connection_mask)) then

       ! Identify cells connected to the ocean.  This includes:
       ! (1) cells that are already ocean
       ! (2) floating cells with the fill color, meaning they are connected to the ocean via other floating cells

       ocean_connection_mask(:,:) = 0

       do j = 1, ny
          do i = 1, nx
             if (ocean_mask(i,j) == 1 .or. color(i,j) == fill_color) then
                ocean_connection_mask(i,j) = 1
             endif
          enddo
       enddo

       call parallel_halo(ocean_connection_mask)

       if (verbose_lake .and. this_rank == rtest) then
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
          print*, 'ocean_connection_mask, rank =', this_rank
          do j = jtest+3, jtest-3, -1
             write(6,'(i6)',advance='no') j
             do i = itest-3, itest+3
                write(6,'(i10)',advance='no') ocean_connection_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

    endif  ! present(ocean_connection_mask)

  end subroutine glissade_lake_mask

!****************************************************************************

  subroutine glissade_extend_mask(nx,       ny,  &
                                  input_mask,      &
                                  extended_mask)

    ! Compute a mask that includes
    ! (1) cells with input_mask = 1
    ! (2) cells that are adjacent to cells with input_mask = 1
    ! For now, assume that both edge and diagonal neighbors are adjacent.
    ! If needed, could add an option to choose only edge neighbors.

    integer, intent(in) ::   &
         nx,  ny                !> number of grid cells in each direction

    integer, dimension(nx,ny), intent(in) :: &
         input_mask             !> input mask to be extended

    integer, dimension(nx,ny), intent(out) :: &
         extended_mask          !> input mask extended by adding neighbor cells

    ! local variables

    integer :: i, j

    extended_mask(:,:) = 0

    do j = 2, ny-1
       do i = 2, nx-1
          if (input_mask(i,j) == 1) then
             extended_mask(i,j) = 1
          elseif (input_mask(i,j-1)   == 1 .or. input_mask(i,j+1)   == 1   .or.  &
                  input_mask(i-1,j+1) == 1 .or. input_mask(i+1,j+1) == 1 .or.  &
                  input_mask(i-1,j)   == 1 .or. input_mask(i+1,j)   == 1 .or.  &
                  input_mask(i-1,j-1) == 1 .or. input_mask(i+1,j-1) == 1) then
             extended_mask(i,j) = 1
          endif
       enddo
    enddo

    call parallel_halo(extended_mask)

  end subroutine glissade_extend_mask

!****************************************************************************

  recursive subroutine glissade_fill(nx,  ny,         &
                                     i,   j,          &
                                     color,           &
                                     fill_mask)

    ! Given a cell (i,j), determine whether it should be given the fill color
    !  and recursively fill neighbor cells.
    ! This subroutine differs from the subroutine above in that cell (i,j) is filled and
    !  the subroutine is called recursively only if fill_mask = 1.
    ! In the subroutine above, cell (i,j) can be filled when active_ice_mask = 0,
    !  but the subroutine is called recursively only if active_ice_mask = 1.

    integer, intent(in) :: nx, ny                       !> domain size
    integer, intent(in) :: i, j                         !> horizontal indices of current cell

    integer, dimension(nx,ny), intent(inout) :: &
         color                                          !> color (initial, fill or boundary)

    integer, dimension(nx,ny), intent(in) :: &
         fill_mask                                      !> = 1 if the cell satisfies the fill criterion

    if (color(i,j) /= fill_color .and. color(i,j) /= boundary_color .and. fill_mask(i,j) == 1) then

       ! assign the fill color to this cell
       color(i,j) = fill_color

       ! recursively call this subroutine for each neighbor to see if it should be filled
       !TODO - May want to rewrite this to avoid recursion, which can crash the code when
       !       the recursion stack is very large on fine grids.
       if (i > 1)  call glissade_fill(nx,    ny,  &
                                      i-1,   j,   &
                                      color, fill_mask)

       if (i < nx) call glissade_fill(nx,    ny,  &
                                      i+1,   j,   &
                                      color, fill_mask)

       if (j > 1)  call glissade_fill(nx,    ny, &
                                      i,     j-1, &
                                      color, fill_mask)

       if (j < ny) call glissade_fill(nx,    ny,  &
                                      i,     j+1, &
                                      color, fill_mask)

    endif   ! not fill color or boundary color

  end subroutine glissade_fill

!****************************************************************************

  recursive subroutine glissade_fill_with_buffer(&
       nx,  ny,         &
       i,   j,          &
       color,           &
       fill_mask)

    ! Given a domain with an initial color, a boundary color and a fill color,
    !  recursively assign the fill color to all cells that are connected to cells
    !  with the fill color.  The connection must pass through cell edges (not vertices).
    ! Note: "with_buffer" refers to the idea that we fill not only cells with fill_mask = 1,
    !       but also cells with fill_mask = 0, provided the cell does not already have
    !       the fill color or boundary color.  But only cells with fill_mask = 1 result
    !       in additional filling.
    !       This logic is used, for example, with iceberg removal, where fill_mask = active_ice_mask.
    !       In this case we fill not only active cells, but also a buffer layer of inactive cells.
    !       Thus, both active cells and inactive buffer cells are filled and are spared from removal.

    integer, intent(in) :: nx, ny                       !> domain size
    integer, intent(in) :: i, j                         !> horizontal indices of current cell

    integer, dimension(nx,ny), intent(inout) :: &
         color                                          !> color (initial, fill or boundary)

    integer, dimension(nx,ny), intent(in) :: &
         fill_mask                                      !> = 1 if the cell satisfies the fill criterion

    if (color(i,j) /= fill_color .and. color(i,j) /= boundary_color) then

       ! assign the fill color to this cell
       color(i,j) = fill_color

       ! If fill_mask = 1, then fill this cell but do not call the subroutine recursively

       if (fill_mask(i,j) == 0) return   ! skip the recursion

       ! recursively call this subroutine for each neighbor to see if it should be filled
       !TODO - May want to rewrite this to avoid recursion, which can crash the code when
       !       the recursion stack is very large on fine grids.
       if (i > 1)  call glissade_fill_with_buffer(nx,    ny,  &
                                                  i-1,   j,   &
                                                  color, fill_mask)

       if (i < nx) call glissade_fill_with_buffer(nx,    ny,  &
                                                  i+1,   j,   &
                                                  color, fill_mask)

       if (j > 1)  call glissade_fill_with_buffer(nx,    ny, &
                                                  i,     j-1, &
                                                  color, fill_mask)

       if (j < ny) call glissade_fill_with_buffer(nx,    ny,  &
                                                  i,     j+1, &
                                                  color, fill_mask)

    endif   ! not fill color or boundary color

  end subroutine glissade_fill_with_buffer

!****************************************************************************

  end module glissade_masks

!****************************************************************************

