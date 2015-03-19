!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_outputs.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glad_outputs

  ! This module defines routines for computing the outputs that CISM sends to a climate
  ! model.

  implicit none
  private

  public :: set_output_fields  ! set all fields output to a climate model

contains

  subroutine set_output_fields(instance, &
       ice_covered, topo, rofi, rofl, hflx, &
       ice_sheet_grid_mask, icemask_coupled_fluxes)

    ! Arguments ----------------------------------------------------------------------------

    type(glint_instance), intent(in) :: instance
    real(dp),dimension(:,:),intent(out) :: ice_covered  ! whether each grid cell is ice-covered [0,1]
    real(dp),dimension(:,:),intent(out) :: topo         ! output surface elevation (m)
    real(dp),dimension(:,:),intent(out) :: hflx         ! output heat flux (W/m^2, positive down)
    real(dp),dimension(:,:),intent(out) :: rofi         ! output ice runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: rofl         ! output liquid runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: ice_sheet_grid_mask !mask of ice sheet grid coverage
    real(dp),dimension(:,:),intent(out) :: icemask_coupled_fluxes !mask of ice sheet grid coverage where we are potentially sending non-zero fluxes

    ! Internal variables ----------------------------------------------------------------------

    integer :: nxl, nyl  ! local grid dimensions
    integer :: i, j      ! indices

    ! Begin subroutine code -------------------------------------------------------------------
    
    ! Initialize arrays. This shouldn't be necessary (because, if the below code is
    ! correct, all points should be explicitly assigned some value), but adds some safety
    ! in case any bugs creep into the below code.
    ice_covered(:,:) = 0.d0
    topo(:,:) = 0.d0
    hflx(:,:) = 0.d0
    rofi(:,:) = 0.d0
    rofl(:,:) = 0.d0
    ice_sheet_grid_mask(:,:) = 0.d0
    icemask_coupled_fluxes(:,:) = 0.d0
    
    nxl = instance%lgrid%size%pt(1)
    nyl = instance%lgrid%size%pt(2)

    do j = 1, nyl
       do i = 1, nxl
          if (is_in_active_grid(instance%model%geometry, i, j)) then
             ice_sheet_grid_mask(i,j) = 1.d0

             if (is_ice_covered(instance%model%geometry, i, j)) then
                ice_covered(i,j) = 1.d0
             else
                ice_covered(i,j) = 0.d0
             end if

             ! Note that we use the same method for computing topo whether this point is
             ! ice-covered or ice-free. This is in contrast to the method for computing
             ! ice-free topo in glint_upscaling_gcm.
             topo(i,j) = thk0 * instance%model%geometry%usrf(i,j)
             
          else
             ! Note that this logic implies that if (in theory) we had an ice-covered
             ! point outside the "active grid", it will get classified as ice-free for
             ! these purposes. This mimics the logic currently in glint_upscaling_gcm.
             ice_sheet_grid_mask(i,j) = 0.d0
             ice_covered(i,j) = 0.d0
             topo(i,j) = 0.d0
          end if

       end do
    end do

    ! TODO(wjs, 2015-03-18) Set hflx, rofi & rofl. Note that these need some
    ! time-averaging, as is done in glint_upscale.F90: glint_upscaling_gcm. For now, I am
    ! simply setting them to 0.
    hflx(:,:) = 0.d0
    rofi(:,:) = 0.d0
    rofl(:,:) = 0.d0

    if (instance%zero_gcm_fluxes == ZERO_GCM_FLUXES_TRUE) then
       icemask_coupled_fluxes(:,:) = 0.d0
       hflx(:,:) = 0.d0
       rofi(:,:) = 0.d0
       rofl(:,:) = 0.d0
    else
       icemask_coupled_fluxes(:,:) = ice_sheet_grid_mask(:,:)
    end if
    
  end subroutine set_output_fields

  
  !===================================================================

  logical function is_in_active_grid(geometry, i, j)
    ! Return true if the given point is inside the "active grid". The active grid includes
    ! any point that can receive a positive surface mass balance, which includes any
    ! point classified as land or ice sheet.
    type(glide_geometry), intent(in) :: geometry
    integer, intent(in) :: i, j  ! point of interest

    real(dp) :: usrf     ! surface elevation (m)

    ! TODO(wjs, 2015-03-18) Could the logic here be replaced by the use of some existing
    ! mask? For now I am simply re-implementing the logic that was in glint.

    usrf = thk0 * instance%model%geometry%usrf(i,j)

    if (usrf > 0.d0) then
       ! points not at sea level are assumed to be land or ice sheet
       is_in_active_grid = .true.
    else
       is_in_active_grid = .false.
    end if
  
  !===================================================================

  logical function is_ice_covered(geometry, i, j)
    ! Return true if the given point is ice-covered
    type(glide_geometry), intent(in) :: geometry
    integer, intent(in) :: i, j  ! point of interest

    real(dp) :: thck     ! ice thickness (m)

    ! TODO(wjs, 2015-03-18) The logic here should probably be replaced by the use of some
    ! existing mask. For now I am simply re-implementing the logic that was in glint.

    thck = thk0 * instance%model%geometry%thck(i,j)

    if (thck > min_thck) then
       is_ice_covered = .true.
    else
       is_ice_covered = .false.
    end if

  end function is_ice_covered
