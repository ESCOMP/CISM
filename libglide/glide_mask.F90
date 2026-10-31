!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_mask.F90 - part of the Community Ice Sheet Model (CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glide_mask

    ! masking ice thicknesses

    use glimmer_global, only : dp

    implicit none

contains

!TODO - Write a new subroutine (in addition to glide_set_mask) to compute mask for staggered grid?
!       This subroutine is now called from glissade_velo_driver with stagthck and stagtopg
!       as input arguments.

  subroutine glide_set_mask(numerics, thck, topg, ewn, nsn, eus, mask)

    use glide_types
    use glimmer_physcon, only : rhoi, rhoo
    implicit none

    type(glide_numerics), intent(in) :: numerics !Numerical parameters structure
    real(dp), dimension(:,:), intent(in) :: thck ! Ice thickness
    real(dp), dimension(:,:), intent(in) :: topg ! Bedrock topography (not lower surface!)
    integer, intent(in) :: ewn, nsn              ! Grid size
    real(dp), intent(in) :: eus                  ! Sea level
    integer, dimension(:,:), intent(inout) :: mask   ! Output mask

    ! local variables
    integer ew,ns

    !Note - This array may not be needed, at least in parallel.

    ! Create an array to "fake" the boundaries of the mask so that boundary
    ! finding can work even on the boundaries of the real mask.

    integer, dimension(0:ewn+1,0:nsn+1) :: maskWithBounds;

    mask = 0

!Note - This mask is confusing.  Wondering if we should replace it by a series of logical masks.

! Would need the following:
! glide_mask_has_ice = 1
! glide_mask_thin_ice = 3
! glide_mask_ocean = 4 (below sea level, with or without ice)
! glide_mask_land = 8 (complement of glide_mask_ocean)
! glide_mask_grounding_line = 16 (could define in terms of margin and has ice?)
! glide_mask_margin = 32 (has_ice + at least one neighbor with no ice)
! glide_mask_dirichlet_bc = 64
! glide_mask_comp_domain_bnd = 128 (no longer needed with new global BC?)
! glide_no_ice (complement of glide_has_ice)
! glide_is_thin
! glide_is_ocean (ocean + no_ice; change to glide_ocean_icefree or remove?)
! glide_is_land (land + no_ice; change to glide_land_icefree or remove?)
! glide_is_ground (land + has_ice)
! glide_is_float (ocean + has_ice)
! glide_is_grounding_line (just inside or just outside? Used only in glide_ground)
! glide_is_margin
! glide_is_land_margin (margin + land + has_ice)
! glide_is_calving (margin + ocean + has_ice; change the name to is_marine_margin?)
! glide_is_marine_ice_edge (margin + (float or GL); may not be needed)
! glide_is_dirichlet_boundary 
! glide_is_comp_domain_bnd (may not be needed with new global BC?)
! 
!       If we keep the present structure, could change glide_is_land to glide_icefree_land,
!                                                      glide_is_ocean to glide_icefree_ocean
!       Could get by with fewer masks in the code by removing some combinations
!       Could remove *BITS

    !Identify points with any ice
    where (thck > 0.d0)
        mask = ior(mask, GLIDE_MASK_HAS_ICE)  ! GLIDE_MASK_HAS_ICE = 1; see glide_mask.inc
    endwhere

    !Identify points where the ice is below the ice dynamics limit
    where (thck > 0.d0 .and. thck < numerics%thklim)
        mask = ior(mask, GLIDE_MASK_THIN_ICE)  ! GLIDE_MASK_THIN_ICE = 3
    endwhere

    !Identify points where the ice is floating or where there is open ocean
    where (topg - eus < (-rhoi/rhoo) * thck)
        mask = ior(mask, GLIDE_MASK_OCEAN)   ! GLIDE_MASK_OCEAN = 8
    elsewhere
        mask = ior(mask, GLIDE_MASK_LAND)    ! GLIDE_MASK_LAND = 4
    endwhere

    !TODO - Replace the following with a halo call for 'mask', with appropriate global BC?
    
    maskWithBounds = 0
    maskWithBounds(1:ewn, 1:nsn) = MASK
    maskWithBounds(0,1:nsn) = mask(1,:)
    maskWithBounds(1:ewn,0) = mask(:,1)
    maskWithBounds(ewn+1,1:nsn) = mask(ewn,:)
    maskWithBounds(1:ewn,nsn+1) = mask(:,nsn)
    maskWithBounds(0,0) = mask(1,1)
    maskWithBounds(ewn+1,nsn+1) = mask(ewn,nsn)
    maskWithBounds(0,nsn+1) = mask(1,nsn)
    maskWithBounds(ewn+1,0) = mask(ewn,1)

    ! finding boundaries

    !Note: If halo cells are present, maskWithBounds array may not be needed; can replace with mask array.
    !      Not sure what happens here when we're computing a mask on the velocity grid.

    do ns = 1,nsn
       do ew = 1,ewn
          !Find the grounding line
          if (GLIDE_IS_GROUND(MASK(ew,ns))) then    ! land + has_ice
             if (GLIDE_IS_FLOAT(maskWithBounds(ew-1,ns)) .or. &
                  GLIDE_IS_FLOAT(maskWithBounds(ew+1,ns)) .or. &
                  GLIDE_IS_FLOAT(maskWithBounds(ew,ns-1)) .or. & 
                  GLIDE_IS_FLOAT(maskWithBounds(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_GROUNDING_LINE)
             end if
          end if

          ! Ice margin
          ! *tb* A point is now masked even if it touches the ocean on one corner.
          if ( GLIDE_HAS_ICE(mask(ew, ns)) .and. &
              (GLIDE_NO_ICE(maskWithBounds(ew-1,ns))   .or. GLIDE_NO_ICE(maskWithBounds(ew+1,ns))   .or. &
               GLIDE_NO_ICE(maskWithBounds(ew,ns-1))   .or. GLIDE_NO_ICE(maskWithBounds(ew,ns+1))   .or. &
               GLIDE_NO_ICE(maskWithBounds(ew-1,ns-1)) .or. GLIDE_NO_ICE(maskWithBounds(ew-1,ns+1)) .or. &
               GLIDE_NO_ICE(maskWithBounds(ew+1,ns-1)) .or. GLIDE_NO_ICE(maskWithBounds(ew+1,ns+1)))) then
             MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_MARGIN)
          end if

       end do
    end do

  end subroutine glide_set_mask

  subroutine augment_kinbc_mask(mask, kinbcmask)

    !  Augments the Glide mask with the location of kinematic (dirichlet) boundary
    !   conditions.  These locations cannot be determined by the model a priori, and
    !   must be specified through a field in a NetCDF file.
    integer, dimension(:,:), target :: mask
    integer, dimension(:,:) :: kinbcmask

    integer, dimension(:,:), pointer :: maskp

    !Because the kinematic boundary conditions are specified on the staggered grid,
    !there may be a size mismatch here depending on whether we are computing a mask
    !for the staggered grid.
    if (size(mask, 1) /= size(kinbcmask, 1)) then
        maskp => mask(1:size(mask,1) - 1, 1:size(mask,2) - 1)
    else
        maskp => mask
    end if

    where (kinbcmask /= 0)
        maskp = ior(maskp, GLIDE_MASK_DIRICHLET_BC)
    endwhere
  end subroutine augment_kinbc_mask

end module glide_mask
