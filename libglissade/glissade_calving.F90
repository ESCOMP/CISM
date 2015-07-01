!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_calving.F90 - part of the Community Ice Sheet Model (CISM)  
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
!CALVING TODO:
! (1) Test Glide v. Glissade masks (then remove Glide masks)
! (2) Scale eus
! (3) Transport damage tracer

!!#ifdef HAVE_CONFIG_H
!!#include "config.inc"
!!#endif
#include "glide_mask.inc"

module glissade_calving

  use glide_types
  use glimmer_global, only: dp
  use parallel

  implicit none

  !WHL - debug
  logical, parameter :: verbose_calving = .true.
  integer, parameter :: jtest = 3

contains
!-------------------------------------------------------------------------------  

  subroutine glissade_calve_ice(which_calving,                  &
                                thck,         relx,             &    
                                topg,          &
                                thkmask,       &  !TODO - Remove Glide thkmask
                                marine_limit, calving_fraction, &    
                                eus,          calving_thck)

    ! Calve ice according to one of several alternative methods

    use glissade_masks, only: glissade_get_masks
    use glimmer_paramets, only: thk0

    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer,  intent(in)                    :: which_calving     !> calving option
    real(dp), dimension(:,:), intent(inout) :: thck              !> ice thickness
    real(dp), dimension(:,:), intent(in)    :: relx              !> relaxed bedrock topography
    real(dp), dimension(:,:), intent(in)    :: topg              !> present bedrock topography
    integer,  dimension(:,:), intent(in)    :: thkmask           !> grid type mask
    real(dp), intent(in)                    :: marine_limit      !> lower limit on topography elevation for ice to be present
    real(dp), intent(in)                    :: calving_fraction  !> fraction of ice lost when calving; used with which_ho_calving = 2
    real(dp), intent(in)                    :: eus               !> eustatic sea level
    real(dp), dimension(:,:), intent(out)   :: calving_thck      !> thickness lost due to calving in each grid cell

    integer :: nx, ny   ! grid dimensions
    integer :: i, j

    integer,  dimension(:,:), allocatable   ::  &
         ice_mask,      &  !> = 1 where ice is present, else = 0
         floating_mask, &  !> = 1 where ice is floating, else = 0
         ocean_mask        !> = 1 where topg is below sea level and ice is absent, else = 0

    !---------------------------------------------------------------------
   
    calving_thck(:,:) = 0.d0

    nx = size(thck,1)
    ny = size(thck,2)

    allocate (ice_mask(nx,ny))
    allocate (floating_mask(nx,ny))
    allocate (ocean_mask(nx,ny))

    ! get masks
    ! Note: Use thickness limit of 0.0 instead of thklim, so as to remove
    !       all ice that meets the calving criteria.

    call glissade_get_masks(nx,          ny,         &
                            thck,        topg,       &
                            eus,         0.0d0,      &    
                            ice_mask,    floating_mask, &
                            ocean_mask)

    !WHL - debug
    if (verbose_calving) then
       j = jtest
       print*, 'Starting glissade_calve_ice: j =', j
       print*, 'i, relx, topg, thck, usfc:'
       do i = nx-11, nx-1
          print*, i, relx(i,j)*thk0, topg(i,j)*thk0, thck(i,j)*thk0, (topg(i,j) + thck(i,j))*thk0
       enddo
    endif

    select case (which_calving)

    case(CALVING_NONE)           ! do nothing

        
    case(CALVING_FLOAT_ZERO)     ! set thickness to zero if ice is floating

       where (floating_mask == 1)
          calving_thck = thck
          thck = 0.0d0
       end where

    case(CALVING_FLOAT_FRACTION) ! remove a specified fraction of ice when floating

       !WHL - Changed definition of calving fraction; now it is the fraction lost
       !      rather than the fraction remaining

       ! Glide calving mask: (1) Ice is present, (2) ocean cell (topg < 0), (3) ice-ocean margin
       ! Equivalent glissade mask: Floating ice is present and borders an ocean cell

       do j = 2, ny-1
          do i = 2, nx-1
             if (floating_mask(i,j) == 1) then
                if (ocean_mask(i-1,j) == 1 .or. ocean_mask(i+1,j) == 1 .or.  &
                    ocean_mask(i,j-1) == 1 .or. ocean_mask(i,j+1) == 1) then 
                   calving_thck(i,j) = calving_fraction * thck(i,j)
                   thck(i,j) = thck(i,j) - calving_thck(i,j)
                endif
             endif
          enddo
       enddo
      
    case(CALVING_RELX_THRESHOLD)   ! set thickness to zero if relaxed bedrock is below a given level

       !WHL - Not sure why this one does not check to see if we are at the ice-ocean margin
       !      (like TOPG_THRESHOLD below)
  
       where (relx <= marine_limit+eus)
          calving_thck = thck
          thck = 0.0d0
       end where
       
    case(CALVING_TOPG_THRESHOLD)   ! set thickness to zero at marine edge if present bedrock is below a given level

       !TODO - Replace Glide marine ice mask
       ! Glide marine_ice_edge mask: (1) Ice is present (but not necessarily floating) 
       !                             (2) at the ocean margin
       !                             (3) topg below prescribed level

       do j = 2, ny-1
          do i = 2, nx-1
             if (topg(i,j) < marine_limit+eus) then
                if (ocean_mask(i-1,j) == 1 .or. ocean_mask(i+1,j) == 1 .or.  &
                    ocean_mask(i,j-1) == 1 .or. ocean_mask(i,j+1) == 1) then 
                   calving_thck(i,j) = thck(i,j)
                   thck(i,j) = 0.0d0
                endif
             endif
          enddo
       enddo

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
             calving_thck = thck
             thck = 0.0d0
          end where
       elseif (eus*thk0 <= -80.d0) then
          where (relx*thk0 <= (2.d0*eus*thk0 - 0.25d0*(eus*thk0 + 80.d0)**2.d0))
             calving_thck = thck
             thck = 0.0d0
          end where
       end if
       
    case(CALVING_DAMAGE)   ! remove ice that is at or near the marine boundary and is sufficiently damaged

       !TODO - Call a suitable subroutine developed by Jeremy B. et al.

    end select

    !WHL - debug
    if (verbose_calving) then
       j = jtest
       print*, 'Done in glissade_calve_ice: j =', j
       print*, 'i, relx, topg, thck, usfc:'
       do i = nx-11, nx-1
          print*, i, relx(i,j)*thk0, topg(i,j)*thk0, thck(i,j)*thk0, (topg(i,j) + thck(i,j))*thk0
       enddo
    endif

  end subroutine glissade_calve_ice

!---------------------------------------------------------------------------

end module glissade_calving

!---------------------------------------------------------------------------
