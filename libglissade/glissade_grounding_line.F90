!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grounding_line.F90 - part of the Community Ice Sheet Model (CISM)  
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
!
! This module contains routines for computing grounding-line fields and diagnostics
! for the Glissade solver.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_grounding_line

    use glimmer_global, only: dp
    use glimmer_physcon, only: rhoi, rhoo
    use glide_types  ! grounding line options
    use parallel

    implicit none

    ! All subroutines in this module are public

  contains

!****************************************************************************

  subroutine glissade_grounded_fraction(nx,            ny,                      &
                                        thck,          topg,                    &
                                        eus,           ice_mask,                &
                                        floating_mask, land_mask,               &
                                        whichground,   whichflotation_function, &
                                        f_ground,      f_flotation)

    !----------------------------------------------------------------
    ! Compute fraction of ice that is grounded.
    ! This fraction is computed at vertices based on the thickness and
    !  topography of the four neighboring cell centers.
    !
    ! There are three options for computing the grounded fraction, based on the value of whichground:
    ! (0) HO_GROUND_NO_GLP: f_ground = 1 for all vertices with grounded, ice-covered neighbor cells
    ! (1) HO_GROUND_GLP: 0 <= f_ground <= 1 based on grounding-line parameterization
    ! (2) HO_GROUND_ALL: f_ground = 1 for all vertices with ice-covered neighbor cells
    !
    ! Notes on whichground:
    !       Both (0) and (1) rely on computing a flotation function at cell centers.
    !       - Method (0) stipulates that a vertex is grounded if any cell neighbor is grounded,
    !         based on the value of the flotation function.
    !       - Method (1) interpolates the flotation function over the bounding box of each
    !         vertex and analytically integrates to compute the grounded and floating fractions.
    !
    ! In addition, there are three options for the flotation function that is interpolated 
    ! from cell centers to vertices as part of the computation of f_ground:
    ! (0) HO_FLOTATION_FUNCTION_PATTYN: f_flotation = (-rhow*b/rhoi*H) = f_pattyn; <=1 for grounded, > 1 for floating
    ! (1) HO_FLOTATION_FUNCTION_INVERSE_PATTYN: f_flotation = (rhoi*H)/(-rhow*b) = 1/f_pattyn; >=1 for grounded, < 1 for floating
    ! (2) HO_FLOTATION_FUNCTION_LINEAR: f_flotation = -rhow*b - rhoi*H = ocean cavity thickness; <=0 for grounded, > 0 for floating
    ! These apply to whichground = 0 or 1; they are irrelevant for whichground = 2.
    !
    ! Notes on whichflotation_function:
    !       Results can be sensitive to the choice of flotation function.
    !       For instance, one such function, used by Pattyn et al. (2006) and Gladstone (2010), is 
    !           f = (-rhow*b) / (rhoi*H)
    !       This function generally works well, but consider a case where one floating cell
    !        (in a fjord, say) is surrounded by grounded cells.  We can have f ~ 10 for the
    !        floating cell and f ~ 0.5 for the grounded cells, so the average > 1 and the vertex
    !        is deemed to be floating with beta = 0.  This can lead to excessive velocities.
    !
    !       Another choice is the inverse of the first:
    !           f = (rhoi*H) / (-rhow*b)
    !       This function is negative for land-based ice (b > 0) and grows without bound as b -> 0.
    !        So it needs to be capped at a large positive value if b > 0 or is small and negative.
    !        This function has the virtue that strongly grounded cells have f >> 1, while floating
    !        cells have 0 < f < 1. So if there are 1 or 2 strongly grounded cells at a vertex, 
    !        then the vertex is deemed to be grounded (or mostly grounded), and excessive velocites 
    !        are less likely.
    !
    !       A third choice (suggested by Xylar Asay-Davis) is the following:
    !           f = -rhoi*b - rhoi*H
    !       If positive, this function is equal to the thickness of the ocean cavity beneath floating ice.
    !       If negative, this function implies that the ice is grounded.
    !       
    !       Currently (as of Sept. 2015), the inverse Pattyn function (1) is the default, but the user can
    !        change this by setting which_ho_flotation_function in the config file.
    !
    !       For the GLP, f_flotation is first computed in all active ice-covered cells.
    !       Then it is extrapolated to adjacent ice-free cells.  Thus, f_flotation has a physically
    !        sensible value (either computed directly, or extrapolated from a neighbor) in all four 
    !        cells surrounding each active vertex. (By definition, an active vertex is a vertex with
    !        at least one active ice-covered neighbor.) Thus, we can interpolate f_flotation
    !        within the bounding box around each active vertex to compute f_ground at the vertex.
    !----------------------------------------------------------------
    
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny                ! number of grid cells in each direction

    ! Default dimensions are meters, but this subroutine will work for
    ! any units as long as thck and topg have the same units.

    real(dp), dimension(nx,ny), intent(in) ::  &
       thck,                 &! ice thickness (m)
       topg                   ! elevation of topography (m)

    real(dp), intent(in) :: &
       eus                    ! eustatic sea level (= 0 by default)

    integer, dimension(nx,ny), intent(in) ::   &
       ice_mask,            & ! = 1 if ice is present (thck > thklim), else = 0
       floating_mask,       & ! = 1 if ice is present (thck > thklim) and floating, else = 0
       land_mask              ! = 1 if topg is at or above sea level

    ! see comments above for more information about these options
    integer, intent(in) ::     &
       whichground,            &! option for computing f_ground
       whichflotation_function  ! option for computing f_flotation

    real(dp), dimension(nx-1,ny-1), intent(out) ::  &
       f_ground               ! grounded ice fraction at vertex, 0 <= f_ground <= 1
                              ! set to special value where vmask = 0

    real(dp), dimension(nx,ny), intent(out) :: &
       f_flotation            ! flotation function
                              ! originally f_pattyn = (-rhow*b) / (rhoi*thck), but added two more options
                              ! see comments above

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
           
    integer :: i, j, ii, jj

    integer, dimension(nx-1,ny-1) ::   &
       vmask                     ! = 1 for vertices neighboring at least one cell where ice is present, else = 0

    real(dp), dimension(nx,ny) :: &
       f_flotation_extrap        ! f_flotation, extrapolated to cells without active ice

    real(dp) :: a, b, c, d       ! coefficients in bilinear interpolation
                                 ! f(x,y) = a + b*x + c*y + d*x*y

    real(dp) :: f1, f2, f3, f4   ! f_flotation at different cell centers

    real(dp) ::  &
       var,                 &! combination of f_flotation terms that determines regions to be integrated
       f_flotation_vertex    ! f_flotation interpolated to vertex

    integer :: nfloat     ! number of grounded vertices of a cell (0 to 4)

    logical, dimension(nx,ny) :: &
       cground,          & ! true if a cell is land and/or has grounded ice, else = false
       cfloat              ! true if flotation condition is satisfied at cell center, else = false
                           ! Note: cground and cfloat are not exact complements; both are false for ice-free ocean

    logical, dimension(2,2) ::   &
       logvar              ! set locally to float or .not.float, depending on nfloat

    real(dp) ::   &
       f_corner, &         ! fractional area in a corner region of the cell
       f_corner1, f_corner2,  &
       f_trapezoid         ! fractional area in a trapezoidal region of the cell

    logical :: adjacent   ! true if two grounded vertices are adjacent (rather than opposite)

    real(dp), parameter :: &
       eps06 = 1.d-06     ! small number

    !WHL - debug
    integer, parameter :: it = 1, jt = 1, rtest = -9999

    logical :: filled     ! true if f_flotation has been filled by extrapolation

    !TODO - Test MISMIP sensitivity to the value of this cap
    ! This needs to be big enough that one land-based cell at a vertex is sufficient to ground the others
    real(dp), parameter :: grounding_factor_max = 10.d0  

    !----------------------------------------------------------------
    ! Compute ice mask at vertices (= 1 if any surrounding cells have ice)
    !----------------------------------------------------------------

    do j = 1, ny-1
       do i = 1, nx-1
          if (ice_mask(i,j+1)==1 .or. ice_mask(i+1,j+1)==1 .or.   &
              ice_mask(i,j)  ==1 .or. ice_mask(i+1,j)  ==1 ) then
             vmask(i,j) = 1
          else
             vmask(i,j) = 0
          endif
       enddo
    enddo

    ! Compute flotation function at cell centers
    ! For diagnostic purposes, f_flotation is always computed, although it affects f_ground
    !  (and thus the velocities) only when running with a GLP.

    ! Note: All cells (including ice-free cells) are initialized to be strongly grounded on land
    !        and almost floating (i.e., barely grounded) in the ocean. Since no cells are initialized
    !        to be floating, all cells are initialized to have cfloat = F in the GLP.
    !       However, these initial values are never used in the calculation of f_ground.
    !       All values used in the f_ground calculation come from cells with ice, or cells adjacent to cells with ice.

    if (whichflotation_function == HO_FLOTATION_FUNCTION_PATTYN) then  ! grounded if f_flotation <= 1, else floating

       ! initialize to a large negative value (strongly grounded) on land, and 1.0 (almost floating) in the ocean
       do j = 1, ny
          do i = 1, nx
             if (land_mask(i,j) == 1) then
                f_flotation(i,j) = -grounding_factor_max
             else    ! marine cell
                f_flotation(i,j) = 1.0d0
             endif
          enddo
       enddo

       ! set f_flotation in ice-covered cells
       do j = 1, ny
          do i = 1, nx
             if (ice_mask(i,j) == 1) then
                f_flotation(i,j) = -rhoo*(topg(i,j) - eus) / (rhoi*thck(i,j))
             endif
          enddo
       enddo

    elseif (whichflotation_function == HO_FLOTATION_FUNCTION_INVERSE_PATTYN) then ! grounded if f_flotation >= 1, else floating

       ! initialize to a large positive value (strongly grounded) on land, and 1.0 (almost floating) in the ocean
       do j = 1, ny
          do i = 1, nx
             if (land_mask(i,j) == 1) then
                f_flotation(i,j) = grounding_factor_max
             else    ! marine cell
                f_flotation(i,j) = 1.0d0
             endif
          enddo
       enddo

       ! set f_flotation in ice-covered cells
       do j = 1, ny
          do i = 1, nx
             if (ice_mask(i,j) == 1) then  ! ice is present
                if (land_mask(i,j) == 1) then
                   ! nothing else to do; f_flotation(i,j) set to grounding_factor_max above
                else    ! marine cell
                   f_flotation(i,j) = rhoi*thck(i,j) / (-rhoo*(topg(i,j) - eus))
                   f_flotation(i,j) = min(f_flotation(i,j), grounding_factor_max)
                endif
             endif  ! ice_mask
          enddo
       enddo

    elseif (whichflotation_function == HO_FLOTATION_FUNCTION_LINEAR) then ! grounded if f_flotation <= 0, else floating

       ! If > 0, f_flotation is the thickness of the ocean cavity beneath the ice shelf.
       ! This function (unlike PATTYN and INVERSE_PATTYN) is linear in both thck and topg.

       ! initialize to a large negative value (strongly grounded) on land, and 0.0 (almost floating) in the ocean
       do j = 1, ny
          do i = 1, nx
             if (land_mask(i,j) == 1) then
                f_flotation(i,j) = -grounding_factor_max
             else    ! marine cell
                f_flotation(i,j) = 0.0d0
             endif
          enddo
       enddo

       ! set f_flotation in ice-covered cells
       do j = 1, ny
          do i = 1, nx
             if (ice_mask(i,j) == 1) then
                f_flotation(i,j) = -rhoo*(topg(i,j) - eus) - rhoi*thck(i,j)
             endif
          enddo
       enddo

    endif  ! whichflotation_function

    ! initialize f_ground
    f_ground(:,:) = 0.0d0

    ! Compute f_ground according to the value of whichground

    select case(whichground)

    case(HO_GROUND_NO_GLP)   ! default: no grounding-line parameterization
                             ! f_ground = 1 at vertex if any neighbor cell is land or has grounded ice

       ! compute a mask that is true for cells that are land and/or have grounded ice
       do j = 1, ny
          do i = 1, nx
             if ((ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) .or. land_mask(i,j) == 1) then
                cground(i,j) = .true.
             else
                cground(i,j) = .false.
             endif
          enddo
       enddo

       ! vertices are grounded if any neighbor cell is land and/or has grounded ice, else are floating

       do j = 1, ny-1
          do i = 1, nx-1
             if (vmask(i,j) == 1) then
                if (cground(i,j+1) .or. cground(i+1,j+1) .or. cground(i,j) .or. cground(i+1,j)) then
                   f_ground(i,j) = 1.d0
                else
                   f_ground(i,j) = 0.d0
                endif
             endif
           enddo
        enddo

    case(HO_GROUND_ALL)

       ! all vertices with ice-covered neighbors are assumed grounded, regardless of thck and topg

       do j = 1, ny-1
          do i = 1, nx-1
             if (vmask(i,j) == 1) then
                f_ground(i,j) = 1.d0
             endif
          enddo
       enddo

    case(HO_GROUND_GLP)      ! grounding-line parameterization

       ! In ice-free cells, fill in f_flotation by extrapolation.
       ! Take the minimum (i.e., most grounded) value from adjacent ice-filled neighbors, using
       !  edge neighbors (if possible) or corner neighbors (if there are no ice-filled edge neighbors).
       ! The reason for this fairly intricate calculation is to make sure that each vertex with vmask = 1
       !  (i.e., with at least one ice-filled neighbor cell) has physically sensible values
       !  of f_flotation in all four neighbor cells, for purposes of interpolation.

       ! In order for the 'min' logic to work, temporarily reverse the sign for the inverse Pattyn function
       if (whichflotation_function == HO_FLOTATION_FUNCTION_INVERSE_PATTYN) then
          f_flotation(:,:) = -f_flotation(:,:)
       endif

       f_flotation_extrap(:,:) = f_flotation(:,:)

       do j = 2, ny-1
          do i = 2, nx-1
             if (ice_mask(i,j) == 0) then

                filled = .false.

                ! loop over edge neighbors
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if ((ii==i .or. jj==j) .and. ice_mask(ii,jj) == 1) then  ! edge neighbor with ice
                         if (.not.filled) then
                            filled = .true.
                            f_flotation_extrap(i,j) = f_flotation(ii,jj)
                         else
                            f_flotation_extrap(i,j) = min(f_flotation_extrap(i,j), f_flotation(ii,jj))
                         endif
                      endif
                   enddo   ! ii
                enddo   ! jj

                ! loop over corner neighbors if necessary
                if (.not.filled) then
                   do jj = j-1, j+1
                      do ii = i-1, i+1
                         if ((abs(ii-i)==1 .and. abs(jj-j)==1) .and. ice_mask(ii,jj)==1) then ! corner neighbor with ice
                            if (.not.filled) then
                               filled = .true.
                               f_flotation_extrap(i,j) = f_flotation(ii,jj)
                            else
                               f_flotation_extrap(i,j) = min(f_flotation_extrap(i,j), f_flotation(ii,jj))
                            endif
                         endif
                      enddo
                   enddo
                endif   ! not filled

             endif   ! ice_mask = 0
          enddo   ! i
       enddo  !j

       ! Reverse the temporary sign change for inverse Pattyn
       if (whichflotation_function == HO_FLOTATION_FUNCTION_INVERSE_PATTYN) then
          f_flotation_extrap(:,:) = -f_flotation_extrap(:,:)
       endif

       ! halo update
       call parallel_halo(f_flotation_extrap)

       ! copy the extrapolated array to the main f_flotation array
       f_flotation(:,:) = f_flotation_extrap(:,:)

       ! Identify cells that contain floating ice.
       ! Note: Ice-free cells adjacent to floating cells (and not adjacent to grounded cells)
       !        will have cfloat = .true. from the extrapolation above.

       if (whichflotation_function == HO_FLOTATION_FUNCTION_PATTYN) then

          ! grounded if f_flotation <= 1, else floating

          do j = 1, ny
             do i = 1, nx
                if (f_flotation(i,j) > 1.d0) then
                   cfloat(i,j) = .true.
                else
                   cfloat(i,j) = .false.
                endif
             enddo
          enddo
          
       elseif (whichflotation_function == HO_FLOTATION_FUNCTION_INVERSE_PATTYN) then
          
          ! grounded if f_flotation >= 1, else floating

          do j = 1, ny
             do i = 1, nx
                if (f_flotation(i,j) < 1.d0) then
                   cfloat(i,j) = .true.
                else
                   cfloat(i,j) = .false.
                endif
             enddo
          enddo
          
       elseif (whichflotation_function == HO_FLOTATION_FUNCTION_LINEAR) then

          ! grounded if f_flotation <= 0, else floating

          do j = 1, ny
             do i = 1, nx
                if (f_flotation(i,j) > 0.d0) then
                   cfloat(i,j) = .true.
                else
                   cfloat(i,j) = .false.
                endif
             enddo
          enddo
          
       endif
       
       !WHL - debug
       if (this_rank == rtest) then
          i = it; j = jt
          print*, 'i, j =', i, j
          print*, 'f_flotation(i:i+1,j+1):', f_flotation(i:i+1,j+1)
          print*, 'f_flotation(i:i+1,j)  :', f_flotation(i:i+1,j)
          print*, 'cfloat(i:i+1,j+1):', cfloat(i:i+1,j+1)
          print*, 'cfloat(i:i+1,j)  :', cfloat(i:i+1,j)
       endif

       ! Loop over vertices, computing f_ground for each vertex with vmask = 1.
       ! Note: All vertices with vmask = 1 (i.e., all vertices with at least one ice-covered neighbor)
       !        are surrounded by four cells with physically meaningful values of f_flotation and cfloat.
       !       That is, f_flotation in each surrounding cell is derived from the cell itself (if the cell
       !        is ice-covered) or from an adjacent ice-covered cell (if the cell is ice-free).
       !       Thus, f_flotation can be interpolated from four physically sensible values.
       !TODO - Test this algorithm for the LINEAR case. Not sure it works when the boundary
       !       between grounded and floating ice is 0 rather than 1.

       do j = 1, ny-1
          do i = 1, nx-1

             if (vmask(i,j) == 1) then  ! ice is present in at least one neighboring cell

                   ! First count the number of floating cells surrounding this vertex

                   nfloat = 0
                   if (cfloat(i,j))     nfloat = nfloat + 1
                   if (cfloat(i+1,j))   nfloat = nfloat + 1
                   if (cfloat(i+1,j+1)) nfloat = nfloat + 1
                   if (cfloat(i,j+1))   nfloat = nfloat + 1

                   !WHL - debug
                   if (i==it .and. j==jt .and. this_rank == rtest) then
                      print*, ' '
                      print*, 'nfloat =', nfloat
                   endif

                   ! Given nfloat, compute f_ground for each vertex
                   ! First the easy cases...
                
                   if (nfloat == 0) then

                      f_ground(i,j) = 1.d0    ! fully grounded

                   elseif (nfloat == 4) then

                      f_ground(i,j) = 0.d0    ! fully floating

                   ! For the other cases the grounding line runs through the rectangular region 
                   !  around this vertex.
                   ! Using the values at the 4 neighboring cells, we approximate f_flotation(x,y) as
                   !  a bilinear function f(x,y) = a + bx + cy + dxy over the region.
                   ! To find f_ground, we integrate over the region with f(x,y) <= 1
                   !  (or alternatively, we find f_float = 1 - f_ground by integrating
                   !  over the region with f(x,y) > 1).
                   !  
                   ! There are 3 patterns to consider:
                   ! (1) nfloat = 1 or nfloat = 3 (one cell neighbor is not like the others)
                   ! (2) nfloat = 2, and adjacent cells are floating
                   ! (3) nfloat = 2, and diagonally opposite cells are floating

                   elseif (nfloat == 1 .or. nfloat == 3) then
 
                      if (nfloat==1) then
                         logvar(1:2,1:2) = cfloat(i:i+1,j:j+1)
                      else  ! nfloat = 3
                         logvar(1:2,1:2) = .not.cfloat(i:i+1,j:j+1)
                      endif
                      
                      ! Identify the cell that is not like the others
                      ! (i.e., the only floating cell if nfloat = 1, or the only
                      !  grounded cell if nfloat = 3)
                      !
                      ! Diagrams below are for the case nfloat = 1.
                      ! If nfloat = 3, the F and G labels are switched.

                      if (logvar(1,1)) then             ! no rotation
                         f1 = f_flotation(i,j)          !   G-----G
                         f2 = f_flotation(i+1,j)        !   |     |
                         f3 = f_flotation(i+1,j+1)      !   |     |
                         f4 = f_flotation(i,j+1)        !   F-----G

                      elseif (logvar(2,1)) then         ! rotate by 90 degrees
                         f4 = f_flotation(i,j)          !   G-----G
                         f1 = f_flotation(i+1,j)        !   |     |
                         f2 = f_flotation(i+1,j+1)      !   |     |
                         f3 = f_flotation(i,j+1)        !   G-----F

                      elseif (logvar(2,2)) then         ! rotate by 180 degrees
                         f3 = f_flotation(i,j)          !   G-----F
                         f4 = f_flotation(i+1,j)        !   |     |
                         f1 = f_flotation(i+1,j+1)      !   |     |
                         f2 = f_flotation(i,j+1)        !   G-----G

                      elseif (logvar(1,2)) then         ! rotate by 270 degrees
                         f2 = f_flotation(i,j)          !   F-----G
                         f3 = f_flotation(i+1,j)        !   |     |
                         f4 = f_flotation(i+1,j+1)      !   |     |
                         f1 = f_flotation(i,j+1)        !   G-----G
                      endif
                      
                      ! Compute coefficients in f(x,y) = a + b*x + c*y + d*x*y
                      ! Note: x is to the right and y is up if the southwest cell is not like the others.
                      !       For the other cases we solve the same problem with x and y rotated.
                      !       The rotations are handled by rotating f1, f2, f3 and f4 above.

                      a = f1
                      b = f2 - f1
                      c = f4 - f1
                      d = f1 + f3 - f2 - f4

                      !WHL - debug
                      if (i==it .and. j==jt .and. this_rank == rtest) then
                         print*, 'f1, f2, f3, f4 =', f1, f2, f3, f4
                         print*, 'a, b, c, d =', a, b, c, d
                      endif

                      ! Compute the fractional area of the corner region 
                      ! (floating if nfloat = 1, grounded if nfloat = 3)
                      !
                      ! Here are the relevant integrals:
                      !
                      ! (1) d /= 0:
                      !     integral_0^x0 {y(x) dx}, where x0   = (1-a)/b
                      !                                    y(x) = (1 - (a+b*x)) / (c+d*x)
                      !     = [bc - ad + d) ln(1 + d(1-a)/(bc)) - (1-a)d] / d^2
                      !
                      ! (2) d = 0:
                      !     integral_0^x0 {y(x) dx}, where x0   = (1-a)/b
                      !                                    y(x) = (1 - (a+b*x)) / c
                      !     = (a-1)(a-1) / (2bc)
                      !
                      ! Note: We cannot have bc = 0, because f_flotation varies in both x and y

                      if (abs(d*(1.d0-a)/(b*c)) > eps06) then
                         f_corner = ((b*c - a*d + d) * log(1.d0 + d*(1.d0 - a)/(b*c)) - (1.d0 - a)*d) / (d*d)
                      else
                         f_corner = (a - 1.d0)*(a - 1.d0) / (2.d0*b*c)
                      endif

                      if (nfloat==1) then  ! f_corner is the floating area
                         f_ground(i,j) = 1.d0 - f_corner
                      else                 ! f_corner is the grounded area
                         f_ground(i,j) = f_corner
                      endif

                      !WHL - debug
                      if (i==it .and. j==jt .and. this_rank == rtest) then
                         print*, 'f_corner =', f_corner
                         print*, 'f_ground =', f_ground(i,j)
                      endif

                   elseif (nfloat == 2) then

                      ! first the 4 cases where the 2 grounded cells are adjacent
                      ! We integrate over the trapezoid in the floating part of the cell

                      if (cfloat(i,j) .and. cfloat(i+1,j)) then  ! no rotation
                         adjacent = .true.              !   G-----G
                         f1 = f_flotation(i,j)          !   |     |
                         f2 = f_flotation(i+1,j)        !   |     |
                         f3 = f_flotation(i+1,j+1)      !   |     |
                         f4 = f_flotation(i,j+1)        !   F-----F

                      elseif (cfloat(i+1,j) .and. cfloat(i+1,j+1)) then  ! rotate by 90 degrees
                         adjacent = .true.              !   G-----F
                         f4 = f_flotation(i,j)          !   |     |
                         f1 = f_flotation(i+1,j)        !   |     |
                         f2 = f_flotation(i+1,j+1)      !   |     |
                         f3 = f_flotation(i,j+1)        !   G-----F

                      elseif (cfloat(i+1,j+1) .and. cfloat(i,j+1)) then  ! rotate by 180 degrees
                         adjacent = .true.              !   F-----F
                         f3 = f_flotation(i,j)          !   |     |
                         f4 = f_flotation(i+1,j)        !   |     |
                         f1 = f_flotation(i+1,j+1)      !   |     |
                         f2 = f_flotation(i,j+1)        !   G-----G

                      elseif (cfloat(i,j+1) .and. cfloat(i,j)) then   ! rotate by 270 degrees
                         adjacent = .true.              !   F-----G
                         f2 = f_flotation(i,j)          !   |     |
                         f3 = f_flotation(i+1,j)        !   |     |
                         f4 = f_flotation(i+1,j+1)      !   |     |
                         f1 = f_flotation(i,j+1)        !   F-----G

                      else   ! the 2 grounded cells are diagonally opposite

                         adjacent = .false.

                         ! We will integrate assuming the two corner regions lie in the lower left
                         ! and upper right, i.e. one of these patterns:
                         !
                         !   F-----G       G-----F
                         !   |     |       |     |
                         !   |  F  |       |  G  |
                         !   |     |       |     |
                         !   G-----F       F-----G
                         !
                         ! Two other patterns are possible, with corner regions in the lower right
                         ! and upper left; these require a rotation before integrating: 
                         
                         !   G-----F       F-----G
                         !   |     |       |     |
                         !   |  F  |       |  G  |
                         !   |     |       |     |
                         !   F-----G       G-----F
                         !   
                         var = f_flotation(i+1,j)*f_flotation(i,j+1) - f_flotation(i,j)*f_flotation(i+1,j+1)   &
                             + f_flotation(i,j) + f_flotation(i+1,j+1) - f_flotation(i+1,j) - f_flotation(i,j+1)

                         if (var >= 0.d0) then   ! we have one of the top two patterns
                            f1 = f_flotation(i,j)
                            f2 = f_flotation(i+1,j)
                            f3 = f_flotation(i+1,j+1)
                            f4 = f_flotation(i,j+1)
                         else   ! we have one of the bottom two patterns; rotate coordinates by 90 degrees
                            f4 = f_flotation(i,j)
                            f1 = f_flotation(i+1,j)
                            f2 = f_flotation(i+1,j+1)
                            f3 = f_flotation(i,j+1)
                         endif

                      endif  ! grounded cells are adjacent

                      ! Compute coefficients in f(x,y) = a + b*x + c*y + d*x*y
                      a = f1
                      b = f2 - f1
                      c = f4 - f1
                      d = f1 + f3 - f2 - f4

                      ! Integrate the corner areas

                      !WHL - debug
                      if (i==it .and. j==jt .and. this_rank == rtest) then
                         print*, 'adjacent =', adjacent
                         print*, 'f1, f2, f3, f4 =', f1, f2, f3, f4
                         print*, 'a, b, c, d =', a, b, c, d
                      endif

                      if (adjacent) then

                         ! Compute the area of the floating part of the cell
                         ! Here are the relevant integrals:
                         !
                         ! (1) d /= 0:
                         !     integral_0^1 {y(x) dx}, where y(x) = (1 - (a+b*x)) / (c+d*x)
                         !                                  
                         !     = [bc - ad + d) ln(1 + d/c) - bd] / d^2
                         !
                         ! (2) d = 0:
                         !     integral_0^1 {y(x) dx}, where y(x) = (1 - (a+b*x)) / c
                         !                                 
                         !     = -(2a + b - 2) / (2c)
                         !
                         ! Note: We cannot have c = 0, because the passage of the GL
                         !       through the region from left to right implies variation in y.
                         
                         if (abs(d/c) > eps06) then   ! the usual case
                            f_trapezoid = ((b*c - a*d + d) * log(1.d0 + d/c) - b*d) / (d*d)
                         else
                            f_trapezoid = -(2.d0*a + b - 2.d0) / (2.d0*c)
                         endif

                         f_ground(i,j) = 1.d0 - f_trapezoid

                         !WHL - debug
                         if (i==it .and. j==jt .and. this_rank == rtest) then
                            print*, 'f_trapezoid =', f_trapezoid
                            print*, 'f_ground =', f_ground(i,j)
                         endif

                      else   ! grounded vertices are diagonally opposite

                         ! bug check: make sure some signs are positive as required by the formulas
                         if (b*c - d*(1.d0-a) < 0.d0) then
                            print*, 'Grounding line error: bc - d(1-a) < 0'
                            stop
                         elseif (b*c < 0.d0) then
                            print*, 'Grounding line error: bc < 0'
                            stop
                         elseif ((b+d)*(c+d) < 0.d0) then
                            print*, 'Grounding line error: (b+d)(c+d) < 0'
                            stop
                         endif

                         ! Compute the combined areas of the two corner regions.
                         ! For the lower left region, the integral is the same as above
                         ! (for the case nfloat = 1 or nfloat = 3, with d /= 0).
                         ! For the upper right region, here is the integral:
                         !
                         !     integral_x1^1 {(1-y(x)) dx}, where x1  = (1-a-c)/(b+d)
                         !                                       y(x) = (1 - (a+b*x)) / (c+d*x)
                         !     = {(bc - ad + d) ln[(bc + d(1-a))/((b+d)(c+d))] + d(a + b + c + d - 1)} / d^2
                         !
                         ! The above integral is valid only if (bc + d(1-a)) > 0.
                         ! If this quantity = 0, then the grounding line lies along two lines,
                         ! x0 = (1-a)/b and y0 = (1-a)/c.
                         ! The lower left area is x0*y0 = (1-a)^2 / (bc).
                         ! The upper right area is (1-x0)*(1-y0) = (a+b-1)(a+c-1) / (bc)
                         !
                         ! Note that this pattern is not possible with d = 0

                         !WHL - debug
                         if (i==it .and. j==jt .and. this_rank == rtest) then
                            print*, 'Pattern 3: i, j, bc + d(1-a) =', i, j, b*c + d*(1.d0-a)
                         endif

                         if (abs(d*(1.d0-a)/(b*c)) > eps06) then  ! the usual case
                            f_corner1 = ((b*c - a*d + d) * log(1.d0 + d*(1.d0-a)/(b*c)) - (1.d0-a)*d) / (d*d)
                            f_corner2 = ((b*c - a*d + d) * log((b*c + d*(1.d0-a))/((b+d)*(c+d)))  &
                                      + d*(a + b + c + d - 1)) / (d*d)
                         else 
                            f_corner1 = (1.d0 - a)*(1.d0 - a) / (b*c)
                            f_corner2 = (a + b - 1.d0)*(a + c - 1.d0) / (b*c)
                         endif
                         
                         ! Determine whether the central point (1/2,1/2) is grounded or floating.
                         ! Then compute the grounded area.
                         ! If the central point is floating, the corner regions are grounded;
                         ! if the central point is grounded, the corner regions are floating.

                         f_flotation_vertex = a + 0.5d0*b + 0.5d0*c + 0.25d0*d
                         if (f_flotation_vertex > 1.d0) then  ! the central point is floating; corners are grounded
                            f_ground(i,j) = f_corner1 + f_corner2
                         else                              ! the central point is grounded; corners are floating
                            f_ground(i,j) = 1.d0 - (f_corner1 + f_corner2)
                         endif

                         !WHL - debug
                         if (i==it .and. j==jt .and. this_rank == rtest) then
                            print*, 'f_flotation_v =', f_flotation_vertex
                            print*, 'f_corner1 =', f_corner1
                            print*, 'f_corner2 =', f_corner2
                            print*, 'f_ground =', f_ground(i,j)
                         endif

                      endif  ! adjacent or opposite

                   endif     ! nfloat

             endif        ! vmask = 1
          enddo           ! i
       enddo              ! j

    end select

  end subroutine glissade_grounded_fraction

!=======================================================================

  subroutine glissade_grounding_line_flux(nx,                       ny,            &
                                          dx,                       dy,            &
                                          sigma,                                   &
                                          thck,                                    &
                                          uvel,                     vvel,          &
                                          ice_mask,                 floating_mask, &
                                          ocean_mask,                              &
                                          gl_flux_east,             gl_flux_north, &
                                          gl_flux                                   )

    ! Computes northward and eastward land ice fluxes at grounding lines,
    !  and a cell-based grounding-line flux field.
    ! Note: Since the GL thicknesses are approximated, the GL fluxes will not exactly 
    !        match the fluxes computed by the transport scheme.
    !       Also, the GL fluxes do not include thinning/calving of grounded marine cliffs.

    use parallel, only: nhalo
    use glimmer_paramets, only: thk0, vel0, len0

    implicit none

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::                     &
        nx, ny                                   !> horizontal grid dimensions

    real(dp), intent(in) ::                    &
        dx, dy                                   !> horizontal grid spacing

    real(dp), dimension(:), intent(in) ::      &
        sigma                                    !> vertical sigma coordinate

    real(dp), dimension(nx,ny), intent(in) ::  &
        thck                                     !> ice thickness

    real(dp), dimension(:,:,:), intent(in) ::  &
        uvel, vvel                               !> ice velocity in x and y directions

    integer, dimension(nx,ny), intent(in) ::  &
        ice_mask,                              & !> = 1 where ice is present, else = 0
        floating_mask,                         & !> = 1 where ice is present and floating, else = 0
        ocean_mask                               !> = 1 for ice-free ocean, else = 0

    ! Note: gl_flux_east and gl_flux_north are directional 
    !       (positive for eastward/northward, negative for westward/southward)
    !       gl_flux is a cell-based quantity based on flux magnitudes on each edge
    !       (so gl_flux >= 0)

    real(dp), dimension(:,:), intent(out) ::   &
        gl_flux_east,                          & !> grounding line flux on east edges
        gl_flux_north,                         & !> grounding line flux on north edges
        gl_flux                                  !> grounding line flux per grid cell


    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer  :: i,j,k                                     !> local cell indices
    integer  :: upn                                       !> vertical grid dimension
    real(dp), dimension(:), allocatable :: uavg, vavg     !> local horizontal velocity averages
    real(dp) :: thck_gl                                   !> GL thickness derived from topg_gl

    upn = size(sigma)

    allocate(uavg(upn), vavg(upn))

    ! Initialize
    gl_flux_east(:,:)  = 0.d0
    gl_flux_north(:,:) = 0.d0
    gl_flux(:,:)       = 0.d0

    ! Compute grounding line fluxes on east and north edges.
    ! Look for edges with a grounded cell on one side and a floating cell on the other.

    do j = nhalo+1, ny-nhalo
        do i = nhalo+1, nx-nhalo

            ! check east edge
           if ( (   (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) .and.   &  ! (i,j) grounded
                (ocean_mask(i+1,j) == 1 .or.  floating_mask(i+1,j) == 1) )     &  ! (i+1,j) floating or ocean
                                        .or.                                   &
                ( (ice_mask(i+1,j) == 1 .and. floating_mask(i+1,j) == 0) .and. &  ! (i+1,j) grounded
                  (ocean_mask(i,j) == 1  .or. floating_mask(i,j) == 1) ) ) then   ! (i,j) floating or ocean

                uavg(:) = (uvel(:,i,j) + uvel(:,i,j-1)) / 2.d0
                if (ice_mask(i,j) == 1 .and. ice_mask(i+1,j) == 1) then
                   ! set GL thickness to the average thickness of the two cells
                   thck_gl = (thck(i,j) + thck(i+1,j)) / 2.d0
                else
                   ! set GL thickness to the thickness of the ice-filled cell
                   thck_gl = max(thck(i,j), thck(i+1,j))
                endif

                do k = 1, upn-1
                    gl_flux_east(i,j) = gl_flux_east(i,j) &
                                        + thck_gl * (sigma(k+1) - sigma(k)) * (uavg(k) + uavg(k+1))/2.d0
                enddo
            endif

            ! check north edge
           if ( (   (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) .and.   &  ! (i,j) grounded
                (ocean_mask(i,j+1) == 1 .or.  floating_mask(i,j+1) == 1) )     &  ! (i,j+1) floating or ocean
                                        .or.                                   &
                ( (ice_mask(i,j+1) == 1 .and. floating_mask(i,j+1) == 0) .and. &  ! (i,j+1) grounded
                  (ocean_mask(i,j) == 1  .or. floating_mask(i,j) == 1) ) ) then   ! (i,j) floating or ocean

                vavg(:) = (vvel(:,i-1,j) + vvel(:,i,j)) / 2.d0
                if (ice_mask(i,j) == 1 .and. ice_mask(i,j+1) == 1) then
                   ! set GL thickness to the average thickness of the two cells
                   thck_gl = (thck(i,j) + thck(i,j+1)) / 2.d0
                else
                   ! set GL thickness to the thickness of the ice-filled cell
                   thck_gl = max(thck(i,j), thck(i,j+1))
                endif

                do k = 1, upn-1
                    gl_flux_north(i,j) = gl_flux_north(i,j) &
                                        + thck_gl * (sigma(k+1) - sigma(k)) * (vavg(k) + vavg(k+1))/2.d0
                enddo
             endif

        enddo   ! i
    enddo   ! j

    ! Compute mass flux through grounding line in each cell.
    ! Only a grounded cell can lose mass. We need to check the direction of the fluxes.

    do j = nhalo+1,ny-nhalo
        do i = nhalo+1,nx-nhalo

            ! Check the sign for east-west flow and assign the flux accordingly
            if (gl_flux_east(i,j) < 0.d0) then
                ! The ice is flowing westward and the flux belongs to the right adjacent cell
                gl_flux(i+1,j) = gl_flux(i+1,j) - gl_flux_east(i,j)
            else
                ! The ice is flowing eastward and the flux belongs to this cell
                gl_flux(i,j) = gl_flux(i,j) + gl_flux_east(i,j)
            endif

            ! Check the sign for north-south flow and assign the flux accordingly
            if (gl_flux_north(i,j) < 0.d0) then
                ! The ice is flowing southward and the flux belongs to the top adjacent cell
                gl_flux(i,j+1) = gl_flux(i,j+1) - gl_flux_north(i,j)
            else
                ! The ice is flowing northward and the flux belongs to this cell
                gl_flux(i,j) = gl_flux(i,j) + gl_flux_north(i,j)
            endif

        enddo   ! i
    enddo   ! j

    ! Convert from model units to kg/m/s
    gl_flux_east  = gl_flux_east  * rhoi*thk0*vel0
    gl_flux_north = gl_flux_north * rhoi*thk0*vel0
    gl_flux       = gl_flux       * rhoi*thk0*vel0

    deallocate(uavg, vavg)

  end subroutine glissade_grounding_line_flux

!****************************************************************************

end module glissade_grounding_line

!****************************************************************************

