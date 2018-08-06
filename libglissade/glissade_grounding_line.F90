!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grounding_line.F90 - part of the Community Ice Sheet Model (CISM)  
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

    private
    public :: glissade_grounded_fraction, glissade_grounding_line_flux

    logical, parameter :: verbose_gl = .false.

  contains

!****************************************************************************

  subroutine glissade_grounded_fraction(nx,              ny,                       &
                                        itest, jtest,    rtest,                       &
                                        thck,            topg,                        &
                                        eus,             ice_mask,                    &
                                        floating_mask,   land_mask,                   &
                                        which_ho_ground, which_ho_flotation_function, &
                                        f_flotation,                                  &
                                        f_ground,        f_ground_cell)

    use glissade_grid_operators, only : glissade_unstagger

    !----------------------------------------------------------------
    ! Compute fraction of ice that is grounded.
    ! This fraction is computed at vertices based on the thickness and topography of the four neighboring cell centers.
    ! More generally, it can be computed at any point centered in a square or rectangle,
    !  given a flotation function (based on thickness and topography) at the corners of the square or rectangle.
    !
    ! There are three options for computing the grounded fraction, based on the value of which_ho_ground:
    ! (0) HO_GROUND_NO_GLP: f_ground = 1 for vertices with grounded and/or land-based neighbor cells
    !                       f_ground = 0 for vertices with floating neighbors only
    ! (1) HO_GROUND_GLP: 0 <= f_ground <= 1 based on grounding-line parameterization
    !        A flotation function is interpolated over the staggered cell around each vertex
    !        and analytically integrated to compute the grounded and floating fractions.
    !        The f_ground_cell is compute by averaging the values computed at vertices.
    ! (2) HO_GROUND_GLP_QUADRANTS: similar to (2), but f_flotation is interpolated over cell quadrants,
    !        and then f_ground is computed in unstaggered and staggered cells by summing over quadrants.
    !
    ! For options (1) and (2), there are three options for the flotation function:
    ! (0) HO_FLOTATION_FUNCTION_PATTYN: f_flotation = (-rhoo*b)/(rhoi*H) - 1 = f_pattyn - 1
    !     Here, f_pattyn = (-rhoo*b)/(rhoi*H) as in Pattyn et al. (2006).
    ! (1) HO_FLOTATION_FUNCTION_INVERSE_PATTYN: f_flotation = 1 - (rhoi*H)/(-rhoo*b) = 1 - 1/f_pattyn
    ! (2) HO_FLOTATION_FUNCTION_LINEAR: f_flotation = -b - (rhoi/rhoo)*H = ocean cavity thickness
    !     This function was suggested by Xylar Asay-Davis and is linear in both b and H.
    ! All three functions are defined such that f <=0 for grounded ice and f > 0 for floating ice.
    ! For each option, land-based cells are assigned a large negative value, so that any vertices
    !  with land-based neighbors are strongly grounded.
    !
    ! We first compute f_flotation in all active ice-covered cells.
    ! Then f_flotation is extrapolated to ice-free neighbors.  Thus, f_flotation has a physically
    !   meaningful value (either computed directly, or extrapolated from a neighbor) in all four
    !   cells surrounding each active vertex. (By definition, an active vertex is a vertex with
    !   at least one active ice-covered neighbor.) Thus, we can interpolate f_flotation
    !   within the staggered cell around each active vertex to compute f_ground at the vertex.
    !
    ! The linear function (2) is the default.
    !
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny                ! number of grid cells in each direction

    integer, intent(in) ::   &
       itest, jtest, rtest    ! coordinates of diagnostic point

    ! Default dimensions are meters.
    ! The subroutine will work for other units as long as thck, topg and eus have the same units,
    !  but f_flotation_land_linear will be scaled incorrectly.

    real(dp), dimension(nx,ny), intent(in) ::  &
       thck,                & ! ice thickness (m)
       topg                   ! elevation of topography (m)

    real(dp), intent(in) :: &
       eus                    ! eustatic sea level (= 0 by default)

    integer, dimension(nx,ny), intent(in) ::   &
       ice_mask,            & ! = 1 if ice is present (thck > thklim), else = 0
       floating_mask,       & ! = 1 if ice is present (thck > thklim) and floating, else = 0
       land_mask              ! = 1 if topg is at or above sea level

    ! see comments above for more information about these options
    integer, intent(in) ::  &
       which_ho_ground,            &! option for computing f_ground
       which_ho_flotation_function  ! option for computing f_flotation

    real(dp), dimension(nx,ny), intent(out) :: &
       f_flotation            ! flotation function; see comments above

    real(dp), dimension(nx-1,ny-1), intent(out) ::  &
       f_ground               ! grounded ice fraction at vertex, 0 <= f_ground <= 1

    real(dp), dimension(nx,ny), intent(out) ::  &
       f_ground_cell          ! grounded ice fraction in cell, 0 <= f_ground <= 1

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, ii, jj

    integer, dimension(nx-1,ny-1) ::   &
         vmask                     ! = 1 for vertices neighboring at least one cell where ice is present, else = 0

    real(dp), dimension(4,nx,ny) :: &
         f_ground_quadrant         ! f_ground for the 4 cell quadrants around each vertex

    real(dp), dimension(nx,ny) :: &
         f_flotation_extrap        ! f_flotation, extrapolated to cells without active ice

    real(dp), dimension(4) :: &
         f_flotation_vector        ! 1D array containing f_flotation for 4 corners of a quadrant

    real(dp) ::  &
         topg_eus_diff      ! topg - eus, limited to be >= f_flotation_land_topg_min

    logical, dimension(nx,ny) :: &
         cground            ! true if a cell is land and/or has grounded ice, else = false

    real(dp), parameter :: &
         f_flotation_land_topg_min = 1.0d0   ! min value of (topg - eus) in f_flotation expression for land cells (m)

    logical :: filled       ! true if f_flotation has been filled by extrapolation

    !TODO - Test sensitivity to these values
    ! These are set to large negative values, so vertices with land-based neighbors are strongly grounded.
    real(dp), parameter :: f_flotation_land_pattyn = -10.d0          ! unitless

    !----------------------------------------------------------------
    ! Compute ice mask at vertices (= 1 if any surrounding cells have ice or are land)
    !----------------------------------------------------------------

    do j = 1, ny-1
       do i = 1, nx-1
          if (ice_mask(i,j+1)==1  .or. ice_mask(i+1,j+1)==1  .or.   &
              ice_mask(i,j)  ==1  .or. ice_mask(i+1,j)  ==1  .or.   &
              land_mask(i,j+1)==1 .or. land_mask(i+1,j+1)==1 .or.   &
              land_mask(i,j)  ==1 .or. land_mask(i+1,j+1)==1) then
             vmask(i,j) = 1
          else
             vmask(i,j) = 0
          endif
       enddo
    enddo

    ! Compute flotation function at cell centers.
    ! For diagnostic purposes, f_flotation is always computed, although it affects f_ground
    !  (and thus the velocities) only when running with a GLP.
    ! Note: f_flotation is set to an arbitrary large negative value for land-based cells.
    !       Ice-free ocean cells have f_flotation = 0.
    ! Note: Values from ice-free cells are not used in the calculation of f_ground.
    !       All values used in the f_ground calculation come from cells with ice,
    !        or are extrapolated from cells with ice.

    if (which_ho_flotation_function == HO_FLOTATION_FUNCTION_PATTYN) then

       ! subtract 1 from (-rhoo*b)/(rhoi*H) so that f > 0 for floating ice, f <= 0 for grounded ice
       do j = 1, ny
          do i = 1, nx
             if (land_mask(i,j) == 1) then  ! topg - eus >= 0
                f_flotation(i,j) = f_flotation_land_pattyn
             elseif (ice_mask(i,j) == 1) then
                f_flotation(i,j) = -rhoo*(topg(i,j) - eus) / (rhoi*thck(i,j)) - 1.0d0
             else  ! ice-free ocean
                f_flotation(i,j) = 0.0d0
             endif
          enddo
       enddo

    elseif (which_ho_flotation_function == HO_FLOTATION_FUNCTION_INVERSE_PATTYN) then

       ! subtract (rhoi*H)/(-rhoo*b) from 1 so that f > 0 for floating ice, f <= 0 for grounded ice
       do j = 1, ny
          do i = 1, nx
             if (land_mask(i,j) == 1) then  ! topg - eus >= 0
                f_flotation(i,j) = f_flotation_land_pattyn
             elseif (ice_mask(i,j) == 1) then
                f_flotation(i,j) = 1.0d0 - rhoi*thck(i,j) / (-rhoo*(topg(i,j) - eus))
                ! Cap at a large minimum value
                f_flotation(i,j) = max(f_flotation(i,j), f_flotation_land_pattyn)
             else  ! ice-free ocean
                f_flotation(i,j) = 0.0d0
             endif
          enddo
       enddo

    elseif (which_ho_flotation_function == HO_FLOTATION_FUNCTION_LINEAR) then

       ! If > 0, f_flotation is the thickness of the ocean cavity beneath the ice shelf.
       ! This function (unlike PATTYN and INVERSE_PATTYN) is linear in both thck and topg.

       do j = 1, ny
          do i = 1, nx
             if (land_mask(i,j) == 1) then
                ! Assign a minimum value to (topg - eus) so that f_flotation is nonzero on land
                topg_eus_diff = max((topg(i,j) - eus), f_flotation_land_topg_min)
                f_flotation(i,j) = -topg_eus_diff
             elseif (ice_mask(i,j) == 1) then
                f_flotation(i,j) = -(topg(i,j) - eus) - (rhoi/rhoo)*thck(i,j)
             else  ! ice-free ocean
                f_flotation(i,j) = 0.0d0
             endif
          enddo
       enddo

    endif  ! which_ho_flotation_function

    ! initialize f_ground arrays
    f_ground(:,:) = 0.0d0
    f_ground_cell(:,:) = 0.0d0

    ! Compute f_ground according to the value of which_ho_ground

    if (which_ho_ground == HO_GROUND_NO_GLP) then

       ! default: no sub-grid grounding-line parameterization
       ! f_ground = 1 at a vertex if any neighbor cell is land or has grounded ice

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

       ! f_ground_cell follows directly from this mask

       where (cground)
          f_ground_cell = 1.0d0
       elsewhere
          f_ground_cell = 0.0d0
       endwhere

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

    else   ! grounding-line parameterization (HO_GROUND_GLP, HO_GROUND_GLP_QUADRANTS)

       ! In ice-free ocean cells, fill in f_flotation by extrapolation.
       ! Take the minimum (i.e., most grounded) value from adjacent ice-filled neighbors, using
       !  edge neighbors (if possible) or corner neighbors (if there are no ice-filled edge neighbors).
       ! The reason for this fairly intricate calculation is to make sure that each vertex with vmask = 1
       !  (i.e., with at least one ice-filled or land-based neighbor cell) has physically sensible values
       !  of f_flotation in all four neighbor cells, for purposes of interpolation.

       f_flotation_extrap(:,:) = f_flotation(:,:)

       do j = 2, ny-1
          do i = 2, nx-1
             if (ice_mask(i,j) == 0 .and. land_mask(i,j) == 0) then

                filled = .false.

                ! loop over edge neighbors
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if ((ii == i .or. jj == j) .and. &
                           (ice_mask(ii,jj) == 1 .or. land_mask(ii,jj) == 1)) then  ! edge neighbor with ice or land
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
                         if ((abs(ii-i) == 1 .and. abs(jj-j) == 1) .and. &
                             (ice_mask(ii,jj) == 1 .or. land_mask(ii,jj) == 1)) then ! corner neighbor with ice or land
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

       ! halo update
       call parallel_halo(f_flotation_extrap)

       ! copy the extrapolated array to the main f_flotation array
       f_flotation(:,:) = f_flotation_extrap(:,:)

       if (verbose_gl .and. this_rank == rtest) then
          i = itest; j = jtest
          print*, ' '
          print*, 'rank, i, j =', this_rank, i, j
          print*, 'f_flotation(i:i+1,j+1):', f_flotation(i:i+1,j+1)
          print*, 'f_flotation(i:i+1,j)  :', f_flotation(i:i+1,j)
       endif

    endif   ! which_ho_ground

    ! Given f_flotation in all cells, we have the information needed to compute f_ground
    !  by either of two methods:
    ! For which_ho_ground = 1, we interpolate f_flotation over the staggered cell around each vertex,
    !  then analytically evaluate the grounded fraction in the staggered cell.
    !  To obtain f_ground_cell, we average from vertices using glissade_unstagger.
    ! For which_ho_ground = 2, we interpolate f_flotation over the four quadrants in each cell,
    !  then analytically evaluate the grounded fraction in the quadrant, and finally
    !  sum over quadrants to obtain f_ground at vertices (staggered grid) and in cells (unstaggered grid).

    if (which_ho_ground == HO_GROUND_GLP) then

       ! Loop over vertices, computing f_ground for each vertex with vmask = 1.
       ! Note: All vertices with vmask = 1 (i.e., all vertices with at least one ice-covered neighbor)
       !        are surrounded by four cells with physically meaningful values of f_flotation.

       do j = 1, ny-1
          do i = 1, nx-1

             if (vmask(i,j) == 1) then  ! at least one neighboring cell is ice-covered or land-based

                f_flotation_vector(1) = f_flotation(i,j)
                f_flotation_vector(2) = f_flotation(i+1,j)
                f_flotation_vector(3) = f_flotation(i+1,j+1)
                f_flotation_vector(4) = f_flotation(i,j+1)

                call compute_grounded_fraction(i,     j,     1,     this_rank, &  ! '1' is the quadrant index
                                               itest, jtest, rtest,            &
                                               f_flotation_vector,             &
                                               f_ground(i,j))

             endif        ! vmask = 1
          enddo           ! i
       enddo              ! j

       ! Average f_ground from vertices to cells (includes a halo update)
       call glissade_unstagger(nx,       ny,  &
                               f_ground, f_ground_cell)

       ! Set f_ground_cell = 1 on land
       where (land_mask == 1)
          f_ground_cell = 1.0d0
       endwhere

    elseif (which_ho_ground == HO_GROUND_GLP_QUADRANTS) then

       !----------------------------------------------------------------
       ! Compute f_ground_vertex by finding f_ground for 4 neighboring quadrants
       !----------------------------------------------------------------

       ! Loop over vertices, computing f_ground for the four quadrants around each vertex with vmask = 1.
       ! Note: All vertices with vmask = 1 (i.e., all vertices with at least one ice-covered neighbor)
       !        are surrounded by four cells with physically meaningful values of f_flotation.

       do j = 1, ny-1
          do i = 1, nx-1

             if (vmask(i,j) == 1) then  ! at least one neighboring cell is ice-covered or land-based

                ! quadrant 1, southwest of the vertex
                f_flotation_vector(1) =          f_flotation(i,j)
                f_flotation_vector(2) = 0.5d0 * (f_flotation(i,j) + f_flotation(i+1,j))
                f_flotation_vector(3) = 0.25d0 * (f_flotation(i,j+1) + f_flotation(i+1,j+1)  &
                                                + f_flotation(i,j)   + f_flotation(i+1,j))
                f_flotation_vector(4) = 0.5d0 * (f_flotation(i,j) + f_flotation(i,j+1))

                call compute_grounded_fraction(i,     j,     1,      this_rank, &  ! '1' is the quadrant index
                                               itest, jtest, rtest,             &
                                               f_flotation_vector,              &
                                               f_ground_quadrant(1,i,j))

                ! quadrant 2, southeast of the vertex
                f_flotation_vector(1) = 0.5d0 * (f_flotation(i+1,j) + f_flotation(i,j))
                f_flotation_vector(2) =          f_flotation(i+1,j)
                f_flotation_vector(3) = 0.5d0 * (f_flotation(i+1,j) + f_flotation(i+1,j+1))
                f_flotation_vector(4) = 0.25d0 * (f_flotation(i,j+1) + f_flotation(i+1,j+1)  &
                                                + f_flotation(i,j)   + f_flotation(i+1,j))

                call compute_grounded_fraction(i,     j,     1,      this_rank, &  ! '1' is the quadrant index
                                               itest, jtest, rtest,             &
                                               f_flotation_vector,              &
                                               f_ground_quadrant(2,i,j))

                ! quadrant 3, northeast of the vertex
                f_flotation_vector(1) = 0.25d0 * (f_flotation(i,j+1) + f_flotation(i+1,j+1)  &
                                                + f_flotation(i,j)   + f_flotation(i+1,j))
                f_flotation_vector(2) = 0.5d0 * (f_flotation(i+1,j+1) + f_flotation(i,j+1))
                f_flotation_vector(3) =          f_flotation(i+1,j+1)
                f_flotation_vector(4) = 0.5d0 * (f_flotation(i+1,j+1) + f_flotation(i+1,j))

                call compute_grounded_fraction(i,     j,     1,      this_rank, &  ! '1' is the quadrant index
                                               itest, jtest, rtest,             &
                                               f_flotation_vector,              &
                                               f_ground_quadrant(3,i,j))

                ! quadrant 4, northwest of the vertex
                f_flotation_vector(1) = 0.5d0 * (f_flotation(i,j+1) + f_flotation(i,j))
                f_flotation_vector(2) = 0.25d0 * (f_flotation(i,j+1) + f_flotation(i+1,j+1)  &
                                                + f_flotation(i,j)   + f_flotation(i+1,j))
                f_flotation_vector(3) = 0.5d0 * (f_flotation(i,j+1) + f_flotation(i+1,j+1))
                f_flotation_vector(4) =          f_flotation(i,j+1)

                call compute_grounded_fraction(i,     j,     1,      this_rank, &  ! '1' is the quadrant index
                                               itest, jtest, rtest,             &
                                               f_flotation_vector,              &
                                               f_ground_quadrant(4,i,j))

                ! average values for the 4 quadrants surrounding the vertex
                f_ground(i,j) = 0.25d0 * (f_ground_quadrant(1,i,j) + f_ground_quadrant(2,i,j)  &
                                        + f_ground_quadrant(3,i,j) + f_ground_quadrant(4,i,j))

                !WHL - debug
                if (verbose_gl .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                   print*, ' '
                   print*, 'f_ground at vertex, r, i, j =', this_rank, i, j
                   print*, 'Quadrant 1:', f_ground_quadrant(1,i,j)
                   print*, 'Quadrant 2:', f_ground_quadrant(2,i,j)
                   print*, 'Quadrant 3:', f_ground_quadrant(3,i,j)
                   print*, 'Quadrant 4:', f_ground_quadrant(4,i,j)
                   print*, 'Average   :', f_ground(i,j)
                endif

             endif        ! vmask = 1
          enddo           ! i
       enddo              ! j

       !----------------------------------------------------------------
       ! Compute f_ground_cell by summing over f_ground_quadrant computed above
       !----------------------------------------------------------------

       ! Loop over cells, summing over f_ground_quadrant
       ! Note: Above, we computed 4 values of f_ground_quadrant for all vertices with at least one
       !       ice-covered or land-based neighbor.
       !       Thus, if a vertex has vmask = 0, it is surrounded by ice-free ocean with f_ground = 0.
       !       So in summing f_ground over the 4 quadrants of each cell, we can ignore quadrants
       !        associated with vertices that have vmask = 0.

       do j = 2, ny-1
          do i = 2, nx-1

             ! SW quadrant: quadrant 3 of vertex (i-1,j-1)
             if (vmask(i-1,j-1) == 1) then
                f_ground_cell(i,j) = f_ground_cell(i,j) + f_ground_quadrant(3,i-1,j-1)
             endif

             ! SE quadrant: quadrant 4 of vertex (i,j-1)
             if (vmask(i,j-1) == 1) then
                f_ground_cell(i,j) = f_ground_cell(i,j) + f_ground_quadrant(4,i,j-1)
             endif

             ! NE quadrant: quadrant 1 of vertex (i,j)
             if (vmask(i,j) == 1) then
                f_ground_cell(i,j) = f_ground_cell(i,j) + f_ground_quadrant(1,i,j)
             endif

             ! NW quadrant: quadrant 2 of vertex (i-1,j)
             if (vmask(i-1,j) == 1) then
                f_ground_cell(i,j) = f_ground_cell(i,j) + f_ground_quadrant(2,i-1,j)
             endif

             ! divide by 4 to get average value of f_ground in the cell
             f_ground_cell(i,j) = 0.25d0 * f_ground_cell(i,j)

             !WHL - debug
             if (verbose_gl .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'f_ground_vertex, r, i, j =', this_rank, i, j
                print*, 'Quadrant 1:', f_ground_quadrant(3,i-1,j-1)
                print*, 'Quadrant 2:', f_ground_quadrant(4,i,j-1)
                print*, 'Quadrant 3:', f_ground_quadrant(1,i,j)
                print*, 'Quadrant 4:', f_ground_quadrant(2,i-1,j)
                print*, 'Average   :', f_ground_cell(i,j)
             endif

          enddo
       enddo

       ! Set f_ground_cell = 1 on land
       where (land_mask == 1)
          f_ground_cell = 1.0d0
       endwhere

       call parallel_halo(f_ground_cell)

    endif   ! which_ho_ground

  end subroutine glissade_grounded_fraction

!=======================================================================

  subroutine compute_grounded_fraction(i,     j,     q,       rank, &
                                       itest, jtest, rtest,         &
                                       f_flotation,  f_ground)

    ! Compute the grounded fraction in a cell quadrant, given f_flotation
    !  at the four corners of the quadrant.
    ! The quadrant can be a grid cell, a staggered cell  around a vertex,
    !  one-quarter of a grid cell, or any other rectangular region.

    integer, intent(in) ::   &
         i, j, q, rank                          !> cell and quadrant indices for this call (diagnostic only)

    integer, intent(in) ::   &
         itest, jtest, rtest                    !> indices of diagnostic cell

    real(dp), dimension(4), intent(in) ::  &
         f_flotation                            !> flotation function at 4 corners of quadrant;
                                                !> should be ordered counterclockwise around the quadrant 

    real(dp), intent(out) ::  &
         f_ground                               !> fraction of grounded ice

    !----------------------------------------------------------------
    ! local variables
    !----------------------------------------------------------------

    integer :: nc  ! corner counting index

    real(dp) :: a, b, c, d       ! coefficients in bilinear interpolation
                                 ! f(x,y) = a + b*x + c*y + d*x*y

    real(dp) :: f1, f2, f3, f4   ! f_flotation at different corners

    real(dp) ::  &
         var                     ! combination of f_flotation terms that determines regions to be integrated

    logical, dimension(4) ::   &
         cfloat                  ! true if flotation condition is satisfied at a corner, else = false

    integer :: nfloat            ! number of grounded vertices of the quadrant (0 to 4)

    logical, dimension(4) ::   &
         logvar                  ! set locally to float or .not.float, depending on nfloat

    real(dp) ::   &
         f_corner, &             ! fractional area in a corner region of the quadrant
         f_corner1, f_corner2, &
         f_trapezoid             ! fractional area in a trapezoidal region of the quadrant

    logical :: adjacent          ! true if two grounded corners are adjacent (rather than opposite)

    logical :: rotated           ! true if a pattern is rotated (used when 2 non-adjacent cells are G, and the other 2 are F)

    real(dp), parameter :: &
         eps06 = 1.d-06          ! small number

    ! Note: By convention, corners are numbered from 1 to 4 proceeding CCW from the southwest corner.
    !       The algorithm will work, however, for any CCW ordering of the input f_flotation array.
    !
    !       4-----3
    !       |     |
    !       |     |
    !       1-----2

    ! Identify corners where the ice is floating (based on f_flotation > 0).
    ! Count the number of floating corners.
    ! Note: Ice-free cells adjacent to floating cells (and not adjacent to grounded cells)
    !        will have f_flotation > 0 as a result of the f_flotation extrapolation above.

    nfloat = 0
    do nc = 1, 4
       if (f_flotation(nc) > 0.0d0) then
          cfloat(nc) = .true.
          nfloat = nfloat + 1
       else
          cfloat(nc) = .false.
       endif
    enddo

    !WHL - debug
    if (verbose_gl .and. i == itest .and. j==jtest .and. rank == rtest) then
       print*, ' '
       print*, 'rank, i, j =', rank, i, j
       print*, 'f_flotation(4:3):', f_flotation(4), f_flotation(3)
       print*, 'f_flotation(1:2):', f_flotation(1), f_flotation(2)
       print*, 'nfloat =', nfloat
    endif

    ! Given nfloat, compute f_ground for each vertex
    ! First the easy cases...

    if (nfloat == 0) then

       f_ground = 1.0d0    ! fully grounded

    elseif (nfloat == 4) then

       f_ground = 0.0d0    ! fully floating

       ! For the other cases the grounding line runs through the rectangular region
       !  defined by the 4 corners.
       ! Using the values at each corner, we approximate f_flotation(x,y) as
       !  a bilinear function f(x,y) = a + bx + cy + dxy over the region.
       ! Assume the grounding line is defined by f(x,y) = 0, with f <= 0 for grounded ice
       !  and f > 0 for floating ice.
       ! To find f_ground, we integrate over the region with f(x,y) <= 0
       !  (or alternatively, we find f_float = 1 - f_ground by integrating
       !  over the region with f(x,y) > 0).
       !  
       ! There are 3 patterns to consider:
       ! (1) nfloat = 1 or nfloat = 3 (one corner is not like the others)
       ! (2) nfloat = 2, and adjacent corners are floating
       ! (3) nfloat = 2, and diagonally opposite corners are floating

    elseif (nfloat == 1 .or. nfloat == 3) then

       if (nfloat==1) then
          logvar(1:4) = cfloat(1:4)
       else  ! nfloat = 3
          logvar(1:4) = .not.cfloat(1:4)
       endif

       ! Identify the corner that is not like the others
       ! (i.e., the only floating corner if nfloat = 1, or the only grounded corner if nfloat = 3)
       !
       ! Diagrams below are for the case nfloat = 1.
       ! If nfloat = 3, the F and G labels are switched.

       if (logvar(1)) then             ! no rotation
          f1 = f_flotation(1)          !   G-----G
          f2 = f_flotation(2)          !   |     |
          f3 = f_flotation(3)          !   |     |
          f4 = f_flotation(4)          !   F-----G

       elseif (logvar(2)) then         ! rotate by 90 degrees
          f4 = f_flotation(1)          !   G-----G
          f1 = f_flotation(2)          !   |     |
          f2 = f_flotation(3)          !   |     |
          f3 = f_flotation(4)          !   G-----F

       elseif (logvar(3)) then         ! rotate by 180 degrees
          f3 = f_flotation(1)          !   G-----F
          f4 = f_flotation(2)          !   |     |
          f1 = f_flotation(3)          !   |     |
          f2 = f_flotation(4)          !   G-----G

       elseif (logvar(4)) then         ! rotate by 270 degrees
          f2 = f_flotation(1)          !   F-----G
          f3 = f_flotation(2)          !   |     |
          f4 = f_flotation(3)          !   |     |
          f1 = f_flotation(4)          !   G-----G
       endif

       ! Compute coefficients in f(x,y) = a + b*x + c*y + d*x*y
       ! Note: x is to the right and y is up if the southwest corner is not like the others.
       !       For the other cases we solve the same problem with x and y rotated.
       !       The rotations are handled by rotating f1, f2, f3 and f4 above.

       a = f1
       b = f2 - f1
       c = f4 - f1
       d = f1 + f3 - f2 - f4

       !WHL - debug
       if (verbose_gl .and. i==itest .and. j==jtest .and. rank == rtest) then
          print*, 'f1, f2, f3, f4 =', f1, f2, f3, f4
          print*, 'a, b, c, d =', a, b, c, d
       endif

       ! Compute the fractional area of the corner region
       ! (floating if nfloat = 1, grounded if nfloat = 3)
       !
       ! Here are the relevant integrals:
       !
       ! (1) d /= 0:
       !     integral_0^x0 {y(x) dx}, where x0   = -a/b
       !                                    y(x) = -(a+bx) / (c+dx)
       !     = [(bc - ad) ln[1 - (ad)/(bc)] + ad] / d^2
       !
       ! (2) d = 0:
       !     integral_0^x0 {y(x) dx}, where x0   = -a/b
       !                                    y(x) = -(a+bx) / c
       !     = a^2 / (2bc)
       !
       ! Note: We cannot have bc = 0, because f_flotation varies in both x and y.
       !       The above rotations ensure that we always take the log of a positive number.
       ! Note: This expression will give a NaN if f_flotation = 0 for land cells.
       !       Thus, f_flotation must be < 0 for land, even if topg - eus = 0.
       if (abs((a*d)/(b*c)) > eps06) then
          f_corner = ((b*c - a*d) * log(abs(1.d0 - (a*d)/(b*c))) + a*d) / (d*d)
       else
          f_corner = (a*a) / (2.0d0*b*c)
       endif

       if (nfloat==1) then  ! f_corner is the floating area
          f_ground = 1.0d0 - f_corner
       else                 ! f_corner is the grounded area
          f_ground = f_corner
       endif

       !WHL - debug
       if (verbose_gl .and. i==itest .and. j==jtest .and. rank == rtest) then
          print*, 'f_corner =', f_corner
          print*, 'f_ground =', f_ground
       endif

    elseif (nfloat == 2) then

       ! first the 4 cases where the 2 grounded corners are adjacent
       ! We integrate over the trapezoid in the floating part of the quadrant

       if (cfloat(1) .and. cfloat(2)) then  ! no rotation
          adjacent = .true.              !   G-----G
          f1 = f_flotation(1)            !   |     |
          f2 = f_flotation(2)            !   |     |
          f3 = f_flotation(3)            !   |     |
          f4 = f_flotation(4)            !   F-----F

       elseif (cfloat(2) .and. cfloat(3)) then  ! rotate by 90 degrees
          adjacent = .true.              !   G-----F
          f4 = f_flotation(1)            !   |     |
          f1 = f_flotation(2)            !   |     |
          f2 = f_flotation(3)            !   |     |
          f3 = f_flotation(4)            !   G-----F

       elseif (cfloat(3) .and. cfloat(4)) then  ! rotate by 180 degrees
          adjacent = .true.              !   F-----F
          f3 = f_flotation(1)            !   |     |
          f4 = f_flotation(2)            !   |     |
          f1 = f_flotation(3)            !   |     |
          f2 = f_flotation(4)            !   G-----G

       elseif (cfloat(4) .and. cfloat(1)) then   ! rotate by 270 degrees
          adjacent = .true.              !   F-----G
          f2 = f_flotation(1)            !   |     |
          f3 = f_flotation(2)            !   |     |
          f4 = f_flotation(3)            !   |     |
          f1 = f_flotation(4)            !   F-----G

       else   ! the 2 grounded corners are diagonally opposite

          adjacent = .false.

          ! We will integrate assuming the two corner regions lie in the lower left
          ! and upper right, i.e. one of these patterns:
          !
          !   F---\-G       G---\-F
          !   |    \|       |    \|
          !   |\    \       \     \
          !   | \   |       |\    |
          !   G--\--F       F-\---G
          !
          ! Two other patterns are possible, with corner regions in the lower right
          ! and upper left; these require a rotation before integrating: 
          !
          !   G-/---F       F-/---G
          !   |/    |       |/    |
          !   /    /|       /    /|
          !   |   / |       |   / |
          !   F--/--G       G--/--F

          var = f_flotation(2)*f_flotation(4) - f_flotation(1)*f_flotation(3)

          if (var >= 0.d0) then   ! we have one of the top two patterns
             f1 = f_flotation(1)
             f2 = f_flotation(2)
             f3 = f_flotation(3)
             f4 = f_flotation(4)
             rotated = .false.
          else   ! we have one of the bottom two patterns; rotate coordinates by 90 degrees
             f4 = f_flotation(1)
             f1 = f_flotation(2)
             f2 = f_flotation(3)
             f3 = f_flotation(4)
             rotated = .true.
          endif

       endif  ! grounded corners are adjacent

       ! Compute coefficients in f(x,y) = a + b*x + c*y + d*x*y
       a = f1
       b = f2 - f1
       c = f4 - f1
       d = f1 + f3 - f2 - f4

       ! Integrate the corner areas

       !WHL - debug
       if (verbose_gl .and. i==itest .and. j==jtest .and. rank == rtest) then
          print*, 'adjacent =', adjacent
          print*, 'f1, f2, f3, f4 =', f1, f2, f3, f4
          print*, 'a, b, c, d =', a, b, c, d
       endif

       if (adjacent) then

          ! Compute the area of the floating part of the quadrant.
          ! Here are the relevant integrals:
          !
          ! (1) d /= 0:
          !     integral_0^1 {y(x) dx}, where y(x) = -(a+bx) / (c+dx)
          !     = [(bc - ad) ln(1 + d/c) - bd] / d^2
          !
          ! (2) d = 0:
          !     integral_0^1 {y(x) dx}, where y(x) = -(a+bx) / c
          !     = -(2a + b) / (2c)
          !
          ! Note: We cannot have c = 0, because the passage of the GL
          !       through the region from left to right implies variation in y.
          !       The above rotations ensure that we always take the log of a positive number

          if (abs(d/c) > eps06) then   ! the usual case
             f_trapezoid = ((b*c - a*d) * log(1.d0 + d/c) - b*d) / (d*d)
          else
             f_trapezoid = -(2.d0*a + b) / (2.d0*c)
          endif

          f_ground = 1.d0 - f_trapezoid

          !WHL - debug
          if (verbose_gl .and. i==itest .and. j==jtest .and. rank == rtest) then
             print*, 'f_trapezoid =', f_trapezoid
             print*, 'f_ground =', f_ground
          endif

       else   ! grounded corners are diagonally opposite

          ! bug check: make sure some signs are positive as required by the formulas
          if (b*c - a*d < 0.d0) then
             print*, 'Grounding line error: bc - ad < 0'
             print*, 'rank, i, j, q =', rank, i, j, q
             stop
          elseif ((b+d)*(c+d) < 0.d0) then
             print*, 'Grounding line error: (b+d)(c+d) < 0'
             print*, 'rank, i, j, q =', rank, i, j, q
             stop
          endif

          ! Compute the combined areas of the two corner regions.
          ! For the lower left region, the integral is the same as above
          ! (for the case nfloat = 1 or nfloat = 3, with d /= 0).
          ! For the upper right region, here is the integral:
          !
          !     integral_x1^1 {(1-y(x)) dx}, where x1  = -(a+c)/(b+d)
          !                                       y(x) = -(a+b*x) / (c+d*x)
          !     = {(bc - ad) ln[(bc - ad)/((b+d)(c+d))] + d(a + b + c + d)} / d^2
          !
          ! The above rotations ensure that bc - ad >= 0.
          ! But the integral is valid only if bc - ad > 0.
          ! If bc - ad = 0, then the grounding line lies along two lines,
          !  x0 = -a/b and y0 = -a/c.
          ! The lower left area is x0*y0 = a^2 / (bc).
          ! The upper right area is (1-x0)*(1-y0) = (a+b)(a+c) / (bc)
          !
          ! Note that this pattern is not possible with d = 0.

          !WHL - debug
          if (verbose_gl .and. i==itest .and. j==jtest .and. rank == rtest) then
             print*, 'Pattern 3: i, j, bc - ad =', i, j, b*c - a*d
          endif

          if (abs(b*c - a*d) > eps06) then  ! the usual case
             f_corner1 = ((b*c - a*d) * log(1.d0 - (a*d)/(b*c)) + a*d) / (d*d)
             f_corner2 = ((b*c - a*d) * log((b*c - a*d)/((b+d)*(c+d)))  &
                  + d*(a + b + c + d)) / (d*d)
          else 
             f_corner1 = a*a / (b*c)
             f_corner2 = (a + b)*(a + c) / (b*c)
          endif

          ! Determine whether the two corner regions are grounded or floating,
          ! and set f_ground accordingly.

          if (f_flotation(1) > 0.0d0) then  ! SW corner is floating
             if (rotated) then   ! corner regions (with areas f_corner1 and f_corner2) are grounded
                f_ground = f_corner1 + f_corner2
             else   ! corner regions are floating
                f_ground = 1.d0 - (f_corner1 + f_corner2)
             endif
          else  ! SW corner is grounded
             if (rotated) then   ! corner regions are floating
                f_ground = 1.d0 - (f_corner1 + f_corner2)
             else   ! corner regions are grounded
                f_ground = f_corner1 + f_corner2
             endif
          endif

          !WHL - debug
          if (verbose_gl .and. i==itest .and. j==jtest .and. rank == rtest) then
             print*, 'f_corner1 =', f_corner1
             print*, 'f_corner2 =', f_corner2
             print*, 'f_ground =', f_ground
          endif

       endif  ! adjacent or opposite

    endif     ! nfloat

  end subroutine compute_grounded_fraction

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

