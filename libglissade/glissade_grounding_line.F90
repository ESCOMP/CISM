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
    use glimmer_log
    use parallel

    implicit none

    ! All subroutines in this module are public

    logical, parameter :: verbose_gl = .false.

  contains

!****************************************************************************

  subroutine glissade_grounded_fraction(nx,            ny,                      &
                                        itest, jtest,  rtest,                   &
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
    ! (0) HO_GROUND_NO_GLP: f_ground = 1 for vertices with grounded and/or land-based neighbor cells
    !                       f_ground = 0 for vertices with floating neighbors only
    ! (1) HO_GROUND_GLP: 0 <= f_ground <= 1 based on grounding-line parameterization
    !        A flotation function is interpolated over the bounding box of each vertex
    !        and analytically integrated to compute the grounded and floating fractions.
    ! (2) HO_GROUND_ALL: f_ground = 1 for all vertices with ice-covered neighbor cells
    !
    ! There are three options for the flotation function:
    ! (0) HO_FLOTATION_FUNCTION_PATTYN: f_flotation = (-rhoo*b)/(rhoi*H) - 1 = f_pattyn - 1
    !     Here, f_pattyn = (-rhoo*b)/(rhoi*H) as in Pattyn et al. (2006).
    ! (1) HO_FLOTATION_FUNCTION_INVERSE_PATTYN: f_flotation = 1 - (rhoi*H)/(-rhoo*b) = 1 - 1/f_pattyn
    ! (2) HO_FLOTATION_FUNCTION_LINEAR: f_flotation = -b - (rhoi/rhoo)*H = ocean cavity thickness
    !     This function was suggested by Xylar Asay-Davis and is linear in both b and H.
    ! All three functions are defined such that f <=0 for grounded ice and f > 0 for floating ice.
    ! For each option, land-based cells are assigned a large negative value, so that any vertices
    !  with land-based neighbors are strongly grounded.
    !
    ! NOTE: Option 1 = HO_GROUND_GLP is not suppported for this release.
    !       The flotation function does not enter f_ground calculations for options 0 and 2,
    !        but still is computed for diagnostic purposes.
    !
    ! We first compute f_flotation in all active ice-covered cells.
    ! Then f_flotation is extrapolated to ice-free neighbors.  Thus, f_flotation has a physically
    !   meaningful value (either computed directly, or extrapolated from a neighbor) in all four
    !   cells surrounding each active vertex. (By definition, an active vertex is a vertex with
    !   at least one active ice-covered neighbor.) Thus, we can interpolate f_flotation
    !   within the bounding box around each active vertex to compute f_ground at the vertex.
    !
    ! Until recently, the inverse Pattyn function (1) was the default.
    ! As of Nov. 2017, the linear function (2) is the default.
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
    ! The subroutine should work for other units as long as thck, topg and eus have the same units.

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

    real(dp), dimension(nx,ny), intent(out) :: &
       f_flotation            ! flotation function; see comments above

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
           
    integer :: i, j, ii, jj

    integer, dimension(nx-1,ny-1) ::   &
       vmask                     ! = 1 for vertices neighboring at least one cell where ice is present, else = 0

    real(dp) ::  &
       topg_eus_diff       ! topg - eus, limited to be >= f_flotation_land_topg_min

    logical, dimension(nx,ny) :: &
       cground             ! true if a cell is land and/or has grounded ice, else = false

    real(dp), parameter :: &
       f_flotation_land_topg_min = 1.0d0   ! min value of (topg - eus) in f_flotation expression for land cells (m)

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

    if (whichflotation_function == HO_FLOTATION_FUNCTION_PATTYN) then

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

    elseif (whichflotation_function == HO_FLOTATION_FUNCTION_INVERSE_PATTYN) then ! grounded if f_flotation >= 1, else floating

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

    elseif (whichflotation_function == HO_FLOTATION_FUNCTION_LINEAR) then

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

    endif  ! whichflotation_function

    ! initialize f_ground
    f_ground(:,:) = 0.0d0

    ! Compute f_ground according to the value of whichground

    select case(whichground)

    case(HO_GROUND_NO_GLP)   ! default: no grounding-line parameterization
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

       ! all vertices with ice-covered or land-based neighbors are assumed grounded, regardless of thck and topg

       do j = 1, ny-1
          do i = 1, nx-1
             if (vmask(i,j) == 1) then
                f_ground(i,j) = 1.d0
             endif
          enddo
       enddo

    case(HO_GROUND_GLP)      ! grounding-line parameterization

       call write_log('The GLP option for which_ho_ground is not supported for this release', GM_FATAL)

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

