!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grid_operators.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains various grid operators for the Glissade dycore, including routines 
! for computing gradients and interpolating between staggered and unstaggered grids.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glissade_grid_operators

    use glimmer_global, only: dp
    use glimmer_log
    use glide_types  ! HO_GRADIENT_MARGIN_*
    use parallel

    implicit none

    private
    public :: glissade_stagger, glissade_unstagger,    &
              glissade_centered_gradient, glissade_upstream_gradient,    &
              glissade_gradient_at_edges,  &
              glissade_surface_elevation_gradient,  &
              glissade_vertical_average

    logical, parameter :: verbose_gradient = .false.

contains

!----------------------------------------------------------------------------

  subroutine glissade_stagger(nx,           ny,        &
                              var,          stagvar,   &
                              ice_mask,     stagger_margin_in)

    !----------------------------------------------------------------
    ! Given a variable on the unstaggered grid (dimension nx, ny), interpolate
    ! to find values on the staggered grid (dimension nx-1, ny-1).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) ::    &
       var                      ! unstaggered field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    integer, dimension(nx,ny), intent(in), optional ::        &
       ice_mask                 ! = 1 where values are included in the average, else = 0
                                ! Typically ice_mask = 1 where ice is present (or thck > thklim), else = 0
                                ! Note: ice_mask is not needed if stagger_margin = 0

    integer, intent(in), optional ::   &
       stagger_margin_in        ! 0 = use all values when interpolating
                                !   may be appropriate when computing stagusrf and stagthck on land
                                ! 1 = use only values where ice_mask = 1
                                !   preferable for tracers (e.g., temperature, flwa) and ocean margins

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sumvar, summask
    integer :: stagger_margin

    if (present(stagger_margin_in)) then
       stagger_margin = stagger_margin_in
    else
       stagger_margin = 0  ! default is to average over all cells, including those where ice is absent
    endif

    if (stagger_margin == 1 .and. .not.present(ice_mask)) then
       call write_log('Must pass in ice_mask to compute staggered field with stagger_margin = 1', GM_FATAL)
    endif

    stagvar(:,:) = 0.d0

    if (stagger_margin == 0) then

       ! Average over all four neighboring cells

       do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          stagvar(i,j) = (var(i,j+1) + var(i+1,j+1) + var(i,j) + var(i+1,j)) / 4.d0
       enddo
       enddo  

    elseif (stagger_margin == 1) then

       ! Average over cells where ice_mask = 1

       do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          sumvar = ice_mask(i,j+1)*var(i,j+1) + ice_mask(i+1,j+1)*var(i+1,j+1)  &
                 + ice_mask(i,j)  *var(i,j)   + ice_mask(i+1,j)  *var(i+1,j)
          summask = real(ice_mask(i,j+1) + ice_mask(i+1,j+1) + ice_mask(i,j) + ice_mask(i+1,j), dp)
          if (summask > 0.d0) stagvar(i,j) = sumvar / summask
       enddo
       enddo  

    endif

  end subroutine glissade_stagger

!----------------------------------------------------------------------------

  subroutine glissade_unstagger(nx,           ny,          &
                                stagvar,      unstagvar,   &
                                vmask,        stagger_margin_in)

    !----------------------------------------------------------------
    ! Given a variable on the staggered grid (dimension nx-1, ny-1), interpolate
    ! to find values on the staggered grid (dimension nx, ny).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx-1,ny-1), intent(in) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    real(dp), dimension(nx,ny), intent(out) ::    &
       unstagvar                ! unstaggered field, defined at cell centers

    integer, dimension(nx-1,ny-1), intent(in), optional  ::        &
       vmask                    ! = 1 for vertices where the value is used in the average, else = 0
                                ! Note: The user needs to compute this mask in the calling subroutine.
                                !       It will likely be based on the scalar ice mask, but the details are left open.

    integer, intent(in), optional ::   &
       stagger_margin_in        ! 0 = use all values when interpolating
                                ! 1 = use only values where vmask = 1

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sumvar, summask
    integer :: stagger_margin

    if (present(stagger_margin_in)) then
       stagger_margin = stagger_margin_in
    else
       stagger_margin = 0  ! default is to average over all cells, including those where ice is absent
    endif

    if (stagger_margin == 1 .and. .not.present(vmask)) then
       call write_log('Must pass in vmask to compute unstaggered field with stagger_margin = 1', GM_FATAL)
    endif

    unstagvar(:,:) = 0.d0

    if (stagger_margin == 0) then

       ! Average over all four neighboring cells

       do j = 2, ny-1   ! loop does not include outer row of cells
       do i = 2, nx-1
          unstagvar(i,j) = (stagvar(i,j) + stagvar(i-1,j) + stagvar(i,j-1) + stagvar(i-1,j-1)) / 4.d0
       enddo
       enddo  

    elseif (stagger_margin == 1) then

       ! Average over vertices with vmask = 1

       do j = 2, ny-1   ! loop does not include outer row of cells
       do i = 2, nx-1
          sumvar = vmask(i-1,j)  *stagvar(i-1,j)   + vmask(i,j)  *stagvar(i,j)  &
                 + vmask(i-1,j-1)*stagvar(i-1,j-1) + vmask(i,j-1)*stagvar(i,j-1)  
          summask = real(vmask(i-1,j) + vmask(i,j) + vmask(i-1,j-1) + vmask(i,j-1), dp)
          if (summask > 0.d0) unstagvar(i,j) = sumvar / summask
       enddo
       enddo  

    endif

    ! Fill in halo values
    call parallel_halo(unstagvar)

  end subroutine glissade_unstagger

!****************************************************************************

  subroutine glissade_centered_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        field,                   &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        gradient_margin_in,      &
                                        usrf,                    &
                                        floating_mask,           &
                                        land_mask,               &
                                        max_slope)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient is evaluated at the four neighboring vertices and is second-order accurate.
    !
    ! The gradient at a given vertex is constructed from gradients at adjacent edges.
    ! If one or more edge gradients is masked out, then df_fx and df_dy are constructed from the others.
    ! If all edge gradients adjacent to a vertex are masked out, then the gradient is set to zero.
    !
    ! Edge gradients are computed in the standard way, taking the difference between
    !  the values in two adjacent cells and dividing by the distance.
    ! At the ice margin, where one cell adjacent to a given edge is ice-free,
    !  edge gradients may be masked in the following ways:
    !
    ! HO_GRADIENT_MARGIN_ALL = 0: Values in both adjacent cells are used to compute the gradient,
    !  including values in ice-free cells.  In other words, there is no masking of edges.
    !  This convention is used by Glide. It works well at land-terminating margins, but performs poorly
    !  for ice shelves with a sudden drop in ice thickness and surface elevation at the margin.
    !
    ! HO_GRADIENT_MARGIN_GROUNDED_ICE = 1: The gradient is computed at edges where either
    ! (1) Both adjacent cells are ice-covered.
    ! (2) One cell is ice-covered and grounded, and lies above the other cell.
    !
    ! The edge is masked out where a floating cell is adjacent to an ice-free ocean cell,
    ! or where an ice-covered cell lies below an ice-free land cell (i.e., a nunatak).
    ! At land-terminating margins the gradient is nonzero (except for the nunatak case),
    ! and at marine-terminating margins the gradient is zero (unless the ice-covered cell is grounded).
    !
    ! This option used to be the default. However, it can overestimate driving forces at grounded marine
    !  margins where spreading is driven by lateral forces even when the surface gradient is zero.
    !
    ! HO_GRADIENT_MARGIN_ICE_ONLY = 2: Only values in ice-covered cells (i.e., cells with thck > thklim)
    !  are used to compute gradients.  If one or both adjacent cells is ice-free, the edge is masked out.
    !  This option works well at shelf margins but less well for land margins (e.g., the Halfar test case).
    !
    ! HO_GRADIENT_MARGIN_ICE_OVER_LAND = 3: The gradient is computed at edges where either
    ! (1) Both adjacent cells are ice-covered.
    ! (2) One cell is ice-covered and lies above ice-free land.
    !
    ! Method 3 does not compute a gradient at edges where ice lies above ocean cells.
    ! At these edges, lateral forces drive spreading in the absence of a surface gradient.
    ! This method works well at both land- and marine-terminating margins.
    ! It is the default for higher-order simulations.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       field                    ! input scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components of input field, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: Compute edge gradient when either cell is ice-covered
                                ! 1: Compute edge gradient for grounded ice above ice-free land or ocean
                                ! 2: Compute edge gradient only when both cells have ice
                                ! 3: Compute edge gradient for ice-covered cell above ice-free land (not ocean)

    real(dp), dimension(nx,ny), intent(in), optional ::       &
       usrf                     ! ice surface elevation)

    integer, dimension(nx,ny), intent(in), optional ::        &
       floating_mask,         & ! = 1 for cells where ice is present and floating, else = 0
       land_mask                ! = 1 for land cells, else = 0

    real(dp), intent(in), optional :: &
       max_slope                ! maximum slope allowed for surface gradient computations (unitless)

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: gradient_margin

    integer :: i, j

    logical, dimension(nx-1,ny) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny-1) ::    &
       edge_mask_y               ! edge mask for computing df/dy

    real(dp) :: df_dx_north, df_dx_south  ! df_dx at neighboring edges
    real(dp) :: df_dy_east, df_dy_west    ! df_dx at neighboring edges

    !WHL - debug
    real(dp) :: dfdx, dfdy
    integer :: edge_count

    !--------------------------------------------------------
    !   Gradient at vertex(i,j) is based on f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)
    !--------------------------------------------------------

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_GROUNDED_ICE
    endif

    ! Create masks identifying edges that will be used in gradient computations

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       edge_mask_x(:,:) = .true.       ! true for all edges
       edge_mask_y(:,:) = .true.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_GROUNDED_ICE) then

       if (present(floating_mask) .and. present(land_mask) .and. present(usrf)) then

          call glissade_edgemask_gradient_margin_grounded_ice(nx,          ny,         &
                                                              ice_mask,                &
                                                              floating_mask,           &
                                                              land_mask,               &
                                                              usrf,                    &
                                                              edge_mask_x, edge_mask_y)
       else
          call write_log('Must pass in floating_mask and usrf to use this gradient_margin option', GM_FATAL)
       endif   ! present(floating_mask), etc.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       ! mask for east and west cell edges
       do j = 1, ny
          do i = 1, nx-1
             if (ice_mask(i,j)==1  .and. ice_mask(i+1,j)==1) then
                edge_mask_x(i,j) = .true.
             else
                edge_mask_x(i,j) = .false.
             endif
          enddo
       enddo
       
       ! mask for north and south edges
       do j = 1, ny-1
          do i = 1, nx
             if (ice_mask(i,j)==1  .and. ice_mask(i,j+1)==1) then
                edge_mask_y(i,j) = .true.
             else
                edge_mask_y(i,j) = .false.
             endif
          enddo
       enddo

    endif  ! gradient_margin

    !WHL - debug - Count number of edges with gradients exceeding max_slope
    if (present(max_slope)) then
       edge_count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             dfdx = (field(i+1,j) - field(i,j)) / dx
             if (abs(dfdx) > max_slope .and. edge_mask_x(i,j)) then
                edge_count = edge_count + 1
             endif
             dfdy = (field(i,j+1) - field(i,j)) / dy
             if (abs(dfdy) > max_slope .and. edge_mask_y(i,j)) then
                edge_count = edge_count + 1
             endif
          enddo
       enddo
       edge_count = parallel_reduce_sum(edge_count)
       if (main_task) then
!          print*, 'Number of edges:', (nx-2*nhalo)*(ny-2*nhalo)*2
!          print*, 'Limit slope: edge_count =', edge_count
       endif
    endif

    ! compute gradient at vertices by averaging gradient at adjacent edges
    ! ignore edges with edge_mask = 0

    do j = 1, ny-1
       do i = 1, nx-1

          ! df/dx
          df_dx_north = (field(i+1,j+1) - field(i,j+1))/ dx
          df_dx_south = (field(i+1,j)   - field(i,j))  / dx

          if (edge_mask_x(i,j) .and. edge_mask_x(i,j+1)) then
             df_dx(i,j) = (df_dx_north + df_dx_south) / 2.d0
          elseif (edge_mask_x(i,j)) then
             df_dx(i,j) = df_dx_south
          elseif (edge_mask_x(i,j+1)) then
             df_dx(i,j) = df_dx_north
          else
             df_dx(i,j) = 0.d0
          endif

          ! Optionally, limit df_dx

          if (present(max_slope)) then
             if (df_dx(i,j) > 0.0d0) then
                df_dx(i,j) = min(df_dx(i,j), max_slope)
             else
                df_dx(i,j) = max(df_dx(i,j), -max_slope)
             endif
          endif

          ! df/dy
          df_dy_east = (field(i+1,j+1) - field(i+1,j))/ dy
          df_dy_west = (field(i,j+1)   - field(i,j))  / dy

          if (edge_mask_y(i,j) .and. edge_mask_y(i+1,j)) then
             df_dy(i,j) = (df_dy_east + df_dy_west) / 2.d0
          elseif (edge_mask_y(i,j)) then
             df_dy(i,j) = df_dy_west
          elseif (edge_mask_y(i+1,j)) then
             df_dy(i,j) = df_dy_east
          else
             df_dy(i,j) = 0.d0
          endif

          ! Optionally, limit df_dy

          if (present(max_slope)) then
             if (df_dy(i,j) > 0.0d0) then
                df_dy(i,j) = min(df_dy(i,j), max_slope)
             else
                df_dy(i,j) = max(df_dy(i,j), -max_slope)
             endif
          endif

       enddo  ! i
    enddo     ! j
   
    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Centered gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
!!          do i = 1, nx-1
          do i = 1, nx/2
             write(6,'(f9.6)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       
       print*, ' '
       print*, 'df_dy:'
       do j = ny-1, 1, -1
!!          do i = 1, nx-1
          do i = 1, nx/2
             write(6,'(f9.6)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_centered_gradient

!****************************************************************************

  subroutine glissade_upstream_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        field,                   &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        usrf,                    &
                                        gradient_margin_in,      &
                                        accuracy_flag_in,        &
                                        floating_mask,           &
                                        land_mask,               &
                                        max_slope)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    !  compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient can be evaluated at one upstream edge (for first-order accuracy) 
    !  or at two upstream edges (for second-order accuracy).
    ! The reason to take a one-sided gradient is to damp checkerboard noise
    !  that often arises with a centered gradient.
    !
    ! Note: Upstream is defined by the direction of higher surface elevation.
    !  For df_dx, the edge gradients are upstream in the y direction,
    !  and for df_dy, the edge gradients are upstream in the x direction.
    !
    ! See comments in subroutine glissade_centered_gradient about the 
    !  various values of gradient_margin. 
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       field                    ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    real(dp), dimension(nx,ny), intent(in) ::       &
       usrf                     ! ice surface elevation (required to determine upstream direction)

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells); see details above
                                !    If one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: use values in ice-covered cells only
                                !    If one or more values is masked out, construct df_fx and df_dy from the others

    integer, intent(in), optional ::    &
       accuracy_flag_in         ! = 1 for 1st order, 2 for 2nd order

    integer, dimension(nx,ny), intent(in), optional ::        &
       floating_mask,        &  ! = 1 where ice is present and floating, else = 0
       land_mask                ! = 1 for land cells, else = 0
                                ! floating and land masks required for gradient_margin = HO_GRADIENT_MARGIN_GROUNDED_ICE

    real(dp), intent(in), optional :: &
       max_slope               ! maximum slope allowed for surface gradient computations (unitless)

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: gradient_margin, accuracy_flag
    integer :: i, j
    integer :: summask
    real(dp) :: sum1, sum2

    logical, dimension(nx-1,ny) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny-1) ::    &
       edge_mask_y               ! edge mask for computing df/dy

    real(dp) :: df_dx_north, df_dx_north2
    real(dp) :: df_dx_south, df_dx_south2
    real(dp) :: df_dy_east, df_dy_east2
    real(dp) :: df_dy_west, df_dy_west2

    !--------------------------------------------------------
    !   First-order upstream gradient at vertex(i,j) is based on two points out of f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)
    !
    !   Second-order gradient is based on four points in the upstream direction
    !--------------------------------------------------------

    if (present(accuracy_flag_in)) then
       accuracy_flag = accuracy_flag_in
    else
       accuracy_flag = 2   ! default to second-order
    endif

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_GROUNDED_ICE
    endif

    ! Set integer edge mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       edge_mask_x(:,:) = .true.       ! true for all edges
       edge_mask_y(:,:) = .true.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_GROUNDED_ICE) then

       if (present(floating_mask) .and. present(land_mask)) then

          call glissade_edgemask_gradient_margin_grounded_ice(nx,          ny,         &
                                                              ice_mask,                &
                                                              floating_mask,           &
                                                              land_mask,               &
                                                              usrf,                    &
                                                              edge_mask_x, edge_mask_y)
       else
          call write_log('Must pass in floating_mask to use this gradient_margin option', GM_FATAL)
       endif   ! present(floating_mask)

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       ! mask for east and west cell edges
       do j = 1, ny
          do i = 1, nx-1
             if (ice_mask(i,j)==1  .and. ice_mask(i+1,j)==1) then
                edge_mask_x(i,j) = .true.
             else
                edge_mask_x(i,j) = .false.
             endif
          enddo
       enddo
       
       ! mask for north and south cell edges
       do j = 1, ny-1
          do i = 1, nx
             if (ice_mask(i,j)==1  .and. ice_mask(i,j+1)==1) then
                edge_mask_y(i,j) = .true.
             else
                edge_mask_y(i,j) = .false.
             endif
          enddo
       enddo
       
    endif  ! gradient_margin

    if (accuracy_flag == 1) then   ! first-order accurate

       do j = 1, ny-1
          do i = 1, nx-1

             if (edge_mask_x(i,j) .or. edge_mask_x(i,j+1)) then

                ! Compute df/dx by taking upstream gradient
                df_dx_north = (field(i+1,j+1) - field(i,j+1)) / dx
                df_dx_south = (field(i+1,j) - field(i,j)) / dx

                sum1 = usrf(i+1,j+1) + usrf(i,j+1)
                sum2 = usrf(i+1,j) + usrf(i,j)

                if (sum1 > sum2) then   ! north is upstream; use north edge gradient if possible

                   if (edge_mask_x(i,j+1)) then
                      df_dx(i,j) = df_dx_north
                   else
                      df_dx(i,j) = df_dx_south
                   endif

                else  ! south is upstream; use south edge gradient if possible

                   if (edge_mask_x(i,j)) then
                      df_dx(i,j) = df_dx_south
                   else
                      df_dx(i,j) = df_dx_north
                   endif
                   
                endif   ! sum1 > sum2

             else    ! both adjacent edge masks = F; punt

                df_dx(i,j) = 0.d0

             endif   ! adjacent edge_mask = T

             ! Optionally, limit df_dx

             if (present(max_slope)) then
                if (df_dx(i,j) > 0.0d0) then
                   df_dx(i,j) = min(df_dx(i,j), max_slope)
                else
                   df_dx(i,j) = max(df_dx(i,j), -max_slope)
                endif
             endif

             if (edge_mask_y(i,j) .or. edge_mask_y(i+1,j)) then

                ! Compute df/dy by taking upstream gradient
                df_dy_east = (field(i+1,j+1) - field(i+1,j)) / dy
                df_dy_west = (field(i,j+1) - field(i,j)) / dy

                sum1 = usrf(i+1,j+1) + usrf(i+1,j)
                sum2 = usrf(i,j+1) + usrf(i,j)
             
                if (sum1 > sum2) then   ! east is upstream; use east edge gradient if possible

                   if (edge_mask_y(i+1,j)) then
                      df_dy(i,j) = df_dy_east
                   else
                      df_dy(i,j) = df_dy_west
                   endif

                else  ! west is upstream; use west edge gradient if possible

                   if (edge_mask_y(i,j)) then
                      df_dy(i,j) = df_dy_west
                   else
                      df_dy(i,j) = df_dy_east
                   endif

                endif   ! sum1 > sum2

             else    ! both adjacent edge masks = F; punt

                df_dy(i,j) = 0.d0

             endif   ! adjacent edge mask = T

             ! Optionally, limit df_dy

             if (present(max_slope)) then
                if (df_dy(i,j) > 0.0d0) then
                   df_dy(i,j) = min(df_dy(i,j), max_slope)
                else
                   df_dy(i,j) = max(df_dy(i,j), -max_slope)
                endif
             endif

          enddo
       enddo

    else    ! second-order accurate

       do j = 2, ny-2   ! loop does not include all of halo
          do i = 2, nx-2

             if (edge_mask_x(i,j) .or. edge_mask_x(i,j+1)) then

                ! Compute df_dx by taking upstream gradient

                df_dx_north2 = (field(i+1,j+2) - field(i,j+2)) / dx
                df_dx_north  = (field(i+1,j+1) - field(i,j+1)) / dx
                df_dx_south  = (field(i+1,j)   - field(i,j))   / dx
                df_dx_south2 = (field(i+1,j-1) - field(i,j-1)) / dx

                sum1 = usrf(i+1,j+1) + usrf(i,j+1) + usrf(i+1,j+2) + usrf(i,j+2)
                sum2 = usrf(i+1,j) + usrf(i,j) + usrf(i+1,j-1) + usrf(i,j-1)

                if (sum1 > sum2) then   ! north is upstream; use north edge gradients if possible

                   if (edge_mask_x(i,j+1) .and. edge_mask_x(i,j+2)) then
                      df_dx(i,j) = 1.5d0 * df_dx_north - 0.5d0 * df_dx_north2
                   elseif (edge_mask_x(i,j+1)) then   ! revert to first order
                      df_dx(i,j) = df_dx_north
                   else      ! first-order downstream
                      df_dx(i,j) = df_dx_south
                   endif
                                            
                else    ! south is upstream; use south edge gradients if possible

                   if (edge_mask_x(i,j) .and. edge_mask_x(i,j-1)) then
                      df_dx(i,j) = 1.5d0 * df_dx_south - 0.5d0 * df_dx_south2
                   elseif (edge_mask_x(i,j)) then   ! revert to first order
                      df_dx(i,j) = df_dx_south
                   else      ! first-order downstream
                      df_dx(i,j) = df_dx_north
                   endif

                endif   ! sum1 > sum2

             else   ! both adjacent edge masks = F; punt

                df_dx(i,j) = 0.d0

             endif  ! adjacent edge mask = T

             ! Optionally, limit df_dx

             if (present(max_slope)) then
                if (df_dx(i,j) > 0.0d0) then
                   df_dx(i,j) = min(df_dx(i,j), max_slope)
                else
                   df_dx(i,j) = max(df_dx(i,j), -max_slope)
                endif
             endif

             if (edge_mask_y(i,j) .or. edge_mask_y(i+1,j)) then

                ! Compute df_dy by taking upstream gradient

                df_dy_east2 = (field(i+2,j+1) - field(i+2,j)) / dy
                df_dy_east  = (field(i+1,j+1) - field(i+1,j)) / dy
                df_dy_west  = (field(i,j+1)   - field(i,j))   / dy
                df_dy_west2 = (field(i-1,j+1) - field(i-1,j)) / dy

                ! determine upstream direction

                sum1 = usrf(i+1,j+1) + usrf(i+1,j) + usrf(i+2,j+1) + usrf(i+2,j)
                sum2 = usrf(i,j+1) + usrf(i,j) + usrf(i-1,j+1) + usrf(i-1,j)
             
                if (sum1 > sum2) then  ! east is upstream; use east edge gradients if possible

                   if (edge_mask_y(i+1,j) .and. edge_mask_y(i+2,j)) then
                      df_dy(i,j) = 1.5d0 * df_dy_east - 0.5d0 * df_dy_east2
                   elseif (edge_mask_y(i+1,j)) then   ! revert to first order
                      df_dy(i,j) = df_dy_east
                   else    ! first-order downstream
                      df_dy(i,j) = df_dy_west
                   endif

                else   ! west is upstream; use west edge gradients if possible

                   if (edge_mask_y(i,j) .and. edge_mask_y(i-1,j)) then
                      df_dy(i,j) = 1.5d0 * df_dy_west - 0.5d0 * df_dy_west2
                   elseif (edge_mask_y(i,j)) then   ! revert to first order
                      df_dy(i,j) = df_dy_west
                   else    ! first_order downstream
                      df_dy(i,j) = df_dy_east
                   endif

                endif   ! sum1 > sum2

             else       ! both adjacent edge masks = F; punt

                df_dy(i,j) = 0.d0

             endif      ! adjacent edge mask = T

             ! Optionally, limit df_dy

             if (present(max_slope)) then
                if (df_dy(i,j) > 0.0d0) then
                   df_dy(i,j) = min(df_dy(i,j), max_slope)
                else
                   df_dy(i,j) = max(df_dy(i,j), -max_slope)
                endif
             endif

          enddo     ! i
       enddo     ! j

       ! fill in halo values
       call staggered_parallel_halo(df_dx)
       call staggered_parallel_halo(df_dy)

    endif   ! first or second order accurate

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'upstream df_dx:'
       do j = ny-2, 2, -1
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'upstream df_dy:'
       do j = ny-2, 2, -1
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo

    endif

  end subroutine glissade_upstream_gradient

!****************************************************************************

  subroutine glissade_gradient_at_edges(nx,           ny,        &
                                        dx,           dy,        &
                                        field,                   &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        gradient_margin_in,      &
                                        usrf,                    &
                                        floating_mask,           &
                                        land_mask,               &
                                        max_slope)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) at cell edges (i.e., the C grid):
    ! df_dx at the midpoint of the east edge and df_dy at the midpoint of
    ! the north edge.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       field                    ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny), intent(out) ::    &
       df_dx                    ! x gradient component, defined on east cell edge

    real(dp), dimension(nx,ny-1), intent(out) ::    &
       df_dy                    ! y gradient component, defined on north cell edge

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells); see details above
                                !    If one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: use values in ice-covered cells only
                                !    If one or more values is masked out, construct df_fx and df_dy from the others

    real(dp), dimension(nx,ny), intent(in), optional ::       &
       usrf                     ! ice surface elevation

    integer, dimension(nx,ny), intent(in), optional ::        &
       floating_mask,         & ! = 1 where ice is present and floating, else = 0
       land_mask                ! = 1 for land cells, else = 0
                                ! floating and land masks required for gradient_margin = HO_GRADIENT_MARGIN_GROUNDED_ICE

    real(dp), intent(in), optional :: &
       max_slope                ! maximum slope allowed for surface gradient computations (unitless)

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: gradient_margin
    integer :: i, j

    logical, dimension(nx-1,ny) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny-1) ::    &
       edge_mask_y               ! edge mask for computing df/dy

    real(dp) ::  &
       edge_thck_upper, edge_factor

    !--------------------------------------------------------
    !   Gradient at east edge(i,j) is based on f(i:i+1,j)
    !   Gradient at north edge(i,j) is based on f(i,j:j+1)
    !
    !   |             |
    !   |   (i,j+1)   |
    !   |             |
    !   |             |
    !   ----df_dy------------------
    !   |             |  
    !   |             |
    !   |   (i,j)   df_dx   (i+1,j)
    !   |             |
    !   |             |
    !   |--------------
    !
    !--------------------------------------------------------

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_GROUNDED_ICE
    endif

    ! Set integer edge mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       edge_mask_x(:,:) = .true.       ! true for all edges
       edge_mask_y(:,:) = .true.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_GROUNDED_ICE) then

       if (present(floating_mask) .and. present(land_mask) .and. present(usrf)) then

          call glissade_edgemask_gradient_margin_grounded_ice(nx,          ny,         &
                                                              ice_mask,                &
                                                              floating_mask,           &
                                                              land_mask,               &
                                                              usrf,                    &
                                                              edge_mask_x, edge_mask_y)
       else
          call write_log('Must pass in floating mask and usrf to use this gradient_margin option', GM_FATAL)
       endif   ! present(floating_mask), etc.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       ! mask for east and west cell edges
       do j = 1, ny
          do i = 1, nx-1
             if (ice_mask(i,j)==1  .and. ice_mask(i+1,j)==1) then
                edge_mask_x(i,j) = .true.
             else
                edge_mask_x(i,j) = .false.
             endif
          enddo
       enddo
       
       ! mask for north and south cell edges
       do j = 1, ny-1
          do i = 1, nx
             if (ice_mask(i,j)==1  .and. ice_mask(i,j+1)==1) then
                edge_mask_y(i,j) = .true.
             else
                edge_mask_y(i,j) = .false.
             endif
          enddo
       enddo
       
    endif  ! gradient_margin

    ! Compute the gradients where edge_mask = .true.

    ! df_dx
    do j = 1, ny
       do i = 1, nx-1
          if (edge_mask_x(i,j)) then
             df_dx(i,j) = (field(i+1,j) - field(i,j)) / dx
          else
             df_dx(i,j) = 0.d0
          endif
       enddo    ! i
    enddo       ! j

    ! df_dy
    do j = 1, ny-1
       do i = 1, nx
          if (edge_mask_y(i,j)) then
             df_dy(i,j) = (field(i,j+1) - field(i,j)) / dy
          else
             df_dy(i,j) = 0.d0
          endif
       enddo    ! i
    enddo       ! j

    if (present(max_slope)) then

       ! Optionally, limit df_dx
       do j = 1, ny
          do i = 1, nx-1
             if (df_dx(i,j) > 0.0d0) then
                df_dx(i,j) = min(df_dx(i,j), max_slope)
             else
                df_dx(i,j) = max(df_dx(i,j), -max_slope)
             endif
          enddo
       enddo
       
       ! Optionally, limit df_dy
       do j = 1, ny-1
          do i = 1, nx
             if (df_dy(i,j) > 0.0d0) then
                df_dy(i,j) = min(df_dy(i,j), max_slope)
             else
                df_dy(i,j) = max(df_dy(i,j), -max_slope)
             endif
          enddo
       enddo

    endif  ! present(max_slope)

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Edge gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
          do i = 1, nx
             write(6,'(f8.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'df_dy:'
       do j = ny, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_gradient_at_edges

!****************************************************************************

  subroutine glissade_edgemask_gradient_margin_grounded_ice(nx,          ny,         &
                                                            ice_mask,                &
                                                            floating_mask,           &
                                                            land_mask,               &
                                                            usrf,                    &
                                                            edge_mask_x, edge_mask_y)
    
    !----------------------------------------------------------------
    ! Compute edge masks required for option gradient_margin = HO_GRADIENT_MARGIN_GROUNDED_ICE.
    !
    ! The mask is set to true at all edges where either
    ! (1) Both adjacent cells are ice-covered.
    ! (2) One cell is ice-covered and grounded, and lies above the other cell.
    !
    ! The mask is set to false where a floating cell is adjacent to an ice-free ocean cell,
    ! or where an ice-covered land cell lies below an ice-free land cell (i.e., a nunatak).
    !
    ! The intent is to give a reasonable gradient at both land-terminating and marine-terminating margins.
    ! At land-terminating margins the gradient is nonzero (except for the nunatak case),
    ! and at marine-terminating margins the gradient is zero (unless the ice-covered cell is grounded).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions
   
    ! Note: land_mask is not currently part of the logic, but might be at some point
    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask,              & ! = 1 where ice is present, else = 0
       floating_mask,         & ! = 1 where ice is present and floating, else = 0
       land_mask                ! = 1 for land cells, else = 0

    real(dp), dimension(nx,ny), intent(in)  ::       &
       usrf                     ! ice surface elevation

    logical, dimension(nx-1,ny), intent(out) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny-1), intent(out) ::    &
       edge_mask_y               ! edge mask for computing df/dy

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    ! compute mask for east and west cell edges
    do j = 1, ny
       do i = 1, nx-1

          ! check for ice in both adjacent cells, or grounded ice in one cell lying above the other cell
          if ( (ice_mask(i,j)==1   .and.  ice_mask(i+1,j)==1)     .or.  &
               (ice_mask(i,j)==1   .and.  floating_mask(i,j)==0   .and. usrf(i,j)   > usrf(i+1,j)) .or.  &
               (ice_mask(i+1,j)==1 .and.  floating_mask(i+1,j)==0 .and. usrf(i+1,j) > usrf(i,j)) ) then

             edge_mask_x(i,j) = .true.

          else

             edge_mask_x(i,j) = .false.

          endif
       enddo
    enddo
    
    ! compute mask for north and south cell edges
    do j = 1, ny-1
       do i = 1, nx

          ! check for ice in both adjacent cells, or grounded ice in one cell lying above the other cell
          if ( (ice_mask(i,j)==1   .and.  ice_mask(i,j+1)==1)      .or.  &
               (ice_mask(i,j)==1   .and.  floating_mask(i,j)==0   .and. usrf(i,j)   > usrf(i,j+1)) .or.  &
               (ice_mask(i,j+1)==1 .and.  floating_mask(i,j+1)==0 .and. usrf(i,j+1) > usrf(i,j)) ) then

             edge_mask_y(i,j) = .true.

          else

             edge_mask_y(i,j) = .false.

          endif
       enddo
    enddo

  end subroutine glissade_edgemask_gradient_margin_grounded_ice

!****************************************************************************

  subroutine glissade_surface_elevation_gradient(nx,          ny,         &
                                                 dx,          dy,         &
                                                 active_ice_mask,         &
                                                 land_mask,               &
                                                 usrf,        thck,       &
                                                 topg,        eus,        &
                                                 thklim,                  &
                                                 thck_gradient_ramp,      &
                                                 ds_dx,       ds_dy,      &
                                                 max_slope)
    
    !----------------------------------------------------------------
    ! Compute surface elevation gradients for option gradient_margin = HO_GRADIENT_MARGIN_ICE_OVER_LAND.
    !
    ! The gradient is nonzero at edges where either
    ! (1) Both adjacent cells have active ice.
    ! (2) An ice-covered cell lies above ice-free land.
    !
    ! The gradient is set to zero where
    ! (1) An ice-covered cell (grounded or floating) lies above ice-free ocean.
    !     Note: Inactive calving-front cells are treated as ice-free ocean.
    ! (2) An ice-covered land cell lies below an ice-free land cell (i.e., a nunatak).
    !
    ! The intent is to give a reasonable gradient at both land-terminating and marine-terminating margins.
    ! At land-terminating margins the gradient is nonzero (except for nunataks), and at marine-terminating
    !  margins the gradient is zero.
    !
    ! At some edges, gradients are reduced to inhibit oscillations:
    ! - If the ice in the upper-lying cell is thin (thklim < thck < thklim + thck_gradient_ramp), the gradient is reduced.
    !   This prevents large changes in the gradient as the ice thickness oscillates about thklim.
    ! - If grounded marine-based ice lies above ice-free ocean, the surface elevation difference is replaced
    !   by thickness above flotation (thck - h_flotation), which approaches zero when the ice is nearly afloat. 
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
         nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
         dx, dy                   ! horizontal grid size

    integer, dimension(nx,ny), intent(in) ::        &
         active_ice_mask,       & ! = 1 where active ice is present, else = 0
         land_mask                ! = 1 for land cells, else = 0

    real(dp), dimension(nx,ny), intent(in) ::       &
         thck,                  & ! ice thickness
         usrf,                  & ! ice surface elevation
         topg                     ! bed elevation 

    real(dp), intent(in) ::       &
         eus,                   & ! eustatic sea level
         thklim,                & ! minimum thickness for active ice
         thck_gradient_ramp

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
         ds_dx,                & ! ds/dx at vertices
         ds_dy                   ! ds/dy at vertices

    real(dp), intent(in), optional ::       &
         max_slope               ! maximum slope allowed for surface gradient computations (unitless)

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j
    integer :: iu, il, ju, jl     ! indices of cells with higher and lower elevations

    real(dp), dimension(nx,ny) ::    &
         ds_dx_edge,            & ! ds/dx on east cell edges
         ds_dy_edge               ! ds/dy on north cell edges

    real(dp) ::  &
         edge_thck_upper,       & ! thickness of higher-elevation cell
         edge_thck_lower,       & ! thickness of lower-elevation cell
         edge_factor,           & ! gradient-weighting factor in range [0,1]
         sign_factor              ! sign factor, +1 or -1

    ! initialize
    ds_dx_edge(:,:) = 0.0d0
    ds_dy_edge(:,:) = 0.0d0

    ds_dx(:,:) = 0.0d0
    ds_dy(:,:) = 0.0d0

    ! compute ds_dx on east edges
    do j = 1, ny
       do i = 1, nx-1

          ! determine which cell is upper and which is lower
          if (usrf(i,j) > usrf(i+1,j)) then
             iu = i
             il = i+1
             sign_factor = -1.0d0
          else
             iu = i+1
             il = i
             sign_factor = 1.0d0
          endif

          if (land_mask(iu,j) == 1) then
             ! Compute a factor that reduces the gradient if ice in the upper cell is thin and land-based.
             ! This inhibits oscillations in the gradient when the thickness in the upper cell is close to thklim.
             edge_thck_upper = thck(iu,j)
             edge_factor = min(1.0d0, (edge_thck_upper - thklim)/thck_gradient_ramp)
             edge_factor = max(edge_factor, 0.0d0)
          else
             edge_factor = 1.0d0
          endif

          if (active_ice_mask(iu,j) == 1 .and. active_ice_mask(il,j) == 1) then  ! both cells have active ice

             ! compute the gradient in the usual way
             ds_dx_edge(i,j) = edge_factor * sign_factor * (usrf(iu,j) - usrf(il,j)) / dx

          elseif (active_ice_mask(iu,j) == 1 .and. land_mask(il,j) == 1) then

             ! upper cell has active ice and ice-free lower cell is land; compute the gradient in the usual way
             ds_dx_edge(i,j) = edge_factor * sign_factor * (usrf(iu,j) - usrf(il,j)) / dx

          endif  ! both cells have ice

       enddo   ! i
    enddo   ! j

    ! compute ds_dy on north edges
    do j = 1, ny-1
       do i = 1, nx

          ! determine which cell is upper and which is lower
          if (usrf(i,j) > usrf(i,j+1)) then
             ju = j
             jl = j+1
             sign_factor = -1.0d0
          else
             ju = j+1
             jl = j
             sign_factor = 1.0d0
          endif

          if (land_mask(i,ju) == 1) then
             ! Compute a factor that reduces the gradient if ice in the upper cell is thin and land-based.
             ! This inhibits oscillations in the gradient when the thickness in the upper cell is close to thklim.
             edge_thck_upper = thck(i,ju)
             edge_factor = min(1.0d0, (edge_thck_upper - thklim)/thck_gradient_ramp)
             edge_factor = max(edge_factor, 0.0d0)
          else
             edge_factor = 1.0d0
          endif

          if (active_ice_mask(i,ju)==1 .and. active_ice_mask(i,jl)==1) then  ! both cells have ice

             ! compute the gradient in the usual way
             ds_dy_edge(i,j) = edge_factor * sign_factor * (usrf(i,ju) - usrf(i,jl)) / dy

          elseif (active_ice_mask(i,ju) == 1 .and. land_mask(i,jl) == 1) then

             ! upper cell has active ice and ice-free lower cell is land; compute the gradient in the usual way
             ds_dy_edge(i,j) = edge_factor * sign_factor * (usrf(i,ju) - usrf(i,jl)) / dy

          endif  ! both cells have ice

       enddo  ! i
    enddo   ! j

    ! Average the edge gradients to vertices

    do j = 1, ny-1
       do i = 1, nx-1
          ds_dx(i,j) = 0.5d0 * (ds_dx_edge(i,j) + ds_dx_edge(i,j+1))
          ds_dy(i,j) = 0.5d0 * (ds_dy_edge(i,j) + ds_dy_edge(i+1,j))
       enddo
    enddo

    !WHL - This is an lternate method of averaging edge gradients to vertices, following glissade_centered_gradient.
    !       To be commented out.
    !      For the dome problem, the differences between gradient_margin methods 1 and 3
    !       come from replacing these lines with the ds_dx and ds_dy calculation above.
    do j = 1, ny-1
       do i = 1, nx-1
          if (ds_dx_edge(i,j) /= 0.0d0 .and. ds_dx_edge(i,j+1) /= 0.0d0) then
             ds_dx(i,j) = 0.5d0 * (ds_dx_edge(i,j) + ds_dx_edge(i,j+1))
          elseif (ds_dx_edge(i,j) /= 0.0d0) then
             ds_dx(i,j) = ds_dx_edge(i,j)
          elseif (ds_dx_edge(i,j+1) /= 0.0d0) then
             ds_dx(i,j) = ds_dx_edge(i,j+1)
          else
             ds_dx(i,j) = 0.d0
          endif
       enddo
    enddo

    do j = 1, ny-1
       do i = 1, nx-1
          if (ds_dy_edge(i,j) /= 0.0d0 .and. ds_dy_edge(i+1,j) /= 0.0d0) then
             ds_dy(i,j) = 0.5d0 * (ds_dy_edge(i,j) + ds_dy_edge(i+1,j))
          elseif (ds_dy_edge(i,j) /= 0.0d0) then
             ds_dy(i,j) = ds_dy_edge(i,j)
          elseif (ds_dy_edge(i+1,j) /= 0.0d0) then
             ds_dy(i,j) = ds_dy_edge(i+1,j)
          else
             ds_dy(i,j) = 0.d0
          endif
       enddo
    enddo


    ! Optionally, limit ds/dx and ds/dy

    if (present(max_slope)) then

       do j = 1, ny-1
          do i = 1, nx-1
             if (ds_dx(i,j) > 0.0d0) then
                ds_dx(i,j) = min(ds_dx(i,j), max_slope)
             else
                ds_dx(i,j) = max(ds_dx(i,j), -max_slope)
             endif
          enddo
       enddo

       do j = 1, ny-1
          do i = 1, nx-1
             if (ds_dy(i,j) > 0.0d0) then
                ds_dy(i,j) = min(ds_dy(i,j), max_slope)
             else
                ds_dy(i,j) = max(ds_dy(i,j), -max_slope)
             endif
          enddo
       enddo

    endif   ! present(max_slope)

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Hybrid gradient:'
       print*, ' '
       print*, 'ds_dx:'
       do j = ny-1, 1, -1
!!          do i = 1, nx-1
          do i = 1, nx/2
             write(6,'(f9.6)',advance='no') ds_dx(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'ds_dy:'
       do j = ny-1, 1, -1
!!          do i = 1, nx-1
          do i = 1, nx/2
             write(6,'(f9.6)',advance='no') ds_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_surface_elevation_gradient

!****************************************************************************

  subroutine glissade_vertical_average(nx,         ny,        &
                                       nz,         sigma,     &
                                       mask,                  &
                                       var,        var_2d)

    !----------------------------------------------------------------
    ! Compute the vertical average of a given variable.
    ! Note: It is assumed that the variable is defined at layer midpoints,
    !       and hence has vertical dimension (nz-1).
    ! Note: This subroutine will work for variables on the staggered
    !       horizontal grid if stagthck is passed in place of thck.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,            &     ! horizontal grid dimensions
       nz                       ! number of vertical levels

    real(dp), dimension(nz), intent(in) ::    &
       sigma                    ! sigma vertical coordinate

    logical, dimension(nx, ny), intent(in) ::    &
       mask                     ! compute var_2d where mask = .true.

    real(dp), dimension(nz-1,nx, ny), intent(in) ::    &
       var                      ! 3D field to be averaged vertically

    real(dp), dimension(nx, ny), intent(out) ::    &
       var_2d                   ! 2D vertically averaged field

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j, k

    do j = 1, ny
       do i = 1, nx

          var_2d(i,j) = 0.d0

          if (mask(i,j)) then
             do k = 1, nz-1
                var_2d(i,j) = var_2d(i,j) + var(k,i,j) * (sigma(k+1) - sigma(k))
             enddo
          endif

       enddo
    enddo

  end subroutine glissade_vertical_average

!****************************************************************************

  end module glissade_grid_operators

!****************************************************************************
