!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grid_operators.F90 - part of the Community Ice Sheet Model (CISM)  
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
    use glide_types
    use parallel

    implicit none

    private
    public :: glissade_stagger, glissade_unstagger, glissade_stagger_real_mask, &
              glissade_stagger_edge,                &
              glissade_gradient, glissade_gradient_at_edges, &
              glissade_surface_elevation_gradient,  &
              glissade_slope_angle,                 &
              glissade_laplacian_smoother,          &
              glissade_vertical_average,            &
              glissade_vertical_interpolate

    logical, parameter :: verbose_gradient = .false.

contains

!----------------------------------------------------------------------------

  subroutine glissade_stagger(nx,           ny,        &
                              var,          stagvar,   &
                              ice_mask,     stagger_margin_in)

    !TODO - Drop the stagger_margin_in argument, and use the optional mask to
    !        determine where to ignore values when interpolating?

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

    ! Note: In most cases it would be safe here to do a staggered parallel halo update of stagvar.
    !       However, if running a problem with nonzero periodic_offset_ew/ns (e.g., the ISMIP-HOM and stream cases),
    !        a halo update would copy incorrect values to vertices on the west and south boundaries
    !        of the global domain. So there is no halo update here.
    !       To ensure correct halo values, the user should either pass in a 'var' field that has already
    !        been updated in halo cells, or do a halo update of 'stagvar' upon return.

  end subroutine glissade_stagger

!----------------------------------------------------------------------------

  subroutine glissade_stagger_real_mask(nx,           ny,        &
                                        var,          stagvar,   &
                                        stagger_rmask)

    !----------------------------------------------------------------
    ! Given a variable on the unstaggered grid (dimension nx, ny), interpolate
    !  to find values on the staggered grid (dimension nx-1, ny-1).
    ! The value in each cell is weighted by a real mask, typically in range [0,1]
    !----------------------------------------------------------------
    !TODO - Combine the two stagger routines using optional arguments.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) ::    &
       var                      ! unstaggered field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    real(dp), dimension(nx,ny), intent(in), optional ::        &
       stagger_rmask            ! real-valued mask that determines how to weight cells in the staggering
                                ! 1 => full weight, 0 => no weight, 0 < mask > 1 => intermediate weight

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp), dimension(nx,ny) :: rmask   ! real-valued mask for weighting
    real(dp) :: sumvar, summask
    integer :: stagger_margin

    if (present(stagger_rmask)) then
       rmask = stagger_rmask
    else
       rmask = 1.0d0   ! default is full weighting of all cells
    endif

    stagvar(:,:) = 0.d0

    do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          sumvar = rmask(i,j+1)*var(i,j+1) + rmask(i+1,j+1)*var(i+1,j+1)  &
                 + rmask(i,j)  *var(i,j)   + rmask(i+1,j)  *var(i+1,j)
          summask = rmask(i,j+1) + rmask(i+1,j+1) + rmask(i,j) + rmask(i+1,j)
          if (summask > 0.d0) stagvar(i,j) = sumvar / summask
       enddo
    enddo

    ! Note: In most cases it would be safe here to do a staggered parallel halo update of stagvar.
    !       However, if running a problem with nonzero periodic_offset_ew/ns (e.g., the ISMIP-HOM and stream cases),
    !        a halo update would copy incorrect values to vertices on the west and south boundaries
    !        of the global domain. So there is no halo update here.
    !       To ensure correct halo values, the user should either pass in a 'var' field that has already
    !        been updated in halo cells, or do a halo update of 'stagvar' upon return.

  end subroutine glissade_stagger_real_mask

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

!----------------------------------------------------------------------------

  subroutine glissade_stagger_edge(nx,           ny,        &
                                   var,                     &
                                   var_east,     var_north, &
                                   stagger_mask)

    !----------------------------------------------------------------
    ! Given a variable on the unstaggered grid (dimension nx, ny), interpolate
    !  to find values at east and north cell edges.
    ! If a mask is present, then cells with mask = 0 are ignored in the average.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
         nx, ny                 ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) ::    &
         var                    ! unstaggered field, defined at cell centers

    real(dp), dimension(nx,ny), intent(out) ::    &
         var_east,            & ! field interpolated to east edges
         var_north              ! field interpolated to north edges

    integer, dimension(nx,ny), intent(in), optional ::        &
       stagger_mask            ! integer mask that determines how to weight cells in the interpolation
                               ! 1 => full weight, 0 => no weight

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    integer, dimension(nx,ny) :: mask   ! mask for weighting
    real(dp) :: sumvar, summask

    if (present(stagger_mask)) then
       mask = stagger_mask
    else
       mask = 1   ! default is full weighting of all cells
    endif

    ! initialize
    var_east(:,:) = 0.0d0
    var_north(:,:) = 0.0d0

    ! east edges
    do j = 1, ny
       do i = 1, nx-1
          sumvar = mask(i,j)*var(i,j) + mask(i+1,j)*var(i+1,j)
          summask = mask(i,j) + mask(i+1,j)
          if (summask > 0.0d0) var_east(i,j) = sumvar / summask
       enddo
    enddo

    ! north edges
    do j = 1, ny-1
       do i = 1, nx
          sumvar = mask(i,j)*var(i,j) + mask(i,j+1)*var(i,j+1)
          summask = mask(i,j) + mask(i,j+1)
          if (summask > 0.0d0) var_north(i,j) = sumvar / summask
       enddo
    enddo

    call parallel_halo(var_east)
    call parallel_halo(var_north)

  end subroutine glissade_stagger_edge

!****************************************************************************

  subroutine glissade_gradient(nx,           ny,        &
                               dx,           dy,        &
                               field,                   &
                               df_dx,        df_dy,     &
                               ice_mask,                &
                               gradient_margin_in)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient is evaluated at the four neighboring vertices and is second-order accurate.
    !
    ! The gradient at a given vertex is taken as the average of the edge gradients
    !  on either side of the vertex.
    ! If gradient_margin_in = 0, then gradients are computed at all edges, even if
    !  one or both cells is ice-free (cf. stagger_margin_in above).
    ! If gradient_margin_in = 1, then gradients are computed only for edges with
    !   ice-covered cells (ice_mask = 1) in each adjacent cells. Other edges are ignored.
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

    integer, dimension(nx,ny), intent(in), optional ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: Compute edge gradient when either cell is ice-covered
                                ! 1: Compute edge gradient only when both cells have ice

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: gradient_margin

    integer :: i, j

    logical, dimension(nx,ny) ::    &
       edge_mask_x,            & ! edge mask for computing df/dx
       edge_mask_y               ! edge mask for computing df/dy

    real(dp) :: df_dx_north, df_dx_south  ! df_dx at neighboring edges
    real(dp) :: df_dy_east, df_dy_west    ! df_dx at neighboring edges

    ! Initialize

    df_dx(:,:) = 0.0d0
    df_dy(:,:) = 0.0d0

    edge_mask_x(:,:) = .false.
    edge_mask_y(:,:) = .false.

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = 0  ! default is to average over all cells, including those where ice is absent
    endif

    if (gradient_margin == 1 .and. .not.present(ice_mask)) then
       call write_log('Must pass in ice_mask to compute gradient with gradient_margin = 1', GM_FATAL)
    endif

    !--------------------------------------------------------
    !   Gradient at vertex(i,j) is based on f(i:i+1,j:j+1)
    !
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)
    !--------------------------------------------------------

    ! Create masks identifying edges that will be used in gradient computations

    if (gradient_margin == 0) then

       edge_mask_x(:,:) = .true.       ! true for all edges
       edge_mask_y(:,:) = .true.

    elseif (gradient_margin == 1) then

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

    ! Compute gradient at vertices by averaging gradient at adjacent edges.
    ! Ignore edges with edge_mask = F

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

  end subroutine glissade_gradient

!****************************************************************************

  subroutine glissade_gradient_at_edges(nx,           ny,        &
                                        dx,           dy,        &
                                        field,                   &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        gradient_margin_in,      &
                                        usrf,                    &
                                        land_mask,               &
                                        max_slope)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) at cell edges (i.e., the C grid):
    ! df_dx at the midpoint of the east edge and df_dy at the midpoint of
    ! the north edge.
    !
    ! This subroutine is called by the glissade SIA solver.
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
       land_mask                ! = 1 for land cells, else = 0
                                ! floating and land masks required for gradient_margin = HO_GRADIENT_MARGIN_HYBRID

    real(dp), intent(in), optional :: &
       max_slope                ! maximum slope allowed for surface gradient computations (unitless)

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: gradient_margin
    integer :: i, j

    logical, dimension(nx,ny) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny) ::    &
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
       gradient_margin = HO_GRADIENT_MARGIN_HYBRID
    endif

    ! Set logical edge mask based on gradient_margin.

    edge_mask_x(:,:) = .false.
    edge_mask_y(:,:) = .false.

    if (gradient_margin == HO_GRADIENT_MARGIN_LAND) then

       edge_mask_x(:,:) = .true.       ! true for all edges
       edge_mask_y(:,:) = .true.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_HYBRID) then

       if (present(land_mask) .and. present(usrf)) then

          call glissade_edgemask_gradient_margin_hybrid(nx,          ny,         &
                                                        ice_mask,                &
                                                        land_mask,               &
                                                        usrf,                    &
                                                        edge_mask_x, edge_mask_y)
       else
          call write_log('Must pass in land_mask and usrf to use this gradient_margin option', GM_FATAL)
       endif

    elseif (gradient_margin == HO_GRADIENT_MARGIN_MARINE) then

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

  subroutine glissade_edgemask_gradient_margin_hybrid(nx,          ny,         &
                                                      ice_mask,                &
                                                      land_mask,               &
                                                      usrf,                    &
                                                      edge_mask_x, edge_mask_y)
    
    !----------------------------------------------------------------
    ! Compute edge masks required for option gradient_margin = HO_GRADIENT_MARGIN_HYBRID
    ! Called from subroutine glissade_gradient_at_edges.
    !
    ! The mask is set to true at all edges where either
    ! (1) Both adjacent cells are ice-covered.
    ! (2) One cell is ice-covered, and lies above ice-free land.
    !
    ! This method sets the gradient to zero at edges where
    ! (1) An ice-covered cell (grounded or floating) lies above ice-free ocean.
    !     Note: Inactive calving-front cells are treated as ice-free ocean.
    ! (2) An ice-covered land cell lies below an ice-free land cell (i.e., a nunatak).
    !
    ! This method aims to give a reasonable gradient at both land-terminating and marine-terminating margins.
    ! At land-terminating margins the gradient is nonzero (except for nunataks), and at marine-terminating
    !  margins the gradient is zero.
    !
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions
   
    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask,              & ! = 1 where ice is present, else = 0
       land_mask                ! = 1 for land cells, else = 0

    real(dp), dimension(nx,ny), intent(in)  ::       &
       usrf                     ! ice surface elevation

    logical, dimension(nx,ny), intent(out) ::    &
       edge_mask_x,           & ! edge mask for computing df/dx
       edge_mask_y              ! edge mask for computing df/dy

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    ! initialize
    edge_mask_x(:,:) = .false.
    edge_mask_y(:,:) = .false.

    ! compute mask for east and west cell edges

    do j = 1, ny
       do i = 1, nx-1

          if (ice_mask(i,j) == 1 .and. ice_mask(i+1,j) == 1) then  ! both cells have ice

             edge_mask_x(i,j) = .true.

          elseif ( (ice_mask(i,j)==1   .and.  land_mask(i+1,j)==1 .and. usrf(i,j)   > usrf(i+1,j)) .or.  &
                   (ice_mask(i+1,j)==1 .and.  land_mask(i,j)==1   .and. usrf(i+1,j) > usrf(i,j)) ) then

             ! ice-covered cell lies above ice-free land
             edge_mask_x(i,j) = .true.

          endif

       enddo
    enddo
    
    ! compute mask for north and south cell edges

    do j = 1, ny-1
       do i = 1, nx

          if (ice_mask(i,j) == 1 .and. ice_mask(i,j+1) == 1) then  ! both cells have ice

             edge_mask_y(i,j) = .true.

          elseif ( (ice_mask(i,j)==1   .and.  land_mask(i,j+1)==1 .and. usrf(i,j)   > usrf(i,j+1)) .or.  &
                   (ice_mask(i,j+1)==1 .and.  land_mask(i,j)==1   .and. usrf(i,j+1) > usrf(i,j)) ) then

             ! ice-covered cell lies above ice-free land
             edge_mask_y(i,j) = .true.

          endif

       enddo
    enddo

  end subroutine glissade_edgemask_gradient_margin_hybrid

!****************************************************************************

  subroutine glissade_surface_elevation_gradient_edges(&
       nx,           ny,          &
       dx,           dy,          &
       itest, jtest, rtest,       &
       active_ice_mask,           &
       land_mask,                 &
       usrf,         thck,        &
       topg,         eus,         &
       thklim,                    &
       thck_gradient_ramp,        &
       ds_dx_east,   ds_dy_east,  &
       ds_dx_north,  ds_dy_north, &
       ho_gradient_margin,        &
       max_slope)

    !----------------------------------------------------------------
    ! Compute surface elevation gradients at east and north cell edges
    !  for different ho_gradient_margin options.
    !
    ! The ds/dx gradient at east edges and the ds_dy gradients at north edges
    !  are computed in the standard way, by taking the difference of the values
    !  in two adjacent cells and dividing by the distance.
    ! The dx_dx gradient at north edges is found by averaging ds_dx from nearby east edges,
    !  and ds_dy at east edges is found by averaging ds_dy from nearby north edges.
    !
    ! At the ice margin, where one or both cells adjacent to a given edge may be ice-free,
    !  edge gradients may be masked in the following ways:
    !
    ! HO_GRADIENT_MARGIN_LAND = 0: Values in both adjacent cells are used to compute the gradient,
    !  including values in ice-free cells.  In other words, there is no masking of edges.
    !  This convention is used by Glide. It works well at land-terminating margins, but performs poorly
    !  for ice shelves with a sharp drop in ice thickness and surface elevation at the margin.
    !
    ! HO_GRADIENT_MARGIN_HYBRID = 1: The gradient is computed at edges where either
    ! (1) Both adjacent cells are ice-covered.
    ! (2) One cell is ice-covered (land or marine-based) and lies above ice-free land.
    !
    ! This method sets the gradient to zero at edges where
    ! (1) An ice-covered cell (grounded or floating) lies above ice-free ocean.
    !     Note: Inactive calving-front cells are treated as ice-free ocean.
    ! (2) An ice-covered land cell lies below an ice-free land cell (i.e., a nunatak).

    ! The aim is to give a reasonable gradient at both land-terminating and marine-terminating margins.
    ! At land-terminating margins the gradient is nonzero (except for nunataks), and at marine-terminating
    !  margins the gradient is zero.
    !
    ! HO_GRADIENT_MARGIN_MARINE = 2: Only values in ice-covered cells (i.e., cells with thck > thklim)
    !  are used to compute gradients.  If one or both adjacent cells is ice-free, the edge is masked out.
    !  This option works well at shelf margins but less well for land margins (e.g., the Halfar test case).
    !
    ! The HO_GRADIENT option is not used.  All gradients are centered; none are upstream-biased.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
         nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
         dx, dy                   ! horizontal grid size

    integer, intent(in) ::      &
         itest, jtest, rtest      ! coordinates of diagnostic point

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

    real(dp), dimension(nx,ny), intent(out) ::    &
         ds_dx_east,  ds_dy_east,  & ! ds/dx at east edges
         ds_dx_north, ds_dy_north    ! ds/dx at north edges

    integer, intent(in) ::      &
         ho_gradient_margin       ! option for computing gradients at ice sheet margin
                                  ! see comments above

    real(dp), intent(in), optional ::       &
         max_slope                ! maximum slope allowed for surface gradient computations (unitless)

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j
    integer :: iu, il, ju, jl     ! indices of cells with higher and lower elevations

    real(dp) ::  &
         edge_thck_upper,       & ! thickness of higher-elevation cell
         edge_thck_lower,       & ! thickness of lower-elevation cell
         edge_factor,           & ! gradient-weighting factor in range [0,1]
         sign_factor              ! sign factor, +1 or -1

    ! initialize

    ds_dx_east(:,:) = 0.0d0
    ds_dy_east(:,:) = 0.0d0

    ds_dx_north(:,:) = 0.0d0
    ds_dy_north(:,:) = 0.0d0

    if (ho_gradient_margin == HO_GRADIENT_MARGIN_LAND) then

       ! Compute ds_dx and ds_dy on all edges, whether or not adjacent cells are ice-covered
       do j = 1, ny
          do i = 1, nx-1
             ds_dx_east(i,j) = (usrf(i+1,j) - usrf(i,j)) / dx
          enddo
       enddo

       do j = 1, ny-1
          do i = 1, nx
             ds_dy_north(i,j) = (usrf(i,j+1) - usrf(i,j)) / dy
          enddo
       enddo

    elseif (ho_gradient_margin == HO_GRADIENT_MARGIN_MARINE) then

       ! Compute ds_dx and ds_dy only on edges with active ice in each adjacent cell
       do j = 1, ny
          do i = 1, nx-1
             if (active_ice_mask(i,j) == 1 .and. active_ice_mask(i+1,j) == 1) then
                ds_dx_east(i,j) = (usrf(i+1,j) - usrf(i,j)) / dx
             endif
          enddo
       enddo

       do j = 1, ny-1
          do i = 1, nx
             if (active_ice_mask(i,j) == 1 .and. active_ice_mask(i,j+1) == 1) then
                ds_dy_east(i,j) = (usrf(i,j+1) - usrf(i,j)) / dy
             endif
          enddo
       enddo

    elseif (ho_gradient_margin == HO_GRADIENT_MARGIN_HYBRID) then

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

                ! compute the gradient
                ds_dx_east(i,j) = edge_factor * sign_factor * (usrf(iu,j) - usrf(il,j)) / dx

             elseif (active_ice_mask(iu,j) == 1 .and. land_mask(il,j) == 1) then

                ! upper cell has active ice, and ice-free lower cell is land; compute the gradient
                ds_dx_east(i,j) = edge_factor * sign_factor * (usrf(iu,j) - usrf(il,j)) / dx

             endif

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

                ! compute the gradient
                ds_dy_north(i,j) = edge_factor * sign_factor * (usrf(i,ju) - usrf(i,jl)) / dy

             elseif (active_ice_mask(i,ju) == 1 .and. land_mask(i,jl) == 1) then

                ! upper cell has active ice, and ice-free lower cell is land; compute the gradient
                ds_dy_north(i,j) = edge_factor * sign_factor * (usrf(i,ju) - usrf(i,jl)) / dy

             endif  ! both cells have ice

          enddo  ! i
       enddo   ! j

    endif   ! ho_gradient_margin

    ! halo updates
    call parallel_halo(ds_dx_east)
    call parallel_halo(ds_dy_north)

    ! Optionally, limit ds_dx_east and ds_dy_north

    if (present(max_slope)) then

       do j = 1, ny
          do i = 1, nx
             if (ds_dx_east(i,j) > 0.0d0) then
                ds_dx_east(i,j) = min(ds_dx_east(i,j), max_slope)
             else
                ds_dx_east(i,j) = max(ds_dx_east(i,j), -max_slope)
             endif
          enddo
       enddo

       do j = 1, ny
          do i = 1, nx
             if (ds_dy_north(i,j) > 0.0d0) then
                ds_dy_north(i,j) = min(ds_dy_north(i,j), max_slope)
             else
                ds_dy_north(i,j) = max(ds_dy_north(i,j), -max_slope)
             endif
          enddo
       enddo

    endif   ! present(max_slope)

    ! Average ds_dx to north edges
    do j = 1, ny-1
       do i = 2, nx
          ds_dx_north(i,j) = 0.25d0 * &
               (ds_dx_east(i-1,j) + ds_dx_east(i,j) + ds_dx_east(i-1,j+1) + ds_dx_east(i,j+1))
       enddo
    enddo

    ! Average ds_dy to east edges
    do j = 2, ny
       do i = 1, nx-1
          ds_dy_east(i,j) = 0.25d0 * &
               (ds_dy_north(i,j-1) + ds_dy_north(i,j) + ds_dy_north(i+1,j-1) + ds_dy_north(i+1,j))
       enddo
    enddo


    ! halo updates
    call parallel_halo(ds_dx_north)
    call parallel_halo(ds_dy_east)


    if (verbose_gradient .and. this_rank==rtest) then
       print*, ' '
       print*, 'Gradients, i, j, task =', itest, jtest, rtest
       print*, ' '
       print*, 'ds_dx_east:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f9.6)',advance='no') ds_dx_east(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'ds_dx_north:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f9.6)',advance='no') ds_dx_north(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'ds_dy_north:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f9.6)',advance='no') ds_dy_north(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'ds_dy_east:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f9.6)',advance='no') ds_dy_east(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_surface_elevation_gradient_edges

!****************************************************************************

  subroutine glissade_surface_elevation_gradient(nx,           ny,        &
                                                 dx,           dy,        &
                                                 itest, jtest, rtest,     &
                                                 active_ice_mask,         &
                                                 land_mask,               &
                                                 usrf,         thck,      &
                                                 topg,         eus,       &
                                                 thklim,                  &
                                                 thck_gradient_ramp,      &
                                                 ds_dx,        ds_dy,     &
                                                 ho_gradient,             &
                                                 ho_gradient_margin,      &
                                                 max_slope)

    !----------------------------------------------------------------
    ! Compute surface elevation gradients for different ho_gradient and ho_gradient_margin options.
    !
    ! The gradient at a given vertex is constructed from gradients at adjacent edges.
    ! Edge gradients are computed in the standard way, taking the difference between
    !  the values in two adjacent cells and dividing by the distance.
    ! At the ice margin, where one or both cells adjacent to a given edge may be ice-free,
    !  edge gradients may be masked in the following ways:
    !
    ! HO_GRADIENT_MARGIN_LAND = 0: Values in both adjacent cells are used to compute the gradient,
    !  including values in ice-free cells.  In other words, there is no masking of edges.
    !  This convention is used by Glide. It works well at land-terminating margins, but performs poorly
    !  for ice shelves with a sharp drop in ice thickness and surface elevation at the margin.
    !
    ! HO_GRADIENT_MARGIN_HYBRID = 1: The gradient is computed at edges where either
    !
    ! (1) Both adjacent cells are ice-covered.
    ! (2) One cell is ice-covered (land or marine-based) and lies above ice-free land.
    !
    ! This method sets the gradient to zero at edges where
    ! (1) An ice-covered cell (grounded or floating) lies above ice-free ocean.
    !     Note: Inactive calving-front cells are treated as ice-free ocean.
    ! (2) An ice-covered land cell lies below an ice-free land cell (i.e., a nunatak).

    ! The aim is to give a reasonable gradient at both land-terminating and marine-terminating margins.
    ! At land-terminating margins the gradient is nonzero (except for nunataks), and at marine-terminating
    !  margins the gradient is zero.
    !
    ! HO_GRADIENT_MARGIN_MARINE = 2: Only values in ice-covered cells (i.e., cells with thck > thklim)
    !  are used to compute gradients.  If one or both adjacent cells is ice-free, the edge is masked out.
    !  This option works well at shelf margins but less well for land margins (e.g., the Halfar test case).
    !
    ! There are three ways to compute vertex gradients from edge gradients, as determined
    !  by the ho_gradient option:
    !
    ! HO_GRADIENT = 0: Standard centered gradient, obtained by averaging the two nearest edge gradients
    !
    ! HO_GRADIENT = 1: First-order upstream gradient, obtained by setting the vertex value to the
    !                  nearest edge gradient on the higher-elevation side of the vertex
    !
    ! HO_GRADIENT = 2: Second-order upstream gradient, obtained by setting the vertex value to a linear
    !                  combination of the two nearest edge gradients on the higher-elevation side
    !
    ! The centered gradient is usually most accurate, but it can lead to checkerboard noise
    !  in the surface elevation field, because a checkerboard pattern is invisible to the gradient.
    ! The upstream gradients are less accurate but are better at damping checkerboard noise.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
         nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
         dx, dy                   ! horizontal grid size

    integer, intent(in) ::      &
         itest, jtest, rtest      ! coordinates of diagnostic point

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
         ds_dx,                 & ! ds/dx at vertices
         ds_dy                    ! ds/dy at vertices

    integer, intent(in) ::      &
         ho_gradient,           & ! gradient type (centered, 1st-order upstream, or 2nd-order upstream)
         ho_gradient_margin       ! option for computing gradients at ice sheet margin
                                  ! see comments above

    real(dp), intent(in), optional ::       &
         max_slope                ! maximum slope allowed for surface gradient computations (unitless)

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

    real(dp) :: sum1, sum2        ! temporary sums

    ! initialize

    ds_dx_edge(:,:) = 0.0d0
    ds_dy_edge(:,:) = 0.0d0

    ds_dx(:,:) = 0.0d0
    ds_dy(:,:) = 0.0d0

    if (ho_gradient_margin == HO_GRADIENT_MARGIN_LAND) then

       ! Compute ds_dx and ds_dy on all edges, whether or not adjacent cells are ice-covered
       do j = 1, ny
          do i = 1, nx-1
             ds_dx_edge(i,j) = (usrf(i+1,j) - usrf(i,j)) / dx
          enddo
       enddo

       do j = 1, ny-1
          do i = 1, nx
             ds_dy_edge(i,j) = (usrf(i,j+1) - usrf(i,j)) / dy
          enddo
       enddo

    elseif (ho_gradient_margin == HO_GRADIENT_MARGIN_MARINE) then

       ! Compute ds_dx and ds_dy only on edges with active ice in each adjacent cell
       do j = 1, ny
          do i = 1, nx-1
             if (active_ice_mask(i,j) == 1 .and. active_ice_mask(i+1,j) == 1) then
                ds_dx_edge(i,j) = (usrf(i+1,j) - usrf(i,j)) / dx
             endif
          enddo
       enddo

       do j = 1, ny-1
          do i = 1, nx
             if (active_ice_mask(i,j) == 1 .and. active_ice_mask(i,j+1) == 1) then
                ds_dy_edge(i,j) = (usrf(i,j+1) - usrf(i,j)) / dy
             endif
          enddo
       enddo

    elseif (ho_gradient_margin == HO_GRADIENT_MARGIN_HYBRID) then

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

                ! compute the gradient
                ds_dx_edge(i,j) = edge_factor * sign_factor * (usrf(iu,j) - usrf(il,j)) / dx

             elseif (active_ice_mask(iu,j) == 1 .and. land_mask(il,j) == 1) then

                ! upper cell has active ice, and ice-free lower cell is land; compute the gradient
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

                ! compute the gradient
                ds_dy_edge(i,j) = edge_factor * sign_factor * (usrf(i,ju) - usrf(i,jl)) / dy

             elseif (active_ice_mask(i,ju) == 1 .and. land_mask(i,jl) == 1) then

                ! upper cell has active ice, and ice-free lower cell is land; compute the gradient
                ds_dy_edge(i,j) = edge_factor * sign_factor * (usrf(i,ju) - usrf(i,jl)) / dy

             endif  ! both cells have ice

          enddo  ! i
       enddo   ! j

    endif   ! ho_gradient_margin


    ! Average the edge gradients to the vertex, depending on the value of ho_gradient.

    if (ho_gradient == HO_GRADIENT_CENTERED) then

       ! Average the edge gradients to vertices using a centered approximation.
       ! This method is 2nd order accurate but can be subject to B grid noise.
       do j = 1, ny-1
          do i = 1, nx-1
             ds_dx(i,j) = 0.5d0 * (ds_dx_edge(i,j) + ds_dx_edge(i,j+1))
             ds_dy(i,j) = 0.5d0 * (ds_dy_edge(i,j) + ds_dy_edge(i+1,j))
          enddo
       enddo

       !WHL - This is an alternate method of averaging edge gradients to vertices, following glissade_centered_gradient.
       !       It is commented out.
       !      For the dome problem, the differences between gradient_margin methods 1 and 3
       !       come from replacing these lines with the ds_dx and ds_dy calculations above.

!    do j = 1, ny-1
!       do i = 1, nx-1
!          if (ds_dx_edge(i,j) /= 0.0d0 .and. ds_dx_edge(i,j+1) /= 0.0d0) then
!             ds_dx(i,j) = 0.5d0 * (ds_dx_edge(i,j) + ds_dx_edge(i,j+1))
!          elseif (ds_dx_edge(i,j) /= 0.0d0) then
!             ds_dx(i,j) = ds_dx_edge(i,j)
!          elseif (ds_dx_edge(i,j+1) /= 0.0d0) then
!             ds_dx(i,j) = ds_dx_edge(i,j+1)
!          else
!             ds_dx(i,j) = 0.d0
!          endif
!       enddo
!    enddo

!    do j = 1, ny-1
!       do i = 1, nx-1
!          if (ds_dy_edge(i,j) /= 0.0d0 .and. ds_dy_edge(i+1,j) /= 0.0d0) then
!             ds_dy(i,j) = 0.5d0 * (ds_dy_edge(i,j) + ds_dy_edge(i+1,j))
!          elseif (ds_dy_edge(i,j) /= 0.0d0) then
!             ds_dy(i,j) = ds_dy_edge(i,j)
!          elseif (ds_dy_edge(i+1,j) /= 0.0d0) then
!             ds_dy(i,j) = ds_dy_edge(i+1,j)
!          else
!             ds_dy(i,j) = 0.d0
!          endif
!       enddo
!    enddo


    elseif (ho_gradient == HO_GRADIENT_UPSTREAM1) then

       ! Take a one-sided, first-order-accurate upstream gradient.
       ! For ds_dx, use the east/west edge that is upstream in the y direction
       ! For ds_dy, use the north/south edge that is upstream in the x direction

       do j = 1, ny-1
          do i = 1, nx-1

             ! Identify the upstream edge in the y direction
             sum1 = usrf(i,j+1) + usrf(i+1,j+1)
             sum2 = usrf(i,j) + usrf(i+1,j)

             if (sum1 > sum2) then  ! north is upstream; use east edge of cell (i,j+1)
                ds_dx(i,j) = ds_dx_edge(i,j+1)
             else                   ! south is upstream; use east edge of cell (i,j)
                ds_dx(i,j) = ds_dx_edge(i,j)
             endif

             ! Identify the upstream edge in the x direction
             sum1 = usrf(i+1,j) + usrf(i+1,j+1)
             sum2 = usrf(i,j) + usrf(i,j+1)

             if (sum1 > sum2) then  ! east is upstream; use north edge of cell (i+1,j)
                ds_dy(i,j) = ds_dy_edge(i+1,j)
             else                   ! west is upstream; use north edge of cell (i,j)
                ds_dy(i,j) = ds_dy_edge(i,j)
             endif

          enddo
       enddo

    elseif (ho_gradient == HO_GRADIENT_UPSTREAM2) then

       ! Take a one-sided, second-order-accurate upstream gradient.

       do j = 2, ny-2
          do i = 2, nx-2

             ! Identify the upstream edge in the y direction
             sum1 = usrf(i,j+1) + usrf(i+1,j+1) + usrf(i,j+2) + usrf(i+1,j+2)
             sum2 = usrf(i,j) + usrf(i+1,j) + usrf(i,j-1) + usrf(i+1,j-1)

             ! Compute df_dx by taking 2nd-order upstream gradient

             if (sum1 > sum2) then   ! north is upstream; use east edge of cells (i,j+1:j+2)
                ds_dx(i,j) = 1.5d0 * ds_dx_edge(i,j+1) - 0.5d0 * ds_dx_edge(i,j+2)
             else                    ! south is upstream; use east edge of cells (i,j-1:j)
                ds_dx(i,j) = 1.5d0 * ds_dx_edge(i,j) - 0.5d0 * ds_dx_edge(i,j-1)
             endif

             ! Identify the upstream edge in the x direction
             sum1 = usrf(i+1,j) + usrf(i+1,j+1) + usrf(i+2,j) + usrf(i+2,j+1)
             sum2 = usrf(i,j) + usrf(i,j+1) + usrf(i-1,j) + usrf(i-1,j+1)

             ! Compute df_dx by taking 2nd-order upstream gradient

             if (sum1 > sum2) then   ! east is upstream; use north edge of cells (i+1:i+2,j)
                ds_dy(i,j) = 1.5d0 * ds_dy_edge(i+1,j) - 0.5d0 * ds_dy_edge(i+2,j)
             else                    ! west is upstream; use north edge of cells (i-1:i,j)
                ds_dy(i,j) = 1.5d0 * ds_dy_edge(i,j) - 0.5d0 * ds_dy_edge(i-1,j)
             endif

          enddo
       enddo

    endif   ! ho_gradient


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

    ! halo update
    call staggered_parallel_halo(ds_dx)
    call staggered_parallel_halo(ds_dy)

    if (verbose_gradient .and. this_rank==rtest) then
       print*, ' '
       print*, 'Hybrid gradient, i, j, task =', itest, jtest, rtest
       print*, ' '
       print*, 'ds_dx_edge:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f9.6)',advance='no') ds_dx_edge(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'ds_dy_edge:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f9.6)',advance='no') ds_dy_edge(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'ds_dx:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f9.6)',advance='no') ds_dx(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'ds_dy:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(f9.6)',advance='no') ds_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_surface_elevation_gradient

!****************************************************************************

  subroutine glissade_slope_angle(nx,         ny,      &
                                  dx,         dy,      &
                                  zsfc,                &
                                  theta_slope,         &
                                  slope_mask_in)

    ! Compute the slope angle between a surface of elevation zsfc and the horizontal.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
         nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
         dx, dy                   ! grid cell length and width
                                  ! assumed to have the same value for each grid cell

    real(dp), dimension(nx,ny), intent(in) ::       &
         zsfc                     ! elevation of the surface whose slope is to be computed

    real(dp), dimension(nx,ny), intent(out) ::       &
         theta_slope              ! angle formed by the surface with the horizontal (x-y) direction

    integer, dimension(nx,ny), intent(in), optional ::  &
         slope_mask_in            ! = 1 for the part of the surface whose slope is computed, else = 0

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer, dimension(nx,ny) ::  &
         slope_mask                 ! local version of slope_mask

    real(dp), dimension(nx-1,ny-1) ::    &
         dz_dx, dz_dy,         &  ! gradient components of zsfc, defined at cell vertices
         stag_slope               ! slope (= magnitude of gradient) at cell vertices

    real(dp), dimension(nx,ny) ::    &
         slope                    ! stag_slope interpolated to cell centers

    ! initialize

    if (present(slope_mask_in)) then
       slope_mask = slope_mask_in
    else
       slope_mask = 1    ! default is to compute the slope everywhere
    endif

    theta_slope = 0.0d0

    ! Compute the x and y components of the surface gradient
    ! Note: With gradient_margin_in = 1, edge gradients are computed only for edges
    !        with slope_mask = 1 on either side.
    !       For instance, if slope_mask = 1 for floating cells only, then dz_dz = dz_dy = 0
    !        for grounded regions.

    call glissade_gradient(nx,         ny,          &
                           dx,         dy,          &
                           zsfc,                    &
                           dz_dx,      dz_dy,       &
                           slope_mask,              &
                           gradient_margin_in = 1)

    ! Compute the magnitude of the gradient.  This is the scalar slope, on the staggered grid.
    stag_slope = sqrt(dz_dx**2 + dz_dy**2)

    ! Interpolate the slope to cell centers

    call glissade_unstagger(nx,          ny,        &
                            stag_slope,  slope)

    ! Compute the slope angle
    theta_slope = atan(slope)

  end subroutine glissade_slope_angle

!****************************************************************************

  subroutine glissade_laplacian_smoother(nx,         ny,            &
                                         var,        var_smooth,    &
                                         smoother_mask,             &
                                         npoints_stencil)

    !----------------------------------------------------------------
    ! Given a 2D field on the ice grid, smooth the field using a 9-point Laplacian stencil.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: nx, ny                             ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) :: var             ! input field, before smoothing

    real(dp), dimension(nx,ny), intent(out) :: var_smooth     ! output field, after smoothing

    real(dp), dimension(nx,ny), intent(in), optional :: &
         smoother_mask                                        ! real mask to weight the cells included in the smoothing

    integer, intent(in), optional :: &
         npoints_stencil                                      ! number of points in stencil, either 5 or 9

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp), dimension(nx,ny) :: rmask   ! real mask set to the optional smoother_mask, else = 1 everywhere

    real(dp) :: sum_mask   ! sum of mask values in the cell and its neighbors, converted to real(dp)

    integer :: npoints  ! set equal to npoints_stencil, else defaults to 5

    integer :: i, j

    if (present(smoother_mask)) then
       rmask = smoother_mask
    else
       rmask = 1.0d0
    endif

    if (present(npoints_stencil)) then
       npoints = npoints_stencil
       if (.not.(npoints == 5 .or. npoints == 9)) then
          call write_log('ERROR, glissade_laplacian_smoother: Must choose 5 or 9 points for the stencil', GM_FATAL)
       endif
    else
       npoints = 5
    endif

    sum_mask = 0.0d0

    if (npoints == 5) then

       do j = 2, ny-1
          do i = 2, nx-1

             sum_mask =  4.0d0 *  rmask(i,j)  &
                       + 1.0d0 * (rmask(i-1,j) + rmask(i+1,j) + rmask(i,j-1) + rmask(i,j+1))

             if (sum_mask > 0.0d0) then

                var_smooth(i,j) = (1.0d0/sum_mask) *  &
                                 ( 4.0d0 *  rmask(i,j)*var(i,j)  &
                                 + 1.0d0 * (rmask(i-1,j)*var(i-1,j) + rmask(i+1,j)*var(i+1,j)   &
                                          + rmask(i,j-1)*var(i,j-1) + rmask(i,j+1)*var(i,j+1)) )
             else
                var_smooth(i,j) = var(i,j)
             endif

          enddo
       enddo

    elseif (npoints == 9) then

       call write_log('Apply Laplacian smoother with 9-point stencil')

       do j = 2, ny-1
          do i = 2, nx-1

             sum_mask =  4.0d0 *  rmask(i,j)  &
                       + 2.0d0 * (rmask(i-1,j) + rmask(i+1,j) + rmask(i,j-1) + rmask(i,j+1))  &
                       + 1.0d0 * (rmask(i-1,j+1) + rmask(i+1,j+1) + rmask(i-1,j-1) + rmask(i+1,j-1))

             if (sum_mask > 0.0d0) then

                var_smooth(i,j) = (1.0d0/sum_mask) *  &
                                 ( 4.0d0 *  rmask(i,j)*var(i,j) &
                                 + 2.0d0 * (rmask(i-1,j)*var(i-1,j) + rmask(i+1,j)*var(i+1,j)   &
                                          + rmask(i,j-1)*var(i,j-1) + rmask(i,j+1)*var(i,j+1))  &
                                 + 1.0d0 * (rmask(i-1,j+1)*var(i-1,j+1) + rmask(i+1,j+1)*var(i+1,j+1)   &
                                          + rmask(i-1,j-1)*var(i-1,j-1) + rmask(i+1,j-1)*var(i+1,j-1)) )
             else
                var_smooth(i,j) = var(i,j)
             endif

          enddo
       enddo

    endif   ! npoints

    !Note: Do not do a halo update, because it is unknown whether we should call parallel_halo or staggered_parallel_halo.
!!    call parallel_halo(var_smooth)

  end subroutine glissade_laplacian_smoother

!****************************************************************************

  subroutine glissade_vertical_average(nx,         ny,        &
                                       nz,         sigma,     &
                                       var,        var_2d,    &
                                       mask)

    !----------------------------------------------------------------
    ! Compute the vertical average of a given variable.
    ! Note: It is assumed that the variable is defined at layer midpoints,
    !       and hence has vertical dimension (nz-1).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,            &     ! horizontal grid dimensions
       nz                       ! number of vertical levels

    real(dp), dimension(nz), intent(in) ::    &
       sigma                    ! sigma vertical coordinate

    real(dp), dimension(nz-1,nx, ny), intent(in) ::    &
       var                      ! 3D field to be averaged vertically

    real(dp), dimension(nx, ny), intent(out) ::    &
       var_2d                   ! 2D vertically averaged field

    logical, dimension(nx, ny), intent(in), optional ::    &
       mask                     ! compute var_2d where mask = .true.


    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j, k

    if (present(mask)) then

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

    else

       do j = 1, ny
          do i = 1, nx
             var_2d(i,j) = 0.d0
             do k = 1, nz-1
                var_2d(i,j) = var_2d(i,j) + var(k,i,j) * (sigma(k+1) - sigma(k))
             enddo
          enddo
       enddo

    endif   ! present(mask)

  end subroutine glissade_vertical_average

!****************************************************************************

  subroutine glissade_vertical_interpolate(nx,       ny,           &
                                           nz,       vert_levels,  &
                                           usrf,                   &
                                           field_3d, field,        &
                                           linear_extrapolate_in)


    ! Given an input field (e.g., surface mass balance), supplied at various vertical levels,
    ! downscale the field to elevation usrf by linear interpolation.

    ! input-output arguments

    integer, intent(in) ::   &
       nx, ny,            &     ! horizontal grid dimensions
       nz                       ! number of vertical levels

    real(dp), dimension(nz), intent(in) :: vert_levels       ! vertical levels at which input data is supplied
    real(dp), dimension(nx,ny), intent(in) :: usrf           ! upper surface elevation
    real(dp), dimension(nz,nx,ny), intent(in) :: field_3d    ! input field supplied at vertical levels
    real(dp), dimension(nx,ny), intent(out) :: field         ! output field downscaled to elevation usrf

    logical, intent(in), optional ::  &
         linear_extrapolate_in    ! if true, then use linear extrapolation outside the given range of vertical levels
                                  ! if false, then simply extrapolate the top and bottom values

    ! local variables

    integer :: i, j, k
    real(dp) :: field_grad        ! vertical gradient of field
    logical :: linear_extrapolate

    if (present(linear_extrapolate_in)) then
       linear_extrapolate = linear_extrapolate_in
    else
       linear_extrapolate = .false.
    endif

    do j = 1, ny
       do i = 1, nx

          if (usrf(i,j) >= vert_levels(nz)) then

             if (linear_extrapolate) then  ! linear extrapolation from top two levels
                field_grad = (field_3d(nz,i,j) - field_3d(nz-1,i,j)) / (vert_levels(nz) - vert_levels(nz-1))
                field(i,j) = field_3d(nz,i,j) + field_grad * (usrf(i,j) - vert_levels(nz))
             else  ! simply extend the top value
                field(i,j) = field_3d(nz,i,j)
             endif

          elseif (usrf(i,j) < vert_levels(1)) then

             if (linear_extrapolate) then  ! linear extrapolation from bottom two levels
                field_grad = (field_3d(2,i,j) - field_3d(1,i,j)) / (vert_levels(2) - vert_levels(1))
                field(i,j) = field_3d(1,i,j) + field_grad * (usrf(i,j) - vert_levels(1))  ! note usrf - vert_levels(1) < 0
             else  ! simply extend the bottom value
                field(i,j) = field_3d(1,i,j)
             endif

          else  ! vert_levels(1) <= usrf <= vert_levels(nz)
                ! linear interpolation between adjacent levels

             do k = 2, nz

                if (usrf(i,j) < vert_levels(k)) then
                   field_grad = (field_3d(k,i,j) - field_3d(k-1,i,j)) / (vert_levels(k) - vert_levels(k-1))
                   field(i,j) = field_3d(k-1,i,j) + field_grad * (usrf(i,j) - vert_levels(k-1))
                endif

             enddo  ! k

          endif  ! usrf v. smb_levels

       enddo   ! i
    enddo   ! j

    call parallel_halo(field)

  end subroutine glissade_vertical_interpolate

!****************************************************************************

  end module glissade_grid_operators

!****************************************************************************
