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
    use glimmer_physcon, only: pi
    use glimmer_log
    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo, &
         parallel_type, parallel_halo, parallel_reduce_sum, parallel_globalindex
  ! Note: Using the glide_diagnostics module creates a circularity.
  !       For that reason, point_diag should be called from a higher level, not from inside this module.
!!    use glide_diagnostics, only: point_diag

    implicit none

    private
    public :: glissade_stagger, glissade_unstagger, glissade_stagger_real_mask, &
              glissade_gradient, glissade_gradient_at_edges, &
              glissade_average_to_edges,            &
              glissade_surface_elevation_gradient,  &
              glissade_laplacian, glissade_laplacian_stagvar, &
              glissade_laplacian_smoother,          &
              glissade_slope_angle, glissade_slope_angle_staggered, &
              glissade_vertical_average,            &
              glissade_vertical_interpolate,        &
              glissade_scalar_extrapolate

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
       stagger_margin = 0  ! default is to average over all vertices, including those where ice is absent
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

    ! Note: Halo update moved to higher level
!    call parallel_halo(unstagvar)

  end subroutine glissade_unstagger

!****************************************************************************

  subroutine glissade_gradient(nx,           ny,        &
                               dx,           dy,        &
                               itest, jtest, rtest,     &
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

    integer, intent(in) ::      &
       itest, jtest, rtest      ! coordinates of diagnostic point

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

  end subroutine glissade_gradient

!****************************************************************************

  subroutine glissade_gradient_at_edges(nx,           ny,        &
                                        dx,           dy,        &
                                        itest, jtest, rtest,     &
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

    integer, intent(in) ::      &
       itest, jtest, rtest      ! coordinates of diagnostic point

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

    ! TODO - Make HO_GRADIENT_MARGIN_LAND the default, since it is simple and requires no optional arguments?
    ! TODO - Make ice_mask an optional argument, = 1 everywhere by default.

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

  end subroutine glissade_gradient_at_edges

!****************************************************************************

  subroutine glissade_average_to_edges(nx,           ny,           &
                                       field,                      &
                                       field_east,   field_north,  &
                                       cell_mask)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! average the field to east and north edges.
    ! Note: The east fields have dimension (nx-1,ny) and the north fields (nx,ny-1).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) ::       &
       field                    ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny), intent(out) ::    &
       field_east               ! field averaged to east edges

    real(dp), dimension(nx,ny-1), intent(out) ::    &
       field_north              ! field averaged to north edges

    integer, dimension(nx,ny), intent(in) ::        &
       cell_mask                ! average at edges only if cell_mask = 1 on either side

    ! Local variables

    integer :: i, j

    field_east(:,:) = 0.0d0
    field_north(:,:) = 0.0d0

    do j = 1, ny
       do i = 1, nx-1
          if (cell_mask(i,j) == 1  .and. cell_mask(i+1,j) == 1) then
             field_east(i,j) = 0.5d0*(field(i,j) + field(i+1,j))
          endif
       enddo
    enddo

    do j = 1, ny-1
       do i = 1, nx
          if (cell_mask(i,j) == 1  .and. cell_mask(i,j+1) == 1) then
             field_north(i,j) = 0.5d0*(field(i,j) + field(i,j+1))
          endif
       enddo
    enddo


  end subroutine glissade_average_to_edges

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

  subroutine glissade_surface_elevation_gradient(nx,           ny,        &
                                                 dx,           dy,        &
                                                 itest, jtest, rtest,     &
                                                 ice_mask,                &
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
         ice_mask,              & ! = 1 where ice is present, else = 0
         land_mask                ! = 1 for land cells, else = 0

    real(dp), dimension(nx,ny), intent(in) ::       &
         thck,                  & ! ice thickness
         usrf,                  & ! ice surface elevation
         topg                     ! bed elevation

    real(dp), intent(in) ::       &
         eus,                   & ! eustatic sea level
         thklim,                & ! minimum thickness for dynamically active ice
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

    logical, parameter :: count_max_slope = .false.
!!    logical, parameter :: count_max_slope = .true.
    integer :: max_slope_count, total_count

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

       ! Compute ds_dx and ds_dy only on edges with ice in each adjacent cell
       do j = 1, ny
          do i = 1, nx-1
             if (ice_mask(i,j) == 1 .and. ice_mask(i+1,j) == 1) then
                ds_dx_edge(i,j) = (usrf(i+1,j) - usrf(i,j)) / dx
             endif
          enddo
       enddo

       do j = 1, ny-1
          do i = 1, nx
             if (ice_mask(i,j) == 1 .and. ice_mask(i,j+1) == 1) then
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

             if (land_mask(iu,j) == 1 .and. thck_gradient_ramp > 0.0d0) then
                ! Compute a factor that reduces the gradient if ice in the upper cell is thin and land-based.
                ! This inhibits oscillations in the gradient when the thickness in the upper cell is close to thklim.
                edge_thck_upper = thck(iu,j)
                edge_factor = min(1.0d0, (edge_thck_upper - thklim)/thck_gradient_ramp)
                edge_factor = max(edge_factor, 0.0d0)
             else
                edge_factor = 1.0d0
             endif

             if (ice_mask(iu,j) == 1 .and. ice_mask(il,j) == 1) then  ! both cells have ice

                ! compute the gradient
                ds_dx_edge(i,j) = edge_factor * sign_factor * (usrf(iu,j) - usrf(il,j)) / dx

             elseif (ice_mask(iu,j) == 1 .and. land_mask(il,j) == 1) then

                ! upper cell has ice, and ice-free lower cell is land; compute the gradient
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

             if (land_mask(i,ju) == 1 .and. thck_gradient_ramp > 0.0d0) then
                ! Compute a factor that reduces the gradient if ice in the upper cell is thin and land-based.
                ! This inhibits oscillations in the gradient when the thickness in the upper cell is close to thklim.
                edge_thck_upper = thck(i,ju)
                edge_factor = min(1.0d0, (edge_thck_upper - thklim)/thck_gradient_ramp)
                edge_factor = max(edge_factor, 0.0d0)
             else
                edge_factor = 1.0d0
             endif

             if (ice_mask(i,ju)==1 .and. ice_mask(i,jl)==1) then  ! both cells have ice

                ! compute the gradient
                ds_dy_edge(i,j) = edge_factor * sign_factor * (usrf(i,ju) - usrf(i,jl)) / dy

             elseif (ice_mask(i,ju) == 1 .and. land_mask(i,jl) == 1) then

                ! upper cell has ice, and ice-free lower cell is land; compute the gradient
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

       if (count_max_slope) then

          max_slope_count = 0
          total_count = 0

          do j = nhalo+1, ny-nhalo
             do i = nhalo+1, nx-nhalo
                if (ice_mask(i,j) == 1 .and. ice_mask(i+1,j) == 1) then
                   total_count = total_count + 1
                   if (abs(ds_dx(i,j)) > max_slope) then
                      max_slope_count = max_slope_count + 1
                   endif
                endif
                if (ice_mask(i,j) == 1 .and. ice_mask(i,j+1) == 1) then
                   total_count = total_count + 1
                   if (abs(ds_dy(i,j)) > max_slope) then
                      max_slope_count = max_slope_count + 1
                   endif
                endif
             enddo
          enddo

          max_slope_count = parallel_reduce_sum(max_slope_count)
          total_count = parallel_reduce_sum(total_count)

          if (main_task) then
             write(6,*) '   max_slope_count, total_count =', max_slope_count, total_count
             write(6,*) '   cell fraction with slope limiting =', real(max_slope_count)/real(total_count)
          endif

       endif   ! count_max_slope

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

    ! Note: Halo updates moved to higher level
!    call staggered_parallel_halo(ds_dx)
!    call staggered_parallel_halo(ds_dy)

  end subroutine glissade_surface_elevation_gradient

!****************************************************************************

  subroutine glissade_slope_angle(&
       nx,           ny,    &
       dx,           dy,    &
       itest, jtest, rtest, &
       zsfc,                &
       theta_slope,         &
       theta_slope_x,       &
       theta_slope_y,       &
       phi_slope,           &
       slope_mask_in)

    ! Compute the slope angle between a surface of elevation zsfc and the horizontal.
    ! Optionally, compute the angle (in the xy plane) of the direction of steepest descent.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
         nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
         dx, dy                   ! grid cell length and width
                                  ! assumed to have the same value for each grid cell

    integer, intent(in) ::      &
         itest, jtest, rtest      ! coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(in) ::       &
         zsfc                     ! elevation of the surface whose slope is to be computed

    real(dp), dimension(nx,ny), intent(out) ::       &
         theta_slope              ! surface elevation angle (radians, 0 to pi/2) relative to the horizontal plane

    real(dp), dimension(nx,ny), intent(out), optional ::       &
         theta_slope_x,         & ! surface elevation angle in the x direction
         theta_slope_y            ! surface elevation angle in the x direction

    real(dp), dimension(nx,ny), intent(out), optional ::       &
         phi_slope                ! direction of steepest descent (radians, 0 to 2*pi), relative to the x axis

    integer, dimension(nx,ny), intent(in), optional ::  &
         slope_mask_in            ! = 1 for the part of the surface whose slope is computed, else = 0

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer, dimension(nx,ny) ::  &
         slope_mask               ! local version of slope_mask

    real(dp), dimension(nx-1,ny-1) ::    &
         stag_dz_dx, stag_dz_dy   ! gradient components of zsfc, defined at cell vertices

    real(dp), dimension(nx,ny) ::    &
         dz_dx, dz_dy,           &! gradient components of zsfc, defined at cell centers
         slope

    integer :: i, j
    real(dp) :: vector_x, vector_y
    character(len=200) :: message

    ! initialize

    if (present(slope_mask_in)) then
       slope_mask = slope_mask_in
    else
       slope_mask = 1    ! default is to compute the slope everywhere
    endif

    theta_slope = 0.0d0

    ! Compute the x and y components of the surface gradient at vertices
    ! Note: With gradient_margin_in = 1, edge gradients are computed only for edges
    !        with slope_mask = 1 on either side.
    !       For instance, if slope_mask = 1 for floating cells only, then dz_dx = dz_dy = 0
    !        for grounded regions.

    call glissade_gradient(nx,           ny,          &
                           dx,           dy,          &
                           itest, jtest, rtest,       &
                           zsfc,                      &
                           stag_dz_dx,   stag_dz_dy,  &
                           slope_mask,                &
                           gradient_margin_in = 1)

    ! Interpolate the gradient components to cell centers

    call glissade_unstagger(nx,          ny,        &
                            stag_dz_dx,  dz_dx)
    call glissade_unstagger(nx,          ny,        &
                            stag_dz_dy,  dz_dy)

    ! Compute the magnitude of the gradient
    slope = sqrt(dz_dx**2 + dz_dy**2)

    ! Compute the slope angle on the staggered grid
    ! Since slope is non-negative, theta_slope is in the range [0,pi/2)
    where (slope >= 0.0d0)
       theta_slope = atan(slope)
    elsewhere
       theta_slope = 0.0d0
    endwhere

    ! Compute slopes in the x and y directions
    if (present(theta_slope_x)) theta_slope_x = atan(abs(dz_dx))
    if (present(theta_slope_y)) theta_slope_y = atan(abs(dz_dy))

    ! Optionally, compute the direction of steepest descent
    if (present(phi_slope)) then
       phi_slope = 0.0d0
       do j = 1, ny-1
          do i = 1, nx-1
             if (slope(i,j) > 0.0d0) then
                vector_x = -dz_dx(i,j)
                vector_y = -dz_dy(i,j)
                if (abs(vector_x) > 0.0d0) then
                   phi_slope(i,j) = atan(vector_y/vector_x)
                elseif (abs(vector_y) > 0.0d0) then
                   if (vector_y > 0.0d0) then
                      phi_slope(i,j) =  pi/2.0d0
                   else
                      phi_slope(i,j) = 3.0d0*pi/2.0d0
                   endif
                endif
             else   ! vector of zero length; default to phi = 0 for flat surfaces
                phi_slope(i,j) = 0.0d0
             endif
             ! The range of atan is (-pi/2, pi/2)
             ! Add pi if the vector lies in quadrant 2 or 3
             if (vector_x < 0.0) phi_slope(i,j) = phi_slope(i,j) + pi
             ! Make sure phi is in the range [0,360)
             if (phi_slope(i,j) < 0.0d0) phi_slope(i,j) = phi_slope(i,j) + 2.0d0*pi
             if (phi_slope(i,j) > 2.0d0*pi) phi_slope(i,j) = phi_slope(i,j) - 2.0d0*pi
             ! bug check
             if (phi_slope(i,j) < 0.0d0 .or. phi_slope(i,j) >= 2.0d0*pi) then
                write(message,*) 'Error, glissade_slope_angle, phi out of range, i, j, phi =', &
                     i, j, phi_slope(i,j), phi_slope(i,j)/pi
                call write_log(message)
             endif
          enddo   ! i
       enddo   ! j
    endif   ! present(phi_slope)

    ! Note: Halo updates of theta_slope and phi_slope moved to higher level

  end subroutine glissade_slope_angle

!****************************************************************************

  subroutine glissade_slope_angle_staggered(&
       nx,           ny,    &
       dx,           dy,    &
       itest, jtest, rtest, &
       zsfc,                &
       theta_slope,         &
       theta_slope_x,       &
       theta_slope_y,       &
       phi_slope,           &
       slope_mask_in)

    ! Compute the slope angle between a surface of elevation zsfc and the horizontal.
    ! Optionally, compute the angle (in the xy plane) of the direction of steepest descent.

    ! This is like the previous subroutine, except that given an elevation field on the
    !  unstaggered grid, it computes theta_slope on the staggered grid.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
         nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
         dx, dy                   ! grid cell length and width
                                  ! assumed to have the same value for each grid cell

    integer, intent(in) ::      &
         itest, jtest, rtest      ! coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(in) ::       &
         zsfc                     ! elevation of the surface whose slope is to be computed

    real(dp), dimension(nx-1,ny-1), intent(out) ::  &
         theta_slope              ! surface elevation angle (radians, 0 to pi/2) relative to the horizontal (x-y) plane

    real(dp), dimension(nx-1,ny-1), intent(out), optional :: &
         theta_slope_x,         & ! surface elevation angle in the x direction
         theta_slope_y            ! surface elevation angle in the x direction

    real(dp), dimension(nx-1,ny-1), intent(out), optional ::  &
         phi_slope                ! direction of steepest descent (radians, 0 to 2*pi)), relative to phi = 0 along the x axis

    integer, dimension(nx,ny), intent(in), optional ::  &
         slope_mask_in            ! = 1 for the part of the surface whose slope is computed, else = 0

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer, dimension(nx,ny) ::  &
         slope_mask               ! local version of slope_mask

    real(dp), dimension(nx-1,ny-1) ::    &
         dz_dx, dz_dy,         &  ! gradient components of zsfc, defined at cell vertices
         slope                    ! slope (= magnitude of gradient) at cell vertices

    integer :: i, j
    real(dp) :: vector_x, vector_y
    character(len=200) :: message

    ! initialize

    if (present(slope_mask_in)) then
       slope_mask = slope_mask_in
    else
       slope_mask = 1    ! default is to compute the slope everywhere
    endif

    ! Compute the x and y components of the surface gradient at cell vertices.
    ! Note: With gradient_margin_in = 1, edge gradients are computed only for edges
    !        with slope_mask = 1 on either side.
    !       For instance, if slope_mask = 1 for floating cells only, then dz_dx = dz_dy = 0
    !        for grounded regions.

    call glissade_gradient(nx,           ny,          &
                           dx,           dy,          &
                           itest, jtest, rtest,       &
                           zsfc,                      &
                           dz_dx,        dz_dy,       &
                           slope_mask,                &
                           gradient_margin_in = 1)

    ! Compute the magnitude of the gradient
    slope = sqrt(dz_dx**2 + dz_dy**2)

    ! Compute the slope angle at vertices
    ! Since slope is non-negative, theta_slope is in the range [0,pi/2)
    where (slope >= 0.0d0)
       theta_slope = atan(slope)
    elsewhere
       theta_slope = 0.0d0
    endwhere

    ! Compute slopes in the x and y directions
    if (present(theta_slope_x)) theta_slope_x = atan(abs(dz_dx))
    if (present(theta_slope_y)) theta_slope_y = atan(abs(dz_dy))


    ! Optionally, compute the direction of steepest descent
    if (present(phi_slope)) then
       phi_slope = 0.0d0
       do j = 1, ny-1
          do i = 1, nx-1
             if (slope(i,j) > 0.0d0) then
                vector_x = -dz_dx(i,j)
                vector_y = -dz_dy(i,j)
                if (abs(vector_x) > 0.0d0) then
                   phi_slope(i,j) = atan(vector_y/vector_x)
                elseif (abs(vector_y) > 0.0d0) then
                   if (vector_y > 0.0d0) then
                      phi_slope(i,j) =  pi/2.0d0
                   else
                      phi_slope(i,j) = 3.0d0*pi/2.0d0
                   endif
                endif
             else   ! vector of zero length; default to phi = 0 for flat surfaces
                phi_slope(i,j) = 0.0d0
             endif
             ! The range of atan is (-pi/2, pi/2)
             ! Add pi if the vector lies in quadrants 2 or 3
             if (vector_x < 0.0) phi_slope(i,j) = phi_slope(i,j) + pi
             ! Make sure phi is in the range [0,360)
             if (phi_slope(i,j) < 0.0) phi_slope(i,j) = phi_slope(i,j) + 2.0d0*pi
             if (phi_slope(i,j) >= 2.0d0*pi) phi_slope(i,j) = phi_slope(i,j) - 2.0d0*pi
             ! bug check
             if (phi_slope(i,j) < 0.0d0 .or. phi_slope(i,j) >= 2.0d0*pi) then
                write(message,*) 'Error, glissade_slope_angle, phi out of range, i, j, phi =', &
                     i, j, phi_slope(i,j), phi_slope(i,j)/pi
                call write_log(message)
             endif
          enddo
       enddo
    endif   ! present(phi_slope)

    ! Note: Halo update of theta_slope and phi_slope moved to higher level

  end subroutine glissade_slope_angle_staggered

!****************************************************************************

  subroutine glissade_laplacian(&
       nx,        ny,          &
       dx,        dy,          &
       field,     del2_field,  &
       del2_mask)

    !----------------------------------------------------------------
    ! Given a scalar variable f of dimension (nx,ny), compute its Laplacian,
    !  d^2f/dx^2 + d^2f/dy^2.
    ! This subroutine is applied to scalar fields located at cell centers,
    !  with dimension (nx,ny).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
         nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::      &
         dx, dy                   ! grid cell length (m)

    real(dp), dimension(nx,ny), intent(in) ::       &
         field                    ! scalar field, defined at cell centers

    real(dp), dimension(nx,ny), intent(out) ::    &
         del2_field               ! Laplacian of the input field

    integer, dimension(nx,ny), intent(in), optional ::        &
         del2_mask                ! = 1 for cells to be included in the Laplacian, else = 0

    ! Local variables

    integer :: i, j

    real(dp) :: &
         df_dx_p, df_dx_m,  &     ! x-derivative terms in the Laplacian
         df_dy_p, df_dy_m         ! y-derivative terms in the Laplacian

    integer, dimension(nx,ny) :: &
         mask                     ! = input mask if present, else defaults to 1 everywhere

    if (present(del2_mask)) then
       mask = del2_mask
    else
       mask = 1
    endif

    del2_field = 0.0d0

    do j = 2, ny-1
       do i = 2, nx-1

          ! x derivative terms
          if (mask(i+1,j) == 1 .and. mask(i,j) == 1) then
             df_dx_p = (field(i+1,j) - field(i,j)) / dx
          else
             df_dx_p = 0.0d0
          endif
          if (mask(i-1,j) == 1 .and. mask(i,j) == 1) then
             df_dx_m = (field(i,j) - field(i-1,j)) / dx
          else
             df_dx_m = 0.0d0
          endif

          ! y derivative terms
          if (mask(i,j+1) == 1 .and. mask(i,j) == 1) then
             df_dy_p = (field(i,j+1) - field(i,j)) / dy
          else
             df_dy_p = 0.0d0
          endif
          if (mask(i,j-1) == 1 .and. mask(i,j) == 1) then
             df_dy_m = (field(i,j) - field(i,j-1)) / dy
          else
             df_dy_m = 0.0d0
          endif

          ! Laplacian
          del2_field(i,j) = (df_dx_p - df_dx_m)/dx + (df_dy_p - df_dy_m)/dy

       enddo
    enddo

  end subroutine glissade_laplacian

!****************************************************************************

  subroutine glissade_laplacian_stagvar(&
       nx,        ny,          &
       dx,        dy,          &
       field,     del2_field,  &
       del2_mask)

    !----------------------------------------------------------------
    ! Given a scalar variable f of dimension (nx,ny), compute its Laplacian,
    !  d^2f/dx^2 + d^2f/dy^2.
    ! This subroutine is like the one above, but it is applied to scalar fields
    !  located at vertices, with dimension (nx-1,ny-1).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
         nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::      &
         dx, dy                   ! grid cell length (m)

    real(dp), dimension(nx-1,ny-1), intent(in) ::       &
         field                    ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
         del2_field               ! Laplacian of the input field

    integer, dimension(nx-1,ny-1), intent(in), optional ::        &
         del2_mask                ! = 1 for cells to be included in the Laplacian, else = 0

    ! Local variables

    integer :: i, j

    real(dp) :: &
         df_dx_p, df_dx_m,  &     ! x-derivative terms in the Laplacian
         df_dy_p, df_dy_m         ! y-derivative terms in the Laplacian

    integer, dimension(nx-1,ny-1) :: &
         mask                     ! = input mask if present, else defaults to 1 everywhere

    if (present(del2_mask)) then
       mask = del2_mask
    else
       mask = 1
    endif

    del2_field = 0.0d0

    do j = 2, ny-2
       do i = 2, nx-2

          ! x derivative terms
          if (mask(i+1,j) == 1 .and. mask(i,j) == 1) then
             df_dx_p = (field(i+1,j) - field(i,j)) / dx
          else
             df_dx_p = 0.0d0
          endif
          if (mask(i-1,j) == 1 .and. mask(i,j) == 1) then
             df_dx_m = (field(i,j) - field(i-1,j)) / dx
          else
             df_dx_m = 0.0d0
          endif

          ! y derivative terms
          if (mask(i,j+1) == 1 .and. mask(i,j) == 1) then
             df_dy_p = (field(i,j+1) - field(i,j)) / dy
          else
             df_dy_p = 0.0d0
          endif
          if (mask(i,j-1) == 1 .and. mask(i,j) == 1) then
             df_dy_m = (field(i,j) - field(i,j-1)) / dy
          else
             df_dy_m = 0.0d0
          endif

          ! Laplacian
          del2_field(i,j) = (df_dx_p - df_dx_m)/dx + (df_dy_p - df_dy_m)/dy

       enddo
    enddo

  end subroutine glissade_laplacian_stagvar

!****************************************************************************

  subroutine glissade_laplacian_smoother(&
       nx,         ny,            &
       var,        var_smooth,    &
       smoother_mask,             &
       npoints_stencil)

    !----------------------------------------------------------------
    ! Given a 2D field on the ice grid, smooth the field using a Laplacian stencil.
    ! Uses a 9-point stencil by default, but optionally can use a 5-point
    !  or 25-point stencil.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: nx, ny                             ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) :: var             ! input field, before smoothing

    real(dp), dimension(nx,ny), intent(out) :: var_smooth     ! output field, after smoothing

    integer, dimension(nx,ny), intent(in), optional :: &
         smoother_mask                                        ! mask to identify the cells included in the smoothing

    integer, intent(in), optional :: &
         npoints_stencil                                      ! number of points in stencil, either 5 or 9

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer, dimension(nx,ny) :: rmask   ! real mask set to the optional smoother_mask, else = 1.0 everywhere

    real(dp) :: sum_mask   ! sum of mask values in the cell and its neighbors, converted to real(dp)

    integer :: npoints  ! set equal to npoints_stencil, else defaults to 9

    integer :: i, j

    if (present(smoother_mask)) then
       rmask = real(smoother_mask, dp)
    else
       rmask = 1.0d0
    endif

    if (present(npoints_stencil)) then
       npoints = npoints_stencil
       if (.not.(npoints == 5 .or. npoints == 9 .or. npoints == 25)) then
          call write_log('ERROR, glissade_laplacian_smoother: Must choose 5, 9 or 25 points for the stencil', GM_FATAL)
       endif
    else
       npoints = 9
    endif

    sum_mask = 0.0d0

    !TODO - Remove the rmask > 0 logic for n = 5 and 9?
    if (npoints == 5) then

       do j = 2, ny-1
          do i = 2, nx-1
             if (rmask(i,j) > 0.0d0) then

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

             endif   ! rmask > 0
          enddo   ! i
       enddo   ! j

    elseif (npoints == 9) then

       do j = 2, ny-1
          do i = 2, nx-1
             if (rmask(i,j) > 0.0d0) then

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

             endif   ! rmask > 0
          enddo   ! i
       enddo   ! j

    elseif (npoints == 25) then

       !Note: not yet tested
       do j = 3, ny-2
          do i = 3, nx-2

             sum_mask =  1.d0 * (rmask(i-2,j-2) + rmask(i-2,j+2) + rmask(i+2,j-2) + rmask(i+2,j+2))  &
                       + 4.d0 * (rmask(i-1,j-2) + rmask(i-2,j-1) + rmask(i-1,j+2) + rmask(i-2,j+1)  &
                               + rmask(i+2,j-1) + rmask(i+1,j-2) + rmask(i+1,j+2) + rmask(i+2,j+1)) &
                       + 6.d0 * (rmask(i-2,j)   + rmask(i,j-2)   + rmask(i+2,j)   + rmask(i,j+2))   &
                       + 16.d0 * (rmask(i-1,j-1) + rmask(i-1,j+1) + rmask(i+1,j-1) + rmask(i+1,j+1)) &
                       + 24.d0 * (rmask(i,j-1)   + rmask(i-1,j)   + rmask(i,j+1)   + rmask(i+1,j))   &
                       + 26.d0 *  rmask(i,j)

             if (sum_mask > 0.0d0) then
                var_smooth(i,j) = (1.d0/sum_mask) * &
                    (1.d0 * (rmask(i-2,j-2)*var(i-2,j-2) + rmask(i-2,j+2)*var(i+2,j+2)   &
                           + rmask(i+2,j-2)*var(i+2,j-2) + rmask(i+2,j+2)*var(i+2,j+2))  &
                   + 4.d0 * (rmask(i-1,j-2)*var(i-1,j-2) + rmask(i-2,j-1)*var(i-2,j-1)   &
                           + rmask(i-1,j+2)*var(i-1,j+2) + rmask(i-2,j+1)*var(i-2,j+1)   &
                           + rmask(i+2,j-1)*var(i+2,j-1) + rmask(i+1,j-2)*var(i+1,j-2)   &
                           + rmask(i+1,j+2)*var(i+1,j+2) + rmask(i+2,j+1)*var(i+2,j+1))  &
                   + 6.d0 * (rmask(i-2,j)*var(i-2,j)     + rmask(i,j-2)*var(i,j-2)       &
                           + rmask(i+2,j)*var(i+2,j)     + rmask(i,j+2)*var(i,j+2))      &
                  + 16.d0 * (rmask(i-1,j-1)*var(i-1,j-1) + rmask(i-1,j+1)*var(i-1,j+1)   &
                           + rmask(i+1,j-1)*var(i+1,j-1) + rmask(i+1,j+1)*var(i+1,j+1))  &
                  + 24.d0 * (rmask(i,j-1)*var(i,j-1)     + rmask(i-1,j)*var(i-1,j)       &
                           + rmask(i,j+1)*var(i,j+1)     + rmask(i+1,j)*var(i+1,j))      &
                  + 36.d0 *  rmask(i,j)*var(i,j))
             else
                var_smooth(i,j) = var(i,j)
             endif

          enddo   ! i
       enddo   ! j

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
    ! interpolate the field to elevation usrf by linear interpolation.

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

    !Note: halo update moved to higher level
!    call parallel_halo(field)

  end subroutine glissade_vertical_interpolate

!-------------------------------------------------------------------------------

  subroutine glissade_scalar_extrapolate(nx,    ny,             &
                                         itest, jtest, rtest,   &
                                         parallel,              &
                                         input_mask,            &
                                         phi_in,                &
                                         output_mask,           &
                                         phi_out,               &
                                         npoints_stencil,       &
                                         apply_smoother,        &
                                         unfilled_value)

    ! Extrapolate an input field phi_in (e.g., ice thickness) from cells with input_mask = 1,
    !  obtaining the filled phi_out in cells with output_mask = 1.
    ! Set phi = phi_unfilled in cells where extrapolation is not possible.
    !
    ! The extrapolation works as follows:
    ! (1) Initialize phi = phi_in where input_mask = 1, and phi = phi_unfilled elsewhere.
    ! (2) For each cell with output_mask = 1 and phi = phi_unfilled, obtain phi as the average value
    !     of neighbors with  phi /= phi_unfilled.
    ! (3) Repeat until no more cells can be filled with extrapolated values.

    ! Input/output arguments

    integer, intent(in) :: nx, ny                  !> horizontal grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest                       !> coordinates of diagnostic point

    type(parallel_type), intent(in) :: parallel    !> info for parallel communication

    integer, dimension(nx,ny), intent(in) :: &
         input_mask,            & !> = 1 for cells where phi_out is initially set to phi_in
         output_mask              !> = 1 for cells to which phi_in is extrapolated; should form a connected region

    real(dp), dimension(nx,ny), intent(in) :: &
         phi_in                   !> input field to be extrapolated

    real(dp), dimension(nx,ny), intent(out) :: &
         phi_out                  !> extrapolated output field

    integer, intent(in), optional :: &
         npoints_stencil          !> number of points in extrapolation stencil, either 5 or 9

    logical, intent(in), optional :: &
         apply_smoother           !> if true, then appply a Laplacian smoother after each extrapolation step

    real(dp), intent(in), optional :: &
         unfilled_value           !> value of phi in unfilled cells, to be replaced where possible by an extrapolated value

    ! Local arguments

    integer :: i, j, iglobal, jglobal, iter

    real(dp), dimension(nx,ny) :: &
         phi                     ! temporary value of phi_out

    integer, dimension(nx,ny) :: &
         filled_mask             ! = 1 for cells that have been filled, else = 0

    integer :: &
         sum_mask                ! sum of filled_mask over neighbor cells

    real(dp) ::  &
         sum_phi                 ! sum of mask*phi over neighbor cells at a given level

    integer :: &
         max_iter,             & ! max(nx,ny) * max(ewtasks, nstasks)
         local_count,          & ! local counter for filled values
         global_count,         & ! global counter for filled values
         global_count_save       ! globalcounter for filled values from previous iteration

    integer :: npoints           ! set to npoints_stencil, else defaults to 5

    logical :: smoother          ! set to apply_smoother, else defaults to false

    real(dp) :: &
         phi_unfilled            ! set to unfilled_value, else defaults to 0.0

    character(len=200) :: message

!    logical, parameter :: verbose_extrapolate = .false.
    logical, parameter :: verbose_extrapolate = .true.

    ! Initialize
    if (present(npoints_stencil)) then
       npoints = npoints_stencil
       if (.not.(npoints == 5 .or. npoints == 9)) then
          call write_log('ERROR, glissade_scalar_extrapolate: Must choose 5 or 9 points for the stencil', GM_FATAL)
       endif
    else
       npoints = 9   ! with 9 points, we typically have less horizontal and vertical striping than with 5 points
    endif

    if (present(apply_smoother)) then
       smoother = apply_smoother
    else
       smoother = .false.
    endif

    if (present(unfilled_value)) then
       phi_unfilled = unfilled_value
    else
       phi_unfilled = 0.0d0
    endif

    where (input_mask == 1)
       phi_out = phi_in
       filled_mask = 1
    elsewhere
       phi_out = phi_unfilled
       filled_mask = 0
    endwhere

    ! Count the number of filled cells

    local_count = 0
    do j = 1+nhalo, ny-nhalo
       do i = 1+nhalo,  nx-nhalo
          local_count = local_count + filled_mask(i,j)
       enddo
    enddo

    global_count_save = parallel_reduce_sum(local_count)

    ! Estimate the max number of iterations
    ! In a worst case, we would start with a few filled cells in one corner of the global domain
    !  and have to extrapolate to the opposite corner

    max_iter = (parallel%ewtasks + parallel%nstasks) * (nx + ny)

    ! Extrapolate the input field horizontally

    do iter = 1, max_iter

       if (verbose_extrapolate) then
          if (this_rank == rtest) write(6,*) 'glissade_scalar_extrapolate, iteration =', iter
       endif

       ! Make a copy of the current output field
       phi(:,:) = phi_out(:,:)

       ! Loop through all locally owned cells, extrapolating phi to cells with one or more
       !  edge neighbors that have been filled.
       ! In the end, all cells with output_mask = 1 should be filled.

       if (npoints == 5) then

          do j = 1+nhalo, ny-nhalo
             do i = 1+nhalo,  nx-nhalo
                if (output_mask(i,j) == 1 .and. filled_mask(i,j) == 0) then   ! not yet filled

                   ! Assign to this cell the average value in edge neighbors that are already filled
                   sum_mask = filled_mask(i-1,j) + filled_mask(i+1,j) + filled_mask(i,j-1) + filled_mask(i,j+1)

                   if (sum_mask > 0) then
                      sum_phi = filled_mask(i-1,j)*phi(i-1,j) + filled_mask(i+1,j)*phi(i+1,j)  &
                              + filled_mask(i,j-1)*phi(i,j-1) + filled_mask(i,j+1)*phi(i,j+1)
                      phi_out(i,j) = sum_phi/sum_mask
                   endif

                endif   ! output_mask = 1 and filled_mask = 0
             enddo   ! i
          enddo   ! j

       else   ! npoints = 9

          do j = 1+nhalo, ny-nhalo
             do i = 1+nhalo,  nx-nhalo
                if (output_mask(i,j) == 1 .and. filled_mask(i,j) == 0) then   ! not yet filled

                   ! Assign to this cell the average value in edge neighbors that are already filled

                   sum_mask = filled_mask(i-1,j+1) + filled_mask(i,j+1) + filled_mask(i+1,j+1) &
                            + filled_mask(i-1,j)   + filled_mask(i,j)   + filled_mask(i+1,j) &
                            + filled_mask(i-1,j-1) + filled_mask(i,j-1) + filled_mask(i+1,j-1)

                   if (sum_mask > 0) then
                      sum_phi = filled_mask(i-1,j+1) * phi(i-1,j+1)  &
                              + filled_mask(i,  j+1) * phi(i,  j+1)  &
                              + filled_mask(i+1,j+1) * phi(i+1,j+1)  &
                              + filled_mask(i-1,j)   * phi(i-1,j)    &
                              + filled_mask(i,  j)   * phi(i,  j)    &
                              + filled_mask(i+1,j)   * phi(i+1,j)    &
                              + filled_mask(i-1,j-1) * phi(i-1,j-1)  &
                              + filled_mask(i,  j-1) * phi(i,  j-1)  &
                              + filled_mask(i+1,j-1) * phi(i+1,j-1)
                      phi_out(i,j) = sum_phi/sum_mask
                   endif

                endif   ! output_mask = 1 and filled_mask = 0
             enddo   ! i
          enddo   ! j

       endif  ! npoints

       call parallel_halo(phi_out, parallel)

       ! Update the filled mask
       where (phi_out == phi_unfilled)
          filled_mask = 0
       elsewhere
          filled_mask = 1
       endwhere

       ! Optionally, apply a Laplacian smoother to the new phi_out
       if (smoother) then
          phi(:,:) = phi_out(:,:)
          call glissade_laplacian_smoother(nx,         ny,         &
                                           phi,        phi_out,    &
                                           filled_mask,            &
                                           npoints_stencil = npoints)
          call parallel_halo(phi_out, parallel)
       endif

       ! Every several iterations, count the number of filled cells.
       ! Exit the do loop when this number is no longer increasing.

       if (iter == 2 .or. iter == 5 .or. mod(iter, 10) == 0) then

          local_count = 0
          do j = 1+nhalo, ny-nhalo
             do i = 1+nhalo,  nx-nhalo
                local_count = local_count + filled_mask(i,j)
             enddo
          enddo

          global_count = parallel_reduce_sum(local_count)

          if (global_count == global_count_save) then
             if (verbose_extrapolate .and. main_task) &
                  write(6,*) 'Extrapolation converged: iter, global_count =', iter, global_count
             exit
          else
             if (verbose_extrapolate .and. main_task) &
                  write(6,*) 'Extrapolation convergence check: iter, global_count =', iter, global_count
             global_count_save = global_count
          endif

       endif   ! time for a convergence check

       if (iter == max_iter) then
          write(6,*) 'iter = max_iter:', max_iter
          call write_log('Extrapolation error; number of filled cells has not plateaued', GM_FATAL)
       endif

    enddo   ! max_iter

    ! Bug check: Check whether all cells with output_mask = 1 have been filled
    ! TODO: Make the warning optional?
    do j = 1, ny
       do i = 1, nx
          if (output_mask(i,j) == 1 .and. filled_mask(i,j) == 0) then
             call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!!             write(6,*) 'i, j, iglobal, jglobal:', i, j, iglobal, jglobal
             write(message,*) &
                  'Extrapolation warning: did not fill cell i, j =', iglobal, jglobal
!!             call write_log(message, GM_FATAL)
             call write_log(message)
          endif
       enddo
    enddo


  end subroutine glissade_scalar_extrapolate

!****************************************************************************

  end module glissade_grid_operators

!****************************************************************************
