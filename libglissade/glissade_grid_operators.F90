!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grid_operators.F90 - part of the Community Ice Sheet Model (CISM)  
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
    public :: glissade_stagger, glissade_unstagger,  &
              glissade_centered_gradient, glissade_upstream_gradient,  &
              glissade_edge_gradient, glissade_gradient_at_grounding_line,  &
              glissade_vertical_average

    logical, parameter :: verbose_gradient = .false.

contains

!----------------------------------------------------------------------------

  subroutine glissade_stagger(nx,           ny,        &
                              var,          stagvar,   &
                              ice_mask,     stagger_margin_in)

    ! Given a variable on the unstaggered grid (dimension nx, ny), interpolate
    ! to find values on the staggered grid (dimension nx-1, ny-1).

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) ::    &
       var                      ! unstaggered field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where values are included in the average, else = 0
                                ! Typically ice_mask = 1 where ice is present (or thck > thklim), else = 0

    integer, intent(in), optional ::   &
       stagger_margin_in        ! 0 = use all values when interpolating (including zeroes where ice is absent)
                                !   may be appropriate when computing stagusrf and stagthck on land
                                ! 1 = use only values where ice is present
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
       stagger_margin = 1  ! default is to average only over the cells with ice present
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

       ! Average over cells with ice present (ice_mask = 1)

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

    ! Given a variable on the staggered grid (dimension nx-1, ny-1), interpolate
    ! to find values on the staggered grid (dimension nx, ny).

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx-1,ny-1), intent(in) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    real(dp), dimension(nx,ny), intent(out) ::    &
       unstagvar                ! unstaggered field, defined at cell centers

    integer, dimension(nx-1,ny-1), intent(in) ::        &
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
       stagger_margin = 1  ! default is to average over cells where vmask = 1
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
                                        f,                       &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        gradient_margin_in,      &
                                        land_mask)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient is evaluated at the four neighboring points and is second-order accurate.
    !
    ! There are several choices for computing gradients at the ice margin:
    ! HO_MARGIN_GRADIENT_ALL = 0: All neighbor values are used to compute the gradient, including 
    !  values in ice-free cells.  This convention is used by Glide, but performs poorly for 
    !  ice shelves with a sudden drop in ice thickness and surface elevation at the margin.
    ! HO_MARGIN_GRADIENT_ICE_LAND = 1: Values in ice-covered and/or land cells are used to compute 
    !  the gradient, but values in ice-free ocean cells are ignored.  Where required values are 
    !  missing, the gradient is set to zero.  This reduces to option (0) for land-based problems 
    !  and (2) for ocean-based problems.
    ! HO_MARGIN_GRADIENT_ICE_ONLY = 2: Only values in ice-covered cells (i.e., cells with thck > thklim) 
    !  are used to compute gradients.  Where required values are missing, the gradient is set to zero.
    !  This option works well at shelf margins but less well for land margins (e.g., the Halfar test case).
    ! Since option (1) generally works well at both land and shelf boundaries, it is the default.
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
       f                        ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell vertices

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, construct df_fx and df_dy from the others

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, dimension(nx,ny), intent(in), optional ::        &
       land_mask                ! = 1 for land cells, else = 0

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer, dimension(nx,ny) :: mask
    integer :: summask, gradient_margin
    integer :: i, j

    !   Gradient at vertex(i,j) is based on f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)


    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    ! Initialize gradients to zero
    df_dx(:,:) = 0.d0
    df_dy(:,:) = 0.d0

    ! Set integer mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       mask(:,:) = 1             ! = 1 for all cells

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask)) then
          mask(:,:) = max(ice_mask(:,:),land_mask(:,:))    ! = 1 if ice_mask = 1 .or. land_mask = 1 
       else
          call write_log('Must pass in land mask to compute centered gradient with gradient_margin = 1', GM_FATAL)
       endif

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

          mask(:,:) = ice_mask(:,:)    ! = 1 for ice-covered cells
    endif

    ! Compute the gradients using info in cells with mask = 1

    do j = 1, ny-1
       do i = 1, nx-1

          summask = mask(i,j) + mask(i+1,j) + mask(i,j+1) + mask(i+1,j+1)

          if (summask == 4) then  ! use info in all four neighbor cells
             df_dx(i,j) = (f(i+1,j) + f(i+1,j+1) - f(i,j) - f(i,j+1)) / (2.d0 * dx)
             df_dy(i,j) = (f(i,j+1) + f(i+1,j+1) - f(i,j) - f(i+1,j)) / (2.d0 * dy)

          else  ! use info only in cells with mask = 1
                ! if info is not available, gradient component = 0

             ! df_dx
             if (mask(i,j)==1 .and. mask(i+1,j)==1) then
                df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
             elseif (mask(i,j+1)==1 .and. mask(i+1,j+1)==1) then
                df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
             endif

             ! df_dy
             if (mask(i,j)==1 .and. mask(i,j+1)==1) then
                df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
             elseif (mask(i+1,j)==1 .and. mask(i+1,j+1)==1) then
                df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
             endif

          endif

       enddo    ! i
    enddo       ! j

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Centered gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       
       print*, ' '
       print*, 'df_dy:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_centered_gradient

!****************************************************************************

  subroutine glissade_upstream_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        f,                       &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        gradient_margin_in,      &
                                        accuracy_flag_in,        &
                                        land_mask)

    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    !  compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient can be evaluated at two upstream points (for first-order accuracy) 
    !  or at four upstream points (for second-order accuracy).
    ! Note: Upstream is defined by the direction of higher surface elevation
    !  rather than the direction the flow is coming from (though these are
    !  usually the same).
    !
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       f                        ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, intent(in), optional ::    &
       accuracy_flag_in         ! = 1 for 1st order, 2 for 2nd order

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, construct df_fx and df_dy from the others

    integer, dimension(nx,ny), intent(in), optional ::        &
       land_mask                ! = 1 for land cells, else = 0

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer, dimension(nx,ny) :: mask
    integer :: i, j
    real(dp) :: sum1, sum2
    integer :: gradient_margin, accuracy_flag, summask

    !   First-order upstream gradient at vertex(i,j) is based on two points out of f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)
    !
    !   Second-order gradient is based on four points in the upstream direction

    if (present(accuracy_flag_in)) then
       accuracy_flag = accuracy_flag_in
    else
       accuracy_flag = 2   ! default to second-order
    endif

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    ! Initialize gradients to zero
    df_dx(:,:) = 0.d0
    df_dy(:,:) = 0.d0

    ! Set integer mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       mask(:,:) = 1             ! = 1 for all cells

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask)) then
          mask(:,:) = max(ice_mask(:,:),land_mask(:,:))    ! = 1 if ice_mask = 1 .or. land_mask = 1 
       else
          call write_log('Must pass in land mask to compute upstream gradient with gradient_margin = 1', GM_FATAL)
       endif

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       mask(:,:) = ice_mask(:,:)    ! = 1 for ice-covered cells

    endif

    if (accuracy_flag == 1) then   ! first-order accurate

       do j = 1, ny-1
          do i = 1, nx-1

             ! Compute gradient only if at least one neighbor is ice-covered
             summask = ice_mask(i,j) + ice_mask(i+1,j) + ice_mask(i,j+1) + ice_mask(i+1,j+1)
           
             if (summask > 0) then

                ! Compute df_dx by taking upstream gradient

                sum1 = f(i+1,j+1) + f(i,j+1)
                sum2 = f(i+1,j) + f(i,j)

                if (sum1 > sum2 .and. mask(i+1,j+1)==1 .and. mask(i,j+1)==1) then
                   df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
                elseif (sum1 <= sum2 .and. mask(i+1,j)==1 .and. mask(i,j)==1) then
                   df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
                endif

                ! Compute df_dy by taking upstream gradient
             
                sum1 = f(i+1,j+1) + f(i+1,j)
                sum2 = f(i,j+1) + f(i,j)
             
                if (sum1 > sum2 .and. mask(i+1,j+1)==1 .and. mask(i+1,j)==1) then
                   df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
                elseif (sum1 <= sum2 .and. mask(i,j+1)==1 .and. mask(i,j)==1) then
                   df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
                else
                   df_dy(i,j) = 0.d0
                endif

             endif  ! summask > 0 (mask = 1 in at least one neighbor cell)

          enddo
       enddo

    else    ! second-order accurate

       do j = 2, ny-2   ! loop does not include all of halo
          do i = 2, nx-2

             ! Compute gradient only if at least one neighbor is ice-covered
             summask = ice_mask(i,j) + ice_mask(i+1,j) + ice_mask(i,j+1) + ice_mask(i+1,j+1)
           
             if (summask > 0) then

                ! Compute df_dx by taking upstream gradient
             
                ! determine upstream direction

                sum1 = f(i+1,j+1) + f(i,j+1) + f(i+1,j+2) + f(i,j+2)
                sum2 = f(i+1,j) + f(i,j) + f(i+1,j-1) + f(i,j-1)

                if (sum1 > sum2) then

                   summask = mask(i+1,j+1) + mask(i,j+1) + mask(i+1,j+2) + mask(i,j+2)

                   if (summask == 4) then ! use info in all four upstream neighbor cells
                      df_dx(i,j) = (1.5d0 * (f(i+1,j+1) - f(i,j+1))     &
                                  - 0.5d0 * (f(i+1,j+2) - f(i,j+2))) / dx
                   elseif (mask(i+1,j+1)==1 .and. mask(i,j+1)==1) then   ! revert to 1st order, using upstream info
                      print*, 'df_dx: i, j, summask =', i, j, summask
                      df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
                   endif

                else  ! sum1 <= sum2

                   summask = mask(i+1,j) + mask(i,j) + mask(i+1,j-1) + mask(i,j-1)

                   if (summask == 4) then ! use info in all four upstream neighbor cells
                      df_dx(i,j) = (1.5d0 * (f(i+1,j)   - f(i,j))     &
                                  - 0.5d0 * (f(i+1,j-1) - f(i,j-1))) / dx
                   elseif (mask(i+1,j)==1 .and. mask(i,j)==1) then   ! revert to 1st order, using upstream info
                      print*, 'df_dx: i, j, summask =', i, j, summask
                      df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
                   endif
                   
                endif   ! sum1 > sum2

                ! Compute df_dy by taking upstream gradient

                ! determine upstream direction

                sum1 = f(i+1,j+1) + f(i+1,j) + f(i+2,j+1) + f(i+2,j)
                sum2 = f(i,j+1) + f(i,j) + f(i-1,j+1) + f(i-1,j)
             
                if (sum1 > sum2) then

                   summask = mask(i+1,j+1) + mask(i+1,j) + mask(i+2,j+1) + mask(i+2,j)

                   if (summask == 4) then ! use info in all four upstream neighbor cells
                      df_dy(i,j) = (1.5d0 * (f(i+1,j+1) - f(i+1,j))     &
                                  - 0.5d0 * (f(i+2,j+1) - f(i+2,j))) / dy
                   elseif (mask(i+1,j+1)==1 .and. mask(i+1,j)==1) then   ! revert to 1st order, using upstream info
                      print*, 'df_dy: i, j, summask =', i, j, summask
                      df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
                   endif

                else   ! sum1 <= sum2

                   summask = mask(i,j+1) + mask(i,j) + mask(i-1,j+1) + mask(i-1,j)
                   
                   if (summask == 4) then ! use info in all four upstream neighbor cells
                      df_dy(i,j) = (1.5d0 * (f(i,j+1)   - f(i,j))     &
                                  - 0.5d0 * (f(i-1,j+1) - f(i-1,j))) / dy
                   elseif (mask(i+1,j+1)==1 .and. mask(i+1,j)==1) then   ! revert to 1st order, using upstream info
                      print*, 'df_dy: i, j, summask =', i, j, summask
                      df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
                   endif

                endif   ! sum1 > sum2

             endif      ! summask > 0 (mask = 1 in at least one neighbor cell)

          enddo     ! i
       enddo        ! j

       ! fill in halo values
       call staggered_parallel_halo(df_dx)
       call staggered_parallel_halo(df_dy)

    endif   ! 1st or 2nd order accurate

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

  subroutine glissade_edge_gradient(nx,           ny,        &
                                    dx,           dy,        &
                                    f,                       &
                                    df_dx,        df_dy,     &
                                    gradient_margin_in,      &
                                    ice_mask,     land_mask)

    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) at cell edges (i.e., the C grid):
    ! df_dx at the midpoint of the east edge and df_dy at the midpoint of
    ! the north edge.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       f                        ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell edges

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, set gradient to zero
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, set gradient to zero

    integer, dimension(nx,ny), intent(in), optional ::        &
       ice_mask,     &          ! = 1 where ice is present, else = 0
       land_mask                ! = 1 for land cells, else = 0

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer, dimension(nx,ny) :: mask
    integer :: gradient_margin
    integer :: i, j

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

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    ! Initialize gradients to zero
    df_dx(:,:) = 0.d0
    df_dy(:,:) = 0.d0

    ! Set integer mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       mask(:,:) = 1             ! = 1 for all cells

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask) .and. present(ice_mask)) then
          mask(:,:) = max(ice_mask(:,:),land_mask(:,:))    ! = 1 if ice_mask = 1 .or. land_mask = 1 
       else
          call write_log('Must pass in land and ice masks to compute edge gradient with gradient_margin = 1', GM_FATAL)
       endif

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       if (present(ice_mask)) then
          mask(:,:) = ice_mask(:,:)    ! = 1 for ice-covered cells
       else
          call write_log('Must pass in ice mask to compute edge gradient with gradient_margin = 2', GM_FATAL)
       endif

    endif

    ! Compute the gradients using info in cells with mask = 1

    do j = 1, ny-1
       do i = 1, nx-1

          ! df_dx

          if (mask(i,j)==1 .and. mask(i+1,j)==1) then
             df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
          endif

          ! df_dy

          if (mask(i,j)==1 .and. mask(i,j+1)==1) then
             df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
          endif

       enddo    ! i
    enddo       ! j

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Edge gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'df_dy:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_edge_gradient

!----------------------------------------------------------------------------
  subroutine glissade_gradient_at_grounding_line(nx,             ny,         &
                                                 dx,             dy,         &
                                                 f_pattyn,       f_ground,   &
                                                 ice_mask,       usrf,       &
!!                                                 dusrf_dx,       dusrf_dy)
                                                 dusrf_dx,       dusrf_dy,   &
                                                 itest,          jtest)

    !WHL - debug - passing in itest and jtest for now

    !----------------------------------------------------------------
    ! Compute a surface elevation gradient for vertices adjacent to the grounding line
    ! (i.e., vertices with 0 < f_ground < 1).
    ! 
    ! Here is the procedure for computing ds/dx adjacent to the GL.
    ! (The procedure is analogous for ds/dy.)
    ! Consider this block of square cells, where each cell has dimension 1 x 1:
    !       ________________________________________________
    !      |           |           |           |           |
    !      |           |           |           |           |
    !      |     *     -     *     -     *     -     *     |
    !      |           |           |           |           |
    !      |___________|___________|___________|___________|
    !      |           |           |           |           |
    !      |           |           |           |           |
    !      |     *     -     *     -     *     -     *     |
    !      |           |           |           |           |
    !      |___________|___________|___________|___________|
    !
    ! Asterisks lie at cell centers, and hyphens lie at the midpoint of vertical edges.
    !
    ! Suppose we want to compute ds/dx at the central vertex (0,0).
    ! A standard centered difference would be computed as
    !  
    !     ds/dx(0,0) = 1/2 * ((s(0.5,0.5) - s(-0.5,0.5)) + (s(0.5,-0.5) - s(-0.5,-0.5)))
    !
    ! In other words, ds/dx is the mean of the edge gradients at (0, ±0.5).
    !
    ! But suppose 0 < f_ground < 1 for the vertex at (0,0). This implies that at least 
    ! one neighbor cell, but not more than three neighbor cells, are floating. (A cell
    ! is floating if f_pattyn = rhow*b/(rhoi*H) > 1.)  We would like to avoid including in the gradient 
    ! any edges that join a floating cell to a grounded cell, because ds/dx is discontinuous
    ! at the GL. Instead, we want to compute gradients based on fully grounded edges 
    ! (i.e., edges with grounded cells, f_pattyn <= 1, on either side) or fully floating edges.
    !
    ! We estimate (ds/dx)_g as follows:
    ! (1) Consider the two edges at (0, ±0.5). If either edge is fully grounded, take
    !     (ds/dx)_g to be the gradient at this edge. By assumption, both edges cannot
    !     be fully grounded.
    ! (2) If neither of the edges at (0, ±0.5) is fully grounded, then consider the
    !     four edges at (±1, ±0.5). If one or more of these edges is fully grounded, take
    !     (ds/dx)_g as the mean of the gradients at these fully grounded edges.
    ! (3) If none of the edges at (±1, ±0.5) is fully grounded, then consider again
    !     the two edges at (0, ±0.5). At least one, and perhaps both, of these edges
    !     is partly grounded or 'mixed' (i.e., it borders one grounded cell). 
    !     Take (ds/dx)_g as the mean of the gradient at the mixed edges.
    !
    ! Similarly, we compute (ds/dx)_f by first searching for fully floating edges 
    ! (steps 1 and 2), and then if necessary using the gradient at the nearest partly 
    ! floating edges.
    !
    ! Finally, set (ds_dx) = f_ground * (ds/dx)_g + (1 - f_ground) * (ds/dx)_f.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 for cells where ice is present, else = 0

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       f_ground                 ! grounded ice fraction at vertex, 0 <= f_ground <= 1
                                ! set to -1 where vmask = 0

    real(dp), dimension(nx,ny), intent(in) :: &
       f_pattyn,              & ! Pattyn flotation function, -rhoo*(topg-eus) / (rhoi*thck)
       usrf                     ! upper surface elevation (m)

    real(dp), dimension(nx-1,ny-1), intent(inout) :: &
       dusrf_dx, dusrf_dy       ! gradient of upper surface elevation (m/m) 

    integer, intent(in) ::   &
       itest, jtest    ! test points
    
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    real(dp) :: ds_dx_g   ! x gradient for grounded ice
    real(dp) :: ds_dx_f   ! x gradient for floating ice
    real(dp) :: ds_dy_g   ! y gradient for grounded ice
    real(dp) :: ds_dy_f   ! y gradient for floating ice
    
    integer :: numedges   ! number of edge gradients to be averaged

    ! Note: Vertical edges with computable gradients have dimension (nx-1,ny), and
    !       horizontal edges with computable gradients have dimension (nx,ny-1).
    !       Here we assign dimensions (nx,ny) to all edge arrays for simplicity.
 
    real(dp), dimension(nx,ny) :: ds_dx_edge   ! ds/dx at vertical edges
    real(dp), dimension(nx,ny) :: ds_dy_edge   ! ds/dy at horizontal edges
    
    integer, dimension(nx,ny) :: grounded_edge_mask   ! = 1 for fully grounded edges, else = 0
    integer, dimension(nx,ny) :: floating_edge_mask   ! = 1 for fully floating edges, else = 0
    integer, dimension(nx,ny) :: mixed_edge_mask      ! = 1 for mixed edges, else = 0

    real(dp), parameter :: eps11 = 1.d-11       ! small number

    !WHL - debug
    logical, parameter :: verbose_gl_gradient = .true.

    !WHL - debug
    if (verbose_gl_gradient) then
       print*, ' '
       print*, 'Starting ds/dx gradient at GL, itest, jtest =', itest, jtest
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(f10.6)',advance='no') dusrf_dx(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Starting ds/dy gradient at GL:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(f10.6)',advance='no') dusrf_dy(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    !WHL - debug - Create some new thickness profiles
!!    i = itest
!!    j = jtest
    ! ground a floating point downstream
!!    f_pattyn(i+1,j) = 1.d0  
!!    f_ground(i,j) = 1.d0

    ! float 3 grounded points upstream
!!    f_pattyn(i-1,j-1:j+1) = 2.d0

    ! Initialize ds/dx and integer masks at vertical edges

    ds_dx_edge(:,:) = 0.d0
    grounded_edge_mask(:,:) = 0
    floating_edge_mask(:,:) = 0
    mixed_edge_mask(:,:) = 0

    ! Characterize vertical edges as fully grounded, fully floating, or mixed,
    ! and compute ds/dx at these edges.  
    ! Vertical edge (i,j) lies to the right of cell center (i,j) and below vertex (i,j). 

    do j = 1, ny     ! loop over vertical edges
       do i = 1, nx-1
          if (ice_mask(i,j) == 1 .and. ice_mask(i+1,j) == 1) then   ! both neighbor cells are ice-covered
             ds_dx_edge(i,j) = (usrf(i+1,j) - usrf(i,j)) / dx
             if (f_pattyn(i,j) <= 1.d0 .and. f_pattyn(i+1,j) <= 1.d0) then    ! edge is fully grounded
                grounded_edge_mask(i,j) = 1
             elseif (f_pattyn(i,j) > 1.d0 .and. f_pattyn(i+1,j) > 1.d0) then  ! edge is fully floating
                floating_edge_mask(i,j) = 1
             else
                mixed_edge_mask(i,j) = 1
             endif
          endif
       enddo
    enddo

    !WHL - debug
    if (verbose_gl_gradient) then
       print*, ' '
       print*, 'f_pattyn at cell center:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(f10.6)',advance='no') f_pattyn(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'f_ground at vertex:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(f10.6)',advance='no') f_ground(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Grounded edge mask for ds/dx:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') grounded_edge_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Floating edge mask for ds/dx:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') floating_edge_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Mixed edge mask for ds/dx:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') mixed_edge_mask(i,j)
          enddo
          write(6,*) ' '
       enddo    
    endif

    ! Compute (ds/dx)_g and (ds/dx)_f at each vertex, then take weighted mean to compute ds/dx.
    do j = 1, ny-1   ! loop over vertices 
       do i = 1, nx-1
          if (f_ground(i,j) > eps11 .and. f_ground(i,j) < (1.d0 - eps11)) then  ! vertex is adjacent to the GL             

             ! Compute the grounded x gradient, (ds/dx)_g, based on gradients at nearby
             ! fully grounded edges (if possible), otherwise based on gradients at mixed edges.

             if (grounded_edge_mask(i,j+1) == 1) then
                ds_dx_g = ds_dx_edge(i,j+1)
             elseif (grounded_edge_mask(i,j) == 1) then
                ds_dx_g = ds_dx_edge(i,j)
             elseif (grounded_edge_mask(i-1,j+1) == 1 .or. grounded_edge_mask(i+1,j+1) == 1  .or.  &
                     grounded_edge_mask(i-1,j)   == 1 .or. grounded_edge_mask(i+1,j)   == 1) then
                numedges = grounded_edge_mask(i-1,j+1) + grounded_edge_mask(i+1,j+1) +   &
                           grounded_edge_mask(i-1,j)   + grounded_edge_mask(i+1,j)
                ds_dx_g = (grounded_edge_mask(i-1,j+1) * ds_dx_edge(i-1,j+1) +  &
                           grounded_edge_mask(i-1,j)   * ds_dx_edge(i-1,j)   +  &
                           grounded_edge_mask(i+1,j+1) * ds_dx_edge(i+1,j+1) +  &
                           grounded_edge_mask(i+1,j)   * ds_dx_edge(i+1,j))     &
                           / numedges
             elseif (mixed_edge_mask(i,j+1) == 1 .or. mixed_edge_mask(i,j) == 1) then
                numedges = mixed_edge_mask(i,j+1) + mixed_edge_mask(i,j)
                ds_dx_g = (mixed_edge_mask(i,j+1) * ds_dx_edge(i,j+1) + &
                           mixed_edge_mask(i,j)   * ds_dx_edge(i,j))   &
                           / numedges
             else   ! punt (TODO - Check whether this ever happens)
                ds_dx_g = 0.d0
             endif

             ! Compute the floating x gradient, (ds/dx)_f, based on gradients at nearby
             ! fully floating edges (if possible), otherwise based on gradients at mixed edges.

             if (floating_edge_mask(i,j+1) == 1) then
                ds_dx_f = ds_dx_edge(i,j+1)
             elseif (floating_edge_mask(i,j) == 1) then
                ds_dx_f = ds_dx_edge(i,j)
             elseif (floating_edge_mask(i-1,j+1) == 1 .or. floating_edge_mask(i+1,j+1) == 1  .or.  &
                     floating_edge_mask(i-1,j)   == 1 .or. floating_edge_mask(i+1,j)   == 1) then
                numedges = floating_edge_mask(i-1,j+1) + floating_edge_mask(i+1,j+1) +   &
                           floating_edge_mask(i-1,j)   + floating_edge_mask(i+1,j)
                ds_dx_f = (floating_edge_mask(i-1,j+1) * ds_dx_edge(i-1,j+1) +  &
                           floating_edge_mask(i-1,j)   * ds_dx_edge(i-1,j)   +  &
                           floating_edge_mask(i+1,j+1) * ds_dx_edge(i+1,j+1) +  &
                           floating_edge_mask(i+1,j)   * ds_dx_edge(i+1,j)) &
                           / numedges
             elseif (mixed_edge_mask(i,j+1) == 1 .or. mixed_edge_mask(i,j) == 1) then
                numedges = mixed_edge_mask(i,j+1) + mixed_edge_mask(i,j)
                ds_dx_f = (mixed_edge_mask(i,j+1) * ds_dx_edge(i,j+1) + &
                           mixed_edge_mask(i,j)   * ds_dx_edge(i,j))   &
                           / numedges
             else   ! punt
                ds_dx_f = 0.d0
             endif

             !WHL - debug
             if (i==itest .and. j==jtest .and. verbose_gl_gradient) then
                print*, 'i, j, f_ground =', i, j, f_ground(i,j)
                print*, 'ds_dx_g, ds_dx_f =', ds_dx_g, ds_dx_f
             endif

             ! Compute ds/dx as a weighted average of (ds/dx)_g and (dx/dx)_f
             dusrf_dx(i,j) = f_ground(i,j) * ds_dx_g + (1.d0 - f_ground(i,j)) * ds_dx_f

          endif   ! vertex is adjacent to the GL
       enddo      ! i
    enddo         ! j

    ! Initialize ds/dy and integer masks at horizontal edges

    ds_dy_edge(:,:) = 0.d0
    grounded_edge_mask(:,:) = 0
    floating_edge_mask(:,:) = 0
    mixed_edge_mask(:,:) = 0

    ! Characterize horizontal edges as fully grounded, fully floating, or mixed,
    ! and compute ds/dy at these edges.  
    ! Horizontal edge (i,j) lies above cell center (i,j) and to the left of vertex (i,j). 

    do j = 1, ny     ! loop over horizontal edges
       do i = 1, nx-1
          if (ice_mask(i,j) == 1 .and. ice_mask(i,j+1) == 1) then   ! both neighbor cells are ice-covered
             ds_dy_edge(i,j) = (usrf(i,j+1) - usrf(i,j)) / dy
             if (f_pattyn(i,j) <= 1.d0 .and. f_pattyn(i,j+1) <= 1.d0) then    ! edge is fully grounded
                grounded_edge_mask(i,j) = 1
             elseif (f_pattyn(i,j) > 1.d0 .and. f_pattyn(i,j+1) > 1.d0) then  ! edge is fully floating
                floating_edge_mask(i,j) = 1
             else
                mixed_edge_mask(i,j) = 1
             endif
          endif
       enddo
    enddo

    !WHL - debug
    if (verbose_gl_gradient) then
       print*, ' '
       print*, 'f_pattyn at cell center:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(f10.6)',advance='no') f_pattyn(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'f_ground at vertex:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(f10.6)',advance='no') f_ground(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Grounded edge mask for ds/dy:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') grounded_edge_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Floating edge mask for ds/dy:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') floating_edge_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Mixed edge mask for ds/dy:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(i10)',advance='no') mixed_edge_mask(i,j)
          enddo
          write(6,*) ' '
       enddo    
    endif

    ! Compute (ds/dy)_g and (ds/dy)_f at each vertex, then take weighted mean to compute ds/dy.
    do j = 1, ny-1   ! loop over vertices
       do i = 1, nx-1
          if (f_ground(i,j) > eps11 .and. f_ground(i,j) < (1.d0 - eps11)) then  ! vertex is adjacent to the GL             

             ! Compute the grounded y gradient, (ds/dy)_g, based on gradients at nearby
             ! fully grounded edges (if possible), otherwise based on gradients at mixed edges.

             if (grounded_edge_mask(i+1,j) == 1) then
                ds_dy_g = ds_dy_edge(i+1,j)
             elseif (grounded_edge_mask(i,j) == 1) then
                ds_dy_g = ds_dy_edge(i,j)
             elseif (grounded_edge_mask(i,j+1) == 1 .or. grounded_edge_mask(i+1,j+1) == 1  .or.  &
                     grounded_edge_mask(i,j-1) == 1 .or. grounded_edge_mask(i+1,j-1) == 1) then
                numedges = grounded_edge_mask(i,j+1) + grounded_edge_mask(i+1,j+1) +   &
                           grounded_edge_mask(i,j-1) + grounded_edge_mask(i+1,j-1)
                ds_dy_g = (grounded_edge_mask(i,j+1)   * ds_dy_edge(i,j+1)   +  &
                           grounded_edge_mask(i,j-1)   * ds_dy_edge(i,j-1)   +  &
                           grounded_edge_mask(i+1,j+1) * ds_dy_edge(i+1,j+1) +  &
                           grounded_edge_mask(i+1,j-1) * ds_dy_edge(i+1,j-1))   &
                           / numedges
             elseif (mixed_edge_mask(i+1,j) == 1 .or. mixed_edge_mask(i,j) == 1) then
                numedges = mixed_edge_mask(i+1,j) + mixed_edge_mask(i,j)
                ds_dy_g = (mixed_edge_mask(i+1,j) * ds_dy_edge(i+1,j) + &
                           mixed_edge_mask(i,j)   * ds_dy_edge(i,j))   &
                           / numedges
             else   ! punt (TODO - Check whether this ever happens)
                ds_dy_g = 0.d0
             endif

             ! Compute the floating y gradient, (ds/dy)_f, based on gradients at nearby
             ! fully floating edges (if possible), otherwise based on gradients at mixed edges.

             if (floating_edge_mask(i+1,j) == 1) then
                ds_dy_f = ds_dy_edge(i+1,j)
             elseif (floating_edge_mask(i,j) == 1) then
                ds_dy_f = ds_dy_edge(i,j)
             elseif (floating_edge_mask(i,j+1) == 1 .or. floating_edge_mask(i+1,j+1) == 1  .or.  &
                     floating_edge_mask(i,j-1) == 1 .or. floating_edge_mask(i+1,j-1) == 1) then
                numedges = floating_edge_mask(i,j+1) + floating_edge_mask(i+1,j+1) +   &
                           floating_edge_mask(i,j-1) + floating_edge_mask(i+1,j-1)
                ds_dy_f = (floating_edge_mask(i,j+1)   * ds_dy_edge(i,j+1)   +  &
                           floating_edge_mask(i,j-1)   * ds_dy_edge(i,j-1)   +  &
                           floating_edge_mask(i+1,j+1) * ds_dy_edge(i+1,j+1) +  &
                           floating_edge_mask(i+1,j-1) * ds_dy_edge(i+1,j-1))   &
                           / numedges
             elseif (mixed_edge_mask(i+1,j) == 1 .or. mixed_edge_mask(i,j) == 1) then
                numedges = mixed_edge_mask(i+1,j) + mixed_edge_mask(i,j)
                ds_dy_f = (mixed_edge_mask(i+1,j) * ds_dy_edge(i+1,j) + &
                           mixed_edge_mask(i,j)   * ds_dy_edge(i,j))   &
                           / numedges
             else   ! punt (TODO - Check whether this ever happens)
                ds_dy_f = 0.d0
             endif

             ! Compute ds/dy as a weighted average of (ds/dx)_g and (dx/dx)_f
             dusrf_dy(i,j) = f_ground(i,j) * ds_dy_g + (1.d0 - f_ground(i,j)) * ds_dy_f

          endif   ! vertex is adjacent to the GL
       enddo      ! i
    enddo         ! j

    !WHL - debug
    if (verbose_gl_gradient) then
       print*, ' '
       print*, 'New ds/dx gradient at GL, itest, jtest =', itest, jtest
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(f10.6)',advance='no') dusrf_dx(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'New ds/dy gradient at GL:'
       do j = jtest+2, jtest-2, -1
          do i = itest-3, itest+3
             write(6,'(f10.6)',advance='no') dusrf_dy(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_gradient_at_grounding_line

!----------------------------------------------------------------------------

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
