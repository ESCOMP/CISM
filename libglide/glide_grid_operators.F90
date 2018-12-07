!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_grid_operators.F90 - part of the Community Ice Sheet Model (CISM)  
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

! Various grid operators for the Glide dycore, including routines for computing gradients
!  and switching between staggered and unstaggered grids.
! Also includes some functions previously used by the Glam dycore.

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_nan.inc"
#include "glide_mask.inc"

module glide_grid_operators

    use glimmer_global, only : dp
    implicit none

contains

!----------------------------------------------------------------------------

  subroutine glide_geometry_derivs(model)

! Compute geometric quantities needed by the glide dycore:
! stagthck (given thck), along with the gradients
! dusrfdew/dns (given usrf) and dthckdew/dns (given thck).

    use glide_types, only: glide_global_type

    implicit none

    type(glide_global_type), intent(inout) :: model

    ! Interpolate ice thickness to velocity points
 
    call stagvarb(model%geometry% thck, &
                  model%geomderv% stagthck,      &
                  model%general%  ewn,           &
                  model%general%  nsn)

    ! Compute EW and NS gradients in usrf and thck

    call geomders(model%numerics, &
                  model%geometry% usrf, &
                  model%geomderv% stagthck,&
                  model%geomderv% dusrfdew, &
                  model%geomderv% dusrfdns)

    call geomders(model%numerics, &
                  model%geometry% thck, &
                  model%geomderv% stagthck,&
                  model%geomderv% dthckdew, &
                  model%geomderv% dthckdns)

    !NOTE: The following commented-out code is included in stagthickness.
!     where (model%geomderv%stagthck == 0.d0)
!            model%geomderv%dusrfdew = 0.d0
!            model%geomderv%dusrfdns = 0.d0
!            model%geomderv%dthckdew = 0.d0
!            model%geomderv%dthckdns = 0.d0
!     endwhere

  end subroutine glide_geometry_derivs

!---------------------------------------------------------------

  subroutine stagvarb(ipvr,opvr,ewn,nsn)

  ! Interpolate a scalar variable such as ice thickness from cell centers to cell corners.
  ! Note:  zero-thickness values are included in the average.

    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:)  :: ipvr
    
    integer, intent(in) :: ewn,nsn

    opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                             ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0

  end subroutine stagvarb

!----------------------------------------------------------------------------

  subroutine stagvarb_3d(ipvr, opvr, ewn, nsn, upn)
    real(dp), intent(in), dimension(:,:,:) :: ipvr
    real(dp), intent(out), dimension(:,:,:) :: opvr
    integer, intent(in) :: ewn, nsn, upn
    integer :: k

    do k = 1, upn
        call stagvarb(ipvr(k,:,:), opvr(k,:,:), ewn, nsn)
    end do

  end subroutine stagvarb_3d

!----------------------------------------------------------------------------

  subroutine stagvarb_mask(ipvr,opvr,ewn,nsn,geometry_mask)

    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:)  :: ipvr
    
    integer, intent(in) :: ewn,nsn
    integer, intent(in), dimension(:,:) :: geometry_mask
    integer :: ew,ns,n
    real(dp) :: tot

        opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                                 ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0

        do ns = 1,nsn-1
            do ew = 1,ewn-1

                !If any of our staggering points are shelf front, ignore zeros when staggering
                if (any(GLIDE_NO_ICE(geometry_mask(ew:ew+1, ns:ns+1)))) then
                    n = 0
                    tot = 0
    
                    if (GLIDE_HAS_ICE(geometry_mask(ew,ns))) then
                        tot = tot + ipvr(ew,ns)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew+1,ns))) then
                        tot = tot + ipvr(ew+1,ns)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew,ns+1))) then
                        tot = tot + ipvr(ew,ns+1)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew+1,ns+1))) then
                        tot = tot + ipvr(ew+1,ns+1)
                        n   = n   + 1
                    end if
                    if (n > 0) then
                        opvr(ew,ns) = tot/n
                    else
                        opvr(ew,ns) = 0
                    end if
                
                !Standard Staggering
                else
                        opvr(ew,ns) = (ipvr(ew+1,ns) + ipvr(ew,ns+1) + &
                                      ipvr(ew+1,ns+1) + ipvr(ew,ns)) / 4.0d0
                end if
  
        end do
    end do

  end subroutine stagvarb_mask

!----------------------------------------------------------------------------

  subroutine stagvarb_3d_mask(ipvr, opvr, ewn, nsn, upn, geometry_mask)
    real(dp), intent(in), dimension(:,:,:) :: ipvr
    real(dp), intent(out), dimension(:,:,:) :: opvr
    integer, intent(in) :: ewn, nsn, upn
    integer, intent(in), dimension(:,:) :: geometry_mask
    integer :: k

    do k = 1, upn
        call stagvarb_mask(ipvr(k,:,:), opvr(k,:,:), ewn, nsn, geometry_mask)
    end do

  end subroutine stagvarb_3d_mask

!----------------------------------------------------------------------------

  subroutine geomders(numerics,ipvr,stagthck,opvrew,opvrns)

    use glide_types, only: glide_numerics

    implicit none 

    type(glide_numerics) :: numerics
    real(dp), intent(out), dimension(:,:) :: opvrew, opvrns
    real(dp), intent(in), dimension(:,:) :: ipvr, stagthck

    real(dp) :: dew2, dns2 
    integer :: ew,ns,ewn,nsn

    ! Obviously we don't need to do this every time,
    ! but will do so for the moment.
    dew2 = 1.d0/(2.0d0 * numerics%dew)
    dns2 = 1.d0/(2.0d0 * numerics%dns)
    ewn=size(ipvr,1)
    nsn=size(ipvr,2)

    do ns=1,nsn-1
       do ew = 1,ewn-1
          if (stagthck(ew,ns) /= 0.0d0) then
             opvrew(ew,ns) = (ipvr(ew+1,ns+1)+ipvr(ew+1,ns)-ipvr(ew,ns)-ipvr(ew,ns+1)) * dew2
             opvrns(ew,ns) = (ipvr(ew+1,ns+1)+ipvr(ew,ns+1)-ipvr(ew,ns)-ipvr(ew+1,ns)) * dns2
          else
             opvrew(ew,ns) = 0.
             opvrns(ew,ns) = 0.
          end if
       end do
    end do
    
  end subroutine geomders

!----------------------------------------------------------------------------

!> NOTE: The rest of this module contains various functions for computing derivatives.
!> These functions previously were used by the Glam dycore.
!> They are not used by Glissade, but are retained in case they are useful.
!> Also, subroutine df_field_2d is called from Glide's basal water routing scheme.

    !------------------------------------------------------------------
    !First Derivative Estimates, Second Order, 2D
    !------------------------------------------------------------------

    !> Computes derivative fields of the given function.
    subroutine df_field_2d(f,  &
                           deltax,      deltay,      &
                           out_dfdx,    out_dfdy,    &
                           direction_x, direction_y)

    ! Note: This subroutine is used to find the gradient of hydropotential
    !       in Glide's basal water routing scheme.
    ! TODO: Write a parallel version of this subroutine?

      use parallel
        implicit none
        real(dp), dimension(:, :), intent(in) :: f
        real(dp), intent(in) :: deltax, deltay
        real(dp), dimension(:, :), intent(out) :: out_dfdx, out_dfdy
        real(dp), dimension(:, :), intent(in), optional  :: direction_x, direction_y
        
        logical :: upwind !Whether or not directions for upwinding were provided

        integer :: grad_x, grad_y !Whether to upwind or downwind at the current point

        integer :: nx, ny, x, y

        !Get the size of the field we're working with
        nx = size(f, 1)
        ny = size(f, 2)
        
        upwind = present(direction_x) .and. present(direction_y)

        !For now, we'll use the function calls defined above.
        !Later on we might want to refactor?

!LOOP: all scalar points (uses upwinding and downwinding to avoid stepping out of bounds)
        do x = 1, nx
            do y = 1, ny
                grad_x = 0
                grad_y = 0
                if (upwind) then
                    if (direction_x(x,y) < 0.d0 .and. x > 2) then ! Upstream case
                        grad_x = -1
                    else if(direction_x(x,y) > 0.d0 .and. x < nx - 1) then ! Downstream case
                        grad_x = 1
                    end if

                    if (direction_y(x,y) < 0.d0 .and. y > 2) then !Upstream case
                        grad_y = -1
                    else if(direction_y(x,y) > 0.d0 .and. y < ny - 1) then !Downstream case
                        grad_y = 1
                    end if
                end if
  
                !For each of the variables in x, y, check whether or not
                !we need to use an upwinding or downwinding differentiation
                !scheme.
                if (x == 1 .or. grad_x > 0) then
                    out_dfdx(x, y) = dfdx_2d_downwind(f, x, y, deltax)
                else if (x == nx .or. grad_x < 0) then
                    out_dfdx(x, y) = dfdx_2d_upwind(f, x, y, deltax)
                else
                    out_dfdx(x, y) = dfdx_2d(f, x, y, deltax)
                end if
                        
                if (y == 1 .or. grad_y > 0) then
                    out_dfdy(x, y) = dfdy_2d_downwind(f, x, y, deltay)
                elseif (y == ny .or. grad_y < 0) then
                    out_dfdy(x, y) = dfdy_2d_upwind(f, x, y, deltay)
                else
                    out_dfdy(x, y) = dfdy_2d(f, x, y, deltay)
                end if
                        
            end do  
        end do

!NOTE:  If halo updates are needed, they should be done at a higher level.

!!        call parallel_halo(out_dfdx)
!!        call parallel_halo(out_dfdy)
        
    end subroutine df_field_2d

!----------------------------------------------------------------------------

    !> Computes derivative fields of the given function.  Places the result
    !> on a staggered grid.  If periodic in one dimension is set, that 
    !> dimension for derivatives must be the same size as the value's dimension.
    !> Otherwise, it should be one less

    subroutine df_field_2d_staggered(f,                  &
                                     deltax,   deltay,   &
                                     out_dfdx, out_dfdy, &
                                     thck,     thklim )

        implicit none
        real(dp), dimension(:, :), intent(in) :: f, thck    ! unstaggered grid
        real(dp), intent(in) :: deltax, deltay, thklim
        real(dp), dimension(:, :), intent(out) :: out_dfdx, out_dfdy  ! staggered grid
        
        integer :: nx, ny, x, y
        
        !Get the size of the field we're working with
        nx = size(f, 1)
        ny = size(f, 2)

        ! intialize to zeros
        out_dfdx = 0.d0
        out_dfdy = 0.d0
        
        ! *SFP* old subroutine calls, commented out below but still available, 
        ! use centered diffs on normal thck / surf grids but do nothing special at lateral
        ! boundaries where centered diffs might give unreasonable values (e.g., due to jumping
        ! from a region of non-zero to zero thickness / elevation). New calls access new 
        ! subroutines which attempt to correct for this if/when possible using approx., first-order
        ! accurate one-sided diffs.

        !Note - Can remove thck and thklim from argument list if not using new calls

        do x = 1, nx - 1  ! We go to nx - 1 because we're using a staggered grid
            do y = 1, ny - 1
                out_dfdx(x,y) = dfdx_2d_stag(f, x, y, deltax) !*SFP* old call
                out_dfdy(x,y) = dfdy_2d_stag(f, x, y, deltay) !*SFP* old call
!                out_dfdx(x,y) = dfdx_2d_stag_os(f, x, y, deltax, thck, thklim )
!                out_dfdy(x,y) = dfdy_2d_stag_os(f, x, y, deltay, thck, thklim )
            end do
        end do

        end subroutine df_field_2d_staggered

!----------------------------------------------------------------------------

    !> Computes derivative with respect to x at a given point.
    !> Applies periodic boundary conditions if needed.

    function dfdx_2d(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_2d
        
        dfdx_2d = (-.5d0/delta)*f(i-1, j) + (.5d0/delta)*f(i+1, j)
        !write(*,*), i, j, f(i,j), ip1, im1, delta, dfdx_2d
    end function dfdx_2d

!----------------------------------------------------------------------------

    !> Computes derivative with respect to y at a given point

    function dfdy_2d(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_2d
        
        integer :: jp1, jm1
        jp1 = j + 1
        jm1 = j - 1
        if (jp1 == size(f, 2)+1) jp1 = 2
        if (jm1 == 0) jm1 = size(f, 2)-1
        
        dfdy_2d = (-.5d0/delta)*f(i, j-1) + (.5d0/delta)*f(i, j+1)
    end function dfdy_2d

!----------------------------------------------------------------------------

    !> Computes derivative with respect to x at the equivalent
    !> point on a staggered grid.

    function dfdx_2d_stag(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_2d_stag
        dfdx_2d_stag = (f(i+1, j) + f(i+1, j+1) - f(i, j) - f(i, j+1))/(2.d0*delta) 
    end function dfdx_2d_stag

!----------------------------------------------------------------------------

    !> Computes derivative with respect to y at the equivalent
    !> point on a staggered grid.

    function dfdy_2d_stag(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_2d_stag
        dfdy_2d_stag = (f(i, j+1) + f(i+1, j+1) - f(i,j) - f(i+1, j))/(2.d0*delta)
    end function dfdy_2d_stag

!----------------------------------------------------------------------------

    !> Computes derivative with respect to x at the given point
    !> using an upwind method (suitable for maximum boundaries)

    function dfdx_2d_upwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_2d_upwind
        dfdx_2d_upwind = (.5d0 * f(i-2,j) - 2.d0 * f(i-1, j) + 1.5d0 * f(i, j))/delta
    end function dfdx_2d_upwind

!----------------------------------------------------------------------------

    !> Computes derivative with respect to y at the given point
    !> using an upwind method (suitable for maximum boundaries)

    function dfdy_2d_upwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_2d_upwind
        dfdy_2d_upwind = (.5d0 * f(i,j-2) - 2.d0 * f(i, j-1) + 1.5d0 * f(i, j))/delta
    end function dfdy_2d_upwind

!----------------------------------------------------------------------------

    !> Computes derivative with respect to x at the given point
    !> using a downwind method (suitable for minimum boundaries)

    function dfdx_2d_downwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_2d_downwind
        dfdx_2d_downwind = (-1.5d0 * f(i, j) + 2.d0 * f(i+1, j) - .5d0 * f(i+2, j))/delta
    end function dfdx_2d_downwind 

!----------------------------------------------------------------------------

    !> Computes derivative with respect to y at the given point
    !> using a downwind method (suitable for minimum boundaries)

    function dfdy_2d_downwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_2d_downwind
        dfdy_2d_downwind = (-1.5d0 * f(i, j) + 2.d0 * f(i, j+1) - .5d0 * f(i, j+2))/delta
    end function dfdy_2d_downwind

!----------------------------------------------------------------------------

    !------------------------------------------------------------------
    !First Derivative Estimates, Second Order, 3D
    !------------------------------------------------------------------

    !> Computes derivative with respect to x at a given point

    function dfdx_3d(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_3d
        dfdx_3d = (-.5d0/delta)*f(k, i-1, j)  + (.5d0/delta)*f(k, i+1, j)
    end function dfdx_3d

!----------------------------------------------------------------------------

    !> Computes derivative with respect to y at a given point

    function dfdy_3d(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_3d
        dfdy_3d = (-.5d0/delta)*f(k, i, j-1) + (.5d0/delta)*f(k, i, j+1)
    end function dfdy_3d
    
!----------------------------------------------------------------------------

    !> Computes derivative with respect to z at a given point
    !> where the Z axis uses an irregular grid defined by \ittext{deltas}.
    !> This derivative is given by the formula:

    function dfdz_3d_irregular(f, i, j, k, dz)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), dimension(:), intent(in) :: dz
        real(dp) :: dfdz_3d_irregular

        dfdz_3d_irregular = f(k-1,i,j)*(dz(k) - dz(k+1))/((dz(k) - dz(k-1))*(dz(k+1)-dz(k-1))) + &
                            f(k,  i,j)*(dz(k+1)-2.d0*dz(k)+dz(k-1))/((dz(k)-dz(k-1))*(dz(k+1)-dz(k))) + &
                            f(k+1,i,j)*(dz(k)-dz(k-1))/((dz(k+1)-dz(k))*(dz(K+1)-dz(k-1)))
    end function
    
!----------------------------------------------------------------------------

    !> Computes derivative with respect to z at a given point using an upwinding
    !> scheme.  The Z axis uses an irregular grid defined by \iittext{deltas}.

    function dfdz_3d_upwind_irregular(f, i, j, k, deltas)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), dimension(:), intent(in) :: deltas
        real(dp) :: dfdz_3d_upwind_irregular
        real(dp) :: zkMinusZkm1, zkMinusZkm2, zkm1MinusZkm2
        zkMinusZkm1 = deltas(k) - deltas(k-1)
        zkMinusZkm2 = deltas(k) - deltas(k-2)
        zkm1MinusZkm2 = deltas(k-1) - deltas(k-2)
        
        dfdz_3d_upwind_irregular = f(k-2, i, j) * zkMinusZkm1 / (zkm1MinusZkm2 * zkMinusZkm2) - &
                                   f(k-1, i, j) * zkMinusZkm2 / (zkMinusZkm1 * zkm1MinusZkm2) + &
                                   f(k,   i, j) * (2.d0*deltas(k) - deltas(k-1) - deltas(k-2)) / (zkMinusZkm1 * zkMinusZkm2)
    end function
    
!----------------------------------------------------------------------------

    !> Computes derivative with respect to z at a given point using a downwinding
    !> scheme.  The Z axis uses an irregular grid defined by \iittext{deltas}.

    function dfdz_3d_downwind_irregular(f, i, j, k, deltas)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), dimension(:), intent(in) :: deltas
        real(dp) :: dfdz_3d_downwind_irregular
        real(dp) :: zkp1MinusZk, zkp2MinusZk, zkp2MinusZkp1
        zkp1MinusZk = deltas(k+1) - deltas(k)
        zkp2MinusZk = deltas(k+2) - deltas(k)
        zkp2MinusZkp1 = deltas(k+2) - deltas(k+1)
        
        dfdz_3d_downwind_irregular =f(k,   i, j) * (-zkp1MinusZk - zkp2MinusZk)/(zkp1MinusZk * zkp2MinusZk) + &
                                    f(k+1, i, j) * zkp2MinusZk / (zkp2MinusZkp1 * zkp1MinusZk) - &
                                    f(k+2, i, j) * zkp1MinusZk / (zkp2MinusZkp1 * zkp2MinusZk)
    end function
    
!----------------------------------------------------------------------------

    !> Computes derivative with respect to x at the equivalent
    !> point on a staggered grid.

    function dfdx_3d_stag(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_3d_stag
        dfdx_3d_stag = (f(k, i+1, j) + f(k, i+1, j+1) - f(k, i, j) - f(k, i, j+1))/(2.d0*delta) 
    end function dfdx_3d_stag

!----------------------------------------------------------------------------

    !> Computes derivative with respect to y at the equivalent
    !> point on a staggered grid.

    function dfdy_3d_stag(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_3d_stag
        dfdy_3d_stag = (f(k, i, j+1) + f(k, i+1, j+1) - f(k, i, j) - f(k, i+1, j))/(2.d0*delta)
    end function dfdy_3d_stag

!----------------------------------------------------------------------------

    !> Computes derivative with respect to x at the given point
    !> using an upwind method (suitable for maximum boundaries)

    function dfdx_3d_upwind(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_3d_upwind
        dfdx_3d_upwind = (.5d0 * f(k, i-2, j) - 2.d0 * f(k, i-1, j) + 1.5d0 * f(k, i, j))/delta
    end function dfdx_3d_upwind

!----------------------------------------------------------------------------

    !> Computes derivative with respect to y at the given point
    !> using an upwind method (suitable for maximum boundaries)

    function dfdy_3d_upwind(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_3d_upwind
        dfdy_3d_upwind = (.5d0 * f(k, i, j-2) - 2.d0 * f(k, i, j-1) + 1.5d0 * f(k, i, j))/delta
    end function dfdy_3d_upwind

!----------------------------------------------------------------------------

    !> Computes derivative with respect to x at the given point
    !> using a downwind method (suitable for minimum boundaries)

    function dfdx_3d_downwind(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j, k
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_3d_downwind
        dfdx_3d_downwind = (-1.5d0 * f(k, i, j) + 2.d0 * f(k, i+1, j) - .5d0 * f(k, i+2, j))/delta
    end function dfdx_3d_downwind 

!----------------------------------------------------------------------------

    !> Computes derivative with respect to y at the given point
    !> using a downwind method (suitable for minimum boundaries)

    function dfdy_3d_downwind(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_3d_downwind
        dfdy_3d_downwind = (-1.5d0 * f(k, i, j) + 2.d0 * f(k, i, j+1) - .5d0 * f(k, i, j+2))/delta
    end function dfdy_3d_downwind

!----------------------------------------------------------------------------

end module glide_grid_operators

!----------------------------------------------------------------------------
