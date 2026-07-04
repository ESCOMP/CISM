!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_plume.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains subroutines for a sub-ice-shelf plume model.
!
! Author: William Lipscomb
!         NSF National Center for Atmospheric Research
!         Climate and Global Dynamics Laboratory
!         Boulder, CO 80303
!         USA
!         <lipscomb@ucar.edu>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
  module glissade_plume

    use glimmer_global, only: dp
    use glimmer_physcon, only: rhoi, rhow, lhci, cpw, pi
    use glimmer_paramets, only: iulog
    use glimmer_log
    use glimmer_utils, only: point_diag
    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo, lhalo, uhalo, &
         parallel_halo, parallel_reduce_max, parallel_global_sum, parallel_globalindex

    implicit none
    save
    private

    public :: glissade_plume_driver, verbose_plume

    logical, parameter :: verbose_plume = .true.

    ! prescribed MISOMIP parameters from Table 4 of Asay-Davis et al. (2016)
    ! Note: cpw and lhci are defined in glimmer_physcon
    !TODO - Add these to a derived type or a constants module?
    real(dp), parameter :: &
         lambda1 = -0.0573d0,        & ! liquidus slope (deg/psu)
         lambda2 =  0.0832d0,        & ! liquidus intercept (deg C)
         lambda3 = -7.53d-8,         & ! liquidus pressure coefficient (deg/Pa)
                                       ! Tb = lambda1*Sb + lambda2 + lambda3*pb
         c_drag = 2.5d-3,            & ! ocean drag coefficient (unitless)
         u_tidal = 0.01d0,           & ! tidal velocity (m/s)
         eos_rho_ref = 1027.51d0,    & ! reference density for linear EOS (kg/m^3)
         eos_Tref = -1.0d0,          & ! reference temperature for linear EOS (deg C)
         eos_Sref = 34.2d0,          & ! reference salinity for linear EOS (deg C)
         eos_alpha = 3.733d-5,       & ! thermal expansion coefficient for linear EOS (deg^-1)
         eos_beta = 7.843d-4,        & ! salinity contraction coefficient for linear EOS (psu^-1)
         ! Note: eos_Tref = -1 C and eos_Sref = 34.2 psu are not used
         f_coriolis = -1.405d-4        ! Coriolis parameter (s^-1) at 75 S = 2*omega*sin(75 deg) (prescribed in text)

!=======================================================================

  contains

!=======================================================================

  subroutine glissade_plume_driver(model)

    type(glide_global_type)  :: model

    ! local variables


  end subroutine glissade_plume_driver

!****************************************************

  subroutine compute_plume_velocity(&
       nx,    ny,               &
       itest, jtest, rtest,     &
       edge_mask,               &
       D_plume,                 &
       pgf_x,                   &
       pgf_y,                   &
       latdrag_x,               &
       latdrag_y,               &
       u_plume,                 &
       v_plume,                 &
       converged_velo,          &
       edge_mask_east_reduce_v, &
       edge_mask_north_reduce_u)
    
    ! Compute the velocity on a set of edges (either east or north)

    integer, intent(in) ::  &
         nx,  ny,           & ! number of grid cells in each dimension
         itest, jtest, rtest  ! test cell coordinates (diagnostic only)
    
    ! Used to be intent(in), but now are module variables
!    real(dp), intent(in) ::   &
!         u_tidal,           & ! tidal velocity (m/s)
!         c_drag,            & ! ocean drag coefficient (unitless)   
!         f_coriolis           ! Coriolis parameter (s^-1)
    
    integer, dimension(nx,ny), intent(in) ::   &
         edge_mask            ! = 1 at edges where velocity is computed

    ! Note: The following variables are co-located with the velocity
    real(dp), dimension(nx,ny), intent(in) ::   &
         D_plume,           & ! plume thickness at edges (m)
         pgf_x,             & ! x component of pressure gradient force
         pgf_y,             & ! y component of pressure gradient force
         latdrag_x,         & ! x component of lateral drag
         latdrag_y            ! y component of lateral drag
    
    real(dp), dimension(nx,ny), intent(inout) ::  &
         u_plume,           & ! x component of plume velocity (m/s)
         v_plume              ! x component of plume velocity (m/s)

    logical, dimension(nx,ny), intent(inout) ::  &
         converged_velo        ! true when velocity has converged at an edge, else false

    !TODO - Remove these terms if lateral drag works
    real(dp), dimension(nx,ny), intent(in), optional :: &
         edge_mask_east_reduce_v,  & ! mask for reducing v on east edges adjacent to a wall
         edge_mask_north_reduce_u    ! mask for reducing u on north edges adjacent to a wall

    ! local variables

    real(dp), dimension(nx,ny) ::   &
         f_x,               &  ! pgf_x + latdrag_x 
         f_y                   ! pgf_y + latdrag_y 

    real(dp), dimension(nx,ny) ::  &
         reduce_v,          &  ! local version of edge_mask_east_reduce_v; no reduction by default
         reduce_u              ! local version of edge_mask_north_reduce_u; no reduction by default

    real(dp) :: &
         plume_speed,       & ! plume speed (m/s)
         x_resid, y_resid,  & ! residuals of momentum balance equations (m^2/s^2)
         denom,             & ! denominator
         a_uu, a_uv,        & ! coefficients for Newton solve
         a_vu, a_vv,        & !
         du, dv               ! change in u_plume and v_plume (m/s)
    
    character(len=128) :: message

    real(dp), parameter :: &
         maxresid_force_balance = 1.0d-8 ! max residual allowed in momentum balance equation (m^2/s^2)
    
    logical, parameter :: &
         velo_newton = .true.  ! if true, use Newton's method; if false, use Picard method

    integer :: i, j

    !TODO - Add lateral drag to the equations
    !       Can be handled numerically by combining with pgf in a single force term

    !--------------------------------------------------------------------
    ! Compute the plume velocity.
    ! Assume a balance between the pressure gradient force, basal drag and Coriolis:
    !
    ! pgf_x - c_d*|U|*u + D*f*v = 0
    ! pgf_y - c_d*|U|*v - D*f*u = 0
    !
    !  where pgf_x = g' * D * db/dx (m^2/s^2) 
    !        pgf_y = g' * D * db/dy (m^2/s^2) 
    !            D = plume boundary-layer thickness
    !           g' = reduced gravity = g*(rhoa - rhop)/rhoo
    !         rhoa = ambient ocean density
    !         rhop = plume density
    !         rhoo = reference ocean density
    !            b = elevation of shelf base
    !          c_d = dimensionless ocean drag coefficient
    !            f = Coriolis coefficient
    !          |U| = sqrt(u^2 + v^2 + u_tidal^2)
    !      u_tidal = a small velocity added for regularization
    !
    ! The solution (assuming D is known) is
    !
    !                c_d*|U|*pgf_x + D*f*pgf_y
    !            u = ________________________
    !                 (D*f)^2 + (c_d*|U|)^2
    !
    !                c_d*|U|*pgf_y - D*f*pgf_x 
    !            v = ________________________
    !                 (D*f)^2 + (c_d*|U|)^2
    !
    ! Since |U| is a function of u and v, we iterate to convergence.
    !
    ! The iteration is sped up by using Newton's method.
    ! We write   u = u0 + du
    !            v = v0 + dv
    !          |U| = U0 + d|U|/du * du + d|U|dv * dv
    ! where the partial derivatives are evaluated at (u,v) = (u0,v0).
    ! 
    ! This gives 
    !           du = (a_vv * R_x - a_uv * R_y) / det|A|
    !           dv = (a_uu * R_y - a_vu * R_x) / det|A|
    ! where    
    !          R_x = pgf_x - c_d*U0*u0 + D*f*v0 = x residual
    !          R_y = pgf_y - c_d*U0*v0 - D*f*u0 = y residual
    !
    !                | a_uu   a_uv |
    ! and        A = |             |     
    !                | a_vu   a_vv |
    !
    ! with    a_uu = c_d*(U0 + u0^2/U0)
    !         a_uv = c_d*u0*v0/U0 - D*f) 
    !         a_vu = c_d*u0*v0/U0 + D*f) 
    !         a_vv = c_d*(U0 + v0^2/U0) 
    !
    ! If reduce_u < 1 or reduce_v < 1, then the Coriolis term in these equations
    ! is reduced proportionately, so as to inhibit flow into walls.
    !
    !--------------------------------------------------------------------

    if (present(edge_mask_north_reduce_u)) then
       reduce_u(:,:) = edge_mask_north_reduce_u(:,:)
    else
       reduce_u(:,:) = 1.0d0  ! no reduction
    endif

    if (present(edge_mask_east_reduce_v)) then
       reduce_v(:,:) = edge_mask_east_reduce_v(:,:)
    else
       reduce_v(:,:) = 1.0d0  ! no reduction
    endif

    ! Combine PGF and lateral drag into one term
    f_x(:,:) = pgf_x(:,:) + latdrag_x(:,:)
    f_y(:,:) = pgf_y(:,:) + latdrag_y(:,:)

    ! Loop over edges of locally owned cells
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          if (edge_mask(i,j) == 1 .and. .not.converged_velo(i,j) ) then
       
             ! Compute plume speed based on current u and v
             plume_speed = sqrt(u_plume(i,j)**2 + v_plume(i,j)**2 + u_tidal**2)
       
             ! Compute residual of the momentum balance
!               x_resid = pgf_x - c_drag*plume_speed*u_plume + f_coriolis*D_plume*v_plume
!               y_resid = pgf_y - c_drag*plume_speed*v_plume - f_coriolis*D_plume*u_plume
             x_resid = f_x(i,j) - c_drag*plume_speed*u_plume(i,j) + reduce_v(i,j)*f_coriolis*D_plume(i,j)*v_plume(i,j)
             y_resid = f_y(i,j) - c_drag*plume_speed*v_plume(i,j) - reduce_u(i,j)*f_coriolis*D_plume(i,j)*u_plume(i,j)

             ! check convergence of plume velocity

             if (abs(x_resid) < maxresid_force_balance .and. abs(y_resid) < maxresid_force_balance) then

                converged_velo(i,j) = .true.

                ! diagnostic print
                if (this_rank == rtest .and. i==itest .and. j==jtest) then
                   write(iulog,*) ' '
                   write(iulog,*) 'Velocity converged: u/v_plume (m/s):', u_plume(i,j), v_plume(i,j)
                endif

             endif

             if (.not.converged_velo(i,j)) then

                if (velo_newton) then
          
                   ! compute some coefficients for the Newton solve
                   a_uu = c_drag * (plume_speed + u_plume(i,j)**2/plume_speed)
                   a_vv = c_drag * (plume_speed + v_plume(i,j)**2/plume_speed)
                      
                   a_uv = c_drag * (u_plume(i,j)*v_plume(i,j))/plume_speed - reduce_v(i,j)*D_plume(i,j)*f_coriolis
                   a_vu = c_drag * (u_plume(i,j)*v_plume(i,j))/plume_speed + reduce_u(i,j)*D_plume(i,j)*f_coriolis
                   
                   ! compute du and dv
                   denom = a_uu*a_vv - a_uv*a_vu
                      
                   if (abs(denom) > 0.0d0) then
                      du = (a_vv*x_resid - a_uv*y_resid) / denom
                      dv = (a_uu*y_resid - a_vu*x_resid) / denom
                         
                      u_plume(i,j) = u_plume(i,j) + du
                      v_plume(i,j) = v_plume(i,j) + dv
                      
                   else  ! denom = 0.0
                      write(iulog,*) 'Error, glissade_plume: ill-posed Newton solve for velocity, rank, i, j:', this_rank, i, j
                      write(iulog,*) 'a_uu, a_vv, a_uv, a_vu =', a_uu, a_vv, a_uv, a_vu
                      write(message,*) 'Error, glissade_plume: ill-posed Newton solve for velocity, rank, i, j:', this_rank, i, j
                      call write_log(message, GM_FATAL)
                   endif
                      
                else  ! simpler Picard solve
          
                   denom = (c_drag*plume_speed)**2 + (D_plume(i,j)*f_coriolis)**2
                   u_plume = (c_drag*plume_speed*f_x(i,j) + reduce_v(i,j)*D_plume(i,j)*f_coriolis*f_y(i,j)) / denom
                   v_plume = (c_drag*plume_speed*f_y(i,j) - reduce_u(i,j)*D_plume(i,j)*f_coriolis*f_x(i,j)) / denom
          
                endif  ! Newton or Picard

             endif  ! .not.converged_velo

             if (verbose_plume .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'plume_speed (m/s) =', plume_speed
                write(iulog,*) 'pgf_x, pgf_y:', pgf_x(i,j), pgf_y(i,j)
                write(iulog,*) 'latdrag_x, latdrag_y:', latdrag_x(i,j), latdrag_y(i,j)
                write(iulog,*) 'Dfv, -Dfu:', D_plume(i,j) * f_coriolis * v_plume(i,j), &
                                     -D_plume(i,j) * f_coriolis * u_plume(i,j)
                write(iulog,*) 'dragu, dragv:', c_drag * plume_speed * u_plume(i,j), &
                                         c_drag * plume_speed * v_plume(i,j)
                write(iulog,*) 'x/y residual:', x_resid, y_resid
                write(iulog,*) 'new u/v_plume:', u_plume(i,j), v_plume(i,j)
             endif

          endif  ! edge_mask
       enddo  ! i
    enddo  ! j
    
  end subroutine compute_plume_velocity

!****************************************************

  subroutine compute_plume_melt_rate(&
       nx,         ny,      &
       gammaT,              &
       gammaS,              &
       plume_mask_cell,     &
       pressure,            &
       entrainment,         &
       u_plume_east,        &
       v_plume_north,       &
       T_ambient,           &
       S_ambient,           &
       T_basal,             &
       S_basal,             &
       T_plume,             &
       S_plume,             &
       itest, jtest, rtest, &
       ustar_plume,         &
       bmlt_float)
    
    !--------------------------------------------------------------------
    ! Compute the melt rate at the ice-ocean interface.
    !
    ! There are 5 equations for 5 unknowns: m, Tb, Sb, T and S
    ! where m = melt rate at ice-ocean interface
    !       Tb = potential temperature at ice-ocean interface
    !       Sb = salinity at ice-ocean interface
    !       T  = potential temperature of boundary-layer plume
    !       S  = salinity of boundary-layer plume
    ! 
    ! (1) rhow * m * L  = rhoo * cw * u_fric * gammaT * (T - Tb)
    ! (2) rhow * m * Sb = rhoo * u_fric * gammaS *(S - Sb)
    ! (3) Tb = lambda1*Sb + lambda2 + lambda3*pb 
    ! (4) L * m = -cw * e * (T - Ta)
    ! (5) S * m = -e * (S - Sa)
    !
    ! Eq. 1 and 2 describe heat and salt transfer at the ice-ocean interface.
    ! Eq. 3 is the linearized liquidus relation that determines the potential freezing point.
    ! Eq. 4 and 5 describe heat and salt entrainment from the ambient ocean to the boundary-layer plume,
    !  where Ta and Sa are the potential temperature and salinity of the ambient ocean.
    !
    ! We can rewrite (1) and (2) as
    !
    ! (1)     m = T_factor * (T - Tb)
    ! (2)  Sb*m = S_factor * (S - Sb) 
    !
    ! where T_factor = (rhoo * cw * ufric * gammaT) / (rhow * L)
    !       S_factor = (rhoo * ufric * gammaS) / rhow
    !
    ! Rearrange (4):  T = Ta - (L/(cw*e)) * m
    ! 
    ! Use (3) and (4) to replace T and Tb in (1):
    !
    ! (1')    m = m1*Sb + m2
    ! where  m1 = -T_factor*lambda1/denom
    !        m2 =  T_factor*(Ta - lambda2 - lambda3*p)/denom
    !     denom = 1 + T_factor*L/(cw*e)
    ! 
    ! Use (5) to replace S in (2):
    !
    ! (2')   Sb = S_factor*e*Sa / ((m+S_factor)*(m+e))
    !
    ! Use (2') to replace Sb in (1') to form a cubic equation for m:
    !
    ! (1'')  a*m^3 + b*m^2 + c*m + d = 0
    !
    !   where a = 1
    !         b = S_factor + e - m2
    !         c = S_factor*e - m2*(S_factor + e)
    !         d = -S_factor*e*(m1*Sa + m2)
    !
    ! Use the cubic_solver subroutine to find m.
    !
    ! Given m, back out the other 4 unknowns.
    !--------------------------------------------------------------------
    
    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    integer, intent(in) ::  &
         itest, jtest, rtest    ! test cell coordinates (diagnostic only)

    ! Note: gammaS and gammaT are config parameters and are passed in as arguments.
    !       Other MISOMIP parameters are declared at the top of the module.
    
    real(dp), intent(in) ::  &
         gammaT,              & ! nondimensional heat transfer coefficient
         gammaS                 ! nondimensional salt transfer coefficient
    
    integer, dimension(nx,ny), intent(in) :: &
         plume_mask_cell        ! = 1 for cells where scalar plume variables are computed

    real(dp), dimension(nx,ny), intent(in) :: &
         pressure,            & ! ocean pressure at base of ice (N/m^2)
         entrainment,         & ! entrainment rate of ambient water into plume (m/s)
         u_plume_east,        & ! u_plume on east edges (m/s)
         v_plume_north,       & ! v_plume on north edges (m/s)
         T_ambient,           & ! ambient ocean potential temperature at depth of ice-ocean interface (deg C)
         S_ambient              ! ambient ocean salinity at depth of ice-ocean interface (psu)
    
    real(dp), dimension(nx,ny), intent(out) :: &
         ustar_plume,         & ! plume friction velocity (m/s) on ice grid, output as a diagnostic
         T_basal,             & ! basal ice temperature (deg C)
         S_basal,             & ! basal ice salinity (psu)
         T_plume,             & ! plume temperature (deg C)
         S_plume,             & ! plume salinity (psu)
         bmlt_float             ! melt rate at base of floating ice (m/s)
    
    ! local variables
    
    real(dp) :: &
         u_plume, v_plume,    & ! plume velocity components at cell center (m/s)
         plume_speed,         & ! plume speed at cell center (m/s)
         T_factor, S_factor,  & ! factors in melt-rate equations
         denom,               & ! denominator
         m1, m2,              & ! factors in relation between m and Sb
         ma, mb, mc, md,      & ! coefficients in cubic equation for m
         bmlt_float_avg         ! average value of bmlt_float in main cavity
    
    integer :: i, j
    
    !WHL - debug -  Test cubic solver
!       ma =    2.d0
!       mb =  -30.d0
!       mc =  162.d0
!       md = -350.d0
!       call cubic_solver(ma, mb, mc, md, solution)
!       write(iulog,*) 'Trial cubic solution =', solution
!       write(iulog,*) 'True solution =', (10.d0 + sqrt(108.d0))**(1.d0/3.d0) - (-10.d0 + sqrt(108.d0))**(1.d0/3.d0) + 5.d0


    ! Loop over locally owned cells
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          
          if (plume_mask_cell(i,j) == 1 .and. entrainment(i,j) > 0.0d0) then
             
             ! Interpolate the plume speed to the cell center, and compute the friction velocity ustar.
             
             u_plume = (u_plume_east(i,j) + u_plume_east(i-1,j)) / 2.0d0
             v_plume = (v_plume_north(i,j) + v_plume_north(i,j-1)) / 2.0d0
             plume_speed = sqrt(u_plume**2 + v_plume**2 + u_tidal**2)
             ustar_plume(i,j) = sqrt(c_drag) * plume_speed

             T_factor = (rhoo * cpw * ustar_plume(i,j) * gammaT) / (rhow * lhci)
             S_factor = (rhoo * ustar_plume(i,j) * gammaS) / rhow
             
             denom = 1.d0 + (T_factor*lhci)/(cpw*entrainment(i,j))
             m1 = -lambda1 * T_factor / denom
             m2 = T_factor * (T_ambient(i,j) - lambda2 - lambda3*pressure(i,j)) / denom
             
             ma = 1.d0
             mb = S_factor + entrainment(i,j) - m2
             mc = S_factor*entrainment(i,j) - m2*(S_factor + entrainment(i,j))
             md = -S_factor*entrainment(i,j)*(m1*S_ambient(i,j) + m2)
             
             ! Solve the cubic equation
             call cubic_solver(&
                  ma, mb, mc, md, &
                  bmlt_float(i,j))

             if (verbose_plume .and. this_rank == rtest .and. i==itest .and. j==jtest) then
                write(iulog,*) ' '
                write(iulog,*) 'Melt rate calc: rank, i, j =', rtest, i, j
                write(iulog,*) 'pressure (Pa) =', pressure(i,j)
                write(iulog,*) 'T_factor (m/s/deg), S_factor (m/s)=', T_factor, S_factor
                write(iulog,*) 'entrainment (m/s) =', entrainment(i,j)
                write(iulog,*) 'm1 (m/s/psu) =', m1
                write(iulog,*) 'm2 (m/s) =', m2
                write(iulog,*) 'denom =', denom
                write(iulog,*) 'a, b, c, d =', ma, mb, mc, md
                write(iulog,*) 'residual of cubic solve =', ma*bmlt_float(i,j)**3 + mb*bmlt_float(i,j)**2 + mc*bmlt_float(i,j) + md
             endif
             
             ! Given the melt rate, compute Sb and Tb
!               S_basal(i,j) = (S_factor * entrainment(i,j) * S_ambient(i,j)) /  &
!                               ( (bmlt_float(i,j) + S_factor) * (bmlt_float(i,j) + entrainment(i,j)) )
             S_basal(i,j) = (bmlt_float(i,j) - m2) / m1
             T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)

             ! Given m, compute S and T for plume
             T_plume(i,j) = T_ambient(i,j) - (lhci/(cpw*entrainment(i,j))) * bmlt_float(i,j)
             S_plume(i,j) = S_ambient(i,j) * entrainment(i,j) / (bmlt_float(i,j) + entrainment(i,j))

             !WHL - debug - check for NaNs
             if (T_plume(i,j) /= T_plume(i,j) .or. S_plume(i,j) /= S_plume(i,j) .or. &
                 T_basal(i,j) /= T_basal(i,j) .or. S_basal(i,j) /= S_basal(i,j) .or. &
                 bmlt_float(i,j) /= bmlt_float(i,j)) then
                write(iulog,*) 'Bad values, i, j =', i, j
                write(iulog,*) 'T_plume, S_plume:', T_plume(i,j), S_plume(i,j)
                write(iulog,*) 'T_basal, S_basal:', T_basal(i,j), S_basal(i,j)
                write(iulog,*) 'bmlt_float:', bmlt_float(i,j)
                stop
             endif

          else    ! plume_mask_cell = 0
             
             bmlt_float(i,j) = 0.0d0
             
             S_plume(i,j) = S_ambient(i,j)
             T_plume(i,j) = T_ambient(i,j)
             
             S_basal(i,j) = S_ambient(i,j)
             T_basal(i,j) = lambda1*S_basal(i,j) + lambda2 + lambda3*pressure(i,j)
             
          endif   ! plume_mask_cell and entrainment > 0
          
       enddo   ! i
    enddo   ! j

  end subroutine compute_plume_melt_rate

!****************************************************
    
  !TODO - Move this subroutine to a utility module?
  !TODO - Pass 3 complex roots in and out.
  subroutine cubic_solver(&
       a, b, c, d, &
       x1,         &
       x2_r, x2_i, &
       x3_r, x3_i)

    !------------------------------------------------
    ! Find the real root of a cubic equation:
    !
    !    ax^3 + bx^2 + cx = d = 0
    !
    ! Do this by making the substitution
    !
    !    x = y - b/(3a)
    !
    ! to convert to a depressed cubic:
    !
    !    y^3 + py + q = 0
    !
    ! where p = (1/a) * (c - b^2/(3a))
    !       q = (1/a) * (d + 2b^3/(27a^2) - bc/(3a))
    !
    !------------------------------------------------

    real(dp), intent(in) ::  &
         a, b, c, d       ! coefficients of cubic equation
                          ! assumed to be real

    real(dp), intent(out) ::  &
         x1               ! real solution of cubic equation

    real(dp), intent(out), optional ::  &
         x2_r, x2_i,    & ! other solutions of cubic equation
         x3_r, x3_i       ! could be either real or complex

    real(dp) :: &
         p, q             ! coefficients of depressed cubic

    real(dp) :: &
         Delta            ! discriminant

    real(dp) :: &
         y1,            & ! solutions of depressed cubic
         y2_r, y2_i,    & !
         y3_r, y3_i

    real(dp) :: &
         u, v,          & ! some intermediate factors
         fu, fv,        &
         phi

    real(dp), parameter :: &
         p333 = 1.d0/3.d0

    !WHL - debug
    logical, parameter :: verbose_cubic = .false.

    ! compute coefficients of depressed cubic, y^3 + py + q = 0

    p = (3.d0*c/a - (b/a)**2) / 3.d0
    q = (2.d0*(b/a)**3 - 9.d0*b*c/(a*a) + 27.d0*d/a) / 27.d0

    ! compute the discriminant
    Delta = (p/3.d0)**3 + (q/2.d0)**2

    if (verbose_cubic) then
       write(iulog,*) 'Delta =', Delta
       if (Delta > 0.d0) then
          write(iulog,*) 'One real root, 2 complex conjugate'
       elseif (Delta == 0.d0) then
          write(iulog,*) 'Three real roots of which at least two are equal'
       elseif (Delta < 0.d0) then
          write(iulog,*) 'Three distinct real roots'
       endif
    endif

    if (Delta >= 0.d0) then   

       if (Delta > 0.d0) then    ! one real root, two complex roots
          fu = -q/2.d0 + sqrt(Delta)
          fv = -q/2.d0 - sqrt(Delta)
       else  ! Delta = 0; three real roots of which at least two are equal
          fu = -q/2.d0
          fv = fu
       endif
 
       ! some logic to avoid taking cube roots of negative numbers
       if (fu >= 0.d0) then
          u = fu**p333
       else
          u = -(-fu)**p333
       endif

       if (fv >= 0.d0) then
          v = fv**p333
       else
          v = -(-fv)**p333
       endif

       ! form solutions of depressed cubic
       y1 = u + v       ! real
       y2_r = -(u+v)/2.d0
       y2_i =  (u-v)*sqrt(3.d0)/2.d0
       y3_r = -(u+v)/2.d0
       y3_r = -(u-v)*sqrt(3.d0)/2.d0

       if (verbose_cubic) then
          write(iulog,*) 'a, b, c, d:', a, b, c, d
          write(iulog,*) 'p, q:', p, q
          write(iulog,*) 'y1 =', y1
          write(iulog,*) 'x1 =', x1
       endif

    else  ! Delta < 0; three distinct real roots
          ! use a trigonometric formulation

       phi = acos(-q/(2.d0*sqrt(abs(p)**3/27.d0)))

       y1 =    2.d0 * sqrt(abs(p)/3.d0) * cos(phi/3.d0)
       y2_r = -2.d0 * sqrt(abs(p)/3.d0) * cos((phi+pi)/3.d0)
       y2_i =  0.d0
       y3_r = -2.d0 * sqrt(abs(p)/3.d0) * cos((phi-pi)/3.d0)
       y3_i =  0.d0

       if (verbose_cubic) then
          write(iulog,*) 'a, b, c, d:', a, b, c, d
          write(iulog,*) 'p, q:', p, q
          write(iulog,*) 'y1, y2, y3 =', y1, y2_r, y3_r
          write(iulog,*) 'b/3a =', b/(3.d0*a)
          write(iulog,*) 'x1 =', y1 - b/(3.d0*a)
          write(iulog,*) 'x2 =', y2_r - b/(3.d0*a)
          write(iulog,*) 'x3 =', y3_r - b/(3.d0*a)
       endif

    endif

    ! Recover the solutions
    ! Mostly likely we are only interested in x1, but compute the others if requested

    x1 = y1 - b/(3.d0*a)

    if (present(x2_r) .and. present(x2_i) .and. present(x3_r) .and. present(x3_i)) then
       x2_r = y2_r - b/(3.d0*a)
       x2_i = y2_i
       x3_r = y3_r - b/(3.d0*a)
       x3_i = y3_i
    endif

  end subroutine cubic_solver

!****************************************************

  !TODO - Move this subroutine to glissade_grid_operations?

  subroutine compute_edge_gradients(&
       nx,              ny,          &
       dx,              dy,          &
       global_bndy_east,             &
       global_bndy_west,             &
       global_bndy_north,            &
       global_bndy_south,            &
       plume_mask_cell,              &
       floating_mask,                &
       lsrf,                         &
       field,                        &
       df_dx_east,      df_dy_east,  &
       df_dx_north,     df_dy_north)
   
    ! Compute the gradients of a scalar field on east and north cell edges.
    ! The procedure for east edges as follows:
    ! (1) Initialize all gradients to zero.
    ! (2) If the plume exists on both sides of an east edge, compute df/dx in the standard way.
    !     Similarly, if the plume exists on both sides of a north edge, compute df/dy in the standard way.
    ! (3) If the edge has a plume cell on one side and floating ice or open water on the other,
    !     and it is not a global boundary edge, then extrapolate the gradient from an adjacent edge.
    ! (4) Compute df/dy on east edges by averaging from adjacent north edges, and compute
    !     df/dx on north edges by extrapolating from adjacent east edges.
    
    integer, intent(in) ::  &
         nx,     ny             ! number of grid cells in each dimension

    real(dp), intent(in) ::  &
         dx,     dy             ! grid cell size (m)
    
    integer, dimension(nx,ny), intent(in) ::  &
         global_bndy_east,   & ! = 1 for edges at global boundaries, else = 0
         global_bndy_west,   &
         global_bndy_north,  &
         global_bndy_south,  &
         plume_mask_cell,    & ! = 1 for cells where scalar plume variables are computed
         floating_mask         ! = 1 where ice is present and floating, else = 0
    
    real(dp), dimension(nx,ny), intent(in) ::  &
         lsrf                  ! lower ice surface (m); used to diagnose open ocean
    
    
    real(dp), dimension(nx,ny), intent(in) :: &
         field                 ! scalar field
    
    real(dp), dimension(nx,ny), intent(out) :: &
         df_dx_east,  df_dy_east,   &  ! gradient components on east edges
         df_dx_north, df_dy_north      ! gradient component on north edges
    
    ! local variables

    integer :: i, j

    ! initialize
    df_dx_east(:,:) = 0.0d0
    df_dy_east(:,:) = 0.0d0
   
    df_dx_north(:,:) = 0.0d0
    df_dy_north(:,:) = 0.0d0
   
    ! Compute gradients at edges with plume cells on each side

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! east edges
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 1) then
             df_dx_east(i,j) = (field(i+1,j) - field(i,j)) / dx
          endif

          ! north edges
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 1) then
             df_dy_north(i,j) = (field(i,j+1) - field(i,j)) / dy
          endif

       enddo
    enddo

    ! Set gradients at edges that have a plume cell on one side and floating ice or water on the other.
    ! Extrapolate the gradient from the nearest neighbor edge.
    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! east edges
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i+1,j) == 0 .and. global_bndy_east(i,j) == 0) then
             if (lsrf(i+1,j) == 0.0d0 .or. floating_mask(i+1,j) == 1) then
                df_dx_east(i,j) = df_dx_east(i-1,j)
             endif
          endif
          if (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i+1,j) == 1 .and. global_bndy_west(i,j) == 0) then
             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                df_dx_east(i,j) = df_dx_east(i+1,j)
             endif
          endif

          ! north edges
          if (plume_mask_cell(i,j) == 1 .and. plume_mask_cell(i,j+1) == 0 .and. global_bndy_north(i,j) == 0) then
             if (lsrf(i,j+1) == 0.0d0 .or. floating_mask(i,j+1) == 1) then
                df_dy_north(i,j) = df_dy_north(i,j-1)
             endif
          endif
          if (plume_mask_cell(i,j) == 0 .and. plume_mask_cell(i,j+1) == 1 .and. global_bndy_south(i,j) == 0) then
             if (lsrf(i,j) == 0.0d0 .or. floating_mask(i,j) == 1) then
                df_dy_north(i,j) = df_dy_north(i,j+1)
             endif
          endif

       enddo
    enddo

    ! Average over 4 neighboring edges to estimate the y derivative on east edges and the x derivative on north edges.

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo

          ! y derivative on east edges
          df_dy_east(i,j) = 0.25d0 * (df_dy_north(i,j)   + df_dy_north(i+1,j)  &
                                    + df_dy_north(i,j-1) + df_dy_north(i+1,j-1))

          ! x derivative on north edges
          df_dx_north(i,j) = 0.25d0 * (df_dx_east(i-1,j+1) + df_dx_east(i,j+1)  &
                                     + df_dx_east(i-1,j)   + df_dx_east(i,j))

       enddo
    enddo

    !TODO - Add a halo update for parallel runs

  end subroutine compute_edge_gradients

!****************************************************

  end module glissade_plume

!****************************************************
