!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo.F90 - part of the Community Ice Sheet Model (CISM)  
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

!TODO - Are all these includes needed?
#ifdef HAVE_CONFIG_H 
#include "config.inc" 
#endif

module glissade_velo

    ! Driver for Glissade velocity solvers

    implicit none
    
contains
        
    subroutine glissade_velo_driver(model)

      ! Glissade higher-order velocity driver

      use glimmer_log
      use glimmer_paramets, only: iulog
      use glimmer_physcon, only: scyr
      use glissade_types
      use glissade_velo_higher, only: glissade_velo_higher_solve
      use glissade_velo_sia, only: glissade_velo_sia_solve
      use cism_profile, only: t_startf, t_stopf

      type(glissade_global_type),intent(inout) :: model

      real(dp), dimension(:,:,:), allocatable :: &
           uvel_sia, vvel_sia        ! temporary SIA velocity

      integer :: i, j, k
      integer :: ewn, nsn, upn
      integer :: itest, jtest, rtest
      integer :: whichbtrc_sav

      logical :: verbose_velo = .false.

      ewn = model%general%ewn
      nsn = model%general%nsn
      upn = model%general%upn

      ! get coordinates of diagnostic point
      rtest = -999
      itest = 1
      jtest = 1
      if (this_rank == model%numerics%rdiag_local) then
         rtest = model%numerics%rdiag_local
         itest = model%numerics%idiag_local
         jtest = model%numerics%jdiag_local
      endif

      !-------------------------------------------------------------------
      ! Call the velocity solver.
      ! The standard glissade higher-order solver is glissade_velo_higher_solve.
      ! There is also a local shallow-ice solver, glissade_velo_sia_solve.
      !-------------------------------------------------------------------
      
      if (model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
 
         call glissade_velo_sia_solve(model, ewn, nsn, upn)

      elseif (model%options%which_ho_approx == HO_APPROX_HYBRID) then

         !-------------------------------------------------------------------
         ! compute the SIA part of the velocity, assuming no basal sliding
         !-------------------------------------------------------------------

         ! make sure basal sliding is turned off
         whichbtrc_sav = model%options%whichbtrc
         model%options%whichbtrc = BTRC_ZERO

         call glissade_velo_sia_solve(model, ewn, nsn, upn)

         ! restore the original value of whichbtrc, just in case
         model%options%whichbtrc = whichbtrc_sav

         ! save the result
         allocate(uvel_sia(upn,ewn,nsn))
         allocate(vvel_sia(upn,ewn,nsn))
         uvel_sia = model%velocity%uvel
         vvel_sia = model%velocity%vvel

         !WHL: Do uvel and vvel have dimensions of m/s at this point?
         if (verbose_velo .and. this_rank == rtest) then
            i = itest
            j = jtest
            write(iulog,*) ' '
            write(iulog,*) 'SIA part of uvel, vvel (m/yr): r, i, j =', rtest, itest, jtest
            write(iulog,*) ' '
            do k = 1, upn
               write(iulog,*) k, model%velocity%uvel(k,i,j)*scyr, &
                             model%velocity%vvel(k,i,j)*scyr
            enddo
         endif

         !-------------------------------------------------------------------
         ! compute the basal velocity using the SSA solver
         !-------------------------------------------------------------------

         ! temporarily set the approximation to SSA
         model%options%which_ho_approx = HO_APPROX_SSA

         call t_startf('glissade_velo_higher_solver')
         call glissade_velo_higher_solve(model, ewn, nsn, upn)
         call t_stopf('glissade_velo_higher_solver')

         if (verbose_velo .and. this_rank == rtest) then
            i = itest
            j = jtest
            write(iulog,*) ' '
            write(iulog,*) 'SSA part of uvel, vvel (m/yr): r, i, j =', rtest, itest, jtest
            write(iulog,*) ' '
            do k = 1, upn
               write(iulog,*) k, model%velocity%uvel(k,i,j)*scyr, &
                             model%velocity%vvel(k,i,j)*scyr
            enddo
         endif

         !-------------------------------------------------------------------
         ! Add the result to the SIA velocity found above
         !-------------------------------------------------------------------

         model%velocity%uvel = model%velocity%uvel + uvel_sia
         model%velocity%vvel = model%velocity%vvel + vvel_sia

         ! restore the approximation option and clean up
         model%options%which_ho_approx = HO_APPROX_HYBRID

         deallocate(uvel_sia)
         deallocate(vvel_sia)

      else   ! standard higher-order solve
             ! can be BP, L1L2, SSA or SIA, depending on model%options%which_ho_approx

         ! Note: The geometry fields (thck, topg, and usrf) must be updated in halos
         !        before calling glissade_velo_higher_solve.
         !       These updates are done in subroutine glissade_diagnostic_variable_solve
         !        in module glissade.F90.

         call t_startf('glissade_velo_higher_solver')
         call glissade_velo_higher_solve(model, ewn, nsn, upn)
         call t_stopf('glissade_velo_higher_solver')

      endif   ! which_ho_approx

      ! optional diagnostics
      if (verbose_velo .and. this_rank == rtest) then
         i = itest
         j = jtest
         write(iulog,*) ' '
         write(iulog,*) 'uvel, vvel (m/yr): r, i, j =', rtest, itest, jtest
         write(iulog,*) ' '
         do k = 1, upn
            write(iulog,*) k, model%velocity%uvel(k,i,j)*scyr, &
                          model%velocity%vvel(k,i,j)*scyr
         enddo
      endif

    end subroutine glissade_velo_driver

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine glissade_wvelintg(uvel,  vvel,      &
                                 geomderv, numerics, &
                                 wgrd,  thck,  stagthck, &
                                 bmlt,  wvel)

      ! Calculate the vertical velocity field, wvel, by integrating
      !
      !    w(sigma) = -int_1^sigma [dU/dsigma (sigma) . (grad(s) - sigma*grad(H))
      !                             + H div(U(sigma))] dsigma + w(1)
      !
      !  (equation 13 in Payne and Dongelmans), where the thickness exceeds thklim.
      ! For Glissade, wvel is a diagnostic output only.
      ! Adapted from subroutine wvelintg in glide_velo.F90.

      use glimmer_global, only: dp
      use glissade_types, only: glissade_geomderv_type, glissade_numerics_type

      real(dp), dimension(:,:,:), intent(in) :: &
           uvel, vvel                ! horizontal velocity components on the staggered grid

      type(glissade_geomderv_type), intent(in) :: &
           geomderv                  ! horizontal derivatives of thickness and upper surface

      type(glissade_numerics_type), intent(in) :: &
           numerics                  ! numerical parameters, including sigma values

      real(dp), dimension(:,:), intent(in) :: &
           wgrd,                   & ! grid vertical velocity at the lowest model level
           thck,                   & ! ice thickness (m)
           stagthck,               & ! ice thickness on the staggered grid (m)
           bmlt                      ! basal melt rate; see Payne and Dongelmans, equation 14

      real(dp), dimension(:,:,:), intent(out) :: &
           wvel                      ! vertical velocity

      ! local variables

      real(dp) :: dew16, dns16       ! inverse of grid spacings multiplied by 16
      real(dp), dimension(6) :: cons ! temporary local values of derivatives
      integer :: ns, ew, up          ! loop indices
      integer :: nsn, ewn, upn       ! domain sizes

      real(dp), dimension(size(uvel,1)) :: &
           dupsw,                  & ! sigma(up+1) - sigma(up)
           depthw,                 & ! (sigma(up+1) + sigma(up)) / 2
           suvel, svvel              ! sums of uvel and vvel over the 4 vertices of a cell

      upn = size(uvel,1) ; ewn = size(uvel,2) ; nsn = size(uvel,3)

      dupsw  = (/ (numerics%sigma(up+1) - numerics%sigma(up), up=1,upn-1), 0.0d0 /)
      depthw = (/ ((numerics%sigma(up+1) + numerics%sigma(up)) / 2.0d0, up=1,upn-1), 0.0d0 /)

      dew16 = 1.d0/(16.0d0 * numerics%dew)
      dns16 = 1.d0/(16.0d0 * numerics%dns)

      do ns = 2, nsn
         do ew = 2, ewn
            if (thck(ew,ns) > numerics%thklim) then

               ! bottom boundary condition
               wvel(upn,ew,ns) = wgrd(ew,ns) - bmlt(ew,ns)

               ! local values of thickness and surface elevation derivatives
               cons(1) = sum(geomderv%dusrfdew(ew-1:ew,ns-1:ns)) / 16.0d0
               cons(2) = sum(geomderv%dthckdew(ew-1:ew,ns-1:ns)) / 16.0d0
               cons(3) = sum(geomderv%dusrfdns(ew-1:ew,ns-1:ns)) / 16.0d0
               cons(4) = sum(geomderv%dthckdns(ew-1:ew,ns-1:ns)) / 16.0d0
               cons(5) = sum(stagthck(ew-1:ew,ns-1:ns))
               cons(6) = cons(5)*dns16
               cons(5) = cons(5)*dew16

               suvel(:) = uvel(:,ew-1,ns-1) + uvel(:,ew-1,ns) + uvel(:,ew,ns-1) + uvel(:,ew,ns)
               svvel(:) = vvel(:,ew-1,ns-1) + vvel(:,ew-1,ns) + vvel(:,ew,ns-1) + vvel(:,ew,ns)

               ! integrate upward from the bed
               do up = upn-1, 1, -1
                  wvel(up,ew,ns) = wvel(up+1,ew,ns) &
                       - dupsw(up) * cons(5) * (sum(uvel(up:up+1,ew,ns-1:ns))  - sum(uvel(up:up+1,ew-1,ns-1:ns))) &
                       - dupsw(up) * cons(6) * (sum(vvel(up:up+1,ew-1:ew,ns))  - sum(vvel(up:up+1,ew-1:ew,ns-1))) &
                       - (suvel(up+1) - suvel(up)) * (cons(1) - depthw(up) * cons(2)) &
                       - (svvel(up+1) - svvel(up)) * (cons(3) - depthw(up) * cons(4))
               enddo

            else   ! too thin; set wvel = 0

               wvel(:,ew,ns) = 0.0d0

            endif
         enddo
      enddo

    end subroutine glissade_wvelintg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module glissade_velo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
