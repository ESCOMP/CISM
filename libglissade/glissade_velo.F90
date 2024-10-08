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
#include "glide_nan.inc"
#include "glide_mask.inc"

module glissade_velo

    ! Driver for Glissade velocity solvers

    implicit none
    
contains
        
    subroutine glissade_velo_driver(model)

      ! Glissade higher-order velocity driver

      use glimmer_log
      use glimmer_paramets, only: vel0
      use glimmer_physcon, only: scyr
      use glide_types
      use glissade_velo_higher, only: glissade_velo_higher_solve
      use glissade_velo_sia, only: glissade_velo_sia_solve
      use glissade_velo_noism, only: glissade_velo_noism_solve
      use glide_mask
      use profile, only: t_startf, t_stopf

      type(glide_global_type),intent(inout) :: model

      real(dp), dimension(:,:,:), allocatable :: &
           uvel_sia, vvel_sia        ! temporary SIA velocity

      integer :: i, j, k
      integer :: ewn, nsn, upn
      integer :: itest, jtest, rtest
      integer :: whichbtrc_sav

      logical, parameter :: verbose_velo = .false.

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
      
      if (model%options%which_ho_approx == HO_APPROX_NOISM) then
 
         call glissade_velo_noism_solve(model, ewn, nsn, upn)

      elseif (model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
 
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

         if (verbose_velo .and. this_rank == rtest) then
            i = itest
            j = jtest
            print*, ' '
            print*, 'SIA part of uvel, vvel (m/yr): r, i, j =', rtest, itest, jtest
            print*, ' '
            do k = 1, upn
               print*, k, model%velocity%uvel(k,i,j)*(vel0*scyr), &
                          model%velocity%vvel(k,i,j)*(vel0*scyr)
            enddo
         endif

         !-------------------------------------------------------------------
         ! compute the basal velocity using the SSA solver
         !-------------------------------------------------------------------

         ! temporarily set the approximation to SSA
         model%options%which_ho_approx = HO_APPROX_SSA

         ! Compute mask for staggered grid. This is needed as an input to calcbeta
         ! (which used to be called here but now is called from glissade_velo_higher_solve).
         ! TODO - Remove the use of stagmask in the Glissade solver?

         call glide_set_mask(model%numerics,                                     &
                             model%geomderv%stagthck, model%geomderv%stagtopg,   &
                             model%general%ewn-1,     model%general%nsn-1,       &
                             model%climate%eus,       model%geometry%stagmask)

         call t_startf('glissade_velo_higher_solver')
         call glissade_velo_higher_solve(model, ewn, nsn, upn)
         call t_stopf('glissade_velo_higher_solver')

         if (verbose_velo .and. this_rank == rtest) then
            i = itest
            j = jtest
            print*, ' '
            print*, 'SSA part of uvel, vvel (m/yr): r, i, j =', rtest, itest, jtest
            print*, ' '
            do k = 1, upn
               print*, k, model%velocity%uvel(k,i,j)*(vel0*scyr), &
                          model%velocity%vvel(k,i,j)*(vel0*scyr)
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

         !-------------------------------------------------------------------
         ! Compute mask for staggered grid. This is needed as an input to calcbeta
         ! (which used to be called here but now is called from glissade_velo_higher_solve).
         ! TODO - Remove the use of stagmask in the Glissade solver?
         !-------------------------------------------------------------------

         call glide_set_mask(model%numerics,                                     &
                             model%geomderv%stagthck, model%geomderv%stagtopg,   &
                             model%general%ewn-1,     model%general%nsn-1,       &
                             model%climate%eus,       model%geometry%stagmask)

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
         print*, ' '
         print*, 'uvel, vvel (m/yr): r, i, j =', rtest, itest, jtest
         print*, ' '
         do k = 1, upn
            print*, k, model%velocity%uvel(k,i,j)*(vel0*scyr), &
                       model%velocity%vvel(k,i,j)*(vel0*scyr)
         enddo
      endif

    end subroutine glissade_velo_driver

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module glissade_velo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
