
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_paramets.F90 - part of the Community Ice Sheet Model (CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> model scaling constants
module glimmer_paramets

  use glimmer_global, only : dp
  use glimmer_physcon, only : scyr

  implicit none
  save

!WHL - logical parameter for code testing
!      If oldglide = T, the glide dycore will reproduce
!      (within roundoff error) the results
!      of Glimmer 1.0.18 for the dome and EISMINT-2 test cases.

  !TODO - Remove oldglide parameter when comparisons to old Glide are no longer desired
  logical, parameter :: oldglide = .false.
!  logical, parameter :: oldglide = .true.

!TODO - redundant output units (stdout and glimmer_unit) 
!           It is redundant to define both stdout (which is public) and 
!            glimmer_unit (which is private to glimmer_log.F90).
!           However, it is sometimes convenient to write to stdout in Glimmer
!            without calling write_log.  
!           May want to delete this later (and declare stdout in glc_constants 
!            for CESM runs).

  integer :: stdout = 6

! logical flag to turn on special DEBUG output (related to test points), false by default

  logical :: GLC_DEBUG = .false.

!TODO: Redefine scaling parameters to have SI or similar units?
!      Considered removing these parameters from the code, but may be too much work.
!
!          Note: If tau0 is redefined in the code in terms of rhoi and grav,
!          then all the scaling parameters could be written in terms of scyr.
!
!          See comments below for details.


  ! Unphysical value used for initializing certain variables (e.g., temperature) so we can tell
  !  later if they were read from an input file or otherwise computed correctly
  real(dp), parameter :: unphys_val = -99999.d0

  real(dp), parameter :: netcdf_fill_value = 9.96921d+36

  ! Other numerical constants
  real(dp), parameter :: eps08 = 1.0d-08  ! small number, useful for some thresholds
  real(dp), parameter :: eps10 = 1.0d-10  ! another small number
  real(dp), parameter :: eps11 = 1.0d-11  ! another small number

  ! scaling parameters

  ! The fundamental scaling parameters are thk0, len0, and vel0. The others are derived from these.

!SCALING - DFM, 2, Oct 2012 - made scaled vs. unscaled values for thk0, len0, 
! and vel0 switchable by the reconstituted NO_RESCALE compilation flag. 
! (necessary to be compatible with alternate dycores) 

#ifndef NO_RESCALE
! The following are the old Glimmer scaling parameters.
  real(dp), parameter :: thk0 = 2000.0d0        ! m 
  real(dp), parameter :: len0 = 200.0d3         ! m 
  real(dp), parameter :: vel0 = 500.d0 / scyr   ! m yr^{-1} converted to S.I. units
!!  real(dp), parameter :: vis0 = 5.70d-18 / scyr  ! yr^{-1} Pa^{-3} converted to S.I. units
#else
! (no rescaling)
  real(dp), parameter :: thk0 = 1.d0        ! no scaling of thickness
  real(dp), parameter :: len0 = 1.d0        ! no scaling of length
  real(dp), parameter :: vel0 = 1.d0 / scyr ! yr * s^{-1}  
! end (no rescaling)
#endif

  !Note: Both the SIA and HO solvers fail unless tim0 = len0/vel0. Not sure if this can be changed.
  !      With the revised scaling, tim0 = scyr.
  real(dp), parameter :: tim0 = len0 / vel0          ! s
  real(dp), parameter :: acc0 = thk0 * vel0 / len0   ! m s^{-1}

!Note - With thk0 = 1, can replace tau0 by rhoi*grav in code and remove stress scaling.
!       Similarly can redefine vis0 and evs0

!Note - The constants rhoi_glam and grav_glam are declared here as parameters because the parameters
!        tau0, evs0 and vis0 depend on them.
!       The values of rhoi and grav used elsewhere in the code are declared in glimmer_physcon but are
!        not parameters, because they can be overridden by the user in the config file.
  real(dp), parameter :: rhoi_glam = 910.0d0          ! kg m^{-3}
  real(dp), parameter :: grav_glam = 9.81d0           ! m s^{-2}

  ! GLAM scaling parameters; units are correct if thk0 has units of meters
  integer, parameter :: gn = 3                              ! Glen flow exponent; fixed at 3 for purposes of setting vis0
  real(dp), parameter :: tau0 = rhoi_glam*grav_glam*thk0    ! stress scale in GLAM ( Pa )  
  real(dp), parameter :: evs0 = tau0 / (vel0/len0)          ! eff. visc. scale in GLAM ( Pa s )
  real(dp), parameter :: vis0 = tau0**(-gn) * (vel0/len0)   ! rate factor scale in GLAM ( Pa^-3 s^-1 )

!SCALING - This is the scaling we would use if we had velocity in m/yr and thk0 = len0 = 1.
!  real(dp), parameter :: thk0 = 1.d0
!  real(dp), parameter :: len0 = 1.d0
!  real(dp), parameter :: vel0 = 1.d0 / scyr
!  real(dp), parameter :: tim0 = scyr
!  real(dp), parameter :: acc0 = 1.d0 / scyr
!  real(dp), parameter :: tau0 = rhoi_glam*grav_glam
!  real(dp), parameter :: evs0 = tau0*scyr
!  real(dp), parameter :: vis0 = tau0**(-gn) / scyr

!WHL - Here I am defining some new constants that have the same values as thk0, len0, etc. in old Glimmer.
!      I am giving the new constants new names to minimize confusion.
!      These are used in only a few places.  For instance, we have this in glide_thck:
!
!          residual = maxval(abs(model%geometry%thck-model%thckwk%oldthck2))
!
!      In old Glimmer, thk0 = 2000 m and thck = O(1)
!      In new CISM, thk0 = 1 and thck = true thickness in meters
!      With thk0 = 1, we need to divide the rhs by 2000 m to reproduce the results of old Glimmer.
!      The following code satisfies either of the two conventions:
!
!          residual = maxval( abs(model%geometry%thck-model%thckwk%oldthck2) * (thk0/thk_scale) )

  real(dp), parameter :: thk_scale = 2000.0d0        ! m
  real(dp), parameter :: len_scale = 200.0d3         ! m
  real(dp), parameter :: vel_scale = 500.0 / scyr    ! m yr^{-1} converted to S.I. units
  real(dp), parameter :: tau_scale = rhoi_glam*grav_glam*thk_scale       ! stress scale in GLAM ( Pa )  
  real(dp), parameter :: vis_scale = tau_scale**(-gn) * (vel_scale/len_scale)  ! rate factor scale in GLAM ( Pa^-3 s^-1 )
  real(dp), parameter :: evs_scale = tau_scale / (vel_scale/len_scale)   ! eff. visc. scale in GLAM ( Pa s )

end module glimmer_paramets
