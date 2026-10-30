
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

  implicit none
  save

  !TODO - Clean up and remove most of the parameters.
  !       Maybe put the rest in a single module, cism_constants.

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

  !WHL - Changing to iulog, following CTSM
!!  integer :: stdout = 6
  integer :: iulog = 6

! logical flag to turn on special DEBUG output (related to test points), false by default

  logical :: GLC_DEBUG = .false.

  ! Unphysical value used for initializing certain variables (e.g., temperature) so we can tell
  !  later if they were read from an input file or otherwise computed correctly
  real(dp), parameter :: unphys_val = -99999.d0
  real(dp), parameter :: netcdf_fill_value = 9.96921d+36

  ! Other numerical constants
  real(dp), parameter :: eps08 = 1.0d-08  ! small number, useful for some thresholds
  real(dp), parameter :: eps10 = 1.0d-10  ! another small number
  real(dp), parameter :: eps11 = 1.0d-11  ! another small number

  ! scaling parameters
  !WHL - Removed the old Glimmer scaling parameters, April 2025
  !      These included thk0, len0, vel0, tim0, acc0, vis0, tau0, and evs0.
  !      A few modules keep local parameters with the old values (e.g., thk_scale = 2000 m
  !      in glide_thck, velo_scale = 500 m/yr in glide_temp and glissade_velo_higher).

end module glimmer_paramets
