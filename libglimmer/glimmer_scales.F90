!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_scales.F90 - part of the Community Ice Sheet Model (CISM)  
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

! This module holds scales for various fields
!TODO - Remove this subroutine and module; no longer used.
!       glide_vars.def now uses scyr as its main scaling parameter

module glimmer_scales

  use glimmer_global, only : dp

  implicit none

contains

  subroutine glimmer_init_scales

    ! set scale factors for I/O (can't have non-integer powers)

    use glimmer_physcon, only : scyr
    implicit none

    !WHL - Removed the old Glimmer scale factors, April 2025
    !      These included scale_uvel, scale_uflx, scale_diff, scale_acab, scale_wvel,
    !       scale_btrc, scale_beta, scale_flwa, scale_tau, scale_efvs, and scale_resid.
    !      Once thk0, vel0, tim0, etc. were removed the code (effectively setting them to 1.0),
    !       the only scale factor remaining is scyr, which is used, for example, to convert
    !       velocity from m/yr in I/O to m/s in the code, and similarly for other quantities
    !       expressed in terms of yr or yr^(-1) in I/O and s or s^(-1) in the code.

 
  end subroutine glimmer_init_scales

end module glimmer_scales
