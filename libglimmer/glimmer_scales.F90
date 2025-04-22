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

module glimmer_scales

  use glimmer_global, only : dp

  implicit none

!!  real(dp) :: scale_uvel, scale_uflx, scale_diffu, scale_acab, scale_wvel, scale_btrc 
!!  real(dp) :: scale_beta, scale_flwa, scale_tau, scale_efvs, scale_resid

contains

  !TODO - Remove this subroutine and module; no longer used
  !       glide_vars.def now uses scyr as its main scaling parameter
  subroutine glimmer_init_scales

    ! set scale factors for I/O (can't have non-integer powers)

    use glimmer_physcon, only : scyr
!!    use glimmer_paramets, only : thk0, tim0, vel0, vis0, len0, acc0, tau0, evs0
    implicit none

!!  !WHL - Reset these scales such that 'scyr' is the only scaling factor; remove thk0, vel0, etc.
    !      With this scaling, velocity has units of m/yr in I/O, but m/s in the code,
    !      and similarly for other variables

!!!    scale_scyr = scyr
!!!    scale_inv_scyr = 1.0d0/scyr

!!    scale_uvel  = scyr * vel0                     ! uvel, vvel, ubas, vbas, etc.
!    scale_uvel  = scyr                               ! uvel, vvel, ubas, vbas, etc.
!!    scale_uflx  = scyr * vel0 * thk0              ! uflx, vflx
!    scale_uflx  = scyr                             ! uflx, vflx
!!    scale_diffu = scyr * vel0 * len0              ! diffu
!    scale_diffu = scyr                           ! diffu
!!    scale_acab = scyr * thk0 / tim0
!    scale_acab  = scyr                            ! acab, bmlt
!!    scale_wvel  = scyr * thk0 / tim0              ! wvel, wgrd
!    scale_wvel  = scyr                            ! wvel, wgrd
!!    scale_btrc  = scyr * vel0 * len0 / (thk0**2)  ! btrc, soft
!    scale_btrc  = scyr                           ! btrc, soft
    
!!    scale_beta  = tau0 / vel0 / scyr              ! units: Pa * sec/m * yr/sec = Pa * yr/m 
!    scale_beta  = 1.0d0 / scyr                     ! units: Pa * sec/m * yr/sec = Pa * yr/m 
!!                                                ! NOTE: on i/o, beta has units of Pa yr/m. Since vel0 has units of m/s, 
                                                  ! the first two terms on the RHS have units of Pa s/m. Thus, the final 
                                                  ! division by scyr here converts s/m to yr/m. All together, the 3 terms 
                                                  ! on the RHS scale on i/o by Pa yr/m (thus, making dimensionless on input, 
                                                  ! assuming the units on input are Pa yr/m, and also converting to Pa yr/m on output)

!!    scale_flwa  = scyr * vis0                     ! flwa
!    scale_flwa  = scyr                            ! flwa
!!    scale_tau   = tau0                            ! tauf, tauxz, btractx
!    scale_tau   = 1.0d0                            ! tauf, tauxz, btractx
!!    scale_efvs  = evs0 / scyr                     ! efvs
!    scale_efvs  = 1.0d0 / scyr                     ! efvs
!!    scale_resid=  tau0 / len0                     ! resid_u, resid_v
!    scale_resid=  1.0d0                            ! resid_u, resid_v
 
  end subroutine glimmer_init_scales

end module glimmer_scales
