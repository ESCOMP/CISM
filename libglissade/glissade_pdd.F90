!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_pdd.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains drivers for handling the positive degree day (PDD) calculations.
!
! Author: Fairuz Ishraque
!         <ifairuz@ucar.edu>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module glissade_pdd
    use glimmer_global, only: dp
    use glimmer_paramets, only: iulog
    use glimmer_log
    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo, lhalo, uhalo, &
         parallel_halo, parallel_reduce_max, parallel_reduce_sum, parallel_globalindex

    implicit none
    save
    private

    use glimmer_physcon, only: pi, scyr, rhow, rhoi
    ! Then derive what you need:
	real(dp), parameter :: seconds_per_day = 86400.0d0
	real(dp)            :: days_per_year   ! not a parameter -- derived at runtime
	real(dp)            :: dt_series, h_days

    days_per_year = scyr / seconds_per_day   ! = 365.0 always in CISM
	dt_series     = scyr / real(model%pdd%n_series, dp)
	h_days        = dt_series / seconds_per_day

!=======================================================================

contains

!=======================================================================
	subroutine glissade_pdd_init(model)
    !-------------------------------------------------------------------
	! Initialize the PDD scheme
	! Set parameters, allocate arrays, precomputes sinusoid amplitude/phase.
	! Call once during model init, after geometry halo updates.
	!--------------------------------------------------------------------
		
	end subroutine glissade_pdd_init

    !-----------------------------------------------------------------------
	pure function calov_greve_integrand(sigma, TacC) result(val)
	!-----------------------------------------------------------------------
	! Compute the integrand of the Calov-Greve (2005) PDD expectation integral.
	!
	! Equation (6) in Calov & Greve (2005): for a temperature drawn from a
	! Gaussian distribution N(TacC, sigma^2), this gives the expected number
	! of positive degree days per unit time.
	!
	! If sigma = 0, reduces to the deterministic limit: max(TacC, 0).
	!
	! Mirrors CalovGreveIntegrand() in PISM localMassBalance.cc.
	!
	! Arguments:
	!   sigma - std dev of temperature variability (deg C)
	!   TacC  - temperature above PDD threshold (deg C)
	! Returns:
	!   val   - integrand (deg C); multiply by h_days (day) to get [K day]
	!-----------------------------------------------------------------------
		
		real(dp), intent(in) :: sigma, TacC
		real(dp)             :: val
		real(dp)             :: Z
		
		if (sigma == 0.0d0) then
		   val = max(TacC, 0.0d0)
		   return
		end if
		
		Z   = TacC / (sqrt(2.0d0) * sigma)
		val = (sigma / sqrt(2.0d0 * pi)) * exp(-Z * Z) &
			+ (TacC / 2.0d0) * erfc(-Z)
		
	end function calov_greve_integrand

end module glissade_pdd