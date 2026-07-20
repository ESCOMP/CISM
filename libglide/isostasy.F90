!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   isostasy.F90 - part of the Community Ice Sheet Model (CISM)  
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

module isostasy

  !-------------------------------------------------------------------------
  ! Some notes on the isostasy calculation (WHL, May 2017; updated July 2026):
  !
  ! The isostasy calculation has been parallelized since the original Glimmer release,
  !  but otherwise the physical is similar. The most common configuration is ELRA =
  !  elastic lithosphere, relaxing asthenosphere.
  !
  ! The following config settings are relevant to the isostasy.
  ! All of these are set in the [isostasy] section unless otherwise specified.
  ! (1) To run with isostasy, set isostasy = 1 in the [options] section.
  !     The default is 0 (no isostasy).
  ! (2) There are two lithosphere options:
  !     * Local lithosphere: lithosphere = 0
  !     * Elastic lithosphere: lithosphere = 1; this is the default
  !     The parameter load_update_interval determines how often the elastic load is updated.
  !     The default is 10 yr.  As long as the load is not recomputed too often, the cost of isostasy
  !     should be minimal compared to the whole simulation.
  !     The parameter flexural_rigidity controls the elastic rigidity; the default is 0.24e25 N m.
  ! (3) There are three asthenosphere options:
  !     * Fluid asthenosphere: asthenosphere = 0
  !     * Relaxing asthenosphere with a constant relaxation factor: asthenosphere = 1; this is the default.
  !     The parameter tau_relax_const is the relaxation time scale for asthenosphere = 1; the default is 3000 yr.
  !     Note: A third option (asthenosphere = 2, with a laterally varying time scale)
  !           is supported for the Glissade dycore.
  ! (4) The which_relaxed parameter determines how the relaxed topography (relx) is computed.
  !     This is the topography we would have eventually (after the asthenosphere fully relaxes) with zero load.
  !     The asthenosphere calculation continually adjusts the topography toward topg = relx - load.
  !     There are three options:
  !     * which_relaxed = 0, the default. Both topg and relx, if present, are read from an input file.
  !       If relx is missing from the input file, the model sets relx = 0.
  !     * which_relaxed = 1. The model sets relx to the input topg. This is appropriate if the model
  !       is initializaed with no ice load, but for an existing ice sheet will be incorrect.
  !     * which_relaxed = 2. The input 'topg' field is interpreted as the equilibrium topography,
  !       given the input load. The relaxed topography is computed at initialization as relx = topg + load.
  !       This setting is appropriate if the topography has had time to adjust fully since the last major change
  !       in load, or if ongoing isostatic adjustment is small compared to the adjustment to be simulated.
  ! Note: Some of these options may not work correctly with the older Glide dycore.
  !-------------------------------------------------------------------------

  !> calculate isostatic adjustment due to changing surface loads

  use glimmer_global, only : dp
  use glimmer_paramets, only: iulog
  use glimmer_physcon, only: scyr
  use glimmer_utils, only: point_diag
  use cism_parallel, only: main_task, this_rank

  implicit none

  private :: relaxing_mantle

  logical, parameter :: verbose_isostasy = .true.

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine init_isostasy(model)

    !> initialise isostasy calculations
    use glide_types
    use glimmer_physcon,  only: scyr
    use isostasy_elastic, only: init_elastic

    implicit none

    type(glide_global_type) :: model

    if (model%isostasy%lithosphere == LITHOSPHERE_ELASTIC) then

       call init_elastic(model%isostasy%rbel, model%numerics%dew)

    end if

    !-----------------------------------------------------------------
    ! Based on the update period, determine how frequently the lithosphere load should be updated.
    ! The load is updated every nlith timesteps.
    ! An integer is used instead of a real number to decide when to update, in order to avoid roundoff issues.
    ! NOTE: The ratio isostasy%period/tinc is rounded to the nearest integer.
    !       Use numerics%tinc because it has units of years (like isostasy%period), whereas numerics%dt has model timeunits.
    !-----------------------------------------------------------------

    if (model%isostasy%load_update_interval > 0.0d0) then
       model%isostasy%nlith = nint(model%isostasy%load_update_interval / model%numerics%tinc)
    else
       model%isostasy%nlith = 0  ! never update
    endif

    model%isostasy%tau_relax_const = model%isostasy%tau_relax_const * scyr

   end subroutine init_isostasy

!-------------------------------------------------------------------------
  
  subroutine isos_icewaterload(model)

    !> calculate surface load factors due to water and ice distribution

    use glimmer_physcon
    use glide_types
    implicit none

    type(glide_global_type) :: model

    real(dp) :: ice_mass, water_depth, water_mass
    integer :: ew,ns
  
     do ns=1,model%general%nsn
       do ew=1,model%general%ewn
          ice_mass = rhoi * model%geometry%thck(ew,ns)

          if (model%geometry%topg(ew,ns) - model%climate%eus < 0.d0) then   ! check if we are below sea level

             water_depth = model%climate%eus - model%geometry%topg(ew,ns)
             water_mass = rhoo * water_depth

             ! Just the water load due to changes in sea-level
             model%isostasy%load_factors(ew,ns) = rhoo* model%climate%eus/rhom

             ! Check if ice is not floating
             if ( ice_mass > water_mass ) then
                model%isostasy%load_factors(ew,ns) = model%isostasy%load_factors(ew,ns) + (ice_mass - water_mass)/rhom
             end if

          else                                       ! bedrock is above sea level

             model%isostasy%load_factors(ew,ns) = ice_mass/rhom

          end if

       end do
    end do

  end subroutine isos_icewaterload

!-------------------------------------------------------------------------

  subroutine isos_compute(model)

    !> calculate isostatic adjustment due to changing surface loads

    use glide_types
    implicit none

    type(glide_global_type) :: model

    ! update load if it is time to do so
    if (model%isostasy%new_load) then

       call isos_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)

       ! update bedrock if the mantle is fluid (non-viscous)
       if (model%isostasy%asthenosphere == ASTHENOSPHERE_FLUID) then
          model%geometry%topg = model%isostasy%relx - model%isostasy%load
       end if

       model%isostasy%new_load = .false.

    end if

    ! update bedrock if the mantle is relaxing
    if (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_CONST) then
       call relaxing_mantle(model)
    end if

  end subroutine isos_compute

!-------------------------------------------------------------------------

  subroutine isos_lithosphere(model,load,load_factors)

    use glide_types
    use isostasy_elastic, only: calc_elastic
    implicit none

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(out) :: load !> loading effect due to load_factors
    real(dp), dimension(:,:), intent(in)  :: load_factors !> load mass divided by mantle density

    integer :: itest, jtest, rtest
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local
    rtest = model%numerics%rdiag_local

    if (model%isostasy%lithosphere == LITHOSPHERE_LOCAL) then

       load = load_factors

    else if (model%isostasy%lithosphere == LITHOSPHERE_ELASTIC) then

       if (verbose_isostasy) then
             if (main_task) then
                write(iulog,*) 'Update lithospheric load: time, tstep_count, nlith =', &
                     model%numerics%time, model%numerics%tstep_count, model%isostasy%nlith
             endif
          call point_diag(load_factors, 'input load_factors', itest, jtest, rtest, 7, 7)
       endif

       call calc_elastic(&
            model%isostasy%rbel,  &
            load_factors,         &
            load,                 &
            model%parallel,       &
            model%numerics%idiag, &
            model%numerics%jdiag, &
            itest, jtest, rtest)

       if (verbose_isostasy) then
          call point_diag(load, 'load after calc_elastic', itest, jtest, rtest, 7, 7)
       endif

    end if

  end subroutine isos_lithosphere

!-------------------------------------------------------------------------

  subroutine isos_relaxed(model)

    !> Calculate the relaxed topography, assuming the isostatic depression
    !> is the equilibrium state for the current topography.

    use glide_types
    implicit none
    type(glide_global_type) :: model

    ! Calculate the load
    call isos_icewaterload(model)

    ! Apply lithosphere model
    call isos_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)

    ! Add to present topography to get relaxed topography
    model%isostasy%relx = model%geometry%topg + model%isostasy%load

  end subroutine isos_relaxed

!-------------------------------------------------------------------------
! private subroutines
!-------------------------------------------------------------------------

  subroutine relaxing_mantle(model)

    !> approximate mantle with a relaxing half-space: dh/dt=-1/tau*(w-h)
    use glide_types
    implicit none
    type(glide_global_type) :: model
    
    integer :: ew,ns
    real(dp) :: ft1, ft2

    integer :: itest, jtest, rtest
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local
    rtest = model%numerics%rdiag_local

    ft1 = exp(-model%numerics%dt/model%isostasy%tau_relax_const)
    ft2 = 1.d0 - ft1

    if (verbose_isostasy) then
       if (this_rank == rtest) then
          write(iulog,*) 'relaxing_mantle, time (yr) =', model%numerics%time
          write(iulog,*) 'tau, dt/tau, relative change =', &
               model%isostasy%tau_relax_const, model%numerics%dt/model%isostasy%tau_relax_const, ft2
       endif
       call point_diag(model%isostasy%relx, 'relx', itest, jtest, rtest, 7, 7)
       call point_diag(model%isostasy%relx - model%isostasy%load, 'relx - load', itest, jtest, rtest, 7, 7)
       call point_diag(model%geometry%topg, 'topg before relaxation', itest, jtest, rtest, 7, 7)
    endif


    do ns=1,model%general%nsn
       do ew=1,model%general%ewn
          model%geometry%topg(ew,ns) = ft2 * (model%isostasy%relx(ew,ns) - model%isostasy%load(ew,ns)) &
                                     + ft1 *  model%geometry%topg(ew,ns)
       end do
    end do

    if (verbose_isostasy) then
       call point_diag(model%geometry%topg, 'topg after relaxation', itest, jtest, rtest, 7, 7)
    endif

  end subroutine relaxing_mantle

!-------------------------------------------------------------------------

end module isostasy

!-------------------------------------------------------------------------
