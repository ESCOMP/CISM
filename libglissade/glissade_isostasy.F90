!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_isostasy.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glissade_isostasy

  !-------------------------------------------------------------------------
  ! Some notes on the glissade_isostasy module (WHL, July 2026):
  !
  ! This module, glissade_isostasy, supersedes the isostasy module in libglide.
  ! Most of the code below is very similar to that in isostasy.F90, with subroutine
  !  names changed to avoid conflicts.
  ! I left isostasy.F90 in place for now (in libglide) but may remove it later.
  !
  ! Changes from isostasy.F90:
  ! * There are just two public subroutines: glissade_isostasy_init and glissade_isostasy_solve.
  !   The rest are private.
  ! * Calculations specific to the elastic lithosphere and relaxing asthenosphere
  !   have been move to glissade_isostasy_elra.

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
  !     * Relaxing asthenosphere with a laterally varying relaxation factor: asthenosphere = 2.
  !     The parameter tau_relax_const is the relaxation time scale for asthenosphere = 1; the default is 3000 yr.
  !     The 2D field tau_relax sets the relaxation time scale for asthenosphere = 2; it is read from an input file.
  ! (4) The which_relaxed parameter determines how the relaxed topography (relx) is computed.
  !     This is the topography we would have eventually (after the asthenosphere fully relaxes) with zero load.
  !     The asthenosphere calculation continually adjusts the topography toward topg = relx - load.
  !     There are three options:
  !     * which_relaxed = 0, the default. Both topg and relx, if present, are read from an input file.
  !       If relx is missing from the input file, the model sets relx = 0.
  !     * which_relaxed = 1. The model sets relx to the input topg. This is appropriate if the model
  !       is initializaed with no ice load, but for an existing ice sheet will be incorrect.
  !     * which_relaxed = 2. The input 'topg' field is interpreted as the equilibrium topography,
  !       given the input load. The relaxed topography is computed at initialization as relx = topg + load;
  !       it retains this value on restart. This setting is appropriate if the topography has had time
  !       to adjust fully since the last major change in load, or if ongoing isostatic adjustment is small
  !       compared to the adjustment to be simulated.
  !-------------------------------------------------------------------------

  ! Calculate isostatic adjustment due to changing surface loads
  ! Note: This module currently supports an ELRA scheme (elastic lithosphere, relaxing asthenosphere).
  !       At some point, it could wrap more complex schemes such as those in the FastIsostasy model.

  use glimmer_global, only : dp
  use glimmer_paramets, only: iulog
  use glimmer_physcon, only: scyr
  use glimmer_utils, only: point_diag
  use glimmer_log
  use cism_parallel, only: this_rank, main_task

  implicit none

  private
  public :: glissade_isostasy_init, glissade_isostasy_solve, verbose_isostasy

  logical, parameter :: verbose_isostasy = .true.

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine glissade_isostasy_init(model)

    !> initialise isostasy calculations
    use glide_types
    use glissade_isostasy_elastic, only: glissade_init_elastic
    use cism_parallel, only: parallel_is_zero

    implicit none

    type(glide_global_type) :: model

    if (model%options%isostasy == ISOSTASY_COMPUTE) then

       if (model%isostasy%lithosphere == LITHOSPHERE_ELASTIC) then
          ! initialize the elastic lithosphere
          call glissade_init_elastic(model%isostasy%rbel, model%numerics%dew)
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

       ! convert asthenosphere relaxation timescale from yr to s
       if (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_CONST) then
          if (model%isostasy%tau_relax_const == 0.0d0) then
             call write_log('tau_relax_const must be nonzero with this asthenosphere option', GM_FATAL)
          endif
          model%isostasy%tau_relax_const = model%isostasy%tau_relax_const * scyr
       elseif (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_LATVAR) then
          if (parallel_is_zero(model%isostasy%tau_relax)) then
             call write_log('tau_relax must be nonzero with this asthenosphere option', GM_FATAL)
          endif
       endif

    endif   ! isostasy_compute

    ! Handle relaxed topography

    select case(model%isostasy%which_relaxed)

    case(RELAXED_TOPO_DEFAULT)
       ! relx, if present in the input file, is read in directly and is distinct from topg;
       ! nothing to do here

    case(RELAXED_TOPO_INPUT)
       ! supplied input topography is relaxed; set relx = topg
       model%isostasy%relx = model%geometry%topg

    case(RELAXED_TOPO_COMPUTE)

       ! supplied topography is in equilibrium with the load;
       ! compute relx based on topg and load

       if (model%options%is_restart == STANDARD_RESTART) then
           ! relx should have been read from the restart file
           if (parallel_is_zero(model%isostasy%relx)) then
              call write_log ('Failed to read relx on restart with which_relaxed = RELAXED_TOPO_COMPUTE', &
                   GM_FATAL)
           endif
        else
           ! Since relx will be computed as topg + load, it should not be present in the input file
           ! Note: For a hybrid restart with 'relx' present in the input restart file,
           !       the user should set which_relaxed = RELAXED_TOPO_STANDARD instead.
           if (.not.parallel_is_zero(model%isostasy%relx)) then
              call write_log ('Do not set which_relaxed = RELAXED_TOPO_COMPUTE if relx is in the input file')
              call write_log ('Either remove relx or set which_relaxed = RELAXED_TOPO_STANDARD', GM_FATAL)
           endif
           ! Compute the load, then compute relx = topg + load
           call isostasy_relaxed(model)
        endif

    end select

  end subroutine glissade_isostasy_init

!-------------------------------------------------------------------------
  
  subroutine glissade_isostasy_solve(model)

    ! ------------------------------------------------------------------------ 
    ! Calculate isostatic adjustment
    ! ------------------------------------------------------------------------ 

    ! ------------------------------------------------------------------------
    ! Note: glissade_isostasy_solve is called near the beginning of glissade_tstep,
    !       just after the previous velocity solve. Following are some old comments
    !       on the question of when to compute isostasy.
    !
    ! Matt Hoffman wrote:
    ! Consider for a forward Euler time step:
    ! With a relaxing mantle model, topg is a prognostic (time-evolving) variable:
    !      topg1 = f(topg0, thk0, ...)
    ! However, for a fluid mantle where the adjustment is instantaneous, topg is a diagnostic variable
    ! (comparable to calculating floatation height of ice in the ocean):
    !      topg1 = f(thk1)
    ! In either case, the topg update should be separate from the thickness evolution (because thk1 = f(thk0, vel0=g(topg0,...)).
    ! However, if the isostasy calculation needs topg0, the icewaterload call should be made BEFORE thck is updated.
    ! If the isostasy calculation needs topg1, the icewaterload call should be made AFTER thck is updated.
    ! Also, we should think about when marinlim, usrf, lsrf, derivatives should be calculated relative to the topg update via isostasy.
    !
    ! WHL writes (May 2017):
    ! When isostasy is turned on, it is usually run with a relaxing mantle.
    ! With the call moved to the start of the time step, both the icewaterload call (if needed) and
    !  the relaxation are done before the ice thickness update. So we have
    !       topg1 = f(topg0, thk0, ...)
    !  followed by
    !       thk1  = f(thk0, vel0=g(topg0,...)
    ! In other words, thk1 is independent of topg1. I think this is what is desired.
    ! 
    ! ------------------------------------------------------------------------

    use glide_types
    use glissade_masks, only: glissade_marine_connection_mask
    use cism_parallel, only: parallel_type, parallel_halo, parallel_halo_extrapolate

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! --- Local variables ---

    type(parallel_type) :: parallel   ! info for parallel communication

    !WHL - debug - for isostasy hack
!    integer :: itest, jtest, rtest
!    itest = model%numerics%idiag_local
!    jtest = model%numerics%jdiag_local
!    rtest = model%numerics%rdiag_local
    
    parallel = model%parallel

    ! ------------------------------------------------------------------------
    ! update the ice/water load at the prescribed interval
    ! ------------------------------------------------------------------------

    if (model%options%isostasy == ISOSTASY_COMPUTE) then

       if (model%isostasy%nlith > 0) then
          if (mod(model%numerics%tstep_count, model%isostasy%nlith) == 0) then

             ! isostasy hack:
!             if (this_rank == rtest) write(iulog,*) 'Isostasy hack: Reduce thck by 20 m'
!             model%geometry%thck = model%geometry%thck - 20.0d0
!             model%geometry%thck = max(model%geometry%thck, 0.0d0)
!             call point_diag(model%geometry%thck, 'adjusted thck', itest, jtest, rtest, 7, 7)

             call isostasy_icewaterload(model)
             model%isostasy%new_load = .true.
          end if
       endif  ! nlith > 0

    end if
   
    ! ------------------------------------------------------------------------ 
    ! Calculate isostatic adjustment
    ! ------------------------------------------------------------------------ 

    if (model%options%isostasy == ISOSTASY_COMPUTE) then

       call compute_isostasy(model)

       ! update topography in halo cells
       ! Note: For outflow BCs, most fields (thck, usrf, temp, etc.) are set to zero in the global halo,
       !        to create ice-free conditions. However, we might not want to set topg = 0 in the global halo,
       !        because then the global halo will be interpreted as ice-free land, whereas we may prefer to
       !        treat it as ice-free ocean. For this reason, topg is extrapolated from adjacent cells.
       !       Similarly, for no_ice BCs, we want to zero out ice state variables adjacent to the global boundary,
       !        but we do not want to zero out the topography.
       ! Note: The topg halo update at initialization has an optional argument periodic_ew,
       !        which is needed for ismip-hom. I doubt ismip-hom will be run with active isostasy,
       !        but the argument is included to be on the safe side.
       !TODO: Do we need similar logic for halo updates of relx?

       if (model%general%global_bc == GLOBAL_BC_OUTFLOW) then
          call parallel_halo_extrapolate(model%geometry%topg, parallel)
       elseif (model%general%global_bc == GLOBAL_BC_NO_ICE) then
          call parallel_halo(model%geometry%topg, parallel, zero_global_boundary_no_ice_bc = .false.)
       else  ! other global BCs, including periodic
          call parallel_halo(model%geometry%topg, parallel, &
                          periodic_offset_ew = model%numerics%periodic_offset_ew, &
                          periodic_offset_ns = model%numerics%periodic_offset_ns)
       endif

       ! update the marine connection mask, which depends on topg

       call glissade_marine_connection_mask(&
            model%general%ewn,          model%general%nsn,          &
            parallel,                                               &
            model%numerics%idiag_local, model%numerics%jdiag_local, &
            model%numerics%rdiag_local,                             &
            model%geometry%thck,        model%geometry%topg,        &
            model%climate%eus,          0.0d0,                      &  ! thklim = 0
            model%geometry%marine_connection_mask)

    end if

  end subroutine glissade_isostasy_solve

!-------------------------------------------------------------------------

  subroutine isostasy_icewaterload(model)

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

  end subroutine isostasy_icewaterload

!-------------------------------------------------------------------------

  subroutine isostasy_relaxed(model)

    !  Calculate the relaxed topography, assuming the isostatic depression
    !   is the equilibrium state for the current topography.
    !  Note: This subroutine is called only at initialization; should not
    !        be called on restart..
    
    use glide_types
    implicit none

    type(glide_global_type) :: model

    ! Calculate the load
    call isostasy_icewaterload(model)

    ! Apply lithosphere model
    call isostasy_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)

    ! Add to present topography to get relaxed topography
    model%isostasy%relx = model%geometry%topg + model%isostasy%load

  end subroutine isostasy_relaxed

!-------------------------------------------------------------------------
  
  subroutine compute_isostasy(model)

    !> calculate isostatic adjustment due to changing surface loads

    use glide_types
    implicit none

    type(glide_global_type) :: model

    ! update load if it is time to do so
    if (model%isostasy%new_load) then

       call isostasy_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)

       ! update bedrock if the mantle is fluid (non-viscous)
       if (model%isostasy%asthenosphere == ASTHENOSPHERE_FLUID) then
          model%geometry%topg = model%isostasy%relx - model%isostasy%load
       end if

       model%isostasy%new_load = .false.

    end if

    ! update bedrock if the mantle is relaxing
    if (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_CONST .or. &
        model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_LATVAR) then
       call relaxing_mantle(model)
    end if

  end subroutine compute_isostasy

!-------------------------------------------------------------------------

  subroutine isostasy_lithosphere(model, load, load_factors)

    ! Update the lithosphere

    use glide_types
    use glissade_isostasy_elastic, only: glissade_calc_elastic
    implicit none

    ! input/output arguments
    !TODO - units of load and load_factors?

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(out) :: load         !> loading effect due to load_factors
    real(dp), dimension(:,:), intent(in)  :: load_factors !> load mass divided by mantle density

    ! local variables

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
       endif
       
       call glissade_calc_elastic(&
            model%isostasy%rbel,  &
            load_factors,         &
            load,                 &
            model%parallel)

       if (verbose_isostasy) then
          call point_diag(load_factors, 'input load_factors', itest, jtest, rtest, 7, 7)
          call point_diag(load, 'load after calc_elastic', itest, jtest, rtest, 7, 7)
       endif

    end if

  end subroutine isostasy_lithosphere

!-------------------------------------------------------------------------

  subroutine relaxing_mantle(model)

    ! Approximate the mantle with a relaxing half-space: dh/dt = -1/tau*(w-h)
    ! The relaxation timescale can be either a constant (tau_relax_const) or a
    !  laterally varying 2D field read in at initialization (tau_relax).
    ! Both dt and tau_relax have units of seconds.

    use glide_types
    implicit none

    type(glide_global_type) :: model
    
    integer :: i, j
    real(dp) :: ft1, ft2
    integer :: itest, jtest, rtest

    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local
    rtest = model%numerics%rdiag_local

    if (verbose_isostasy) then
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'relaxing_mantle, time (yr) =', model%numerics%time
       endif
       call point_diag(model%isostasy%relx, 'relx', itest, jtest, rtest, 7, 7)
       call point_diag(model%isostasy%relx - model%isostasy%load, 'relx - load', itest, jtest, rtest, 7, 7)
       call point_diag(model%geometry%topg, 'topg before relaxation', itest, jtest, rtest, 7, 7)
    endif

    if (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_CONST) then

       ft1 = exp(-model%numerics%dt/model%isostasy%tau_relax_const)
       ft2 = 1.d0 - ft1
       do j = 1, model%general%nsn
          do i = 1, model%general%ewn
             model%geometry%topg(i,j) = ft2 * (model%isostasy%relx(i,j) - model%isostasy%load(i,j)) &
                                      + ft1 *  model%geometry%topg(i,j)
          end do
       end do

       if (verbose_isostasy .and. this_rank == rtest) then
          write(iulog,*) 'tau_relax_const (yr)', model%isostasy%tau_relax_const/scyr
       endif

    elseif (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_LATVAR) then

       do j = 1, model%general%nsn
          do i = 1, model%general%ewn
             if (model%isostasy%tau_relax(i,j) > 0.0d0) then
                ft1 = exp(-model%numerics%dt/model%isostasy%tau_relax(i,j))
             else
                ft1 = exp(-model%numerics%dt/model%isostasy%tau_relax_const)
             endif
             ft2 = 1.d0 - ft1
             model%geometry%topg(i,j) = ft2 * (model%isostasy%relx(i,j) - model%isostasy%load(i,j)) &
                                      + ft1 *  model%geometry%topg(i,j)
          end do
       end do

       if (verbose_isostasy) then
          call point_diag(model%isostasy%tau_relax/scyr, 'tau_relax (yr)', itest, jtest, rtest, 7, 7)
       endif

    endif

    if (verbose_isostasy) then
       call point_diag(model%geometry%topg, 'topg after relaxation', itest, jtest, rtest, 7, 7)
    endif

  end subroutine relaxing_mantle

!-------------------------------------------------------------------------

end module glissade_isostasy

!-------------------------------------------------------------------------
