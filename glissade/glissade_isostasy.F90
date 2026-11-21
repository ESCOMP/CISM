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
  ! This module, glissade_isostasy, supersedes the isostasy module in glide (glide_isostasy).
  ! Most of the code below is very similar to that in glide_isostasy.F90, with subroutine
  !  names changed to avoid conflicts.
  !
  ! Changes from glide_isostasy.F90:
  ! * There are just two public subroutines: glissade_isostasy_init and glissade_isostasy_solve.
  !   The rest are private.
  ! * The elastic lithosphere calculation and the Kelvin functions it uses are private
  !   procedures at the end of this module (formerly the modules glissade_isostasy_elastic
  !   and glissade_isostasy_kelvin).

  !-------------------------------------------------------------------------
  ! Some notes on the isostasy calculation (WHL, May 2017; updated July 2026):
  !
  ! The isostasy calculation has been parallelized since the original Glimmer release,
  !  but otherwise the physical is similar. The most common configuration is ELRA =
  !  elastic lithosphere, relaxing asthenosphere.
  !
  ! The following config settings are relevant to the isostasy.
  ! All of these are set in the [isostasy] section unless otherwise specified.
  ! (1) To run with isostasy, set compute_isostasy = .true. in the [isostasy] section.
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
  ! (4) The input_topg_method option determines how the relaxed topography (relx) is computed.
  !     This is the topography we would have eventually (after the asthenosphere fully relaxes) with zero load.
  !     The asthenosphere calculation continually adjusts the topography toward topg = relx - load.
  !     There are three options:
  !     * input_topg_method = separate_relx, the default. Both topg and relx, if present, are read from an input file.
  !       If relx is missing from the input file, the model sets relx = 0.
  !     * input_topg_method = relx_from_topg. The model sets relx to the input topg. This is appropriate if the model
  !       is initializaed with no ice load, but for an existing ice sheet will be incorrect.
  !     * input_topg_method = relx_from_load. The input 'topg' field is interpreted as the equilibrium topography,
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
  use glimmer_physcon, only: scyr, pi
  use glimmer_utils, only: point_diag
  use glimmer_log
  use cism_parallel, only: this_rank, main_task
  use glissade_types, only: glissade_isos_elastic_type

  implicit none

  private
  public :: glissade_isostasy_init, glissade_isostasy_solve, verbose_isostasy

  logical :: verbose_isostasy = .false.
  logical :: verbose_elastic = .false.  ! if true, print diagnostic messages

  ! Kelvin functions (used for the elastic lithosphere)
  real(kind=dp), parameter :: euler_gamma=0.577215664901532860606512d0 !< Euler's constant
  integer :: j_max = 40 !< maximum number of iterations
  real(kind=dp) :: tolerance = 1.d-10 !< the tolerance

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine glissade_isostasy_init(model)

    !> initialise isostasy calculations
    use glissade_types
    use cism_parallel, only: parallel_is_zero

    implicit none

    type(glissade_global_type) :: model

    if (model%isostasy%compute_isostasy) then

       if (model%isostasy%lithosphere_method == LITHOSPHERE_METHOD_ELASTIC) then
          ! initialize the elastic lithosphere
          call elastic_init(model%isostasy%rbel, model%grid%dew, model%isostasy%flexural_rigidity)
       end if

       !-----------------------------------------------------------------
       ! Based on the update period, determine how frequently the lithosphere load should be updated.
       ! The load is updated every nlith timesteps.
       ! An integer is used instead of a real number to decide when to update, in order to avoid roundoff issues.
       ! NOTE: The ratio isostasy%period/tinc is rounded to the nearest integer.
       !       Use time%tinc because it has units of years (like isostasy%period), whereas time%dt has model timeunits.
       !-----------------------------------------------------------------

       if (model%isostasy%load_update_interval > 0.0d0) then
          model%isostasy%nlith = nint(model%isostasy%load_update_interval / model%time%tinc)
       else
          model%isostasy%nlith = 0  ! never update
       endif

       ! convert asthenosphere relaxation timescale from yr to s
       if (model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_RELAXING) then
          if (model%isostasy%tau_relax_const == 0.0d0) then
             call write_log('tau_relax_const must be nonzero with this asthenosphere option', GM_FATAL)
          endif
          model%isostasy%tau_relax_const = model%isostasy%tau_relax_const * scyr
       elseif (model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_RELAXING_LATVAR) then
          if (parallel_is_zero(model%isostasy%tau_relax)) then
             call write_log('tau_relax must be nonzero with this asthenosphere option', GM_FATAL)
          endif
       endif

    endif   ! isostasy_compute

    ! Handle relaxed topography

    select case(model%isostasy%input_topg_method)

    case(INPUT_TOPG_METHOD_SEPARATE_RELX)
       ! relx, if present in the input file, is read in directly and is distinct from topg;
       ! nothing to do here

    case(INPUT_TOPG_METHOD_RELX_FROM_TOPG)
       ! supplied input topography is relaxed; set relx = topg
       model%isostasy%relx = model%geometry%topg

    case(INPUT_TOPG_METHOD_RELX_FROM_LOAD)

       ! supplied topography is in equilibrium with the load;
       ! compute relx based on topg and load

       if (model%options%restart == STANDARD_RESTART) then
           ! relx should have been read from the restart file
           if (parallel_is_zero(model%isostasy%relx)) then
              call write_log ('Failed to read relx on restart with input_topg_method = relx_from_load', &
                   GM_FATAL)
           endif
        else
           ! Since relx will be computed as topg + load, it should not be present in the input file
           ! Note: For a hybrid restart with 'relx' present in the input restart file,
           !       the user should set input_topg_method = separate_relx instead.
           if (.not.parallel_is_zero(model%isostasy%relx)) then
              call write_log ('Do not set input_topg_method = relx_from_load if relx is in the input file')
              call write_log ('Either remove relx or set input_topg_method = separate_relx', GM_FATAL)
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

    use glissade_types
    use glissade_masks, only: glissade_marine_connection_mask
    use cism_parallel, only: parallel_type, parallel_halo, parallel_halo_extrapolate

    implicit none

    type(glissade_global_type), intent(inout) :: model   ! model instance

    ! --- Local variables ---

    type(parallel_type) :: parallel   ! info for parallel communication

    !WHL - debug - for isostasy hack
!    integer :: itest, jtest, rtest
!    itest = model%diagnostics%idiag_local
!    jtest = model%diagnostics%jdiag_local
!    rtest = model%diagnostics%rdiag_local
    
    parallel = model%parallel

    ! ------------------------------------------------------------------------
    ! update the ice/water load at the prescribed interval
    ! ------------------------------------------------------------------------

    if (model%isostasy%compute_isostasy) then

       if (model%isostasy%nlith > 0) then
          if (mod(model%time%tstep_count, model%isostasy%nlith) == 0) then

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

    if (model%isostasy%compute_isostasy) then

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

       if (model%grid%global_bc == GLOBAL_BC_OUTFLOW) then
          call parallel_halo_extrapolate(model%geometry%topg, parallel)
       elseif (model%grid%global_bc == GLOBAL_BC_NO_ICE) then
          call parallel_halo(model%geometry%topg, parallel, zero_global_boundary_no_ice_bc = .false.)
       else  ! other global BCs, including periodic
          call parallel_halo(model%geometry%topg, parallel, &
                          periodic_offset_ew = model%numerics%periodic_offset_ew, &
                          periodic_offset_ns = model%numerics%periodic_offset_ns)
       endif

       ! update the marine connection mask, which depends on topg

       call glissade_marine_connection_mask(&
            model%grid%ewn,                model%grid%nsn,                &
            parallel,                                                     &
            model%diagnostics%idiag_local, model%diagnostics%jdiag_local, &
            model%diagnostics%rdiag_local,                                &
            model%geometry%thck,           model%geometry%topg,           &
            model%climate%eus,             0.0d0,                         &  ! thck_threshold = 0
            model%geometry%marine_connection_mask)

    end if

  end subroutine glissade_isostasy_solve

!-------------------------------------------------------------------------

  subroutine isostasy_icewaterload(model)

    !> calculate surface load factors due to water and ice distribution

    use glimmer_physcon
    use glissade_types
    implicit none

    type(glissade_global_type) :: model

    real(dp) :: ice_mass, water_depth, water_mass
    integer :: ew,ns
  
     do ns=1,model%grid%nsn
       do ew=1,model%grid%ewn
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
    
    use glissade_types
    implicit none

    type(glissade_global_type) :: model

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

    use glissade_types
    implicit none

    type(glissade_global_type) :: model

    ! update load if it is time to do so
    if (model%isostasy%new_load) then

       call isostasy_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)

       ! update bedrock if the mantle is fluid (non-viscous)
       if (model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_FLUID) then
          model%geometry%topg = model%isostasy%relx - model%isostasy%load
       end if

       model%isostasy%new_load = .false.

    end if

    ! update bedrock if the mantle is relaxing
    if (model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_RELAXING .or. &
        model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_RELAXING_LATVAR) then
       call relaxing_mantle(model)
    end if

  end subroutine compute_isostasy

!-------------------------------------------------------------------------

  subroutine isostasy_lithosphere(model, load, load_factors)

    ! Update the lithosphere

    use glissade_types
    implicit none

    ! input/output arguments
    !TODO - units of load and load_factors?

    type(glissade_global_type) :: model
    real(dp), dimension(:,:), intent(out) :: load         !> loading effect due to load_factors
    real(dp), dimension(:,:), intent(in)  :: load_factors !> load mass divided by mantle density

    ! local variables

    integer :: itest, jtest, rtest

    itest = model%diagnostics%idiag_local
    jtest = model%diagnostics%jdiag_local
    rtest = model%diagnostics%rdiag_local

    if (model%isostasy%lithosphere_method == LITHOSPHERE_METHOD_LOCAL) then

       load = load_factors

    else if (model%isostasy%lithosphere_method == LITHOSPHERE_METHOD_ELASTIC) then

       if (verbose_isostasy) then
          if (main_task) then
             write(iulog,*) 'Update lithospheric load: time, tstep_count, nlith =', &
                  model%time%time, model%time%tstep_count, model%isostasy%nlith
          endif
       endif
       
       call elastic_calc(&
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

    use glissade_types
    implicit none

    type(glissade_global_type) :: model
    
    integer :: i, j
    real(dp) :: ft1, ft2
    integer :: itest, jtest, rtest

    itest = model%diagnostics%idiag_local
    jtest = model%diagnostics%jdiag_local
    rtest = model%diagnostics%rdiag_local

    if (verbose_isostasy) then
       if (this_rank == rtest) then
          write(iulog,*) ' '
          write(iulog,*) 'relaxing_mantle, time (yr) =', model%time%time
       endif
       call point_diag(model%isostasy%relx, 'relx', itest, jtest, rtest, 7, 7)
       call point_diag(model%isostasy%relx - model%isostasy%load, 'relx - load', itest, jtest, rtest, 7, 7)
       call point_diag(model%geometry%topg, 'topg before relaxation', itest, jtest, rtest, 7, 7)
    endif

    if (model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_RELAXING) then

       ft1 = exp(-model%time%dt/model%isostasy%tau_relax_const)
       ft2 = 1.d0 - ft1
       do j = 1, model%grid%nsn
          do i = 1, model%grid%ewn
             model%geometry%topg(i,j) = ft2 * (model%isostasy%relx(i,j) - model%isostasy%load(i,j)) &
                                      + ft1 *  model%geometry%topg(i,j)
          end do
       end do

       if (verbose_isostasy .and. this_rank == rtest) then
          write(iulog,*) 'tau_relax_const (yr)', model%isostasy%tau_relax_const/scyr
       endif

    elseif (model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_RELAXING_LATVAR) then

       do j = 1, model%grid%nsn
          do i = 1, model%grid%ewn
             if (model%isostasy%tau_relax(i,j) > 0.0d0) then
                ft1 = exp(-model%time%dt/model%isostasy%tau_relax(i,j))
             else
                ft1 = exp(-model%time%dt/model%isostasy%tau_relax_const)
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
! Elastic lithosphere
!-------------------------------------------------------------------------

  ! Code for an elastic lithosphere; typically combined with a relaxing asthenosphere
  !  Notes:
  !  * This code is based on the isostasy_elastic module from the old Glide code.
  !    It was copied to glissade, renamed and restructured in July 2026.
  !  * The elastic lithosphere calculation is done on a single task.
  !    When running in parallel runs, data are gathered onto the main task for the computation,
  !    then scattered back to local processors.
  !    This procedure does not scale well, although is is manageable on a 4 km mesh as long as
  !    the update is done at a frequency of once every few decades or less.

  subroutine elastic_init(rbel, deltax, flexural_rigidity)

    !> initialize elastic lithosphere calculations

    use glimmer_physcon, only : pi
    implicit none

    type(glissade_isos_elastic_type) :: rbel          !> structure holding elastic lithosphere data
    real(dp), intent(in) :: deltax      !> grid spacing
    real(dp), intent(in) :: flexural_rigidity   !> flexural rigidity of the lithosphere (N m)

    ! local variables
    real(dp) :: a     ! radius of disk
    real(dp) :: r     ! distance from centre
    integer :: i,j

    real(dp), parameter :: r_lr = 6.0d0   ! influence of disk load at (0,0) is felt within a radius of r_lr*rbel_r

    ! calculate a so that a circle of radius a is equivalent to a square with size deltax
    a = deltax/sqrt(pi)

    rbel%d = flexural_rigidity

    ! initialise w
    call init_rbel(rbel, a)

    ! calculate size of operator
    rbel%wsize = int(r_lr*rbel%lr/deltax)

    ! allocate memory for operator
    allocate(rbel%w(0:rbel%wsize,0:rbel%wsize))

    ! calculating points within disk
    rbel%w(0,0) = rbel_iw(rbel,0.d0)
    r = deltax/rbel%lr
    rbel%w(0,1) = rbel_iw(rbel,r)
    rbel%w(1,0) = rbel%w(0,1)

    ! calculating points outside disk
    do j=0,rbel%wsize
       do i=2,rbel%wsize
          r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
          rbel%w(i,j) = rbel_ow(rbel,r)
       end do
    end do

    do j=2,rbel%wsize
       do i=0,1
          r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
          rbel%w(i,j) = rbel_ow(rbel,r)
       end do
    end do

    i=1
    j=1
    r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
    rbel%w(i,j) = rbel_ow(rbel,r)

#ifdef DEB_REBOUND
    open(1,file='w.dat',status='UNKNOWN')
    do j=0,rbel%wsize
       do i=0,rbel%wsize
          write(1,*) i,j,rbel%w(i,j)
       end do
    end do
    close(1)
#endif

  end subroutine elastic_init

!-------------------------------------------------------------------------

  subroutine elastic_calc(&
       rbel,                 &
       load_factors,         &
       load,                 &
       parallel)

    !> Calculate surface loading using an elastic lithosphere approximation.
    !> The results match those of subroutine calc_elastic in the older Glide isostasy model.
    !> The main difference is that this subroutine uses a global gather and broadcast
    !  to enable each task to compute the load for locally owned cells only.

    use cism_parallel, only: this_rank, main_task, nhalo, parallel_type, &
         gather_var, broadcast, parallel_halo, parallel_globalindex

    implicit none

    ! input-output arguments
    type(glissade_isos_elastic_type) :: rbel                             !> structure holding elastic litho data
    real(dp), dimension(:,:), intent(in)  :: load_factors  !> load mass divided by mantle density
    real(dp), dimension(:,:), intent(out) :: load          !> loading effect due to load_factors

    type(parallel_type), intent(in) :: parallel            !> info for parallel communication

    ! local variables

    integer :: ewn, nsn                   ! grid dimensions on the local task; includes halo cells
    integer :: global_ewn, global_nsn     ! global grid dimensions

    integer :: i, j, n, m
    integer :: ig, jg                     ! global indices

    real(dp), dimension(:,:), allocatable :: &
         load_global,             & !> global version of the output 'load' array
         load_factors_global        !> global version of the input 'load_factors' array

    character(len=100) :: message

    ! initialize

    ewn = size(load,1)
    nsn = size(load,2)
    global_ewn = parallel%global_ewn
    global_nsn = parallel%global_nsn

    load = 0.0d0

    if (verbose_elastic .and. main_task) then
       write(iulog,*) 'In elastic_calc'
       write(iulog,*) 'local ewn/nsn =', ewn, nsn
       write(iulog,*) 'global_ewn/nsn =', global_ewn, global_nsn
    endif

    ! Gather the local load_factors arrays onto the main task
    ! Note: global arrays are allocated in the subroutine
    call gather_var(load_factors, load_factors_global, parallel)

    if (verbose_elastic .and. main_task) then
       if (sum(load_factors_global) > 0.0d0) then
          write(iulog,*) 'my_task, sum(load_factors_global) =', this_rank, sum(load_factors_global)
       endif
       write(iulog,*) 'Allocate load_factors_global'
    endif

    ! allocate load_factors_global on tasks other than main
    if (.not.main_task) then
       if (allocated(load_factors_global)) deallocate(load_factors_global)
       allocate(load_factors_global(global_ewn,global_nsn))
    endif

    if (verbose_elastic .and. main_task) then
       write(iulog,*) 'Broadcast load_factors_global to each task'
    endif

    ! broadcast load_factors_global from main_task to all processors
    call broadcast(load_factors_global)

    if (sum(load_factors_global) == 0.0d0) then
       write(message,*) 'Error, calc_elastic, sum(load_factors_global) = 0, my_task =', this_rank
       call write_log(message, GM_FATAL)
    endif

    if (verbose_elastic .and. main_task) then
       write(iulog,*) 'Compute load locally on each task'
    endif

    ! loop over locally owned cells
    do j = nhalo+1, nsn-nhalo
       do i = nhalo+1, ewn-nhalo
          call parallel_globalindex(i, j, ig, jg, parallel)

          ! Compute load terms by summing over cells in the radius of influence
          do n = max(1,jg-rbel%wsize), min(global_nsn,jg+rbel%wsize)
             do m = max(1,ig-rbel%wsize), min(global_ewn,ig+rbel%wsize)
                load(i,j) = load(i,j) + load_factors_global(m,n) * rbel%w(abs(m-ig),abs(n-jg))
             end do   ! m
          end do   ! n

       enddo   ! i
    enddo   ! j

    ! update halo cells
    call parallel_halo(load, parallel)

    ! deallocate global arrays
    deallocate(load_factors_global)

  end subroutine elastic_calc

!-------------------------------------------------------------------------

  subroutine init_rbel(rbel, a)

    !> initialize elastic lithosphere calculations

    use glimmer_physcon, only: rhom, grav
    implicit none

    type(glissade_isos_elastic_type) :: rbel        !> structure holding elastic litho data
    real(dp), intent(in) :: a         !> radius of disk

    real(dp) :: dummy_a

    call set_kelvin(1.d-10,40)

    rbel%lr = (rbel%d/(rhom*grav))**0.25d0
    rbel%a  = a

    dummy_a = rbel%a/rbel%lr
    
    rbel%c1  =  dummy_a * dker(dummy_a)
    rbel%c2  = -dummy_a * dkei(dummy_a)
    rbel%cd3 =  dummy_a * dber(dummy_a)
    rbel%cd4 = -dummy_a * dbei(dummy_a)
    
  end subroutine init_rbel

!-------------------------------------------------------------------------

  function rbel_ow(rbel,r)

    !> calculate deflection outside disk

    implicit none
    real(dp) :: rbel_ow
    real(dp), intent(in) :: r          !> radius, r should be scaled with lr
    type(glissade_isos_elastic_type) :: rbel     !> structure holding elastic litho data
    
    rbel_ow = rbel%cd3*ker(r) + rbel%cd4*kei(r)

  end function rbel_ow

!-------------------------------------------------------------------------

  function rbel_iw(rbel,r)

    !> calculate deflection inside disk
    implicit none
    real(dp) :: rbel_iw
    real(dp), intent(in) :: r          !> radius, r should be scaled with lr
    type(glissade_isos_elastic_type) :: rbel         !> structure holding elastic litho data
    
    rbel_iw = 1.d0 + rbel%c1*ber(r) + rbel%c2*bei(r)

  end function rbel_iw

!-------------------------------------------------------------------------
! Kelvin functions
!-------------------------------------------------------------------------

  !> Code for calculating zeroth order Kelvin functions and their derivatives.
  !!
  !! \author Magnus Hagdorn
  !! \date June 2000
  !! Added to the Glissade dycore by William Lipscomb, July 2026

  !> set tolerance and maximum number of iterations
  subroutine set_kelvin(tol, jmax)
    implicit none
    real(kind=dp), intent(in) :: tol
    integer, intent(in) :: jmax
    j_max = jmax
    tolerance = tol
  end subroutine set_kelvin

  function ber(x)
    implicit none
    real(kind=dp) :: ber
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_ber
    real(kind=dp) :: factorial
    real(kind=dp) :: sign
    integer :: j

    p_ber = 0.d0
    factorial = 1.d0

    ber = 1.d0
    arg = (x/2.d0)**4
    arg_d = arg
    sign = -1.d0

    j=1
    do while (j < j_max)
       p_ber = ber
       factorial = factorial*2*j*(2*j-1.d0)
       ber = ber + sign*arg_d/(factorial*factorial)
       if (abs(ber-p_ber) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function ber

  function bei(x)
    implicit none
    real(kind=dp) :: bei
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_bei
    real(kind=dp) :: factorial
    real(kind=dp) :: sign
    integer :: j

    p_bei = 1.d12
    factorial = 1.d0

    arg = (x/2.d0)**2
    bei = arg
    arg_d = arg*arg*arg
    arg = arg*arg
    sign = -1.d0

    j=1
    do while (j < j_max)
       p_bei = bei
       factorial = factorial*2*j*(2*j+1.d0)
       bei = bei + sign*arg_d/(factorial*factorial)
       if (abs(bei-p_bei) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function bei

  function ker(x)
    implicit none
    real(kind=dp) :: ker
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_ker
    real(kind=dp) :: factorial
    real(kind=dp) :: phi
    real(kind=dp) :: sign
    integer :: j

    p_ker = 0.d0
    factorial = 1.d0

    arg = (x/2.d0)**4
    arg_d = arg
    sign = -1.d0
    phi = 0.d0
    ker = -(log(x/2.d0)+euler_gamma)*ber(x)+(pi/4.d0)*bei(x)

    j=1
    do while (j < j_max)
       p_ker = ker
       factorial = factorial*2*j*(2*j-1.d0)
       phi = phi + 1.d0/(2.d0*j-1.d0) + 1.d0/(2.d0*j)
       ker = ker + sign*phi*arg_d/(factorial*factorial)
       if (abs(ker-p_ker) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function ker

  function kei(x)
    implicit none
    real(kind=dp) :: kei
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_kei
    real(kind=dp) :: factorial
    real(kind=dp) :: phi
    real(kind=dp) :: sign
    integer :: j

    p_kei = 0.d0
    factorial = 1.d0

    arg = (x/2.d0)**2
    sign = -1.d0
    phi = 1.d0
    kei = -(log(x/2.d0)+euler_gamma)*bei(x)-(pi/4.d0)*ber(x)+arg
    arg_d = arg
    arg = arg*arg
    arg_d = arg_d*arg

    j=1
    do while (j < j_max)
       p_kei = kei
       factorial = factorial*2*j*(2*j+1.d0)
       phi = phi + 1.d0/(2.d0*j+1.d0) + 1.d0/(2.d0*j)
       kei = kei + sign*phi*arg_d/(factorial*factorial)
       if (abs(kei-p_kei) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function kei

  function dber(x)
    implicit none
    real(kind=dp) :: dber
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_dber
    real(kind=dp) :: factorial
    real(kind=dp) :: sign
    integer :: j

    p_dber = 0.d0
    factorial = 1.d0

    dber = 0.d0
    arg = (x/2.d0)**4
    arg_d = (x/2.d0)**3
    sign = -1.d0

    j=1
    do while (j < j_max)
       p_dber = dber
       factorial = factorial*2*j*(2*j-1.d0)
       dber = dber + sign*2.d0*j*arg_d/(factorial*factorial)
       if (abs(dber-p_dber) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function dber

  function dbei(x)
    implicit none
    real(kind=dp) :: dbei
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_dbei
    real(kind=dp) :: factorial
    real(kind=dp) :: sign
    integer :: j

    p_dbei = 1.d12
    factorial = 1.d0

    arg = (x/2.d0)**4
    arg_d = arg*(x/2.d0)
    dbei = (x/2.d0)
    sign = -1.d0

    j=1
    do while (j < j_max)
       p_dbei = dbei
       factorial = factorial*2*j*(2*j+1.d0)
       dbei = dbei + sign*(2.d0*j+1.d0)*arg_d/(factorial*factorial)
       if (abs(dbei-p_dbei) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function dbei

  function dker(x)
    implicit none
    real(kind=dp) :: dker
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_dker
    real(kind=dp) :: factorial
    real(kind=dp) :: phi
    real(kind=dp) :: sign
    integer :: j

    p_dker = 0.d0
    factorial = 1.d0

    arg = (x/2.d0)**4
    arg_d = (x/2.d0)**3
    sign = -1.d0
    phi = 0.d0
    dker = -(log(x/2.d0)+euler_gamma)*dber(x)-ber(x)/x+(pi/4.d0)*dbei(x)

    j=1
    do while (j < j_max)
       p_dker = dker
       factorial = factorial*2*j*(2*j-1.d0)
       phi = phi + 1.d0/(2.d0*j-1.d0) + 1.d0/(2.d0*j)
       dker = dker + sign*phi*2.d0*j*arg_d/(factorial*factorial)
       if (abs(dker-p_dker) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function dker

  function dkei(x)
    implicit none
    real(kind=dp) :: dkei
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_dkei
    real(kind=dp) :: factorial
    real(kind=dp) :: phi
    real(kind=dp) :: sign
    integer :: j

    p_dkei = 0.d0
    factorial = 1.d0

    arg = (x/2.d0)
    sign = -1.d0
    phi = 1.d0
    dkei = -(log(x/2.d0)+euler_gamma)*dbei(x)-bei(x)/x-(pi/4.d0)*dber(x)+arg
    arg_d = arg**5
    arg = arg**4

    j=1
    do while (j < j_max)
       p_dkei = dkei
       factorial = factorial*2*j*(2*j+1.d0)
       phi = phi + 1.d0/(2.d0*j+1.d0) + 1.d0/(2.d0*j)
       dkei = dkei + sign*phi*(2.d0*j+1.d0)*arg_d/(factorial*factorial)
       if (abs(dkei-p_dkei) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function dkei

!-------------------------------------------------------------------------

end module glissade_isostasy

!-------------------------------------------------------------------------
