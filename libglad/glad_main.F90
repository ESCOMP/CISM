!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_main.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!
!   Copyright (C) 2005-2014
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

module glad_main

  ! This module provides an interface to GCMs in the case where fields have already been
  ! downscaled to the ice sheet grid (and the GCM does its own upscaling from the ice
  ! sheet grid to the land grid).
  !
  ! This only provides code for the SMB case, not for the PDD case.

  use glimmer_global, only: dp, fname_length
  use glad_type
  use glad_constants

  use glimmer_paramets, only: stdout, GLC_DEBUG

  implicit none
  private
  
  ! ------------------------------------------------------------
  ! glad_params derived type definition
  ! This is where default values are set.
  ! ------------------------------------------------------------

  type, public :: glad_params 

     !> Derived type containing parameters relevant to all instances of 
     !> the model - i.e. those parameters which pertain to the global model. 

     ! Ice model instances --------------------------------------

     integer                                   :: ninstances = 1       !> Number of ice model instances
     character(fname_length),pointer,dimension(:) :: config_fnames => null()    ! array of config filenames
     type(glad_instance),pointer,dimension(:) :: instances  => null() !> Array of glimmer\_instances

     ! Global model parameters ----------------------------------

     integer  :: tstep_mbal = 1        !> Mass-balance timestep (hours)
     integer  :: start_time            !> Time of first call to glad (hours)
     integer  :: time_step             !> Calling timestep of global model (hours)

     ! Parameters that can be set by the GCM calling Glad

     logical  :: gcm_restart = .false. !> If true, restart the model from a GCM restart file
     character(fname_length) :: gcm_restart_file   !> Name of restart file
     integer  :: gcm_fileunit = 99     !> Fileunit specified by GCM for reading config files

  end type glad_params

  !---------------------------------------------------------------------------------------
  ! Use of the routines here:
  !
  ! In model initialization:
  ! - Call glad_initialize once
  ! - Call glad_initialize_instance once per instance
  ! - Call glad_initialization_wrapup once
  !
  ! In the model run loop:
  ! - Call glad_gcm once per instance
  !---------------------------------------------------------------------------------------
  
  public :: glad_initialize
  public :: glad_initialize_instance
  public :: glad_initialization_wrapup
  
  public :: glad_gcm

  
  !---------------------------------------------------------------------------------------
  ! Some notes on coupling to the Community Earth System Model (CESM).  These may be applicable
  ! for coupling to other GCMs:
  !
  ! When coupled to CESM, Glad receives two fields from the coupler on the ice sheet grid:
  !   qsmb = surface mass balance (kg/m^2/s)
  !   tsfc = surface ground temperature (deg C)
  ! Both qsmb and tsfc are computed in the CESM land model.
  ! Seven fields are returned to CESM on the ice sheet grid:
  !   ice_covered = whether a grid cell is ice-covered [0,1]
  !   topo = surface elevation (m)
  !   hflx = heat flux from the ice interior to the surface (W/m^2)
  !   rofi = ice runoff (i.e., calving) (kg/m^2/s)
  !   rofl = liquid runoff (i.e., basal melting; the land model handles sfc runoff) (kg/m^2/s)
  !   ice_sheet_grid_mask = mask of ice sheet grid coverage
  !   icemask_coupled_fluxes = mask of ice sheet grid coverage where we are potentially
  !     sending non-zero fluxes
  !
  ! Note about ice_sheet_grid_mask and icemask_coupled_fluxes: ice_sheet_grid_mask is
  ! non-zero wherever CISM is operating - i.e., grid cells with icesheet or bare land (but
  ! not ocean). icemask_coupled_fluxes is similar, but is 0 for icesheet instances that
  ! have zero_gcm_fluxes = .true. Thus, icemask_coupled_fluxes can be used to determine
  ! the regions of the world in which CISM is operating and potentially sending non-zero
  ! fluxes to the climate model.
  !
  ! The land model has the option to update its ice coverage and surface elevation, given
  ! the fields returned from Glad.
  !
  ! There are two driver subroutines in this module for CESM coupling: 
  !  initialise_glad_gcm (for initialization) and glad_gcm (for timestepping).
  !
  !---------------------------------------------------------------------------------------
  
contains

  subroutine glad_initialize(params, time_step, paramfile, daysinyear, start_time, &
                             gcm_restart, gcm_restart_file, gcm_debug, gcm_fileunit)

    ! Initialize the model for runs coupled to a GCM. This routine initializes variables
    ! shared between instances. See above for documentation of the full initialization
    ! sequence.

    ! Subroutine argument declarations --------------------------------------------------------
    
    type(glad_params),              intent(inout) :: params      !> parameters to be set
    integer,                           intent(in)  :: time_step   !> Timestep of calling model (hours)
    character(*),dimension(:),         intent(in)  :: paramfile   !> array of configuration filenames.
    integer,                  optional,intent(in)  :: daysinyear  !> Number of days in the year
    integer,                  optional,intent(in)  :: start_time  !> Time of first call to glad (hours)
    logical,                  optional,intent(in)  :: gcm_restart ! logical flag to restart from a GCM restart file
    character(*),             optional,intent(in)  :: gcm_restart_file ! restart filename for a GCM restart
                                                                  ! (currently assumed to be CESM)
    logical,                  optional,intent(in)  :: gcm_debug   ! logical flag from GCM to output debug information
    integer,                  optional,intent(in)  :: gcm_fileunit! fileunit for reading config files
    
    ! Internal variables -----------------------------------------------------------------------

    type(ConfigSection), pointer :: global_config
    
    ! Begin subroutine code --------------------------------------------------------------------


    if (present(gcm_debug)) then
       GLC_DEBUG = gcm_debug
    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Initializing glad'
    end if

    ! Initialise start time and calling model time-step (time_step = integer number of hours)
    ! We ignore t=0 by default 

    params%time_step = time_step

    ! Note: start_time = nhour_glad = 0 for an initial run.
    !       Does this create problems given that Glad convention is to ignore t = 0?

    if (present(start_time)) then
       params%start_time = start_time
    else
       params%start_time = time_step
    end if

    params%gcm_restart = .false.
    if (present(gcm_restart)) then
       params%gcm_restart = gcm_restart
    endif

    params%gcm_restart_file = ''
    if (present(gcm_restart_file)) then
       params%gcm_restart_file = gcm_restart_file
    endif

    params%gcm_fileunit = 99
    if (present(gcm_fileunit)) then
       params%gcm_fileunit = gcm_fileunit
    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'time_step     =', params%time_step
       write(stdout,*) 'start_time    =', params%start_time
    end if

    ! Initialise year-length -------------------------------------------------------------------

    if (present(daysinyear)) then
       call glad_set_year_length(daysinyear)
    end if

    ! ---------------------------------------------------------------
    ! Determine how many instances there are, according to what
    ! configuration files we've been provided with
    ! ---------------------------------------------------------------

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Read paramfile'
       write(stdout,*) 'paramfile =', paramfile
    end if

    if (size(paramfile) == 1) then
       ! Load the configuration file into the linked list
       call ConfigRead(process_path(paramfile(1)), global_config, params%gcm_fileunit)    
       ! Parse the list
       call glad_readconfig(global_config, params%ninstances, params%config_fnames, paramfile)
    else
       params%ninstances = size(paramfile)
       allocate(params%config_fnames(params%ninstances))
       params%config_fnames(:) = paramfile(:)
    end if

    allocate(params%instances(params%ninstances))

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Number of instances =', params%ninstances
    end if

  end subroutine glad_initialize

  subroutine glad_initialize_instance(params,         instance_index,        &
                                      ice_covered,    topo,                  &
                                      rofi,           rofl,           hflx,  &
                                      ice_sheet_grid_mask,                   &
                                      icemask_coupled_fluxes,                &
                                      output_flag)

    ! Initialize one instance in the params structure. See above for documentation of
    ! the full initialization sequence.
    
    ! Subroutine argument declarations --------------------------------------------------------

    type(glad_params),              intent(inout) :: params          !> parameters to be set
    integer,                         intent(in)    :: instance_index  !> index of current ice sheet instance

    real(dp),dimension(:,:),intent(out) :: ice_covered  ! whether each grid cell is ice-covered [0,1]
    real(dp),dimension(:,:),intent(out) :: topo         ! output surface elevation (m)
    real(dp),dimension(:,:),intent(out) :: hflx         ! output heat flux (W/m^2, positive down)
    real(dp),dimension(:,:),intent(out) :: rofi         ! output ice runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: rofl         ! output liquid runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: ice_sheet_grid_mask !mask of ice sheet grid coverage
    real(dp),dimension(:,:),intent(out) :: icemask_coupled_fluxes !mask of ice sheet grid coverage where we are potentially sending non-zero fluxes
    
    logical,                  optional,intent(out) :: output_flag !> Flag to show output set (provided for consistency)
    
    ! Internal variables -----------------------------------------------------------------------

    type(ConfigSection), pointer :: instance_config
    
    ! Begin subroutine code --------------------------------------------------------------------

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Read config file and initialize instance #', instance_index
    end if

    call ConfigRead(process_path(params%config_fnames(instance_index)),&
         instance_config, params%gcm_fileunit)

    call glad_i_initialise_gcm(instance_config,     params%instances(i),     &
                                params%g_grid,                                &
                                params%start_time,   params%time_step,        &
                                params%gcm_restart,  params%gcm_restart_file, &
                                params%gcm_fileunit )

    call set_output_fields(params%instances(instance_index), &
         ice_covered, topo, rofi, rofl, hflx, &
         ice_sheet_grid_mask, icemask_coupled_fluxes)

    if (present(output_flag)) output_flag = .true.
    
  end subroutine glad_initialize_instance

  subroutine glad_initialization_wrapup(params, ice_dt)

    type(glad_params),              intent(inout) :: params      !> parameters to be set
    integer,                  optional,intent(out) :: ice_dt      !> Ice dynamics time-step in hours

    ! Wrapup glad initialization - perform error checks, etc. See above for documentation
    ! of the full initialization sequence

    ! Check that all mass-balance time-steps are the same length and
    ! assign that value to the top-level variable

    params%tstep_mbal = check_mbts(params%instances(:)%mbal_tstep)

    if (present(ice_dt)) then
       ice_dt = check_mbts(params%instances(:)%ice_tstep)
    end if

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'tstep_mbal =', params%tstep_mbal
       write(stdout,*) 'start_time =', params%start_time
       write(stdout,*) 'time_step =',  params%time_step
       if (present(ice_dt)) write(stdout,*) 'ice_dt =', ice_dt
    end if

    ! Check time-steps divide into one another appropriately.

    if (.not. (mod (params%tstep_mbal, params%time_step) == 0)) then
       call write_log('The mass-balance timestep must be an integer multiple of the forcing time-step', &
                       GM_FATAL,__FILE__,__LINE__)
    end if
    

  end subroutine glad_initialization_wrapup
  
  !===================================================================

  subroutine glad_gcm(params,         instance_index, time,  &
                      qsmb,           tsfc,                  &
                      ice_covered,    topo,                  &
                      rofi,           rofl,           hflx,  &
                      ice_sheet_grid_mask,                   &
                      icemask_coupled_fluxes,                &
                      output_flag,    ice_tstep)

    ! Main Glad subroutine for GCM coupling.
    !
    ! It does all necessary temporal averaging, 
    ! and calls the dynamic ice sheet model when required. 
    !
    ! Input fields should be taken as means over the period since the last call.
    ! See the user documentation for more information.

    use glimmer_utils
    use glad_timestep, only: glad_i_tstep_gcm
    use glimmer_log
    use glimmer_paramets, only: scyr
    use parallel, only: main_task, tasks

    implicit none

    ! Subroutine argument declarations -------------------------------------------------------------

    type(glad_params),              intent(inout) :: params          !> parameters for this run
    integer,                         intent(in)    :: instance_index  !> index of current ice sheet instance
    integer,                         intent(in)    :: time            !> Current model time        (hours)

    real(dp),dimension(:,:),intent(in)    :: qsmb          ! input surface mass balance of glacier ice (kg/m^2/s)
    real(dp),dimension(:,:),intent(in)    :: tsfc          ! input surface ground temperature (deg C)

    real(dp),dimension(:,:),intent(out) :: ice_covered  ! whether each grid cell is ice-covered [0,1]
    real(dp),dimension(:,:),intent(out) :: topo         ! output surface elevation (m)
    real(dp),dimension(:,:),intent(out) :: hflx         ! output heat flux (W/m^2, positive down)
    real(dp),dimension(:,:),intent(out) :: rofi         ! output ice runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: rofl         ! output liquid runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: ice_sheet_grid_mask !mask of ice sheet grid coverage
    real(dp),dimension(:,:),intent(out) :: icemask_coupled_fluxes !mask of ice sheet grid coverage where we are potentially sending non-zero fluxes

    logical,optional,intent(out)   :: output_flag     ! Set true if outputs are set
    logical,optional,intent(out)   :: ice_tstep       ! Set when an ice dynamic timestep has been done
                                                      !  and new output is available

    ! Internal variables ----------------------------------------------------------------------------

    logical :: icets
    character(250) :: message

    integer :: av_start_time  ! value of time from the last occasion averaging was restarted (hours)

    ! Check input fields are correct ----------------------------------------------------------------

    ! Reset output flag

    if (present(output_flag)) output_flag = .false.
    if (present(ice_tstep))   ice_tstep = .false.

    ! Accumulate input fields for later averaging

    call accumulate_averages(params%instances(instance_index)%glad_inputs, &
         qsmb = qsmb, tsfc = tsfc, time = time)

    ! ---------------------------------------------------------
    ! If this is a mass balance timestep, prepare global fields, and do a timestep
    ! for each model instance
    ! ---------------------------------------------------------

    av_start_time = get_av_start_time(params%instances(instance_index)%glad_inputs)
    
    if (mod (time - av_start_time, params%time_step) /= 0) then
       
       write(message,*) 'Unexpected calling of GLAD at time ', time
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    
    else if (time - av_start_time + params%time_step > params%tstep_mbal) then

       write(message,*) &
            'Incomplete forcing of GLAD mass-balance time-step detected at time ', time
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
       
    else if (time - av_start_time + params%time_step == params%tstep_mbal) then

       ! Set output_flag

       ! At present, outputs are done for each mass-balance timestep, since
       ! that involved least change to the code. However, it might be good
       ! to change the output to occur with user-specified frequency.

       if (present(output_flag)) output_flag = .true.

       ! Do a timestep for this instance

       if (time == params%instances(instance_index)%next_time) then

          params%instances(instance_index)%next_time = &
               params%instances(instance_index)%next_time + &
               params%instances(instance_index)%mbal_tstep

          ! Calculate averages by dividing by number of steps elapsed
          ! since last model timestep.

          call calculate_averages(params%instances(instance_index)%glad_inputs, &
               qsmb = params%instances(instance_index)%acab, &
               tsfc = params%instances(instance_index)%artm)

          ! Calculate total surface mass balance - multiply by time since last model timestep
          ! Note on units: We want acab to have units of meters w.e. (accumulated over mass balance time step)
          ! Initial units are kg m-2 s-1 = mm s-1
          ! Divide by 1000 to convert from mm to m
          ! Multiply by hours2seconds = 3600 to convert from 1/s to 1/hr.  (tstep_mbal has units of hours)

          !TODO - Modify code so that qsmb and acab are always in kg m-2 s-1 water equivalent?
          params%instances(instance_index)%acab(:,:) = &
               params%instances(instance_index)%acab(:,:) * &
               params%tstep_mbal * hours2seconds / 1000.d0

          if (GLC_DEBUG .and. main_task) write(stdout,*) 'Take a glad time step, instance', instance_index
          call glad_i_tstep_gcm(time,                  &
               params%instances(instance_index),   &
               icets)

          call set_output_fields(params%instances(instance_index), &
               ice_covered, topo, rofi, rofl, hflx, &
               ice_sheet_grid_mask, icemask_coupled_fluxes)
          

          ! Set flag
          if (present(ice_tstep)) then
             ice_tstep = (ice_tstep .or. icets)
          end if

       endif   ! time = next_time

       ! ---------------------------------------------------------
       ! Reset averaging fields, flags and counters
       ! ---------------------------------------------------------

       call reset_glad_input_averages(params%instances(instance_index)%glad_inputs, &
            next_av_start = time + params%time_step)
       
       if (GLC_DEBUG .and. main_task) then
          write(stdout,*) 'Done in glad_gcm'
       endif

   endif    ! time - av_start_time + params%time_step > params%tstep_mbal

  end subroutine glad_gcm

  !===================================================================

  subroutine end_glad(params,close_logfile)

    !> tidy-up operations for Glad
    use glad_initialise
    use glimmer_log
    implicit none

    type(glad_params),intent(inout) :: params          ! parameters for this run
    logical, intent(in), optional    :: close_logfile   ! if true, then close the log file
                                                        ! (GCM may do this elsewhere)                                  
    integer :: i

    ! end individual instances

    do i = 1, params%ninstances
       call glad_i_end(params%instances(i))
    enddo

    if (present(close_logfile)) then
       if (close_logfile) call close_log
    else
       call close_log
    endif

    deallocate(params%config_fnames)
    deallocate(params%instances)

  end subroutine end_glad

  !----------------------------------------------------------------------
  ! PRIVATE INTERNAL GLIMMER SUBROUTINES FOLLOW.............
  !----------------------------------------------------------------------

  !TODO - Move subroutine glad_readconfig to a glad_setup module, in analogy to glide_setup?

  subroutine glad_readconfig(config, ninstances, fnames, infnames)

    !> Determine whether a given config file is a
    !> top-level glad config file, and return parameters
    !> accordingly.

    use glimmer_config
    use glimmer_log
    implicit none

    ! Arguments -------------------------------------------

    type(ConfigSection),      pointer :: config !> structure holding sections of configuration file
    integer,              intent(out) :: ninstances !> Number of instances to create
    character(fname_length),dimension(:),pointer :: fnames !> list of filenames (output)
    character(fname_length),dimension(:) :: infnames !> list of filenames (input)

    ! Internal variables ----------------------------------

    type(ConfigSection), pointer :: section
    character(len=100) :: message
    integer :: i

    if (associated(fnames)) nullify(fnames)

    call GetSection(config,section,'GLAD')
    if (associated(section)) then
       call GetValue(section,'n_instance',ninstances)
       allocate(fnames(ninstances))
       do i=1,ninstances
          call GetSection(section%next,section,'GLAD instance')
          if (.not.associated(section)) then
             write(message,*) 'Must specify ',ninstances,' instance config files'
             call write_log(message,GM_FATAL,__FILE__,__LINE__)
          end if
          call GetValue(section,'name',fnames(i))
       end do
    else
       ninstances=1
       allocate(fnames(1))
       fnames=infnames
    end if

    ! Print some configuration information

!!$    call write_log('GLAD global')
!!$    call write_log('------------')
!!$    write(message,*) 'number of instances :',params%ninstances
!!$    call write_log(message)
!!$    call write_log('')

  end subroutine glad_readconfig


  !========================================================

  integer function check_mbts(timesteps)

    !> Checks to see that all mass-balance time-steps are
    !> the same. Flags a fatal error if not, else assigns that
    !> value to the output

    use glimmer_log

    implicit none

    integer,dimension(:) :: timesteps !> Array of mass-balance timsteps

    integer :: n,i

    n = size(timesteps)
    if (n==0) then
       check_mbts = 0
       return
    endif

    check_mbts = timesteps(1)

    do i = 2,n
       if (timesteps(i) /= check_mbts) then
          call write_log('All instances must have the same mass-balance and ice timesteps', &
               GM_FATAL,__FILE__,__LINE__)
       endif
    enddo

  end function check_mbts

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glad_main

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
