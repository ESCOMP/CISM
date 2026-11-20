!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_setup.F90 - part of the Community Ice Sheet Model (CISM)
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

module glissade_setup

  ! general routines for initialisation, etc, called from top-level glimmer subroutines
  !
  ! Configuration for the Glissade dycore, used in standalone runs and by GLAD.

  use glimmer_global, only: dp
  use glimmer_paramets, only : iulog

  implicit none

  private
  public :: glissade_config, &
            glissade_readconfig, glissade_printconfig, glissade_scale_params, &
            glissade_load_sigma, glissade_read_sigma, glissade_calc_sigma, &
            glissade_get_zocn, glissade_get_zatm

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine glissade_config(model,config)

    ! Read Glissade configuration from file and print it to the log

    use glissade_types
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    use glimmer_filenames

    implicit none

    type(glissade_global_type), intent(inout) :: model  ! model instance
    type(ConfigSection), pointer  :: config          ! structure holding sections of configuration file

    ! read configuration file
    call glissade_readconfig (model,config)

    ! check settings that involve more than one key, and set defaults that depend on other options
    call glissade_check_config(model, config)

    ! read the ocean levels ([ocean_data]), if present, or compute them from dzocn
    call glissade_get_zocn(model,config)

    ! read the atmosphere levels ([climate]), if present
    call glissade_get_zatm(model,config)

    ! write the configuration to the log
    call glissade_printconfig(model,config)

    ! read sigma levels from config file, if present
    ! Note: The sigma levels are written to the log when they are computed (glissade_load_sigma).
    call glissade_read_sigma(model,config)

    !WHL - Moved isostasy configuration to glissade_setup
!    call isos_readconfig(model%isos,config)
!    call isos_printconfig(model%isos)

    ! read mapping from config file
    ! **** Use of dew and dns here is an ugly fudge that
    ! **** allows the use of old [GLINT projection] config section
    ! **** for backwards compatibility. It will be deleted soon.
    ! **** (You have been warned!)
    ! **** N.B. Here, dew and dns are unscaled - i.e. real distances in m

    call glimmap_readconfig(model%projection,   config,   &
                            model%grid%dew,     model%grid%dns)

    ! netCDF I/O ([CF ...] sections)
    call glimmer_nc_readparams(model%funits, model%options%restart, config)

  end subroutine glissade_config

!-------------------------------------------------------------------------

  subroutine glissade_check_config(model, config)

    ! Check settings that involve keys in more than one section, and set the defaults
    ! that depend on another option. Called after glissade_readconfig and before glissade_printconfig.
    ! A dependent default applies only when its key is not in the config file;
    ! an explicit setting that conflicts with the other option stops the run.

    use glissade_types
    use glimmer_config
    use glimmer_log

    implicit none

    type(glissade_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config     !> structure holding sections of configuration file

    ! The local SIA solver cannot use first-order dissipation.
    ! With stokes_approx = local_sia, the default dissipation_method is sia.
    if (model%velocity%stokes_approx == STOKES_APPROX_LOCAL_SIA .and. &
        model%temperature%dissipation_method == DISSIPATION_METHOD_FIRST_ORDER) then
       if (key_given(config, 'temperature', 'dissipation_method')) then
          call write_log('Error, dissipation_method = first_order cannot be used with stokes_approx = local_sia; ' // &
               'use sia or none', GM_FATAL)
       else
          model%temperature%dissipation_method = DISSIPATION_METHOD_SIA
       endif
    endif

    ! Computing on active blocks only requires no-ice boundary conditions.
    ! With compute_blocks = active_only, the default global_bc is no_ice.
    if (model%numerics%compute_blocks == COMPUTE_BLOCKS_ACTIVE_ONLY .and. &
        model%grid%global_bc /= GLOBAL_BC_NO_ICE) then
       if (key_given(config, 'grid', 'global_bc')) then
          call write_log('Error, compute_blocks = active_only requires global_bc = no_ice', GM_FATAL)
       else
          model%grid%global_bc = GLOBAL_BC_NO_ICE
       endif
    endif

  end subroutine glissade_check_config

!-------------------------------------------------------------------------

  logical function key_given(config, secname, key)

    ! True if the config file sets this key in this section

    use glimmer_config
    implicit none
    type(ConfigSection), pointer :: config
    character(len=*), intent(in) :: secname, key

    type(ConfigSection), pointer :: section
    character(len=512) :: val

    call GetSection(config, section, secname)
    key_given = .false.
    if (associated(section)) key_given = ConfigSectionHasValue(section, key, val)

  end function key_given

!-------------------------------------------------------------------------

  subroutine glissade_readconfig(model,config)

    ! read Glissade configuration file
    ! Note: sigma coordinates are handled by a subsequent call to glissade_read_sigma
 
    use glissade_types
    use glimmer_config
    use glimmer_log
    implicit none
    type(glissade_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section

    ! The sections are read in the order in which they are written to the log (glissade_printconfig).
    ! A section is read whenever it is present; if it is absent, its keys keep their default values.
    ! Note: The presence of a section does not turn on its feature; the feature's switch does
    !       (e.g., isostasy, lateral_melt_method, bwat_method, enable_glaciers).
    !       The handlers only read values, so reading a section for a feature that is off
    !       has no effect on the results.
    ! Note: The ocean and atmosphere levels (nzocn, dzocn, zocn, nbasin in [ocean_data]; nzatm, zatm in [climate])
    !       are read by glissade_get_zocn and glissade_get_zatm, and [sigma] by glissade_read_sigma.

    ! read grid size parameters
    call GetSection(config,section,'grid')
    if (associated(section)) then
       call handle_grid(section, model)
    end if

    ! read time parameters
    call GetSection(config,section,'time')
    if (associated(section)) then
       call handle_time(section, model)
    end if

    ! read diagnostic options
    call GetSection(config,section,'diagnostics')
    if (associated(section)) then
       call handle_diagnostics(section, model)
    end if

    ! read options
    call GetSection(config,section,'options')
    if (associated(section)) then
       call handle_options(section, model)
    end if

    ! External dycores (BISICLES, Albany-FELIX) are no longer supported
    call GetSection(config,section,'external_dycore_options')
    if (associated(section)) then
       call write_log('Error, external dycores are no longer supported', GM_FATAL)
    end if

    ! read numerical options and parameters
    call GetSection(config,section,'numerics')
    if (associated(section)) then
       call handle_numerics(section, model)
    end if

    ! read physical constants
    call GetSection(config,section,'constants')
    if (associated(section)) then
       call handle_constants(section, model)
    end if

    ! read options for adjusting the input fields at initialization
    call GetSection(config,section,'initialization')
    if (associated(section)) then
       call handle_initialization(section, model)
    end if

    ! read velocity options (Stokes approximation, viscosity, surface gradient)
    call GetSection(config,section,'velocity')
    if (associated(section)) then
       call handle_velocity(section, model)
    end if

    ! read velocity solver options (nonlinear and linear solvers, matrix assembly)
    call GetSection(config,section,'velocity_solver')
    if (associated(section)) then
       call handle_velocity_solver(section, model)
    end if

    ! read ice advection options (thickness and tracer transport)
    call GetSection(config,section,'advection')
    if (associated(section)) then
       call handle_advection(section, model)
    end if

    ! read temperature options and parameters
    call GetSection(config,section,'temperature')
    if (associated(section)) then
       call handle_temperature(section, model)
    end if

    ! read flow law options and parameters
    call GetSection(config,section,'flow_law')
    if (associated(section)) then
       call handle_flow_law(section, model)
    end if

    ! read basal physics (sliding laws and basal friction parameters)
    call GetSection(config,section,'basal_physics')
    if (associated(section)) then
       call handle_basal_physics(section, model)
    end if

    ! read basal hydrology info
    call GetSection(config,section,'basal_hydro')
    if (associated(section)) then
       call handle_basal_hydro(section, model)
    end if

    ! read grounding-line options
    call GetSection(config,section,'grounding_line')
    if (associated(section)) then
       call handle_grounding_line(section, model)
    end if

    ! read basal melting of floating ice
    call GetSection(config,section,'basal_melt')
    if (associated(section)) then
       call handle_basal_melt(section, model)
    end if

    ! read calving
    call GetSection(config,section,'calving')
    if (associated(section)) then
       call handle_calving(section, model)
    end if

    ! read lateral melt info
    call GetSection(config,section,'lateral_melt')
    if (associated(section)) then
       call handle_lateral_melt(section, model)
    end if

    ! read ice removal (other than calving)
    call GetSection(config,section,'removal')
    if (associated(section)) then
       call handle_removal(section, model)
    end if

    ! read surface climate forcing (SMB, surface temperature)
    call GetSection(config,section,'climate')
    if (associated(section)) then
       call handle_climate(section, model)
    end if

    ! read ocean data (thermal forcing) used by the basal melt parameterizations
    call GetSection(config,section,'ocean_data')
    if (associated(section)) then
       call handle_ocean_data(section, model)
    end if

    ! read geothermal heat flux options
    call GetSection(config,section,'geothermal')
    if (associated(section)) then
       call handle_geothermal(section, model)
    end if

    ! read isostasy options and parameters
    call GetSection(config,section,'isostasy')
    if (associated(section)) then
       call handle_isostasy(section, model)
    end if

    ! read inversion parameters (basal friction and ocean temperature inversions)
    call GetSection(config,section,'inversion')
    if (associated(section)) then
       call handle_inversion(section, model)
    end if

    ! read glacier info
    call GetSection(config,section,'glaciers')
    if (associated(section)) then
       call handle_glaciers(section, model)
    end if

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file.
    ! (Glint restart variables are handled separately by Glint setup routines.)
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13

    call define_glissade_restart_variables(model, model%model_id)

  end subroutine glissade_readconfig

!-------------------------------------------------------------------------

  subroutine glissade_printconfig(model, config)

    !> print model configuration to log, one config section at a time
    !> Order: grid, time and run-level choices; the dynamical core (velocity, advection, temperature);
    !> physics from the ice interior to the bed to floating ice and other ice removal; external forcing
    !> from the atmosphere, ocean and solid Earth; and other sections.
    !> Note: The sigma levels are written later, when they are computed (glissade_load_sigma).
    use glimmer_log
    use glimmer_config
    use glissade_types
    implicit none
    type(glissade_global_type)  :: model !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file

    call write_log_div
    call print_grid(model)
    call print_time(model)
    call print_diagnostics(model)
    call print_options(model)
    call print_numerics(model)
    call print_constants(model)
    call print_initialization(model)
    call print_velocity(model)
    call print_velocity_solver(model)
    call print_advection(model)
    call print_temperature(model)
    call print_flow_law(model)
    call print_basal_physics(model)
    call print_basal_hydro(model)
    call print_grounding_line(model)
    call print_basal_melt(model)
    call print_calving(model)
    call print_lateral_melt(model)
    call print_removal(model)
    call print_climate(model, config)
    call print_ocean_data(model, config)
    call print_geothermal(model)
    call print_isostasy(model)
    call print_inversion(model)
    call print_glaciers(model)

  end subroutine glissade_printconfig

!-------------------------------------------------------------------------
    
  subroutine glissade_scale_params(model)

    !TODO - Most of the scaling involves converting input parameters from yr or yr^{-1} to s or s^{-1}.
    !       Some input parameters are not rescaled here, but are effectively rescaled during calculations
    !        or subroutine calls, e.g. writing (basal_melt%bmlt_float_omega/scyr) when bmlt_float_omega
    !        is used in a calculation with MKS units.
    !       It would be better to have a consistent approach, possibly by removing most of the conversions below.

    !> scale parameters
    use glissade_types
    use glimmer_physcon,  only: scyr

    implicit none

    type(glissade_global_type)  :: model !> model instance

    ! convert dt from yr to s
    model%time%dt     = model%time%tinc * scyr

    ! scale calving parameters
    model%calving%calving_timescale = model%calving%calving_timescale * scyr                ! convert from yr to s
    model%calving%cliff_timescale = model%calving%cliff_timescale * scyr    ! convert from yr to s

    ! scale basal traction parameter
    model%basal_physics%btrc_const = model%basal_physics%btrc_const/scyr

    ! scale inversion parameters
    model%inversion%inversion_babc_timescale = model%inversion%inversion_babc_timescale * scyr    ! convert yr to s
    model%inversion%inversion_deltaT_ocn_timescale = model%inversion%inversion_deltaT_ocn_timescale * scyr   ! yr to s

    ! scale SMB/acab parameterss
    model%climate%overwrite_acab_value = model%climate%overwrite_acab_value/scyr

    ! scale basal_hydro parameters
    model%basal_hydro%sliding_speed_fixed = model%basal_hydro%sliding_speed_fixed / scyr  ! m/yr to m/s
    model%basal_hydro%flwa_basal = model%basal_hydro%flwa_basal / scyr  ! Pa^{-3} yr^{-1} to Pa^{-3} s^{-1}

  end subroutine glissade_scale_params

!-------------------------------------------------------------------------

  subroutine glissade_read_sigma(model,config)

    ! read sigma levels from configuration file, if present
    ! called immediately after glissade_readconfig

    use glissade_types
    use glimmer_config
    use glimmer_log
    implicit none

    type(glissade_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file
        
    ! local variables
    type(ConfigSection), pointer :: section

    ! read sigma levels
    ! NOTE: The [sigma] section is ignored unless model%grid%sigma_method = SIGMA_METHOD_CONFIG

    if (model%grid%sigma_method == SIGMA_METHOD_CONFIG) then
       call GetSection(config,section,'sigma')
       if (associated(section)) then
          call handle_sigma(section, model)
       else
          call write_log('Error, sigma_method = config requires a [sigma] section with sigma_levels', GM_FATAL)
       end if
    endif

  end subroutine glissade_read_sigma

!-------------------------------------------------------------------------

  subroutine glissade_load_sigma(model,unit)

    ! Compute sigma coordinates or read them from a file
    ! Note: This subroutine is called from glissade_initialise.
    !       If sigma levels are provided in the config file, then they are read
    !        in by glissade_read_sigma, and model%grid%sigma_method is set to
    !        SIGMA_METHOD_CONFIG, in which case this subroutine does nothing.

    use glissade_types
    use glimmer_log
    use glimmer_filenames
    use cism_parallel, only: main_task, broadcast

    implicit none

    ! Arguments
    type(glissade_global_type),intent(inout) :: model !> Ice model to use
    integer,               intent(in)    :: unit   !> Logical file unit to use. 
                                                   !> (Must not already be in use)

    ! Internal variables

    integer :: up,upn
    logical :: there
    real(dp) :: level

    ! Beginning of code

    upn=model%grid%upn

    select case(model%grid%sigma_method)

    case(SIGMA_METHOD_GLIDE)   !  compute standard Glide sigma levels

       do up = 1,upn
          level = real(up-1,kind=dp) / real(upn-1,kind=dp)
          model%grid%sigma(up) = glissade_calc_sigma(level, 2.d0)
       end do

       call write_log('Computing Glide sigma levels')

    case(SIGMA_METHOD_FILE)        ! read from external file

       if (main_task) inquire (exist=there, file=process_path(model%grid%sigma_file))
       call broadcast(there)
       if (.not.there) then
          call write_log('Sigma levels file: '//trim(process_path(model%grid%sigma_file))// &
               ' does not exist',GM_FATAL)
       end if
       call write_log('Reading sigma file: '//process_path(model%grid%sigma_file))
       if (main_task) then
          open(unit,file=process_path(model%grid%sigma_file))
          read(unit,'(f9.7)',err=10,end=10) (model%grid%sigma(up), up=1,upn)
          close(unit)
       end if
       call broadcast(model%grid%sigma)

    case(SIGMA_METHOD_CONFIG)          ! read from config file

       ! sigma levels have already been read from glissade_read_sigma

       call write_log('Getting sigma levels from configuration file')

    case(SIGMA_METHOD_UNIFORM)

       do up = 1,upn
          model%grid%sigma(up) = real(up-1,kind=dp) / real(upn-1,kind=dp)
       enddo

       call write_log('Computing evenly spaced sigma levels')

    case(SIGMA_METHOD_PATTYN)

       do up = 1,upn
          if (up == 1) then
             model%grid%sigma(up) = 0.d0
          else if (up == upn) then
             model%grid%sigma(up) = 1.d0
          else
             level = real(up-1,kind=dp) / real(upn-1,kind=dp)
             model%grid%sigma(up) = glissade_calc_sigma_pattyn(level)
          end if
       enddo

       call write_log('Computing Pattyn sigma levels')

    end select

    ! Compute stagsigma (= sigma values at layers midpoints)

    model%grid%stagsigma(1:upn-1) =   &
            (model%grid%sigma(1:upn-1) + model%grid%sigma(2:upn)) / 2.0_dp

    ! Compute stagwbndsigma, adding the boundaries to stagsigma

    model%grid%stagwbndsigma(1:upn-1) = model%grid%stagsigma(1:upn-1)
    model%grid%stagwbndsigma(0) = 0.d0
    model%grid%stagwbndsigma(upn) = 1.d0        

    call print_sigma(model)

    return

10  call write_log('something wrong with sigma coord file',GM_FATAL)
    
  end subroutine glissade_load_sigma

!--------------------------------------------------------------------------------

  function glissade_calc_sigma(x,n)

     implicit none
     real(dp) :: glissade_calc_sigma, x, n
      
     glissade_calc_sigma = (1-(x+1)**(-n)) / (1-2**(-n))

  end function glissade_calc_sigma

!--------------------------------------------------------------------------------

  function glissade_calc_sigma_pattyn(x)

     ! Implements an alternate set of sigma levels that encourages better
     ! convergence for higher-order velocities

     implicit none
     real(dp) :: glissade_calc_sigma_pattyn, x

     glissade_calc_sigma_pattyn =   &
         (-2.5641025641d-4)*(41d0*x)**2+3.5256410256d-2*(41d0*x)-8.0047080075d-13

  end function glissade_calc_sigma_pattyn

!--------------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! grid sizes

  subroutine handle_grid(section, model)
    use glimmer_config
    use glissade_types
    use glimmer_filenames
    implicit none
    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    call GetValue(section,'ewn',model%grid%ewn)
    call GetValue(section,'nsn',model%grid%nsn)
    call GetValue(section,'upn',model%grid%upn)
    call GetValue(section,'dew',model%grid%dew)
    call GetValue(section,'dns',model%grid%dns)
    call GetValue(section,'sigma_file',model%grid%sigma_file)
    call GetValue(section,'nx_block', model%grid%nx_block)
    call GetValue(section,'ny_block', model%grid%ny_block)
    call GetOption(section,'global_bc',global_bc_choices,model%grid%global_bc)

    ! If a sigma file is given, then sigma_method = file.
    ! A warning/error is generated if sigma levels are specified in some other way.
    ! If the file is not in the working directory, look for it in the config file directory.
    if (trim(model%grid%sigma_file) /= '') then
       model%grid%sigma_file = filenames_inputname(model%grid%sigma_file)
       model%grid%sigma_method = SIGMA_METHOD_FILE
    end if

    ! Note: sigma_method is read after sigma_file, so that an explicit sigma_method overrides
    !       the SIGMA_METHOD_FILE setting above (as when the sigma option was read from [options]).
    call GetOption(section,'sigma_method',sigma_method_choices,model%grid%sigma_method)
    call GetLogical(section,'read_lat_lon',model%grid%read_lat_lon)

  end subroutine handle_grid

!--------------------------------------------------------------------------------

  subroutine print_grid(model)

    ! Grid options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    call log_section('grid')
    call log_param('ewn', model%grid%ewn)
    call log_param('nsn', model%grid%nsn)
    call log_param('upn', model%grid%upn)
    call log_param('dew', model%grid%dew, 'm')
    call log_param('dns', model%grid%dns, 'm')
    if (model%numerics%compute_blocks == COMPUTE_BLOCKS_ACTIVE_ONLY .or. &
        model%numerics%compute_blocks == COMPUTE_BLOCKS_INQUIRE) then
       call log_param('nx_block', model%grid%nx_block, 'block size for compute_blocks')
       call log_param('ny_block', model%grid%ny_block, 'block size for compute_blocks')
    endif
    call log_option('global_bc', model%grid%global_bc, global_bc_choices)
    call log_option('sigma_method', model%grid%sigma_method, sigma_method_choices)
    if (len_trim(model%grid%sigma_file) > 0) then
       call log_string('sigma_file', model%grid%sigma_file)
    endif
    call log_flag('read_lat_lon', model%grid%read_lat_lon, &
         'if true, lat and lon are read from the input file and written to restart')

  end subroutine print_grid

!--------------------------------------------------------------------------------


  ! Sigma levels
  subroutine handle_sigma(section, model)

    use glimmer_config
    use glissade_types
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    if (model%grid%sigma_method==SIGMA_METHOD_FILE) then
       call write_log('Sigma levels specified twice - use only'// &
            ' config file or separate file, not both',GM_FATAL)
    else
       call GetValue(section,'sigma_levels',model%grid%sigma,model%grid%upn)
    end if

  end subroutine handle_sigma

!--------------------------------------------------------------------------------

  subroutine print_sigma(model)

    ! Sigma levels (from the [sigma] section, a sigma file, or computed; see sigma_method in [grid])

    use glissade_types
    use glissade_log_config
    implicit none
    type(glissade_global_type)  :: model

    call log_section('sigma')
    call log_list('sigma_levels', model%grid%sigma(1:model%grid%upn), digits=6)

  end subroutine print_sigma

!--------------------------------------------------------------------------------

  ! time
  subroutine handle_time(section, model)
    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    real(dp) :: dt_in      ! dt from the config file, if present
    integer  :: nsteps_in  ! nsteps_per_year from the config file, if present

!       Set dt_input_option to specify whether the dynamic timestep (model%time%tinc)
!        is based on dt (the timestep in years) or nsteps_per_year in the config file.
!       The default (dt_input_option = dt) is to read dt directly,
!        but the nsteps_per_year option allows more flexibility.
!       Typically, we want an integral number of timesteps per year,
!        e.g., if running coupled with the mass balance timestep of 1 year,
!        or if writing output at annual intervals.
!       Note: This used to be called dt_option; renamed to avoid a naming conflict with CESM.
    call GetValue(section,'tstart',model%time%tstart)
    call GetValue(section,'tend',model%time%tend)
    call GetOption(section,'dt_input_option',dt_input_option_choices,model%time%dt_input_option)
    ! Give either dt (dt_input_option = dt) or nsteps_per_year (dt_input_option = nsteps_per_year), not both.
    ! Initialize to impossible values to tell whether each key is in the config file.
    dt_in = -huge(1.0d0)
    nsteps_in = -huge(1)
    call GetValue(section,'dt',dt_in)
    call GetValue(section,'nsteps_per_year',nsteps_in)
    if (model%time%dt_input_option == DT_INPUT_OPTION_NSTEPS_PER_YEAR) then
       if (dt_in /= -huge(1.0d0)) then
          call write_log('Error, dt is set in [time] with dt_input_option = nsteps_per_year; ' // &
               'set nsteps_per_year instead', GM_FATAL)
       endif
    else
       if (nsteps_in /= -huge(1)) then
          call write_log('Error, nsteps_per_year is set in [time] with dt_input_option = dt; set dt instead, '// &
               'or set dt_input_option = nsteps_per_year', GM_FATAL)
       endif
    endif
    if (dt_in /= -huge(1.0d0)) model%time%tinc = dt_in
    if (nsteps_in /= -huge(1)) model%time%nsteps_per_year = nsteps_in
    call GetValue(section,'adaptive_cfl_threshold', model%time%adaptive_cfl_threshold)

    ! If the time step was entered in number of steps per year, then set the timestep in years.
    ! If the time step was entered in years, then set nsteps_per_year.
    if (model%time%dt_input_option == DT_INPUT_OPTION_NSTEPS_PER_YEAR) then
       if (model%time%nsteps_per_year > 0) then
          model%time%tinc = 1.d0 / real(model%time%nsteps_per_year, dp)
       else
          call write_log('Must set nsteps_per_year > 0 with this dt option', GM_FATAL)
       endif
    else   ! input dt (in years) directly
       if (model%time%tinc > 0.0d0) then
          model%time%nsteps_per_year = nint(1.d0/model%time%tinc)
       else
          call write_log('Must set dt > 0.0 with this dt option', GM_FATAL)
       endif
    endif   ! dt_input_option

  end subroutine handle_time
  
!--------------------------------------------------------------------------------

  subroutine print_time(model)

    ! Time options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    call log_section('time')
    call log_param('tstart', model%time%tstart, 'yr')
    call log_param('tend', model%time%tend, 'yr')
    call log_option('dt_input_option', model%time%dt_input_option, dt_input_option_choices)
    if (model%time%dt_input_option == DT_INPUT_OPTION_NSTEPS_PER_YEAR) then
       call log_param('nsteps_per_year', model%time%nsteps_per_year)
       call log_param('dt', model%time%tinc, 'yr; computed as 1/nsteps_per_year')
    else
       call log_param('dt', model%time%tinc, 'yr')
    endif
    call log_param('adaptive_cfl_threshold', model%time%adaptive_cfl_threshold, &
         'advection is subcycled when CFL exceeds this value; 0 = no adaptive subcycling')

  end subroutine print_time

!--------------------------------------------------------------------------------

  subroutine handle_diagnostics(section, model)

    ! Diagnostic options: how often diagnostics are written, the diagnostic point,
    ! and the units and cells used for the global diagnostics.

    use glimmer_config
    use glimmer_log
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    character(len=512) :: val

    !Note: Give either dt_diag (yr) or ndiag (time steps) in the config file, not both.
    !      If dt_diag is given, it is used to compute ndiag. (Output is written every ndiag timesteps.)
    if (ConfigSectionHasValue(section,'dt_diag',val) .and. ConfigSectionHasValue(section,'ndiag',val)) then
       call write_log('Error, dt_diag and ndiag are both set in [diagnostics]; set only one of them', GM_FATAL)
    endif
    call GetValue(section,'idiag',        model%diagnostics%idiag)
    call GetValue(section,'jdiag',        model%diagnostics%jdiag)
    call GetValue(section,'dt_diag',      model%diagnostics%dt_diag)
    call GetValue(section,'ndiag',        model%diagnostics%ndiag)
    call GetOption(section,'dm_dt_diag',   dm_dt_diag_choices, model%diagnostics%dm_dt_diag)
    call GetOption(section,'diag_ice_threshold', diag_ice_threshold_choices, model%diagnostics%diag_ice_threshold)

  end subroutine handle_diagnostics

!--------------------------------------------------------------------------------

  subroutine print_diagnostics(model)

    ! Diagnostic options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    if (model%diagnostics%idiag < 1 .or. model%diagnostics%idiag > model%grid%ewn     &
                                        .or.                                                     &
        model%diagnostics%jdiag < 1 .or. model%diagnostics%jdiag > model%grid%nsn) then
        call write_log('Error, global diagnostic point (idiag, jdiag) is out of bounds', GM_FATAL)
    endif

    call log_section('diagnostics')
    call log_param('idiag', model%diagnostics%idiag, 'global i index of the diagnostic point')
    call log_param('jdiag', model%diagnostics%jdiag, 'global j index of the diagnostic point')
    if (model%diagnostics%dt_diag > 0.d0) then
       call log_param('dt_diag', model%diagnostics%dt_diag, 'yr')
    elseif (model%diagnostics%ndiag > 0) then
       call log_param('ndiag', model%diagnostics%ndiag, 'time steps')
    else
       call log_param('dt_diag', model%diagnostics%dt_diag, 'yr; 0 = no diagnostic output')
    endif
    call log_option('dm_dt_diag', model%diagnostics%dm_dt_diag, dm_dt_diag_choices)
    call log_option('diag_ice_threshold', model%diagnostics%diag_ice_threshold, diag_ice_threshold_choices)

  end subroutine print_diagnostics

!--------------------------------------------------------------------------------

  ! options
  subroutine handle_options(section, model)

    use glimmer_config
    use glissade_types
    use glimmer_log

    implicit none
    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model
    integer :: log_level     ! passed to glimmer_log; not stored in the model

    log_level = GM_levels-GM_ERROR
    call GetValue(section,'log_level',log_level)
    call glimmer_set_msg_level(log_level)

    call GetDycore(section,model%options%dycore)   ! 2 or glissade (see glimmer_config)
    call GetLogical(section,'forcewrite_final', model%options%forcewrite_final)
    call GetOption(section,'restart',restart_choices,model%options%restart)
    call GetLogical(section,'restart_extend_velo',model%options%restart_extend_velo)

  end subroutine handle_options

!--------------------------------------------------------------------------------

  subroutine print_options(model)

    ! General options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    if (model%options%dycore /= DYCORE_GLISSADE) then
       call write_log('Error, the Glissade dycore requires dycore = glissade (or 2)',GM_FATAL)
    end if

    call log_section('options')
    call log_option('dycore', model%options%dycore, dycore_choices)
    call log_option('restart', model%options%restart, restart_choices)
    if (model%options%restart /= NO_RESTART) then
       call log_flag('restart_extend_velo', model%options%restart_extend_velo, &
            'if true, use extended velocity fields for restart')
    endif
    call log_flag('forcewrite_final', model%options%forcewrite_final, &
         'if true, force a write to output files when the run completes')

  end subroutine print_options

!--------------------------------------------------------------------------------

  ! parameters
  subroutine handle_numerics(section, model)

    ! Numerical options and parameters: computational blocks, reproducible sums,
    ! thickness limits, and the periodic offsets used by ISMIP-HOM.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetOption(section,'compute_blocks',    compute_blocks_choices, model%numerics%compute_blocks)
    call GetLogical(section,'reproducible_sums', model%numerics%reproducible_sums)

    ! parameters
    call GetValue(section,'thcklim',            model%numerics%thcklim)
    call GetValue(section,'thcklim_temp',       model%numerics%thcklim_temp)
    call GetValue(section,'thck_gradient_ramp', model%numerics%thck_gradient_ramp)
    call GetValue(section,'periodic_offset_ew', model%numerics%periodic_offset_ew)
    call GetValue(section,'periodic_offset_ns', model%numerics%periodic_offset_ns)

  end subroutine handle_numerics

!--------------------------------------------------------------------------------

  subroutine print_numerics(model)

    ! Numerical options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    if (model%numerics%compute_blocks /= COMPUTE_BLOCKS_ALL) then
       if (model%grid%nx_block < 1 .or. model%grid%ny_block < 1) then
          call write_log('Must set nx_block and ny_block > 0 for this compute_blocks option', GM_FATAL)
       endif
    endif

    !Note: The Glissade dycore is known to crash for thcklim = 0, but has not
    !      been extensively tested for small values of thcklim.
    !      Values smaller than 1 mm may be OK, but no guarantees.
    if (model%numerics%thcklim < 1.d-3) then   ! 1 mm
       call write_log('ice limit (thcklim) is too small for Glissade dycore', GM_FATAL)
    endif

    call log_section('numerics')
    call log_option('compute_blocks', model%numerics%compute_blocks, compute_blocks_choices)
    call log_flag('reproducible_sums', model%numerics%reproducible_sums)
    call log_param('thcklim', model%numerics%thcklim, 'm; thickness limit for dynamically active ice')
    call log_param('thcklim_temp', model%numerics%thcklim_temp, 'm; thickness limit for temperature calculations')
    call log_param('thck_gradient_ramp', model%numerics%thck_gradient_ramp, &
         'm; 0 = no ramp of land-margin thickness gradients')
    if (model%grid%global_bc == GLOBAL_BC_PERIODIC) then
       call log_param('periodic_offset_ew', model%numerics%periodic_offset_ew, 'm')
       call log_param('periodic_offset_ns', model%numerics%periodic_offset_ns, 'm')
    endif

  end subroutine print_numerics

!--------------------------------------------------------------------------------

  subroutine handle_constants(section, model)

    ! Physical constants.
    !Note: The following physical constants have default values in glimmer_physcon.F90.
    !      Some test cases (e.g., MISMIP) specify different values. The default values
    !       can therefore be overridden by the user in the config file.
    !      For coupled CESM runs, however, CISM uses the values in CESM's shr_const_mod,
    !       which cannot be overridden.
    !      These constants are *not* part of the model derived type.
    !      If running multiple instances, the user should either use the default values
    !       or specify identical values in each config file.  Otherwise, the run will use
    !       whatever values are specified in the last config file to be read.

    use glimmer_config
    use glissade_types
    use glimmer_physcon, only: rhoi, rhoo, grav, shci, lhci, trpt
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

#ifndef CCSMCOUPLED
    call GetValue(section,'rhoi', rhoi)
    call GetValue(section,'rhoo', rhoo)
    call GetValue(section,'grav', grav)
    call GetValue(section,'shci', shci)
    call GetValue(section,'lhci', lhci)
    call GetValue(section,'trpt', trpt)
#endif

  end subroutine handle_constants

!--------------------------------------------------------------------------------

  subroutine print_constants(model)

    ! Physical constants
    ! Note: In CESM runs, CISM uses the CESM shared constants, which are written here.

    use glissade_types
    use glimmer_log
    use glissade_log_config
    use glimmer_physcon, only: rhoi, rhoo, lhci, shci, trpt, grav

    implicit none

    type(glissade_global_type)  :: model

    call log_section('constants')
    call log_param('rhoi', rhoi, 'kg/m^3; ice density')
    call log_param('rhoo', rhoo, 'kg/m^3; ocean density')
    call log_param('grav', grav, 'm/s^2; gravitational acceleration')
    call log_param('shci', shci, 'J/kg/K; heat capacity of ice')
    call log_param('lhci', lhci, 'J/kg; latent heat of ice')
    call log_param('trpt', trpt, 'K; triple point of water')

  end subroutine print_constants

!--------------------------------------------------------------------------------

  subroutine handle_initialization(section, model)

    ! Adjustments of the input fields at initialization: thickness, topography and
    ! upper surface.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetLogical(section,'adjust_input_thickness',  model%initialization%adjust_input_thickness)
    call GetLogical(section,'smooth_input_usrf',       model%initialization%smooth_input_usrf)

  end subroutine handle_initialization

!--------------------------------------------------------------------------------

  subroutine print_initialization(model)

    ! Options and parameters for adjusting the input fields at initialization

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    call log_section('initialization')
    call log_flag('adjust_input_thickness', model%initialization%adjust_input_thickness, &
         'if true, input thickness is adjusted based on surface and bed topography')
    call log_flag('smooth_input_usrf', model%initialization%smooth_input_usrf)

  end subroutine print_initialization

!--------------------------------------------------------------------------------

  subroutine handle_velocity(section, model)

    ! Velocity options: Stokes approximation, effective viscosity, and surface gradient.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetOption(section,'stokes_approx',                  stokes_approx_choices, model%velocity%stokes_approx)
    call GetLogical(section,'diva_slope_correction',          model%velocity%diva_slope_correction)
    call GetOption(section,'efvs_method',                    efvs_method_choices, model%velocity%efvs_method)
    call GetOption(section,'surface_gradient_method',        surface_gradient_method_choices, model%velocity%surface_gradient_method)
    call GetOption(section,'surface_gradient_margin_method', surface_gradient_margin_method_choices, model%velocity%surface_gradient_margin_method)

    ! parameters
    call GetValue(section,'efvs_const',    model%velocity%efvs_const)
    call GetValue(section,'effstrain_min', model%velocity%effstrain_min)
    call GetValue(section,'max_slope',     model%velocity%max_slope)

  end subroutine handle_velocity

!--------------------------------------------------------------------------------

  subroutine print_velocity(model)

    ! Velocity options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    if (model%velocity%stokes_approx == STOKES_APPROX_LOCAL_SIA) then
       if (model%grid%global_bc == GLOBAL_BC_NO_PENETRATION) then
          call write_log('Error, cannot use no-penetration BC with local SIA solver', GM_FATAL)
       endif
    endif  ! Glissade local SIA solver

    call log_section('velocity')
    call log_option('stokes_approx', model%velocity%stokes_approx, stokes_approx_choices)
    if (model%velocity%stokes_approx == STOKES_APPROX_DIVA) then
       call log_flag('diva_slope_correction', model%velocity%diva_slope_correction, &
            'if true, DIVA includes a correction for steep slopes')
    endif
    call log_option('efvs_method', model%velocity%efvs_method, efvs_method_choices)
    if (model%velocity%efvs_method == EFVS_METHOD_CONSTANT) then
       call log_param('efvs_const', model%velocity%efvs_const, 'Pa yr')
    elseif (model%velocity%efvs_method == EFVS_METHOD_NONLINEAR) then
       call log_param('effstrain_min', model%velocity%effstrain_min, 'yr^-1')
    end if
    call log_option('surface_gradient_method', model%velocity%surface_gradient_method, surface_gradient_method_choices)
    call log_option('surface_gradient_margin_method', model%velocity%surface_gradient_margin_method, &
         surface_gradient_margin_method_choices)
    call log_param('max_slope', model%velocity%max_slope, 'maximum surface slope')

  end subroutine print_velocity

!--------------------------------------------------------------------------------

  subroutine handle_velocity_solver(section, model)

    ! Velocity solver options: nonlinear iteration, linear solver and preconditioner,
    ! and the assembly of the matrix terms.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! nonlinear iteration
    call GetOption(section,'nonlinear_method',    nonlinear_method_choices, model%velocity_solver%nonlinear_method)
    call GetOption(section,'nonlinear_residual',  nonlinear_residual_choices, model%velocity_solver%nonlinear_residual)
    call GetValue(section,'nonlinear_maxiters',  model%velocity_solver%nonlinear_maxiters)

    ! linear solver
    call GetOption(section,'linear_solver',       linear_solver_choices, model%velocity_solver%linear_solver)
    call GetOption(section,'preconditioner',      preconditioner_choices, model%velocity_solver%preconditioner)
    call GetValue(section,'linear_solve_ncheck', model%velocity_solver%linear_solve_ncheck)
    call GetValue(section,'linear_maxiters',     model%velocity_solver%linear_maxiters)
    call GetValue(section,'linear_tolerance',    model%velocity_solver%linear_tolerance)

    ! matrix assembly
    call GetOption(section,'assemble_beta',       assemble_beta_choices, model%velocity_solver%assemble_beta)
    call GetOption(section,'assemble_taud',       assemble_taud_choices, model%velocity_solver%assemble_taud)
    call GetOption(section,'assemble_bfric',      assemble_bfric_choices, model%velocity_solver%assemble_bfric)
    call GetOption(section,'assemble_lateral',    assemble_lateral_choices, model%velocity_solver%assemble_lateral)

  end subroutine handle_velocity_solver

!--------------------------------------------------------------------------------

  subroutine print_velocity_solver(model)

    ! Velocity solver options

    use glissade_types
    use glimmer_log
    use glissade_log_config
    use glimmer_config, only: option_name
    use cism_parallel, only: tasks

    implicit none

    type(glissade_global_type)  :: model

    ! Forbidden options for running in parallel
    ! Note: The local SIA solver does not use a sparse linear solver, so linear_solver is ignored.
    if (tasks > 1 .and. model%velocity%stokes_approx /= STOKES_APPROX_LOCAL_SIA        .and.  &
                        (model%velocity_solver%linear_solver==LINEAR_SOLVER_SLAP_BICG  .or.  &
                         model%velocity_solver%linear_solver==LINEAR_SOLVER_SLAP_GMRES .or.  &
                         model%velocity_solver%linear_solver==LINEAR_SOLVER_SLAP_PCG) ) then
       call write_log('Error, SLAP solver not supported for more than one processor', GM_FATAL)
    end if

    ! The SIA preconditioner cannot be used for 2D solves
    if ( (model%velocity%stokes_approx == STOKES_APPROX_SSA  .or.  &
          model%velocity%stokes_approx == STOKES_APPROX_L1L2 .or.  &
          model%velocity%stokes_approx == STOKES_APPROX_DIVA .or.  &
          model%velocity%stokes_approx == STOKES_APPROX_HYBRID)    &
                             .and.                            &
          (model%velocity_solver%linear_solver == LINEAR_SOLVER_PCG .or.    &
           model%velocity_solver%linear_solver == LINEAR_SOLVER_PCG_CHRONGEAR) ) then
       if (model%velocity_solver%preconditioner == PRECONDITIONER_SIA) then
          call write_log('Error, preconditioner = sia cannot be used for the 2D solves of stokes_approx = ' // &
               trim(option_name(stokes_approx_choices, model%velocity%stokes_approx)) // &
               '; use diagonal, tridiag_local or tridiag_global', GM_FATAL)
       endif
    endif

    if ((model%velocity_solver%linear_solver == LINEAR_SOLVER_PCG .or.  &
         model%velocity_solver%linear_solver == LINEAR_SOLVER_PCG_CHRONGEAR) ) then
       if (model%numerics%reproducible_sums) then
          if (model%velocity_solver%preconditioner == PRECONDITIONER_TRIDIAG_LOCAL .or. &
              model%velocity_solver%preconditioner == PRECONDITIONER_TRIDIAG_GLOBAL) then
             call write_log ('Tridiagonal preconditioners are not supported with reproducible sums.')
             call write_log ('Please choose a different preconditioner (e.g., diagonal)', GM_FATAL)
          endif
       endif
    endif

    call log_section('velocity_solver')
    call log_option('nonlinear_method', model%velocity_solver%nonlinear_method, nonlinear_method_choices)
    call log_option('nonlinear_residual', model%velocity_solver%nonlinear_residual, nonlinear_residual_choices)
    call log_param('nonlinear_maxiters', model%velocity_solver%nonlinear_maxiters)
    call log_option('linear_solver', model%velocity_solver%linear_solver, linear_solver_choices)
    if ((model%velocity_solver%linear_solver == LINEAR_SOLVER_PCG .or.  &
         model%velocity_solver%linear_solver == LINEAR_SOLVER_PCG_CHRONGEAR) ) then
       call log_option('preconditioner', model%velocity_solver%preconditioner, preconditioner_choices)
    end if
    call log_param('linear_solve_ncheck', model%velocity_solver%linear_solve_ncheck)
    call log_param('linear_maxiters', model%velocity_solver%linear_maxiters)
    call log_param('linear_tolerance', model%velocity_solver%linear_tolerance)
    call log_option('assemble_beta', model%velocity_solver%assemble_beta, assemble_beta_choices)
    call log_option('assemble_taud', model%velocity_solver%assemble_taud, assemble_taud_choices)
    call log_option('assemble_bfric', model%velocity_solver%assemble_bfric, assemble_bfric_choices)
    call log_option('assemble_lateral', model%velocity_solver%assemble_lateral, assemble_lateral_choices)

    ! unsupported resid options
    if (model%velocity_solver%nonlinear_residual == NONLINEAR_RESIDUAL_MAX) then
      call write_log('Residual as max. value of normalized velocity vector update is not currently scientifically supported.  &
           &USE AT YOUR OWN RISK.', GM_WARNING)
    endif
    if (model%velocity_solver%nonlinear_residual == NONLINEAR_RESIDUAL_MAX_NO_UBAS) then
      call write_log('Residual as max. value of normalized velocity vector update with basal velocity omitted is not currently &
           &scientifically supported.  USE AT YOUR OWN RISK.', GM_WARNING)
    endif
    if (model%velocity_solver%nonlinear_residual == NONLINEAR_RESIDUAL_MEAN) then
      call write_log('Residual as mean value of normalized velocity vector update is not currently scientifically supported.  &
           &USE AT YOUR OWN RISK.', GM_WARNING)
    endif

  end subroutine print_velocity_solver

!--------------------------------------------------------------------------------

  subroutine handle_advection(section, model)

    ! Ice advection: the thickness evolution (transport) scheme, vertical remapping,
    ! and the ice age tracer.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    call GetOption(section,'advection_method',     advection_method_choices, model%advection%advection_method)
    call GetOption(section,'vertical_remap_order', vertical_remap_order_choices, model%advection%vertical_remap_order)
    call GetLogical(section,'compute_ice_age',     model%advection%compute_ice_age)

  end subroutine handle_advection

!--------------------------------------------------------------------------------

  subroutine print_advection(model)

    ! Advection options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    call log_section('advection')
    call log_option('advection_method', model%advection%advection_method, advection_method_choices)
    call log_option('vertical_remap_order', model%advection%vertical_remap_order, vertical_remap_order_choices)
    call log_flag('compute_ice_age', model%advection%compute_ice_age)

  end subroutine print_advection

!--------------------------------------------------------------------------------

  subroutine handle_temperature(section, model)

    ! Ice temperature: the temperature (or enthalpy) scheme, its initialization,
    ! the dissipation term, the order of the thermal and transport solves, and
    ! the pressure-melting-point parameters.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetOption(section,'temperature_method', temperature_method_choices, model%temperature%temperature_method)
    call GetOption(section,'temperature_init',   temperature_init_choices,   model%temperature%temperature_init)
    call GetOption(section,'dissipation_method', dissipation_method_choices, model%temperature%dissipation_method)

    ! parameters
    call GetValue(section,'pmp_offset',    model%temperature%pmp_offset)
    call GetValue(section,'pmp_threshold', model%temperature%pmp_threshold)

  end subroutine handle_temperature

!--------------------------------------------------------------------------------

  subroutine print_temperature(model)

    ! Temperature options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    call log_section('temperature')
    call log_option('temperature_method', model%temperature%temperature_method, temperature_method_choices)
    ! Note: If temperature is read from an input or restart file, the temperature_init option is overridden.
    !       The method actually used is written to the log file by glissade_init_temp.
    call log_option('temperature_init', model%temperature%temperature_init, temperature_init_choices, &
         'not used if temperature is read from a file')
    call log_option('dissipation_method', model%temperature%dissipation_method, dissipation_method_choices)
    call log_param('pmp_threshold', model%temperature%pmp_threshold, 'K')
    call log_param('pmp_offset', model%temperature%pmp_offset, 'K')

    ! unsupported temperature options
    if (model%temperature%temperature_method == TEMPERATURE_METHOD_ENTHALPY) then
      call write_log('WARNING: Enthalpy-based formulation for solving temperature evolution is not currently &
           &scientifically supported. USE AT YOUR OWN RISK.', GM_WARNING)
    endif

  end subroutine print_temperature

!--------------------------------------------------------------------------------

  subroutine handle_flow_law(section, model)

    ! Flow law: how the rate factor flwa is computed, its enhancement factors,
    ! and the exponent in Glen's flow law.

    use glimmer_config
    use glissade_types
    use glimmer_physcon, only: n_glen
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetOption(section,'flow_law',             flow_law_choices, model%flow_law%flow_law)
    call GetLogical(section,'damage_flwa_feedback', model%flow_law%damage_flwa_feedback)

    ! parameters
    call GetValue(section,'n_glen',             n_glen)
    call GetValue(section,'flow_factor_ground', model%flow_law%flow_factor_ground)
    call GetValue(section,'flow_factor_float',  model%flow_law%flow_factor_float)
    call GetValue(section,'flwa_const',         model%flow_law%flwa_const)

  end subroutine handle_flow_law

!--------------------------------------------------------------------------------

  subroutine print_flow_law(model)

    ! Flow law options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config
    use glimmer_physcon, only: n_glen

    implicit none

    type(glissade_global_type)  :: model

    call log_section('flow_law')
    call log_option('flow_law', model%flow_law%flow_law, flow_law_choices)
    if (model%flow_law%flow_law == FLOW_LAW_CONSTANT) then
       call log_param('flwa_const', model%flow_law%flwa_const, 'Pa^-n yr^-1; constant flow factor')
    end if
    call log_param('n_glen', n_glen, 'exponent in Glen''s flow law')
    call log_param('flow_factor_ground', model%flow_law%flow_factor_ground, 'grounded ice')
    call log_param('flow_factor_float', model%flow_law%flow_factor_float, 'floating ice')
    call log_flag('damage_flwa_feedback', model%flow_law%damage_flwa_feedback)

  end subroutine print_flow_law

!--------------------------------------------------------------------------------

  subroutine handle_basal_physics(section, model)

    ! Basal physics: the sliding law and its friction parameters, beta limits,
    ! and the effective pressure used by the sliding laws.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetOption(section,'sliding_law',               sliding_law_choices, model%basal_physics%sliding_law)
    call GetLogical(section,'use_c_space_factor',        model%basal_physics%use_c_space_factor)
    call GetOption(section,'beta_limit_method',         beta_limit_method_choices, model%basal_physics%beta_limit_method)
    call GetOption(section,'powerlaw_c_method',         powerlaw_c_method_choices, model%basal_physics%powerlaw_c_method)
    call GetOption(section,'coulomb_c_method',          coulomb_c_method_choices, model%basal_physics%coulomb_c_method)
    call GetLogical(section,'elevation_based_coulomb_c', model%basal_physics%elevation_based_coulomb_c)
    call GetOption(section,'effecpress_method',         effecpress_method_choices, model%basal_physics%effecpress_method)

    ! basal traction for the local SIA solver (Glide-style; inverse sense of beta)
    call GetOption(section,'btrc_method', btrc_method_choices, model%basal_physics%btrc_method)
    call GetValue(section,'btrc_const',  model%basal_physics%btrc_const)

    ! beta parameters
    call GetValue(section,'beta_grounded_min', model%basal_physics%beta_grounded_min)
    call GetValue(section,'beta_const',        model%basal_physics%beta_const)
    call GetValue(section,'beta_small',        model%basal_physics%beta_small)
    call GetValue(section,'beta_large',        model%basal_physics%beta_large)

    ! basal friction parameters
    call GetValue(section,'powerlaw_c_const',        model%basal_physics%powerlaw_c_const)
    call GetValue(section,'powerlaw_c_max',          model%basal_physics%powerlaw_c_max)
    call GetValue(section,'powerlaw_c_min',          model%basal_physics%powerlaw_c_min)
    call GetValue(section,'powerlaw_m',              model%basal_physics%powerlaw_m)
    call GetValue(section,'coulomb_c_const',         model%basal_physics%coulomb_c_const)
    call GetValue(section,'coulomb_c_max',           model%basal_physics%coulomb_c_max)
    call GetValue(section,'coulomb_c_min',           model%basal_physics%coulomb_c_min)
    call GetValue(section,'coulomb_c_hi',            model%basal_physics%coulomb_c_hi)
    call GetValue(section,'coulomb_c_lo',            model%basal_physics%coulomb_c_lo)
    call GetValue(section,'coulomb_c_bed_hi',        model%basal_physics%coulomb_c_bed_hi)
    call GetValue(section,'coulomb_c_bed_lo',        model%basal_physics%coulomb_c_bed_lo)
    call GetValue(section,'schoof_n',                model%basal_physics%schoof_n)
    call GetValue(section,'schoof_gamma',            model%basal_physics%schoof_gamma)
    call GetValue(section,'schoof_p',                model%basal_physics%schoof_p)
    call GetValue(section,'beta_powerlaw_umax',      model%basal_physics%beta_powerlaw_umax)
    call GetValue(section,'zoet_iverson_ut',         model%basal_physics%zoet_iverson_ut)
    call GetValue(section,'zoet_iverson_nmax',       model%basal_physics%zoet_iverson_nmax)
    call GetValue(section,'friction_powerlaw_k',     model%basal_physics%friction_powerlaw_k)
    call GetValue(section,'coulomb_bump_max_slope',  model%basal_physics%coulomb_bump_max_slope)
    call GetValue(section,'coulomb_bump_wavelength', model%basal_physics%coulomb_bump_wavelength)
    call GetValue(section,'flwa_basal_const',        model%basal_physics%flwa_basal_const)

    ! effective pressure: penetration of ocean water pressure
    call GetValue(section,'p_ocean_penetration', model%basal_physics%p_ocean_penetration)

    ! pseudo-plastic parameters
    call GetValue(section,'pseudo_plastic_q',      model%basal_physics%pseudo_plastic_q)
    call GetValue(section,'pseudo_plastic_u0',     model%basal_physics%pseudo_plastic_u0)
    !WHL TODO - During the config reorg, the next four should be removed
    call GetValue(section,'pseudo_plastic_phimin', model%basal_physics%pseudo_plastic_phimin)
    call GetValue(section,'pseudo_plastic_phimax', model%basal_physics%pseudo_plastic_phimax)
    call GetValue(section,'pseudo_plastic_bedmin', model%basal_physics%pseudo_plastic_bedmin)
    call GetValue(section,'pseudo_plastic_bedmax', model%basal_physics%pseudo_plastic_bedmax)

  end subroutine handle_basal_physics

!--------------------------------------------------------------------------------

  subroutine print_basal_physics(model)

    ! Basal physics (sliding law and basal friction) options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model
    character(len=500) :: message

    if (model%basal_physics%use_c_space_factor) then
       if (model%basal_physics%sliding_law /= SLIDING_LAW_COULOMB .and.  &
           model%basal_physics%sliding_law /= SLIDING_LAW_SCHOOF .and. &
           model%basal_physics%sliding_law /= SLIDING_LAW_TSAI) then
          call write_log('Error, C_space_factor not supported for this choice of sliding_law', GM_FATAL)
       endif
    endif

    ! Inversion options

    ! coulomb_c inversion not supported with elevation-based coulomb_c
    if (model%basal_physics%elevation_based_coulomb_c) then
       if (model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION .or. &
           model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN) then
          call write_log('Error, inversion not supported with elevation_based_coulomb_c = T', GM_FATAL)
       endif
    endif

    ! Note: Cp inversion supported for the basic power law and the Schoof and Tsai laws
    if (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION .or. &
        model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN) then
       if (model%basal_physics%sliding_law == SLIDING_LAW_POWERLAW .or. &
           model%basal_physics%sliding_law == SLIDING_LAW_SCHOOF .or.  &
           model%basal_physics%sliding_law == SLIDING_LAW_TSAI) then
          ! inversion for Cp is supported
       else
          call write_log('Error, Cp inversion is not supported for this sliding_law; ' // &
               'it is supported for powerlaw, schoof and tsai', GM_FATAL)
       endif
    endif

    ! Note: Cc inversion supported for the Zoet-Iverson law, pseudoplastic, Schoof and Tsai laws
    if (model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION .or. &
        model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN) then
       if (model%basal_physics%sliding_law == SLIDING_LAW_ZOET_IVERSON .or. &
           model%basal_physics%sliding_law == SLIDING_LAW_PSEUDO_PLASTIC .or. &
           model%basal_physics%sliding_law == SLIDING_LAW_SCHOOF .or.  &
           model%basal_physics%sliding_law == SLIDING_LAW_TSAI) then
          ! inversion for Cc is supported
       else
          call write_log('Error, Cc inversion is not supported for this sliding_law; ' // &
               'it is supported for zoet_iverson, pseudo_plastic, schoof and tsai', GM_FATAL)
       endif
    endif

    ! Note: The Schoof and Tsai laws permit inversion for Cc combined with
    !        powerlaw_c_method = from_coulomb_c.
    !       In these cases we invert for Cc, with Cp coming along for the ride.
    !       There is currently no option to invert for Cc and Cp independently.
    if ( (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION .or. &
          model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN) &
          .and. &
         (model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION .or. &
          model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN) ) then
       call write_log('Simultaneous independent inversion for Cp and Cc is not supported', &
            GM_FATAL)
    endif

    ! Make sure the user isn't trying to invert for powerlaw_c and coulomb_c independently.
    ! We can invert for both at once only if the two coefficients have a specified functional relationship.
    if ( (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION .and. &
          model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION) &
                                        .or. &
         (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN .and. &
          model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN) ) then
       call write_log('Cannot invert for powerlaw_c and coulomb_c independently')
       call write_log('Please set powerlaw_c_method = from_coulomb_c', GM_FATAL)
    endif

    call log_section('basal_physics')
    call log_option('sliding_law', model%basal_physics%sliding_law, sliding_law_choices)
    if (model%velocity%stokes_approx == STOKES_APPROX_LOCAL_SIA) then
       call log_option('btrc_method', model%basal_physics%btrc_method, btrc_method_choices)
    endif
    call log_option('beta_limit_method', model%basal_physics%beta_limit_method, beta_limit_method_choices)
    call log_option('powerlaw_c_method', model%basal_physics%powerlaw_c_method, powerlaw_c_method_choices)
    call log_option('coulomb_c_method', model%basal_physics%coulomb_c_method, coulomb_c_method_choices)
    call log_flag('use_c_space_factor', model%basal_physics%use_c_space_factor, 'if true, beta is multiplied by C_space_factor')
    call log_flag('elevation_based_coulomb_c', model%basal_physics%elevation_based_coulomb_c)
    call log_option('effecpress_method', model%basal_physics%effecpress_method, effecpress_method_choices)

    if (model%basal_physics%btrc_method == BTRC_METHOD_CONSTANT      .or.  &
        model%basal_physics%btrc_method == BTRC_METHOD_CONSTANT_BWAT .or.  &
        model%basal_physics%btrc_method == BTRC_METHOD_CONSTANT_BPMP) then
       call log_param('btrc_const', model%basal_physics%btrc_const, 'm/yr/Pa; basal traction parameter for the local SIA solver')
    end if

    if (model%basal_physics%sliding_law == SLIDING_LAW_BETA_CONSTANT) then
       call log_param('beta_const', model%basal_physics%beta_const, 'Pa yr/m')
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_NO_SLIP) then
       call log_param('beta_large', model%basal_physics%beta_large, 'Pa yr/m')
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_BETA_BPMP) then
       call log_param('beta_large', model%basal_physics%beta_large, 'Pa yr/m')
       call log_param('beta_small', model%basal_physics%beta_small, 'Pa yr/m')
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_PSEUDO_PLASTIC_OLD .or.  &
            model%basal_physics%sliding_law == SLIDING_LAW_PSEUDO_PLASTIC) then
       call log_param('pseudo_plastic_q', model%basal_physics%pseudo_plastic_q)
       call log_param('pseudo_plastic_u0', model%basal_physics%pseudo_plastic_u0)
       if (model%basal_physics%sliding_law == SLIDING_LAW_PSEUDO_PLASTIC_OLD) then
          call log_param('pseudo_plastic_phimin', model%basal_physics%pseudo_plastic_phimin, 'deg')
          call log_param('pseudo_plastic_phimax', model%basal_physics%pseudo_plastic_phimax, 'deg')
          call log_param('pseudo_plastic_bedmin', model%basal_physics%pseudo_plastic_bedmin, 'm')
          call log_param('pseudo_plastic_bedmax', model%basal_physics%pseudo_plastic_bedmax, 'm')
       endif
       if (model%velocity_solver%assemble_beta == ASSEMBLE_BETA_STANDARD) then
          call write_log('Local beta assembly is recommended for the pseudo-plastic sliding law; ' // &
               'set assemble_beta = local', GM_WARNING)
       endif
       if (model%velocity_solver%assemble_bfric == ASSEMBLE_BFRIC_STANDARD) then
          call write_log('Local bfric assembly is recommended for the pseudo-plastic sliding law; ' // &
               'set assemble_bfric = local', GM_WARNING)
       endif
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_ZOET_IVERSON) then
       ! Note: The Zoet-Iverson law typically uses a spatially variable coulomb_c.
       !       If so, the value written here is just the initial value.
       call log_param('coulomb_c_const', model%basal_physics%coulomb_c_const)
       call log_param('coulomb_c_max', model%basal_physics%coulomb_c_max)
       call log_param('coulomb_c_min', model%basal_physics%coulomb_c_min)
       call log_param('powerlaw_m', model%basal_physics%powerlaw_m)
       call log_param('zoet_iverson_ut', model%basal_physics%zoet_iverson_ut, 'm/yr')
       call log_param('zoet_iverson_nmax', model%basal_physics%zoet_iverson_nmax, 'Pa')
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_ISMIP_HOM_C) then
       if (model%grid%ewn /= model%grid%nsn) then
          call write_log('Error, must have ewn = nsn for ISMIP-HOM test C', GM_FATAL)
       endif
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_POWERLAW) then
       call log_param('powerlaw_c_const', model%basal_physics%powerlaw_c_const, 'Pa (m/yr)^(-1/m)')
       call log_param('powerlaw_c_max', model%basal_physics%powerlaw_c_max)
       call log_param('powerlaw_c_min', model%basal_physics%powerlaw_c_min)
       call log_param('powerlaw_m', model%basal_physics%powerlaw_m)
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_COULOMB) then
       call log_param('coulomb_c_const', model%basal_physics%coulomb_c_const)
       call log_param('coulomb_bump_max_slope', model%basal_physics%coulomb_bump_max_slope)
       call log_param('coulomb_bump_wavelength', model%basal_physics%coulomb_bump_wavelength, 'm')
       call log_param('flwa_basal_const', model%basal_physics%flwa_basal_const, 'Pa^-3 yr^-1; flow factor at the bed')
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_SCHOOF) then
       ! Note: The Schoof and Tsai laws often invert for a spatially variable powerlaw_c or coulomb_c.
       !       If so, the values written here are just the initial values.
       call log_param('coulomb_c_const', model%basal_physics%coulomb_c_const)
       if (model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION .or. &
           model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN) then
          call log_param('coulomb_c_max', model%basal_physics%coulomb_c_max)
          call log_param('coulomb_c_min', model%basal_physics%coulomb_c_min)
       endif
       if (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_FROM_COULOMB_C) then
          call log_comment('Cp = schoof_gamma * Cc^schoof_p')
          call log_param('schoof_gamma', model%basal_physics%schoof_gamma)
          call log_param('schoof_p', model%basal_physics%schoof_p)
       else
          call log_param('powerlaw_c_const', model%basal_physics%powerlaw_c_const, 'Pa (m/yr)^(-1/m)')
       endif
       if (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION .or. &
           model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN) then
          call log_param('powerlaw_c_max', model%basal_physics%powerlaw_c_max)
          call log_param('powerlaw_c_min', model%basal_physics%powerlaw_c_min)
       endif
       call log_param('powerlaw_m', model%basal_physics%powerlaw_m)
       call log_param('schoof_n', model%basal_physics%schoof_n)
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_TSAI) then
       call log_param('coulomb_c_const', model%basal_physics%coulomb_c_const)
       if (model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION .or. &
           model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN) then
          call log_param('coulomb_c_max', model%basal_physics%coulomb_c_max)
          call log_param('coulomb_c_min', model%basal_physics%coulomb_c_min)
       endif
       if (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_FROM_COULOMB_C) then
          call log_comment('Cp = schoof_gamma * Cc^schoof_p')
          call log_param('schoof_gamma', model%basal_physics%schoof_gamma)
          call log_param('schoof_p', model%basal_physics%schoof_p)
       else
          call log_param('powerlaw_c_const', model%basal_physics%powerlaw_c_const, 'Pa (m/yr)^(-1/m)')
       endif
       if (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION .or. &
           model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN) then
          call log_param('powerlaw_c_max', model%basal_physics%powerlaw_c_max)
          call log_param('powerlaw_c_min', model%basal_physics%powerlaw_c_min)
       endif
       call log_param('powerlaw_m', model%basal_physics%powerlaw_m)
    elseif (model%basal_physics%sliding_law == SLIDING_LAW_POWERLAW_EFFECPRESS) then
       call write_log('Weertman-style power law higher-order basal boundary condition is not currently scientifically &
            &supported.  USE AT YOUR OWN RISK.', GM_WARNING)
       !TODO - Use powerlaw_c instead of friction_powerlaw_k?  Allow p and q to be set in config file instead of hard-wired?
       call log_param('friction_powerlaw_k', model%basal_physics%friction_powerlaw_k, 'roughness parameter for the power-law friction law')
    endif

    ! Coulomb elevation parameters
    if (model%basal_physics%elevation_based_coulomb_c) then
       call log_param('coulomb_c_hi', model%basal_physics%coulomb_c_hi)
       call log_param('coulomb_c_lo', model%basal_physics%coulomb_c_lo)
       call log_param('coulomb_c_bed_hi', model%basal_physics%coulomb_c_bed_hi, 'm')
       call log_param('coulomb_c_bed_lo', model%basal_physics%coulomb_c_bed_lo, 'm')
    endif

    call log_param('beta_powerlaw_umax', model%basal_physics%beta_powerlaw_umax, &
         'm/yr; max ice speed when evaluating beta(u); 0 = no limit')

    if (model%basal_physics%beta_limit_method == BETA_LIMIT_METHOD_ABSOLUTE) then
       call log_param('beta_grounded_min', model%basal_physics%beta_grounded_min, 'Pa yr/m; min beta for grounded ice')
    endif

    ! ocean_p parameters
    !TODO - Move these to basal_hydro?
    call log_param('p_ocean_penetration', model%basal_physics%p_ocean_penetration, &
         'ocean connection reduces the effective pressure; 0 = off')

  end subroutine print_basal_physics

!--------------------------------------------------------------------------------

  subroutine handle_basal_hydro(section, model)

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    ! basal water option (the switch for this section)
    call GetOption(section, 'bwat_method', bwat_method_choices, model%basal_hydro%bwat_method)

    ! constant basal water
    call GetValue(section, 'bwat_const', model%basal_hydro%bwat_const)

    ! local basal till
    call GetValue(section, 'bwat_till_max', model%basal_hydro%bwat_till_max)
    call GetValue(section, 'till_drainage_rate', model%basal_hydro%till_drainage_rate)

    ! flux routing
    call GetOption(section, 'flux_routing_method', flux_routing_method_choices, model%basal_hydro%flux_routing_method)
    call GetValue(section, 'bwat_source_const', model%basal_hydro%bwat_source_const)
    call GetValue(section, 'btemp_flow_scale', model%basal_hydro%btemp_flow_scale)
    call GetValue(section, 'btemp_freeze_scale', model%basal_hydro%btemp_freeze_scale)

    ! effective pressure options and parameters
    call GetValue(section, 'effecpress_delta',   model%basal_hydro%effecpress_delta)
    call GetValue(section, 'bwat_threshold', model%basal_hydro%bwat_threshold)
    call GetValue(section, 'bpmp_threshold', model%basal_hydro%bpmp_threshold)
    call GetValue(section, 'haf_threshold', model%basal_hydro%haf_threshold)
    call GetOption(section, 'cavity_open_slide', cavity_open_slide_choices, model%basal_hydro%cavity_open_slide)
    call GetOption(section, 'cavity_open_melt', cavity_open_melt_choices, model%basal_hydro%cavity_open_melt)
    call GetValue(section, 'bump_height', model%basal_hydro%bump_height)
    call GetValue(section, 'bump_wavelength', model%basal_hydro%bump_wavelength)
    call GetValue(section, 'sliding_speed_fixed', model%basal_hydro%sliding_speed_fixed)
    call GetValue(section, 'flwa_basal', model%basal_hydro%flwa_basal)
    call GetValue(section, 'cavity_closing_const', model%basal_hydro%cavity_closing_const)

  end subroutine handle_basal_hydro

!--------------------------------------------------------------------------------

  subroutine print_basal_hydro(model)

    ! Basal hydrology options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none
    type(glissade_global_type)  :: model

    ! check the options

    if (model%basal_hydro%bwat_method == BWAT_METHOD_LOCAL_TILL) then
       if (model%basal_physics%effecpress_method /= EFFECPRESS_METHOD_BWAT_BVP) then
          call write_log('Error, bwat_method = local_till requires effecpress_method = bwat_bvp', GM_FATAL)
       endif
    elseif (model%basal_hydro%bwat_method == BWAT_METHOD_FLUX_ROUTING) then
       if (model%numerics%reproducible_sums) then
          if (model%basal_hydro%flux_routing_method /= FLUX_ROUTING_METHOD_D8) then
             call write_log('Error, this flux_routing scheme is not supported with reproducible sums')
             call write_log('Please use D8 flux-routing or turn off reproducible sums', GM_FATAL)
          endif
          if (model%basal_hydro%btemp_freeze_scale > 0.0d0) then
             call write_log('Error, flux-routing does not support refreezing with reproducible sums')
             call write_log('Please set btemp_freeze_scale = 0 or turn off reproducible sums', GM_FATAL)
          endif
       endif  ! reproducible sums
    endif

    ! write the options and parameters

    call log_section('basal_hydro')
    call log_option('bwat_method', model%basal_hydro%bwat_method, bwat_method_choices)

    if (model%basal_hydro%bwat_method == BWAT_METHOD_CONSTANT) then
       call log_param('bwat_const', model%basal_hydro%bwat_const, 'm; constant basal water depth')
    elseif (model%basal_hydro%bwat_method == BWAT_METHOD_LOCAL_TILL) then
       call log_param('bwat_till_max', model%basal_hydro%bwat_till_max, 'm; maximum till water depth')
       call log_param('till_drainage_rate', model%basal_hydro%till_drainage_rate, 'm/yr; till drainage rate')
       call log_param('effecpress_delta', model%basal_hydro%effecpress_delta, &
            'min effective pressure relative to overburden')
    elseif (model%basal_hydro%bwat_method == BWAT_METHOD_FLUX_ROUTING) then
       call log_option('flux_routing_method', model%basal_hydro%flux_routing_method, flux_routing_method_choices)
       call log_param('bwat_source_const', model%basal_hydro%bwat_source_const, 'm/yr; constant melt source at the bed')
       call log_param('btemp_flow_scale', model%basal_hydro%btemp_flow_scale, &
            'deg C; temperature scale for flow around frozen bed; 0 = off')
       call log_param('btemp_freeze_scale', model%basal_hydro%btemp_freeze_scale, &
            'deg C; temperature scale for refreezing at the bed; 0 = off')
       if (model%basal_physics%effecpress_method == EFFECPRESS_METHOD_BWAT) then
          call log_param('effecpress_delta', model%basal_hydro%effecpress_delta, &
               'min effective pressure relative to overburden')
          call log_param('bwat_threshold', model%basal_hydro%bwat_threshold, 'm')
       elseif (model%basal_physics%effecpress_method == EFFECPRESS_METHOD_CAVITY_SHEET) then
          call log_option('cavity_open_slide', model%basal_hydro%cavity_open_slide, cavity_open_slide_choices)
          call log_option('cavity_open_melt', model%basal_hydro%cavity_open_melt, cavity_open_melt_choices)
          if (model%basal_hydro%cavity_open_slide /= CAVITY_OPEN_SLIDE_NONE) then
             call log_param('bump_height', model%basal_hydro%bump_height, 'm')
             call log_param('bump_wavelength', model%basal_hydro%bump_wavelength, 'm')
             if (model%basal_hydro%cavity_open_slide == CAVITY_OPEN_SLIDE_FIXED_UB) then
                call log_param('sliding_speed_fixed', model%basal_hydro%sliding_speed_fixed, 'm/yr')
             endif
          endif   ! opening by sliding
          call log_param('flwa_basal', model%basal_hydro%flwa_basal, 'Pa^-3 yr^-1; flow factor for cavity closing')
          call log_param('cavity_closing_const', model%basal_hydro%cavity_closing_const, 'constant for cavity closing')
       endif   ! cavity-sheet parameterization
    endif   ! bwat_method

    ! other effective pressure options

    if (model%basal_physics%effecpress_method == EFFECPRESS_METHOD_BPMP) then
       call log_param('effecpress_delta', model%basal_hydro%effecpress_delta, &
            'min effective pressure relative to overburden')
       call log_param('bpmp_threshold', model%basal_hydro%bpmp_threshold, 'K')
    endif

    if (model%basal_physics%p_ocean_penetration > 0.0d0) then
       call log_param('haf_threshold', model%basal_hydro%haf_threshold, &
            'm; cap on height above flotation for land-based ice')
    endif

  end subroutine print_basal_hydro

!--------------------------------------------------------------------------------

  subroutine handle_grounding_line(section, model)

    ! Grounding-line options: the grounding-line parameterization (GLP) for basal stress,
    ! the grounded fraction without GLP, subgrid basal melting near the grounding line,
    ! and the flotation function.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    call GetOption(section,'grounding_line_param', grounding_line_param_choices, model%grounding_line%grounding_line_param)
    call GetOption(section,'fground_no_glp',       fground_no_glp_choices, model%grounding_line%fground_no_glp)
    call GetOption(section,'grounding_line_bmlt',  grounding_line_bmlt_choices, model%grounding_line%grounding_line_bmlt)
    call GetOption(section,'flotation_function',   flotation_function_choices, model%grounding_line%flotation_function)

  end subroutine handle_grounding_line

!--------------------------------------------------------------------------------

  subroutine print_grounding_line(model)

    ! Grounding-line options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model
    character(len=500) :: message

    if (model%grounding_line%grounding_line_bmlt == GROUNDING_LINE_BMLT_FLOATING_FRACTION  .or. &
        model%grounding_line%grounding_line_bmlt == GROUNDING_LINE_BMLT_ZERO_GROUNDED) then
       if (model%grounding_line%grounding_line_param /= GROUNDING_LINE_PARAM_DELUXE) then
          call write_log('Error, this setting of grounding_line_bmlt requires grounding_line_param = deluxe', GM_FATAL)
       endif
    endif

    ! write the options

    call log_section('grounding_line')
    call log_option('grounding_line_param', model%grounding_line%grounding_line_param, grounding_line_param_choices)
    if (model%grounding_line%grounding_line_param == GROUNDING_LINE_PARAM_NONE) then
       call log_option('fground_no_glp', model%grounding_line%fground_no_glp, fground_no_glp_choices)
    endif
    call log_option('grounding_line_bmlt', model%grounding_line%grounding_line_bmlt, grounding_line_bmlt_choices)
    call log_option('flotation_function', model%grounding_line%flotation_function, flotation_function_choices)

    if (model%grounding_line%grounding_line_param == GROUNDING_LINE_PARAM_NONE .and. &
        model%grounding_line%flotation_function == FLOTATION_FUNCTION_PATTYN) then
       write(message,*) 'WARNING: Pattyn flotation function with no GLP tends to be unstable'
       call write_log(message, GM_WARNING)
    endif

  end subroutine print_grounding_line

!--------------------------------------------------------------------------------

  subroutine handle_basal_melt(section, model)

    ! Basal melting of floating ice: the bmlt_float scheme and its parameters,
    ! anomaly forcing, and whether basal melt enters the continuity equation.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetOption(section,'bmlt_float_method',                bmlt_float_method_choices, model%basal_melt%bmlt_float_method)
    call GetOption(section,'bmlt_float_thermal_forcing_param', bmlt_float_thermal_forcing_param_choices, model%basal_melt%bmlt_float_thermal_forcing_param)
    call GetLogical(section,'bmlt_float_init',                  model%basal_melt%bmlt_float_init)
    call GetOption(section,'deltaT_ocn_method',                deltaT_ocn_method_choices, model%basal_melt%deltaT_ocn_method)
    call GetLogical(section,'enable_bmlt_anomaly',              model%basal_melt%enable_bmlt_anomaly)
    call GetLogical(section,'basal_mbal_in_continuity', model%basal_melt%basal_mbal_in_continuity)

    ! general parameters
    call GetValue(section,'bmlt_float_factor',      model%basal_melt%bmlt_float_factor)
    call GetValue(section,'bmlt_cavity_h0',         model%basal_melt%bmlt_cavity_h0)
    call GetValue(section,'bmlt_anomaly_tstart',    model%basal_melt%bmlt_anomaly_tstart)
    call GetValue(section,'bmlt_anomaly_timescale', model%basal_melt%bmlt_anomaly_timescale)

    ! MISMIP+ basal melting parameters
    call GetValue(section,'bmlt_float_omega', model%basal_melt%bmlt_float_omega)
    call GetValue(section,'bmlt_float_h0',    model%basal_melt%bmlt_float_h0)
    call GetValue(section,'bmlt_float_z0',    model%basal_melt%bmlt_float_z0)
    call GetValue(section,'bmlt_float_const', model%basal_melt%bmlt_float_const)
    call GetValue(section,'bmlt_float_xlim',  model%basal_melt%bmlt_float_xlim)

    ! depth-dependent basal melting parameters
    call GetValue(section,'bmlt_float_depth_meltmax',  model%basal_melt%bmlt_float_depth_meltmax)
    call GetValue(section,'bmlt_float_depth_frzmax',   model%basal_melt%bmlt_float_depth_frzmax)
    call GetValue(section,'bmlt_float_depth_zmeltmax', model%basal_melt%bmlt_float_depth_zmeltmax)
    call GetValue(section,'bmlt_float_depth_zmelt0',   model%basal_melt%bmlt_float_depth_zmelt0)
    call GetValue(section,'bmlt_float_depth_zfrzmax',  model%basal_melt%bmlt_float_depth_zfrzmax)
    call GetValue(section,'bmlt_float_depth_meltmin',  model%basal_melt%bmlt_float_depth_meltmin)
    call GetValue(section,'bmlt_float_depth_zmeltmin', model%basal_melt%bmlt_float_depth_zmeltmin)
    call GetValue(section,'bmlt_float_gamma0',         model%basal_melt%bmlt_float_gamma0)

  end subroutine handle_basal_melt

!--------------------------------------------------------------------------------

  subroutine print_basal_melt(model)

    ! Basal melting of floating ice: options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model
    character(len=500) :: message

    if (model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_DTHCK_DT) then
       if (model%basal_melt%bmlt_float_thermal_forcing_param /= BMLT_FLOAT_TF_PARAM_ISMIP6_LOCAL .and. &
           model%basal_melt%bmlt_float_thermal_forcing_param /= BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL .and. &
           model%basal_melt%bmlt_float_thermal_forcing_param /= BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL_SLOPE) then
          write(message,*) 'deltaT_ocn_method dthck_dt option supported only for ISMIP6 bmlt_float schemes'
          call write_log(message, GM_FATAL)
       endif
    endif

    ! write the options and parameters

    call log_section('basal_melt')
    call log_option('bmlt_float_method', model%basal_melt%bmlt_float_method, bmlt_float_method_choices)
    if (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_THERMAL_FORCING) then
       call log_option('bmlt_float_thermal_forcing_param', model%basal_melt%bmlt_float_thermal_forcing_param, &
            bmlt_float_thermal_forcing_param_choices)
       call log_param('bmlt_float_gamma0', model%basal_melt%bmlt_float_gamma0, 'm/yr; melt rate coefficient')
       call log_flag('bmlt_float_init', model%basal_melt%bmlt_float_init, &
            'compute bmlt_float at initialization')
    endif
    call log_flag('basal_mbal_in_continuity', model%basal_melt%basal_mbal_in_continuity, &
         'if true, the basal mass balance is applied in the continuity equation')

    if (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_THERMAL_FORCING) then
       call log_option('deltaT_ocn_method', model%basal_melt%deltaT_ocn_method, deltaT_ocn_method_choices)
    elseif (model%basal_melt%deltaT_ocn_method /= DELTAT_OCN_METHOD_NONE) then
       call log_option('deltaT_ocn_method', model%basal_melt%deltaT_ocn_method, deltaT_ocn_method_choices, &
            'ignored; supported only for bmlt_float_method = thermal_forcing')
    endif

    call log_flag('enable_bmlt_anomaly', model%basal_melt%enable_bmlt_anomaly, 'bmlt_float anomaly forcing')
    if (model%basal_melt%enable_bmlt_anomaly) then
       call log_param('bmlt_anomaly_tstart', model%basal_melt%bmlt_anomaly_tstart, 'yr')
       call log_param('bmlt_anomaly_timescale', model%basal_melt%bmlt_anomaly_timescale, 'yr')
    endif

    if (model%basal_melt%bmlt_float_method /= BMLT_FLOAT_METHOD_NONE .and. &
        model%basal_melt%bmlt_float_method /= BMLT_FLOAT_METHOD_MISMIP) then
       call log_param('bmlt_cavity_h0', model%basal_melt%bmlt_cavity_h0, &
            'm; melt is reduced in cavities thinner than this; 0 = off')
    endif

    if (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_FILE) then
       call log_param('bmlt_float_factor', model%basal_melt%bmlt_float_factor, &
            'multiplies bmlt_float from the input file')
    elseif (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_CONSTANT) then
       call log_param('bmlt_float_const', model%basal_melt%bmlt_float_const, 'm/yr')
       call log_param('bmlt_float_xlim', model%basal_melt%bmlt_float_xlim, 'm')
    elseif (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_MISMIP) then
       call log_param('bmlt_float_omega', model%basal_melt%bmlt_float_omega, 'yr^-1')
       call log_param('bmlt_float_h0', model%basal_melt%bmlt_float_h0, 'm')
       call log_param('bmlt_float_z0', model%basal_melt%bmlt_float_z0, 'm')
    elseif (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_DEPTH) then
       call log_param('bmlt_float_depth_meltmax', model%basal_melt%bmlt_float_depth_meltmax, 'm/yr')
       call log_param('bmlt_float_depth_frzmax', model%basal_melt%bmlt_float_depth_frzmax, 'm/yr')
       call log_param('bmlt_float_depth_zmeltmax', model%basal_melt%bmlt_float_depth_zmeltmax, 'm')
       call log_param('bmlt_float_depth_zmelt0', model%basal_melt%bmlt_float_depth_zmelt0, 'm')
       call log_param('bmlt_float_depth_zfrzmax', model%basal_melt%bmlt_float_depth_zfrzmax, 'm')
       call log_param('bmlt_float_depth_meltmin', model%basal_melt%bmlt_float_depth_meltmin, 'm/yr')
       call log_param('bmlt_float_depth_zmeltmin', model%basal_melt%bmlt_float_depth_zmeltmin, 'm')
    endif

  end subroutine print_basal_melt

!--------------------------------------------------------------------------------

  subroutine handle_calving(section, model)

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! calving options
    call GetOption(section,'calving_law',           calving_law_choices, model%calving%calving_law)
    call GetLogical(section,'calve_at_startup', model%calving%calve_at_startup)
    call GetOption(section,'calving_domain',        calving_domain_choices, model%calving%calving_domain)
    call GetLogical(section,'apply_calving_mask',    model%calving%apply_calving_mask)
    call GetLogical(section,'remove_icebergs',       model%calving%remove_icebergs)
    call GetLogical(section,'remove_isthmuses',      model%calving%remove_isthmuses)
    call GetLogical(section,'limit_marine_cliffs',   model%calving%limit_marine_cliffs)
    call GetLogical(section,'subgrid_calving_front', model%calving%subgrid_calving_front)
    call GetOption(section,'calvingmip_domain',     calvingmip_domain_choices, model%calving%calvingmip_domain)

    ! calving parameters
    call GetValue(section,'marine_limit',       model%calving%marine_limit)
    call GetValue(section,'calving_timescale',  model%calving%calving_timescale)
    call GetValue(section,'calving_thck_threshold', model%calving%calving_thck_threshold)
    call GetValue(section,'dthck_dx_cf',        model%calving%dthck_dx_cf)
    call GetValue(section,'thck_effective_min', model%calving%thck_effective_min)
    call GetValue(section,'eigencalving_const', model%calving%eigencalving_const)
    call GetValue(section,'tau_eigenconstant1', model%calving%tau_eigenconstant1)
    call GetValue(section,'tau_eigenconstant2', model%calving%tau_eigenconstant2)
    call GetValue(section,'stress_threshold',   model%calving%stress_threshold)
    call GetValue(section,'lateral_rate_min',   model%calving%lateral_rate_min)
    call GetValue(section,'effec_stress_min',   model%calving%effec_stress_min)
    call GetValue(section,'calving_length_scale', model%calving%calving_length_scale)
    call GetValue(section,'damage_const',       model%calving%damage_const)
    call GetValue(section,'taumax_cliff',       model%calving%taumax_cliff)
    call GetValue(section,'cliff_timescale',    model%calving%cliff_timescale)
    call GetValue(section,'calving_front_x',    model%calving%calving_front_x)
    call GetValue(section,'calving_front_y',    model%calving%calving_front_y)
    call GetValue(section,'calving_front_radius',         model%calving%calving_front_radius)
    call GetValue(section,'cf_advance_retreat_amplitude', model%calving%cf_advance_retreat_amplitude)
    call GetValue(section,'cf_advance_retreat_period',    model%calving%cf_advance_retreat_period)
    call GetValue(section,'f_ground_threshold', model%calving%f_ground_threshold)

  end subroutine handle_calving

!--------------------------------------------------------------------------------

  subroutine print_calving(model)

    ! Calving options and parameters.
    ! Note: remove_isthmuses = .true. requires remove_icebergs = .true. (for stability).

    use glissade_types
    use glimmer_log
    use glissade_log_config
    use glimmer_config, only: option_name

    implicit none

    type(glissade_global_type)  :: model
    character(len=500) :: message

    ! for stability, icebergs must be removed when isthmuses are removed
    if (model%calving%remove_isthmuses .and. .not.model%calving%remove_icebergs) then
       call write_log('Error, remove_isthmuses = .true. requires remove_icebergs = .true. (for stability)', GM_FATAL)
    endif

    if (model%calving%calving_law == CALVING_LAW_ADVANCE_RETREAT_RATE    .or.  &
        model%calving%calving_law == CALVING_LAW_THCK_THRESHOLD     .or.  &
        model%calving%calving_law == CALVING_LAW_STRESS             .or.  &
        model%calving%calving_law == CALVING_LAW_STRESS_STOCHASTIC  .or.  &
        model%calving%calving_law == CALVING_LAW_EIGENCALVING              .or.  &
        model%calving%calving_law == CALVING_LAW_DAMAGE) then
       if (.not. model%calving%subgrid_calving_front) then
          call write_log('Error, calving_law = ' // trim(option_name(calving_law_choices, model%calving%calving_law)) // &
               ' requires subgrid_calving_front = .true.', GM_FATAL)
       endif
    endif

    if (model%calving%calving_law == CALVING_LAW_THCK_THRESHOLD) then
       if (model%calving%calving_thck_threshold <= 0.0d0) then
          write(message,*) 'Error, this calving option needs calving_thck_threshold > 0'
          call write_log(message, GM_FATAL)
       endif
    endif

    if (model%calving%subgrid_calving_front) then
       if (.not.model%calving%remove_icebergs) then
          call write_log('Error, set remove_icebergs = .true. for stability when subgrid_calving_front = .true.', GM_FATAL)
       endif
    endif

    ! write the options and parameters

    call log_section('calving')
    call log_option('calving_law', model%calving%calving_law, calving_law_choices)
    call log_flag('calve_at_startup', model%calving%calve_at_startup, &
         'if true, ice calves at initialization')
    call log_option('calving_domain', model%calving%calving_domain, calving_domain_choices)

    ! parameters of the calving law
    if (model%calving%calving_law == CALVING_LAW_RELX_THRESHOLD .or.  &
            model%calving%calving_law == CALVING_LAW_TOPG_THRESHOLD) then
       call log_param('marine_limit', model%calving%marine_limit, 'm; marine depth limit')
    elseif (model%calving%calving_law == CALVING_LAW_EIGENCALVING) then
       call log_param('eigencalving_const', model%calving%eigencalving_const, 'm')
    elseif (model%calving%calving_law == CALVING_LAW_STRESS) then
       call log_param('tau_eigenconstant1', model%calving%tau_eigenconstant1)
       call log_param('tau_eigenconstant2', model%calving%tau_eigenconstant2)
       call log_param('stress_threshold', model%calving%stress_threshold, 'Pa')
       call log_param('lateral_rate_min', model%calving%lateral_rate_min, 'm/yr')
    elseif (model%calving%calving_law == CALVING_LAW_STRESS_STOCHASTIC) then
       call log_param('tau_eigenconstant1', model%calving%tau_eigenconstant1)
       call log_param('tau_eigenconstant2', model%calving%tau_eigenconstant2)
       call log_param('stress_threshold', model%calving%stress_threshold, 'Pa')
       call log_param('effec_stress_min', model%calving%effec_stress_min)
       call log_param('calving_length_scale', model%calving%calving_length_scale, 'm')
    elseif (model%calving%calving_law == CALVING_LAW_DAMAGE) then
       call log_param('tau_eigenconstant1', model%calving%tau_eigenconstant1)
       call log_param('tau_eigenconstant2', model%calving%tau_eigenconstant2)
       call log_param('damage_const', model%calving%damage_const, 'Pa yr')
    elseif (model%calving%calving_law == CALVING_LAW_ADVANCE_RETREAT_RATE) then
       call log_param('cf_advance_retreat_amplitude', model%calving%cf_advance_retreat_amplitude, 'm/yr')
       call log_param('cf_advance_retreat_period', model%calving%cf_advance_retreat_period, 'yr')
    elseif (model%calving%calving_law == CALVING_LAW_GRID_MASK) then
       call log_param('calving_front_x', model%calving%calving_front_x, 'm; 0 = not used')
       call log_param('calving_front_y', model%calving%calving_front_y, 'm; 0 = not used')
       call log_param('calving_front_radius', model%calving%calving_front_radius, 'm; 0 = not used')
    endif
    if (model%calving%calving_law /= CALVING_LAW_NONE) then
       call log_param('calving_timescale', model%calving%calving_timescale, 'yr; 0 = calve immediately')
    endif
    ! The threshold thickness is used by the thickness-threshold law, by stress-based calving
    !  (if > thcklim), and by mask-based calving without a subgrid calving front when
    !  calving_timescale > 1 yr (the minimum thinning rate is calving_thck_threshold/calving_timescale).
    if (model%calving%calving_law == CALVING_LAW_THCK_THRESHOLD .or. &
        model%calving%calving_law == CALVING_LAW_STRESS .or. &
        ((model%calving%calving_law == CALVING_LAW_GRID_MASK .or. model%calving%apply_calving_mask) .and. &
         .not. model%calving%subgrid_calving_front .and. &
         model%calving%calving_timescale > 1.0d0)) then
       call log_param('calving_thck_threshold', model%calving%calving_thck_threshold, 'm')
    endif

    call log_flag('subgrid_calving_front', model%calving%subgrid_calving_front)
    if (model%calving%subgrid_calving_front) then
       call log_param('dthck_dx_cf', model%calving%dthck_dx_cf)
       call log_param('thck_effective_min', model%calving%thck_effective_min, 'm')
    endif
    if (model%calving%calvingmip_domain /= CALVINGMIP_DOMAIN_NONE) then
       call log_option('calvingmip_domain', model%calving%calvingmip_domain, calvingmip_domain_choices)
    end if

    call log_flag('apply_calving_mask', model%calving%apply_calving_mask)
    call log_flag('remove_icebergs', model%calving%remove_icebergs)
    call log_flag('remove_isthmuses', model%calving%remove_isthmuses)
    if (model%calving%remove_icebergs .or. model%calving%remove_isthmuses) then
       call log_param('f_ground_threshold', model%calving%f_ground_threshold)
    endif
    call log_flag('limit_marine_cliffs', model%calving%limit_marine_cliffs)
    if (model%calving%limit_marine_cliffs) then
       call log_param('taumax_cliff', model%calving%taumax_cliff, 'Pa')
       call log_param('cliff_timescale', model%calving%cliff_timescale, 'yr')
    endif

  end subroutine print_calving

!--------------------------------------------------------------------------------

  subroutine handle_lateral_melt(section, model)

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    call GetOption(section, 'lateral_melt_method', lateral_melt_method_choices, model%lateral_melt%lateral_melt_method)
    call GetValue(section, 'melt_rate_const', model%lateral_melt%melt_rate_const)
    call GetValue(section, 'melt_factor', model%lateral_melt%melt_factor)
    call GetLogical(section, 'runoff_from_ablation', model%lateral_melt%runoff_from_ablation)
    call GetLogical(section, 'thermal_forcing_avg_3d_to_2d', model%lateral_melt%thermal_forcing_avg_3d_to_2d)
    call GetValue(section, 'ztop_tfavg', model%lateral_melt%ztop_tfavg)
    call GetValue(section, 'zbot_tfavg', model%lateral_melt%zbot_tfavg)

  end subroutine handle_lateral_melt

  !--------------------------------------------------------------------------------

  subroutine print_lateral_melt(model)

    ! Lateral melt options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none
    type(glissade_global_type)  :: model

    call log_section('lateral_melt')
    call log_option('lateral_melt_method', model%lateral_melt%lateral_melt_method, lateral_melt_method_choices)

    if (model%lateral_melt%lateral_melt_method == LATERAL_MELT_METHOD_CONSTANT) then
       call log_param('melt_rate_const', model%lateral_melt%melt_rate_const, 'm/yr')
    elseif (model%lateral_melt%lateral_melt_method == LATERAL_MELT_METHOD_ISMIP) then
       call log_param('melt_factor', model%lateral_melt%melt_factor)
       call log_flag('runoff_from_ablation', model%lateral_melt%runoff_from_ablation, &
            'if false, runoff is read from an input file')
       call log_flag('thermal_forcing_avg_3d_to_2d', model%lateral_melt%thermal_forcing_avg_3d_to_2d, &
            'if false, 2D thermal forcing is read from an input file')
       if (model%lateral_melt%thermal_forcing_avg_3d_to_2d) then
          call log_param('ztop_tfavg', model%lateral_melt%ztop_tfavg, 'm; top of the TF averaging range')
          call log_param('zbot_tfavg', model%lateral_melt%zbot_tfavg, 'm; bottom of the TF averaging range')
       endif
    endif

  end subroutine print_lateral_melt

!--------------------------------------------------------------------------------

  subroutine handle_removal(section, model)

    ! Ice removal that is not calving: removal of ice caps, and forced retreat.
    ! The removed ice is tracked separately from the calving flux.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    call GetLogical(section,'remove_ice_caps', model%removal%remove_ice_caps)
    call GetOption(section,'force_retreat',   force_retreat_choices, model%removal%force_retreat)

  end subroutine handle_removal

!--------------------------------------------------------------------------------

  subroutine print_removal(model)

    ! Ice removal options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    call log_section('removal')
    call log_flag('remove_ice_caps', model%removal%remove_ice_caps, 'removed ice goes to the removal flux')
    call log_option('force_retreat', model%removal%force_retreat, force_retreat_choices)
    if (model%removal%force_retreat == FORCE_RETREAT_FLOATING_ICE .and. &
        .not.model%calving%remove_isthmuses) then
       call log_comment('Warning: force_retreat = floating_ice can be unstable when remove_isthmuses = .false.')
    endif

  end subroutine print_removal

!--------------------------------------------------------------------------------

  logical function atm_levels_given(section)

    ! True if [climate] gives the atmosphere levels (nzatm or zatm)

    use glimmer_config
    implicit none
    type(ConfigSection), pointer :: section
    character(len=512) :: val

    atm_levels_given = .false.
    if (.not.associated(section)) return
    atm_levels_given = ConfigSectionHasValue(section,'nzatm',val) .or. ConfigSectionHasValue(section,'zatm',val)

  end function atm_levels_given

!--------------------------------------------------------------------------------

  logical function ocn_levels_given(section)

    ! True if [ocean_data] gives the ocean levels or basins (nzocn, dzocn, zocn or nbasin)

    use glimmer_config
    implicit none
    type(ConfigSection), pointer :: section
    character(len=512) :: val

    ocn_levels_given = .false.
    if (.not.associated(section)) return
    ocn_levels_given = ConfigSectionHasValue(section,'nzocn',val) .or. ConfigSectionHasValue(section,'dzocn',val) .or. &
                       ConfigSectionHasValue(section,'zocn',val)  .or. ConfigSectionHasValue(section,'nbasin',val)

  end function ocn_levels_given

!--------------------------------------------------------------------------------

  subroutine glissade_get_zatm(model,config)

    ! Read the atmosphere levels (nzatm, zatm in [climate]), if present, from the config file.
    ! Called after glissade_readconfig and before glissade_printconfig (see print_climate)

    use glissade_types
    use glimmer_config
    use glimmer_log
    use cism_parallel, only: main_task

    implicit none

    type(glissade_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file

    ! local variables
    type(ConfigSection), pointer :: section
    character(len=512) :: message
    integer :: k

    ! Use the levels from [climate] if nzatm or zatm is given there
    call GetSection(config,section,'climate')

    if (atm_levels_given(section)) then

       call GetValue(section,'nzatm',model%climate%nzatm)

       if (.not.associated(model%climate%zatm)) allocate(model%climate%zatm(model%climate%nzatm))

       if (size(model%climate%zatm) /= model%climate%nzatm) then
          write(message,*) 'size of zatm : ', size(model%climate%zatm)
          call write_log(trim(message))
          call write_log ( 'zatm has the wrong size; must equal nzatm', GM_FATAL)
       end if

       ! Load the zatm levels from the config file.
       ! For example, the config file could have a section like this:
       ! [climate]
       ! nzatm = 10
       ! zatm = 100., 300., 550., 850., 1150., 1450., 1800., 2250., 2750., 3500.

       ! TODO: Is it possible to get the levels from the netCDF file instead?
       ! Note: By convention, zatm increases with increasing k.
       !       If the input zatm levels do not satisfy this criterion, the code aborts.
       call GetValue(section, 'zatm', model%climate%zatm, model%climate%nzatm)
       do k = 2, model%climate%nzatm
          if (model%climate%zatm(k) - model%climate%zatm(k-1) < 1.0d0) then
             write(message,*) 'nzatm, zatm =', model%climate%nzatm, model%climate%zatm(:)
             call write_log(trim(message))
             write(message,*) 'Must have zatm increasing with increasing k in the [climate] section'
             call write_log(message, GM_FATAL)
          endif
       enddo

    elseif (model%climate%smb_input_function == SMB_INPUT_FUNCTION_XYZ .or. &
            model%climate%artm_input_function == ARTM_INPUT_FUNCTION_XYZ) then

       ! Allocate and initialize model%climate%zatm using the standard CESM elevation classes.
       model%climate%nzatm = 10
       if (associated(model%climate%zatm)) deallocate(model%climate%zatm)
       allocate(model%climate%zatm(model%climate%nzatm))
       model%climate%zatm(:) = &
            (/ 100.d0, 300.d0, 550.d0, 850.d0, 1150.d0, 1450.d0, 1800.d0, 2250.d0, 2750.d0, 3500.d0 /)

    else   ! zatm not needed, but allocate and initialize to be safe

       model%climate%nzatm = 1
       if (associated(model%climate%zatm)) deallocate(model%climate%zatm)
       allocate(model%climate%zatm(model%climate%nzatm))
       model%climate%zatm(:) = (/ 0.0d0 /)

    endif   ! atm_levels_given

  end subroutine glissade_get_zatm

!--------------------------------------------------------------------------------

  subroutine handle_climate(section, model)

    ! Surface climate forcing: SMB and surface air temperature (input options, anomalies,
    ! overwrites, and the parameters of the internal SMB functions).

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetOption(section,'smb_input_units',                 smb_input_units_choices, model%climate%smb_input_units)
    call GetOption(section,'smb_input_function',              smb_input_function_choices, model%climate%smb_input_function)
    call GetOption(section,'artm_input_function',             artm_input_function_choices, model%climate%artm_input_function)
    call GetLogical(section,'enable_smb_anomaly',              model%climate%enable_smb_anomaly)
    call GetLogical(section,'enable_artm_anomaly',             model%climate%enable_artm_anomaly)
    call GetLogical(section,'enable_snow_anomaly',             model%climate%enable_snow_anomaly)
    call GetLogical(section,'enable_precip_anomaly',           model%climate%enable_precip_anomaly)
    call GetOption(section,'overwrite_acab',                  overwrite_acab_choices, model%climate%overwrite_acab)
    call GetLogical(section,'enable_acab_dthck_dt_correction', model%climate%enable_acab_dthck_dt_correction)
    call GetLogical(section,'block_inception',                 model%climate%block_inception)

    ! parameters of the internal SMB and artm functions
    call GetValue(section,'t_lapse',            model%climate%t_lapse)
    call GetValue(section,'degree_factor',      model%climate%degree_factor)
    call GetValue(section,'tmlt',               model%climate%tmlt)
    call GetValue(section,'snow_threshold_min', model%climate%snow_threshold_min)
    call GetValue(section,'snow_threshold_max', model%climate%snow_threshold_max)

    ! anomaly and overwrite parameters
    call GetValue(section,'smb_anomaly_tstart',     model%climate%smb_anomaly_tstart)
    call GetValue(section,'smb_anomaly_timescale',  model%climate%smb_anomaly_timescale)
    call GetValue(section,'overwrite_acab_value',   model%climate%overwrite_acab_value)
    call GetValue(section,'overwrite_acab_thck_min', model%climate%overwrite_acab_thck_min)
    call GetValue(section,'artm_anomaly_const',     model%climate%artm_anomaly_const)
    call GetValue(section,'artm_anomaly_tstart',    model%climate%artm_anomaly_tstart)
    call GetValue(section,'artm_anomaly_timescale', model%climate%artm_anomaly_timescale)

  end subroutine handle_climate

!--------------------------------------------------------------------------------

  subroutine print_climate(model, config)

    ! Surface climate (SMB and surface temperature) options and parameters,
    ! and the atmosphere levels (read or set by glissade_get_zatm), written if there is more than one level

    use glissade_types
    use glimmer_config
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: section
    character(len=500) :: message

    ! check the options

    if (model%climate%smb_input_units == SMB_INPUT_UNITS_M_ICE_PER_YR) then
       if (model%climate%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ .or. &
           model%climate%smb_input_function == SMB_INPUT_FUNCTION_XYZ .or. &
           model%climate%enable_smb_anomaly .or. &
           model%glaciers%enable_glaciers) then
          call write_log('Error, this setting requires input SMB in mm/yr w.e.: set smb_input_units = mm_we_per_yr', &
               GM_FATAL)
       endif
    endif

    ! write the options and parameters

    call log_section('climate')
    call log_option('smb_input_units', model%climate%smb_input_units, smb_input_units_choices)
    call log_option('smb_input_function', model%climate%smb_input_function, smb_input_function_choices)
    if (model%climate%smb_input_function == SMB_INPUT_FUNCTION_PDD) then
       call log_param('degree_factor', model%climate%degree_factor, 'mm/yr w.e./deg C')
       call log_param('tmlt', model%climate%tmlt, 'deg C; melt threshold')
       call log_param('snow_threshold_min', model%climate%snow_threshold_min, 'deg C')
       call log_param('snow_threshold_max', model%climate%snow_threshold_max, 'deg C')
    endif
    call log_option('artm_input_function', model%climate%artm_input_function, artm_input_function_choices)
    if (model%climate%artm_input_function == ARTM_INPUT_FUNCTION_XY_LAPSE) then
       call log_param('t_lapse', model%climate%t_lapse, 'deg C/m; artm lapse rate')
    endif

    call log_flag('enable_smb_anomaly', model%climate%enable_smb_anomaly)
    if (model%climate%enable_smb_anomaly) then
       call log_param('smb_anomaly_tstart', model%climate%smb_anomaly_tstart, 'yr')
       call log_param('smb_anomaly_timescale', model%climate%smb_anomaly_timescale, 'yr')
    endif
    call log_flag('enable_artm_anomaly', model%climate%enable_artm_anomaly)
    if (model%climate%enable_artm_anomaly) then
       call log_param('artm_anomaly_const', model%climate%artm_anomaly_const, 'deg C')
       call log_param('artm_anomaly_tstart', model%climate%artm_anomaly_tstart, 'yr')
       call log_param('artm_anomaly_timescale', model%climate%artm_anomaly_timescale, 'yr')
    endif
    call log_flag('enable_snow_anomaly', model%climate%enable_snow_anomaly)
    call log_flag('enable_precip_anomaly', model%climate%enable_precip_anomaly)

    call log_option('overwrite_acab', model%climate%overwrite_acab, overwrite_acab_choices)
    if (model%climate%overwrite_acab /= OVERWRITE_ACAB_NONE) then
       call log_param('overwrite_acab_value', model%climate%overwrite_acab_value, 'm/yr')
       if (model%climate%overwrite_acab == OVERWRITE_ACAB_THCK_MIN) then
          call log_param('overwrite_acab_thck_min', model%climate%overwrite_acab_thck_min, 'm')
       endif
    endif
    call log_flag('enable_acab_dthck_dt_correction', model%climate%enable_acab_dthck_dt_correction, &
         'correct acab based on dthck_dt_obs')
    call log_flag('block_inception', model%climate%block_inception, &
         'block inception outside the main ice sheet')

    if (model%climate%nzatm > 1) then
       call GetSection(config,section,'climate')
       if (.not.atm_levels_given(section)) then
          call log_comment('no nzatm or zatm in [climate]; using the standard CESM elevation classes')
       endif
       call log_param('nzatm', model%climate%nzatm)
       call log_list('zatm', model%climate%zatm, 'm')
    endif

  end subroutine print_climate

!-------------------------------------------------------------------------

  subroutine glissade_get_zocn(model,config)

    ! Read the ocean levels and basins (nzocn, dzocn, zocn, nbasin in [ocean_data]), if present, from the config file.
    ! Called after glissade_readconfig and before glissade_printconfig (see print_ocean_data)

    use glissade_types
    use glimmer_config
    use glimmer_log
    use cism_parallel, only: main_task

    implicit none

    type(glissade_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file

    ! local variables
    type(ConfigSection), pointer :: section
    character(len=512) :: message
    integer :: k

    ! Use the levels from [ocean_data] if nzocn, dzocn, zocn or nbasin is given there
    call GetSection(config,section,'ocean_data')

    if (ocn_levels_given(section)) then

       call GetValue(section,'nzocn',model%ocean_data%nzocn)
       call GetValue(section,'dzocn',model%ocean_data%dzocn)     ! m
       call GetValue(section,'nbasin',model%ocean_data%nbasin)

       if (.not.associated(model%ocean_data%zocn)) &
            allocate(model%ocean_data%zocn(model%ocean_data%nzocn))

       if ( size(model%ocean_data%zocn) /= model%ocean_data%nzocn )then
          write(message,*) 'size of zocn : ', size(model%ocean_data%zocn)
          call write_log(trim(message))
          call write_log ( 'zocn has the wrong size; must equal nzocn', GM_FATAL)
       end if

       ! There are two ways to get zocn levels:
       ! (1) Compute uniform levels based on dzocn.  This is done if dzocn is set to a nonzero value in the config file.
       ! (2) Load zocn levels from the config file.  This is done if dzocn = 0 (the default value).
       ! Note: By convention, k = 1 is the top level, and z becomes more negative with increasing k.
       !       If the input zocn levels do not satisfy this criterion (or if dzocn and zocn are not in the config file),
       !        The code aborts with an error message.

       if (model%ocean_data%dzocn > 0.0d0) then
          do k = 1, model%ocean_data%nzocn
             model%ocean_data%zocn(k) = -model%ocean_data%dzocn * (real(k,dp) - 0.5d0)
          enddo
       else
          call GetValue(section, 'zocn', model%ocean_data%zocn, model%ocean_data%nzocn)
          do k = 2, model%ocean_data%nzocn
             if (model%ocean_data%zocn(k-1) - model%ocean_data%zocn(k) < 1.0d0) then
                write(message,*) 'nzocn, zocn =', model%ocean_data%nzocn, model%ocean_data%zocn(:)
                call write_log(trim(message))
                write(message,*) 'Must have zocn decreasing with increasing k in the [ocean_data] section'
                call write_log(message, GM_FATAL)
             endif
          enddo
       endif

    else    ! no ocean levels in [ocean_data]

       ! Allocate and initialize model%ocean_data%zocn.
       ! Note: Coupled CESM runs require this array to be allocated, even if CISM
       !       is not receiving realistic ocean forcing fields (as of Jan. 2021).
       model%ocean_data%nzocn = get_nzocn(model)
       if (.not.associated(model%ocean_data%zocn)) &
            allocate(model%ocean_data%zocn(model%ocean_data%nzocn))
       model%ocean_data%zocn(:) = 0.0d0

       if (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_THERMAL_FORCING) then
          write(message,*) 'Must set nzocn and zocn (or dzocn) in [ocean_data] to use the thermal forcing option'
          call write_log(message, GM_FATAL)
       endif

    endif   ! ocn_levels_given

  end subroutine glissade_get_zocn

!--------------------------------------------------------------------------------

  subroutine handle_ocean_data(section, model)

    ! Ocean data (thermal forcing) for the thermal-forcing basal melt schemes.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetOption(section,'ocean_data_domain',      ocean_data_domain_choices, model%ocean_data%ocean_data_domain)
    call GetLogical(section,'extrapolate_ocean_data', model%ocean_data%extrapolate_ocean_data)

    ! parameters
    call GetValue(section,'thermal_forcing_basin_min',         model%ocean_data%thermal_forcing_basin_min)
    call GetValue(section,'thermal_forcing_basin_max',         model%ocean_data%thermal_forcing_basin_max)
    call GetValue(section,'thermal_forcing_anomaly',           model%ocean_data%thermal_forcing_anomaly)
    call GetValue(section,'thermal_forcing_anomaly_tstart',    model%ocean_data%thermal_forcing_anomaly_tstart)
    call GetValue(section,'thermal_forcing_anomaly_timescale', model%ocean_data%thermal_forcing_anomaly_timescale)
    call GetValue(section,'thermal_forcing_anomaly_basin',     model%ocean_data%thermal_forcing_anomaly_basin)

  end subroutine handle_ocean_data

!--------------------------------------------------------------------------------

  subroutine print_ocean_data(model, config)

    ! Ocean data (thermal forcing) options and parameters, and the ocean levels and basins.
    ! The options and parameters are written only for the thermal-forcing basal melt scheme, which uses them;
    ! the levels and basins are written when they are given in the config file (see glissade_get_zocn).

    use glissade_types
    use glimmer_config
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: section
    logical :: levels_given

    call GetSection(config,section,'ocean_data')
    levels_given = ocn_levels_given(section)

    if (model%basal_melt%bmlt_float_method /= BMLT_FLOAT_METHOD_THERMAL_FORCING .and. .not.levels_given) return

    call log_section('ocean_data')

    if (levels_given) then
       call log_param('nzocn', model%ocean_data%nzocn)
       if (model%ocean_data%dzocn > 0.0d0) then
          call log_param('dzocn', model%ocean_data%dzocn, 'm; uniform level spacing')
          call log_list('zocn', model%ocean_data%zocn, 'm; computed from dzocn')
       else
          call log_list('zocn', model%ocean_data%zocn, 'm')
       endif
       call log_param('nbasin', model%ocean_data%nbasin)
    endif

    if (model%basal_melt%bmlt_float_method /= BMLT_FLOAT_METHOD_THERMAL_FORCING) return

    call log_option('ocean_data_domain', model%ocean_data%ocean_data_domain, ocean_data_domain_choices)
    if (model%ocean_data%extrapolate_ocean_data) then
       call log_flag('extrapolate_ocean_data', .true., 'CISM extrapolates the ocean data to cavities at run time')
    else
       call log_flag('extrapolate_ocean_data', .false., 'CISM does not extrapolate the ocean data to cavities at run time')
    endif

    if (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_THERMAL_FORCING .and. &
        (model%basal_melt%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL .or. &
         model%basal_melt%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL_SLOPE)) then
       call log_param('thermal_forcing_basin_min', model%ocean_data%thermal_forcing_basin_min, &
            'deg C; min basin-scale TF; 0 = no limit')
       call log_param('thermal_forcing_basin_max', model%ocean_data%thermal_forcing_basin_max, &
            'deg C; max basin-scale TF; 0 = no limit')
    endif   ! nonlocal melt schemes

    call log_param('thermal_forcing_anomaly', model%ocean_data%thermal_forcing_anomaly, 'deg C; 0 = no anomaly')
    if (model%ocean_data%thermal_forcing_anomaly /= 0.0d0) then
       call log_param('thermal_forcing_anomaly_tstart', model%ocean_data%thermal_forcing_anomaly_tstart, 'yr')
       call log_param('thermal_forcing_anomaly_timescale', model%ocean_data%thermal_forcing_anomaly_timescale, 'yr')
       if (model%ocean_data%thermal_forcing_anomaly_basin == 0) then
          call log_param('thermal_forcing_anomaly_basin', model%ocean_data%thermal_forcing_anomaly_basin, &
               'all basins')
       else
          call log_param('thermal_forcing_anomaly_basin', model%ocean_data%thermal_forcing_anomaly_basin)
       endif
    endif

  end subroutine print_ocean_data

!--------------------------------------------------------------------------------

  subroutine handle_geothermal(section, model)

    ! Geothermal heat flux: uniform (geothermal_flux) or read from the input file.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    call GetOption(section,'geothermal_method', geothermal_method_choices, model%geothermal%geothermal_method)
    call GetValue(section,'geothermal_flux',   model%geothermal%geothermal_flux)

  end subroutine handle_geothermal

!--------------------------------------------------------------------------------

  subroutine print_geothermal(model)

    ! Geothermal heat flux option and parameter

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    ! NOTE: Set geothermal_method = file in the config file to read the geothermal heat flux from an input file.
    !       Otherwise it will be overwritten, even if the 'bheatflx' field is present.

    call log_section('geothermal')
    call log_option('geothermal_method', model%geothermal%geothermal_method, geothermal_method_choices)
    call log_param('geothermal_flux', model%geothermal%geothermal_flux, 'W/m^2')

  end subroutine print_geothermal

!--------------------------------------------------------------------------------

  subroutine handle_isostasy(section, model)
    use glimmer_config
    use glissade_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    call GetLogical(section,'compute_isostasy', model%isostasy%compute_isostasy)
    call GetOption(section,'lithosphere_method',lithosphere_method_choices, model%isostasy%lithosphere_method)
    call GetOption(section,'asthenosphere_method',asthenosphere_method_choices, model%isostasy%asthenosphere_method)
    call GetOption(section,'input_topg_method',input_topg_method_choices, model%isostasy%input_topg_method)
    call GetValue(section,'tau_relax_const',model%isostasy%tau_relax_const)
    call GetValue(section,'load_update_interval',model%isostasy%load_update_interval)
    call GetValue(section,'flexural_rigidity',model%isostasy%flexural_rigidity)

  end subroutine handle_isostasy

!--------------------------------------------------------------------------------

  subroutine print_isostasy(model)

    ! Isostasy options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none
    type(glissade_global_type)  :: model

    call log_section('isostasy')
    call log_flag('compute_isostasy', model%isostasy%compute_isostasy)

    if (model%isostasy%compute_isostasy) then
       call log_option('lithosphere_method', model%isostasy%lithosphere_method, lithosphere_method_choices)
       if (model%isostasy%lithosphere_method == LITHOSPHERE_METHOD_ELASTIC) then
          call log_param('flexural_rigidity', model%isostasy%flexural_rigidity, 'N m')
          call log_param('load_update_interval', model%isostasy%load_update_interval, 'yr')
       end if
       call log_option('asthenosphere_method', model%isostasy%asthenosphere_method, asthenosphere_method_choices)
       if (model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_RELAXING) then
          call log_param('tau_relax_const', model%isostasy%tau_relax_const, 'yr; relaxation timescale')
       endif
       call log_option('input_topg_method', model%isostasy%input_topg_method, input_topg_method_choices)
    endif   ! compute isostasy

  end subroutine print_isostasy

!--------------------------------------------------------------------------------

  subroutine handle_inversion(section, model)

    ! Parameters of the inversions for the basal friction parameters (powerlaw_c, coulomb_c)
    ! and for the ocean temperature correction deltaT_ocn.
    ! The inversions themselves are selected in [basal_physics] (powerlaw_c_method,
    ! coulomb_c_method) and [basal_melt] (deltaT_ocn_method).

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! general
    call GetValue(section,'inversion_thck_flotation_buffer', model%inversion%inversion_thck_flotation_buffer)
    call GetValue(section,'inversion_thck_threshold',        model%inversion%inversion_thck_threshold)
    call GetValue(section,'inversion_damping_factor',        model%inversion%inversion_damping_factor)

    ! basal friction inversion
    call GetValue(section,'inversion_babc_timescale',    model%inversion%inversion_babc_timescale)
    call GetValue(section,'inversion_babc_thck_scale',   model%inversion%inversion_babc_thck_scale)
    call GetValue(section,'inversion_babc_length_scale', model%inversion%inversion_babc_length_scale)
    call GetValue(section,'inversion_babc_relax_factor', model%inversion%inversion_babc_relax_factor)

    ! deltaT_ocn inversion
    call GetValue(section,'inversion_deltaT_ocn_timescale',    model%inversion%inversion_deltaT_ocn_timescale)
    call GetValue(section,'inversion_deltaT_ocn_thck_scale',   model%inversion%inversion_deltaT_ocn_thck_scale)
    call GetValue(section,'inversion_deltaT_ocn_length_scale', model%inversion%inversion_deltaT_ocn_length_scale)
    call GetValue(section,'inversion_deltaT_ocn_temp_scale',   model%inversion%inversion_deltaT_ocn_temp_scale)
    call GetValue(section,'inversion_deltaT_basin_relax',      model%inversion%inversion_deltaT_basin_relax)

    ! basin-scale deltaT_ocn inversion
    call GetValue(section,'inversion_basin_flotation_threshold',    model%inversion%inversion_basin_flotation_threshold)
    call GetValue(section,'inversion_basin_mass_correction',        model%inversion%inversion_basin_mass_correction)
    call GetValue(section,'inversion_basin_number_mass_correction', model%inversion%inversion_basin_number_mass_correction)

  end subroutine handle_inversion

!--------------------------------------------------------------------------------

  subroutine print_inversion(model)

    ! Inversion parameters.
    ! The section is written only when an inversion is selected: powerlaw_c_method or
    ! coulomb_c_method in [basal_physics], deltaT_ocn_method in [basal_melt], or
    ! glacier_powerlaw_c_method in [glaciers].

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    logical :: babc_inversion     ! true for an inversion for powerlaw_c or coulomb_c
    logical :: deltaT_inversion   ! true for an inversion for deltaT_ocn
    logical :: glacier_inversion  ! true for an inversion for glacier powerlaw_c

    babc_inversion = (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION .or. &
                      model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN .or. &
                      model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION .or. &
                      model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN)

    deltaT_inversion = (model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_INVERSION .or. &
                        model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_INVERSION_BASIN)

    glacier_inversion = (model%glaciers%enable_glaciers .and. &
                         model%glaciers%glacier_powerlaw_c_method == GLACIER_POWERLAW_C_METHOD_INVERSION)

    if (.not. (babc_inversion .or. deltaT_inversion .or. glacier_inversion)) return

    call log_section('inversion')
    if (babc_inversion .or. deltaT_inversion) then
       call log_param('inversion_thck_flotation_buffer', model%inversion%inversion_thck_flotation_buffer, 'm')
       if (model%inversion%inversion_thck_threshold < model%numerics%thcklim) then
          call log_param('inversion_thck_threshold', model%inversion%inversion_thck_threshold, &
               'm; thcklim (' // real_to_string(model%numerics%thcklim) // ' m) is used, since it is larger')
       else
          call log_param('inversion_thck_threshold', model%inversion%inversion_thck_threshold, 'm')
       endif
       call log_param('inversion_damping_factor', model%inversion%inversion_damping_factor)
    endif

    if (babc_inversion .or. glacier_inversion) then
       call log_param('inversion_babc_thck_scale', model%inversion%inversion_babc_thck_scale, 'm')
       call log_param('inversion_babc_timescale', model%inversion%inversion_babc_timescale, 'yr')
       call log_param('inversion_babc_length_scale', model%inversion%inversion_babc_length_scale, 'm')
       call log_param('inversion_babc_relax_factor', model%inversion%inversion_babc_relax_factor)
    endif

    if (deltaT_inversion) then
       call log_param('inversion_deltaT_ocn_thck_scale', model%inversion%inversion_deltaT_ocn_thck_scale, 'm')
       call log_param('inversion_deltaT_ocn_timescale', model%inversion%inversion_deltaT_ocn_timescale, 'yr')
       call log_param('inversion_deltaT_ocn_temp_scale', model%inversion%inversion_deltaT_ocn_temp_scale, 'deg C')
       if (model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_INVERSION_BASIN) then
          call log_param('inversion_basin_flotation_threshold', model%inversion%inversion_basin_flotation_threshold, 'm')
          call log_param('inversion_deltaT_basin_relax', model%inversion%inversion_deltaT_basin_relax, &
               'deg C; value toward which deltaT_basin relaxes')
          call log_param('inversion_basin_number_mass_correction', model%inversion%inversion_basin_number_mass_correction, &
               '0 = no mass correction')
          call log_param('inversion_basin_mass_correction', model%inversion%inversion_basin_mass_correction, 'Gt')
       else   ! 2D deltaT_ocn inversion
          call log_param('inversion_deltaT_ocn_length_scale', model%inversion%inversion_deltaT_ocn_length_scale, 'm')
       endif
    endif

  end subroutine print_inversion

!--------------------------------------------------------------------------------

  subroutine handle_glaciers(section, model)

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    call GetLogical(section,'enable_glaciers',         model%glaciers%enable_glaciers)
    call GetOption(section,'mu_star_method',          mu_star_method_choices, model%glaciers%mu_star_method)
    call GetOption(section,'alpha_snow_method',       alpha_snow_method_choices, model%glaciers%alpha_snow_method)
    call GetOption(section,'glacier_powerlaw_c_method', glacier_powerlaw_c_method_choices, model%glaciers%glacier_powerlaw_c_method)
    call GetOption(section,'snow_method',             snow_method_choices, model%glaciers%snow_method)
    call GetLogical(section,'scale_area',              model%glaciers%scale_area)
    call GetValue(section,'length_scale_factor',     model%glaciers%length_scale_factor)
    call GetValue(section,'glacier_tmlt',            model%glaciers%glacier_tmlt)
    call GetValue(section,'mu_star_const',           model%glaciers%mu_star_const)
    call GetValue(section,'mu_star_min',             model%glaciers%mu_star_min)
    call GetValue(section,'mu_star_max',             model%glaciers%mu_star_max)
    call GetValue(section,'alpha_snow_const',        model%glaciers%alpha_snow_const)
    call GetValue(section,'alpha_snow_min',          model%glaciers%alpha_snow_min)
    call GetValue(section,'alpha_snow_max',          model%glaciers%alpha_snow_max)
    call GetValue(section,'beta_artm_max',           model%glaciers%beta_artm_max)
    call GetValue(section,'beta_artm_increment',     model%glaciers%beta_artm_increment)
    call GetValue(section,'glacier_snow_threshold_min', model%glaciers%glacier_snow_threshold_min)
    call GetValue(section,'glacier_snow_threshold_max', model%glaciers%glacier_snow_threshold_max)
    call GetValue(section,'baseline_date',           model%glaciers%baseline_date)
    call GetValue(section,'rgi_date',                model%glaciers%rgi_date)
    call GetValue(section,'recent_date',             model%glaciers%recent_date)
    call GetValue(section,'glacier_diag_thck_min',   model%glaciers%glacier_diag_thck_min)
    call GetLogical(section,'redistribute_advanced_ice',  model%glaciers%redistribute_advanced_ice)
    call GetValue(section,'thinning_rate_advanced_ice', model%glaciers%thinning_rate_advanced_ice)
    call GetValue(section,'smb_weight_advanced_ice', model%glaciers%smb_weight_advanced_ice)

  end subroutine handle_glaciers

!--------------------------------------------------------------------------------

  subroutine print_glaciers(model)

    ! Glacier options and parameters.
    ! Note: length_scale_factor is used only with scale_area = .true. (see glissade_glacier_init).
    ! The parameters of the glacier powerlaw_c inversion are written in [inversion].

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none
    type(glissade_global_type)  :: model

    if (model%glaciers%mu_star_method /= MU_STAR_METHOD_INVERSION) then
       if (model%glaciers%alpha_snow_method == ALPHA_SNOW_METHOD_INVERSION) then
          call write_log('Error, must invert for mu_star if inverting for alpha_snow', GM_FATAL)
       elseif (model%glaciers%glacier_powerlaw_c_method == GLACIER_POWERLAW_C_METHOD_INVERSION) then
          call write_log('Error, must invert for mu_star if inverting for powerlaw_c', GM_FATAL)
       endif
    endif

    ! write the options and parameters

    call log_section('glaciers')
    call log_flag('enable_glaciers', model%glaciers%enable_glaciers, 'glacier tracking and tuning')
    call log_option('mu_star_method', model%glaciers%mu_star_method, mu_star_method_choices)
    call log_option('alpha_snow_method', model%glaciers%alpha_snow_method, alpha_snow_method_choices)
    call log_option('glacier_powerlaw_c_method', model%glaciers%glacier_powerlaw_c_method, glacier_powerlaw_c_method_choices)
    call log_option('snow_method', model%glaciers%snow_method, snow_method_choices)
    if (model%glaciers%snow_method == SNOW_METHOD_PRECIP_ARTM) then
       call log_param('glacier_snow_threshold_min', model%glaciers%glacier_snow_threshold_min, 'deg C')
       call log_param('glacier_snow_threshold_max', model%glaciers%glacier_snow_threshold_max, 'deg C')
    endif
    call log_param('glacier_tmlt', model%glaciers%glacier_tmlt, 'deg C')

    call log_flag('scale_area', model%glaciers%scale_area, 'scale glacier area based on latitude')
    if (model%glaciers%scale_area) then
       call log_param('length_scale_factor', model%glaciers%length_scale_factor, 'scales dew and dns')
    elseif (model%glaciers%length_scale_factor /= 1.0d0) then
       call log_param('length_scale_factor', model%glaciers%length_scale_factor, &
            'ignored; used only with scale_area = .true.')
    endif

    call log_param('glacier_diag_thck_min', model%glaciers%glacier_diag_thck_min, 'm')

    call log_flag('redistribute_advanced_ice', model%glaciers%redistribute_advanced_ice, &
         'redistribute advanced ice in the accumulation zone')
    if (model%glaciers%redistribute_advanced_ice) then
       call log_param('thinning_rate_advanced_ice', model%glaciers%thinning_rate_advanced_ice, 'm/yr')
    endif

    if (model%glaciers%mu_star_method == MU_STAR_METHOD_INVERSION) then
       call log_param('smb_weight_advanced_ice', model%glaciers%smb_weight_advanced_ice, &
            'SMB weight in the advanced ablation zone')
       call log_param('mu_star_const', model%glaciers%mu_star_const, 'mm/yr w.e./deg C')
       call log_param('mu_star_min', model%glaciers%mu_star_min, 'mm/yr w.e./deg C')
       call log_param('mu_star_max', model%glaciers%mu_star_max, 'mm/yr w.e./deg C')
       if (model%glaciers%alpha_snow_method == ALPHA_SNOW_METHOD_INVERSION) then
          call log_param('alpha_snow_const', model%glaciers%alpha_snow_const)
          call log_param('alpha_snow_min', model%glaciers%alpha_snow_min)
          call log_param('alpha_snow_max', model%glaciers%alpha_snow_max)
          call log_param('beta_artm_max', model%glaciers%beta_artm_max, 'deg C')
          call log_param('beta_artm_increment', model%glaciers%beta_artm_increment, 'deg C')
          call log_param('baseline_date', model%glaciers%baseline_date, 'yr; glaciers assumed in balance')
          call log_param('rgi_date', model%glaciers%rgi_date, 'yr; date of RGI observations')
          call log_param('recent_date', model%glaciers%recent_date, 'yr; date of recent SMB observations')
       endif
    endif

  end subroutine print_glaciers

!--------------------------------------------------------------------------------

  subroutine define_glissade_restart_variables(model, model_id)

    !> This subroutine analyzes the glissade options input by the user in the config file
    !> and determines which variables are necessary for an exact restart.  MJH 1/11/2013

    ! Please comment thoroughly the reasons why a particular variable needs to be a restart variable for a given config.
    ! Note: This subroutine assumes that any restart variables you add are loadable.
    !       Check glissade_vars.def to make sure any added variables have 'load: 1'

    use glissade_types
    use glissade_io, only: glissade_add_to_restart_variable_list

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glissade_global_type), intent (in) :: model  !> Derived type holding all model info

    integer, intent(in) :: model_id  !> identifier of this ice sheet instance (1 - N, where N is the total number of ice sheet models in this run)

    !------------------------------------------------------------------------------------

    !This was the restart list as of 1/11/13 using the old hot=1 system in glide_vars.def:
    !restart_variable_list=' lat  relx  tauf  thk  thkmask  topg  bheatflx  bmlt_ground  bwat  uvel  vvel  wgrd  flwa  temp  litho_temp  age '

    ! Start with a few variables that we always want - prognostic variables and b.c.
    ! topg - needed to reconstruct all other geometry fields
    ! thk - prognostic variable
    ! temp - prognostic variable
    ! Note: the conversion from temp/flwa to tempstag/flwastag (if necessary) happens in glissade_io.F90
    ! bheatflx, artm, acab - boundary conditions.  Of course if these fields are 0 they don't need 
    !        to be in the restart file, but without adding a check for that we cannot assume any of them are.
    !        There are some options where artm would not be needed.  Logic could be added to make that distinction.
    !        Note that bheatflx may not be an input variable but can also be assigned as a parameter in the config file!
    call glissade_add_to_restart_variable_list('topg thk temp bheatflx artm', model_id)

    ! add the SMB variable, based on model%climate%smb_input_units
    ! Note: SMB can be input in one of several functional forms:
    !       - SMB(x,y)
    !       - SMB(x,y) at reference elevation, plus vertical gradient dSMB/dz(x,y)
    !       - SMB(x,y,z), where z is defined by reference levels
    ! Note: If the SMB field is 'acab', it is assumed to have units of m/y ice
    !       If the SMB field is 'smb', it is assumed to have units of mm/y w.e.

    select case(model%climate%smb_input_function)

       case(SMB_INPUT_FUNCTION_XY)

          select case (model%climate%smb_input_units)
          case (SMB_INPUT_UNITS_M_ICE_PER_YR)
             call glissade_add_to_restart_variable_list('acab', model_id)
          case (SMB_INPUT_UNITS_MM_WE_PER_YR)
             call glissade_add_to_restart_variable_list('smb', model_id)
          end select  ! smb_input_units

       case(SMB_INPUT_FUNCTION_XY_GRADZ)
          call glissade_add_to_restart_variable_list('smb_ref', model_id)
          call glissade_add_to_restart_variable_list('smb_gradz', model_id)
          call glissade_add_to_restart_variable_list('usrf_ref', model_id)

       case(SMB_INPUT_FUNCTION_XYZ)
          call glissade_add_to_restart_variable_list('smb_3d', model_id)

       ! Note: For SMB_INPUT_FUNCTION_PDD, precip and artm are read in as forcing
       !       and are not needed for exact restart

    end select  ! smb_input_function

    ! Similarly for surface temperature (artm), based on model%options%artm_input
    ! Note: These options share usrf_ref and zatm with the SMB options above.

    select case(model%climate%artm_input_function)

       case(ARTM_INPUT_FUNCTION_XY_GRADZ)
          call glissade_add_to_restart_variable_list('artm_ref', model_id)
          call glissade_add_to_restart_variable_list('artm_gradz', model_id)
          if (model%climate%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
             ! usrf_ref was added to restart above; nothing to do here
          else
             call glissade_add_to_restart_variable_list('usrf_ref', model_id)
          endif

       case(ARTM_INPUT_FUNCTION_XYZ)
          call glissade_add_to_restart_variable_list('artm_3d', model_id)
          if (model%climate%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then
             ! smb_levels was added to restart above; nothing to do here
          else
             call glissade_add_to_restart_variable_list('smb_levels', model_id)
          endif

       case(ARTM_INPUT_FUNCTION_XY_LAPSE)
          call glissade_add_to_restart_variable_list('artm_ref', model_id)
          ! Note: Instead of artm_gradz, there is a uniform lapse rate
          if (model%climate%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
             ! usrf_ref was added to restart above; nothing to do here
          else
             call glissade_add_to_restart_variable_list('usrf_ref', model_id)
          endif

    end select  ! artm_input_function

    ! Add anomaly forcing variables
    !TODO - Remove these? Anomaly forcing is typically in a forcing file, not the main input file.

    if (model%climate%enable_smb_anomaly) then
          call glissade_add_to_restart_variable_list('smb_anomaly', model_id)
    endif

    if (model%climate%enable_artm_anomaly) then
       call glissade_add_to_restart_variable_list('artm_anomaly', model_id)
    endif

    if (model%basal_melt%enable_bmlt_anomaly) then
       call glissade_add_to_restart_variable_list('bmlt_float_anomaly', model_id)
    endif

    if (model%grid%read_lat_lon) then
       ! If lat and lon are to be read from the input file, they should be written to the restart file
       call glissade_add_to_restart_variable_list('lat lon', model_id)
    endif

    select case (model%basal_melt%bmlt_float_method)

       ! If bmlt_float is read from an external file at startup, then it needs to be in the restart file
       case (BMLT_FLOAT_METHOD_FILE)
          call glissade_add_to_restart_variable_list('bmlt_float_external', model_id)

       ! If prescribing a warm ocean mask for depth-dependent melting, this needs to be read on restart
       case (BMLT_FLOAT_METHOD_DEPTH)
          call glissade_add_to_restart_variable_list('warm_ocean_mask', model_id)

       case (BMLT_FLOAT_METHOD_THERMAL_FORCING)
          ! Need the latest value of the thermal forcing field.
          ! This could be either the baseline value (if not updating during runtime), or a value read from a forcing file.
          ! If the latter, this field may not be needed, but include to be on the safe side, in case the forcing file
          !  is not read at restart.
          call glissade_add_to_restart_variable_list('thermal_forcing', model_id)

    end select  ! bmlt_float_method

    ! If using an ISMIP6 basin-based melt parameterization, and/or inverting for
    !  basin-scale quantities, we need a 2D field of basin numbers.
    ! Not strictly needed for the ISMIP6_LOCAL option, but included for diagnostics.
    ! The thermal-forcing parameterization matters only when bmlt_float_method uses thermal forcing.
    if ( (model%basal_melt%bmlt_float_method == BMLT_FLOAT_METHOD_THERMAL_FORCING .and.  &
          (model%basal_melt%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_PARAM_ISMIP6_LOCAL .or.  &
           model%basal_melt%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL .or. &
           model%basal_melt%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL_SLOPE) ) .or. &
         model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_INVERSION_BASIN) then
       call glissade_add_to_restart_variable_list('basin_number', model_id)
    endif

    ! add restart variables that depend on the velocity solver and other options
     ! beta - b.c. needed for runs with sliding - could add logic to only include in that case
     ! flwa is not needed for glissade.

      call glissade_add_to_restart_variable_list('kinbcmask bfricflx dissip', model_id)

     ! uvel,vvel: These are needed for an exact restart because we can only recalculate
     !            them to within the picard/jfnk convergence tolerance.
     ! uvel/vvel_extend - These are identical to uvel and vvel, except that the mesh includes
     !                     points along the north and east boundaries of the domain.
     !                    CISM requires these fields for exact restart if the boundary velocities are nonzero,
     !                     as in MISMIP test problems with periodic BCs.
     !                    To output these fields, the user must set restart_extend_velo = .true. in the config file.
     ! Note: It never hurts to write uvel/vvel_extend in place of uvel/vvel. But for most cases where restart
     !       is required (e.g., whole-ice-sheet simulations), velocities are zero along the boundaries
     !       and uvel/vvel are sufficient.

     if (model%options%restart_extend_velo) then
        call glissade_add_to_restart_variable_list('uvel_extend vvel_extend', model_id)
     else
        call glissade_add_to_restart_variable_list('uvel vvel', model_id)
     endif

     ! Glissade approximation options
     select case (model%velocity%stokes_approx)

        case (STOKES_APPROX_DIVA)
           ! DIVA requires the 2D velocity, basal traction and effective viscosity for exact restart.
           ! Since the 2D velocity and basal traction are located on the staggered grid, these fields
           !  must be written to and read from the extended grid if velocites are nonzero at the boundaries.
           !
           ! Note: The 2D velocity is needed if the DIVA scheme solves for the mean velocity.
           !       If DIVA is configured to solve for the velocity at a specific level (e.g., the surface),
           !       then the 2D velocity could instead be copied from the 3D velocity array.
           ! Note: In addition to uvel/vvel_2D, DIVA requires the full 3D velocity field for exact restart,
           !       because horizontal transport is done before updating the velocity.
           
           if (model%options%restart_extend_velo) then
              call glissade_add_to_restart_variable_list('uvel_2d_extend vvel_2d_extend btractx_extend btracty_extend efvs', &
                   model_id)
           else
              call glissade_add_to_restart_variable_list('uvel_2d vvel_2d btractx btracty efvs', model_id)
           endif
           
        case default
           ! Other approximations (including SSA and L1L2) use the 3D uvel and vvel to initialize the velocity

     end select   ! stokes_approx
        
     ! mask used to limit computation to active blocks
     if (model%numerics%compute_blocks == COMPUTE_BLOCKS_ACTIVE_ONLY) then
        call glissade_add_to_restart_variable_list('ice_domain_mask', model_id)
     endif

     ! basal water option for Glissade
     select case (model%basal_hydro%bwat_method)
     case (BWAT_METHOD_NONE, BWAT_METHOD_CONSTANT)
        ! no restart variables needed
     case (BWAT_METHOD_FLUX_ROUTING)
        ! need to know the heat flux associated with refreezing meltwater
        call glissade_add_to_restart_variable_list('bhydroflx', model_id)
     case default
        ! restart needs to know bwat value
        call glissade_add_to_restart_variable_list('bwat', model_id)
     end select

     ! calving options for Glissade

     !TODO: CALVING_LAW_GRID_MASK and apply_calving_mask are redundant; remove one option?
     !      The advantage of apply_calving_mask is that it can be combined with other calving_law options.
     if (model%calving%calving_law == CALVING_LAW_GRID_MASK .or. model%calving%apply_calving_mask .or. &
          model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_INVERSION) then
        if (.not. model%calving%subgrid_calving_front) then
           call glissade_add_to_restart_variable_list('calving_mask', model_id)
        else   ! using a subgrid CF scheme
           call glissade_add_to_restart_variable_list('subgrid_calving_mask', model_id)
        endif
     endif

     ! The eigencalving calculation requires the product of eigenvalues of the horizontal strain rate tensor,
     !  which depends on the stress tensor, which is computed by the HO solver.
     ! On restart, the correct stress and strain rate tensors are not available, so we read in the eigenproduct.
     if (model%calving%calving_law == CALVING_LAW_EIGENCALVING) then
        call glissade_add_to_restart_variable_list('eps_eigen1', model_id)
        call glissade_add_to_restart_variable_list('eps_eigen2', model_id)
     elseif (model%calving%calving_law == CALVING_LAW_STRESS .or. &
             model%calving%calving_law == CALVING_LAW_STRESS_STOCHASTIC .or. &
             model%calving%calving_law == CALVING_LAW_DAMAGE) then
        call glissade_add_to_restart_variable_list('tau_eigen1', model_id)
        call glissade_add_to_restart_variable_list('tau_eigen2', model_id)
     endif

     ! If forcing ice retreat, then we need ice_fraction_retreat_mask (which specifies the cells where retreat is forced)
     !  and reference_thck (which sets up an upper thickness limit for partly retreating cells)
     if (model%removal%force_retreat /= FORCE_RETREAT_NONE) then
        call glissade_add_to_restart_variable_list('ice_fraction_retreat_mask', model_id)
        call glissade_add_to_restart_variable_list('reference_thck', model_id)
     endif

     ! lateral melt options
     if (model%lateral_melt%lateral_melt_method == LATERAL_MELT_METHOD_ISMIP) then
        if (model%lateral_melt%runoff_from_ablation) then
           ! Need input SMB; this is already in the restart file
        else   ! lateral melt parameterization needs runoff_applied (= subglacial discharge)
           call glissade_add_to_restart_variable_list('runoff_applied', model_id)
        endif
        if (model%lateral_melt%thermal_forcing_avg_3d_to_2d) then  ! need 3d thermal_forcing
           call glissade_add_to_restart_variable_list('thermal_forcing', model_id)
        else   ! need vertical-average thermal forcing
           call glissade_add_to_restart_variable_list('thermal_forcing_applied', model_id)
        endif
     endif

     ! other Glissade options

     ! If overwriting acab in certain grid cells, than overwrite_acab_mask needs to be in the restart file.
     ! This mask is read in at model initialization, or is set based on the input acab or ice thickness.
     if (model%climate%overwrite_acab /= OVERWRITE_ACAB_NONE) then
        call glissade_add_to_restart_variable_list('overwrite_acab_mask', model_id)
     endif

    ! ==== Other options ====

    ! internal water option (for enthalpy scheme)
    select case (model%temperature%temperature_method)
      case (TEMPERATURE_METHOD_ENTHALPY)
         ! restart needs to know internal water fraction
         call glissade_add_to_restart_variable_list('waterfrac', model_id)
      case default
         ! no restart variables needed
    end select

    ! basal sliding option
    select case (model%basal_physics%sliding_law)
      case (SLIDING_LAW_COULOMB, SLIDING_LAW_SCHOOF, SLIDING_LAW_TSAI)
         ! Note: These options compute beta internally, so it does not need to be in the restart file.
         if (model%basal_physics%use_c_space_factor) then
            ! c_space_factor needs to be in the restart file
            call glissade_add_to_restart_variable_list('c_space_factor', model_id)
         endif
      case default
         ! Other HO basal boundary conditions may need the external beta field (although there are a few that don't)
         ! Note: If using beta from an external file, then 'beta' here needs to be the fixed, external field,
         !       and not the internal beta field that may have been weighted by the grounded fraction or otherwise adjusted.
         call glissade_add_to_restart_variable_list('beta', model_id)
    end select

    ! basal friction inversion options

    if (model%basal_physics%powerlaw_c_method /= POWERLAW_C_METHOD_CONSTANT .and. &
        model%basal_physics%powerlaw_c_method /= POWERLAW_C_METHOD_FROM_COULOMB_C) then
       call glissade_add_to_restart_variable_list('powerlaw_c', model_id)
    endif

    if (model%basal_physics%coulomb_c_method /= COULOMB_C_METHOD_CONSTANT) then
       call glissade_add_to_restart_variable_list('coulomb_c', model_id)
    endif

    ! If using the basin-scale inversion option for powerlaw_c or coulomb_c, we need a target thickness for grounded ice
    if (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN .or. &
        model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN) then
       call glissade_add_to_restart_variable_list('grounded_thck_target', model_id)
    endif

    ! inversion options for ocean temperature corrections
    if (model%basal_melt%deltaT_ocn_method /= DELTAT_OCN_METHOD_NONE) then
       call glissade_add_to_restart_variable_list('deltaT_ocn', model_id)
    endif

    ! If using the basin-scale inversion option for deltaT_ocn, we need a target thickness for floating ice
    if (model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_INVERSION_BASIN) then
       call glissade_add_to_restart_variable_list('floating_thck_target', model_id)
    endif

    ! fields needed for inversion options that try to match local thickness, elevation, and/or surface speed
    ! Note: If usrf_obs is supplied, thck_obs will be computed at initialization.
    ! Note: Some options do not need usrf_obs and/or velo_sfc_obs as targets,
    !       but they can still be useful diagnostics.
    if (model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION       .or. &
        model%basal_physics%coulomb_c_method  == COULOMB_C_METHOD_INVERSION        .or. &
        model%basal_melt%deltaT_ocn_method    == DELTAT_OCN_METHOD_INVERSION       .or. &
        model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN .or. &
        model%basal_physics%coulomb_c_method  == COULOMB_C_METHOD_INVERSION_BASIN  .or. &
        model%basal_melt%deltaT_ocn_method    == DELTAT_OCN_METHOD_INVERSION_BASIN) then
       call glissade_add_to_restart_variable_list('usrf_obs', model_id)
       call glissade_add_to_restart_variable_list('velo_sfc_obs', model_id)
    endif

    ! fields needed for inversion options that try to match observed dthck_dt
    if (model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_DTHCK_DT .or. &
        model%climate%enable_acab_dthck_dt_correction) then
       call glissade_add_to_restart_variable_list('dthk_dt_obs', model_id)
       call glissade_add_to_restart_variable_list('dthk_dt_obs_basin', model_id)
    endif

    if (model%isostasy%compute_isostasy) then
         ! restart needs to know relaxed topography (the topography to which the mantle would relax with no load)
         ! MJH: I suspect that relx is only needed when asthenosphere=1 (relaxing mantle), but I'm not sure -
         !      this should be tested when isostasy implementation is finalized/tested.
         ! WHL: Looking at subroutine isos_compute, I think relx is needed also when asthenosphere = 0 (fluid mantle).
         !      In this case, topg is set instantaneously to relx - load.
         call glissade_add_to_restart_variable_list('relx', model_id)
         !WHL - The load field also needs to be in the restart file. The reason is that the load is updated
         !      at a period set by isostasy%period. If we restart between two updates, we need to use the most
         !      recently computed load. If we recompute the load right after restarting, the restart may not be exact.
         call glissade_add_to_restart_variable_list('load', model_id)
         if (model%isostasy%asthenosphere_method == ASTHENOSPHERE_METHOD_RELAXING_LATVAR) then
            ! The relaxation timescale is a 2d field read at initialization and again on restart
            call glissade_add_to_restart_variable_list('tau_relax', model_id)
         endif
    endif

    if (model%advection%compute_ice_age) then
       call glissade_add_to_restart_variable_list('ice_age', model_id)
    endif

    if (model%glaciers%enable_glaciers) then
       ! some fields related to glacier indexing
       !TODO - Do we need all the SMB masks?
       call glissade_add_to_restart_variable_list('rgi_glacier_id', model_id)
       call glissade_add_to_restart_variable_list('cism_glacier_id', model_id)
       call glissade_add_to_restart_variable_list('cism_glacier_id_init', model_id)
       call glissade_add_to_restart_variable_list('cism_glacier_id_baseline', model_id)
       call glissade_add_to_restart_variable_list('smb_glacier_id', model_id)
       call glissade_add_to_restart_variable_list('smb_glacier_id_init', model_id)
       call glissade_add_to_restart_variable_list('smb_glacier_id_baseline', model_id)
       call glissade_add_to_restart_variable_list('cism_to_rgi_glacier_id', model_id)
       ! SMB is computed at the end of each year to apply during the next year
       call glissade_add_to_restart_variable_list('smb', model_id)
       call glissade_add_to_restart_variable_list('smb_rgi', model_id)
       call glissade_add_to_restart_variable_list('smb_recent', model_id)
       ! mu_star, alpha_snow, and beta_artm are inversion parameters
       call glissade_add_to_restart_variable_list('glacier_mu_star', model_id)
       call glissade_add_to_restart_variable_list('glacier_alpha_snow', model_id)
       call glissade_add_to_restart_variable_list('glacier_beta_artm', model_id)
       ! smb_obs and usrf_obs are used to invert for mu_star
       call glissade_add_to_restart_variable_list('glacier_smb_obs', model_id)
       call glissade_add_to_restart_variable_list('usrf_obs', model_id)
       ! powerlaw_c is used for power law sliding
       call glissade_add_to_restart_variable_list('powerlaw_c', model_id)
       !TODO: Are area_init and volume_init needed in the restart file?
       !      These could be computed based on cism_glacier_id_init and usrf_obs.
       call glissade_add_to_restart_variable_list('glacier_volume_init', model_id)
       call glissade_add_to_restart_variable_list('glacier_area_init', model_id)
    endif

  end subroutine define_glissade_restart_variables

!--------------------------------------------------------------------------------

end module glissade_setup

!--------------------------------------------------------------------------------
