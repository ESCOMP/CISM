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

  subroutine glissade_config(model,config,fileunit)

    ! Read Glissade configuration from file and print it to the log

    use glissade_types
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    use glimmer_filenames

    implicit none

    type(glissade_global_type), intent(inout) :: model  ! model instance
    type(ConfigSection), pointer  :: config          ! structure holding sections of configuration file
    integer, intent(in), optional :: fileunit        ! fileunit for reading config file 

    type(ConfigSection), pointer :: ncconfig
    integer :: unit

    unit = 99
    if (present(fileunit)) then
       unit = fileunit
    endif

    ! read configuration file
    call glissade_readconfig (model,config)

    ! optionally, read zocn levels from config file, if present, or compute them from dzocn
    call glissade_get_zocn(model,config)

    ! optionally, read zatm levels from config file, if present
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
                            model%numerics%dew, model%numerics%dns)

    ! netCDF I/O
    if (trim(model%funits%ncfile) == '') then
       ncconfig => config
    else
       call ConfigRead(process_path(model%funits%ncfile), ncconfig, unit)
    end if

    call glimmer_nc_readparams(model%funits, model%options%is_restart, ncconfig)

  end subroutine glissade_config

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
    ! Note: The [grid_ocn] and [grid_atm] sections are read by glissade_get_zocn and glissade_get_zatm,
    !       and [sigma] by glissade_read_sigma.

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
    call print_grid_atm(model, config)
    call print_climate(model)
    call print_grid_ocn(model, config)
    call print_ocean_data(model)
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

    model%numerics%dttem = model%numerics%ntem * model%numerics%tinc   

    ! convert dt and dttem from yr to s
    model%numerics%dt     = model%numerics%tinc * scyr
    model%numerics%dttem  = model%numerics%dttem * scyr

    ! allow for subcycling of ice transport
    model%numerics%dt_transport = model%numerics%dt / real(model%numerics%subcyc, dp)

    ! scale calving parameters
    model%calving%timescale = model%calving%timescale * scyr                ! convert from yr to s
    model%calving%cliff_timescale = model%calving%cliff_timescale * scyr    ! convert from yr to s

    ! scale basal traction parameter
    model%paramets%btrc_const = model%paramets%btrc_const/scyr

    ! scale inversion parameters
    model%inversion%babc_timescale = model%inversion%babc_timescale * scyr    ! convert yr to s
    model%inversion%deltaT_ocn_timescale = model%inversion%deltaT_ocn_timescale * scyr   ! yr to s

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
    ! NOTE: The [sigma] section is ignored unless model%options%which_sigma = SIGMA_CONFIG

    if (model%options%which_sigma == SIGMA_CONFIG) then
       call GetSection(config,section,'sigma')
       if (associated(section)) then
          call handle_sigma(section, model)
       else
          model%options%which_sigma = SIGMA_COMPUTE_GLIDE  ! default to standard sigma levels
          call write_log('No [sigma] section present; will compute standard Glide sigma levels')
       end if
    endif

  end subroutine glissade_read_sigma

!-------------------------------------------------------------------------

  subroutine glissade_load_sigma(model,unit)

    ! Compute sigma coordinates or read them from a file
    ! Note: This subroutine is called from glissade_initialise.
    !       If sigma levels are provided in the config file, then they are read
    !        in by glissade_read_sigma, and model%options%which_sigma is set to
    !        SIGMA_CONFIG, in which case this subroutine does nothing.

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

    upn=model%general%upn

    select case(model%options%which_sigma)

    case(SIGMA_COMPUTE_GLIDE)   !  compute standard Glide sigma levels

       do up = 1,upn
          level = real(up-1,kind=dp) / real(upn-1,kind=dp)
          model%numerics%sigma(up) = glissade_calc_sigma(level, 2.d0)
       end do

       call write_log('Computing Glide sigma levels')

    case(SIGMA_EXTERNAL)        ! read from external file

       if (main_task) inquire (exist=there, file=process_path(model%funits%sigfile))
       call broadcast(there)
       if (.not.there) then
          call write_log('Sigma levels file: '//trim(process_path(model%funits%sigfile))// &
               ' does not exist',GM_FATAL)
       end if
       call write_log('Reading sigma file: '//process_path(model%funits%sigfile))
       if (main_task) then
          open(unit,file=process_path(model%funits%sigfile))
          read(unit,'(f9.7)',err=10,end=10) (model%numerics%sigma(up), up=1,upn)
          close(unit)
       end if
       call broadcast(model%numerics%sigma)

    case(SIGMA_CONFIG)          ! read from config file

       ! sigma levels have already been read from glissade_read_sigma

       call write_log('Getting sigma levels from configuration file')

    case(SIGMA_COMPUTE_EVEN)

       do up = 1,upn
          model%numerics%sigma(up) = real(up-1,kind=dp) / real(upn-1,kind=dp)
       enddo

       call write_log('Computing evenly spaced sigma levels')

    case(SIGMA_COMPUTE_PATTYN)

       do up = 1,upn
          if (up == 1) then
             model%numerics%sigma(up) = 0.d0
          else if (up == upn) then
             model%numerics%sigma(up) = 1.d0
          else
             level = real(up-1,kind=dp) / real(upn-1,kind=dp)
             model%numerics%sigma(up) = glissade_calc_sigma_pattyn(level)
          end if
       enddo

       call write_log('Computing Pattyn sigma levels')

    end select

    ! Compute stagsigma (= sigma values at layers midpoints)

    model%numerics%stagsigma(1:upn-1) =   &
            (model%numerics%sigma(1:upn-1) + model%numerics%sigma(2:upn)) / 2.0_dp

    ! Compute stagwbndsigma, adding the boundaries to stagsigma

    model%numerics%stagwbndsigma(1:upn-1) = model%numerics%stagsigma(1:upn-1)
    model%numerics%stagwbndsigma(0) = 0.d0
    model%numerics%stagwbndsigma(upn) = 1.d0        

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

    call GetValue(section,'ewn',model%general%ewn)
    call GetValue(section,'nsn',model%general%nsn)
    call GetValue(section,'upn',model%general%upn)
    call GetValue(section,'dew',model%numerics%dew)
    call GetValue(section,'dns',model%numerics%dns)
    call GetValue(section,'sigma_file',model%funits%sigfile)
    call GetValue(section,'nx_block', model%general%nx_block)
    call GetValue(section,'ny_block', model%general%ny_block)
    call GetValue(section,'global_bc',model%general%global_bc)

    ! We set this flag to one to indicate we've got a sigfile name.
    ! A warning/error is generated if sigma levels are specified in some other way
    ! and mangle the name
    if (trim(model%funits%sigfile) /= '') then
       model%funits%sigfile = filenames_inputname(model%funits%sigfile)
       model%options%which_sigma = SIGMA_EXTERNAL
    end if

    ! Note: sigma_method is read after sigma_file, so that an explicit sigma_method overrides
    !       the SIGMA_EXTERNAL setting above (as when the sigma option was read from [options]).
    call GetValue(section,'sigma_method',model%options%which_sigma)
    call GetValue(section,'read_lat_lon',model%options%read_lat_lon)

  end subroutine handle_grid

!--------------------------------------------------------------------------------

  subroutine print_grid(model)

    ! Grid options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    character(len=*), dimension(0:3), parameter :: global_bc = (/ &
         'periodic boundaries          ', &
         'outflow at boundaries        ', &
         'no penetration at boundaries ', &
         'no ice at boundaries         ' /)

    character(len=*), dimension(0:4), parameter :: sigma_method = (/ &
         'standard Glide sigma levels   ', &
         'sigma levels from sigma_file  ', &
         'sigma levels from [sigma]     ', &
         'evenly spaced sigma levels    ', &
         'Pattyn sigma levels           ' /)

    if (model%general%global_bc < 0 .or. model%general%global_bc >= size(global_bc)) then
       call write_log('Error: Invalid option for global_bc', GM_FATAL)
    endif

    call log_section('grid')
    call log_param('ewn', model%general%ewn)
    call log_param('nsn', model%general%nsn)
    call log_param('upn', model%general%upn)
    call log_param('dew', model%numerics%dew, 'm')
    call log_param('dns', model%numerics%dns, 'm')
    if (model%options%compute_blocks >= ACTIVE_BLOCKS_ONLY) then
       call log_param('nx_block', model%general%nx_block, 'block size for compute_blocks')
       call log_param('ny_block', model%general%ny_block, 'block size for compute_blocks')
    endif
    call log_option('global_bc', model%general%global_bc, global_bc(model%general%global_bc))
    if (model%options%which_sigma >= 0 .and. model%options%which_sigma < size(sigma_method)) then
       call log_option('sigma_method', model%options%which_sigma, sigma_method(model%options%which_sigma))
    else
       call log_option('sigma_method', model%options%which_sigma)
    endif
    if (len_trim(model%funits%sigfile) > 0) then
       call log_string('sigma_file', model%funits%sigfile)
    endif
    call log_flag('read_lat_lon', model%options%read_lat_lon, &
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

    if (model%options%which_sigma==SIGMA_EXTERNAL) then
       call write_log('Sigma levels specified twice - use only'// &
            ' config file or separate file, not both',GM_FATAL)
    else
       call GetValue(section,'sigma_levels',model%numerics%sigma,model%general%upn)
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
    call log_list('sigma_levels', model%numerics%sigma(1:model%general%upn), digits=6)

  end subroutine print_sigma

!--------------------------------------------------------------------------------

  ! time
  subroutine handle_time(section, model)
    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

!       Set dt_input_option to specify whether the dynamic timestep (model%numerics%tinc)
!        is based on dt (the timestep in years) or nsteps_per_year in the config file.
!       The default (dt_input_option = DT_IN_YEARS) is to read dt directly,
!        but the nsteps_per_year option allows more flexibility.
!       Typically, we want an integral number of timesteps per year,
!        e.g., if running coupled with the mass balance timestep of 1 year,
!        or if writing output at annual intervals.
!       Note: This used to be called dt_option; renamed to avoid a naming conflict with CESM.
    call GetValue(section,'tstart',model%numerics%tstart)
    call GetValue(section,'tend',model%numerics%tend)
    call GetValue(section,'dt_input_option',model%options%dt_input_option)
    call GetValue(section,'dt',model%numerics%tinc)
    call GetValue(section,'nsteps_per_year',model%numerics%nsteps_per_year)
    call GetValue(section,'subcyc',model%numerics%subcyc)
    call GetValue(section,'adaptive_cfl_threshold', model%numerics%adaptive_cfl_threshold)
    call GetValue(section,'ntem',model%numerics%ntem)

    ! If the time step was entered in number of steps per year, then set the timestep in years.
    ! If the time step was entered in years, then set nsteps_per_year.
    if (model%options%dt_input_option == DT_STEPS_PER_YEAR) then
       if (model%numerics%nsteps_per_year > 0) then
          model%numerics%tinc = 1.d0 / real(model%numerics%nsteps_per_year, dp)
       else
          call write_log('Must set nsteps_per_year > 0 with this dt option', GM_FATAL)
       endif
    else   ! input dt (in years) directly
       if (model%numerics%tinc > 0.0d0) then
          model%numerics%nsteps_per_year = nint(1.d0/model%numerics%tinc)
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

    character(len=100) :: note

    character(len=*), dimension(0:1), parameter :: dt_input_option = (/ &
         'dt (yr) set in the config file      ', &
         'nsteps_per_year set in config file  ' /)

    if (model%options%dt_input_option < 0 .or. model%options%dt_input_option >= size(dt_input_option)) then
       call write_log('Error, dt_input_option out of range', GM_FATAL)
    end if

    if ( (model%numerics%ntem < 1.0d0) .or. &
       (floor(model%numerics%ntem) /= model%numerics%ntem) ) then
       call write_log('ntem is a multiplier on the basic time step.  It should be a positive integer.  Aborting.',GM_FATAL)
    endif

    ! Note: With adaptive CFL subcycling, subcyc is set to 1; subcycling is adaptive only.
    note = ''
    if (model%numerics%adaptive_cfl_threshold > 0.0d0) then
       if (model%numerics%subcyc /= 1) then
          model%numerics%subcyc = 1
          note = 'reset to 1: subcycling is adaptive only when adaptive_cfl_threshold > 0'
       endif
    endif

    call log_section('time')
    call log_param('tstart', model%numerics%tstart, 'yr')
    call log_param('tend', model%numerics%tend, 'yr')
    call log_option('dt_input_option', model%options%dt_input_option, &
         dt_input_option(model%options%dt_input_option))
    if (model%options%dt_input_option == DT_STEPS_PER_YEAR) then
       call log_param('nsteps_per_year', model%numerics%nsteps_per_year)
       call log_param('dt', model%numerics%tinc, 'yr; computed as 1/nsteps_per_year')
    else
       call log_param('dt', model%numerics%tinc, 'yr')
    endif
    call log_param('ntem', model%numerics%ntem, 'thermal time step / dt')
    call log_param('adaptive_cfl_threshold', model%numerics%adaptive_cfl_threshold, &
         'advection is subcycled when CFL exceeds this value; 0 = no adaptive subcycling')
    call log_param('subcyc', model%numerics%subcyc, note)

  end subroutine print_time

!--------------------------------------------------------------------------------

  subroutine handle_diagnostics(section, model)

    ! Diagnostic options: how often diagnostics are written, the diagnostic point,
    ! and the units and cells used for the global diagnostics.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    !Note: Either dt_diag or ndiag can be specified in the config file.
    !      If dt_diag is specified, it is used to compute ndiag. (Output is written every ndiag timesteps.)
    call GetValue(section,'dt_diag',      model%numerics%dt_diag)
    call GetValue(section,'ndiag',        model%numerics%ndiag)
    call GetValue(section,'idiag',        model%numerics%idiag)
    call GetValue(section,'jdiag',        model%numerics%jdiag)
    call GetValue(section,'dm_dt_diag',   model%options%dm_dt_diag)
    call GetValue(section,'diag_minthck', model%options%diag_minthck)

  end subroutine handle_diagnostics

!--------------------------------------------------------------------------------

  subroutine print_diagnostics(model)

    ! Diagnostic options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    character(len=*), dimension(0:1), parameter :: dm_dt_diag = (/ &
         'write dmass/dt diagnostic in units of kg/s ',  &
         'write dmass/dt diagnostic in units of Gt/yr'/)

    character(len=*), dimension(0:1), parameter :: diag_minthck = (/ &
         'include cells with H > 0 in global diagnostics     ', &
         'include cells with H > thklim in global diagnostics'/)

    if (model%options%dm_dt_diag < 0 .or. model%options%dm_dt_diag >= size(dm_dt_diag)) then
       call write_log('Error, dm_dt_diag out of range',GM_FATAL)
    end if

    if (model%options%diag_minthck < 0 .or. model%options%diag_minthck >= size(diag_minthck)) then
       call write_log('Error, diag_minthck out of range',GM_FATAL)
    end if

    if (model%numerics%idiag < 1 .or. model%numerics%idiag > model%general%ewn     &
                                        .or.                                                     &
        model%numerics%jdiag < 1 .or. model%numerics%jdiag > model%general%nsn) then
        call write_log('Error, global diagnostic point (idiag, jdiag) is out of bounds', GM_FATAL)
    endif

    call log_section('diagnostics')
    if (model%numerics%dt_diag > 0.d0) then
       call log_param('dt_diag', model%numerics%dt_diag, 'yr')
    elseif (model%numerics%ndiag > 0) then
       call log_param('ndiag', model%numerics%ndiag, 'time steps')
    else
       call log_param('dt_diag', model%numerics%dt_diag, 'yr; 0 = no diagnostic output')
    endif
    call log_param('idiag', model%numerics%idiag, 'global i index of the diagnostic point')
    call log_param('jdiag', model%numerics%jdiag, 'global j index of the diagnostic point')
    call log_option('dm_dt_diag', model%options%dm_dt_diag, dm_dt_diag(model%options%dm_dt_diag))
    call log_option('diag_minthck', model%options%diag_minthck, diag_minthck(model%options%diag_minthck))

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
    integer :: loglevel

    loglevel = GM_levels-GM_ERROR
    call GetValue(section,'log_level',loglevel)
    call glimmer_set_msg_level(loglevel)

    call GetDycore(section,model%options%whichdycore)   ! 2 or glissade (see glimmer_config)
    call GetValue(section,'ioparams',model%funits%ncfile)
    call GetValue(section,'forcewrite_final', model%options%forcewrite_final)
    call GetValue(section,'restart',model%options%is_restart)
    call GetValue(section,'restart_extend_velo',model%options%restart_extend_velo)

  end subroutine handle_options

!--------------------------------------------------------------------------------

  subroutine print_options(model)

    ! General options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    character(len=*), dimension(0:2), parameter :: restart = (/ &
         'not a restart                               ', &
         'standard restart                            ', &
         'hybrid restart (model time is reset)        ' /)

    if (model%options%whichdycore /= DYCORE_GLISSADE) then
       call write_log('Error, the Glissade dycore requires dycore = 2 or glissade',GM_FATAL)
    end if
    if (model%options%is_restart < 0 .or. model%options%is_restart >= size(restart)) then
       call write_log('Error, restart option out of range', GM_FATAL)
    end if

    call log_section('options')
    call log_option('dycore', model%options%whichdycore, 'glissade')
    if (len_trim(model%funits%ncfile) > 0) then
       call log_string('ioparams', model%funits%ncfile)
    endif
    call log_option('restart', model%options%is_restart, restart(model%options%is_restart))
    if (model%options%is_restart /= 0) then
       call log_option('restart_extend_velo', model%options%restart_extend_velo, &
            'if 1, use extended velocity fields for restart')
    endif
    call log_flag('forcewrite_final', model%options%forcewrite_final, &
         'if true, force a write to output files when the run completes')

  end subroutine print_options

!--------------------------------------------------------------------------------

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
    call GetValue(section,'compute_blocks',    model%options%compute_blocks)
    call GetValue(section,'reproducible_sums', model%options%reproducible_sums)

    ! parameters
    !TODO - Change 'ice_limit' to 'thklim'?
    call GetValue(section,'ice_limit',          model%numerics%thklim)
    call GetValue(section,'ice_limit_temp',     model%numerics%thklim_temp)
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

    character(len=*), dimension(0:2), parameter :: compute_blocks = (/ &
         'compute on all blocks              ', &
         'compute on active blocks only      ', &
         'inquire the number of active blocks'/)

    if (model%options%compute_blocks < 0 .or. model%options%compute_blocks >= size(compute_blocks)) then
       call write_log('Error, compute_blocks out of range',GM_FATAL)
    end if
    if (model%options%compute_blocks /= 0) then
       if (model%general%nx_block < 1 .or. model%general%ny_block < 1) then
          call write_log('Must set nx_block and ny_block > 0 for this compute_blocks option', GM_FATAL)
       endif
    endif

    !Note: The Glissade dycore is known to crash for thklim = 0, but has not
    !      been extensively tested for small values of thklim.
    !      Values smaller than 1 mm may be OK, but no guarantees.
    if (model%numerics%thklim < 1.d-3) then   ! 1 mm
       call write_log('ice limit (thklim) is too small for Glissade dycore', GM_FATAL)
    endif

    call log_section('numerics')
    call log_option('compute_blocks', model%options%compute_blocks, &
         compute_blocks(model%options%compute_blocks))
    call log_flag('reproducible_sums', model%options%reproducible_sums)
    call log_param('ice_limit', model%numerics%thklim, 'm; thickness limit for dynamically active ice')
    call log_param('ice_limit_temp', model%numerics%thklim_temp, 'm; thickness limit for temperature calculations')
    call log_param('thck_gradient_ramp', model%numerics%thck_gradient_ramp, &
         'm; 0 = no ramp of land-margin thickness gradients')
    if (model%general%global_bc == GLOBAL_BC_PERIODIC) then
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
    call GetValue(section,'adjust_input_thickness',  model%options%adjust_input_thickness)
    call GetValue(section,'smooth_input_topography', model%options%smooth_input_topography)
    call GetValue(section,'smooth_input_usrf',       model%options%smooth_input_usrf)
    call GetValue(section,'adjust_input_topography', model%options%adjust_input_topography)

    ! parameters to adjust input topography
    call GetValue(section,'adjust_topg_xmin',       model%paramets%adjust_topg_xmin)
    call GetValue(section,'adjust_topg_xmax',       model%paramets%adjust_topg_xmax)
    call GetValue(section,'adjust_topg_ymin',       model%paramets%adjust_topg_ymin)
    call GetValue(section,'adjust_topg_ymax',       model%paramets%adjust_topg_ymax)
    call GetValue(section,'adjust_topg_no_adjust',  model%paramets%adjust_topg_no_adjust)
    call GetValue(section,'adjust_topg_max_adjust', model%paramets%adjust_topg_max_adjust)
    call GetValue(section,'adjust_topg_delta',      model%paramets%adjust_topg_delta)

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
    call log_flag('adjust_input_thickness', model%options%adjust_input_thickness, &
         'if true, input thickness is adjusted based on surface and bed topography')
    call log_flag('smooth_input_usrf', model%options%smooth_input_usrf)
    call log_flag('smooth_input_topography', model%options%smooth_input_topography)
    call log_flag('adjust_input_topography', model%options%adjust_input_topography, &
         'if true, input topography is adjusted in a selected region')
    if (model%options%adjust_input_topography) then
       call log_param('adjust_topg_xmin', model%paramets%adjust_topg_xmin, 'm')
       call log_param('adjust_topg_xmax', model%paramets%adjust_topg_xmax, 'm')
       call log_param('adjust_topg_ymin', model%paramets%adjust_topg_ymin, 'm')
       call log_param('adjust_topg_ymax', model%paramets%adjust_topg_ymax, 'm')
       call log_param('adjust_topg_no_adjust', model%paramets%adjust_topg_no_adjust, 'm')
       call log_param('adjust_topg_max_adjust', model%paramets%adjust_topg_max_adjust, 'm')
       call log_param('adjust_topg_delta', model%paramets%adjust_topg_delta, 'm')
    endif

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
    call GetValue(section,'stokes_approx',                  model%options%which_ho_approx)
    call GetValue(section,'diva_slope_correction',          model%options%diva_slope_correction)
    call GetValue(section,'efvs_method',                    model%options%which_ho_efvs)
    call GetValue(section,'surface_gradient_method',        model%options%which_ho_gradient)
    call GetValue(section,'surface_gradient_margin_method', model%options%which_ho_gradient_margin)

    ! parameters
    call GetValue(section,'efvs_constant', model%paramets%efvs_constant)
    call GetValue(section,'effstrain_min', model%paramets%effstrain_min)
    call GetValue(section,'max_slope',     model%paramets%max_slope)

  end subroutine handle_velocity

!--------------------------------------------------------------------------------

  subroutine print_velocity(model)

    ! Velocity options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    character(len=*), dimension(0:2), parameter :: efvs_method = (/ &
         'constant value                 ', & 
         'multiple of flow factor        ', &
         'nonlinear, from eff strain rate' /)

    character(len=*), dimension(-1:5), parameter :: stokes_approx = (/ &
         'SIA only (glissade_velo_sia)                     ', &
         'SIA only (glissade_velo_higher)                  ', &
         'SSA only (glissade_velo_higher)                  ', &
         'Blatter-Pattyn HO (glissade_velo_higher)         ', &
         'Depth-integrated L1L2 (glissade_velo_higher)     ', &
         'Depth-integrated viscosity (glissade_velo_higher)', &
         'Hybrid SIA/SSA                                   ' /)

    character(len=*), dimension(0:2), parameter :: surface_gradient_method = (/ &
         'centered gradient (glissade)             ', &
         'first-order upstream gradient (glissade) ', &
         'second-order upstream gradient (glissade)' /)

    character(len=*), dimension(0:2), parameter :: surface_gradient_margin_method = (/ &
         'compute edge gradient when either cell has ice         ', &
         'compute edge gradient when ice lies above ice-free land', & 
         'compute edge gradient when both cells have ice         ' /)

    if (model%options%which_ho_approx < -1 .or. model%options%which_ho_approx >= size(stokes_approx)-1) then
       call write_log('Error, Stokes approximation out of range for glissade dycore', GM_FATAL)
    end if
    if (model%options%which_ho_efvs < 0 .or. model%options%which_ho_efvs >= size(efvs_method)) then
       call write_log('Error, HO effective viscosity input out of range', GM_FATAL)
    end if
    if (model%options%which_ho_gradient < 0 .or. model%options%which_ho_gradient >= size(surface_gradient_method)) then
       call write_log('Error, gradient option out of range for glissade dycore', GM_FATAL)
    end if
    if (model%options%which_ho_gradient_margin < 0 .or. &
        model%options%which_ho_gradient_margin >= size(surface_gradient_margin_method)) then
       call write_log('Error, gradient margin option out of range for glissade dycore', GM_FATAL)
    end if

    if (model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
       if (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then
          call write_log('Error, cannot use no-penetration BC with local SIA solver', GM_FATAL)
       endif
    endif  ! Glissade local SIA solver

    call log_section('velocity')
    call log_option('stokes_approx', model%options%which_ho_approx, &
         stokes_approx(model%options%which_ho_approx))
    if (model%options%which_ho_approx == HO_APPROX_DIVA) then
       call log_flag('diva_slope_correction', model%options%diva_slope_correction, &
            'if true, DIVA includes a correction for steep slopes')
    endif
    call log_option('efvs_method', model%options%which_ho_efvs, efvs_method(model%options%which_ho_efvs))
    if (model%options%which_ho_efvs == HO_EFVS_CONSTANT) then
       call log_param('efvs_constant', model%paramets%efvs_constant, 'Pa yr')
    elseif (model%options%which_ho_efvs == HO_EFVS_NONLINEAR) then
       call log_param('effstrain_min', model%paramets%effstrain_min, 'yr^-1')
    end if
    call log_option('surface_gradient_method', model%options%which_ho_gradient, &
         surface_gradient_method(model%options%which_ho_gradient))
    call log_option('surface_gradient_margin_method', model%options%which_ho_gradient_margin, &
         surface_gradient_margin_method(model%options%which_ho_gradient_margin))
    call log_param('max_slope', model%paramets%max_slope, 'maximum surface slope')

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
    call GetValue(section,'nonlinear_method',    model%options%which_ho_nonlinear)
    call GetValue(section,'nonlinear_residual',  model%options%which_ho_resid)
    call GetValue(section,'nonlinear_maxiters',  model%options%glissade_maxiter)

    ! linear solver
    call GetValue(section,'linear_solver',       model%options%which_ho_sparse)
    call GetValue(section,'preconditioner',      model%options%which_ho_precond)
    call GetValue(section,'linear_solve_ncheck', model%options%linear_solve_ncheck)
    call GetValue(section,'linear_maxiters',     model%options%linear_maxiters)
    call GetValue(section,'linear_tolerance',    model%options%linear_tolerance)

    ! matrix assembly
    call GetValue(section,'assemble_beta',       model%options%which_ho_assemble_beta)
    call GetValue(section,'assemble_taud',       model%options%which_ho_assemble_taud)
    call GetValue(section,'assemble_bfric',      model%options%which_ho_assemble_bfric)
    call GetValue(section,'assemble_lateral',    model%options%which_ho_assemble_lateral)

  end subroutine handle_velocity_solver

!--------------------------------------------------------------------------------

  subroutine print_velocity_solver(model)

    ! Velocity solver options

    use glissade_types
    use glimmer_log
    use glissade_log_config
    use cism_parallel, only: tasks

    implicit none

    type(glissade_global_type)  :: model

    character(len=100) :: note

    character(len=*), dimension(0:1), parameter :: nonlinear_method = (/ &
         'use standard Picard iteration          ', &
         'use Picard iteration with acceleration '/)

    character(len=*), dimension(0:4), parameter :: nonlinear_residual = (/ &
         'max value                   ', &
         'max value ignoring ubas     ', &
         'mean value                  ', &
         'L2 norm of Ax-b = resid     ', &
         'relative L2 norm, |Ax-b|/|b|' /)

    character(len=*), dimension(-1:3), parameter :: linear_solver = (/ &
         'PCG with incomplete Cholesky preconditioner', &
         'BiCG with LU preconditioner                ', &
         'GMRES with LU preconditioner               ', &
         'Native PCG solver, standard                ', &
         'Native PCG solver, Chronopoulos-Gear       '/)

    character(len=*), dimension(0:4), parameter :: preconditioner = (/ &
         'No preconditioner (native PCG)                ', &
         'Diagonal preconditioner (native PCG)          ', &
         'SIA preconditioner (native PCG)               ', &
         'Local tridiagonal preconditioner (native PCG) ', &
         'Global tridiagonal preconditioner (native PCG)' /)

    character(len=*), dimension(0:1), parameter :: assemble_beta = (/ &
         'standard finite-element assembly (glissade dycore) ', &
         'use local beta at each vertex (glissade dycore)    '  /)

    character(len=*), dimension(0:1), parameter :: assemble_taud = (/ &
         'standard finite-element assembly (glissade dycore)       ', &
         'use local driving stress at each vertex (glissade dycore)'  /)

    character(len=*), dimension(0:1), parameter :: assemble_bfric = (/ &
         'standard finite-element assembly (glissade dycore)       ', &
         'use local basal friction at each vertex (glissade dycore)'  /)

    character(len=*), dimension(0:1), parameter :: assemble_lateral = (/ &
         'standard finite-element assembly (glissade dycore)         ', &
         'use local thck and usrf on each cell face (glissade dycore)'  /)

    if (model%options%which_ho_nonlinear < 0 .or. model%options%which_ho_nonlinear >= size(nonlinear_method)) then
       call write_log('Error, HO nonlinear solution input out of range', GM_FATAL)
    end if
    if (model%options%which_ho_resid < 0 .or. model%options%which_ho_resid >= size(nonlinear_residual)) then
       call write_log('Error, HO residual input out of range', GM_FATAL)
    end if
    if (model%options%which_ho_sparse < lbound(linear_solver,1) .or. &
        model%options%which_ho_sparse > ubound(linear_solver,1)) then
       call write_log('Error, HO sparse solver input out of range', GM_FATAL)
    end if
    if (model%options%which_ho_precond < 0 .or. model%options%which_ho_precond >= size(preconditioner)) then
       call write_log('Error, glissade preconditioner out of range', GM_FATAL)
    end if
    if (model%options%which_ho_assemble_beta < 0 .or. &
        model%options%which_ho_assemble_beta >= size(assemble_beta)) then
       call write_log('Error, beta assembly option out of range for glissade dycore', GM_FATAL)
    end if
    if (model%options%which_ho_assemble_taud < 0 .or. &
        model%options%which_ho_assemble_taud >= size(assemble_taud)) then
       call write_log('Error, driving-stress assembly option out of range for glissade dycore', GM_FATAL)
    end if
    if (model%options%which_ho_assemble_bfric < 0 .or. &
        model%options%which_ho_assemble_bfric >= size(assemble_bfric)) then
       call write_log('Error, basal-friction assembly option out of range for glissade dycore', GM_FATAL)
    end if
    if (model%options%which_ho_assemble_lateral < 0 .or. &
        model%options%which_ho_assemble_lateral >= size(assemble_lateral)) then
       call write_log('Error, lateral-stress assembly option out of range for glissade dycore', GM_FATAL)
    end if

    ! Forbidden options for running in parallel
    ! Note: The local SIA solver does not use a sparse linear solver, so which_ho_sparse is ignored.
    if (tasks > 1 .and. model%options%which_ho_approx /= HO_APPROX_LOCAL_SIA .and.  &
                        (model%options%which_ho_sparse==HO_SPARSE_BICG  .or.  &
                         model%options%which_ho_sparse==HO_SPARSE_GMRES .or.  &
                         model%options%which_ho_sparse==HO_SPARSE_PCG_INCH) ) then
       call write_log('Error, SLAP solver not supported for more than one processor', GM_FATAL)
    end if

    ! The SIA preconditioner cannot be used for 2D solves
    note = ''
    if ( (model%options%which_ho_approx == HO_APPROX_SSA  .or.  &
          model%options%which_ho_approx == HO_APPROX_L1L2 .or.  &
          model%options%which_ho_approx == HO_APPROX_DIVA .or.  &
          model%options%which_ho_approx == HO_APPROX_HYBRID)    &
                             .and.                            &
          (model%options%which_ho_sparse == HO_SPARSE_PCG_STANDARD .or.    &
           model%options%which_ho_sparse == HO_SPARSE_PCG_CHRONGEAR) ) then
       if (model%options%which_ho_precond == HO_PRECOND_SIA) then
          !TODO - Change default to tridiagonal if it turns out to be faster than diagonal
          model%options%which_ho_precond = HO_PRECOND_DIAG
          note = '; reset from SIA, which cannot be used for 2D solves'
       endif
    endif

    if ((model%options%which_ho_sparse == HO_SPARSE_PCG_STANDARD .or.  &
         model%options%which_ho_sparse == HO_SPARSE_PCG_CHRONGEAR) ) then
       if (model%options%reproducible_sums) then
          if (model%options%which_ho_precond == HO_PRECOND_TRIDIAG_LOCAL .or. &
              model%options%which_ho_precond == HO_PRECOND_TRIDIAG_GLOBAL) then
             call write_log ('Tridiagonal preconditioners are not supported with reproducible sums.')
             call write_log ('Please choose a different preconditioner (e.g., diagonal)', GM_FATAL)
          endif
       endif
    endif

    call log_section('velocity_solver')
    call log_option('nonlinear_method', model%options%which_ho_nonlinear, &
         nonlinear_method(model%options%which_ho_nonlinear))
    call log_option('nonlinear_residual', model%options%which_ho_resid, &
         nonlinear_residual(model%options%which_ho_resid))
    call log_param('nonlinear_maxiters', model%options%glissade_maxiter)
    call log_option('linear_solver', model%options%which_ho_sparse, &
         linear_solver(model%options%which_ho_sparse))
    if ((model%options%which_ho_sparse == HO_SPARSE_PCG_STANDARD .or.  &
         model%options%which_ho_sparse == HO_SPARSE_PCG_CHRONGEAR) ) then
       call log_option('preconditioner', model%options%which_ho_precond, &
            trim(preconditioner(model%options%which_ho_precond)) // trim(note))
    end if
    call log_param('linear_solve_ncheck', model%options%linear_solve_ncheck)
    call log_param('linear_maxiters', model%options%linear_maxiters)
    call log_param('linear_tolerance', model%options%linear_tolerance)
    call log_option('assemble_beta', model%options%which_ho_assemble_beta, &
         assemble_beta(model%options%which_ho_assemble_beta))
    call log_option('assemble_taud', model%options%which_ho_assemble_taud, &
         assemble_taud(model%options%which_ho_assemble_taud))
    call log_option('assemble_bfric', model%options%which_ho_assemble_bfric, &
         assemble_bfric(model%options%which_ho_assemble_bfric))
    call log_option('assemble_lateral', model%options%which_ho_assemble_lateral, &
         assemble_lateral(model%options%which_ho_assemble_lateral))

    ! unsupported resid options
    if (model%options%which_ho_resid == HO_RESID_MAXU) then
      call write_log('Residual as max. value of normalized velocity vector update is not currently scientifically supported.  &
           &USE AT YOUR OWN RISK.', GM_WARNING)
    endif
    if (model%options%which_ho_resid == HO_RESID_MAXU_NO_UBAS) then
      call write_log('Residual as max. value of normalized velocity vector update with basal velocity omitted is not currently &
           &scientifically supported.  USE AT YOUR OWN RISK.', GM_WARNING)
    endif
    if (model%options%which_ho_resid == HO_RESID_MEANU) then
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

    call GetValue(section,'advection_method',     model%options%whichevol)
    call GetValue(section,'vertical_remap_order', model%options%which_ho_vertical_remap)
    call GetValue(section,'compute_ice_age',      model%options%which_ho_ice_age)

  end subroutine handle_advection

!--------------------------------------------------------------------------------

  subroutine print_advection(model)

    ! Advection options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    !WHL TODO - Renumber during the config reorg
    character(len=*), dimension(3:5), parameter :: advection_method = (/ &
         'incremental remapping                 ', &   
         '1st order upwind                      ', &
         'thickness fixed at initial value      '  /)

    character(len=*), dimension(0:1), parameter :: vertical_remap_order = (/ &
         'first-order accurate  ', &
         'second-order accurate ' /)

    character(len=*), dimension(0:1), parameter :: compute_ice_age = (/ &
         'ice age computation off', &
         'ice age computation on ' /)

    !WHL TODO: Renumber during the config reorg and restore the following commented-out line
!!    if (model%options%whichevol < 0 .or. model%options%whichevol >= size(advection_method)) then
    if (model%options%whichevol < EVOL_INC_REMAP .or. model%options%whichevol > EVOL_NO_THICKNESS) then
       call write_log('Error, advection_method out of range',GM_FATAL)
    end if
    if (model%options%which_ho_vertical_remap < 0 .or. &
        model%options%which_ho_vertical_remap >= size(vertical_remap_order)) then
       call write_log('Error, vertical remap option out of range for glissade dycore', GM_FATAL)
    end if
    if (model%options%which_ho_ice_age < 0 .or. model%options%which_ho_ice_age >= size(compute_ice_age)) then
       call write_log('Error, compute_ice_age option out of range for glissade dycore', GM_FATAL)
    end if

    call log_section('advection')
    call log_option('advection_method', model%options%whichevol, advection_method(model%options%whichevol))
    call log_option('vertical_remap_order', model%options%which_ho_vertical_remap, &
         vertical_remap_order(model%options%which_ho_vertical_remap))
    call log_option('compute_ice_age', model%options%which_ho_ice_age, &
         compute_ice_age(model%options%which_ho_ice_age))

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
    call GetValue(section,'temperature_method', model%options%whichtemp)
    call GetValue(section,'temp_init',          model%options%temp_init)
    call GetValue(section,'dissipation_method', model%options%which_ho_disp)
    call GetValue(section,'thermal_timestep',   model%options%which_ho_thermal_timestep)

    ! parameters
    call GetValue(section,'pmp_offset',    model%temper%pmp_offset)
    call GetValue(section,'pmp_threshold', model%temper%pmp_threshold)

  end subroutine handle_temperature

!--------------------------------------------------------------------------------

  subroutine print_temperature(model)

    ! Temperature options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    character(len=100) :: note

    character(len=*), dimension(0:3), parameter :: temperature_method = (/ &
         'isothermal            ', &
         'prognostic temperature', &
         'constant in time      ', &
         'prognostic enthalpy   ' /)

    character(len=*), dimension(0:4), parameter :: temp_init = (/ &
         'set to 0 C                  ', &
         'set to surface air temp     ', &
         'linear vertical profile     ', &
         'advective-diffusive balance ',&
         'temp from external file     ' /)

    character(len=*), dimension(-1:1), parameter :: dissipation_method = (/ &
         'no dissipation                    ', &
         '0-order SIA                       ', &
         'first-order model (Blatter-Pattyn)' /)

    character(len=*), dimension(0:2), parameter :: thermal_timestep = (/ &
         'vertical thermal solve before transport    ', &
         'vertical thermal solve after transport     ', &
         'vertical thermal solve split into two parts' /)

    !NOTE : Old option 3 (TEMP_REMAP_ADV) has been removed.
    ! If this has been set, then change to option 1 (TEMP_PROGNOSTIC), which applies to any dycore.

    if (model%options%whichtemp < 0 .or. model%options%whichtemp >= size(temperature_method)) then
       call write_log('Error, temperature_method out of range',GM_FATAL)
    end if
    if (model%options%temp_init < 0 .or. model%options%temp_init >= size(temp_init)) then
       call write_log('Error, temp_init option out of range',GM_FATAL)
    end if
    if (model%options%which_ho_disp < -1 .or. model%options%which_ho_disp >= size(dissipation_method)-1) then
       call write_log('Error, HO dissipation input out of range', GM_FATAL)
    end if
    if (model%options%which_ho_thermal_timestep < 0 .or. &
         model%options%which_ho_thermal_timestep >= size(thermal_timestep)) then
       call write_log('Error, HO thermal timestep input out of range', GM_FATAL)
    end if

    ! The local SIA solver cannot use first-order dissipation
    note = ''
    if (model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
       if (model%options%which_ho_disp == HO_DISP_FIRSTORDER ) then
          model%options%which_ho_disp = HO_DISP_SIA
          note = '; reset from first-order, which cannot be used with the local SIA solver'
       endif
    endif  ! Glissade local SIA solver

    call log_section('temperature')
    call log_option('temperature_method', model%options%whichtemp, temperature_method(model%options%whichtemp))
    ! Note: If temperature is read from an input or restart file, the temp_init option is overridden.
    !       The method actually used is written to the log file by glissade_init_temp.
    call log_option('temp_init', model%options%temp_init, &
         trim(temp_init(model%options%temp_init)) // '; not used if temperature is read from a file')
    call log_option('dissipation_method', model%options%which_ho_disp, &
         trim(dissipation_method(model%options%which_ho_disp)) // trim(note))
    call log_option('thermal_timestep', model%options%which_ho_thermal_timestep, &
         thermal_timestep(model%options%which_ho_thermal_timestep))
    call log_param('pmp_threshold', model%temper%pmp_threshold, 'K')
    call log_param('pmp_offset', model%temper%pmp_offset, 'K')

    ! unsupported temperature options
    if (model%options%whichtemp == TEMP_ENTHALPY) then
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
    call GetValue(section,'flow_law',             model%options%whichflwa)
    call GetValue(section,'damage_flwa_feedback', model%options%damage_flwa_feedback)

    ! parameters
    call GetValue(section,'n_glen',             n_glen)
    call GetValue(section,'flow_factor_ground', model%paramets%flow_enhancement_factor_ground)
    call GetValue(section,'flow_factor_float',  model%paramets%flow_enhancement_factor_float)
    !TODO - Change default_flwa to flwa_constant?  Would have to change config files.
    call GetValue(section,'default_flwa',       model%paramets%default_flwa)

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

    character(len=*), dimension(0:2), parameter :: flow_law = (/ &
         'uniform factor flwa         ', &
         'Paterson and Budd (T = -5 C)', &
         'Paterson and Budd           ' /)

    if (model%options%whichflwa < 0 .or. model%options%whichflwa >= size(flow_law)) then
       call write_log('Error, flow_law out of range',GM_FATAL)
    end if

    call log_section('flow_law')
    call log_option('flow_law', model%options%whichflwa, flow_law(model%options%whichflwa))
    if (model%options%whichflwa == FLWA_CONST_FLWA) then
       call log_param('default_flwa', model%paramets%default_flwa, 'Pa^-n yr^-1; constant flow factor')
    end if
    call log_param('n_glen', n_glen, 'exponent in Glen''s flow law')
    call log_param('flow_factor_ground', model%paramets%flow_enhancement_factor_ground, 'grounded ice')
    call log_param('flow_factor_float', model%paramets%flow_enhancement_factor_float, 'floating ice')
    call log_flag('damage_flwa_feedback', model%options%damage_flwa_feedback)

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
    call GetValue(section,'sliding_law',               model%options%which_ho_babc)
    call GetValue(section,'use_c_space_factor',        model%options%use_c_space_factor)
    call GetValue(section,'beta_limit_method',         model%options%which_ho_beta_limit)
    call GetValue(section,'powerlaw_c_method',         model%options%which_ho_powerlaw_c)
    call GetValue(section,'coulomb_c_method',          model%options%which_ho_coulomb_c)
    call GetValue(section,'elevation_based_coulomb_c', model%options%elevation_based_coulomb_c)
    call GetValue(section,'effecpress_method',         model%options%which_ho_effecpress)

    ! basal traction for the local SIA solver (Glide-style; inverse sense of beta)
    call GetValue(section,'btrc_method', model%options%whichbtrc)
    call GetValue(section,'btrc_const',  model%paramets%btrc_const)

    ! beta parameters
    call GetValue(section,'beta_grounded_min', model%basal_physics%beta_grounded_min)
    call GetValue(section,'beta_const',        model%basal_physics%ho_beta_const)
    call GetValue(section,'beta_small',        model%basal_physics%ho_beta_small)
    call GetValue(section,'beta_large',        model%basal_physics%ho_beta_large)

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

    ! Note: Glide also supports 'linear function of bmlt' (4) and 'tanh function of bwat' (5)
    character(len=*), dimension(0:3), parameter :: btrc_method = (/ &
         'no basal sliding       ', &
         'constant basal traction', &
         'constant where bwat > 0', &
         'constant where T = Tpmp' /)

    !WHL TODO - Remove option 3 (and maybe others) during the config reorg
    character(len=*), dimension(0:14), parameter :: sliding_law = (/ &
         'constant beta                                    ', &
         'beta depends on basal temp (melting or frozen)   ', &
         'pseudo-plastic sliding law, new C_c options      ', &
         'pseudo-plastic sliding law, old tan(phi) options ', &
         'no slip (using large B^2)                        ', &
         'beta from external file                          ', &
         'no slip (Dirichlet implementation)               ', &
         'Zoet-Iverson sliding law                         ', &
         'beta as in ISMIP-HOM test C                      ', &
         'power law                                        ', &
         'Coulomb friction law w/ effec press              ', &
         'Schoof friction law                              ', &
         'min of Coulomb stress and power-law stress (Tsai)', &
         'power law using effective pressure               ', &
         'till yield stress (Picard)                       ' /)

    character(len=*), dimension(0:1), parameter :: beta_limit_method = (/ &
         'absolute beta limit based on beta_grounded_min   ', &
         'beta is limited, then scaled by f_ground_cell    ' /)

    character(len=*), dimension(0:4), parameter :: powerlaw_c_method = (/ &
         'spatially uniform friction parameter Cp     ', &
         'invert for 2D friction parameter Cp         ', &
         'read friction parameter Cp from file        ', &
         'invert for basin-scale friction parameter Cp', &
         'friction parameter Cp is a function of Cc   ' /)

    character(len=*), dimension(0:3), parameter :: coulomb_c_method = (/ &
         'spatially uniform friction parameter Cc     ', &
         'invert for 2D friction parameter Cc         ', &
         'read friction parameter Cc from file        ', &
         'invert for basin-scale friction parameter Cc' /)

    character(len=*), dimension(0:4), parameter :: effecpress_method = (/ &
         'full overburden pressure                             ', &
         'reduced effecpress near pressure melting point       ', &
         'macroporous sheet; effecpress depends on bwat        ', &
         'cavity-sheet; effecpress depends on opening, closing ', &
         'reduced effecpress where bwat > 0 (B/vP)             ' /)

    if (model%options%whichbtrc < 0 .or. model%options%whichbtrc >= size(btrc_method)) then
       call write_log('Error, btrc_method out of range',GM_FATAL)
    end if
    if (model%options%which_ho_babc < 0 .or. model%options%which_ho_babc >= size(sliding_law)) then
       call write_log('Error, HO basal BC input out of range', GM_FATAL)
    end if
    if (model%options%which_ho_beta_limit < 0 .or. model%options%which_ho_beta_limit >= size(beta_limit_method)) then
       call write_log('Error, HO beta limit input out of range', GM_FATAL)
    end if
    if (model%options%which_ho_powerlaw_c < 0 .or. model%options%which_ho_powerlaw_c >= size(powerlaw_c_method)) then
       call write_log('Error, HO powerlaw_c input out of range', GM_FATAL)
    end if
    if (model%options%which_ho_coulomb_c < 0 .or. model%options%which_ho_coulomb_c >= size(coulomb_c_method)) then
       call write_log('Error, coulomb_c_method out of range', GM_FATAL)
    end if
    if (model%options%which_ho_effecpress < 0 .or. model%options%which_ho_effecpress >= size(effecpress_method)) then
       call write_log('Error, HO effective pressure input out of range', GM_FATAL)
    end if

    if (model%options%use_c_space_factor) then
       if (model%options%which_ho_babc /= HO_BABC_COULOMB_FRICTION .and.  &
           model%options%which_ho_babc /= HO_BABC_SCHOOF .and. &
           model%options%which_ho_babc /= HO_BABC_TSAI) then
          call write_log('Error, C_space_factor not supported for this choice of sliding_law', GM_FATAL)
       endif
    endif

    ! Inversion options

    ! coulomb_c inversion not supported with elevation-based coulomb_c
    if (model%options%elevation_based_coulomb_c) then
       if (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION .or. &
           model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION_BASIN) then
          call write_log('Error, inversion not supported with elevation_based_coulomb_c = T', GM_FATAL)
       endif
    endif

    ! Note: Cp inversion supported for the basic power law and the Schoof and Tsai laws
    if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or. &
        model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION_BASIN) then
       if (model%options%which_ho_babc == HO_BABC_POWERLAW .or. &
           model%options%which_ho_babc == HO_BABC_SCHOOF .or.  &
           model%options%which_ho_babc == HO_BABC_TSAI) then
          ! inversion for Cp is supported
       else
          call write_log('Error, Cp inversion is not supported for this basal BC option')
          write(message,*) 'Cp inversion is supported for these options: ', &
               HO_BABC_POWERLAW, HO_BABC_SCHOOF, HO_BABC_TSAI
          call write_log(message, GM_FATAL)
       endif
    endif

    ! Note: Cc inversion supported for the Zoet-Iverson law, pseudoplastic, Schoof and Tsai laws
    if (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION .or. &
        model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION_BASIN) then
       if (model%options%which_ho_babc == HO_BABC_ZOET_IVERSON .or. &
           model%options%which_ho_babc == HO_BABC_PSEUDO_PLASTIC .or. &
           model%options%which_ho_babc == HO_BABC_SCHOOF .or.  &
           model%options%which_ho_babc == HO_BABC_TSAI) then
          ! inversion for Cc is supported
       else
          call write_log('Error, Cc inversion is not supported for this basal BC option')
          write(message,*) 'Cc inversion is supported for these options: ', &
               HO_BABC_ZOET_IVERSON, HO_BABC_PSEUDO_PLASTIC, HO_BABC_SCHOOF, HO_BABC_TSAI
          call write_log(message, GM_FATAL)
       endif
    endif

    ! Note: The Schoof and Tsai laws permit inversion for Cc combined with
    !        which_ho_powerlaw_c = HO_POWERLAW_C_FUNCTION_COULOMB_C.
    !       In these cases we invert for Cc, with Cp coming along for the ride.
    !       There is currently no option to invert for Cc and Cp independently.
    if ( (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or. &
          model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION_BASIN) &
          .and. &
         (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION .or. &
          model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION_BASIN) ) then
       call write_log('Simultaneous independent inversion for Cp and Cc is not supported', &
            GM_FATAL)
    endif

    ! Make sure the user isn't trying to invert for powerlaw_c and coulomb_c independently.
    ! We can invert for both at once only if the two coefficients have a specified functional relationship.
    if ( (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .and. &
          model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION) &
                                        .or. &
         (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION_BASIN .and. &
          model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION_BASIN) ) then
       call write_log('Cannot invert for powerlaw_c and coulomb_c independently')
       write(message,*) 'Please set powerlaw_c_method =', HO_POWERLAW_C_FUNCTION_COULOMB_C
       call write_log(message, GM_FATAL)
    endif

    call log_section('basal_physics')
    call log_option('sliding_law', model%options%which_ho_babc, sliding_law(model%options%which_ho_babc))
    if (model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
       call log_option('btrc_method', model%options%whichbtrc, btrc_method(model%options%whichbtrc))
    endif
    call log_option('beta_limit_method', model%options%which_ho_beta_limit, &
         beta_limit_method(model%options%which_ho_beta_limit))
    call log_option('powerlaw_c_method', model%options%which_ho_powerlaw_c, &
         powerlaw_c_method(model%options%which_ho_powerlaw_c))
    call log_option('coulomb_c_method', model%options%which_ho_coulomb_c, &
         coulomb_c_method(model%options%which_ho_coulomb_c))
    call log_flag('use_c_space_factor', model%options%use_c_space_factor, 'if true, beta is multiplied by C_space_factor')
    call log_flag('elevation_based_coulomb_c', model%options%elevation_based_coulomb_c)
    call log_option('effecpress_method', model%options%which_ho_effecpress, &
         effecpress_method(model%options%which_ho_effecpress))

    if (model%options%whichbtrc == BTRC_CONSTANT      .or.  &
        model%options%whichbtrc == BTRC_CONSTANT_BWAT .or.  &
        model%options%whichbtrc == BTRC_CONSTANT_BPMP) then
       call log_param('btrc_const', model%paramets%btrc_const, 'm/yr/Pa; basal traction parameter for the local SIA solver')
    end if

    if (model%options%which_ho_babc == HO_BABC_BETA_CONSTANT) then
       call log_param('beta_const', model%basal_physics%ho_beta_const, 'Pa yr/m')
    elseif (model%options%which_ho_babc == HO_BABC_BETA_LARGE) then
       call log_param('beta_large', model%basal_physics%ho_beta_large, 'Pa yr/m')
    elseif (model%options%which_ho_babc == HO_BABC_BETA_BPMP) then
       call log_param('beta_large', model%basal_physics%ho_beta_large, 'Pa yr/m')
       call log_param('beta_small', model%basal_physics%ho_beta_small, 'Pa yr/m')
    elseif (model%options%which_ho_babc == HO_BABC_PSEUDO_PLASTIC_OLD .or.  &
            model%options%which_ho_babc == HO_BABC_PSEUDO_PLASTIC) then
       call log_param('pseudo_plastic_q', model%basal_physics%pseudo_plastic_q)
       call log_param('pseudo_plastic_u0', model%basal_physics%pseudo_plastic_u0)
       if (model%options%which_ho_babc == HO_BABC_PSEUDO_PLASTIC_OLD) then
          call log_param('pseudo_plastic_phimin', model%basal_physics%pseudo_plastic_phimin, 'deg')
          call log_param('pseudo_plastic_phimax', model%basal_physics%pseudo_plastic_phimax, 'deg')
          call log_param('pseudo_plastic_bedmin', model%basal_physics%pseudo_plastic_bedmin, 'm')
          call log_param('pseudo_plastic_bedmax', model%basal_physics%pseudo_plastic_bedmax, 'm')
       endif
       if (model%options%which_ho_assemble_beta == HO_ASSEMBLE_BETA_STANDARD) then
          call write_log('WARNING: local beta assembly is recommended for the pseudo-plastic sliding law')
          write(message,*) 'Set assemble_beta =', HO_ASSEMBLE_BETA_LOCAL
          call write_log(message, GM_WARNING)
       endif
       if (model%options%which_ho_assemble_bfric == HO_ASSEMBLE_BFRIC_STANDARD) then
          call write_log('WARNING: local bfric assembly is recommended for the pseudo-plastic sliding law')
          write(message,*) 'Set assemble_bfric =', HO_ASSEMBLE_BFRIC_LOCAL
          call write_log(message, GM_WARNING)
       endif
       if (model%options%which_ho_thermal_timestep == HO_THERMAL_BEFORE_TRANSPORT) then
          call write_log('WARNING: Best to do vertical thermal solve after transport for the pseudo-plastic sliding law')
          write(message,*) 'Set thermal_timestep to one of the following values:', &
               HO_THERMAL_AFTER_TRANSPORT, HO_THERMAL_SPLIT_TIMESTEP
          call write_log(message, GM_WARNING)
       endif
    elseif (model%options%which_ho_babc == HO_BABC_ZOET_IVERSON) then
       ! Note: The Zoet-Iverson law typically uses a spatially variable coulomb_c.
       !       If so, the value written here is just the initial value.
       call log_param('coulomb_c_const', model%basal_physics%coulomb_c_const)
       call log_param('coulomb_c_max', model%basal_physics%coulomb_c_max)
       call log_param('coulomb_c_min', model%basal_physics%coulomb_c_min)
       call log_param('powerlaw_m', model%basal_physics%powerlaw_m)
       call log_param('zoet_iverson_ut', model%basal_physics%zoet_iverson_ut, 'm/yr')
       call log_param('zoet_iverson_nmax', model%basal_physics%zoet_iverson_nmax, 'Pa')
    elseif (model%options%which_ho_babc == HO_BABC_ISHOMC) then
       if (model%general%ewn /= model%general%nsn) then
          call write_log('Error, must have ewn = nsn for ISMIP-HOM test C', GM_FATAL)
       endif
    elseif (model%options%which_ho_babc == HO_BABC_POWERLAW) then
       call log_param('powerlaw_c_const', model%basal_physics%powerlaw_c_const, 'Pa (m/yr)^(-1/m)')
       call log_param('powerlaw_c_max', model%basal_physics%powerlaw_c_max)
       call log_param('powerlaw_c_min', model%basal_physics%powerlaw_c_min)
       call log_param('powerlaw_m', model%basal_physics%powerlaw_m)
    elseif (model%options%which_ho_babc == HO_BABC_COULOMB_FRICTION) then
       call log_param('coulomb_c_const', model%basal_physics%coulomb_c_const)
       call log_param('coulomb_bump_max_slope', model%basal_physics%coulomb_bump_max_slope)
       call log_param('coulomb_bump_wavelength', model%basal_physics%coulomb_bump_wavelength, 'm')
    elseif (model%options%which_ho_babc == HO_BABC_SCHOOF) then
       ! Note: The Schoof and Tsai laws often invert for a spatially variable powerlaw_c or coulomb_c.
       !       If so, the values written here are just the initial values.
       call log_param('coulomb_c_const', model%basal_physics%coulomb_c_const)
       if (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION .or. &
           model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION_BASIN) then
          call log_param('coulomb_c_max', model%basal_physics%coulomb_c_max)
          call log_param('coulomb_c_min', model%basal_physics%coulomb_c_min)
       endif
       if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_FUNCTION_COULOMB_C) then
          call log_comment('Cp = schoof_gamma * Cc^schoof_p')
          call log_param('schoof_gamma', model%basal_physics%schoof_gamma)
          call log_param('schoof_p', model%basal_physics%schoof_p)
       else
          call log_param('powerlaw_c_const', model%basal_physics%powerlaw_c_const, 'Pa (m/yr)^(-1/m)')
       endif
       if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or. &
           model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION_BASIN) then
          call log_param('powerlaw_c_max', model%basal_physics%powerlaw_c_max)
          call log_param('powerlaw_c_min', model%basal_physics%powerlaw_c_min)
       endif
       call log_param('powerlaw_m', model%basal_physics%powerlaw_m)
       call log_param('schoof_n', model%basal_physics%schoof_n)
    elseif (model%options%which_ho_babc == HO_BABC_TSAI) then
       call log_param('coulomb_c_const', model%basal_physics%coulomb_c_const)
       if (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION .or. &
           model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION_BASIN) then
          call log_param('coulomb_c_max', model%basal_physics%coulomb_c_max)
          call log_param('coulomb_c_min', model%basal_physics%coulomb_c_min)
       endif
       if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_FUNCTION_COULOMB_C) then
          call log_comment('Cp = schoof_gamma * Cc^schoof_p')
          call log_param('schoof_gamma', model%basal_physics%schoof_gamma)
          call log_param('schoof_p', model%basal_physics%schoof_p)
       else
          call log_param('powerlaw_c_const', model%basal_physics%powerlaw_c_const, 'Pa (m/yr)^(-1/m)')
       endif
       if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or. &
           model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION_BASIN) then
          call log_param('powerlaw_c_max', model%basal_physics%powerlaw_c_max)
          call log_param('powerlaw_c_min', model%basal_physics%powerlaw_c_min)
       endif
       call log_param('powerlaw_m', model%basal_physics%powerlaw_m)
    elseif (model%options%which_ho_babc == HO_BABC_POWERLAW_EFFECPRESS) then
       call write_log('Weertman-style power law higher-order basal boundary condition is not currently scientifically &
            &supported.  USE AT YOUR OWN RISK.', GM_WARNING)
       !TODO - Use powerlaw_c instead of friction_powerlaw_k?  Allow p and q to be set in config file instead of hard-wired?
       call log_param('friction_powerlaw_k', model%basal_physics%friction_powerlaw_k, 'roughness parameter for the power-law friction law')
    endif

    ! Coulomb elevation parameters
    if (model%options%elevation_based_coulomb_c) then
       call log_param('coulomb_c_hi', model%basal_physics%coulomb_c_hi)
       call log_param('coulomb_c_lo', model%basal_physics%coulomb_c_lo)
       call log_param('coulomb_c_bed_hi', model%basal_physics%coulomb_c_bed_hi, 'm')
       call log_param('coulomb_c_bed_lo', model%basal_physics%coulomb_c_bed_lo, 'm')
    endif

    call log_param('beta_powerlaw_umax', model%basal_physics%beta_powerlaw_umax, &
         'm/yr; max ice speed when evaluating beta(u); 0 = no limit')

    if (model%options%which_ho_beta_limit == HO_BETA_LIMIT_ABSOLUTE) then
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
    call GetValue(section, 'bwat_method', model%options%which_ho_bwat)

    ! constant basal water
    call GetValue(section, 'const_bwat', model%basal_hydro%const_bwat)

    ! local basal till
    call GetValue(section, 'bwat_till_max', model%basal_hydro%bwat_till_max)
    call GetValue(section, 'c_drainage', model%basal_hydro%c_drainage)

    ! flux routing
    call GetValue(section, 'flux_routing_method', model%basal_hydro%ho_flux_routing_scheme)
    call GetValue(section, 'const_source', model%basal_hydro%const_source)
    call GetValue(section, 'btemp_flow_scale', model%basal_hydro%btemp_flow_scale)
    call GetValue(section, 'btemp_freeze_scale', model%basal_hydro%btemp_freeze_scale)

    ! effective pressure options and parameters
    call GetValue(section, 'effecpress_delta',   model%basal_hydro%effecpress_delta)
    call GetValue(section, 'bwat_threshold', model%basal_hydro%bwat_threshold)
    call GetValue(section, 'bpmp_threshold', model%basal_hydro%bpmp_threshold)
    call GetValue(section, 'haf_threshold', model%basal_hydro%haf_threshold)
    call GetValue(section, 'cavity_open_slide', model%basal_hydro%cavity_open_slide)
    call GetValue(section, 'cavity_open_melt', model%basal_hydro%cavity_open_melt)
    call GetValue(section, 'bump_height', model%basal_hydro%bump_height)
    call GetValue(section, 'bump_wavelength', model%basal_hydro%bump_wavelength)
    call GetValue(section, 'sliding_speed_fixed', model%basal_hydro%sliding_speed_fixed)
    call GetValue(section, 'flwa_basal', model%basal_hydro%flwa_basal)
    call GetValue(section, 'c_close', model%basal_hydro%c_close)

  end subroutine handle_basal_hydro

!--------------------------------------------------------------------------------

  subroutine print_basal_hydro(model)

    ! Basal hydrology options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none
    type(glissade_global_type)  :: model
    character(len=100) :: message

    character(len=*), dimension(0:3), parameter :: bwat_method = (/ &
         'zero basal water depth                           ', &
         'constant basal water depth                       ', &
         'basal water depth computed from local till model ', &
         'steady-state flux routing                        ' /)

    character(len=*), dimension(0:2), parameter :: flux_routing_method = (/ &
         'D8; route flux to lowest-elevation neighbor      ', &
         'Dinf; route flux to two lower-elevation neighbors', &
         'FD8; route flux to all lower-elevation neighbors ' /)

    character(len=*), dimension(0:2), parameter :: cavity_open_slide = (/ &
         'no opening by sliding                           ', &
         'opening by sliding with fixed u_b               ', &
         'opening by sliding with dynamic u_b             ' /)

    character(len=*), dimension(0:3), parameter :: cavity_open_melt = (/ &
         'no opening by melting                           ', &
         'opening by melting based on bmlt_ground         ', &
         'opening by melting based on cavity dissipation  ', &
         'opening based on bmlt_ground plus dissipation   ' /)

    ! check the options

    if (model%options%which_ho_bwat < 0 .or. model%options%which_ho_bwat >= size(bwat_method)) then
       call write_log('Error, HO basal water input out of range', GM_FATAL)
    end if

    if (model%options%which_ho_bwat == HO_BWAT_LOCAL_TILL) then
       if (model%options%which_ho_effecpress /= HO_EFFECPRESS_BWAT_BVP) then
          write(message,*) 'local till option must use effecpress_method =', HO_EFFECPRESS_BWAT_BVP
          call write_log(message, GM_FATAL)
       endif
    elseif (model%options%which_ho_bwat == HO_BWAT_FLUX_ROUTING) then
       if (model%basal_hydro%ho_flux_routing_scheme < 0.or. &
            model%basal_hydro%ho_flux_routing_scheme >= size(flux_routing_method)) then
          call write_log('Error, HO flux routing scheme out of range', GM_FATAL)
       end if
       if (model%options%reproducible_sums) then
          if (model%basal_hydro%ho_flux_routing_scheme /= HO_FLUX_ROUTING_D8) then
             call write_log('Error, this flux_routing scheme is not supported with reproducible sums')
             call write_log('Please use D8 flux-routing or turn off reproducible sums', GM_FATAL)
          endif
          if (model%basal_hydro%btemp_freeze_scale > 0.0d0) then
             call write_log('Error, flux-routing does not support refreezing with reproducible sums')
             call write_log('Please set btemp_freeze_scale = 0 or turn off reproducible sums', GM_FATAL)
          endif
       endif  ! reproducible sums
       if (model%options%which_ho_effecpress == HO_EFFECPRESS_CAVITY_SHEET) then
          if (model%basal_hydro%cavity_open_slide < 0 .or. &
               model%basal_hydro%cavity_open_slide >= size(cavity_open_slide)) then
             call write_log('Error, cavity_open_slide out of range', GM_FATAL)
          endif
          if (model%basal_hydro%cavity_open_melt < 0 .or. &
               model%basal_hydro%cavity_open_melt >= size(cavity_open_melt)) then
             call write_log('Error, cavity_open_melt out of range', GM_FATAL)
          endif
       endif
    endif

    ! write the options and parameters

    call log_section('basal_hydro')
    call log_option('bwat_method', model%options%which_ho_bwat, bwat_method(model%options%which_ho_bwat))

    if (model%options%which_ho_bwat == HO_BWAT_CONSTANT) then
       call log_param('const_bwat', model%basal_hydro%const_bwat, 'm; constant basal water depth')
    elseif (model%options%which_ho_bwat == HO_BWAT_LOCAL_TILL) then
       call log_param('bwat_till_max', model%basal_hydro%bwat_till_max, 'm; maximum till water depth')
       call log_param('c_drainage', model%basal_hydro%c_drainage, 'm/yr; till drainage rate')
       call log_param('effecpress_delta', model%basal_hydro%effecpress_delta, &
            'min effective pressure relative to overburden')
    elseif (model%options%which_ho_bwat == HO_BWAT_FLUX_ROUTING) then
       call log_option('flux_routing_method', model%basal_hydro%ho_flux_routing_scheme, &
            flux_routing_method(model%basal_hydro%ho_flux_routing_scheme))
       call log_param('const_source', model%basal_hydro%const_source, 'm/yr; constant melt source at the bed')
       call log_param('btemp_flow_scale', model%basal_hydro%btemp_flow_scale, &
            'deg C; temperature scale for flow around frozen bed; 0 = off')
       call log_param('btemp_freeze_scale', model%basal_hydro%btemp_freeze_scale, &
            'deg C; temperature scale for refreezing at the bed; 0 = off')
       if (model%options%which_ho_effecpress == HO_EFFECPRESS_BWAT) then
          call log_param('effecpress_delta', model%basal_hydro%effecpress_delta, &
               'min effective pressure relative to overburden')
          call log_param('bwat_threshold', model%basal_hydro%bwat_threshold, 'm')
       elseif (model%options%which_ho_effecpress == HO_EFFECPRESS_CAVITY_SHEET) then
          call log_option('cavity_open_slide', model%basal_hydro%cavity_open_slide, &
               cavity_open_slide(model%basal_hydro%cavity_open_slide))
          call log_option('cavity_open_melt', model%basal_hydro%cavity_open_melt, &
               cavity_open_melt(model%basal_hydro%cavity_open_melt))
          if (model%basal_hydro%cavity_open_slide /= CAVITY_OPEN_SLIDE_NONE) then
             call log_param('bump_height', model%basal_hydro%bump_height, 'm')
             call log_param('bump_wavelength', model%basal_hydro%bump_wavelength, 'm')
             if (model%basal_hydro%cavity_open_slide == CAVITY_OPEN_SLIDE_FIXED_UB) then
                call log_param('sliding_speed_fixed', model%basal_hydro%sliding_speed_fixed, 'm/yr')
             endif
          endif   ! opening by sliding
          call log_param('flwa_basal', model%basal_hydro%flwa_basal, 'Pa^-3 yr^-1; flow factor for cavity closing')
          call log_param('c_close', model%basal_hydro%c_close, 'constant for cavity closing')
       endif   ! cavity-sheet parameterization
    endif   ! which_ho_bwat

    ! other effective pressure options

    if (model%options%which_ho_effecpress == HO_EFFECPRESS_BPMP) then
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

    call GetValue(section,'grounding_line_param', model%options%which_ho_ground)
    call GetValue(section,'fground_no_glp',       model%options%which_ho_fground_no_glp)
    call GetValue(section,'grounding_line_bmlt',  model%options%which_ho_ground_bmlt)
    call GetValue(section,'flotation_function',   model%options%which_ho_flotation_function)

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

    character(len=*), dimension(0:2), parameter :: grounding_line_param = (/ &
         'f_ground = 0 or 1; no GLP  (glissade dycore)               ', &
         'GLP for basal friction with 0 <= f_ground <= 1 for vertices', &
         'deluxe GLP, 0 <= f_ground <= 1 for both vertices and cells ' /)

    character(len=*), dimension(0:1), parameter :: fground_no_glp = (/ &
         'f_ground = 1 at vertex if any neighbor cell is grounded   ', &
         'f_ground = 0 or 1 at vertex based on staggered f_flotation' /)

    character(len=*), dimension(0:2), parameter :: grounding_line_bmlt = (/ &
         'no GLP for bmlt_float                        ', &
         'weigh bmlt_float by floating fraction of cell', &
         'set bmlt_float = 0 in partly grounded cells  ' /)

    character(len=*), dimension(0:2), parameter :: flotation_function = (/ &
         'f_pattyn = (-rhoo*b)/(rhoi*H)     ', &
         '1/fpattyn = (rhoi*H)/(-rhoo*b)    ', &
         'linear = -b - (rhoi/rhoo)*H       ' /)

    ! check the options

    if (model%options%which_ho_ground < 0 .or. model%options%which_ho_ground >= size(grounding_line_param)) then
       call write_log('Error, grounding_line_param out of range for glissade dycore', GM_FATAL)
    end if

    if (model%options%which_ho_ground == HO_GROUND_NO_GLP) then
       if (model%options%which_ho_fground_no_glp < 0 .or. &
           model%options%which_ho_fground_no_glp >= size(fground_no_glp)) then
          call write_log('Error, fground_no_glp out of range for glissade dycore', GM_FATAL)
       end if
    endif

    if (model%options%which_ho_ground_bmlt < 0 .or. &
        model%options%which_ho_ground_bmlt >= size(grounding_line_bmlt)) then
       call write_log('Error, grounding_line_bmlt out of range for glissade dycore', GM_FATAL)
    end if

    if (model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_FLOATING_FRAC  .or. &
        model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_ZERO_GROUNDED) then
       if (model%options%which_ho_ground /= HO_GROUND_GLP_DELUXE) then
          write(message,'(a,i2)') 'This setting of grounding_line_bmlt requires grounding_line_param = ', &
               HO_GROUND_GLP_DELUXE
          call write_log(message)
          call write_log('Incompatible grounding_line_bmlt and grounding_line_param', GM_FATAL)
       endif
    endif

    if (model%options%which_ho_flotation_function < 0 .or. &
         model%options%which_ho_flotation_function >= size(flotation_function)) then
       call write_log('Error, flotation_function option out of range for glissade dycore', GM_FATAL)
    endif

    ! write the options

    call log_section('grounding_line')
    call log_option('grounding_line_param', model%options%which_ho_ground, &
         grounding_line_param(model%options%which_ho_ground))
    if (model%options%which_ho_ground == HO_GROUND_NO_GLP) then
       call log_option('fground_no_glp', model%options%which_ho_fground_no_glp, &
            fground_no_glp(model%options%which_ho_fground_no_glp))
    endif
    call log_option('grounding_line_bmlt', model%options%which_ho_ground_bmlt, &
         grounding_line_bmlt(model%options%which_ho_ground_bmlt))
    call log_option('flotation_function', model%options%which_ho_flotation_function, &
         flotation_function(model%options%which_ho_flotation_function))

    if (model%options%which_ho_ground == HO_GROUND_NO_GLP .and. &
        model%options%which_ho_flotation_function == HO_FLOTATION_FUNCTION_PATTYN) then
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
    call GetValue(section,'bmlt_float_method',                model%options%whichbmlt_float)
    call GetValue(section,'bmlt_float_thermal_forcing_param', model%options%bmlt_float_thermal_forcing_param)
    call GetValue(section,'bmlt_float_init',                  model%options%bmlt_float_init)
    call GetValue(section,'deltaT_ocn_method',                model%options%which_ho_deltaT_ocn)
    call GetValue(section,'enable_bmlt_anomaly',              model%options%enable_bmlt_anomaly)
    call GetValue(section,'basal_mass_balance',               model%options%basal_mbal)

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

    character(len=*), dimension(0:1), parameter :: basal_mass_balance = (/ &
         'not in continuity eqn', &
         'in continuity eqn    ' /)

    character(len=*), dimension(0:6), parameter :: bmlt_float_method = (/ &
         'none                                  ', &
         'MISMIP+ melt rate profile             ', &
         'constant melt rate                    ', &
         'depth-dependent melt rate             ', &
         'melt rate from external file          ', &
         'unused (formerly MISOMIP plume model) ', &
         'melt rate from thermal forcing        ' /)

    character(len=*), dimension(0:3), parameter :: bmlt_float_thermal_forcing_param = (/ &
         'quadratic function of thermal forcing     ', &
         'ISMIP6 local quadratic                    ', &
         'ISMIP6 nonlocal quadratic                 ', &
         'ISMIP6 nonlocal quadratic, slope-dependent' /)

    character(len=*), dimension(0:5), parameter :: deltaT_ocn_method = (/ &
         'deltaT_ocn = 0                              ', &
         'invert for 2D deltaT_ocn based on thck      ', &
         'read deltaT_ocn from external file          ', &
         'invert for basin-scale deltaT_ocn           ', &
         'calibrate deltaT_basin to match melt target ', &
         'set deltaT_ocn based on observed dthck_dt   ' /)

    ! check the options

    if (model%options%whichbmlt_float < 0 .or. model%options%whichbmlt_float >= size(bmlt_float_method)) then
       call write_log('Error, bmlt_float_method out of range',GM_FATAL)
    end if

    if (model%options%whichbmlt_float == 5) then   ! formerly the MISOMIP plume option, now removed
       call write_log('Error, bmlt_float_method = 5 (MISOMIP plume model) is no longer supported', GM_FATAL)
    endif

    if (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then
       if (model%options%bmlt_float_thermal_forcing_param < 0 .or. &
           model%options%bmlt_float_thermal_forcing_param >= size(bmlt_float_thermal_forcing_param)) then
          call write_log('Error, bmlt_float_thermal_forcing_param out of range', GM_FATAL)
       end if
    endif

    if (model%options%basal_mbal < 0 .or. model%options%basal_mbal >= size(basal_mass_balance)) then
       call write_log('Error, basal_mass_balance out of range',GM_FATAL)
    end if

    if (model%options%which_ho_deltaT_ocn < 0 .or. &
         model%options%which_ho_deltaT_ocn >= size(deltaT_ocn_method)) then
       call write_log('Error, deltaT_ocn_method out of range', GM_FATAL)
    end if

    if (model%options%which_ho_deltaT_ocn == HO_DELTAT_OCN_DTHCK_DT) then
       if (model%options%bmlt_float_thermal_forcing_param /= BMLT_FLOAT_TF_ISMIP6_LOCAL .and. &
           model%options%bmlt_float_thermal_forcing_param /= BMLT_FLOAT_TF_ISMIP6_NONLOCAL .and. &
           model%options%bmlt_float_thermal_forcing_param /= BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then
          write(message,*) 'deltaT_ocn_method dthck_dt option supported only for ISMIP6 bmlt_float schemes'
          call write_log(message, GM_FATAL)
       endif
    endif

    ! write the options and parameters

    call log_section('basal_melt')
    call log_option('bmlt_float_method', model%options%whichbmlt_float, &
         bmlt_float_method(model%options%whichbmlt_float))
    if (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then
       call log_option('bmlt_float_thermal_forcing_param', model%options%bmlt_float_thermal_forcing_param, &
            bmlt_float_thermal_forcing_param(model%options%bmlt_float_thermal_forcing_param))
       call log_flag('bmlt_float_init', model%options%bmlt_float_init, &
            'compute bmlt_float at initialization')
    endif
    call log_option('basal_mass_balance', model%options%basal_mbal, &
         basal_mass_balance(model%options%basal_mbal))

    if (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then
       call log_option('deltaT_ocn_method', model%options%which_ho_deltaT_ocn, &
            deltaT_ocn_method(model%options%which_ho_deltaT_ocn))
    elseif (model%options%which_ho_deltaT_ocn /= HO_DELTAT_OCN_NONE) then
       write(message,'(a,i0)') 'ignored; supported only for bmlt_float_method = ', BMLT_FLOAT_THERMAL_FORCING
       call log_option('deltaT_ocn_method', model%options%which_ho_deltaT_ocn, trim(message))
    endif

    call log_flag('enable_bmlt_anomaly', model%options%enable_bmlt_anomaly, 'bmlt_float anomaly forcing')
    if (model%options%enable_bmlt_anomaly) then
       call log_param('bmlt_anomaly_tstart', model%basal_melt%bmlt_anomaly_tstart, 'yr')
       call log_param('bmlt_anomaly_timescale', model%basal_melt%bmlt_anomaly_timescale, 'yr')
    endif

    if (model%options%whichbmlt_float /= BMLT_FLOAT_NONE .and. &
        model%options%whichbmlt_float /= BMLT_FLOAT_MISMIP) then
       call log_param('bmlt_cavity_h0', model%basal_melt%bmlt_cavity_h0, &
            'm; melt is reduced in cavities thinner than this; 0 = off')
    endif

    if (model%options%whichbmlt_float == BMLT_FLOAT_EXTERNAL) then
       call log_param('bmlt_float_factor', model%basal_melt%bmlt_float_factor, &
            'multiplies bmlt_float from the input file')
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_CONSTANT) then
       call log_param('bmlt_float_const', model%basal_melt%bmlt_float_const, 'm/yr')
       call log_param('bmlt_float_xlim', model%basal_melt%bmlt_float_xlim, 'm')
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_MISMIP) then
       call log_param('bmlt_float_omega', model%basal_melt%bmlt_float_omega, 'yr^-1')
       call log_param('bmlt_float_h0', model%basal_melt%bmlt_float_h0, 'm')
       call log_param('bmlt_float_z0', model%basal_melt%bmlt_float_z0, 'm')
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_DEPTH) then
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
    call GetValue(section,'calving_law',           model%options%whichcalving)
    call GetValue(section,'calving_init',          model%options%calving_init)
    call GetValue(section,'calving_domain',        model%options%calving_domain)
    call GetValue(section,'apply_calving_mask',    model%options%apply_calving_mask)
    call GetValue(section,'remove_icebergs',       model%options%remove_icebergs)
    call GetValue(section,'remove_isthmuses',      model%options%remove_isthmuses)
    call GetValue(section,'limit_marine_cliffs',   model%options%limit_marine_cliffs)
    call GetValue(section,'subgrid_calving_front', model%options%which_ho_calving_front)
    call GetValue(section,'calvingmip_domain',     model%options%which_ho_calvingmip_domain)

    ! calving parameters
    call GetValue(section,'marine_limit',       model%calving%marine_limit)
    call GetValue(section,'calving_fraction',   model%calving%calving_fraction)
    call GetValue(section,'calving_timescale',  model%calving%timescale)
    call GetValue(section,'calving_minthck',    model%calving%minthck)
    call GetValue(section,'dthck_dx_cf',        model%calving%dthck_dx_cf)
    call GetValue(section,'thck_effective_min', model%calving%thck_effective_min)
    call GetValue(section,'eigenconstant',      model%calving%eigenconstant)
    call GetValue(section,'tau_eigenconstant1', model%calving%tau_eigenconstant1)
    call GetValue(section,'tau_eigenconstant2', model%calving%tau_eigenconstant2)
    call GetValue(section,'stress_threshold',   model%calving%stress_threshold)
    call GetValue(section,'lateral_rate_min',   model%calving%lateral_rate_min)
    call GetValue(section,'effec_stress_min',   model%calving%effec_stress_min)
    call GetValue(section,'calving_length_scale', model%calving%length_scale)
    call GetValue(section,'damage_constant',    model%calving%damage_constant)
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
    ! Note: If remove_isthmuses = .true., remove_icebergs is set to .true. here (for stability).

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model
    character(len=500) :: message
    logical :: icebergs_reset   ! true if remove_icebergs is reset to .true.

    character(len=*), dimension(0:11), parameter :: calving_law = (/ &
         'no calving law                    ', &
         'remove all floating ice           ', &
         'remove fraction of floating ice   ', &
         'relaxed bedrock threshold         ', &
         'present bedrock threshold         ', &
         'calving based on grid location    ', &
         'prescribe CF advance/retreat rate ', &
         'ice thickness threshold           ', &
         'deterministic stress-based calving', &
         'stochastic stress-based calving   ', &
         'damage-based calving scheme       ', &
         'strain-rate (eigencalving)        ' /)

    character(len=*), dimension(0:1), parameter :: calving_init = (/ &
         'no calving at initialization    ', &
         'ice calves at initialization    ' /)

    character(len=*), dimension(0:1), parameter :: calving_domain = (/ &
         'calving only at the ocean edge             ',  &
         'calving in all cells where criterion is met'/)

    character(len=*), dimension(0:1), parameter :: subgrid_calving_front = (/ &
         'no subgrid calving front parameterization', &
         'subgrid calving front                    ' /)

    character(len=*), dimension(0:2), parameter :: calvingmip_domain = (/ &
         'none       ', &
         'circular   ', &
         'Thule      '  /)

    ! check the options

    if (model%options%whichcalving < 0 .or. model%options%whichcalving >= size(calving_law)) then
       call write_log('Error, calving_law out of range',GM_FATAL)
    end if

    if (model%options%calving_init < 0 .or. model%options%calving_init >= size(calving_init)) then
       call write_log('Error, calving_init out of range',GM_FATAL)
    end if

    if (model%options%calving_domain < 0 .or. model%options%calving_domain >= size(calving_domain)) then
       call write_log('Error, calving_domain out of range',GM_FATAL)
    end if

    if (model%options%which_ho_calving_front < 0 .or. &
         model%options%which_ho_calving_front >= size(subgrid_calving_front)) then
       call write_log('Error, calving front option out of range for glissade dycore', GM_FATAL)
    end if

    if (model%options%which_ho_calvingmip_domain < 0 .or. &
         model%options%which_ho_calvingmip_domain >= size(calvingmip_domain)) then
       call write_log('Error, calvingMIP domain option out of range for glissade dycore', GM_FATAL)
    end if

    ! for stability, remove icebergs when removing isthmuses
    icebergs_reset = .false.
    if (model%options%remove_isthmuses .and. .not.model%options%remove_icebergs) then
       model%options%remove_icebergs = .true.
       icebergs_reset = .true.
    endif

    if (model%options%whichcalving == CF_ADVANCE_RETREAT_RATE    .or.  &
        model%options%whichcalving == CALVING_THCK_THRESHOLD     .or.  &
        model%options%whichcalving == CALVING_STRESS             .or.  &
        model%options%whichcalving == CALVING_STRESS_STOCHASTIC  .or.  &
        model%options%whichcalving == EIGEN_CALVING              .or.  &
        model%options%whichcalving == CALVING_DAMAGE) then
       if (model%options%which_ho_calving_front == HO_CALVING_FRONT_NO_SUBGRID) then
          write(message,*) &
               'Calving option ', model%options%whichcalving, ' requires a subgrid calving front'
          call write_log(message, GM_FATAL)
       endif
    endif

    if (model%options%whichcalving == CALVING_THCK_THRESHOLD) then
       if (model%calving%minthck <= 0.0d0) then
          write(message,*) 'Error, this calving option needs calving_minthck > 0'
          call write_log(message, GM_FATAL)
       endif
    endif

    if (model%options%which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then
       if (.not.model%options%remove_icebergs) then
          write(message,*) 'Set remove_icebergs = T for stability when using subgrid calving_front scheme'
          call write_log(message, GM_FATAL)
       endif
    endif

    ! write the options and parameters

    call log_section('calving')
    call log_option('calving_law', model%options%whichcalving, calving_law(model%options%whichcalving))
    call log_option('calving_init', model%options%calving_init, calving_init(model%options%calving_init))
    call log_option('calving_domain', model%options%calving_domain, calving_domain(model%options%calving_domain))

    ! parameters of the calving law
    if (model%options%whichcalving == CALVING_FLOAT_FRACTION) then
       call log_param('calving_fraction', model%calving%calving_fraction, 'ice fraction lost in calving')
    elseif (model%options%whichcalving == CALVING_RELX_THRESHOLD .or.  &
            model%options%whichcalving == CALVING_TOPG_THRESHOLD) then
       call log_param('marine_limit', model%calving%marine_limit, 'm; marine depth limit')
    elseif (model%options%whichcalving == CALVING_THCK_THRESHOLD) then
       call log_param('calving_minthck', model%calving%minthck, 'm')
    elseif (model%options%whichcalving == EIGEN_CALVING) then
       call log_param('eigenconstant', model%calving%eigenconstant, 'm')
    elseif (model%options%whichcalving == CALVING_STRESS) then
       call log_param('tau_eigenconstant1', model%calving%tau_eigenconstant1)
       call log_param('tau_eigenconstant2', model%calving%tau_eigenconstant2)
       call log_param('stress_threshold', model%calving%stress_threshold, 'Pa')
       call log_param('lateral_rate_min', model%calving%lateral_rate_min, 'm/yr')
    elseif (model%options%whichcalving == CALVING_STRESS_STOCHASTIC) then
       call log_param('tau_eigenconstant1', model%calving%tau_eigenconstant1)
       call log_param('tau_eigenconstant2', model%calving%tau_eigenconstant2)
       call log_param('stress_threshold', model%calving%stress_threshold, 'Pa')
       call log_param('effec_stress_min', model%calving%effec_stress_min)
       call log_param('calving_length_scale', model%calving%length_scale, 'm')
    elseif (model%options%whichcalving == CALVING_DAMAGE) then
       call log_param('tau_eigenconstant1', model%calving%tau_eigenconstant1)
       call log_param('tau_eigenconstant2', model%calving%tau_eigenconstant2)
    elseif (model%options%whichcalving == CF_ADVANCE_RETREAT_RATE) then
       call log_param('cf_advance_retreat_amplitude', model%calving%cf_advance_retreat_amplitude, 'm/yr')
       call log_param('cf_advance_retreat_period', model%calving%cf_advance_retreat_period, 'yr')
    elseif (model%options%whichcalving == CALVING_GRID_MASK) then
       call log_param('calving_front_x', model%calving%calving_front_x, 'm; 0 = not used')
       call log_param('calving_front_y', model%calving%calving_front_y, 'm; 0 = not used')
       call log_param('calving_front_radius', model%calving%calving_front_radius, 'm; 0 = not used')
    endif
    if (model%options%whichcalving /= CALVING_NONE) then
       call log_param('calving_timescale', model%calving%timescale, 'yr; 0 = calve immediately')
    endif

    call log_option('subgrid_calving_front', model%options%which_ho_calving_front, &
         subgrid_calving_front(model%options%which_ho_calving_front))
    if (model%options%which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then
       call log_param('dthck_dx_cf', model%calving%dthck_dx_cf)
       call log_param('thck_effective_min', model%calving%thck_effective_min, 'm')
    endif
    if (model%options%which_ho_calvingmip_domain /= HO_CALVINGMIP_DOMAIN_NONE) then
       call log_option('calvingmip_domain', model%options%which_ho_calvingmip_domain, &
            calvingmip_domain(model%options%which_ho_calvingmip_domain))
    end if

    call log_flag('apply_calving_mask', model%options%apply_calving_mask)
    if (icebergs_reset) then
       call log_flag('remove_icebergs', model%options%remove_icebergs, &
            'reset from .false. for stability when remove_isthmuses = .true.')
    else
       call log_flag('remove_icebergs', model%options%remove_icebergs)
    endif
    call log_flag('remove_isthmuses', model%options%remove_isthmuses)
    if (model%options%remove_icebergs .or. model%options%remove_isthmuses) then
       call log_param('f_ground_threshold', model%calving%f_ground_threshold)
    endif
    call log_flag('limit_marine_cliffs', model%options%limit_marine_cliffs)
    if (model%options%limit_marine_cliffs) then
       call log_param('taumax_cliff', model%calving%taumax_cliff, 'Pa')
       call log_param('cliff_timescale', model%calving%cliff_timescale, 'yr')
    endif

    if (model%options%whichcalving == CALVING_FLOAT_FRACTION) then
       write(message,*) 'WARNING: calving float fraction option deprecated with Glissade_dycore; set calving_timescale instead'
       call write_log(message, GM_WARNING)
    endif

  end subroutine print_calving

!--------------------------------------------------------------------------------

  subroutine handle_lateral_melt(section, model)

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    call GetValue(section, 'lateral_melt_method', model%options%which_lateral_melt)
    call GetValue(section, 'melt_rate_const', model%lateral_melt%melt_rate_const)
    call GetValue(section, 'melt_factor', model%lateral_melt%melt_factor)
    call GetValue(section, 'runoff_from_ablation', model%lateral_melt%runoff_from_ablation)
    call GetValue(section, 'thermal_forcing_avg_3d_to_2d', model%lateral_melt%thermal_forcing_avg_3d_to_2d)
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

    character(len=*), dimension(0:2), parameter :: lateral_melt_method = (/ &
         'no lateral melt at marine margin         ', &
         'constant lateral melt rate               ', &
         'ISMIP lateral melt from TF and discharge ' /)

    if (model%options%which_lateral_melt < 0 .or. model%options%which_lateral_melt >= size(lateral_melt_method)) then
       call write_log('Error, lateral melt option out of range', GM_FATAL)
    endif

    call log_section('lateral_melt')
    call log_option('lateral_melt_method', model%options%which_lateral_melt, &
         lateral_melt_method(model%options%which_lateral_melt))

    if (model%options%which_lateral_melt == LATERAL_MELT_CONSTANT) then
       call log_param('melt_rate_const', model%lateral_melt%melt_rate_const, 'm/yr')
    elseif (model%options%which_lateral_melt == LATERAL_MELT_ISMIP) then
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

    call GetValue(section,'remove_ice_caps', model%options%remove_ice_caps)
    call GetValue(section,'force_retreat',   model%options%force_retreat)

  end subroutine handle_removal

!--------------------------------------------------------------------------------

  subroutine print_removal(model)

    ! Ice removal options

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    character(len=*), dimension(0:2), parameter :: force_retreat = (/ &
         'no forced retreat                                  ', &
         'force retreat of all ice with retreat mask         ', &
         'force retreat of floating ice with retreat mask    ' /)

    if (model%options%force_retreat < 0 .or. model%options%force_retreat >= size(force_retreat)) then
       call write_log('Error, force_retreat out of range', GM_FATAL)
    end if

    call log_section('removal')
    call log_flag('remove_ice_caps', model%options%remove_ice_caps, 'removed ice goes to the removal flux')
    call log_option('force_retreat', model%options%force_retreat, force_retreat(model%options%force_retreat))
    if (model%options%force_retreat == FORCE_RETREAT_FLOATING_ICE .and. &
        .not.model%options%remove_isthmuses) then
       call log_comment('Warning: force_retreat = 2 can be unstable when remove_isthmuses = .false.')
    endif

  end subroutine print_removal

!--------------------------------------------------------------------------------

  subroutine glissade_get_zatm(model,config)

    ! Read atmosphere grid information, if present, from the config file.
    ! Called after glissade_readconfig and before glissade_printconfig (see print_grid_atm)

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

    ! Check for section [grid_atm]
    call GetSection(config,section,'grid_atm')

    if (associated(section)) then

       call GetValue(section,'nzatm',model%climate%nzatm)

       if (.not.associated(model%climate%zatm)) allocate(model%climate%zatm(model%climate%nzatm))

       if (size(model%climate%zatm) /= model%climate%nzatm) then
          write(message,*) 'size of zatm : ', size(model%climate%zatm)
          call write_log(trim(message))
          call write_log ( 'zatm has the wrong size; must equal nzatm', GM_FATAL)
       end if

       ! Load the zatm levels from the config file.
       ! For example, the config file could have a section like this:
       ! [grid_atm]
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
             write(message,*) 'Must have zatm increasing with increasing k in the [grid_atm] section'
             call write_log(message, GM_FATAL)
          endif
       enddo

    elseif (model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ .or. &
            model%options%artm_input_function == ARTM_INPUT_FUNCTION_XYZ) then

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

    endif   ! associated(section)

  end subroutine glissade_get_zatm

!--------------------------------------------------------------------------------

  subroutine print_grid_atm(model, config)

    ! Atmosphere grid ([grid_atm]; read or set by glissade_get_zatm), written if there is more than one level

    use glissade_types
    use glimmer_config
    use glissade_log_config
    implicit none
    type(glissade_global_type)   :: model
    type(ConfigSection), pointer :: config

    type(ConfigSection), pointer :: section

    if (model%climate%nzatm <= 1) return

    call GetSection(config,section,'grid_atm')
    call log_section('grid_atm')
    if (.not.associated(section)) then
       call log_comment('no [grid_atm] section; using the standard CESM elevation classes')
    endif
    call log_param('nzatm', model%climate%nzatm)
    call log_list('zatm', model%climate%zatm, 'm')

  end subroutine print_grid_atm

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
    call GetValue(section,'smb_input',                       model%options%smb_input)
    call GetValue(section,'smb_input_function',              model%options%smb_input_function)
    call GetValue(section,'artm_input_function',             model%options%artm_input_function)
    call GetValue(section,'enable_smb_anomaly',              model%options%enable_smb_anomaly)
    call GetValue(section,'enable_artm_anomaly',             model%options%enable_artm_anomaly)
    call GetValue(section,'enable_snow_anomaly',             model%options%enable_snow_anomaly)
    call GetValue(section,'enable_precip_anomaly',           model%options%enable_precip_anomaly)
    call GetValue(section,'overwrite_acab',                  model%options%overwrite_acab)
    call GetValue(section,'enable_acab_dthck_dt_correction', model%options%enable_acab_dthck_dt_correction)
    call GetValue(section,'block_inception',                 model%options%block_inception)

    ! parameters of the internal SMB and artm functions
    call GetValue(section,'t_lapse',            model%climate%t_lapse)
    call GetValue(section,'degree_factor',      model%climate%degree_factor)
    call GetValue(section,'tmlt',               model%climate%tmlt)
    call GetValue(section,'snow_threshold_min', model%climate%snow_threshold_min)
    call GetValue(section,'snow_threshold_max', model%climate%snow_threshold_max)
    call GetValue(section,'smb_factor',         model%climate%smb_factor)

    ! anomaly and overwrite parameters
    call GetValue(section,'smb_anomaly_tstart',     model%climate%smb_anomaly_tstart)
    call GetValue(section,'smb_anomaly_timescale',  model%climate%smb_anomaly_timescale)
    call GetValue(section,'overwrite_acab_value',   model%climate%overwrite_acab_value)
    call GetValue(section,'overwrite_acab_minthck', model%climate%overwrite_acab_minthck)
    call GetValue(section,'artm_anomaly_const',     model%climate%artm_anomaly_const)
    call GetValue(section,'artm_anomaly_tstart',    model%climate%artm_anomaly_tstart)
    call GetValue(section,'artm_anomaly_timescale', model%climate%artm_anomaly_timescale)

  end subroutine handle_climate

!--------------------------------------------------------------------------------

  subroutine print_climate(model)

    ! Surface climate (SMB and surface temperature) options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model
    character(len=500) :: message

    character(len=*), dimension(0:1), parameter :: smb_input = (/ &
         'SMB input in units of m/yr ice  ', &
         'SMB input in units of mm/yr w.e.' /)

    character(len=*), dimension(0:3), parameter :: smb_input_function = (/ &
         'SMB input as function of (x,y)              ', &
         'SMB and d(SMB)/dz input as function of (x,y)', &
         'SMB input as function of (x,y,z)            ', &
         'SMB computed with positive-degree scheme    ' /)

    character(len=*), dimension(0:3), parameter :: artm_input_function = (/ &
         'artm input as function of (x,y)               ', &
         'artm and d(artm)/dz input as function of (x,y)', &
         'artm input as function of (x,y,z)             ', &
         'artm input as function of (x,y) w/ lapse rate ' /)

    character(len=*), dimension(0:3), parameter :: overwrite_acab = (/ &
         'do not overwrite acab anywhere            ', &
         'overwrite acab where input acab = 0       ', &
         'overwrite acab where input thck <= minthck', &
         'overwrite acab based on input mask        ' /)

    ! check the options

    if (model%options%smb_input < 0 .or. model%options%smb_input >= size(smb_input)) then
       call write_log('Error, smb_input option out of range',GM_FATAL)
    end if

    if (model%options%smb_input_function < 0 .or. model%options%smb_input_function >= size(smb_input_function)) then
       call write_log('Error, smb_input_function option out of range',GM_FATAL)
    end if

    if (model%options%smb_input == SMB_INPUT_MYR_ICE) then
       if (model%options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ .or. &
           model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ .or. &
           model%climate%smb_factor /= 1.0d0 .or. &
           model%options%enable_smb_anomaly .or. &
           model%options%enable_glaciers) then
          write(message,*) 'Error, input SMB must have units of mm/yr w.e., smb_input =', SMB_INPUT_MMYR_WE
          call write_log(message, GM_FATAL)
       endif
    endif

    if (model%options%artm_input_function < 0 .or. model%options%artm_input_function >= size(artm_input_function)) then
       call write_log('Error, artm_input_function option out of range',GM_FATAL)
    end if

    if (model%options%overwrite_acab < 0 .or. model%options%overwrite_acab >= size(overwrite_acab)) then
       call write_log('Error, overwrite_acab option out of range',GM_FATAL)
    end if

    ! write the options and parameters

    call log_section('climate')
    call log_option('smb_input', model%options%smb_input, smb_input(model%options%smb_input))
    call log_option('smb_input_function', model%options%smb_input_function, &
         smb_input_function(model%options%smb_input_function))
    if (model%options%smb_input_function == SMB_INPUT_FUNCTION_PDD) then
       call log_param('degree_factor', model%climate%degree_factor, 'mm/yr w.e./deg C')
       call log_param('tmlt', model%climate%tmlt, 'deg C; melt threshold')
       call log_param('snow_threshold_min', model%climate%snow_threshold_min, 'deg C')
       call log_param('snow_threshold_max', model%climate%snow_threshold_max, 'deg C')
    endif
    call log_option('artm_input_function', model%options%artm_input_function, &
         artm_input_function(model%options%artm_input_function))
    if (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_LAPSE) then
       call log_param('t_lapse', model%climate%t_lapse, 'deg C/m; artm lapse rate')
    endif
    call log_param('smb_factor', model%climate%smb_factor, 'multiplies the input SMB')

    call log_flag('enable_smb_anomaly', model%options%enable_smb_anomaly)
    if (model%options%enable_smb_anomaly) then
       call log_param('smb_anomaly_tstart', model%climate%smb_anomaly_tstart, 'yr')
       call log_param('smb_anomaly_timescale', model%climate%smb_anomaly_timescale, 'yr')
    endif
    call log_flag('enable_artm_anomaly', model%options%enable_artm_anomaly)
    if (model%options%enable_artm_anomaly) then
       call log_param('artm_anomaly_const', model%climate%artm_anomaly_const, 'deg C')
       call log_param('artm_anomaly_tstart', model%climate%artm_anomaly_tstart, 'yr')
       call log_param('artm_anomaly_timescale', model%climate%artm_anomaly_timescale, 'yr')
    endif
    call log_flag('enable_snow_anomaly', model%options%enable_snow_anomaly)
    call log_flag('enable_precip_anomaly', model%options%enable_precip_anomaly)

    call log_option('overwrite_acab', model%options%overwrite_acab, overwrite_acab(model%options%overwrite_acab))
    if (model%options%overwrite_acab /= OVERWRITE_ACAB_NONE) then
       call log_param('overwrite_acab_value', model%climate%overwrite_acab_value, 'm/yr')
       if (model%options%overwrite_acab == OVERWRITE_ACAB_THCKMIN) then
          call log_param('overwrite_acab_minthck', model%climate%overwrite_acab_minthck, 'm')
       endif
    endif
    call log_flag('enable_acab_dthck_dt_correction', model%options%enable_acab_dthck_dt_correction, &
         'correct acab based on dthck_dt_obs')
    call log_flag('block_inception', model%options%block_inception, &
         'block inception outside the main ice sheet')

  end subroutine print_climate

!-------------------------------------------------------------------------

  subroutine glissade_get_zocn(model,config)

    ! Read ocean grid information, if present, from the config file.
    ! Called after glissade_readconfig and before glissade_printconfig (see print_grid_ocn)

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

    ! Check for section [grid_ocn]
    call GetSection(config,section,'grid_ocn')

    if (associated(section)) then

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
                write(message,*) 'Must have zocn decreasing with increasing k in the [grid_ocn] section'
                call write_log(message, GM_FATAL)
             endif
          enddo
       endif

    else    ! no 'grid_ocn' section

       ! Allocate and initialize model%ocean_data%zocn.
       ! Note: Coupled CESM runs require this array to be allocated, even if CISM
       !       is not receiving realistic ocean forcing fields (as of Jan. 2021).
       model%ocean_data%nzocn = get_nzocn(model)
       if (.not.associated(model%ocean_data%zocn)) &
            allocate(model%ocean_data%zocn(model%ocean_data%nzocn))
       model%ocean_data%zocn(:) = 0.0d0

       if (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then
          write(message,*) 'Must have a [grid_ocn] section to use the thermal forcing option'
          call write_log(message, GM_FATAL)
       endif

    endif   ! associated(section)

  end subroutine glissade_get_zocn

!--------------------------------------------------------------------------------

  subroutine print_grid_ocn(model, config)

    ! Ocean grid ([grid_ocn]; read by glissade_get_zocn), written if the section is present

    use glissade_types
    use glimmer_config
    use glissade_log_config
    implicit none
    type(glissade_global_type)   :: model
    type(ConfigSection), pointer :: config

    type(ConfigSection), pointer :: section

    call GetSection(config,section,'grid_ocn')
    if (.not.associated(section)) return

    call log_section('grid_ocn')
    call log_param('nzocn', model%ocean_data%nzocn)
    if (model%ocean_data%dzocn > 0.0d0) then
       call log_param('dzocn', model%ocean_data%dzocn, 'm; uniform level spacing')
       call log_list('zocn', model%ocean_data%zocn, 'm; computed from dzocn')
    else
       call log_list('zocn', model%ocean_data%zocn, 'm')
    endif
    call log_param('nbasin', model%ocean_data%nbasin)

  end subroutine print_grid_ocn

!--------------------------------------------------------------------------------

  subroutine handle_ocean_data(section, model)

    ! Ocean data (thermal forcing) for the thermal-forcing basal melt schemes.

    use glimmer_config
    use glissade_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glissade_global_type)   :: model

    ! options
    call GetValue(section,'ocean_data_domain',      model%options%ocean_data_domain)
    call GetValue(section,'ocean_data_extrapolate', model%options%ocean_data_extrapolate)

    ! parameters
    call GetValue(section,'gamma0',                            model%ocean_data%gamma0)
    call GetValue(section,'thermal_forcing_basin_min',         model%ocean_data%thermal_forcing_basin_min)
    call GetValue(section,'thermal_forcing_basin_max',         model%ocean_data%thermal_forcing_basin_max)
    call GetValue(section,'thermal_forcing_anomaly',           model%ocean_data%thermal_forcing_anomaly)
    call GetValue(section,'thermal_forcing_anomaly_tstart',    model%ocean_data%thermal_forcing_anomaly_tstart)
    call GetValue(section,'thermal_forcing_anomaly_timescale', model%ocean_data%thermal_forcing_anomaly_timescale)
    call GetValue(section,'thermal_forcing_anomaly_basin',     model%ocean_data%thermal_forcing_anomaly_basin)

  end subroutine handle_ocean_data

!--------------------------------------------------------------------------------

  subroutine print_ocean_data(model)

    ! Ocean data (thermal forcing) options and parameters.
    ! The section is written only for the thermal-forcing basal melt scheme, which uses it.

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    character(len=*), dimension(0:2), parameter :: ocean_data_domain = (/ &
         'ocean data computed internally by CISM', &
         'ocean data read from external file    ', &
         'ocean data from coupler via Glad      ' /)

    character(len=*), dimension(0:1), parameter :: ocean_data_extrapolate = (/ &
         'ocean data not extrapolated to cavities', &
         'ocean data extrapolated to cavities    ' /)

    if (model%options%whichbmlt_float /= BMLT_FLOAT_THERMAL_FORCING) return

    ! check the options

    if (model%options%ocean_data_domain < 0 .or. &
        model%options%ocean_data_domain >= size(ocean_data_domain)) then
       call write_log('Error, ocean_data_domain out of range', GM_FATAL)
    end if

    if (model%options%ocean_data_extrapolate < 0 .or. &
        model%options%ocean_data_extrapolate >= size(ocean_data_extrapolate)) then
       call write_log('Error, ocean_data_extrapolate out of range', GM_FATAL)
    end if

    ! write the options and parameters

    call log_section('ocean_data')
    call log_option('ocean_data_domain', model%options%ocean_data_domain, &
         ocean_data_domain(model%options%ocean_data_domain))
    call log_option('ocean_data_extrapolate', model%options%ocean_data_extrapolate, &
         ocean_data_extrapolate(model%options%ocean_data_extrapolate))
    call log_param('gamma0', model%ocean_data%gamma0, 'm/yr; melt rate coefficient')

    if (model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or. &
        model%options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then
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

    call GetValue(section,'geothermal_method', model%options%gthf)
    call GetValue(section,'geothermal_flux',   model%paramets%geot)

  end subroutine handle_geothermal

!--------------------------------------------------------------------------------

  subroutine print_geothermal(model)

    ! Geothermal heat flux option and parameter

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none

    type(glissade_global_type)  :: model

    ! NOTE: Set geothermal_method = 1 in the config file to read the geothermal heat flux from an input file.
    !       Otherwise it will be overwritten, even if the 'bheatflx' field is present.

    character(len=*), dimension(0:1), parameter :: geothermal_method = (/ &
         'uniform geothermal flux         ', &
         'read flux from file, if present ' /)

    if (model%options%gthf < 0 .or. model%options%gthf >= size(geothermal_method)) then
       call write_log('Error, geothermal flux option out of range',GM_FATAL)
    end if

    call log_section('geothermal')
    call log_option('geothermal_method', model%options%gthf, geothermal_method(model%options%gthf))
    call log_param('geothermal_flux', model%paramets%geot, 'W/m^2')

  end subroutine print_geothermal

!--------------------------------------------------------------------------------

  subroutine handle_isostasy(section, model)
    use glimmer_config
    use glissade_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glissade_global_type)  :: model

    call GetValue(section,'isostasy',model%options%isostasy)
    call GetValue(section,'lithosphere_method',model%isostasy%lithosphere)
    call GetValue(section,'asthenosphere_method',model%isostasy%asthenosphere)
    call GetValue(section,'input_topg_method',model%isostasy%which_relaxed)
    call GetValue(section,'tau_relax_const',model%isostasy%tau_relax_const)
    call GetValue(section,'load_update_interval',model%isostasy%load_update_interval)
    call GetValue(section,'flexural_rigidity',model%isostasy%rbel%d)

  end subroutine handle_isostasy

!--------------------------------------------------------------------------------

  subroutine print_isostasy(model)

    ! Isostasy options and parameters

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none
    type(glissade_global_type)  :: model

    character(len=*), dimension(0:1), parameter :: isostasy = (/ &
         'no isostasy calculation         ', &
         'compute isostasy with model     ' /)

    character(len=*), dimension(0:1), parameter :: lithosphere_method = (/ &
         'local lithosphere         ', &
         'elastic lithosphere       ' /)

    character(len=*), dimension(0:2), parameter :: asthenosphere_method = (/ &
         'fluid asthenosphere                ', &
         'relaxing with constant timescale   ', &
         'relaxing with lat varying timescale' /)

    character(len=*), dimension(0:2), parameter :: input_topg_method = (/ &
         'read both topg and relx from input file', &
         'set relx to input topg                 ', &
         'compute relx = input topg + load       ' /)

    ! check the options

    if (model%options%isostasy < 0 .or. model%options%isostasy >= size(isostasy)) then
       call write_log('Error, isostasy option out of range',GM_FATAL)
    end if

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       if (model%isostasy%lithosphere < 0 .or. model%isostasy%lithosphere >= size(lithosphere_method)) then
          call write_log('Error, lithosphere_method out of range', GM_FATAL)
       endif
       if (model%isostasy%asthenosphere < 0 .or. model%isostasy%asthenosphere >= size(asthenosphere_method)) then
          call write_log('Error, asthenosphere_method out of range', GM_FATAL)
       endif
       if (model%isostasy%which_relaxed < 0 .or. model%isostasy%which_relaxed >= size(input_topg_method)) then
          call write_log('Error, input_topg_method out of range', GM_FATAL)
       endif
    endif

    ! write the options and parameters

    call log_section('isostasy')
    call log_option('isostasy', model%options%isostasy, isostasy(model%options%isostasy))

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call log_option('lithosphere_method', model%isostasy%lithosphere, &
            lithosphere_method(model%isostasy%lithosphere))
       if (model%isostasy%lithosphere == LITHOSPHERE_ELASTIC) then
          call log_param('flexural_rigidity', model%isostasy%rbel%d, 'N m')
          call log_param('load_update_interval', model%isostasy%load_update_interval, 'yr')
       end if
       call log_option('asthenosphere_method', model%isostasy%asthenosphere, &
            asthenosphere_method(model%isostasy%asthenosphere))
       if (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_CONST) then
          call log_param('tau_relax_const', model%isostasy%tau_relax_const, 'yr; relaxation timescale')
       endif
       call log_option('input_topg_method', model%isostasy%which_relaxed, &
            input_topg_method(model%isostasy%which_relaxed))
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
    call GetValue(section,'inversion_thck_flotation_buffer', model%inversion%thck_flotation_buffer)
    call GetValue(section,'inversion_thck_threshold',        model%inversion%thck_threshold)
    call GetValue(section,'inversion_damping_factor',        model%inversion%damping_factor)

    ! basal friction inversion
    call GetValue(section,'inversion_babc_timescale',    model%inversion%babc_timescale)
    call GetValue(section,'inversion_babc_thck_scale',   model%inversion%babc_thck_scale)
    call GetValue(section,'inversion_babc_length_scale', model%inversion%babc_length_scale)
    call GetValue(section,'inversion_babc_relax_factor', model%inversion%babc_relax_factor)

    ! deltaT_ocn inversion
    call GetValue(section,'inversion_deltaT_ocn_timescale',    model%inversion%deltaT_ocn_timescale)
    call GetValue(section,'inversion_deltaT_ocn_thck_scale',   model%inversion%deltaT_ocn_thck_scale)
    call GetValue(section,'inversion_deltaT_ocn_length_scale', model%inversion%deltaT_ocn_length_scale)
    call GetValue(section,'inversion_deltaT_ocn_temp_scale',   model%inversion%deltaT_ocn_temp_scale)
    call GetValue(section,'inversion_deltaT_basin_relax',      model%inversion%deltaT_basin_relax)

    ! basin-scale deltaT_ocn inversion
    call GetValue(section,'inversion_basin_flotation_threshold',    model%inversion%basin_flotation_threshold)
    call GetValue(section,'inversion_basin_mass_correction',        model%inversion%basin_mass_correction)
    call GetValue(section,'inversion_basin_number_mass_correction', model%inversion%basin_number_mass_correction)

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

    babc_inversion = (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or. &
                      model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION_BASIN .or. &
                      model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION .or. &
                      model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION_BASIN)

    deltaT_inversion = (model%options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION .or. &
                        model%options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION_BASIN)

    glacier_inversion = (model%options%enable_glaciers .and. &
                         model%glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION)

    if (.not. (babc_inversion .or. deltaT_inversion .or. glacier_inversion)) return

    call log_section('inversion')
    if (babc_inversion .or. deltaT_inversion) then
       call log_param('inversion_thck_flotation_buffer', model%inversion%thck_flotation_buffer, 'm')
       call log_param('inversion_thck_threshold', model%inversion%thck_threshold, 'm')
       call log_param('inversion_damping_factor', model%inversion%damping_factor)
    endif

    if (babc_inversion .or. glacier_inversion) then
       call log_param('inversion_babc_thck_scale', model%inversion%babc_thck_scale, 'm')
       call log_param('inversion_babc_timescale', model%inversion%babc_timescale, 'yr')
       call log_param('inversion_babc_length_scale', model%inversion%babc_length_scale, 'm')
       call log_param('inversion_babc_relax_factor', model%inversion%babc_relax_factor)
    endif

    if (deltaT_inversion) then
       call log_param('inversion_deltaT_ocn_thck_scale', model%inversion%deltaT_ocn_thck_scale, 'm')
       call log_param('inversion_deltaT_ocn_timescale', model%inversion%deltaT_ocn_timescale, 'yr')
       call log_param('inversion_deltaT_ocn_temp_scale', model%inversion%deltaT_ocn_temp_scale, 'deg C')
       if (model%options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION_BASIN) then
          call log_param('inversion_basin_flotation_threshold', model%inversion%basin_flotation_threshold, 'm')
          call log_param('inversion_deltaT_basin_relax', model%inversion%deltaT_basin_relax, &
               'deg C; value toward which deltaT_basin relaxes')
          call log_param('inversion_basin_number_mass_correction', model%inversion%basin_number_mass_correction, &
               '0 = no mass correction')
          call log_param('inversion_basin_mass_correction', model%inversion%basin_mass_correction, 'Gt')
       else   ! 2D deltaT_ocn inversion
          call log_param('inversion_deltaT_ocn_length_scale', model%inversion%deltaT_ocn_length_scale, 'm')
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

    call GetValue(section,'enable_glaciers',         model%options%enable_glaciers)
    call GetValue(section,'mu_star_method',          model%glacier%set_mu_star)
    call GetValue(section,'alpha_snow_method',       model%glacier%set_alpha_snow)
    call GetValue(section,'glacier_powerlaw_c_method', model%glacier%set_powerlaw_c)
    call GetValue(section,'snow_method',             model%glacier%snow_calc)
    call GetValue(section,'scale_area',              model%glacier%scale_area)
    call GetValue(section,'length_scale_factor',     model%glacier%length_scale_factor)
    call GetValue(section,'glacier_tmlt',            model%glacier%tmlt)
    call GetValue(section,'mu_star_const',           model%glacier%mu_star_const)
    call GetValue(section,'mu_star_min',             model%glacier%mu_star_min)
    call GetValue(section,'mu_star_max',             model%glacier%mu_star_max)
    call GetValue(section,'alpha_snow_const',        model%glacier%alpha_snow_const)
    call GetValue(section,'alpha_snow_min',          model%glacier%alpha_snow_min)
    call GetValue(section,'alpha_snow_max',          model%glacier%alpha_snow_max)
    call GetValue(section,'beta_artm_max',           model%glacier%beta_artm_max)
    call GetValue(section,'beta_artm_increment',     model%glacier%beta_artm_increment)
    call GetValue(section,'glacier_snow_threshold_min', model%glacier%snow_threshold_min)
    call GetValue(section,'glacier_snow_threshold_max', model%glacier%snow_threshold_max)
    call GetValue(section,'baseline_date',           model%glacier%baseline_date)
    call GetValue(section,'rgi_date',                model%glacier%rgi_date)
    call GetValue(section,'recent_date',             model%glacier%recent_date)
    call GetValue(section,'diagnostic_minthck',      model%glacier%diagnostic_minthck)
    call GetValue(section,'redistribute_advanced_ice',  model%glacier%redistribute_advanced_ice)
    call GetValue(section,'thinning_rate_advanced_ice', model%glacier%thinning_rate_advanced_ice)
    call GetValue(section,'smb_weight_advanced_ice', model%glacier%smb_weight_advanced_ice)

  end subroutine handle_glaciers

!--------------------------------------------------------------------------------

  subroutine print_glaciers(model)

    ! Glacier options and parameters.
    ! Note: If scale_area = .false., length_scale_factor is reset to 1 here.
    ! The parameters of the glacier powerlaw_c inversion are written in [inversion].

    use glissade_types
    use glimmer_log
    use glissade_log_config

    implicit none
    type(glissade_global_type)  :: model
    real(dp) :: length_scale_factor_in   ! length_scale_factor before a reset

    character(len=*), dimension(0:2), parameter :: mu_star_method = (/ &
         'spatially uniform glacier parameter mu_star', &
         'glacier-specific mu_star found by inversion', &
         'glacier-specific mu_star read from file    ' /)

    character(len=*), dimension(0:2), parameter :: alpha_snow_method = (/ &
         'spatially uniform glacier parameter alpha_snow', &
         'glacier-specific alpha_snow found by inversion', &
         'glacier-specific alpha_snow read from file    ' /)

    character(len=*), dimension(0:2), parameter :: glacier_powerlaw_c_method = (/ &
         'spatially uniform glacier parameter Cp', &
         'Cp for glaciers found by inversion    ', &
         'Cp for glaciers read from file        ' /)

    character(len=*), dimension(0:1), parameter :: snow_method = (/ &
         'read in snowfall rate directly            ', &
         'compute snowfall rate from precip and artm' /)

    if (.not.model%options%enable_glaciers) then
       call log_section('glaciers')
       call log_flag('enable_glaciers', model%options%enable_glaciers)
       return
    endif

    ! check the options

    if (model%glacier%set_mu_star < 0 .or. &
        model%glacier%set_mu_star >= size(mu_star_method)) then
       call write_log('Error, glacier_set_mu_star option out of range', GM_FATAL)
    end if

    if (model%glacier%set_alpha_snow < 0 .or. &
        model%glacier%set_alpha_snow >= size(alpha_snow_method)) then
       call write_log('Error, glacier_set_alpha_snow option out of range', GM_FATAL)
    end if

    if (model%glacier%set_powerlaw_c < 0 .or. &
        model%glacier%set_powerlaw_c >= size(glacier_powerlaw_c_method)) then
       call write_log('Error, glacier_set_powerlaw_c option out of range', GM_FATAL)
    end if

    if (model%glacier%snow_calc < 0 .or. &
        model%glacier%snow_calc >= size(snow_method)) then
       call write_log('Error, glacier_snow_calc option out of range', GM_FATAL)
    end if

    if (model%glacier%set_mu_star /= GLACIER_MU_STAR_INVERSION) then
       if (model%glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
          call write_log('Error, must invert for mu_star if inverting for alpha_snow', GM_FATAL)
       elseif (model%glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          call write_log('Error, must invert for mu_star if inverting for powerlaw_c', GM_FATAL)
       endif
    endif

    ! length_scale_factor is used only with scale_area
    length_scale_factor_in = model%glacier%length_scale_factor
    if (.not.model%glacier%scale_area) model%glacier%length_scale_factor = 1.0d0

    ! write the options and parameters

    call log_section('glaciers')
    call log_flag('enable_glaciers', model%options%enable_glaciers, 'glacier tracking and tuning')
    call log_option('mu_star_method', model%glacier%set_mu_star, mu_star_method(model%glacier%set_mu_star))
    call log_option('alpha_snow_method', model%glacier%set_alpha_snow, &
         alpha_snow_method(model%glacier%set_alpha_snow))
    call log_option('glacier_powerlaw_c_method', model%glacier%set_powerlaw_c, &
         glacier_powerlaw_c_method(model%glacier%set_powerlaw_c))
    call log_option('snow_method', model%glacier%snow_calc, snow_method(model%glacier%snow_calc))
    if (model%glacier%snow_calc == GLACIER_SNOW_CALC_PRECIP_ARTM) then
       call log_param('glacier_snow_threshold_min', model%glacier%snow_threshold_min, 'deg C')
       call log_param('glacier_snow_threshold_max', model%glacier%snow_threshold_max, 'deg C')
    endif
    call log_param('glacier_tmlt', model%glacier%tmlt, 'deg C')

    call log_flag('scale_area', model%glacier%scale_area, 'scale glacier area based on latitude')
    if (model%glacier%scale_area) then
       call log_param('length_scale_factor', model%glacier%length_scale_factor, 'scales dew and dns')
    elseif (length_scale_factor_in /= 1.0d0) then
       call log_param('length_scale_factor', model%glacier%length_scale_factor, &
            'reset from ' // real_to_string(length_scale_factor_in) // '; used only with scale_area = .true.')
    endif

    call log_param('diagnostic_minthck', model%glacier%diagnostic_minthck, 'm')

    call log_flag('redistribute_advanced_ice', model%glacier%redistribute_advanced_ice, &
         'redistribute advanced ice in the accumulation zone')
    if (model%glacier%redistribute_advanced_ice) then
       call log_param('thinning_rate_advanced_ice', model%glacier%thinning_rate_advanced_ice, 'm/yr')
    endif

    if (model%glacier%set_mu_star == GLACIER_MU_STAR_INVERSION) then
       call log_param('smb_weight_advanced_ice', model%glacier%smb_weight_advanced_ice, &
            'SMB weight in the advanced ablation zone')
       call log_param('mu_star_const', model%glacier%mu_star_const, 'mm/yr w.e./deg C')
       call log_param('mu_star_min', model%glacier%mu_star_min, 'mm/yr w.e./deg C')
       call log_param('mu_star_max', model%glacier%mu_star_max, 'mm/yr w.e./deg C')
       if (model%glacier%set_alpha_snow == GLACIER_ALPHA_SNOW_INVERSION) then
          call log_param('alpha_snow_const', model%glacier%alpha_snow_const)
          call log_param('alpha_snow_min', model%glacier%alpha_snow_min)
          call log_param('alpha_snow_max', model%glacier%alpha_snow_max)
          call log_param('beta_artm_max', model%glacier%beta_artm_max, 'deg C')
          call log_param('beta_artm_increment', model%glacier%beta_artm_increment, 'deg C')
          call log_param('baseline_date', model%glacier%baseline_date, 'yr; glaciers assumed in balance')
          call log_param('rgi_date', model%glacier%rgi_date, 'yr; date of RGI observations')
          call log_param('recent_date', model%glacier%recent_date, 'yr; date of recent SMB observations')
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
    ! Internal variables
    !------------------------------------------------------------------------------------
    type(glissade_options_type) :: options  !> Derived type holding all model options

    ! Copy model%options to options to save typing below
    ! Note: Originally, only model%options was passed in, but passing in the full model derived type
    !       allows the restart logic to be based on parameter values also.

    options = model%options

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

    ! add the SMB variable, based on options%smb_input
    ! Note: SMB can be input in one of several functional forms:
    !       - SMB(x,y)
    !       - SMB(x,y) at reference elevation, plus vertical gradient dSMB/dz(x,y)
    !       - SMB(x,y,z), where z is defined by reference levels
    ! Note: If the SMB field is 'acab', it is assumed to have units of m/y ice
    !       If the SMB field is 'smb', it is assumed to have units of mm/y w.e.

    select case(options%smb_input_function)

       case(SMB_INPUT_FUNCTION_XY)

          select case (options%smb_input)
          case (SMB_INPUT_MYR_ICE)
             call glissade_add_to_restart_variable_list('acab', model_id)
          case (SMB_INPUT_MMYR_WE)
             call glissade_add_to_restart_variable_list('smb', model_id)
          end select  ! smb_input

       case(SMB_INPUT_FUNCTION_XY_GRADZ)
          call glissade_add_to_restart_variable_list('smb_ref', model_id)
          call glissade_add_to_restart_variable_list('smb_gradz', model_id)
          call glissade_add_to_restart_variable_list('usrf_ref', model_id)

       case(SMB_INPUT_FUNCTION_XYZ)
          call glissade_add_to_restart_variable_list('smb_3d', model_id)

       ! Note: For SMB_INPUT_FUNCTION_PDD, precip and artm are read in as forcing
       !       and are not needed for exact restart

    end select  ! smb_input_function

    ! Similarly for surface temperature (artm), based on options%artm_input
    ! Note: These options share usrf_ref and zatm with the SMB options above.

    select case(options%artm_input_function)

       case(ARTM_INPUT_FUNCTION_XY_GRADZ)
          call glissade_add_to_restart_variable_list('artm_ref', model_id)
          call glissade_add_to_restart_variable_list('artm_gradz', model_id)
          if (options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
             ! usrf_ref was added to restart above; nothing to do here
          else
             call glissade_add_to_restart_variable_list('usrf_ref', model_id)
          endif

       case(ARTM_INPUT_FUNCTION_XYZ)
          call glissade_add_to_restart_variable_list('artm_3d', model_id)
          if (options%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then
             ! smb_levels was added to restart above; nothing to do here
          else
             call glissade_add_to_restart_variable_list('smb_levels', model_id)
          endif

       case(ARTM_INPUT_FUNCTION_XY_LAPSE)
          call glissade_add_to_restart_variable_list('artm_ref', model_id)
          ! Note: Instead of artm_gradz, there is a uniform lapse rate
          if (options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
             ! usrf_ref was added to restart above; nothing to do here
          else
             call glissade_add_to_restart_variable_list('usrf_ref', model_id)
          endif

    end select  ! artm_input_function

    ! Add anomaly forcing variables
    !TODO - Remove these? Anomaly forcing is typically in a forcing file, not the main input file.

    if (options%enable_smb_anomaly) then
          call glissade_add_to_restart_variable_list('smb_anomaly', model_id)
    endif

    if (options%enable_artm_anomaly) then
       call glissade_add_to_restart_variable_list('artm_anomaly', model_id)
    endif

    if (options%enable_bmlt_anomaly) then
       call glissade_add_to_restart_variable_list('bmlt_float_anomaly', model_id)
    endif

    if (options%read_lat_lon) then
       ! If lat and lon are to be read from the input file, they should be written to the restart file
       call glissade_add_to_restart_variable_list('lat lon', model_id)
    endif

    select case (options%whichbmlt_float)

       ! If bmlt_float is read from an external file at startup, then it needs to be in the restart file
       case (BMLT_FLOAT_EXTERNAL)
          call glissade_add_to_restart_variable_list('bmlt_float_external', model_id)

       ! If prescribing a warm ocean mask for depth-dependent melting, this needs to be read on restart
       case (BMLT_FLOAT_DEPTH)
          call glissade_add_to_restart_variable_list('warm_ocean_mask', model_id)

       case (BMLT_FLOAT_THERMAL_FORCING)
          ! Need the latest value of the thermal forcing field.
          ! This could be either the baseline value (if not updating during runtime), or a value read from a forcing file.
          ! If the latter, this field may not be needed, but include to be on the safe side, in case the forcing file
          !  is not read at restart.
          call glissade_add_to_restart_variable_list('thermal_forcing', model_id)

    end select  ! whichbmlt_float

    ! If using an ISMIP6 basin-based melt parameterization, and/or inverting for
    !  basin-scale quantities, we need a 2D field of basin numbers.
    ! Not strictly needed for the ISMIP6_LOCAL option, but included for diagnostics
    if (options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL .or.  &
        options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or. &
        options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE .or. &
        options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION_BASIN) then
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
     !                    To output these fields, the user must set restart_extend_velo = 1 in the config file.
     ! Note: It never hurts to write uvel/vvel_extend in place of uvel/vvel. But for most cases where restart
     !       is required (e.g., whole-ice-sheet simulations), velocities are zero along the boundaries
     !       and uvel/vvel are sufficient.

     if (options%restart_extend_velo == RESTART_EXTEND_VELO_TRUE) then
        call glissade_add_to_restart_variable_list('uvel_extend vvel_extend', model_id)
     else
        call glissade_add_to_restart_variable_list('uvel vvel', model_id)
     endif

     ! Glissade approximation options
     select case (options%which_ho_approx)

        case (HO_APPROX_DIVA)
           ! DIVA requires the 2D velocity, basal traction and effective viscosity for exact restart.
           ! Since the 2D velocity and basal traction are located on the staggered grid, these fields
           !  must be written to and read from the extended grid if velocites are nonzero at the boundaries.
           !
           ! Note: The 2D velocity is needed if the DIVA scheme solves for the mean velocity.
           !       If DIVA is configured to solve for the velocity at a specific level (e.g., the surface),
           !       then the 2D velocity could instead be copied from the 3D velocity array.
           ! Note: In addition to uvel/vvel_2D, DIVA requires the full 3D velocity field for exact restart,
           !       because horizontal transport is done before updating the velocity.
           
           if (options%restart_extend_velo == RESTART_EXTEND_VELO_TRUE) then
              call glissade_add_to_restart_variable_list('uvel_2d_extend vvel_2d_extend btractx_extend btracty_extend efvs', &
                   model_id)
           else
              call glissade_add_to_restart_variable_list('uvel_2d vvel_2d btractx btracty efvs', model_id)
           endif
           
        case default
           ! Other approximations (including SSA and L1L2) use the 3D uvel and vvel to initialize the velocity

     end select   ! which_ho_approx
        
     ! mask used to limit computation to active blocks
     if (options%compute_blocks == ACTIVE_BLOCKS_ONLY) then
        call glissade_add_to_restart_variable_list('ice_domain_mask', model_id)
     endif

     ! basal water option for Glissade
     select case (options%which_ho_bwat)
     case (HO_BWAT_NONE, HO_BWAT_CONSTANT)
        ! no restart variables needed
     case (HO_BWAT_FLUX_ROUTING)
        ! need to know the heat flux associated with refreezing meltwater
        call glissade_add_to_restart_variable_list('bhydroflx', model_id)
     case default
        ! restart needs to know bwat value
        call glissade_add_to_restart_variable_list('bwat', model_id)
     end select

     ! calving options for Glissade

     !TODO: CALVING_GRID_MASK and apply_calving_mask are redundant; remove one option?
     !      The advantage of apply_calving_mask is that it can be combined with other whichcalving options.
     if (options%whichcalving == CALVING_GRID_MASK .or. options%apply_calving_mask .or. &
          options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION) then
        if (options%which_ho_calving_front == HO_CALVING_FRONT_NO_SUBGRID) then
           call glissade_add_to_restart_variable_list('calving_mask', model_id)
        else   ! using a subgrid CF scheme
           call glissade_add_to_restart_variable_list('subgrid_calving_mask', model_id)
        endif
     endif

     ! The eigencalving calculation requires the product of eigenvalues of the horizontal strain rate tensor,
     !  which depends on the stress tensor, which is computed by the HO solver.
     ! On restart, the correct stress and strain rate tensors are not available, so we read in the eigenproduct.
     if (options%whichcalving == EIGEN_CALVING) then
        call glissade_add_to_restart_variable_list('eps_eigen1', model_id)
        call glissade_add_to_restart_variable_list('eps_eigen2', model_id)
     elseif (options%whichcalving == CALVING_STRESS .or. &
             options%whichcalving == CALVING_STRESS_STOCHASTIC .or. &
             options%whichcalving == CALVING_DAMAGE) then
        call glissade_add_to_restart_variable_list('tau_eigen1', model_id)
        call glissade_add_to_restart_variable_list('tau_eigen2', model_id)
     endif

     ! If forcing ice retreat, then we need ice_fraction_retreat_mask (which specifies the cells where retreat is forced)
     !  and reference_thck (which sets up an upper thickness limit for partly retreating cells)
     if (options%force_retreat /= FORCE_RETREAT_NONE) then
        call glissade_add_to_restart_variable_list('ice_fraction_retreat_mask', model_id)
        call glissade_add_to_restart_variable_list('reference_thck', model_id)
     endif

     ! lateral melt options
     if (options%which_lateral_melt == LATERAL_MELT_ISMIP) then
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
     if (options%overwrite_acab /= 0) then
        call glissade_add_to_restart_variable_list('overwrite_acab_mask', model_id)
     endif

    ! ==== Other options ====

    ! internal water option (for enthalpy scheme)
    select case (options%whichtemp)
      case (TEMP_ENTHALPY)
         ! restart needs to know internal water fraction
         call glissade_add_to_restart_variable_list('waterfrac', model_id)
      case default
         ! no restart variables needed
    end select

    ! basal sliding option
    select case (options%which_ho_babc)
      case (HO_BABC_COULOMB_FRICTION, HO_BABC_SCHOOF, HO_BABC_TSAI)
         ! Note: These options compute beta internally, so it does not need to be in the restart file.
         if (options%use_c_space_factor) then
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

    if (options%which_ho_powerlaw_c /= HO_POWERLAW_C_CONSTANT .and. &
        options%which_ho_powerlaw_c /= HO_POWERLAW_C_FUNCTION_COULOMB_C) then
       call glissade_add_to_restart_variable_list('powerlaw_c', model_id)
    endif

    if (options%which_ho_coulomb_c /= HO_COULOMB_C_CONSTANT) then
       call glissade_add_to_restart_variable_list('coulomb_c', model_id)
    endif

    ! If using the basin-scale inversion option for powerlaw_c or coulomb_c, we need a target thickness for grounded ice
    if (options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION_BASIN .or. &
        options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION_BASIN) then
       call glissade_add_to_restart_variable_list('grounded_thck_target', model_id)
    endif

    ! inversion options for ocean temperature corrections
    if (options%which_ho_deltaT_ocn /= HO_DELTAT_OCN_NONE) then
       call glissade_add_to_restart_variable_list('deltaT_ocn', model_id)
    endif

    ! If using the basin-scale inversion option for deltaT_ocn, we need a target thickness for floating ice
    if (options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION_BASIN) then
       call glissade_add_to_restart_variable_list('floating_thck_target', model_id)
    endif

    ! fields needed for inversion options that try to match local thickness, elevation, and/or surface speed
    ! Note: If usrf_obs is supplied, thck_obs will be computed at initialization.
    ! Note: Some options do not need usrf_obs and/or velo_sfc_obs as targets,
    !       but they can still be useful diagnostics.
    if (options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION .or. &
        options%which_ho_coulomb_c  == HO_COULOMB_C_INVERSION  .or. &
        options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION .or. &
        options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION_BASIN .or. &
        options%which_ho_coulomb_c  == HO_COULOMB_C_INVERSION_BASIN  .or. &
        options%which_ho_deltaT_ocn == HO_DELTAT_OCN_INVERSION_BASIN) then
       call glissade_add_to_restart_variable_list('usrf_obs', model_id)
       call glissade_add_to_restart_variable_list('velo_sfc_obs', model_id)
    endif

    ! fields needed for inversion options that try to match observed dthck_dt
    if (options%which_ho_deltaT_ocn == HO_DELTAT_OCN_DTHCK_DT .or. &
        options%enable_acab_dthck_dt_correction) then
       call glissade_add_to_restart_variable_list('dthk_dt_obs', model_id)
       call glissade_add_to_restart_variable_list('dthk_dt_obs_basin', model_id)
    endif

    select case (options%isostasy)
      case(ISOSTASY_COMPUTE)
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
         if (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING_LATVAR) then
            ! The relaxation timescale is a 2d field read at initialization and again on restart
            call glissade_add_to_restart_variable_list('tau_relax', model_id)
         endif
         case default
         ! no new restart variables needed
    end select

    select case (options%which_ho_ice_age)
       case(HO_ICE_AGE_COMPUTE)
          call glissade_add_to_restart_variable_list('ice_age', model_id)
       case default
          ! no restart variables needed
    end select

    if (model%options%enable_glaciers) then
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
