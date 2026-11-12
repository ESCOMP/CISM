!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_setup.F90 - part of the Community Ice Sheet Model (CISM)  
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

#include "glide_mask.inc"

module glide_setup

  ! general routines for initialisation, etc, called from top-level glimmer subroutines

  use glimmer_global, only: dp
  use glimmer_paramets, only : iulog

  implicit none

  private
  public :: glide_readconfig, glide_printconfig, glide_scale_params, &
            glide_load_sigma, glide_read_sigma, glide_calc_sigma

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine glide_readconfig(model,config)

    ! read Glide configuration file
    ! Note: sigma coordinates are handled by a subsequent call to glide_read_sigma
 
    use glide_types
    use glimmer_config
    use glimmer_log
    implicit none
    type(glide_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section

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

    ! read parameters
    !TODO - Put some sets of related parameters in their own section
    call GetSection(config,section,'parameters')
    if (associated(section)) then
       call handle_parameters(section, model)
    end if

    ! read geothermal heat flux
    ! NOTE: The [GTHF] section is ignored unless model%options%gthf = GTHF_COMPUTE
    if (model%options%gthf == GTHF_COMPUTE) then
       call GetSection(config,section,'GTHF')
       if (associated(section)) then
          call handle_gthf(section, model)
       end if
    endif

    ! read isostasy
    ! NOTE: The [isostasy] section is ignored unless model%options%isostasy = ISOSTASY_COMPUTE
    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call GetSection(config,section,'isostasy')
       if (associated(section)) then
          call handle_isostasy(section, model)
       end if
    endif

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file.
    ! (Glint restart variables are handled separately by Glint setup routines.)
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13

    call define_glide_restart_variables(model, model%model_id)

  end subroutine glide_readconfig

!-------------------------------------------------------------------------

  subroutine glide_printconfig(model)

    !> print model configuration to log
    use glimmer_log
    use glide_types
    implicit none
    type(glide_global_type)  :: model !> model instance

    call write_log_div
    call print_grid(model)
    call print_time(model)
    call print_options(model)
    call print_parameters(model)
    call print_gthf(model)
    call print_isostasy(model)

  end subroutine glide_printconfig

!-------------------------------------------------------------------------
    
  subroutine glide_scale_params(model)
    !> scale parameters
    use glide_types
    use glimmer_physcon,  only: scyr

    implicit none

    type(glide_global_type)  :: model !> model instance

    model%numerics%dttem = model%numerics%ntem * model%numerics%tinc   

    ! convert dt and dttem from yr to s
    model%numerics%dt     = model%numerics%tinc * scyr
    model%numerics%dttem  = model%numerics%dttem * scyr

    ! scale glide basal traction parameters
    model%paramets%btrac_const = model%paramets%btrac_const/scyr
    model%paramets%btrac_max   = model%paramets%btrac_max/scyr

  end subroutine glide_scale_params

!-------------------------------------------------------------------------

  subroutine glide_read_sigma(model,config)

    ! read sigma levels from configuration file, if present
    ! called immediately after glide_readconfig

    use glide_types
    use glimmer_config
    use glimmer_log
    implicit none

    type(glide_global_type) :: model        !> model instance
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

  end subroutine glide_read_sigma

!-------------------------------------------------------------------------

  subroutine glide_load_sigma(model,unit)

    ! Compute sigma coordinates or read them from a file
    ! Note: This subroutine is called from glide_initialise or glissade_initialise.
    !       If sigma levels are provided in the config file, then they are read
    !        in by glide_read_sigma, and model%options%which_sigma is set to
    !        SIGMA_CONFIG, in which case this subroutine does nothing.

    use glide_types
    use glimmer_log
    use glimmer_filenames
    use cism_parallel, only: main_task, broadcast

    implicit none

    ! Arguments
    type(glide_global_type),intent(inout) :: model !> Ice model to use
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
          model%numerics%sigma(up) = glide_calc_sigma(level, 2.d0)
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

       ! sigma levels have already been read from glide_read_sigma

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
             model%numerics%sigma(up) = glide_calc_sigma_pattyn(level)
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
    
  end subroutine glide_load_sigma

!--------------------------------------------------------------------------------

  function glide_calc_sigma(x,n)

     implicit none
     real(dp) :: glide_calc_sigma, x, n
      
     glide_calc_sigma = (1-(x+1)**(-n)) / (1-2**(-n))

  end function glide_calc_sigma

!--------------------------------------------------------------------------------

  function glide_calc_sigma_pattyn(x)

     ! Implements an alternate set of sigma levels that encourages better
     ! convergence for higher-order velocities

     implicit none
     real(dp) :: glide_calc_sigma_pattyn, x

     glide_calc_sigma_pattyn =   &
         (-2.5641025641d-4)*(41d0*x)**2+3.5256410256d-2*(41d0*x)-8.0047080075d-13

  end function glide_calc_sigma_pattyn

!-------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! grid sizes

  subroutine handle_grid(section, model)
    use glimmer_config
    use glide_types
    use glimmer_filenames
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'ewn',model%general%ewn)
    call GetValue(section,'nsn',model%general%nsn)
    call GetValue(section,'upn',model%general%upn)
    call GetValue(section,'dew',model%numerics%dew)
    call GetValue(section,'dns',model%numerics%dns)
    call GetValue(section,'sigma_file',model%funits%sigfile)

    ! We set this flag to one to indicate we've got a sigfile name.
    ! A warning/error is generated if sigma levels are specified in some other way
    ! and mangle the name
    if (trim(model%funits%sigfile) /= '') then
       model%funits%sigfile = filenames_inputname(model%funits%sigfile)
       model%options%which_sigma = SIGMA_EXTERNAL
    end if

  end subroutine handle_grid

!--------------------------------------------------------------------------------

  subroutine print_grid(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=512) :: message

    call write_log('Grid specification')
    call write_log('------------------')
    write(message,*) 'ewn             : ',model%general%ewn
    call write_log(trim(message))
    write(message,*) 'nsn             : ',model%general%nsn
    call write_log(trim(message))
    write(message,*) 'upn             : ',model%general%upn
    call write_log(trim(message))
    write(message,*) 'EW grid spacing : ',model%numerics%dew
    call write_log(trim(message))
    write(message,*) 'NS grid spacing : ',model%numerics%dns
    call write_log(trim(message))

    write(message,*) 'sigma file      : ',trim(model%funits%sigfile)
    call write_log(trim(message))
    call write_log('')

  end subroutine print_grid

!--------------------------------------------------------------------------------

  ! time
  subroutine handle_time(section, model)
    use glimmer_config
    use glide_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

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
    call GetValue(section,'ntem',model%numerics%ntem)
    call GetValue(section,'profile',model%numerics%profile_period)

    call GetValue(section,'idiag',model%numerics%idiag)
    call GetValue(section,'jdiag',model%numerics%jdiag)

    !Note: Either dt_diag or ndiag can be specified in the config file.
    !      If dt_diag is specified, it is used to compute ndiag. (Output is written every ndiag timesteps.)
    call GetValue(section,'dt_diag',model%numerics%dt_diag)
    call GetValue(section,'ndiag',model%numerics%ndiag)

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
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log('Time steps')
    call write_log('----------')
    write(message,*) 'start time (yr)     : ',model%numerics%tstart
    call write_log(message)
    write(message,*) 'end time (yr)       : ',model%numerics%tend
    call write_log(message)
    write(message,*) 'time step (yr)      : ',model%numerics%tinc
    call write_log(message)
    if (model%numerics%nsteps_per_year > 0) then
       write(message,*) 'nsteps per year     : ',model%numerics%nsteps_per_year
       call write_log(message)
    endif
    write(message,*) 'thermal dt factor   : ',model%numerics%ntem
    call write_log(message)
    if ( (model%numerics%ntem < 1.0d0) .or. & 
       (floor(model%numerics%ntem) /= model%numerics%ntem) ) then
       call write_log('ntem is a multiplier on the basic time step.  It should be a positive integer.  Aborting.',GM_FATAL)
    endif

    write(message,*) 'profile frequency   : ',model%numerics%profile_period
    call write_log(message)

    if (model%numerics%dt_diag > 0.d0) then
       write(message,*) 'diagnostic interval (years):',model%numerics%dt_diag
       call write_log(message)
    elseif (model%numerics%ndiag > 0) then
       write(message,*) 'diagnostic interval (steps):',model%numerics%ndiag
       call write_log(message)
    endif

     !WHL - Written to log in glide_init_diag
!    write(message,*) 'idiag               : ',model%numerics%idiag
!    call write_log(message)
!    write(message,*) 'jdiag               : ',model%numerics%jdiag
!    call write_log(message)

    call write_log('')

  end subroutine print_time

!--------------------------------------------------------------------------------

  ! options
  subroutine handle_options(section, model)

    use glimmer_config
    use glide_types

    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'dycore',model%options%whichdycore)
    call GetValue(section,'evolution',model%options%whichevol)
    call GetValue(section,'temperature',model%options%whichtemp)
    call GetValue(section,'temp_init',model%options%temp_init)
    call GetValue(section,'flow_law',model%options%whichflwa)
    call GetValue(section,'slip_coeff',model%options%whichbtrc)
    call GetValue(section,'basal_water',model%options%whichbwat)
    call GetValue(section,'basal_mass_balance',model%options%basal_mbal)
    call GetValue(section,'smb_input',model%options%smb_input)
    call GetValue(section,'gthf',model%options%gthf)
    call GetValue(section,'isostasy',model%options%isostasy)
    call GetValue(section,'marine_margin',model%options%whichcalving)
    call GetValue(section,'calving_init',model%options%calving_init)
    call GetValue(section,'dm_dt_diag',model%options%dm_dt_diag)
    call GetValue(section,'diag_minthck',model%options%diag_minthck)
    call GetValue(section,'vertical_integration',model%options%whichwvel)
    call GetValue(section,'periodic_ew',model%options%periodic_ew)
    call GetValue(section,'sigma',model%options%which_sigma)
    call GetValue(section,'ioparams',model%funits%ncfile)
    call GetValue(section,'forcewrite_final', model%options%forcewrite_final)
    call GetValue(section,'restart',model%options%is_restart)

  end subroutine handle_options

!--------------------------------------------------------------------------------

  subroutine print_options(model)

    use glide_types
    use glimmer_log
    use glimmer_physcon, only: scyr
    use cism_parallel, only: tasks

    implicit none

    type(glide_global_type)  :: model
    character(len=500) :: message

    ! basic Glide/Glimmer options

    ! Glide is the only dycore option here (dycore = 0)
    character(len=*), parameter :: dycore = 'glide              '  ! Glimmer SIA

    character(len=*), dimension(0:2), parameter :: evolution = (/ &
         'pseudo-diffusion         ', &
         'ADI scheme               ', &
         'iterated diffusion       ' /)

    character(len=*), dimension(0:2), parameter :: temperature = (/ &
         'isothermal            ', &
         'prognostic temperature', &
         'constant in time      ' /)

    character(len=*), dimension(0:2), parameter :: temp_init = (/ &
         'set to 0 C                  ', &
         'set to surface air temp     ', &
         'linear vertical profile     ' /)

    character(len=*), dimension(0:2), parameter :: flow_law = (/ &
         'uniform factor flwa         ', &
         'Paterson and Budd (T = -5 C)', &
         'Paterson and Budd           ' /)

    !TODO - Rename slip_coeff to which_btrc?
    character(len=*), dimension(0:5), parameter :: slip_coeff = (/ &
         'no basal sliding       ', &
         'constant basal traction', &
         'constant where bwat > 0', &
         'constant where T = Tpmp', &
         'linear function of bmlt', &
         'tanh function of bwat  ' /)

    character(len=*), dimension(0:3), parameter :: basal_water = (/ &
         'none                     ', &
         'local water balance      ', &
         'local + steady-state flux', &
         'Constant value (= 10 m)  ' /)

    character(len=*), dimension(0:1), parameter :: b_mbal = (/ &
         'not in continuity eqn', &
         'in continuity eqn    ' /)

    character(len=*), dimension(0:1), parameter :: smb_input = (/ &
         'SMB input in units of m/yr ice  ', &
         'SMB input in units of mm/yr w.e.' /)

    ! NOTE: Set gthf = 1 in the config file to read the geothermal heat flux from an input file.
    !       Otherwise it will be overwritten, even if the 'bheatflx' field is present.

    character(len=*), dimension(0:2), parameter :: gthf = (/ &
         'uniform geothermal flux         ', &
         'read flux from file, if present ', &
         'compute flux from diffusion eqn ' /)

    character(len=*), dimension(0:1), parameter :: isostasy = (/ &
         'no isostasy calculation         ', &
         'compute isostasy with model     ' /)

    character(len=*), dimension(0:5), parameter :: marine_margin = (/ &
         'no calving law                    ', &
         'remove all floating ice           ', &
         'remove fraction of floating ice   ', &
         'relaxed bedrock threshold         ', &
         'present bedrock threshold         ', &
         'Huybrechts calving                '/)

    character(len=*), dimension(0:1), parameter :: init_calving = (/ &
         'no calving at initialization    ', &
         'ice calves at initialization    ' /)

    character(len=*), dimension(0:1), parameter :: dm_dt_diag = (/ &
         'write dmass/dt diagnostic in units of kg/s ',  &
         'write dmass/dt diagnostic in units of Gt/yr'/)

    character(len=*), dimension(0:1), parameter :: diag_minthck = (/ &
         'include cells with H > 0 in global diagnostics     ', &
         'include cells with H > thklim in global diagnostics'/)

    character(len=*), dimension(0:1), parameter :: vertical_integration = (/ &
         'standard     ', &
         'obey upper BC' /)

    ! higher-order options

    ! basal water options

    call write_log('Dycore options')
    call write_log('-------------')

    write(message,*) 'I/O parameter file      : ',trim(model%funits%ncfile)
    call write_log(message)

    if (model%options%whichdycore /= DYCORE_GLIDE) then
       call write_log('Error, the Glide dycore requires dycore = 0',GM_FATAL)
    end if
    write(message,*) 'Dycore                  : ',model%options%whichdycore,dycore
    call write_log(message)

    if (tasks > 1) then
       call write_log('Error, Glide dycore not supported for runs with more than one processor', GM_FATAL)
    end if

    if (model%options%whichevol == EVOL_ADI) then
       call write_log('WARNING: exact restarts are not currently possible with ADI evolution', GM_WARNING)
    endif

    if (model%options%whichtemp < 0 .or. model%options%whichtemp >= size(temperature)) then
       call write_log('Error, temperature option out of range',GM_FATAL)
    end if
    write(message,*) 'temperature calculation : ',model%options%whichtemp,temperature(model%options%whichtemp)
    call write_log(message)

    if (model%options%temp_init < 0 .or. model%options%temp_init >= size(temp_init)) then
       call write_log('Error, temp_init option out of range',GM_FATAL)
    end if
    ! Note: If reading temperature from an input or restart file, the temp_init option is overridden,
    !        in which case it could be confusing here to write the option to the log file.
    !       The method actually used is written to the log file by glide_init_temp. 

    if (model%options%whichflwa < 0 .or. model%options%whichflwa >= size(flow_law)) then
       call write_log('Error, flow_law out of range',GM_FATAL)
    end if
    write(message,*) 'flow law                : ',model%options%whichflwa,flow_law(model%options%whichflwa)
    call write_log(message)

    if (model%options%whichbwat < 0 .or. model%options%whichbwat >= size(basal_water)) then
       call write_log('Error, basal_water out of range',GM_FATAL)
    end if
    write(message,*) 'basal_water             : ',model%options%whichbwat,basal_water(model%options%whichbwat)
    call write_log(message)

    ! unsupported basal_water options
    if (model%options%whichbwat == BWATER_FLUX) then
      call write_log('WARNING: Steady state routing basal_water option is not scientifically supported.  &
           &USE AT YOUR OWN RISK.', GM_WARNING)
    endif

    if (model%options%whichcalving < 0 .or. model%options%whichcalving >= size(marine_margin)) then
       call write_log('Error, marine_margin out of range',GM_FATAL)
    end if
    write(message,*) 'marine_margin           : ', model%options%whichcalving, marine_margin(model%options%whichcalving)
    call write_log(message)

    if (model%options%calving_init < 0 .or. model%options%calving_init >= size(init_calving)) then
       call write_log('Error, calving_init out of range',GM_FATAL)
    end if
    write(message,*) 'calving_init            : ', model%options%calving_init, init_calving(model%options%calving_init)
    call write_log(message)

    if (model%options%whichbtrc < 0 .or. model%options%whichbtrc >= size(slip_coeff)) then
       call write_log('Error, slip_coeff out of range',GM_FATAL)
    end if

    write(message,*) 'slip_coeff              : ', model%options%whichbtrc, slip_coeff(model%options%whichbtrc)
    call write_log(message)

    if (model%options%whichevol < 0 .or. model%options%whichevol >= size(evolution)) then
       call write_log('Error, evolution out of range',GM_FATAL)
    end if

    write(message,*) 'evolution               : ', model%options%whichevol, evolution(model%options%whichevol)
    call write_log(message)

    if (model%options%dm_dt_diag < 0 .or. model%options%dm_dt_diag >= size(dm_dt_diag)) then
       call write_log('Error, dm_dt_diag out of range',GM_FATAL)
    end if

    if (model%options%diag_minthck < 0 .or. model%options%diag_minthck >= size(diag_minthck)) then
       call write_log('Error, diag_minthck out of range',GM_FATAL)
    end if

    write(message,*) 'minthck for diagnostics : ',model%options%diag_minthck, &
         diag_minthck(model%options%diag_minthck)
    call write_log(message)

    if (model%options%whichwvel < 0 .or. model%options%whichwvel >= size(vertical_integration)) then
       call write_log('Error, vertical_integration out of range',GM_FATAL)
    end if

    write(message,*) 'vertical_integration    : ',model%options%whichwvel,vertical_integration(model%options%whichwvel)
    call write_log(message)

    if (model%options%basal_mbal < 0 .or. model%options%basal_mbal >= size(b_mbal)) then
       call write_log('Error, basal_mass_balance out of range',GM_FATAL)
    end if

    write(message,*) 'basal mass balance      : ',model%options%basal_mbal,b_mbal(model%options%basal_mbal)
    call write_log(message)

    if (model%options%smb_input < 0 .or. model%options%smb_input >= size(smb_input)) then
       call write_log('Error, smb_input option out of range',GM_FATAL)
    end if

    write(message,*) 'smb input units         : ',model%options%smb_input,smb_input(model%options%smb_input)
    call write_log(message)

    if (model%options%smb_input == SMB_INPUT_MYR_ICE) then
    endif

    if (model%options%gthf < 0 .or. model%options%gthf >= size(gthf)) then
       call write_log('Error, geothermal flux option out of range',GM_FATAL)
    end if

    write(message,*) 'geothermal heat flux    : ',model%options%gthf,gthf(model%options%gthf)
    call write_log(message)

    if (model%options%isostasy < 0 .or. model%options%isostasy >= size(isostasy)) then
       call write_log('Error, isostasy option out of range',GM_FATAL)
    end if

    write(message,*) 'isostasy                : ',model%options%isostasy,isostasy(model%options%isostasy)
    call write_log(message)

    if (model%options%periodic_ew) then
       if (model%options%whichevol == EVOL_ADI) then
          call write_log('Periodic boundary conditions not implemented in ADI scheme',GM_FATAL)
       end if
       call write_log('Periodic EW lateral boundary condition')
       call write_log('  Slightly cheated with how temperature is implemented.',GM_WARNING)
    end if

    if (model%options%forcewrite_final) then
       call write_log('Force write to output files when the run completes')
    endif

    if (model%options%is_restart == STANDARD_RESTART) then
       call write_log('Restarting model from a previous run')
    elseif (model%options%is_restart == HYBRID_RESTART) then
       call write_log('Hybrid restart from a previous run')
    end if

  end subroutine print_options

!--------------------------------------------------------------------------------

  ! parameters
  subroutine handle_parameters(section, model)

    use glimmer_config
    use glide_types
    use glimmer_log
    use glimmer_physcon, only: rhoi, rhoo, grav, shci, lhci, trpt, n_glen

    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model
    real(dp), pointer, dimension(:) :: tempvar => NULL()
    integer :: loglevel

    !Note: The following physical constants have default values in glimmer_physcon.F90.
    !      Some test cases might specify different values. The default values
    !       can therefore be overridden by the user in the config file.
    !      Glide is not run in CESM, but this file is compiled there. With CCSMCOUPLED,
    !       these constants are shr_const_mod parameters and cannot be passed to GetValue.
#ifndef CCSMCOUPLED
    call GetValue(section,'rhoi', rhoi)
    call GetValue(section,'rhoo', rhoo)
    call GetValue(section,'grav', grav)
    call GetValue(section,'shci', shci)
    call GetValue(section,'lhci', lhci)
    call GetValue(section,'trpt', trpt)
#endif
    call GetValue(section,'n_glen', n_glen)

    loglevel = GM_levels-GM_ERROR
    call GetValue(section,'log_level',loglevel)
    call glimmer_set_msg_level(loglevel)

    !TODO - Change 'ice_limit' to 'thklim'?
    call GetValue(section,'ice_limit',          model%numerics%thklim)
    call GetValue(section,'geothermal',         model%paramets%geot)
    call GetValue(section,'flow_factor_ground', model%paramets%flow_enhancement_factor_ground)
    call GetValue(section,'default_flwa',       model%paramets%default_flwa)
    call GetValue(section,'hydro_time',         model%paramets%hydtim)

    ! SMB and BMB parameters

    ! calving parameters
    call GetValue(section,'marine_limit',       model%calving%marine_limit)
    call GetValue(section,'calving_fraction',   model%calving%calving_fraction)

    ! NOTE: bpar is used only for BTRC_TANH_BWAT
    !       btrac_max and btrac_slope are used (with btrac_const) for BTRC_LINEAR_BMLT
    !       btrac_const is used for several options
    call GetValue(section,'basal_tract_const', model%paramets%btrac_const)
    call GetValue(section,'basal_tract_max',   model%paramets%btrac_max)
    call GetValue(section,'basal_tract_slope', model%paramets%btrac_slope)
    call GetValue(section,'basal_tract_tanh',  tempvar, 5)
    if (associated(tempvar)) then
       model%paramets%bpar(:) = tempvar(:)
       deallocate(tempvar)
    end if

  end subroutine handle_parameters

!--------------------------------------------------------------------------------

  subroutine print_parameters(model)

    use glimmer_physcon, only: rhoi, rhoo, lhci, shci, trpt, grav, n_glen
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log(' ')
    call write_log('Parameters')
    call write_log('----------')

    write(message,*) 'thickness limit for dynamically active ice (m) : ', model%numerics%thklim
    call write_log(message)

    ! calving options

    if (model%options%whichcalving == CALVING_FLOAT_FRACTION) then
       write(message,*) 'ice fraction lost in calving  : ', model%calving%calving_fraction
       call write_log(message)
    end if

    if (model%options%whichcalving == CALVING_RELX_THRESHOLD .or.  &
        model%options%whichcalving == CALVING_TOPG_THRESHOLD) then
       write(message,*) 'marine depth limit (m)        : ', model%calving%marine_limit
       call write_log(message)
    endif

    write(message,*) 'ice density (kg/m^3)          : ', rhoi
    call write_log(message)

    write(message,*) 'ocean density (kg/m^3)        : ', rhoo
    call write_log(message)

    write(message,*) 'gravitational accel (m/s^2)   : ', grav
    call write_log(message)

    write(message,*) 'heat capacity of ice (J/kg/K) : ', shci
    call write_log(message)

    write(message,*) 'latent heat of ice (J/kg)     : ', lhci
    call write_log(message)

    write(message,*) 'triple point of water (K)     : ', trpt
    call write_log(message)

    write(message,*) 'Glen flow law exponent        : ', n_glen
    call write_log(message)

    write(message,*) 'geothermal flux  (W/m^2)      : ', model%paramets%geot
    call write_log(message)

    write(message,*) 'flow factor (grounded ice)    : ', model%paramets%flow_enhancement_factor_ground
    call write_log(message)

    write(message,*) 'basal hydro time constant (yr): ', model%paramets%hydtim
    call write_log(message)

    if (model%options%whichflwa == FLWA_CONST_FLWA) then
       write(message,*) 'constant flow factor (Pa^-n yr^-1) :', model%paramets%default_flwa
       call write_log(message)
    end if

    if (model%options%whichbtrc == BTRC_CONSTANT      .or.  &
        model%options%whichbtrc == BTRC_CONSTANT_BWAT .or.  &
        model%options%whichbtrc == BTRC_LINEAR_BMLT   .or.  &
        model%options%whichbtrc == BTRC_CONSTANT_BPMP) then
       write(message,*) 'basal traction param (m/yr/Pa)      : ', model%paramets%btrac_const
       call write_log(message)
    end if

    if (model%options%whichbtrc == BTRC_TANH_BWAT) then
       write(message,*) 'basal traction tanh factors: ',model%paramets%bpar(1)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(2)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(3)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(4)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(5)
       call write_log(message)
    end if

    if (model%options%whichbtrc == BTRC_LINEAR_BMLT) then
       write(message,*) 'basal traction max            : ',model%paramets%btrac_max
       call write_log(message)
       write(message,*) 'basal traction slope          : ',model%paramets%btrac_slope
       call write_log(message)
    end if

    if (model%numerics%idiag < 1 .or. model%numerics%idiag > model%general%ewn     &
                                        .or.                                                     &
        model%numerics%jdiag < 1 .or. model%numerics%jdiag > model%general%nsn) then
        call write_log('Error, global diagnostic point (idiag, jdiag) is out of bounds', GM_FATAL)
    endif

  end subroutine print_parameters

!--------------------------------------------------------------------------------

  ! Sigma levels
  subroutine handle_sigma(section, model)

    use glimmer_config
    use glide_types
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    if (model%options%which_sigma==SIGMA_EXTERNAL) then
       call write_log('Sigma levels specified twice - use only'// &
            ' config file or separate file, not both',GM_FATAL)
    else
       call GetValue(section,'sigma_levels',model%numerics%sigma,model%general%upn)
    end if

  end subroutine handle_sigma

!--------------------------------------------------------------------------------

  subroutine print_sigma(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message,temp
    integer :: i

    call write_log('Sigma levels:')
    call write_log('------------------')
    message=''
    do i=1,model%general%upn
       write(temp,'(f6.3)') model%numerics%sigma(i)
       message=trim(message)//trim(temp)
    enddo
    call write_log(trim(message))
    call write_log('')
    
  end subroutine print_sigma

!--------------------------------------------------------------------------------

  ! geothermal heat flux calculations
  subroutine handle_gthf(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'num_dim',model%lithot%num_dim)
    call GetValue(section,'nlayer',model%lithot%nlayer)
    call GetValue(section,'surft',model%lithot%surft)
    call GetValue(section,'rock_base',model%lithot%rock_base)
    call GetValue(section,'numt',model%lithot%numt)
    call GetValue(section,'rho',model%lithot%rho_r)
    call GetValue(section,'shc',model%lithot%shc_r)
    call GetValue(section,'con',model%lithot%con_r)
  end subroutine handle_gthf

!--------------------------------------------------------------------------------

  subroutine print_gthf(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message
    
    if (model%options%gthf == GTHF_COMPUTE) then
       call write_log('Geothermal heat flux configuration')
       call write_log('----------------------------------')
       if (model%lithot%num_dim==1) then
          call write_log('solve 1D diffusion equation')
       else if (model%lithot%num_dim==3) then          
          call write_log('solve 3D diffusion equation')
       else
          call write_log('Wrong number of dimensions.',GM_FATAL,__FILE__,__LINE__)
       end if
       write(message,*) 'number of layers                     : ',model%lithot%nlayer
       call write_log(message)
       write(message,*) 'initial surface temperature          : ',model%lithot%surft
       call write_log(message)
       write(message,*) 'rock base                            : ',model%lithot%rock_base
       call write_log(message)
       write(message,*) 'density of rock layer                : ',model%lithot%rho_r
       call write_log(message)
       write(message,*) 'specific heat capacity of rock layer : ',model%lithot%shc_r
       call write_log(message)
       write(message,*) 'thermal conductivity of rock layer   : ',model%lithot%con_r
       call write_log(message)
       write(message,*) 'number of time steps for spin-up     : ',model%lithot%numt
       call write_log(message)
       call write_log('')
    end if
  end subroutine print_gthf

!--------------------------------------------------------------------------------

  subroutine handle_isostasy(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'lithosphere',model%isostasy%lithosphere)
    call GetValue(section,'asthenosphere',model%isostasy%asthenosphere)
    call GetValue(section,'which_relaxed',model%isostasy%which_relaxed)
    call GetValue(section,'tau_relax_const',model%isostasy%tau_relax_const)
    call GetValue(section,'load_update_interval',model%isostasy%load_update_interval)
    call GetValue(section,'flexural_rigidity',model%isostasy%rbel%d)

  end subroutine handle_isostasy

!--------------------------------------------------------------------------------

  subroutine print_isostasy(model)

    use glide_types
    use glimmer_log

    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    character(len=*), dimension(0:1), parameter :: lithosphere = (/ &
         'local lithosphere         ', &
         'elastic lithosphere       ' /)

    character(len=*), dimension(0:1), parameter :: asthenosphere = (/ &
         'fluid asthenosphere                ', &
         'relaxing with constant timescale   ' /)

    character(len=*), dimension(0:2), parameter :: which_relaxed = (/ &
         'read both topg and relx from input file', &
         'set relx to input topg                 ', &
         'compute relx = input topg + load       ' /)
    
    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call write_log(' ')
       call write_log('Isostasy')
       call write_log('--------')

       if (model%isostasy%lithosphere < 0 .or. model%isostasy%lithosphere >= size(lithosphere)) then
          call write_log('Error, lithosphere option out of range', GM_FATAL)
       else
          write(message,*) 'lithosphere                 : ',model%isostasy%lithosphere,  &
               lithosphere(model%isostasy%lithosphere)
          call write_log(message)
       endif

       if (model%isostasy%lithosphere==LITHOSPHERE_ELASTIC) then
          write(message,*) 'flexural rigidity (N m)     : ', model%isostasy%rbel%d
          call write_log(message)
          write(message,*) 'load update interval (yr)   : ', model%isostasy%load_update_interval
          call write_log(message)
       end if

       if (model%isostasy%asthenosphere < 0 .or. model%isostasy%asthenosphere >= size(asthenosphere)) then
          call write_log('Error, asthenosphere option out of range', GM_FATAL)
       else
          write(message,*) 'asthenosphere               : ',model%isostasy%asthenosphere,  &
               asthenosphere(model%isostasy%asthenosphere)
          call write_log(message)
       endif

       if (model%isostasy%which_relaxed < 0 .or. model%isostasy%which_relaxed >= size(which_relaxed)) then
          call write_log('Error, which_relaxed option out of range', GM_FATAL)
       else
          write(message,*) 'which_relaxed               : ',model%isostasy%which_relaxed,  &
               which_relaxed(model%isostasy%which_relaxed)
          call write_log(message)
       endif

       if (model%isostasy%asthenosphere==ASTHENOSPHERE_RELAXING) then
          write(message,*) 'relaxation constant (yr)    : ', model%isostasy%tau_relax_const
          call write_log(message)
       endif

    endif   ! compute isostasy

  end subroutine print_isostasy

!--------------------------------------------------------------------------------

  subroutine define_glide_restart_variables(model, model_id)

    !> This subroutine analyzes the glide/glissade options input by the user in the config file
    !> and determines which variables are necessary for an exact restart.  MJH 1/11/2013

    ! Please comment thoroughly the reasons why a particular variable needs to be a restart variable for a given config.
    ! Note: This subroutine assumes that any restart variables you add are loadable.
    !       Check glide_vars.def to make sure any added variables have 'load: 1'

    use glide_types
    use glide_io, only: glide_add_to_restart_variable_list

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glide_global_type), intent (in) :: model  !> Derived type holding all model info

    integer, intent(in) :: model_id  !> identifier of this ice sheet instance (1 - N, where N is the total number of ice sheet models in this run)

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------
    type(glide_options) :: options  !> Derived type holding all model options

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
    ! Note: the conversion from temp/flwa to tempstag/flwastag (if necessary) happens in glide_io.F90
    ! bheatflx, artm, acab - boundary conditions.  Of course if these fields are 0 they don't need 
    !        to be in the restart file, but without adding a check for that we cannot assume any of them are.
    !        There are some options where artm would not be needed.  Logic could be added to make that distinction.
    !        Note that bheatflx may not be an input variable but can also be assigned as a parameter in the config file!
    call glide_add_to_restart_variable_list('topg thk temp bheatflx artm', model_id)

    ! add the SMB variable, based on options%smb_input
    ! Note: If the SMB field is 'acab', it is assumed to have units of m/y ice
    !       If the SMB field is 'smb', it is assumed to have units of mm/y w.e.
    select case (options%smb_input)
    case (SMB_INPUT_MYR_ICE)
       call glide_add_to_restart_variable_list('acab', model_id)
    case (SMB_INPUT_MMYR_WE)
       call glide_add_to_restart_variable_list('smb', model_id)
    end select  ! smb_input

    ! add Glide velocity and thermal restart variables
    ! wgrd & wvel - temp driver calculates weff = f(wgrd, wvel) so both are needed by temp code.
    !               It looks possible to calculate wvel on a restart from wgrd because wvel does not 
    !               appear to require a time derivative (see subroutine wvelintg).  
    !               wgrd does require time derivatives and therefore should be
    !               calculated at the end of each time step and stored as a restart variable
    !               so that the time derivatives do not need to be restart variables.
    !               For now I am calculating wvel at the same time (end of glide time step) 
    !               and then saving both as restart variables.  This has the advantage of
    !               them being on consistent time levels in the output file.  
    !               (If we waited to calculate wvel in the temp driver, we would not need to
    !                add it as a restart variable, been then in the output wgrd and wvel would
    !                be based on different time levels.)
    ! flwa - in principal this could be reconstructed from temp.  However in the current 
    !        implementation of glide the flwa calculation occurs after temp evolution but 
    !        before thk evolution.  This means flwa is calculated from the current temp and 
    !        the old thk.  The old thk is not available on a restart (just the current thk).
    !        (thk is needed to calculate flwa for 1) a mask for where ice is, 2) correction for pmp.)
    call glide_add_to_restart_variable_list('thkmask wgrd wvel flwa uvel vvel', model_id)

    ! basal water option for Glide
    select case (options%whichbwat)
    case (BWATER_NONE, BWATER_CONST)
       ! no restart variables needed
    case default
       ! restart needs to know bwat value
       call glide_add_to_restart_variable_list('bwat', model_id)
    end select

    ! geothermal heat flux option
    select case (options%gthf)
      case(GTHF_COMPUTE)
         ! restart needs to know lithosphere temperature
         call glide_add_to_restart_variable_list('litho_temp', model_id)
      case default
         ! no restart variables needed
    end select

    select case (options%isostasy)
      case(ISOSTASY_COMPUTE)
         ! restart needs to know relaxed topography (the topography to which the mantle would relax with no load)
         ! MJH: I suspect that relx is only needed when asthenosphere=1 (relaxing mantle), but I'm not sure -
         !      this should be tested when isostasy implementation is finalized/tested.
         ! WHL: Looking at subroutine isos_compute, I think relx is needed also when asthenosphere = 0 (fluid mantle).
         !      In this case, topg is set instantaneously to relx - load.
         call glide_add_to_restart_variable_list('relx', model_id)
         !WHL - The load field also needs to be in the restart file. The reason is that the load is updated
         !      at a period set by isostasy%period. If we restart between two updates, we need to use the most
         !      recently computed load. If we recompute the load right after restarting, the restart may not be exact.
         call glide_add_to_restart_variable_list('load', model_id)
      case default
         ! no new restart variables needed
    end select

  end subroutine define_glide_restart_variables

!--------------------------------------------------------------------------------

end module glide_setup

!--------------------------------------------------------------------------------
