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

  implicit none

  private
  public :: glide_readconfig, glide_printconfig, glide_scale_params, &
            glide_load_sigma, glide_read_sigma, glide_calc_sigma, glide_get_zocn

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine glide_readconfig(model,config)

    ! read Glide configuration file
    ! Note: sigma coordinates are handled by a subsequent call to glide_read_sigma
 
    use glide_types
    use glimmer_config
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

    ! read options for higher-order computation
    call GetSection(config,section,'ho_options')
    if (associated(section)) then
        call handle_ho_options(section, model)
    end if

    ! read options for computation using an external dycore -- Doug Ranken 04/20/12
    call GetSection(config,section,'external_dycore_options')
    if (associated(section)) then
        call handle_dycore_options(section, model)
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

    ! read glacier info
    if (model%options%enable_glaciers) then
       call GetSection(config,section,'glaciers')
       if (associated(section)) then
          call handle_glaciers(section, model)
       end if
    endif

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file.
    ! (Glint restart variables are handled separately by Glint setup routines.)
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13

    call define_glide_restart_variables(model)

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
    call print_glaciers(model)

  end subroutine glide_printconfig

!-------------------------------------------------------------------------
    
  subroutine glide_scale_params(model)
    !> scale parameters
    use glide_types
    use glimmer_physcon,  only: scyr
    use glimmer_paramets, only: thk0, tim0, len0, vel0, vis0, acc0, tau0

    implicit none

    type(glide_global_type)  :: model !> model instance

    model%numerics%dttem = model%numerics%ntem * model%numerics%tinc   

    ! convert dt and dttem from yr to scaled time units
    model%numerics%dt     = model%numerics%tinc * scyr / tim0   
    model%numerics%dttem  = model%numerics%dttem * scyr / tim0   

    ! allow for subcycling of ice transport
    model%numerics%dt_transport = model%numerics%dt / real(model%numerics%subcyc, dp)

    model%numerics%thklim = model%numerics%thklim / thk0
    model%numerics%thklim_temp = model%numerics%thklim_temp / thk0
    model%numerics%thck_gradient_ramp = model%numerics%thck_gradient_ramp / thk0

    model%numerics%dew = model%numerics%dew / len0
    model%numerics%dns = model%numerics%dns / len0

    ! scale calving parameters
    model%calving%marine_limit = model%calving%marine_limit / thk0
    model%calving%timescale = model%calving%timescale * scyr                ! convert from yr to s
    model%calving%cliff_timescale = model%calving%cliff_timescale * scyr    ! convert from yr to s
    model%calving%eigencalving_constant = model%calving%eigencalving_constant / scyr    ! convert from m/yr/Pa to m/s/Pa
    model%calving%damage_constant = model%calving%damage_constant / scyr    ! convert from yr^{-1} to s^{-1}
    model%calving%lateral_rate_max = model%calving%lateral_rate_max / scyr  ! convert from m/yr to m/s

    ! scale periodic offsets for ISMIP-HOM
    model%numerics%periodic_offset_ew = model%numerics%periodic_offset_ew / thk0
    model%numerics%periodic_offset_ns = model%numerics%periodic_offset_ns / thk0

    ! scale glide basal traction parameters
    model%velowk%trc0   = vel0 * len0 / (thk0**2)
    model%velowk%btrac_const = model%paramets%btrac_const/model%velowk%trc0/scyr
    model%velowk%btrac_max   = model%paramets%btrac_max / model%velowk%trc0/scyr    
    model%velowk%btrac_slope = model%paramets%btrac_slope*acc0/model%velowk%trc0

    ! scale basal melting parameters (1/yr -> 1/s, or m/yr -> m/s)
    model%basal_melt%bmlt_float_omega = model%basal_melt%bmlt_float_omega / scyr
    model%basal_melt%bmlt_float_const = model%basal_melt%bmlt_float_const / scyr
    model%basal_melt%bmlt_float_depth_meltmax = model%basal_melt%bmlt_float_depth_meltmax / scyr
    model%basal_melt%bmlt_float_depth_frzmax = model%basal_melt%bmlt_float_depth_frzmax / scyr
    model%basal_melt%bmlt_float_depth_meltmin = model%basal_melt%bmlt_float_depth_meltmin / scyr

    ! scale basal inversion parameters
    model%inversion%babc_timescale = model%inversion%babc_timescale * scyr    ! convert yr to s
    model%inversion%thck_threshold = model%inversion%thck_threshold / thk0
    model%inversion%thck_flotation_buffer = model%inversion%thck_flotation_buffer / thk0
    model%inversion%dbmlt_dtemp_scale = model%inversion%dbmlt_dtemp_scale / scyr   ! m/yr/degC to m/s/degC
    model%inversion%bmlt_basin_timescale = model%inversion%bmlt_basin_timescale * scyr   ! yr to s

    ! scale SMB/acab parameters
    model%climate%overwrite_acab_value = model%climate%overwrite_acab_value*tim0/(scyr*thk0)
    model%climate%overwrite_acab_minthck = model%climate%overwrite_acab_minthck / thk0

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

  subroutine glide_get_zocn(model,config)

    ! Read ocean grid information, if present, from the config file.
    ! Called after glide_readconfig

    use glide_types
    use glimmer_config
    use glimmer_log
    use cism_parallel, only: main_task

    implicit none

    type(glide_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file

    ! local variables
    type(ConfigSection), pointer :: section
    character(len=512) :: message
    character(len=16) :: message_tmp
    integer :: k

    ! Check for section [grid_ocn]
    call GetSection(config,section,'grid_ocn')

    if (associated(section)) then

       call GetValue(section,'nzocn',model%ocean_data%nzocn)
       call GetValue(section,'dzocn',model%ocean_data%dzocn)     ! m
       call GetValue(section,'nbasin',model%ocean_data%nbasin)

       call write_log(' ')
       write(message,*) 'number of ocean levels        : ',model%ocean_data%nzocn
       call write_log(trim(message))

       if (.not.associated(model%ocean_data%zocn)) &
            allocate(model%ocean_data%zocn(model%ocean_data%nzocn))

       if ( size(model%ocean_data%zocn) /= model%ocean_data%nzocn )then
          write(message,*) 'size of zocn : ',size(model%ocean_data%zocn)
          call write_log(trim(message))
          call write_log (' zocn is allocated to be too small', GM_FATAL)
       end if

       ! There are two ways to get zocn levels:
       ! (1) Compute uniform levels based on dzocn.  This is done if dzocn is set to a nonzero value in the config file.
       ! (2) Load zocn levels from the config file.  This is done if dzocn = 0 (the default value).
       ! Note: By convention, k = 1 is the top level, and z becomes more negative with increasing k.
       !       If the input zocn levels do not satisfy this criterion (or if dzocn and zocn are not in the config file),
       !        The code aborts with an error message.

       if (model%ocean_data%dzocn > 0.0d0) then
          call write_log (' Computing zocn levels based on dzocn')
          write(message,*) 'dz (m) for ocean levels       : ',model%ocean_data%dzocn
          call write_log(trim(message))
          do k = 1, model%ocean_data%nzocn
             model%ocean_data%zocn(k) = -model%ocean_data%dzocn * (real(k,dp) - 0.5d0)
          enddo
       else
          call write_log (' Reading zocn levels from config file')
          call GetValue(section,'zocn',model%ocean_data%zocn, model%ocean_data%nzocn)
          do k = 2, model%ocean_data%nzocn
             if (model%ocean_data%zocn(k-1) - model%ocean_data%zocn(k) < 1.0d0) then
                write(message,*) 'nzocn, zocn =', model%ocean_data%nzocn, model%ocean_data%zocn(:)
                call write_log(trim(message))
                write(message,*) 'Must have zocn decreasing with increasing k in the [grid_ocn] section'
                call write_log(message, GM_FATAL)
             endif
          enddo
       endif

       if (model%ocean_data%nzocn > 1) then
          call write_log('')
          call write_log('zocn levels (m):')
          call write_log('------------------')
          message = ''
          do k = 1, model%ocean_data%nzocn
             write(message_tmp,'(f8.1)') model%ocean_data%zocn(k)
             message = trim(message)//trim(message_tmp)
          enddo
          call write_log(trim(message))
          call write_log('')
       endif

       if (model%ocean_data%nbasin >= 1) then
          call write_log('')
          write(message,*) 'number of ocean basins: ', model%ocean_data%nbasin
          call write_log(trim(message))
       else
          call write_log('No ocean basins')
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

  end subroutine glide_get_zocn

!--------------------------------------------------------------------------------

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

    if (model%general%nx_block > 0 .and. model%general%ny_block> 0) then
       write(message,*) 'nx_block        :', model%general%nx_block
       call write_log(trim(message))
       write(message,*) 'ny_block        :', model%general%ny_block
       call write_log(trim(message))
    endif

    if (model%general%global_bc==GLOBAL_BC_PERIODIC) then
       write(message,*) 'Periodic global boundary conditions'
       call write_log(trim(message))
    elseif (model%general%global_bc==GLOBAL_BC_OUTFLOW) then
       write(message,*) 'Outflow global boundary conditions; scalars in global halo will be set to zero'
       call write_log(trim(message))
    elseif (model%general%global_bc==GLOBAL_BC_NO_PENETRATION) then
       write(message,*) 'No-penetration global boundary conditions; outflow set to zero at global boundaries'
       call write_log(trim(message))
    elseif (model%general%global_bc==GLOBAL_BC_NO_ICE) then
       write(message,*) 'No-ice global boundary conditions; scalars set to zero adjacent to global boundaries'
       call write_log(trim(message))
    else
       write(message,*) ('Error: Invalid option for global_bc')
       call write_log(trim(message), GM_FATAL)
    endif

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
    call GetValue(section,'subcyc',model%numerics%subcyc)
    call GetValue(section,'adaptive_cfl_threshold', model%numerics%adaptive_cfl_threshold)
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

    if (model%options%whichdycore == DYCORE_GLIDE) then  ! Glide option only
       write(message,*) 'profile frequency   : ',model%numerics%profile_period
       call write_log(message)
    endif

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
    call GetValue(section,'bmlt_float',model%options%whichbmlt_float)
    call GetValue(section,'bmlt_float_thermal_forcing_param',model%options%bmlt_float_thermal_forcing_param)
    call GetValue(section,'ocean_data_domain',model%options%ocean_data_domain)
    call GetValue(section,'ocean_data_extrapolate',model%options%ocean_data_extrapolate)
    call GetValue(section,'enable_bmlt_anomaly',model%options%enable_bmlt_anomaly)
    call GetValue(section,'basal_mass_balance',model%options%basal_mbal)
    call GetValue(section,'smb_input',model%options%smb_input)
    call GetValue(section,'smb_input_function',model%options%smb_input_function)
    call GetValue(section,'artm_input_function',model%options%artm_input_function)
    call GetValue(section,'nlev_smb',model%climate%nlev_smb)
    call GetValue(section,'enable_acab_anomaly',model%options%enable_acab_anomaly)
    call GetValue(section,'enable_artm_anomaly',model%options%enable_artm_anomaly)
    call GetValue(section,'overwrite_acab',model%options%overwrite_acab)
    call GetValue(section,'gthf',model%options%gthf)
    call GetValue(section,'isostasy',model%options%isostasy)
    call GetValue(section,'marine_margin',model%options%whichcalving)
    call GetValue(section,'calving_init',model%options%calving_init)
    call GetValue(section,'calving_domain',model%options%calving_domain)
    call GetValue(section,'apply_calving_mask', model%options%apply_calving_mask)
    call GetValue(section,'remove_icebergs', model%options%remove_icebergs)
    call GetValue(section,'remove_isthmuses', model%options%remove_isthmuses)
    call GetValue(section,'expand_calving_mask', model%options%expand_calving_mask)
    call GetValue(section,'limit_marine_cliffs', model%options%limit_marine_cliffs)
    call GetValue(section,'cull_calving_front', model%options%cull_calving_front)
    call GetValue(section,'adjust_input_thickness', model%options%adjust_input_thickness)
    call GetValue(section,'smooth_input_topography', model%options%smooth_input_topography)
    call GetValue(section,'smooth_input_usrf', model%options%smooth_input_usrf)
    call GetValue(section,'adjust_input_topography', model%options%adjust_input_topography)
    call GetValue(section,'read_lat_lon',model%options%read_lat_lon)
    call GetValue(section,'dm_dt_diag',model%options%dm_dt_diag)
    call GetValue(section,'diag_minthck',model%options%diag_minthck)
    call GetValue(section,'vertical_integration',model%options%whichwvel)
    call GetValue(section,'periodic_ew',model%options%periodic_ew)
    call GetValue(section,'sigma',model%options%which_sigma)
    call GetValue(section,'ioparams',model%funits%ncfile)

    !Note: Previously, the terms 'hotstart' and 'restart' were both supported in the config file.
    !      Going forward, only 'restart' is supported.
    call GetValue(section,'restart',model%options%is_restart)

    call GetValue(section,'restart_extend_velo',model%options%restart_extend_velo)

  end subroutine handle_options

!--------------------------------------------------------------------------------

  !Higher order options
  subroutine handle_ho_options(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type) :: model
    
    call GetValue(section, 'compute_blocks',              model%options%compute_blocks)
    call GetValue(section, 'which_ho_efvs',               model%options%which_ho_efvs)
    call GetValue(section, 'which_ho_disp',               model%options%which_ho_disp)
    call GetValue(section, 'which_ho_thermal_timestep',   model%options%which_ho_thermal_timestep)
    call GetValue(section, 'which_ho_babc',               model%options%which_ho_babc)
    call GetValue(section, 'use_c_space_factor',          model%options%use_c_space_factor)
    call GetValue(section, 'which_ho_beta_limit',         model%options%which_ho_beta_limit)
    call GetValue(section, 'which_ho_powerlaw_c',         model%options%which_ho_powerlaw_c)
    call GetValue(section, 'which_ho_coulomb_c',          model%options%which_ho_coulomb_c)
    call GetValue(section, 'which_ho_bmlt_basin_inversion', model%options%which_ho_bmlt_basin_inversion)
    call GetValue(section, 'which_ho_bwat',               model%options%which_ho_bwat)
    call GetValue(section, 'ho_flux_routing_scheme',      model%options%ho_flux_routing_scheme)
    call GetValue(section, 'which_ho_effecpress',         model%options%which_ho_effecpress)
    call GetValue(section, 'which_ho_resid',              model%options%which_ho_resid)
    call GetValue(section, 'which_ho_nonlinear',          model%options%which_ho_nonlinear)
    call GetValue(section, 'which_ho_sparse',             model%options%which_ho_sparse)
    call GetValue(section, 'which_ho_approx',             model%options%which_ho_approx)
    call GetValue(section, 'which_ho_precond',            model%options%which_ho_precond)
    call GetValue(section, 'which_ho_gradient',           model%options%which_ho_gradient)
    call GetValue(section, 'which_ho_gradient_margin',    model%options%which_ho_gradient_margin)
    call GetValue(section, 'which_ho_vertical_remap',     model%options%which_ho_vertical_remap)
    call GetValue(section, 'which_ho_assemble_beta',      model%options%which_ho_assemble_beta)
    call GetValue(section, 'which_ho_assemble_taud',      model%options%which_ho_assemble_taud)
    call GetValue(section, 'which_ho_assemble_bfric',     model%options%which_ho_assemble_bfric)
    call GetValue(section, 'which_ho_assemble_lateral',   model%options%which_ho_assemble_lateral)
    call GetValue(section, 'which_ho_calving_front',      model%options%which_ho_calving_front)
    call GetValue(section, 'which_ho_ground',             model%options%which_ho_ground)
    call GetValue(section, 'which_ho_fground_no_glp',     model%options%which_ho_fground_no_glp)
    call GetValue(section, 'which_ho_ground_bmlt',        model%options%which_ho_ground_bmlt)
    call GetValue(section, 'which_ho_flotation_function', model%options%which_ho_flotation_function)
    call GetValue(section, 'block_inception',             model%options%block_inception)
    call GetValue(section, 'remove_ice_caps',             model%options%remove_ice_caps)
    call GetValue(section, 'force_retreat',               model%options%force_retreat)
    call GetValue(section, 'which_ho_ice_age',            model%options%which_ho_ice_age)
    call GetValue(section, 'enable_glaciers',             model%options%enable_glaciers)
    call GetValue(section, 'glissade_maxiter',            model%options%glissade_maxiter)
    call GetValue(section, 'linear_solve_ncheck',         model%options%linear_solve_ncheck)
    call GetValue(section, 'linear_maxiters',             model%options%linear_maxiters)
    call GetValue(section, 'linear_tolerance',            model%options%linear_tolerance)

  end subroutine handle_ho_options

!--------------------------------------------------------------------------------

  ! Handles external dycore options -- Doug Ranken 03/26/12
  subroutine handle_dycore_options(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type) :: model
    
    call GetValue(section, 'external_dycore_type', model%options%external_dycore_type)
    call GetValue(section, 'dycore_input_file',  model%options%dycore_input_file)
    if (model%options%external_dycore_type .eq. 1) model%options%whichdycore = 4 ! DYCORE_BISICLES
    if (model%options%external_dycore_type .eq. 2) model%options%whichdycore = 3 ! DYCORE_ALBANYFELIX

    print *,"In handle_dycore_options, external dycore type, input file = ", &
             model%options%external_dycore_type,model%options%dycore_input_file 
    ! print *,"In handle_dycore_options, whichdycore = ",model%options%whichdycore
  end subroutine handle_dycore_options

!--------------------------------------------------------------------------------

  subroutine print_options(model)

    use glide_types
    use glimmer_log
    use cism_parallel, only: tasks

    implicit none

    type(glide_global_type)  :: model
    character(len=500) :: message

    ! basic Glide/Glimmer options

    character(len=*), dimension(0:4), parameter :: dycore = (/ &
         'glide              ', &  ! Glimmer SIA
         'glam               ', &  ! Payne-Price finite difference; no longer supported
         'glissade           ', &  ! prototype finite element
         'albany-felix       ', &  ! External Albany-FELIX finite element
         'bisicles           ' /)  ! External BISICLES-Chombo FVM 

    character(len=*), dimension(0:2), parameter :: compute_blocks = (/ &
         'compute on all blocks              ', &
         'compute on active blocks only      ', &
         'inquire the number of active blocks'/)

    character(len=*), dimension(0:5), parameter :: evolution = (/ &
         'pseudo-diffusion                      ', &
         'ADI scheme                            ', &
         'iterated diffusion                    ', &
         'incremental remapping                 ', &   
         '1st order upwind                      ', &
         'thickness fixed at initial value      '  /)

    character(len=*), dimension(0:3), parameter :: temperature = (/ &
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

    character(len=*), dimension(0:6), parameter :: which_bmlt_float = (/ &
         'none                                  ', &
         'MISMIP+ melt rate profile             ', &
         'constant melt rate                    ', &
         'depth-dependent melt rate             ', &
         'melt rate from external file          ', &
         'melt rate from MISOMIP T/S profile    ', &   ! not supported
         'melt rate from thermal forcing        ' /)

    character(len=*), dimension(0:3), parameter :: bmlt_float_thermal_forcing_param = (/ &
         'quadratic function of thermal forcing     ', &
         'ISMIP6 local quadratic                    ', &
         'ISMIP6 nonlocal quadratic                 ', &
         'ISMIP6 nonlocal quadratic, slope-dependent' /)

    character(len=*), dimension(0:2), parameter :: ocean_data_domain = (/ &
         'ocean data computed internally by CISM', &
         'ocean data read from external file    ', &
         'ocean data from coupler via Glad      ' /)

    character(len=*), dimension(0:1), parameter :: ocean_data_extrapolate = (/ &
         'ocean data not extrapolated to cavities', &
         'ocean data extrapolated to cavities    ' /)

    character(len=*), dimension(0:1), parameter :: smb_input = (/ &
         'SMB input in units of m/yr ice  ', &
         'SMB input in units of mm/yr w.e.' /)

    character(len=*), dimension(0:2), parameter :: smb_input_function = (/ &
         'SMB input as function of (x,y)              ', &
         'SMB and d(SMB)/dz input as function of (x,y)', &
         'SMB input as function of (x,y,z)            ' /)

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

    ! NOTE: Set gthf = 1 in the config file to read the geothermal heat flux from an input file.
    !       Otherwise it will be overwritten, even if the 'bheatflx' field is present.

    character(len=*), dimension(0:2), parameter :: gthf = (/ &
         'uniform geothermal flux         ', &
         'read flux from file, if present ', &
         'compute flux from diffusion eqn ' /)

    ! NOTE: This option has replaced the old do_isos option
    character(len=*), dimension(0:1), parameter :: isostasy = (/ &
         'no isostasy calculation         ', &
         'compute isostasy with model     ' /)

    !TODO - Change 'marine_margin' to 'calving'?  Would have to modify many config files
    character(len=*), dimension(0:9), parameter :: marine_margin = (/ &
         'no calving law                  ', &
         'remove all floating ice         ', &
         'remove fraction of floating ice ', &
         'relaxed bedrock threshold       ', &
         'present bedrock threshold       ', &
         'calving based on grid location  ', &
         'ice thickness threshold         ', &
         'eigencalving scheme             ', & 
         'damage-based calving scheme     ', & 
         'Huybrechts grounding-line scheme' /) 

    character(len=*), dimension(0:1), parameter :: init_calving = (/ &
         'no calving at initialization    ', &
         'ice calves at initialization    ' /)

    character(len=*), dimension(0:1), parameter :: domain_calving = (/ &
         'calving only at the ocean edge             ',  &
         'calving in all cells where criterion is met'/)

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

    character(len=*), dimension(0:2), parameter :: ho_whichefvs = (/ &
         'constant value                 ', & 
         'multiple of flow factor        ', &
         'nonlinear, from eff strain rate' /)

    character(len=*), dimension(-1:1), parameter :: ho_whichdisp = (/ &
         'no dissipation                    ', &
         '0-order SIA                       ', &
         'first-order model (Blatter-Pattyn)' /)

    character(len=*), dimension(0:2), parameter :: ho_whichthermal_timestep = (/ &
         'vertical thermal solve before transport    ', &
         'vertical thermal solve after transport     ', &
         'vertical thermal solve split into two parts' /)

    character(len=*), dimension(0:15), parameter :: ho_whichbabc = (/ &
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
         'Coulomb friction law w/ effec press, const flwa_b', &
         'min of Coulomb stress and power-law stress (Tsai)', &
         'power law using effective pressure               ', &
         'simple pattern of beta                           ', &
         'till yield stress (Picard)                       ' /)

    character(len=*), dimension(0:1), parameter :: ho_whichbeta_limit = (/ &
         'absolute beta limit based on beta_grounded_min   ', &
         'beta is limited, then scaled by f_ground_cell    ' /)

    character(len=*), dimension(0:2), parameter :: ho_powerlaw_c = (/ &
         'spatially uniform friction parameter Cp ', &
         'friction parameter Cp found by inversion', &
         'friction parameter Cp read from file    ' /)

    character(len=*), dimension(0:3), parameter :: ho_coulomb_c = (/ &
         'spatially uniform friction parameter Cc ', &
         'friction parameter Cc found by inversion', &
         'friction parameter Cc read from file    ', &
         'Cc is a function of bed elevation       ' /)

    character(len=*), dimension(0:2), parameter :: ho_bmlt_basin_whichinversion = (/ &
         'no inversion for basin-based basal melting parameters      ', &
         'invert for basin-based basal melting parameters            ', &
         'apply basin basal melting parameters from earlier inversion' /)

    character(len=*), dimension(0:3), parameter :: ho_whichbwat = (/ &
         'zero basal water depth                          ', &
         'constant basal water depth                      ', &
         'basal water depth computed from local till model', &
         'steady-state water routing with flux calculation' /)

    character(len=*), dimension(0:2), parameter :: ho_flux_routing_scheme = (/ &
         'D8; route flux to lowest-elevation neighbor      ', &
         'Dinf; route flux to two lower-elevation neighbors', &
         'FD8; route flux to all lower-elevation neighbors ' /)

    character(len=*), dimension(0:4), parameter :: ho_whicheffecpress = (/ &
         'full overburden pressure                             ', &
         'reduced effecpress near pressure melting point       ', &
         'reduced effecpress where bwat > 0 (ramp)             ', &
         'reduced effecpress where bwatflx > 0                 ', &
         'reduced effecpress where bwat > 0 (B/vP)             '/)

    character(len=*), dimension(0:1), parameter :: which_ho_nonlinear = (/ &
         'use standard Picard iteration          ', &
         'use Picard iteration with acceleration '/)

    character(len=*), dimension(0:4), parameter :: ho_whichresid = (/ &
         'max value                   ', &
         'max value ignoring ubas     ', &
         'mean value                  ', &
         'L2 norm of Ax-b = resid     ', &
         'relative L2 norm, |Ax-b|/|b|' /)

    character(len=*), dimension(-1:4), parameter :: ho_whichsparse = (/ &
         'PCG with incomplete Cholesky preconditioner', &
         'BiCG with LU preconditioner                ', &
         'GMRES with LU preconditioner               ', &
         'Native PCG solver, standard                ', &
         'Native PCG solver, Chronopoulos-Gear       ', &
         'Trilinos interface                         '/)

    character(len=*), dimension(-1:5), parameter :: ho_whichapprox = (/ &
         'SIA only (glissade_velo_sia)                     ', &
         'SIA only (glissade_velo_higher)                  ', &
         'SSA only (glissade_velo_higher)                  ', &
         'Blatter-Pattyn HO (glissade_velo_higher)         ', &
         'Depth-integrated L1L2 (glissade_velo_higher)     ', &
         'Depth-integrated viscosity (glissade_velo_higher)', &
         'Hybrid SIA/SSA                                   ' /)

    character(len=*), dimension(0:4), parameter :: ho_whichprecond = (/ &
         'No preconditioner (native PCG)                ', &
         'Diagonal preconditioner (native PCG)          ', &
         'SIA preconditioner (native PCG)               ', &
         'Local tridiagonal preconditioner (native PCG) ', &
         'Global tridiagonal preconditioner (native PCG)' /)

    character(len=*), dimension(0:2), parameter :: ho_whichgradient = (/ &
         'centered gradient (glissade)             ', &
         'first-order upstream gradient (glissade) ', &
         'second-order upstream gradient (glissade)' /)

    character(len=*), dimension(0:2), parameter :: ho_whichgradient_margin = (/ &
         'compute edge gradient when either cell has ice         ', &
         'compute edge gradient when ice lies above ice-free land', & 
         'compute edge gradient when both cells have ice         ' /)

    character(len=*), dimension(0:1), parameter :: ho_whichvertical_remap = (/ &
         'first-order accurate  ', &
         'second-order accurate ' /)

    character(len=*), dimension(0:1), parameter :: ho_whichassemble_beta = (/ &
         'standard finite-element assembly (glissade dycore) ', &
         'use local beta at each vertex (glissade dycore)    '  /)

    character(len=*), dimension(0:1), parameter :: ho_whichassemble_taud = (/ &
         'standard finite-element assembly (glissade dycore)       ', &
         'use local driving stress at each vertex (glissade dycore)'  /)

    character(len=*), dimension(0:1), parameter :: ho_whichassemble_bfric = (/ &
         'standard finite-element assembly (glissade dycore)       ', &
         'use local basal friction at each vertex (glissade dycore)'  /)

    character(len=*), dimension(0:1), parameter :: ho_whichassemble_lateral = (/ &
         'standard finite-element assembly (glissade dycore)         ', &
         'use local thck and usrf on each cell face (glissade dycore)'  /)

    character(len=*), dimension(0:1), parameter :: ho_whichcalving_front = (/ &
         'no subgrid calving front parameterization      ', &
         'subgrid calving front scheme; inactive CF cells' /)

    character(len=*), dimension(0:2), parameter :: ho_whichground = (/ &
         'f_ground = 0 or 1; no GLP  (glissade dycore)               ', &
         'GLP for basal friction with 0 <= f_ground <= 1 for vertices', &
         'deluxe GLP, 0 <= f_ground <= 1 for both vertices and cells ' /)

    character(len=*), dimension(0:1), parameter :: ho_whichfground_no_glp = (/ &
         'f_ground = 1 at vertex if any neighbor cell is grounded   ', &
         'f_ground = 0 or 1 at vertex based on staggered f_flotation' /)

    character(len=*), dimension(0:2), parameter :: ho_whichground_bmlt = (/ &
         'no GLP for bmlt_float                        ', &
         'weigh bmlt_float by floating fraction of cell', &
         'set bmlt_float = 0 in partly grounded cells  ' /)

    character(len=*), dimension(0:4), parameter :: ho_whichflotation_function = (/ &
         'f_pattyn = (-rhoo*b)/(rhoi*H)       ', &
         '1/fpattyn = (rhoi*H)/(-rhoo*b)      ', &
         'linear = -b - (rhoi/rhoo)*H         ', &
         'modified linear = -b - (rhoi/rhoo)*H', &
         'modified linear; topg stdev addition'/)

    character(len=*), dimension(0:1), parameter :: ho_whichice_age = (/ &
         'ice age computation off', &
         'ice age computation on ' /)

    call write_log('Dycore options')
    call write_log('-------------')

    write(message,*) 'I/O parameter file      : ',trim(model%funits%ncfile)
    call write_log(message)

    if (model%options%whichdycore < 0 .or. model%options%whichdycore >= size(dycore)) then
       call write_log('Error, dycore option out of range',GM_FATAL)
    end if
    write(message,*) 'Dycore                  : ',model%options%whichdycore,dycore(model%options%whichdycore)
    call write_log(message)

    ! unsupported dycore options
    if (model%options%whichdycore == DYCORE_GLAM) then
      call write_log('ERROR: The Glam dycore is no longer supported', GM_FATAL)
    endif
    if (model%options%whichdycore == DYCORE_ALBANYFELIX) then
      call write_log('WARNING: Albany-FELIX dycore is not currently scientifically supported.  USE AT YOUR OWN RISK.', GM_WARNING)
    endif
    if (model%options%whichdycore == DYCORE_BISICLES) then
      call write_log('WARNING: BISICLES dycore is not currently scientifically supported.  USE AT YOUR OWN RISK.', GM_WARNING)
    endif

    ! Forbidden options associated with the Glide dycore
    if (model%options%whichdycore == DYCORE_GLIDE) then

       if (model%options%whichevol == EVOL_INC_REMAP     .or.  &
           model%options%whichevol == EVOL_UPWIND        .or.  &
           model%options%whichevol == EVOL_NO_THICKNESS) then
          call write_log('Error, Glissade thickness evolution options cannot be used with Glide dycore', GM_FATAL)
       endif

       if (model%options%whichtemp == TEMP_ENTHALPY) then
          call write_log('Error, Enthalpy scheme cannot be used with Glide dycore', GM_FATAL)
       endif

       if (tasks > 1) then
          call write_log('Error, Glide dycore not supported for runs with more than one processor', GM_FATAL)
       end if

       if (model%options%whichevol == EVOL_ADI) then
          call write_log('WARNING: exact restarts are not currently possible with ADI evolution', GM_WARNING)
       endif

       if (model%options%overwrite_acab /= OVERWRITE_ACAB_NONE) then
          call write_log('WARNING: overwrite_acab option will be ignored for Glide dycore', GM_WARNING)
       endif

    else   ! forbidden evolution options with dycores other than Glide

       if (model%options%whichevol == EVOL_PSEUDO_DIFF .or.  &
           model%options%whichevol == EVOL_ADI         .or.  &
           model%options%whichevol == EVOL_DIFFUSION) then
          call write_log('Error, Glide thickness evolution options cannot be used with Glissade dycore', GM_FATAL)
       endif

    endif

    ! Forbidden options for running in parallel
    if (tasks > 1 .and. (model%options%which_ho_sparse==HO_SPARSE_BICG  .or.  &
                         model%options%which_ho_sparse==HO_SPARSE_GMRES .or.  &
                         model%options%which_ho_sparse==HO_SPARSE_PCG_INCH) ) then
       call write_log('Error, SLAP solver not supported for more than one processor', GM_FATAL)
    end if

    if (tasks > 1 .and. model%options%whichbwat==BWATER_FLUX) then
       call write_log('Error, flux-based basal water option not yet supported for more than one processor', GM_FATAL)
    endif

    ! Forbidden options associated with Glissade dycore
   
    if (model%options%whichdycore == DYCORE_GLISSADE) then 

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
             call write_log('Warning, cannot use SIA preconditioning for 2D solve')
             call write_log('Defaulting to diagonal preconditioner')
          endif
       endif

       if (model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
          
          if (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then
             call write_log('Error, cannot use no-penetration BC with local SIA solver', GM_FATAL)
          endif

          if (model%options%which_ho_disp == HO_DISP_FIRSTORDER ) then
             model%options%which_ho_disp = HO_DISP_SIA
             call write_log('Warning: Cannot use first-order dissipation with local SIA solver')
             call write_log('Defaulting to SIA dissipation')
          endif

       endif  ! Glissade local SIA solver

    endif

    if (model%options%whichdycore /= DYCORE_GLISSADE) then 

       if (model%options%which_ho_sparse == HO_SPARSE_PCG_STANDARD .or.   &
           model%options%which_ho_sparse == HO_SPARSE_PCG_CHRONGEAR) then
          call write_log('Error, native PCG solver requires glissade dycore', GM_FATAL)
       endif

       if (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then
          call write_log('Error, no-penetration BC requires glissade dycore', GM_FATAL)
       endif

    endif

    ! Config specific to Albany-Felix dycore   
    if (model%options%whichdycore == DYCORE_ALBANYFELIX) then
       call write_log('WARNING: Albany-FELIX dycore requires external libraries, and it is still in development!!!', GM_WARNING)
    endif

    !NOTE : Old option 3 (TEMP_REMAP_ADV) has been removed.
    ! If this has been set, then change to option 1 (TEMP_PROGNOSTIC), which applies to any dycore.

    if (model%options%whichtemp < 0 .or. model%options%whichtemp >= size(temperature)) then
       call write_log('Error, temperature option out of range',GM_FATAL)
    end if
    write(message,*) 'temperature calculation : ',model%options%whichtemp,temperature(model%options%whichtemp)
    call write_log(message)

    ! unsupported temperature options
    if (model%options%whichtemp == TEMP_ENTHALPY) then
      call write_log('WARNING: Enthalpy-based formulation for solving temperature evolution is not currently &
           &scientifically supported. USE AT YOUR OWN RISK.', GM_WARNING)
    endif

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
      call write_log('WARNING: Steady state routing basal_water option is not currently scientifically supported.  &
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

    if (model%options%calving_domain < 0 .or. model%options%calving_domain >= size(domain_calving)) then
       call write_log('Error, calving_domain out of range',GM_FATAL)
    end if
    write(message,*) 'calving_domain          : ', model%options%calving_domain, domain_calving(model%options%calving_domain)
    call write_log(message)
    
    ! dycore-dependent options; most of these are supported for Glissade only

    if (model%options%whichdycore == DYCORE_GLISSADE) then

       if (model%options%apply_calving_mask) then
          call write_log(' A calving mask will be applied')
       endif

       if (model%options%remove_icebergs) then
          call write_log(' Icebergs will be removed')
       else
          call write_log(' Icebergs will not be removed')
       endif

       if (model%options%remove_isthmuses) then
          if (.not.model%options%remove_icebergs) then
             model%options%remove_icebergs = .true.
             write(message,*) ' Setting remove_icebergs = T for stability when remove_isthmuses = T'
             call write_log(message)
          endif
          call write_log(' Isthmuses will be removed')
       endif
       
       if (model%options%expand_calving_mask) then
          if (model%options%whichcalving == CALVING_GRID_MASK .or. model%options%apply_calving_mask) then
             call write_log(' The calving mask will be expanded to include floating ice in select basins')
          else
             call write_log(' Not using a calving_mask; expand_calving_mask = T will be ignored')
          endif
       endif

       if (model%options%limit_marine_cliffs) then
          call write_log(' The thickness of marine ice cliffs will be limited')
          call write_log(message)
       else
          call write_log(' The thickness of marine ice cliffs will not be limited')
       endif

       if (model%options%cull_calving_front) then
          write(message,*) ' Calving-front cells will be culled', model%calving%ncull_calving_front, 'times at initialization'
          call write_log(message)
       else
          call write_log(' Calving-front cells will not be culled at initialization')
       endif

       if (model%options%whichcalving == CALVING_FLOAT_FRACTION) then
          write(message,*) 'WARNING: calving float fraction option deprecated with Glissade_dycore; set calving_timescale instead'
          call write_log(message, GM_WARNING)
       endif

       if (model%options%adjust_input_thickness) then
          write(message,*) ' Input ice thickness will be adjusted based on surface and bed topography'
          call write_log(message)
       endif

       if (model%options%smooth_input_usrf) then
          write(message,*) ' Input usrf will be smoothed'
          call write_log(message)
       endif

       if (model%options%smooth_input_topography) then
          write(message,*) ' Input topography will be smoothed'
          call write_log(message)
       endif

       if (model%options%adjust_input_topography) then
          write(message,*) ' Input topography in a selected region will be adjusted'
          call write_log(message)
       endif

       if (model%options%read_lat_lon) then
          write(message,*) ' Lat and lon fields will be read from input files and written to restart'
          call write_log(message)
       endif

    else   ! not Glissade

       if (model%options%whichcalving == CALVING_THCK_THRESHOLD) then
          call write_log('Error, calving thickness threshold option is supported for Glissade dycore only', GM_FATAL)
       endif
       if (model%options%whichcalving == EIGENCALVING) then
          call write_log('Error, eigencalving option is supported for Glissade dycore only', GM_FATAL)
       endif
       if (model%options%whichcalving == CALVING_GRID_MASK) then
          call write_log('Error, calving grid mask option is supported for Glissade dycore only', GM_FATAL)
       endif
       if (model%options%whichcalving == CALVING_DAMAGE) then
          call write_log('Error, calving damage option is supported for Glissade dycore only', GM_FATAL)
       endif
       if (model%options%calving_domain /= CALVING_DOMAIN_OCEAN_EDGE) then
          write(message,*) 'WARNING: calving domain can be selected for Glissade dycore only; user selection ignored'
          call write_log(message, GM_WARNING)
       endif
       if (model%calving%timescale > 0.0d0) then
          write(message,*) 'WARNING: calving timescale option supported for Glissade dycore only; user selection ignored'
          call write_log(message, GM_WARNING)
       endif
       if (model%options%adjust_input_thickness) then
          write(message,*) 'WARNING: adjust_input_thickness supported for Glissade dycore only; user selection ignored'
          call write_log(message, GM_WARNING)
       endif
       if (model%options%smooth_input_topography) then
          write(message,*) 'WARNING: smooth_input_topography supported for Glissade dycore only; user selection ignored'
          call write_log(message, GM_WARNING)
       endif
       if (model%options%adjust_input_topography) then
          write(message,*) 'WARNING: adjust_input_topography supported for Glissade dycore only; user selection ignored'
          call write_log(message, GM_WARNING)
       endif

    endif  ! dycore

    if (model%options%whichbtrc < 0 .or. model%options%whichbtrc >= size(slip_coeff)) then
       call write_log('Error, slip_coeff out of range',GM_FATAL)
    end if

    if (model%options%whichdycore == DYCORE_GLIDE .or.  &
         (model%options%whichdycore == DYCORE_GLISSADE .and. model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) ) then
       write(message,*) 'slip_coeff              : ', model%options%whichbtrc, slip_coeff(model%options%whichbtrc)
       call write_log(message)

       !Note: Not all basal traction options are supported for the Glissade SIA solver
       if (model%options%whichdycore == DYCORE_GLISSADE .and. model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
          if (model%options%whichbtrc > BTRC_CONSTANT_BPMP) then
             call write_log('Error, slip_coeff out of range for Glissade dycore',GM_FATAL)
          end if
       endif
    endif

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

    if (model%options%whichwvel /= VERTINT_STANDARD .and. model%options%whichdycore /= DYCORE_GLIDE) then
       call write_log('Error, only standard vertical velocity calculation is supported for higher-order dycores.',GM_FATAL)
    end if

    write(message,*) 'vertical_integration    : ',model%options%whichwvel,vertical_integration(model%options%whichwvel)
    call write_log(message)

    if (model%options%whichbmlt_float < 0 .or. model%options%whichbmlt_float >= size(which_bmlt_float)) then
       call write_log('Error, bmlt_float out of range',GM_FATAL)
    end if

    write(message,*) 'basal melt, floating ice: ',model%options%whichbmlt_float, which_bmlt_float(model%options%whichbmlt_float)
    call write_log(message)

    if (model%options%whichbmlt_float == BMLT_FLOAT_MISOMIP) then
       call write_log('Error, BMLT_FLOAT_MISOMIP option is not supported', GM_FATAL)
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then
       write(message,*) 'melt parameterization   : ', model%options%bmlt_float_thermal_forcing_param, &
            bmlt_float_thermal_forcing_param(model%options%bmlt_float_thermal_forcing_param)
       call write_log(message)
       write(message,*) 'ocean data domain       : ', model%options%ocean_data_domain, &
            ocean_data_domain(model%options%ocean_data_domain)
       call write_log(message)
       write(message,*) 'ocean data extrapolate  : ', model%options%ocean_data_extrapolate, &
            ocean_data_extrapolate(model%options%ocean_data_extrapolate)
       call write_log(message)
    endif

    if (model%options%basal_mbal < 0 .or. model%options%basal_mbal >= size(b_mbal)) then
       call write_log('Error, basal_mass_balance out of range',GM_FATAL)
    end if

    if (model%options%enable_bmlt_anomaly) then
       call write_log('bmlt_float anomaly forcing is enabled')
    endif

    write(message,*) 'basal mass balance      : ',model%options%basal_mbal,b_mbal(model%options%basal_mbal)
    call write_log(message)

    if (model%options%smb_input < 0 .or. model%options%smb_input >= size(smb_input)) then
       call write_log('Error, smb_input option out of range',GM_FATAL)
    end if

    write(message,*) 'smb input units         : ',model%options%smb_input,smb_input(model%options%smb_input)
    call write_log(message)

    if (model%options%smb_input_function < 0 .or. model%options%smb_input_function >= size(smb_input_function)) then
       call write_log('Error, smb_input_function option out of range',GM_FATAL)
    end if

    write(message,*) 'smb input function      : ', &
         model%options%smb_input_function, smb_input_function(model%options%smb_input_function)
    call write_log(message)
    if (model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then
       write(message,*) 'number of SMB levels    : ', model%climate%nlev_smb
       call write_log(message)
       if (model%climate%nlev_smb < 2) then
          call write_log('Error, must have nlev_smb >= 2 for this input function', GM_FATAL)
       endif
    endif

    if (model%options%artm_input_function < 0 .or. model%options%artm_input_function >= size(artm_input_function)) then
       call write_log('Error, artm_input_function option out of range',GM_FATAL)
    end if

    write(message,*) 'artm input function     : ', &
         model%options%artm_input_function, artm_input_function(model%options%artm_input_function)
    call write_log(message)
    if (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XYZ) then
       write(message,*) 'number of artm levels   : ', model%climate%nlev_smb
       call write_log(message)
       if (model%climate%nlev_smb < 2) then
          call write_log('Error, must have nlev_smb >= 2 for this input function', GM_FATAL)
       endif
    endif

    if (model%options%enable_acab_anomaly) then
       call write_log('acab/SMB anomaly forcing is enabled')
    endif

    if (model%options%enable_artm_anomaly) then
       call write_log('artm anomaly forcing is enabled')
    endif

    if (model%options%overwrite_acab < 0 .or. model%options%overwrite_acab >= size(overwrite_acab)) then
       call write_log('Error, overwrite_acab option out of range',GM_FATAL)
    end if

    write(message,*) 'overwrite_acab          : ',model%options%overwrite_acab,overwrite_acab(model%options%overwrite_acab)
    call write_log(message)

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

    if (model%options%is_restart == RESTART_TRUE) then
       call write_log('Restarting model from a previous run')
       if (model%options%restart_extend_velo == RESTART_EXTEND_VELO_TRUE) then
          call write_log('Using extended velocity fields for restart')
       endif
    end if

    !HO options

    if (model%options%whichdycore /= DYCORE_GLIDE) then   ! glissade higher-order

       call write_log(' ')
       call write_log('Higher-order options:')
       call write_log('----------')

       if (model%options%compute_blocks /= 0) then
          write(message,*) 'compute_blocks          : ', model%options%compute_blocks, &
               compute_blocks(model%options%compute_blocks)
          call write_log(message)
          if (model%general%nx_block < 1 .or. model%general%ny_block < 1) then
             write(message,*) 'Must set nx_block and ny_block > 0 for this compute_blocks option'
             call write_log(message, GM_FATAL)
          endif
       endif
       if (model%options%compute_blocks < 0 .or. model%options%compute_blocks >= size(compute_blocks)) then
          call write_log('Error, compute_blocks out of range',GM_FATAL)
       end if

       write(message,*) 'ho_whichefvs            : ',model%options%which_ho_efvs,  &
                         ho_whichefvs(model%options%which_ho_efvs)
       call write_log(message)
       if (model%options%which_ho_efvs < 0 .or. model%options%which_ho_efvs >= size(ho_whichefvs)) then
          call write_log('Error, HO effective viscosity input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichdisp            : ',model%options%which_ho_disp,  &
                         ho_whichdisp(model%options%which_ho_disp)
       call write_log(message)
       if (model%options%which_ho_disp < -1 .or. model%options%which_ho_disp >= size(ho_whichdisp)-1) then
          call write_log('Error, HO dissipation input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichthermal_timestep: ',model%options%which_ho_thermal_timestep,  &
                         ho_whichthermal_timestep(model%options%which_ho_thermal_timestep)
       call write_log(message)
       if (model%options%which_ho_thermal_timestep < 0 .or. &
            model%options%which_ho_thermal_timestep >= size(ho_whichthermal_timestep)) then
          call write_log('Error, HO thermal timestep input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichbabc            : ',model%options%which_ho_babc,  &
                         ho_whichbabc(model%options%which_ho_babc)
       call write_log(message)
       if (model%options%which_ho_babc < 0 .or. model%options%which_ho_babc >= size(ho_whichbabc)) then
          call write_log('Error, HO basal BC input out of range', GM_FATAL)
       end if

       if (model%options%use_c_space_factor) then
          if (model%options%which_ho_babc == HO_BABC_COULOMB_FRICTION .or.  &
              model%options%which_ho_babc == HO_BABC_COULOMB_POWERLAW_SCHOOF .or. &
              model%options%which_ho_babc == HO_BABC_COULOMB_POWERLAW_TSAI) then
             write(message,*) 'Multiplying beta by C_space_factor'
             call write_log(message)
          else
             call write_log('Error, C_space_factor not supported for this choice of which_ho_babc', GM_FATAL)
          endif
       endif

       write(message,*) 'ho_whichbeta_limit      : ',model%options%which_ho_beta_limit,  &
                         ho_whichbeta_limit(model%options%which_ho_beta_limit)
       call write_log(message)
       if (model%options%which_ho_beta_limit < 0 .or. model%options%which_ho_beta_limit >= size(ho_whichbeta_limit)) then
          call write_log('Error, HO beta limit input out of range', GM_FATAL)
       end if

       ! basal friction options

       write(message,*) 'ho_powerlaw_c           : ',model%options%which_ho_powerlaw_c,  &
                         ho_powerlaw_c(model%options%which_ho_powerlaw_c)
       call write_log(message)
       if (model%options%which_ho_powerlaw_c < 0 .or. model%options%which_ho_beta_limit >= size(ho_powerlaw_c)) then
          call write_log('Error, HO powerlaw_c input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_coulomb_c            : ',model%options%which_ho_coulomb_c,  &
                         ho_coulomb_c(model%options%which_ho_coulomb_c)
       call write_log(message)
       if (model%options%which_ho_coulomb_c < 0 .or. model%options%which_ho_beta_limit >= size(ho_coulomb_c)) then
          call write_log('Error, HO coulomb_c input out of range', GM_FATAL)
       end if

       ! Inversion options

       ! Note: Inversion for Cp is currently supported for the Schoof sliding law, Tsai law, and basic power law
       if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION) then
          if (model%options%which_ho_babc == HO_BABC_COULOMB_POWERLAW_SCHOOF .or.  &
              model%options%which_ho_babc == HO_BABC_COULOMB_POWERLAW_TSAI .or.  &
              model%options%which_ho_babc == HO_BABC_POWERLAW) then
             ! inversion for Cp is supported
          else
             call write_log('Error, Cp inversion is not supported for this basal BC option')
             write(message,*) 'Cp inversion is supported only for these options: ', &
                  HO_BABC_COULOMB_POWERLAW_SCHOOF, HO_BABC_COULOMB_POWERLAW_TSAI, HO_BABC_POWERLAW
             call write_log(message, GM_FATAL)
          endif
       endif

       ! Note: Inversion for Cc is currently supported only for the Zoet-Iverson law
       if (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION) then
          if (model%options%which_ho_babc == HO_BABC_ZOET_IVERSON) then
             ! inversion for Cc is supported
          else
             call write_log('Error, Cc inversion is not supported for this basal BC option')
             write(message,*) 'Cc inversion is supported only for these options: ', &
                  HO_BABC_ZOET_IVERSON
             call write_log(message, GM_FATAL)
          endif
       endif

       if (model%options%which_ho_bmlt_basin_inversion /= HO_BMLT_BASIN_INVERSION_NONE) then
          write(message,*) 'ho_bmlt_basin_whichinversion : ',model%options%which_ho_bmlt_basin_inversion,  &
                            ho_bmlt_basin_whichinversion(model%options%which_ho_bmlt_basin_inversion)
          call write_log(message)
          if (model%options%whichbmlt_float /= BMLT_FLOAT_THERMAL_FORCING) then
             call write_log('Error, bmlt_basin inversion is not supported for this bmlt_float option')
             write(message,*) 'bmlt_basin inversion is supported only for bmlt_float = ', BMLT_FLOAT_THERMAL_FORCING
             call write_log(message, GM_FATAL)
          endif
       endif

       if (model%options%which_ho_bmlt_basin_inversion < 0 .or. &
            model%options%which_ho_bmlt_basin_inversion >= size(ho_bmlt_basin_whichinversion)) then
          call write_log('Error, bmlt_basin inversion input out of range', GM_FATAL)
       end if

       ! basal water options

       write(message,*) 'ho_whichbwat            : ',model%options%which_ho_bwat,  &
                         ho_whichbwat(model%options%which_ho_bwat)
       call write_log(message)
       if (model%options%which_ho_bwat < 0 .or. model%options%which_ho_bwat >= size(ho_whichbwat)) then
          call write_log('Error, HO basal water input out of range', GM_FATAL)
       end if

    if (model%options%which_ho_bwat == HO_BWAT_CONSTANT) then
       write(message,*) 'constant basal water depth (m): ', model%basal_hydro%const_bwat
       call write_log(message)
    elseif (model%options%which_ho_bwat == HO_BWAT_LOCAL_TILL) then
       write(message,*) 'maximum till water depth (m)  : ', model%basal_hydro%bwat_till_max
       call write_log(message)
       write(message,*) 'till drainage rate (m/yr)     : ', model%basal_hydro%c_drainage
       call write_log(message)
    elseif (model%options%which_ho_bwat == HO_BWAT_FLUX_ROUTING) then
       if (model%options%ho_flux_routing_scheme < 0.or. &
           model%options%ho_flux_routing_scheme >= size(ho_flux_routing_scheme)) then
          call write_log('Error, HO flux routing scheme out of range', GM_FATAL)
       end if
       write(message,*) 'ho_flux_routing_scheme  : ',model%options%ho_flux_routing_scheme,  &
            ho_flux_routing_scheme(model%options%ho_flux_routing_scheme)
       call write_log(message)
    endif

       write(message,*) 'ho_whicheffecpress      : ',model%options%which_ho_effecpress,  &
                         ho_whicheffecpress(model%options%which_ho_effecpress)
       call write_log(message)
       if (model%options%which_ho_effecpress < 0 .or. model%options%which_ho_effecpress >= size(ho_whicheffecpress)) then
          call write_log('Error, HO effective pressure input out of range', GM_FATAL)
       end if

       write(message,*) 'which_ho_nonlinear      : ',model%options%which_ho_nonlinear,  &
                         which_ho_nonlinear(model%options%which_ho_nonlinear)
       call write_log(message)
       if (model%options%which_ho_nonlinear < 0 .or. model%options%which_ho_nonlinear >= size(which_ho_nonlinear)) then
          call write_log('Error, HO nonlinear solution input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichresid           : ',model%options%which_ho_resid,  &
                         ho_whichresid(model%options%which_ho_resid)
       call write_log(message)
       if (model%options%which_ho_resid < 0 .or. model%options%which_ho_resid >= size(ho_whichresid)) then
          call write_log('Error, HO residual input out of range', GM_FATAL)
       end if
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

       write(message,*) 'ho_whichsparse          : ',model%options%which_ho_sparse,  &
                         ho_whichsparse(model%options%which_ho_sparse)
       call write_log(message)
       if (model%options%which_ho_sparse < -1 .or. model%options%which_ho_sparse >= size(ho_whichsparse)) then
          call write_log('Error, HO sparse solver input out of range', GM_FATAL)
       end if

       if (model%options%whichdycore == DYCORE_GLISSADE) then

          write(message,*) 'ho_whichapprox          : ',model%options%which_ho_approx,  &
                            ho_whichapprox(model%options%which_ho_approx)
          call write_log(message)
          if (model%options%which_ho_approx < -1 .or. model%options%which_ho_approx >= size(ho_whichapprox)-1) then
             call write_log('Error, Stokes approximation out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichgradient        : ',model%options%which_ho_gradient,  &
                            ho_whichgradient(model%options%which_ho_gradient)
          call write_log(message)
          if (model%options%which_ho_gradient < 0 .or. model%options%which_ho_gradient >= size(ho_whichgradient)) then
             call write_log('Error, gradient option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichgradient_margin : ',model%options%which_ho_gradient_margin,  &
                            ho_whichgradient_margin(model%options%which_ho_gradient_margin)
          call write_log(message)
          if (model%options%which_ho_gradient_margin < 0 .or. &
              model%options%which_ho_gradient_margin >= size(ho_whichgradient_margin)) then
             call write_log('Error, gradient margin option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichvertical_remap  : ',model%options%which_ho_vertical_remap,  &
                            ho_whichvertical_remap(model%options%which_ho_vertical_remap)
          call write_log(message)
          if (model%options%which_ho_vertical_remap < 0 .or. &
              model%options%which_ho_vertical_remap >= size(ho_whichvertical_remap)) then
             call write_log('Error, vertical remap option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichassemble_beta   : ',model%options%which_ho_assemble_beta,  &
                            ho_whichassemble_beta(model%options%which_ho_assemble_beta)
          call write_log(message)
          if (model%options%which_ho_assemble_beta < 0 .or. &
              model%options%which_ho_assemble_beta >= size(ho_whichassemble_beta)) then
             call write_log('Error, beta assembly option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichassemble_taud   : ',model%options%which_ho_assemble_taud,  &
                            ho_whichassemble_taud(model%options%which_ho_assemble_taud)
          call write_log(message)
          if (model%options%which_ho_assemble_taud < 0 .or. &
              model%options%which_ho_assemble_taud >= size(ho_whichassemble_taud)) then
             call write_log('Error, driving-stress assembly option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichassemble_bfric  : ',model%options%which_ho_assemble_bfric,  &
                            ho_whichassemble_bfric(model%options%which_ho_assemble_bfric)
          call write_log(message)
          if (model%options%which_ho_assemble_bfric < 0 .or. &
              model%options%which_ho_assemble_bfric >= size(ho_whichassemble_bfric)) then
             call write_log('Error, basal-friction assembly option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichassemble_lateral  : ',model%options%which_ho_assemble_lateral,  &
                            ho_whichassemble_lateral(model%options%which_ho_assemble_lateral)
          call write_log(message)
          if (model%options%which_ho_assemble_lateral < 0 .or. &
              model%options%which_ho_assemble_lateral >= size(ho_whichassemble_lateral)) then
             call write_log('Error, lateral-stress assembly option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichcalving_front   : ',model%options%which_ho_calving_front,  &
                            ho_whichcalving_front(model%options%which_ho_calving_front)
          call write_log(message)
          if (model%options%which_ho_calving_front < 0 .or. &
               model%options%which_ho_calving_front >= size(ho_whichcalving_front)) then
             call write_log('Error, calving front option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichground          : ',model%options%which_ho_ground,  &
                            ho_whichground(model%options%which_ho_ground)
          call write_log(message)
          if (model%options%which_ho_ground < 0 .or. model%options%which_ho_ground >= size(ho_whichground)) then
             call write_log('Error, ho_ground option out of range for glissade dycore', GM_FATAL)
          end if

          if (model%options%which_ho_ground == HO_GROUND_NO_GLP) then
             write(message,*) 'ho_whichfground_no_glp  : ', model%options%which_ho_fground_no_glp,  &
                               ho_whichfground_no_glp(model%options%which_ho_fground_no_glp)
             call write_log(message)
             if (model%options%which_ho_fground_no_glp < 0 .or. &
                 model%options%which_ho_fground_no_glp >= size(ho_whichfground_no_glp)) then
                call write_log('Error, ho_fground_no_glp option out of range for glissade dycore', GM_FATAL)
             end if
          endif

          write(message,*) 'ho_whichground_bmlt     : ',model%options%which_ho_ground_bmlt,  &
                            ho_whichground_bmlt(model%options%which_ho_ground_bmlt)
          call write_log(message)
          if (model%options%which_ho_ground_bmlt < 0 .or. &
              model%options%which_ho_ground_bmlt >= size(ho_whichground_bmlt)) then
             call write_log('Error, ho_ground_bmlt option out of range for glissade dycore', GM_FATAL)
          end if

          if (model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_FLOATING_FRAC  .or. &
              model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_ZERO_GROUNDED) then
             if (model%options%which_ho_ground /= HO_GROUND_GLP_DELUXE) then
                write(message,*) 'This setting of which_ho_ground_bmlt requires which_ho_ground =', &
                     HO_GROUND_GLP_DELUXE
                call write_log(message)
                write(message,*)  'Defaulting to which_ho_ground_bmlt =', HO_GROUND_BMLT_NO_GLP
                call write_log(message)
                model%options%which_ho_ground_bmlt = HO_GROUND_BMLT_NO_GLP
             endif
          endif

          write(message,*) 'ho_whichflotation_function:',model%options%which_ho_flotation_function,  &
               ho_whichflotation_function(model%options%which_ho_flotation_function)
          call write_log(message)
          if (model%options%which_ho_flotation_function < 0 .or. &
               model%options%which_ho_flotation_function >= size(ho_whichflotation_function)) then
             call write_log('Error, flotation_function option out of range for glissade dycore', GM_FATAL)
          endif

          if (model%options%block_inception) then
             write(message,*) 'Inception outside the main ice sheet will be blocked'
             call write_log(message)
          endif

          if (model%options%remove_ice_caps) then
             write(message,*) 'Ice caps will be removed and added to the calving flux'
             call write_log(message)
          endif

          if (model%options%force_retreat == FORCE_RETREAT_ALL_ICE) then
             write(message,*) 'Ice retreat will be forced using ice_fraction_retreat_mask'
             call write_log(message)
          elseif (model%options%force_retreat == FORCE_RETREAT_FLOATING_ICE) then
             write(message,*) 'Floating ice retreat will be forced using ice_fraction_retreat_mask'
             call write_log(message)
             if (.not.model%options%remove_isthmuses) then
                call write_log('  Warning: Can be unstable when remove_isthmuses = F')
             endif
          endif

          write(message,*) 'ho_whichice_age         : ',model%options%which_ho_ice_age,  &
                            ho_whichice_age(model%options%which_ho_ice_age)
          call write_log(message)
          if (model%options%which_ho_ice_age < 0 .or. model%options%which_ho_ice_age >= size(ho_whichice_age)) then
             call write_log('Error, ice_age option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'glissade_maxiter        : ',model%options%glissade_maxiter
          call write_log(message)

          write(message,*) 'linear_solve_ncheck     : ',model%options%linear_solve_ncheck
          call write_log(message)

          write(message,*) 'linear_maxiters         : ',model%options%linear_maxiters
          call write_log(message)

          write(message,*) 'linear_tolerance        : ',model%options%linear_tolerance
          call write_log(message)

       end if   ! DYCORE_GLISSADE

       if (model%options%whichdycore == DYCORE_GLISSADE .and.   &
           model%options%which_ho_ground == HO_GROUND_NO_GLP .and. &
           model%options%which_ho_flotation_function == HO_FLOTATION_FUNCTION_PATTYN) then
          write(message,*) 'WARNING: Pattyn flotation function with no GLP tends to be unstable'
          call write_log(message, GM_WARNING)
       endif

       if (model%options%whichdycore == DYCORE_GLISSADE .and.   &
           (model%options%which_ho_sparse == HO_SPARSE_PCG_STANDARD .or.  &
            model%options%which_ho_sparse == HO_SPARSE_PCG_CHRONGEAR) ) then 
          write(message,*) 'ho_whichprecond         : ',model%options%which_ho_precond,  &
                            ho_whichprecond(model%options%which_ho_precond)
          call write_log(message)
          if (model%options%which_ho_precond < 0 .or. model%options%which_ho_precond >= size(ho_whichprecond)) then
             call write_log('Error, glissade preconditioner out of range', GM_FATAL)
          end if
       end if

    endif   ! whichdycore

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

    !TODO - Reorganize parameters into sections based on relevant physics
    !Note: The following physical constants have default values in glimmer_physcon.F90.
    !      Some test cases (e.g., MISMIP) specify different values. The default values
    !       can therefore be overridden by the user in the config file.
    !      For coupled CESM runs, however, CISM uses the values in CESM's shr_const_mod,
    !       which cannot be overridden.
    !      These constants are *not* part of the model derived type.
    !      If running multiple instances, the user should either use the default values
    !       or specify identical values in each config file.  Otherwise, the run will use
    !       whatever values are specified in the last config file to be read.
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
    call GetValue(section,'ice_limit_temp',     model%numerics%thklim_temp)
    call GetValue(section,'thck_gradient_ramp', model%numerics%thck_gradient_ramp)
    call GetValue(section,'pmp_offset',         model%temper%pmp_offset)
    call GetValue(section,'pmp_threshold',      model%temper%pmp_threshold)
    call GetValue(section,'geothermal',         model%paramets%geot)
    call GetValue(section,'flow_factor',        model%paramets%flow_enhancement_factor)
    call GetValue(section,'flow_factor_float',  model%paramets%flow_enhancement_factor_float)
    !TODO - Change default_flwa to flwa_constant?  Would have to change config files.
    call GetValue(section,'default_flwa',       model%paramets%default_flwa)
    call GetValue(section,'efvs_constant',      model%paramets%efvs_constant)
    call GetValue(section,'effstrain_min',      model%paramets%effstrain_min)
    call GetValue(section,'hydro_time',         model%paramets%hydtim)
    call GetValue(section,'max_slope',          model%paramets%max_slope)

    ! parameters to adjust external forcing
    call GetValue(section,'t_lapse',            model%climate%t_lapse)
    call GetValue(section,'acab_factor',        model%climate%acab_factor)
    call GetValue(section,'bmlt_float_factor',  model%basal_melt%bmlt_float_factor)

    ! calving parameters
    call GetValue(section,'marine_limit',       model%calving%marine_limit)
    call GetValue(section,'calving_fraction',   model%calving%calving_fraction)
    call GetValue(section,'calving_minthck',    model%calving%minthck)
    call GetValue(section,'lateral_rate_max',   model%calving%lateral_rate_max)
    call GetValue(section,'eigencalving_constant', model%calving%eigencalving_constant)
    call GetValue(section,'eigen2_weight',      model%calving%eigen2_weight)
    call GetValue(section,'damage_constant',    model%calving%damage_constant)
    call GetValue(section,'taumax_cliff',       model%calving%taumax_cliff)
    call GetValue(section,'cliff_timescale',    model%calving%cliff_timescale)
    call GetValue(section,'ncull_calving_front',   model%calving%ncull_calving_front)
    call GetValue(section,'calving_timescale',  model%calving%timescale)
    call GetValue(section,'calving_front_x',    model%calving%calving_front_x)
    call GetValue(section,'calving_front_y',    model%calving%calving_front_y)
    call GetValue(section,'damage_threshold',   model%calving%damage_threshold)

    ! NOTE: bpar is used only for BTRC_TANH_BWAT
    !       btrac_max and btrac_slope are used (with btrac_const) for BTRC_LINEAR_BMLT
    !       btrac_const is used for several options

    call GetValue(section,'basal_tract_const', model%paramets%btrac_const)
    call GetValue(section,'basal_tract_max',   model%paramets%btrac_max)
    call GetValue(section,'basal_tract_slope', model%paramets%btrac_slope)

    !WHL - Changed this so that bpar can be read correctly from config file.
    !      This parameter is now called 'basal_tract_tanh' instead of 'basal_tract'.
    call GetValue(section,'basal_tract_tanh',  tempvar, 5)
    if (associated(tempvar)) then
!!       model%paramets%btrac_const = tempvar(1)  ! old code
       model%paramets%bpar(:) = tempvar(:)
       deallocate(tempvar)
    end if

    call GetValue(section,'beta_grounded_min', model%basal_physics%beta_grounded_min)
    call GetValue(section,'ho_beta_const', model%basal_physics%ho_beta_const)
    call GetValue(section,'ho_beta_small', model%basal_physics%ho_beta_small)
    call GetValue(section,'ho_beta_large', model%basal_physics%ho_beta_large)

    ! basal friction parameters
    call GetValue(section, 'powerlaw_c_const', model%basal_physics%powerlaw_c_const)
    call GetValue(section, 'powerlaw_c_max', model%basal_physics%powerlaw_c_max)
    call GetValue(section, 'powerlaw_c_min', model%basal_physics%powerlaw_c_min)
    call GetValue(section, 'powerlaw_m', model%basal_physics%powerlaw_m)
    call GetValue(section, 'coulomb_c_const', model%basal_physics%coulomb_c_const)
    call GetValue(section, 'coulomb_c_max', model%basal_physics%coulomb_c_max)
    call GetValue(section, 'coulomb_c_min', model%basal_physics%coulomb_c_min)
    call GetValue(section, 'coulomb_c_bedmax', model%basal_physics%coulomb_c_bedmax)
    call GetValue(section, 'coulomb_c_bedmin', model%basal_physics%coulomb_c_bedmin)
    call GetValue(section, 'beta_powerlaw_umax', model%basal_physics%beta_powerlaw_umax)
    call GetValue(section, 'zoet_iversion_ut', model%basal_physics%zoet_iverson_ut)
    call GetValue(section, 'zoet_iversion_nmax', model%basal_physics%zoet_iverson_nmax)
    call GetValue(section, 'friction_powerlaw_k', model%basal_physics%friction_powerlaw_k)
    call GetValue(section, 'flwa_basal', model%basal_physics%flwa_basal)
    call GetValue(section, 'coulomb_bump_max_slope', model%basal_physics%coulomb_bump_max_slope)
    call GetValue(section, 'coulomb_bump_wavelength', model%basal_physics%coulomb_bump_wavelength)

    ! effective pressure parameters
    call GetValue(section, 'p_ocean_penetration', model%basal_physics%p_ocean_penetration)
    call GetValue(section, 'ocean_p_timescale', model%basal_physics%ocean_p_timescale)
    call GetValue(section, 'effecpress_delta', model%basal_physics%effecpress_delta)
    call GetValue(section, 'effecpress_bpmp_threshold', model%basal_physics%effecpress_bpmp_threshold)
    call GetValue(section, 'effecpress_bwat_threshold', model%basal_physics%effecpress_bwat_threshold)
    call GetValue(section, 'effecpress_bwatflx_threshold', model%basal_physics%effecpress_bwatflx_threshold)
    call GetValue(section, 'effecpress_timescale', model%basal_physics%effecpress_timescale)

    ! basal water parameters
    call GetValue(section, 'const_bwat', model%basal_hydro%const_bwat)
    call GetValue(section, 'bwat_till_max', model%basal_hydro%bwat_till_max)
    call GetValue(section, 'c_drainage', model%basal_hydro%c_drainage)

    ! pseudo-plastic parameters
    call GetValue(section, 'pseudo_plastic_q', model%basal_physics%pseudo_plastic_q)
    call GetValue(section, 'pseudo_plastic_u0', model%basal_physics%pseudo_plastic_u0)
    !TODO - next four to be removed in favor of coulomb_c_min, etc.
    call GetValue(section, 'pseudo_plastic_phimin', model%basal_physics%pseudo_plastic_phimin)
    call GetValue(section, 'pseudo_plastic_phimax', model%basal_physics%pseudo_plastic_phimax)
    call GetValue(section, 'pseudo_plastic_bedmin', model%basal_physics%pseudo_plastic_bedmin)
    call GetValue(section, 'pseudo_plastic_bedmax', model%basal_physics%pseudo_plastic_bedmax)

    ! ocean data parameters
    call GetValue(section, 'gamma0', model%ocean_data%gamma0)
    call GetValue(section, 'thermal_forcing_anomaly', model%ocean_data%thermal_forcing_anomaly)
    call GetValue(section, 'thermal_forcing_anomaly_tstart', model%ocean_data%thermal_forcing_anomaly_tstart)
    call GetValue(section, 'thermal_forcing_anomaly_timescale', model%ocean_data%thermal_forcing_anomaly_timescale)
    call GetValue(section, 'thermal_forcing_anomaly_basin', model%ocean_data%thermal_forcing_anomaly_basin)

    ! parameters to adjust input topography
    call GetValue(section, 'adjust_topg_xmin', model%paramets%adjust_topg_xmin)
    call GetValue(section, 'adjust_topg_xmax', model%paramets%adjust_topg_xmax)
    call GetValue(section, 'adjust_topg_ymin', model%paramets%adjust_topg_ymin)
    call GetValue(section, 'adjust_topg_ymax', model%paramets%adjust_topg_ymax)
    call GetValue(section, 'adjust_topg_no_adjust',  model%paramets%adjust_topg_no_adjust)
    call GetValue(section, 'adjust_topg_max_adjust', model%paramets%adjust_topg_max_adjust)
    call GetValue(section, 'adjust_topg_delta',   model%paramets%adjust_topg_delta)

    ! basal inversion parameters
    !TODO - Put inversion parameters in a separate section
    call GetValue(section, 'inversion_thck_flotation_buffer', model%inversion%thck_flotation_buffer)
    call GetValue(section, 'inversion_thck_threshold', model%inversion%thck_threshold)

    call GetValue(section, 'inversion_babc_timescale', model%inversion%babc_timescale)
    call GetValue(section, 'inversion_babc_thck_scale', model%inversion%babc_thck_scale)
    call GetValue(section, 'inversion_babc_velo_scale', model%inversion%babc_velo_scale)

    call GetValue(section, 'inversion_dbmlt_dtemp_scale', model%inversion%dbmlt_dtemp_scale)
    call GetValue(section, 'inversion_bmlt_basin_timescale', model%inversion%bmlt_basin_timescale)
    call GetValue(section, 'inversion_bmlt_basin_flotation_threshold', &
         model%inversion%bmlt_basin_flotation_threshold)
    call GetValue(section, 'inversion_bmlt_basin_mass_correction', &
         model%inversion%bmlt_basin_mass_correction)
    call GetValue(section, 'inversion_bmlt_basin_number_mass_correction', &
         model%inversion%bmlt_basin_number_mass_correction)

    ! ISMIP-HOM parameters
    call GetValue(section,'periodic_offset_ew',model%numerics%periodic_offset_ew)
    call GetValue(section,'periodic_offset_ns',model%numerics%periodic_offset_ns)

    ! parameters for acab/artm anomaly and overwrite options
    call GetValue(section,'acab_anomaly_timescale', model%climate%acab_anomaly_timescale)
    call GetValue(section,'overwrite_acab_value', model%climate%overwrite_acab_value)
    call GetValue(section,'overwrite_acab_minthck', model%climate%overwrite_acab_minthck)
    call GetValue(section,'bmlt_anomaly_timescale', model%basal_melt%bmlt_anomaly_timescale)

    ! parameters for artm anomaly option
    call GetValue(section,'artm_anomaly_timescale', model%climate%artm_anomaly_timescale)

    ! basal melting parameters
    call GetValue(section,'bmlt_cavity_h0', model%basal_melt%bmlt_cavity_h0)

    ! MISMIP+ basal melting parameters
    call GetValue(section,'bmlt_float_omega', model%basal_melt%bmlt_float_omega)
    call GetValue(section,'bmlt_float_h0', model%basal_melt%bmlt_float_h0)
    call GetValue(section,'bmlt_float_z0', model%basal_melt%bmlt_float_z0)
    call GetValue(section,'bmlt_float_const', model%basal_melt%bmlt_float_const)
    call GetValue(section,'bmlt_float_xlim', model%basal_melt%bmlt_float_xlim)

    ! depth-dependent basal melting parameters
    call GetValue(section,'bmlt_float_depth_meltmax', model%basal_melt%bmlt_float_depth_meltmax)
    call GetValue(section,'bmlt_float_depth_frzmax', model%basal_melt%bmlt_float_depth_frzmax)
    call GetValue(section,'bmlt_float_depth_zmeltmax', model%basal_melt%bmlt_float_depth_zmeltmax)
    call GetValue(section,'bmlt_float_depth_zmelt0', model%basal_melt%bmlt_float_depth_zmelt0)
    call GetValue(section,'bmlt_float_depth_zfrzmax', model%basal_melt%bmlt_float_depth_zfrzmax)
    call GetValue(section,'bmlt_float_depth_meltmin', model%basal_melt%bmlt_float_depth_meltmin)
    call GetValue(section,'bmlt_float_depth_zmeltmin', model%basal_melt%bmlt_float_depth_zmeltmin)

    ! MISOMIP plume parameters
    !TODO - Put MISMIP+ and MISOMIP parameters in their own section
    call GetValue(section,'T0',   model%plume%T0)
    call GetValue(section,'Tbot', model%plume%Tbot)
    call GetValue(section,'S0',   model%plume%S0)
    call GetValue(section,'Sbot', model%plume%Sbot)
    call GetValue(section,'gammaT',    model%plume%gammaT)
    call GetValue(section,'gammaS',    model%plume%gammaS)

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

    !Note: The Glissade dycore is known to crash for thklim = 0, but has not
    !      been extensively tested for small values of thklim.
    !      Values smaller than 1 mm may be OK, but no guarantees.
    if (model%options%whichdycore == DYCORE_GLISSADE .and.   &
        model%numerics%thklim < 1.d-3) then   ! 1 mm
       call write_log('ice limit (thklim) is too small for Glissade dycore', GM_FATAL)
    endif

    if (model%options%whichdycore == DYCORE_GLISSADE .and.   &
         model%numerics%adaptive_cfl_threshold > 0.0d0) then
       write(message,*) 'Advection will be subcycled when CFL >', model%numerics%adaptive_cfl_threshold
       call write_log(message)
       if (model%numerics%subcyc /= 1) then  ! set subcyc = 1 with adaptive CFL subcycling
          model%numerics%subcyc = 1
          write(message,*) 'Setting numerics%subcyc = 1; subcycling is adaptive only'
          call write_log(message)
       endif
    endif

    if (model%options%whichdycore /= DYCORE_GLIDE) then
       write(message,*) 'thickness limit for temperature calculations (m) : ', model%numerics%thklim_temp
       call write_log(message)
       if (model%numerics%thck_gradient_ramp > 0.0d0) then
          write(message,*) 'thickness scale for gradient ramp (m):', model%numerics%thck_gradient_ramp
          call write_log(message)
       endif
       write(message,*) 'pmp threshold for temperature (K): ', model%temper%pmp_threshold
       call write_log(message)
    endif

    if (model%climate%acab_factor /= 1.0d0) then
       write(message,*) 'Input acab multiplied by      :', model%climate%acab_factor
       call write_log(message)
    endif

    if (model%options%whichcalving == CALVING_FLOAT_FRACTION) then
       write(message,*) 'ice fraction lost in calving  : ', model%calving%calving_fraction
       call write_log(message)
    end if

    if (model%options%whichcalving == CALVING_RELX_THRESHOLD .or.  &
        model%options%whichcalving == CALVING_TOPG_THRESHOLD) then
       write(message,*) 'marine depth limit (m)        : ', model%calving%marine_limit
       call write_log(message)
    endif

    ! thickness-based calving options
    if (model%options%whichcalving == CALVING_THCK_THRESHOLD .or. &
        model%options%whichcalving == EIGENCALVING           .or. &
        model%options%whichcalving == CALVING_DAMAGE) then

       if (model%calving%timescale <= 0.0d0) then
          write(message,*) 'Must set calving_timescale to a positive nonzero value for this calving option'
          call write_log(message, GM_FATAL)
       endif

       if (model%options%whichcalving == EIGENCALVING .or. &
           model%options%whichcalving == CALVING_DAMAGE) then
          if (model%options%which_ho_calving_front == HO_CALVING_FRONT_NO_SUBGRID) then
             write(message,*) &
                  'Calving option ', model%options%whichcalving, ' requires a subgrid calving front'
             call write_log(message, GM_FATAL)
          endif
       endif

       if (model%options%whichcalving == CALVING_THCK_THRESHOLD) then
          if (model%calving%minthck > 0.0d0) then
             write(message,*) 'calving thickness threshold (m) : ', model%calving%minthck
             call write_log(message)
          else
             write(message,*) 'Will use a 2D calving thickness threshold field'
             call write_log(message)
          endif
       elseif (model%options%whichcalving == EIGENCALVING) then
          if (model%calving%minthck == 0.0d0) then
             write(message,*) 'Error: Eigencalving requires minthck > 0'
             call write_log(message, GM_FATAL)
          else
             write(message,*) 'calving thickness threshold (m) : ', model%calving%minthck
             call write_log(message)
          endif
          write(message,*) 'eigencalving constant (m yr^-1 Pa^-1): ', model%calving%eigencalving_constant
          call write_log(message)
          write(message,*) 'eigenvalue 2 weight (unitless)       : ', model%calving%eigen2_weight
          call write_log(message)
       elseif (model%options%whichcalving == CALVING_DAMAGE) then
          if (model%calving%minthck == 0.0d0) then
             write(message,*) 'Error: Damage-based calving requires minthck > 0'
             call write_log(message, GM_FATAL)
          else
             write(message,*) 'calving thickness threshold (m) : ', model%calving%minthck
             call write_log(message)
          endif
          write(message,*) 'damage constant (yr^-1)              : ', model%calving%damage_constant
          call write_log(message)
          write(message,*) 'damage threshold                     : ', model%calving%damage_threshold
          call write_log(message)
          write(message,*) 'max lateral calving rate (m/yr)      : ', model%calving%lateral_rate_max
          call write_log(message)
       endif
    endif   ! CALVING_THCK_THRESHOLD, EIGENCALVING, CALVING_DAMAGE

    if (model%options%which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then
       if (.not.model%options%remove_icebergs) then
          model%options%remove_icebergs = .true.
          write(message,*) 'Setting remove_icebergs = T for stability when using subgrid calving_front scheme'
          call write_log(message)
       endif
    endif

    if (model%options%limit_marine_cliffs) then
       write(message,*) 'taumax_cliff                  : ', model%calving%taumax_cliff
       call write_log(message)
       write(message,*) 'cliff time scale (yr)       : ', model%calving%cliff_timescale
       call write_log(message)
    endif

    if (model%options%whichcalving == CALVING_GRID_MASK) then
       if (model%calving%calving_front_x > 0.0d0) then
          write(message,*) 'x calving front (m)           : ', model%calving%calving_front_x
          call write_log(message)
       endif
       if (model%calving%calving_front_y > 0.0d0) then
          write(message,*) 'y calving front (m)           : ', model%calving%calving_front_y
          call write_log(message)
       endif
    endif

    if (model%calving%timescale > 0.0d0) then
       write(message,*) 'calving time scale (yr)       : ', model%calving%timescale
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

    write(message,*) 'flow factor (grounded ice)    : ', model%paramets%flow_enhancement_factor
    call write_log(message)

    write(message,*) 'flow factor (floating ice)    : ', model%paramets%flow_enhancement_factor_float
    call write_log(message)

    if (model%options%whichdycore == DYCORE_GLIDE) then
       write(message,*) 'basal hydro time constant (yr): ', model%paramets%hydtim
       call write_log(message)
    endif

    if (model%options%whichdycore == DYCORE_GLISSADE) then
       write(message,*) 'max surface slope             : ', model%paramets%max_slope
       call write_log(message)
    end if       
 
    if (model%options%whichflwa == FLWA_CONST_FLWA) then
       write(message,*) 'constant flow factor (Pa^-n yr^-1) :', model%paramets%default_flwa
       call write_log(message)
    end if

    if (model%options%which_ho_efvs == HO_EFVS_CONSTANT) then
       write(message,*) 'constant effec viscosity (Pa yr)   :', model%paramets%efvs_constant
       call write_log(message)
    elseif (model%options%which_ho_efvs == HO_EFVS_NONLINEAR) then
       write(message,*) 'min effective strain rate (yr^-1)  :', model%paramets%effstrain_min
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

    if (model%options%which_ho_babc == HO_BABC_BETA_CONSTANT) then
       write(message,*) 'uniform beta (Pa yr/m)        : ',model%basal_physics%ho_beta_const
       call write_log(message)
    elseif (model%options%which_ho_babc == HO_BABC_BETA_BPMP) then
       write(message,*) 'large (frozen) beta (Pa yr/m) : ',model%basal_physics%ho_beta_large
       call write_log(message)
       write(message,*) 'small (thawed) beta (Pa yr/m) : ',model%basal_physics%ho_beta_small
       call write_log(message)
    elseif (model%options%which_ho_babc == HO_BABC_PSEUDO_PLASTIC_OLD .or.  &
            model%options%which_ho_babc == HO_BABC_PSEUDO_PLASTIC) then
       write(message,*) 'pseudo-plastic q              : ',model%basal_physics%pseudo_plastic_q
       call write_log(message)
       write(message,*) 'pseudo-plastic u0             : ',model%basal_physics%pseudo_plastic_u0
       call write_log(message)
       if (model%options%which_ho_babc == HO_BABC_PSEUDO_PLASTIC_OLD) then
          write(message,*) 'pseudo-plastic phi_min (deg)  : ',model%basal_physics%pseudo_plastic_phimin
          call write_log(message)
          write(message,*) 'pseudo-plastic phi_max (deg)  : ',model%basal_physics%pseudo_plastic_phimax
          call write_log(message)
          write(message,*) 'pseudo-plastic bed min (m)    : ',model%basal_physics%pseudo_plastic_bedmin
          call write_log(message)
          write(message,*) 'pseudo-plastic bed max (m)    : ',model%basal_physics%pseudo_plastic_bedmax
          call write_log(message)
       endif
       ! Note: For the new Coulomb_C elevation option, phimin/phimax/bedmin/bedmax are written below.
       if (model%options%which_ho_assemble_beta == HO_ASSEMBLE_BETA_STANDARD) then
          call write_log('WARNING: local beta assembly is recommended for the pseudo-plastic sliding law')
          write(message,*) 'Set which_ho_assemble_beta =', HO_ASSEMBLE_BETA_LOCAL
          call write_log(message, GM_WARNING)
       endif
       if (model%options%which_ho_assemble_bfric == HO_ASSEMBLE_BFRIC_STANDARD) then
          call write_log('WARNING: local bfric assembly is recommended for the pseudo-plastic sliding law')
          write(message,*) 'Set which_ho_assemble_bfric =', HO_ASSEMBLE_BFRIC_LOCAL
          call write_log(message, GM_WARNING)
       endif
       if (model%options%which_ho_thermal_timestep == HO_THERMAL_BEFORE_TRANSPORT) then
          call write_log('WARNING: Best to do vertical thermal solve after transport for the pseudo-plastic sliding law')
          write(message,*) 'Set which_ho_thermal_timestep to one of the following values:', &
               HO_THERMAL_AFTER_TRANSPORT, HO_THERMAL_SPLIT_TIMESTEP
          call write_log(message, GM_WARNING)
       endif
    elseif (model%options%which_ho_babc == HO_BABC_ZOET_IVERSON) then
       ! Note: The Zoet-Iverson law typically uses a spatially variable coulomb_c.
       !       If so, the value written here is just the initial value.
       write(message,*) 'Cc for Zoet-Iversion law                     : ', model%basal_physics%coulomb_c_const
       call write_log(message)
       write(message,*) 'm exponent for Zoet-Iverson law              : ', model%basal_physics%powerlaw_m
       call write_log(message)
       write(message,*) 'threshold speed for Zoet-Iverson law (m/yr)  : ', model%basal_physics%zoet_iverson_ut
       call write_log(message)
       write(message,*) 'max effecpress for Zoet-Iverson law (Pa)     : ', model%basal_physics%zoet_iverson_nmax
       call write_log(message)
    elseif (model%options%which_ho_babc == HO_BABC_ISHOMC) then
       if (model%general%ewn /= model%general%nsn) then
          call write_log('Error, must have ewn = nsn for ISMIP-HOM test C', GM_FATAL)
       endif
    elseif (model%options%which_ho_babc == HO_BABC_POWERLAW) then
       write(message,*) 'Cp for power law, Pa (m/yr)^(-1/3)           : ', model%basal_physics%powerlaw_c_const
       call write_log(message)
       write(message,*) 'm exponent for power law                     : ', model%basal_physics%powerlaw_m
       call write_log(message)
    elseif (model%options%which_ho_babc == HO_BABC_COULOMB_FRICTION) then
       write(message,*) 'Cc for Coulomb friction law                  : ', model%basal_physics%coulomb_c_const
       call write_log(message)
       write(message,*) 'bed bump max slope for Coulomb friction law  : ', model%basal_physics%coulomb_bump_max_slope
       call write_log(message)
       write(message,*) 'bed bump wavelength for Coulomb friction law : ', model%basal_physics%coulomb_bump_wavelength
       call write_log(message)
    elseif (model%options%which_ho_babc == HO_BABC_COULOMB_POWERLAW_SCHOOF) then
       ! Note: The Schoof law typically uses a spatially variable powerlaw_c.
       !       If so, the value written here is just the initial value.
       write(message,*) 'Cc for Schoof Coulomb law                    : ', model%basal_physics%coulomb_c_const
       call write_log(message)
       write(message,*) 'Cp for Schoof power law, Pa (m/yr)^(-1/3)    : ', model%basal_physics%powerlaw_c_const
       call write_log(message)
       write(message,*) 'm exponent for Schoof power law              : ', model%basal_physics%powerlaw_m
       call write_log(message)
    elseif (model%options%which_ho_babc == HO_BABC_COULOMB_POWERLAW_TSAI) then
       ! Note: The Tsai law typically uses a spatially variable powerlaw_c. 
       !       If so, the value written here is just the initial value.
       write(message,*) 'Cc for Tsai Coulomb law                      : ', model%basal_physics%coulomb_c_const
       call write_log(message)
       write(message,*) 'Cp for Tsai power law, Pa (m/yr)^(-1/3)      : ', model%basal_physics%powerlaw_c_const
       call write_log(message)
       write(message,*) 'm exponent for Tsai power law                : ', model%basal_physics%powerlaw_m
       call write_log(message)
    elseif (model%options%which_ho_babc == HO_BABC_POWERLAW_EFFECPRESS) then
       call write_log('Weertman-style power law higher-order basal boundary condition is not currently scientifically &
            &supported.  USE AT YOUR OWN RISK.', GM_WARNING)
       !TODO - Use powerlaw_c instead of friction_powerlaw_k?  Allow p and q to be set in config file instead of hard-wired?
       write(message,*) 'roughness parameter, k, for power-law friction law : ',model%basal_physics%friction_powerlaw_k
       call write_log(message)
    endif

    ! Coulomb elevation parameters
    if (model%options%which_ho_coulomb_c == HO_COULOMB_C_ELEVATION) then
       write(message,*) 'coulomb_c_max                                : ',model%basal_physics%coulomb_c_max
       call write_log(message)
       write(message,*) 'coulomb_c_min                                : ',model%basal_physics%coulomb_c_min
       call write_log(message)
       write(message,*) 'coulomb_c_bedmax (m)                         : ',model%basal_physics%coulomb_c_bedmax
       call write_log(message)
       write(message,*) 'coulomb_c_bedmin (m)                         : ',model%basal_physics%coulomb_c_bedmin
       call write_log(message)
    endif

    if (model%options%adjust_input_topography) then
       call write_log('Input topography will be adjusted')
       write(message,*) 'adjust_topg_xmin (m)                         : ', model%paramets%adjust_topg_xmin
       call write_log(message)
       write(message,*) 'adjust_topg_xmax (m)                         : ', model%paramets%adjust_topg_xmax
       call write_log(message)
       write(message,*) 'adjust_topg_ymin (m)                         : ', model%paramets%adjust_topg_ymin
       call write_log(message)
       write(message,*) 'adjust_topg_ymax (m)                         : ', model%paramets%adjust_topg_ymax
       call write_log(message)
       write(message,*) 'adjust_topg_no_adjust (m)                    : ', model%paramets%adjust_topg_no_adjust
       call write_log(message)
       write(message,*) 'adjust_topg_max_adjust (m)                   : ', model%paramets%adjust_topg_max_adjust
       call write_log(message)
       write(message,*) 'adjust_topg_delta (m)                        : ', model%paramets%adjust_topg_delta
       call write_log(message)
    endif

    ! inversion parameters

    if (model%options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION) then
       write(message,*) 'inversion flotation thickness buffer (m)     : ', &
            model%inversion%thck_flotation_buffer
       call write_log(message)
       write(message,*) 'inversion thickness threshold (m)            : ', &
            model%inversion%thck_threshold
       call write_log(message)
       write(message,*) 'powerlaw_c max, Pa (m/yr)^(-1/3)             : ', &
            model%basal_physics%powerlaw_c_max
       call write_log(message)
       write(message,*) 'powerlaw_c min, Pa (m/yr)^(-1/3)             : ', &
            model%basal_physics%powerlaw_c_min
       call write_log(message)
       write(message,*) 'inversion basal friction timescale (yr)      : ', &
            model%inversion%babc_timescale
       call write_log(message)
       if (model%inversion%babc_thck_scale > 0.0d0) then
          write(message,*) 'inversion thickness scale (m)                : ', &
               model%inversion%babc_thck_scale
          call write_log(message)
       endif
       if (model%inversion%babc_velo_scale > 0.0d0) then
          write(message,*) 'inversion velocity scale (m/yr)              : ', &
               model%inversion%babc_velo_scale
          call write_log(message)
       endif
    endif   ! which_ho_powerlaw_c

    if (model%options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION) then
       write(message,*) 'coulomb_c max                                : ', &
            model%basal_physics%coulomb_c_max
       call write_log(message)
       write(message,*) 'coulomb_c min                                : ', &
            model%basal_physics%coulomb_c_min
       call write_log(message)
       write(message,*) 'inversion basal friction timescale (yr)      : ', &
            model%inversion%babc_timescale
       call write_log(message)
       if (model%inversion%babc_thck_scale > 0.0d0) then
          write(message,*) 'inversion thickness scale (m)                : ', &
               model%inversion%babc_thck_scale
          call write_log(message)
       endif
       if (model%inversion%babc_velo_scale > 0.0d0) then
          write(message,*) 'inversion velocity scale (m/yr)              : ', &
               model%inversion%babc_velo_scale
          call write_log(message)
       endif
    endif   ! which_ho_coulomb_c

    if (model%options%which_ho_bmlt_basin_inversion == HO_BMLT_BASIN_INVERSION_COMPUTE) then
       write(message,*) 'timescale (yr) for adjusting deltaT_basin    : ', model%inversion%bmlt_basin_timescale
       call write_log(message)
       write(message,*) 'dbmlt/dtemp scale (m/yr/deg C)               : ', model%inversion%dbmlt_dtemp_scale
       call write_log(message)
       write(message,*) 'Flotation threshold (m) for bmlt_basin inversion: ', &
            model%inversion%bmlt_basin_flotation_threshold
       call write_log(message)
       if (abs(model%inversion%bmlt_basin_mass_correction) > 0.0d0 .and. &
            model%inversion%bmlt_basin_number_mass_correction > 0) then
          write(message,*) 'Inversion mass correction applied to basin # :', &
               model%inversion%bmlt_basin_number_mass_correction
          call write_log(message)
          write(message,*) 'Mass correction (Gt)                         :', &
               model%inversion%bmlt_basin_mass_correction
          call write_log(message)
       endif
    endif

    if (model%basal_physics%beta_powerlaw_umax > 0.0d0) then
       write(message,*) 'max ice speed (m/yr) when evaluating beta(u) : ', model%basal_physics%beta_powerlaw_umax
       call write_log(message)
    endif

    if (model%basal_physics%beta_grounded_min > 0.d0) then
       write(message,*) 'min beta, grounded ice (Pa yr/m)             : ', model%basal_physics%beta_grounded_min
       call write_log(message)
    endif

    ! effective pressure parameters

    if (model%options%which_ho_effecpress == HO_EFFECPRESS_BPMP) then
       write(message,*) 'effective pressure delta      : ', model%basal_physics%effecpress_delta
       call write_log(message)
       write(message,*) 'effecpress bpmp threshold (K) : ', model%basal_physics%effecpress_bpmp_threshold
       call write_log(message)
    elseif (model%options%which_ho_effecpress == HO_EFFECPRESS_BWAT) then
       write(message,*) 'effective pressure delta      : ', model%basal_physics%effecpress_delta
       call write_log(message)
       write(message,*) 'effecpress bwat threshold (m) : ', model%basal_physics%effecpress_bwat_threshold
       call write_log(message)
    elseif (model%options%which_ho_effecpress == HO_EFFECPRESS_BWATFLX) then
       write(message,*) 'effecpress bwatflx threshold (m/yr) : ', model%basal_physics%effecpress_bwatflx_threshold
       call write_log(message)
       write(message,*) 'effecpress timescale (yr)     : ', model%basal_physics%effecpress_timescale
       call write_log(message)
       write(message,*) 'effective pressure delta      : ', model%basal_physics%effecpress_delta
       call write_log(message)
    elseif (model%options%which_ho_effecpress == HO_EFFECPRESS_BWAT_BVP) then
       write(message,*) 'effective pressure delta      : ', model%basal_physics%effecpress_delta
       call write_log(message)
       write(message,*) 'effecpress bwat threshold (m) : ', model%basal_physics%effecpress_bwat_threshold
       call write_log(message)
    endif

    if (model%basal_physics%p_ocean_penetration > 0.0d0) then
       call write_log('Apply ocean connection to reduce effective pressure')
       write(message,*) 'p_ocean_penetration           : ', model%basal_physics%p_ocean_penetration
       call write_log(message)
       if (model%basal_physics%ocean_p_timescale > 0.0d0) then
          write(message,*) 'ocean_p relaxation time (yr)  : ', model%basal_physics%ocean_p_timescale
          call write_log(message)
       endif
    endif

    if (model%numerics%idiag < 1 .or. model%numerics%idiag > model%general%ewn     &
                                        .or.                                                     &
        model%numerics%jdiag < 1 .or. model%numerics%jdiag > model%general%nsn) then
        call write_log('Error, global diagnostic point (idiag, jdiag) is out of bounds', GM_FATAL)
    endif

    ! ISMIP-HOM parameters
    if (model%numerics%periodic_offset_ew /= 0.d0) then
       write(message,*) 'periodic offset_ew (m)  : ',model%numerics%periodic_offset_ew
       call write_log(message)
    endif

    if (model%numerics%periodic_offset_ns /= 0.d0) then
       write(message,*) 'periodic offset_ns (m)  : ',model%numerics%periodic_offset_ns
       call write_log(message)
    endif

    ! initMIP parameters
    if (model%climate%acab_anomaly_timescale > 0.0d0) then
       write(message,*) 'acab_anomaly_timescale (yr): ', model%climate%acab_anomaly_timescale
       call write_log(message)
    endif

    if (model%options%overwrite_acab /= OVERWRITE_ACAB_NONE) then
       write(message,*) 'overwrite_acab_value (m/yr)   : ', model%climate%overwrite_acab_value
       call write_log(message)
       if (model%options%overwrite_acab == OVERWRITE_ACAB_THCKMIN) then
          write(message,*) 'overwrite_acab_minthck (m)    : ', model%climate%overwrite_acab_minthck
          call write_log(message)
       endif
    endif

    ! parameters for artm anomaly option
    if (model%climate%artm_anomaly_timescale > 0.0d0) then
       write(message,*) 'artm_anomaly_timescale (yr): ', model%climate%artm_anomaly_timescale
       call write_log(message)
    endif

    ! lapse rate
    if (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_LAPSE) then
       write(message,*) 'artm lapse rate (deg/m) : ', model%climate%t_lapse
       call write_log(message)
    endif

    if (model%basal_melt%bmlt_anomaly_timescale > 0.0d0) then
       write(message,*) 'bmlt_anomaly_timescale (yr): ', model%basal_melt%bmlt_anomaly_timescale
       call write_log(message)
    endif

    ! parameters for basal melting of floating ice (including MISMIP+ and MISOMIP)

    if (model%basal_melt%bmlt_cavity_h0 > 0.0d0 .and. &
        model%options%whichbmlt_float /= BMLT_FLOAT_MISMIP) then
       write(message,*) 'bmlt_cavity_h0 (m)       :  ', model%basal_melt%bmlt_cavity_h0
       call write_log(message)
    endif

    if (model%options%whichbmlt_float == BMLT_FLOAT_EXTERNAL) then
       if (model%basal_melt%bmlt_float_factor /= 1.0d0) then
          write(message,*) 'Input bmlt_float multiplied by: ', model%basal_melt%bmlt_float_factor
          call write_log(message)
       endif
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_CONSTANT) then
       write(message,*) 'bmlt_float_const (m/yr)  :  ', model%basal_melt%bmlt_float_const
       call write_log(message)
       write(message,*) 'bmlt_float_xlim (m)      :  ', model%basal_melt%bmlt_float_xlim
       call write_log(message)
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_MISMIP) then
       write(message,*) 'bmlt_float_omega (yr^-1) :  ', model%basal_melt%bmlt_float_omega
       call write_log(message)
       write(message,*) 'bmlt_float_h0 (m)        :  ', model%basal_melt%bmlt_float_h0
       call write_log(message)
       write(message,*) 'bmlt_float_z0 (m)        :  ', model%basal_melt%bmlt_float_z0
       call write_log(message)
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_DEPTH) then
       write(message,*) 'bmlt_float_depth_meltmax (m/yr):  ', model%basal_melt%bmlt_float_depth_meltmax
       call write_log(message)
       write(message,*) 'bmlt_float_depth_frzmax (m/yr) :  ', model%basal_melt%bmlt_float_depth_frzmax
       call write_log(message)
       write(message,*) 'bmlt_float_depth_zmeltmax (m)  :  ', model%basal_melt%bmlt_float_depth_zmeltmax
       call write_log(message)
       write(message,*) 'bmlt_float_depth_zmelt0 (m)    :  ', model%basal_melt%bmlt_float_depth_zmelt0
       call write_log(message)
       write(message,*) 'bmlt_float_depth_zfrzmax (m)   :  ', model%basal_melt%bmlt_float_depth_zfrzmax
       call write_log(message)
       write(message,*) 'warm ocean meltmin (m/yr)      :  ', model%basal_melt%bmlt_float_depth_meltmin
       call write_log(message)
       write(message,*) 'warm ocean zmeltmin (m)        :  ', model%basal_melt%bmlt_float_depth_zmeltmin
       call write_log(message)
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_MISOMIP) then
       write(message,*) 'T0 (deg C)               :  ', model%plume%T0
       call write_log(message)
       write(message,*) 'Tbot (deg C)             :  ', model%plume%Tbot
       call write_log(message)
       write(message,*) 'S0 (psu)                 :  ', model%plume%S0
       call write_log(message)
       write(message,*) 'Sbot (deg C)             :  ', model%plume%Sbot
       call write_log(message)
       write(message,*) 'gammaT (nondimensional)  :  ', model%plume%gammaT
       call write_log(message)
       write(message,*) 'gammaS (nondimensional)  :  ', model%plume%gammaS
       call write_log(message)
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_THERMAL_FORCING) then
       write(message,*) 'gamma0 (m/yr)            :  ', model%ocean_data%gamma0
       call write_log(message)
       if (model%ocean_data%thermal_forcing_anomaly /= 0.0d0) then
          write(message,*) 'thermal forcing anomaly (C) :', model%ocean_data%thermal_forcing_anomaly
          call write_log(message)
          write(message,*) 'TF anomaly start time (yr)  :', model%ocean_data%thermal_forcing_anomaly_tstart
          call write_log(message)
          write(message,*) 'TF anomaly timescale (yr)   :', model%ocean_data%thermal_forcing_anomaly_timescale
          call write_log(message)
          if (model%ocean_data%thermal_forcing_anomaly_basin == 0) then
             call write_log('TF anomaly will be applied to all basins')
          else
             write(message,*) 'TF anomaly will be applied to basin', &
                  model%ocean_data%thermal_forcing_anomaly_basin
             call write_log(message)
          endif
       endif
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
    call GetValue(section,'whichrelaxed',model%isostasy%whichrelaxed)
    call GetValue(section,'relaxed_tau',model%isostasy%relaxed_tau)
    call GetValue(section,'lithosphere_period',model%isostasy%period)

    !NOTE: This value used to be in a separate section ('elastic lithosphere')
    call GetValue(section,'flexural_rigidity',model%isostasy%rbel%d)

  end subroutine handle_isostasy

!--------------------------------------------------------------------------------

  subroutine print_isostasy(model)

    use glide_types
    use glimmer_log
    use cism_parallel, only: tasks

    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message
    
    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call write_log('Isostasy')
       call write_log('--------')

       if (model%isostasy%lithosphere==LITHOSPHERE_LOCAL) then
          call write_log('using local lithosphere approximation')
       else if (model%isostasy%lithosphere==LITHOSPHERE_ELASTIC) then
          call write_log('using elastic lithosphere approximation')
          if (tasks > 1) then
             call write_log('Warning, load calculation will be gathered to one processor; does not scale well',GM_WARNING)
          endif
          write(message,*) ' flexural rigidity : ', model%isostasy%rbel%d
          call write_log(message)
          write(message,*) ' lithosphere update period (yr): ', model%isostasy%period
          call write_log(message)
       else
          call write_log('Error, unknown lithosphere option',GM_FATAL)
       end if

       if (model%isostasy%asthenosphere==ASTHENOSPHERE_FLUID) then
          call write_log('using fluid mantle')
       else if (model%isostasy%asthenosphere==ASTHENOSPHERE_RELAXING) then
          call write_log('using relaxing mantle')
          write(message,*) ' characteristic time constant (yr): ', model%isostasy%relaxed_tau
          call write_log(message)
       else
          call write_log('Error, unknown asthenosphere option',GM_FATAL)
       end if

       if (model%isostasy%whichrelaxed==RELAXED_TOPO_DEFAULT) then
          call write_log('reading topg and relx as separate input fields')
       elseif (model%isostasy%whichrelaxed==RELAXED_TOPO_INPUT) then
          call write_log('setting relx to first slice of input topg')
       elseif (model%isostasy%whichrelaxed==RELAXED_TOPO_COMPUTE) then
          call write_log('computing relx, given that input topg is in equilibrium')
       else
          call write_log('Error, unknown whichrelaxed option',GM_FATAL)
       end if

       call write_log('')

    endif   ! compute isostasy

  end subroutine print_isostasy

!--------------------------------------------------------------------------------

  subroutine handle_glaciers(section, model)

    use glimmer_config
    use glide_types
    implicit none

    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'set_mu_star',    model%glacier%set_mu_star)
    call GetValue(section,'set_powerlaw_c', model%glacier%set_powerlaw_c)
    call GetValue(section,'t_mlt',          model%glacier%t_mlt)

  end subroutine handle_glaciers

!--------------------------------------------------------------------------------

  subroutine print_glaciers(model)

    use glide_types
    use glimmer_log

    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    ! glacier inversion options

    character(len=*), dimension(0:2), parameter :: glacier_set_mu_star = (/ &
         'spatially uniform glacier parameter mu_star', &
         'glacier-specific mu_star found by inversion', &
         'glacier-specific mu_star read from file    ' /)

    character(len=*), dimension(0:2), parameter :: glacier_set_powerlaw_c = (/ &
         'spatially uniform glacier parameter Cp', &
         'glacier-specific Cp found by inversion', &
         'glacier-specific Cp read from file    ' /)

    if (model%options%enable_glaciers) then

       call write_log(' ')
       call write_log('Glaciers')
       call write_log('--------')

       call write_log('Glacier tracking and tuning is enabled')

       write(message,*) 'set_mu_star               : ', model%glacier%set_mu_star, &
            glacier_set_mu_star(model%glacier%set_mu_star)
       call write_log(message)
       if (model%glacier%set_mu_star < 0 .or. &
            model%glacier%set_mu_star >= size(glacier_set_mu_star)) then
          call write_log('Error, glacier_set_mu_star option out of range', GM_FATAL)
       end if

       write(message,*) 'set_powerlaw_c            : ', model%glacier%set_powerlaw_c, &
            glacier_set_powerlaw_c(model%glacier%set_powerlaw_c)
       call write_log(message)
       if (model%glacier%set_powerlaw_c < 0 .or. &
            model%glacier%set_powerlaw_c >= size(glacier_set_powerlaw_c)) then
          call write_log('Error, glacier_set_powerlaw_c option out of range', GM_FATAL)
       end if

       write(message,*) 'glacier T_mlt (deg C)     :  ', model%glacier%t_mlt
       call write_log(message)

    endif   ! enable_glaciers

  end subroutine print_glaciers

!--------------------------------------------------------------------------------

  subroutine define_glide_restart_variables(model)

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
    call glide_add_to_restart_variable_list('topg thk temp bheatflx artm')

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
             call glide_add_to_restart_variable_list('acab')
          case (SMB_INPUT_MMYR_WE)
             call glide_add_to_restart_variable_list('smb')
          end select  ! smb_input

       case(SMB_INPUT_FUNCTION_XY_GRADZ)

          select case (options%smb_input)
          case (SMB_INPUT_MYR_ICE)
             call glide_add_to_restart_variable_list('acab_ref')
             call glide_add_to_restart_variable_list('acab_gradz')
          case (SMB_INPUT_MMYR_WE)
             call glide_add_to_restart_variable_list('smb_ref')
             call glide_add_to_restart_variable_list('smb_gradz')
          end select

          call glide_add_to_restart_variable_list('usrf_ref')

       case(SMB_INPUT_FUNCTION_XYZ)

          select case (options%smb_input)
          case (SMB_INPUT_MYR_ICE)
             call glide_add_to_restart_variable_list('acab_3d')
          case (SMB_INPUT_MMYR_WE)
             call glide_add_to_restart_variable_list('smb_3d')
          end select

          call glide_add_to_restart_variable_list('smb_levels')

    end select  ! smb_input_function

    ! Similarly for surface temperature (artm), based on options%artm_input
    ! Note: These options share usrf_ref and smb_levels with the SMB options above.

    select case(options%artm_input_function)

       case(ARTM_INPUT_FUNCTION_XY_GRADZ)
          call glide_add_to_restart_variable_list('artm_ref')
          call glide_add_to_restart_variable_list('artm_gradz')
          if (options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
             ! usrf_ref was added to restart above; nothing to do here
          else
             call glide_add_to_restart_variable_list('usrf_ref')
          endif

       case(ARTM_INPUT_FUNCTION_XYZ)
          call glide_add_to_restart_variable_list('artm_3d')
          if (options%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then
             ! smb_levels was added to restart above; nothing to do here
          else
             call glide_add_to_restart_variable_list('smb_levels')
          endif

       case(ARTM_INPUT_FUNCTION_XY_LAPSE)
          call glide_add_to_restart_variable_list('artm_ref')
          ! Note: Instead of artm_gradz, there is a uniform lapse rate
          if (options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
             ! usrf_ref was added to restart above; nothing to do here
          else
             call glide_add_to_restart_variable_list('usrf_ref')
          endif

    end select  ! artm_input_function

    ! Add anomaly forcing variables

    if (options%enable_acab_anomaly) then
       select case (options%smb_input)
       case (SMB_INPUT_MYR_ICE)
          call glide_add_to_restart_variable_list('acab_anomaly')
       case (SMB_INPUT_MMYR_WE)
          call glide_add_to_restart_variable_list('smb_anomaly')
       end select
    endif

    if (options%enable_artm_anomaly) then
       call glide_add_to_restart_variable_list('artm_anomaly')
    endif

    if (options%enable_bmlt_anomaly) then
       call glide_add_to_restart_variable_list('bmlt_float_anomaly')
    endif

    if (options%read_lat_lon) then
       ! If lat and lon are to be read from the input file, they should be written to the restart file
       call glide_add_to_restart_variable_list('lat lon')
    endif

    select case (options%whichbmlt_float)

       ! If bmlt_float is read from an external file at startup, then it needs to be in the restart file
       case (BMLT_FLOAT_EXTERNAL)
          call glide_add_to_restart_variable_list('bmlt_float_external')

       ! If prescribing a warm ocean mask for depth-dependent melting, this needs to be read on restart
       case (BMLT_FLOAT_DEPTH)
          call glide_add_to_restart_variable_list('warm_ocean_mask')

       case (BMLT_FLOAT_THERMAL_FORCING)

          ! Need the latest value of the thermal forcing field.
          ! This could be either the baseline value (if not updating during runtime), or a value read from a forcing file.
          ! If the latter, this field may not be needed, but include to be on the safe side, in case the forcing file
          !  is not read at restart.
          call glide_add_to_restart_variable_list('thermal_forcing')

          ! If using an ISMIP6 melt parameterization (either local or nonlocal),
          !  we need basin numbers and deltaT values for the parameterization.
          if (options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_LOCAL .or.  &
              options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL .or. &
              options%bmlt_float_thermal_forcing_param == BMLT_FLOAT_TF_ISMIP6_NONLOCAL_SLOPE) then
             call glide_add_to_restart_variable_list('basin_number')
             call glide_add_to_restart_variable_list('deltaT_basin')
          endif

    end select  ! whichbmlt_float

    ! add dycore specific restart variables
    select case (options%whichdycore)

      case (DYCORE_GLIDE)
        ! thkmask - TODO is this needed?
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
        call glide_add_to_restart_variable_list('thkmask wgrd wvel flwa uvel vvel')
    
        ! slip option for SIA
        select case (options%whichbtrc)
          case (0)
            ! no restart variable needed when no-slip is chosen
          case default
            ! when a slip option is chosen, ubas & vbas are needed by the temperature solver
            ! for calculating basal heating prior to the first calculation of velocity.
            ! Rather than recalculate the sliding field on restart, it is easier and 
            ! less error-prone to have them be restart variables.  
            ! This could either be done by making ubas, vbas restart variables or
            ! having them assigned from the bottom level of uvel,vvel on init
            ! Note that btrc and soft are not needed as restart variables because
            ! their current implementation is as a scalar ('basal_tract_const' config parameter).
            ! If they are ever implemented as 2-d fields, then they (probably just one of them)
            ! should become restart variables.
            
            ! Nothing needs to happen because ubas,vbas are assigned from uvel,vel in glide_init_state_diagnostic()
        end select

        ! basal water option for Glide
        select case (options%whichbwat)
        case (BWATER_NONE, BWATER_CONST)
           ! no restart variables needed
        case default
           ! restart needs to know bwat value
           call glide_add_to_restart_variable_list('bwat')
        end select


      case (DYCORE_GLISSADE)
        ! beta - b.c. needed for runs with sliding - could add logic to only include in that case
        ! flwa is not needed for glissade.
        ! TODO not sure if thkmask is needed for HO

        call glide_add_to_restart_variable_list('thkmask kinbcmask bfricflx dissip')

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
           call glide_add_to_restart_variable_list('uvel_extend vvel_extend')
        else
           call glide_add_to_restart_variable_list('uvel vvel')
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
                 call glide_add_to_restart_variable_list('uvel_2d_extend vvel_2d_extend btractx_extend btracty_extend efvs')
              else
                 call glide_add_to_restart_variable_list('uvel_2d vvel_2d btractx btracty efvs')
              endif
              
           case default
              ! Other approximations (including SSA and L1L2) use the 3D uvel and vvel to initialize the velocity

        end select   ! which_ho_approx
           
        ! mask used to limit computation to active blocks
        if (options%compute_blocks == ACTIVE_BLOCKS_ONLY) then
           call glide_add_to_restart_variable_list('ice_domain_mask')
        endif

        ! basal water option for Glissade
        select case (options%which_ho_bwat)
        case (HO_BWAT_NONE, HO_BWAT_CONSTANT)
           ! no restart variables needed
        case default
           ! restart needs to know bwat value
           call glide_add_to_restart_variable_list('bwat')
        end select

        ! grounding-line option for Glissade
        if (options%which_ho_flotation_function == HO_FLOTATION_FUNCTION_LINEAR_STDEV) then
           ! This option uses a correction based on topg_stdev to compute the flotation function.
           call glide_add_to_restart_variable_list('topg_stdev')
        endif

        ! calving options for Glissade

        !TODO: CALVING_GRID_MASK and apply_calving_mask are redundant; remove one option
        if (options%whichcalving == CALVING_GRID_MASK .or. options%apply_calving_mask) then
           call glide_add_to_restart_variable_list('calving_mask')
        endif

        if (options%whichcalving == CALVING_THCK_THRESHOLD) then
           call glide_add_to_restart_variable_list('thck_calving_threshold')
        endif

        ! The eigencalving calculation requires the product of eigenvalues of the horizontal strain rate tensor,
        !  which depends on the stress tensor, which is computed by the HO solver.
        ! On restart, the correct stress and strain rate tensors are not available, so we read in the eigenproduct.
        if (options%whichcalving == EIGENCALVING .or. options%whichcalving == CALVING_DAMAGE) then
           call glide_add_to_restart_variable_list('tau_eigen1')
           call glide_add_to_restart_variable_list('tau_eigen2')
        endif

        ! If forcing ice retreat, then we need ice_fraction_retreat_mask (which specifies the cells where retreat is forced)
        !  and reference_thck (which sets up an upper thickness limit for partly retreating cells)
        if (options%force_retreat /= FORCE_RETREAT_NONE) then
           call glide_add_to_restart_variable_list('ice_fraction_retreat_mask')
           call glide_add_to_restart_variable_list('reference_thck')
        endif

        ! other Glissade options

        ! If overwriting acab in certain grid cells, than overwrite_acab_mask needs to be in the restart file.
        ! This mask is read in at model initialization, or is set based on the input acab or ice thickness.
        if (options%overwrite_acab /= 0) then
           call glide_add_to_restart_variable_list('overwrite_acab_mask')
        endif

      end select ! which_dycore

    ! ==== Other non-dycore specific options ====

    ! internal water option (for enthalpy scheme)
    select case (options%whichtemp)
      case (TEMP_ENTHALPY)
        ! restart needs to know internal water fraction
        call glide_add_to_restart_variable_list('waterfrac')
      case default
        ! no restart variables needed
    end select

    ! basal sliding option
    select case (options%which_ho_babc)
       !WHL - Removed effecpress as a restart variable; it is recomputed with each velocity solve.
!!      case (HO_BABC_POWERLAW, HO_BABC_COULOMB_FRICTION, HO_BABC_COULOMB_POWERLAW_SCHOOF)
!!        ! These friction laws need effective pressure
!!        call glide_add_to_restart_variable_list('effecpress')
!!      case(HO_BABC_COULOMB_POWERLAW_TSAI)
!!        call glide_add_to_restart_variable_list('effecpress')
      case (HO_BABC_COULOMB_FRICTION, HO_BABC_COULOMB_POWERLAW_SCHOOF, HO_BABC_COULOMB_POWERLAW_TSAI)
         ! Note: These options compute beta internally, so it does not need to be in the restart file.
         if (options%use_c_space_factor) then
            ! c_space_factor needs to be in the restart file
            call glide_add_to_restart_variable_list('c_space_factor')
         endif
      case default
        ! Other HO basal boundary conditions may need the external beta field (although there are a few that don't)
        ! Note: If using beta from an external file, then 'beta' here needs to be the fixed, external field,
        !       and not the internal beta field that may have been weighted by the grounded fraction or otherwise adjusted.
        call glide_add_to_restart_variable_list('beta')
    end select

    ! basal friction options

    if (options%which_ho_powerlaw_c == HO_POWERLAW_C_INVERSION) then
       call glide_add_to_restart_variable_list('powerlaw_c')
       call glide_add_to_restart_variable_list('usrf_obs')
    elseif (options%which_ho_powerlaw_c == HO_POWERLAW_C_EXTERNAL) then
       call glide_add_to_restart_variable_list('powerlaw_c')
    endif

    if (options%which_ho_coulomb_c == HO_COULOMB_C_INVERSION) then
       call glide_add_to_restart_variable_list('coulomb_c')
       call glide_add_to_restart_variable_list('usrf_obs')

    elseif (options%which_ho_coulomb_c == HO_COULOMB_C_EXTERNAL) then
       call glide_add_to_restart_variable_list('coulomb_c')
    endif

    ! If inverting for coulomb_c or powerlaw_c based on observed surface speed
    ! (with model%inversion%babc_velo_scale > 0), then write velo_sfc_obs to the restart file.
    if (model%inversion%babc_velo_scale > 0.0d0) then
       call glide_add_to_restart_variable_list('velo_sfc_obs')
    endif

    ! effective pressure options
    ! f_effecpress_bwat represents the reduction of overburden pressure from bwatflx
    if (options%which_ho_effecpress == HO_EFFECPRESS_BWATFLX) then
       call glide_add_to_restart_variable_list('f_effecpress_bwat')
    endif

    ! f_effecpress_ocean_p represents the reduction of overburden pressure when ocean_p > 0
    ! Needs to be saved in case this fraction is relaxed over time toward (1 - Hf/H)^p
    if (model%basal_physics%p_ocean_penetration > 0.0d0) then
       call glide_add_to_restart_variable_list('f_effecpress_ocean_p')
    endif

    ! The bmlt_basin inversion option needs a thickness target for floating ice
    ! Note: deltaT_basin is added to the restart file above.
    if (options%which_ho_bmlt_basin_inversion == HO_BMLT_BASIN_INVERSION_COMPUTE) then
       call glide_add_to_restart_variable_list('floating_thck_target')
    endif

    ! geothermal heat flux option
    select case (options%gthf)
      case(GTHF_COMPUTE)
         ! restart needs to know lithosphere temperature
         call glide_add_to_restart_variable_list('litho_temp')
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
         call glide_add_to_restart_variable_list('relx')
         !WHL - The load field also needs to be in the restart file. The reason is that the load is updated
         !      at a period set by isostasy%period. If we restart between two updates, we need to use the most
         !      recently computed load. If we recompute the load right after restarting, the restart may not be exact.
         call glide_add_to_restart_variable_list('load')
      case default
         ! no new restart variables needed
    end select

    !WHL - added ice_age option
    !      Note: Ice age is a diagnostic field, not part of the prognostic ice state.
    !      Omitting it from restart will only break the diagnostic.
    select case (options%which_ho_ice_age)
       case(HO_ICE_AGE_COMPUTE)
          call glide_add_to_restart_variable_list('ice_age')
       case default
          ! no restart variables needed
    end select

    if (model%options%enable_glaciers) then
       ! Save some arrays related to glacier indexing
       call glide_add_to_restart_variable_list('rgi_glacier_id')
       call glide_add_to_restart_variable_list('cism_glacier_id')
       call glide_add_to_restart_variable_list('cism_glacier_id_init')
       call glide_add_to_restart_variable_list('cism_to_rgi_glacier_id')
       ! Save some arrays used to find the SMB and basal friction
       if (model%glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then
          call glide_add_to_restart_variable_list('usrf_obs')
          call glide_add_to_restart_variable_list('powerlaw_c')
       elseif (model%glacier%set_powerlaw_c == GLACIER_POWERLAW_C_EXTERNAL) then
          call glide_add_to_restart_variable_list('powerlaw_c')
       endif
       !TODO: Are area_target and volume_target needed?
       !      These could be computed based on cism_glacier_id_init and usrf_obs.
       call glide_add_to_restart_variable_list('glacier_volume_target')
       call glide_add_to_restart_variable_list('glacier_area_target')
       ! mu_star is needed only if relaxing toward the desired value;
       !  not needed if computed based on SMB = 0 over the target area
!!       call glide_add_to_restart_variable_list('glacier_mu_star')
    endif

    ! TODO bmlt was set as a restart variable, but I'm not sure when or if it is needed.

    ! TODO age should be a restart variable if it is an input variable.  
    ! Same goes for b.c. (bheatflxm, artm, acab) and any other tracers that get introduced.
    ! These could be included all the time (as I have down above for b.c.), or 
    ! we could add logic to only include them when they were in the input file.
    ! To do this, this subroutine would have to be moved to after where input files are read,
    ! glide_io_readall(), but before the output files are created, glide_io_createall()

    ! TODO lat is only needed for some climate drivers.  It is not needed for cism_driver.
    ! Need to add logic that will add it only when those drivers are used.

  end subroutine define_glide_restart_variables

!--------------------------------------------------------------------------------

end module glide_setup

!--------------------------------------------------------------------------------
