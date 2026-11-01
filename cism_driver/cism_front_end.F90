!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_front_end.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_front_end
  ! The CISM front-end is used to connect both the standalone driver
  ! (cism_driver) or the CISM interface to CESM (cism_cesm_interface),
  ! to the dycore.
  ! There is a separate init/run/finalize routine for each dycore (Glide and Glissade),
  !  since each dycore has its own model instance type.

  use eismint_forcing, only: eismint_climate_type

  implicit none
  private
  public :: cism_init_glide, cism_run_glide, cism_finalize_glide
  public :: cism_init_glissade, cism_run_glissade, cism_finalize_glissade

  ! parameters for the EISMINT climate (used only for EISMINT test cases)
  type(eismint_climate_type), save :: eismint_climate

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Glide
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine cism_init_glide(model)

  use glimmer_global, only: dp
  use glimmer_model_registry, only : register_model
  use glide
  use eismint_forcing, only: eismint_initialise
  use glimmer_log
  use glimmer_config
  use glide_nc_custom, only: glide_nc_fillall
  use glide_io, only: glide_io_writeall
  use glide_diagnostics, only: glide_write_diagnostics
  use profile, only: profile_init, t_startf, t_stopf, t_adj_detailf

  implicit none

  type(glide_global_type) :: model        ! model instance
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=dp) :: time                   ! model time in years

  ! read the configuration
  call cism_front_end_read_config(config)

  ! initialise profiling
  call profile_init(model%profile,'glide.profile')

  call t_startf('cism')
 
  ! initialise GLIDE
  call t_startf('initialization')

  call register_model(model%model_id)

  call glide_config(model,config)

  ! This call is needed only if running the EISMINT test cases
  call eismint_initialise(eismint_climate,config)

  call glide_initialise(model)

  call CheckSections(config)

  ! fill dimension variables on output files
  call glide_nc_fillall(model)

  time = model%numerics%tstart
  model%numerics%time = time

  ! Set EISMINT forcing for initial time
  call set_eismint_forcing_glide(model, time)

  ! read forcing time slice if needed - this will overwrite values from IC file if there is a conflict.
  ! Note: The first 'model' is passed to the argument 'data', which is filled by calling glide_read.
  call glide_read_forcing(model, model)

  ! Optionally, read all the time slices at once from selected forcing files.
  call glide_read_forcing_once(model, model)

  call spinup_lithot(model)

  call t_stopf('initialization')

  !MJH Created this block here to fill out initial state without needing to enter time stepping loop.  This allows
  ! a run with tend=tstart to be run without time-stepping at all.  It requires solving all diagnostic (i.e. not
  ! time depdendent) variables (most important of which is velocity) for the initial state and then writing the 
  ! initial state as time 0 (or more accurately, as time=tstart).  Also, halo updates need to occur after the 
  ! diagnostic variables are calculated.

! ------------- Calculate initial state and output it -----------------

  call t_startf('initial_diag_var_solve')

  if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
     ! disable further profiling in normal usage
     call t_adj_detailf(+10)
  endif

  call glide_init_state_diagnostic(model)

  if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
     ! restore profiling to normal settings
     call t_adj_detailf(-10)
  endif

  call t_stopf('initial_diag_var_solve')

  ! Write initial diagnostic output to log file
  ! Note: tstep_count is set to 0 at model initialization and then is incremented in cism_run_glide
  !  before each call to the dycore.

  call t_startf('initial_write_diagnostics')
  call glide_write_diagnostics(model,        time,       &
                               tstep_count = model%numerics%tstep_count)
  call t_stopf('initial_write_diagnostics')

  ! --- Output the initial state -------------
  ! Note: For a standard restart, the initial state is not output, because this state
  !       should already have been written to the output file when the previous run ended.

  if (model%options%is_restart == NO_RESTART .or. model%options%is_restart == HYBRID_RESTART) then
     call t_startf('initial_io_writeall')
     call glide_io_writeall(model, model, time=time)
     call t_stopf('initial_io_writeall')
  endif

end subroutine cism_init_glide


subroutine cism_run_glide(model)

  use glimmer_global, only: dp
  use glide
  use glide_io, only: glide_io_writeall
  use glide_diagnostics, only: glide_write_diagnostics
  use profile, only: t_startf, t_stopf

  implicit none

  type(glide_global_type) :: model        ! model instance
  real(kind=dp) :: time                   ! model time in years
  real(kind=dp) :: time_eps               ! tolerance within which times are equal 

  time = model%numerics%tstart
  time_eps = model%numerics%tinc/1000.0d0

  ! ------------- Begin time step loop -----------------
 
  ! check if we're doing any evolution
  if (time < model%numerics%tend) then
    do while(time + time_eps < model%numerics%tend)

      ! Forcing from a 'forcing' data file - will read time slice if needed
      ! (See the notes in cism_run_glissade.)
      call t_startf('read_forcing')
      call glide_read_forcing(model, model)

      ! If any forcing data have been read once into Fortran arrays at initialization,
      !  simply copy the data based on the current forcing time.
      call glide_retrieve_forcing(model, model)
      call t_stopf('read_forcing')

      ! Increment time step
      time = time + model%numerics%tinc
      model%numerics%time = time  ! TODO This is redundant with what is happening in glide, but this is needed for forcing to work properly.

      call t_startf('tstep')

      call t_startf('glide_tstep_p1')
      call glide_tstep_p1(model,time)
      call t_stopf('glide_tstep_p1')

      call t_startf('glide_tstep_p2')
      call glide_tstep_p2(model)
      call t_stopf('glide_tstep_p2')

      call t_startf('glide_tstep_p3')
      call glide_tstep_p3(model)
      call t_stopf('glide_tstep_p3')

      call t_stopf('tstep')

      ! write ice sheet diagnostics to log file at desired interval (model%numerics%dt_diag)

      call t_startf('write_diagnostics')
      call glide_write_diagnostics(model,        time,       &
                                   tstep_count = model%numerics%tstep_count)
      call t_stopf('write_diagnostics')

      ! update time from dycore advance
      model%numerics%time = time

      ! --- Set forcing ---
      ! (See the notes in cism_run_glissade.)
      call t_startf('set_forcing')
      call set_eismint_forcing_glide(model, time)
      call t_stopf('set_forcing')

      ! Write to output netCDF files at desired intervals
      call t_startf('io_writeall')
      call glide_io_writeall(model,model)
      call t_stopf('io_writeall')
    end do   ! time < model%numerics%tend
  endif    
  
end subroutine cism_run_glide


subroutine cism_finalize_glide(model)

  use glide_types, only: glide_global_type
  use glide_stop, only: glide_finalise
  use glimmer_log, only: close_log
  use profile, only: t_stopf

  implicit none

  type(glide_global_type) :: model        ! model instance

  call t_stopf('cism')

  ! finalise the model
  call glide_finalise(model)

  call close_log

end subroutine cism_finalize_glide


subroutine set_eismint_forcing_glide(model, time)

  ! Set the surface mass balance and air temperature for EISMINT test cases
  ! Note: These calls do nothing unless an EISMINT case is run.

  use glimmer_global, only: dp
  use glide_types, only: glide_global_type
  use eismint_forcing, only: eismint_massbalance, eismint_surftemp, eismint_exact_forcing

  implicit none

  type(glide_global_type) :: model        ! model instance
  real(dp), intent(in) :: time            ! model time in years

  call eismint_massbalance(eismint_climate, time, model%numerics%dew, model%parallel, &
                           model%climate%acab)
  call eismint_surftemp(eismint_climate, time, model%numerics%dew, model%parallel, &
                        model%geometry%thck, model%climate%artm)

  ! exact verification test (Glide only)
  call eismint_exact_forcing(eismint_climate, model, time)

end subroutine set_eismint_forcing_glide

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Glissade
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine cism_init_glissade(model)

  use glimmer_global, only: dp
  use glimmer_model_registry, only : register_model
  use glissade_types
  use glissade, only: glissade_initialise, glissade_diagnostic_variable_solve
  use glissade_setup, only: glissade_config
  use eismint_forcing, only: eismint_initialise
  use glimmer_log
  use glimmer_config
  use glissade_nc_custom, only: glissade_nc_fillall
  use glissade_io, only: glissade_io_writeall, glissade_read_forcing, glissade_read_forcing_once
  use glissade_diagnostics, only: glissade_write_diagnostics
  use profile, only: profile_init, t_startf, t_stopf, t_adj_detailf

  implicit none

  type(glissade_global_type) :: model     ! model instance
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=dp) :: time                   ! model time in years

  ! read the configuration
  call cism_front_end_read_config(config)

  ! initialise profiling
  call profile_init(model%profile,'glide.profile')

  call t_startf('cism')
 
  call t_startf('initialization')

  call register_model(model%model_id)

  call glissade_config(model,config)

  ! This call is needed only if running the EISMINT test cases
  call eismint_initialise(eismint_climate,config)
  if (eismint_climate%eismint_type == 4) then
     call write_log('The EXACT verification test is supported for the Glide dycore only', GM_FATAL)
  endif

  call glissade_initialise(model)

  call CheckSections(config)

  ! fill dimension variables on output files
  call glissade_nc_fillall(model)

  time = model%numerics%tstart
  model%numerics%time = time    ! MJH added 1/10/13 - the initial diagnostic glissade solve won't know 
                                !                     the correct time on a restart unless we set it here.

  ! Set EISMINT forcing for initial time
  call set_eismint_forcing_glissade(model, time)

  ! read forcing time slice if needed - this will overwrite values from IC file if there is a conflict.
  ! Note: The first 'model' is passed to the argument 'data', which is filled by calling glissade_read.
  call glissade_read_forcing(model, model)

  ! Optionally, read all the time slices at once from selected forcing files.
  call glissade_read_forcing_once(model, model)

  call t_stopf('initialization')

  !MJH Created this block here to fill out initial state without needing to enter time stepping loop.  This allows
  ! a run with tend=tstart to be run without time-stepping at all.  It requires solving all diagnostic (i.e. not
  ! time depdendent) variables (most important of which is velocity) for the initial state and then writing the 
  ! initial state as time 0 (or more accurately, as time=tstart).  Also, halo updates need to occur after the 
  ! diagnostic variables are calculated.

! ------------- Calculate initial state and output it -----------------

  call t_startf('initial_diag_var_solve')

  if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
     ! disable further profiling in normal usage
     call t_adj_detailf(+10)
  endif

  ! solve the remaining diagnostic variables for the initial state
  call glissade_diagnostic_variable_solve(model)  ! velocity, usrf, etc.

  if (model%numerics%tstart < (model%numerics%tend - model%numerics%tinc)) then
     ! restore profiling to normal settings
     call t_adj_detailf(-10)
  endif

  call t_stopf('initial_diag_var_solve')

  ! Write initial diagnostic output to log file
  ! Note: tstep_count is set to 0 at model initialization and then is incremented in cism_run_glissade
  !  before each call to the dycore.

  call t_startf('initial_write_diagnostics')
  call glissade_write_diagnostics(model,        time,       &
                                  tstep_count = model%numerics%tstep_count)
  call t_stopf('initial_write_diagnostics')

  ! --- Output the initial state -------------
  ! Note: For a standard restart, the initial state is not output, because this state
  !       should already have been written to the output file when the previous run ended.

  if (model%options%is_restart == NO_RESTART .or. model%options%is_restart == HYBRID_RESTART) then
     call t_startf('initial_io_writeall')
     call glissade_io_writeall(model, model, time=time)
     call t_stopf('initial_io_writeall')
  endif

end subroutine cism_init_glissade


subroutine cism_run_glissade(model)

  use glimmer_global, only: dp
  use glissade_types
  use glissade, only: glissade_tstep
  use glissade_io, only: glissade_io_writeall, glissade_read_forcing, glissade_retrieve_forcing
  use glissade_diagnostics, only: glissade_write_diagnostics
  use profile, only: t_startf, t_stopf

  implicit none

  type(glissade_global_type) :: model     ! model instance
  real(kind=dp) :: time                   ! model time in years
  real(kind=dp) :: time_eps               ! tolerance within which times are equal 

  time = model%numerics%tstart
  time_eps = model%numerics%tinc/1000.0d0

  ! ------------- Begin time step loop -----------------
 
  ! check if we're doing any evolution
  if (time < model%numerics%tend) then
    do while(time + time_eps < model%numerics%tend)

      !!! SFP moved block of code for applying time dependent forcing read in from netCDF here,
      !!! as opposed to at the end of the time step (commented it out in its original location for now)
      !!! This is a short-term fix. See additional discussion as part of issue #19 (in cism-piscees github repo).

      ! Forcing from a 'forcing' data file - will read time slice if needed
      call t_startf('read_forcing')
      call glissade_read_forcing(model, model)

      ! If any forcing data have been read once into Fortran arrays at initialization,
      !  simply copy the data based on the current forcing time.
      call glissade_retrieve_forcing(model, model)
      call t_stopf('read_forcing')

      ! Increment time step
      time = time + model%numerics%tinc
      model%numerics%time = time  ! TODO This is redundant with what is happening in glissade, but this is needed for forcing to work properly.

      call t_startf('tstep')
      call glissade_tstep(model,time)
      call t_stopf('tstep')

      ! write ice sheet diagnostics to log file at desired interval (model%numerics%dt_diag)

      call t_startf('write_diagnostics')
      call glissade_write_diagnostics(model,        time,       &
                                      tstep_count = model%numerics%tstep_count)
      call t_stopf('write_diagnostics')

      ! update time from dycore advance
      model%numerics%time = time

      ! --- Set forcing ---
      ! Setting forcing at the end of the time step maintains consistency
      ! with a forward Euler time step and ensures consistency of the time stamp
      ! to fields in input and output files.  
      ! For forward Euler time stepping we want S^n+1 = g(S^n, F^n)
      ! where S is the model state, F is forcing, and n, n+1 are time levels
      ! We also want a forcing field in the output file to have a time stamp
      ! that matches its time stamp in the input file or the EISMINT analytic function.
      ! The simplest way to ensure both of these criteria is to set forcing at the
      ! end of each time step.
      ! EISMINT forcing
      ! NOTE: these only do something when an EISMINT case is run
      call t_startf('set_forcing')
      call set_eismint_forcing_glissade(model, time)
      call t_stopf('set_forcing')

      ! Write to output netCDF files at desired intervals
      call t_startf('io_writeall')
      call glissade_io_writeall(model,model)
      call t_stopf('io_writeall')
    end do   ! time < model%numerics%tend
  endif    
  
end subroutine cism_run_glissade


subroutine cism_finalize_glissade(model)

  use glissade_types
  use glissade, only: glissade_finalise
  use glimmer_log, only: close_log
  use profile, only: t_stopf

  implicit none

  type(glissade_global_type) :: model     ! model instance

  call t_stopf('cism')

  ! finalise the model
  call glissade_finalise(model)

  call close_log

end subroutine cism_finalize_glissade


subroutine set_eismint_forcing_glissade(model, time)

  ! Set the surface mass balance and air temperature for EISMINT test cases
  ! Note: These calls do nothing unless an EISMINT case is run.

  use glimmer_global, only: dp
  use glissade_types
  use eismint_forcing, only: eismint_massbalance, eismint_surftemp

  implicit none

  type(glissade_global_type) :: model     ! model instance
  real(dp), intent(in) :: time            ! model time in years

  call eismint_massbalance(eismint_climate, time, model%numerics%dew, model%parallel, &
                           model%climate%acab)
  call eismint_surftemp(eismint_climate, time, model%numerics%dew, model%parallel, &
                        model%geometry%thck, model%climate%artm)

end subroutine set_eismint_forcing_glissade

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Steps shared by both dycores
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine cism_front_end_read_config(config)

  ! Read the command line and the configuration file

  use glimmer_config
  use glimmer_commandline
  use glimmer_filenames, only : filenames_init

  implicit none

  type(ConfigSection), pointer :: config  ! configuration stuff

  call glimmer_GetCommandline()

  ! DMR -- open_log call commented out, since called in gci_init_interface()
  ! start logging
  ! call open_log(unit=50, fname=logname(commandline_configname))
  ! Note: In principle, the standalone CISM driver could support multiple ice sheet instances.
  !       At least for now, assume that there is only one instance, with one config file.

  ! setup paths
  call filenames_init(commandline_configname)

  ! read configuration
  call ConfigRead(commandline_configname,config)

end subroutine cism_front_end_read_config

end module cism_front_end
