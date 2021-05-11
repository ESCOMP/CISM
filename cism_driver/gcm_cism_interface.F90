!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   gcm_cism_interface.F90 - part of the Community Ice Sheet Model (CISM)  
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

! from glide_types.F90:
!  integer, parameter :: DYCORE_GLIDE = 0     ! old shallow-ice dycore from Glimmer
!  integer, parameter :: DYCORE_GLAM = 1      ! Payne-Price finite-difference solver
!  integer, parameter :: DYCORE_GLISSADE = 2  ! finite-element solver
!  integer, parameter :: DYCORE_ALBANYFELIX = 3  ! External Albany-Felix finite-element solver
!  integer, parameter :: DYCORE_BISICLES = 4     ! BISICLES external dycore

module gcm_cism_interface

  use glint_commandline
  use glide
  use cism_front_end

  use glint_example_clim
  use glint_main
  use gcm_to_cism_glint

  ! Note: Options 0 and 1 used to be called GCM_MINIMAL_MODEL AND GCM_DATA_MODEL.
  !       'BASIC' refers to standalone CISM with a single input data file.
  !       'GLINT' refers to runs with the Glint interface, with a climate config file
  !         and one or more ice sheet config files.
  !       The 'GCM' prefix is a misnomer, in that both of these options use standalone CISM,
  !         uncoupled to a global climate model.

  integer, parameter :: GCM_BASIC_MODEL = 0
  integer, parameter :: GCM_GLINT_MODEL = 1
  integer, parameter :: GCM_CESM = 2

contains

subroutine gci_init_interface(which_gcm,g2c)

  use glint_commandline
  use glimmer_config
  use glide
  use glide_types
  use cism_parallel, only: main_task
 
  use cism_front_end 

  integer, intent(in) :: which_gcm
  type(gcm_to_cism_type) :: g2c   ! holds everything

  integer :: whichdycore
  type(ConfigSection), pointer :: config   ! configuration stuff
  type(ConfigSection), pointer :: section  ! pointer to the section to be checked

  character(len=fname_length) :: fname_root ! root of log file name
  integer :: num_icesheet_config            ! number of ice sheet config files
  integer :: i

  ! call parallel_initialise
  
  ! Set the name of the log file
  ! Note: The log file name is constructed from the name of the config file.
  !       If there are multiple config files, the log file name is built by concatenating them.
  call glint_GetCommandline()

  fname_root = trim(commandline_configname(1))
  num_icesheet_config = size(commandline_configname)
  if (num_icesheet_config > 1) then
     do i = 2, num_icesheet_config
        fname_root = trim(fname_root)//'..'//trim(commandline_configname(i))
     enddo
  endif
  call open_log(unit=50, fname=logname(fname_root))

  ! Get the CISM dycore to be used
  ! Note: If there are multiple ice sheet config files, CISM assumes that each will use the same dycore,
  !        and gets the dycore from the first config file in the configname array.
  !       Either all instances are basic or all instances use glint; a mix is not supported.
  call ConfigRead(commandline_configname(1),config)
  call GetSection(config,section,'options')
  call GetValue(section,'dycore',whichdycore)
  if (main_task) print *,'CISM dycore type (0=Glide, 1=Glam, 2=Glissade, 3=AlbanyFelix, 4 = BISICLES) = ', whichdycore  

  ! Check to see if running basic GCM or glint GCM.  Still need to add CESM GCM:
  call GetSection(config,section,'GLINT climate')

  if (associated(section)) then
     g2c%which_gcm = GCM_GLINT_MODEL
  else 
     g2c%which_gcm = GCM_BASIC_MODEL
  end if
  if (main_task) print *,'g2c%which_gcm (0 = basic, 1 = Glint) = ', g2c%which_gcm

  select case (g2c%which_gcm)
    case (GCM_BASIC_MODEL)
      if (main_task) print*, 'call cism_init_dycore'
      call cism_init_dycore(g2c%glide_model)
 
    case (GCM_GLINT_MODEL)
      if (main_task) print*, 'call g2c_glint_init'
      call g2c_glint_init(g2c)

    case (GCM_CESM)
      ! call gcm_glint_GetCommandline_proxy()
      ! call g2c_glint_init(g2c) 

    case default
      if (main_task) print *,"Error -- unknown GCM type."
  end select

end subroutine gci_init_interface   

subroutine gci_run_model(g2c)

  type(gcm_to_cism_type) :: g2c 

  logical :: finished = .false.

  do while (.not. finished)
    select case (g2c%which_gcm)
      case (GCM_BASIC_MODEL)
        ! call gcm_update_model(gcm_model,cism_model)
!        if (main_task) print *,"In gci_run_model, calling cism_run_dycore"
        call cism_run_dycore(g2c%glide_model)

      case (GCM_GLINT_MODEL,GCM_CESM)
!        if (main_task) print *,"In gci_run_model, calling g2c_glint_run"
        call g2c_glint_run(g2c)
        call g2c_glint_climate_time_step(g2c)
      case default
    end select
    finished = (gci_finished(g2c))
  end do
end subroutine gci_run_model


! gci_finished is used to test status of GCM
function gci_finished(g2c) result(finished)

  type(gcm_to_cism_type) :: g2c
  logical :: finished
 
  select case (g2c%which_gcm)
    case (GCM_BASIC_MODEL)
      finished = .true.

    case (GCM_GLINT_MODEL,GCM_CESM)
      call g2c_glint_check_finished(g2c,finished)
    case default
  end select
  !if (main_task) print *,"In gci_finished, finished = ",finished  

end function gci_finished


subroutine gci_finalize_interface(g2c)

  type(gcm_to_cism_type) :: g2c

  select case (g2c%which_gcm)
    case (GCM_BASIC_MODEL)
      call cism_finalize_dycore(g2c%glide_model)
 
    case (GCM_GLINT_MODEL)
      call g2c_glint_end(g2c)

    case (GCM_CESM)
      ! call g2c_glint_end(g2c)
    case default
  end select

end subroutine gci_finalize_interface


end module gcm_cism_interface
