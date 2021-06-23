!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_commandline.F90 - part of the Community Ice Sheet Model (CISM)  
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

#include "config.inc"

#ifdef HAVE_2003ARGS
#define NARGS   command_argument_count()
#define GETARG  get_command_argument
#else
#define NARGS   iargc
#define GETARG  getarg
#endif


module glint_commandline

  use glimmer_global, only:fname_length

  implicit none

  character(len=5000)         :: commandline_history            !< complete command line
  character(len=fname_length), dimension(:), allocatable :: &
                                 commandline_configname         !< name of ice sheet configuration file(s)
  character(len=fname_length) :: commandline_configname_scalar  !< configuration filenames, concatenated into a single string
  character(len=fname_length) :: commandline_resultsname        !< name of results file
  character(len=fname_length) :: commandline_climatename        !< name of climate configuration file

contains

  !> get the command line and parse it
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
  !! Modified by William Lipscomb, March 2021, to allow multiple ice sheet instances,
  !!  each with its own config file.
  !! The command line should be of the form:
  !! > ./cism_driver icesheet1.config icesheet2.config climate.config
  !! There can be any number of ice sheet config files, followed by a climate config file that is listed last.
  !! Note: Although Glint supports multiple ice sheet instances (plus a climate file) on the command line,
  !!       Glimmer supports only a single instance (without a climate file).

  subroutine glint_GetCommandline()

    use cism_parallel, only: main_task
    implicit none

    integer :: numargs, nfiles
    integer :: i
    integer :: num_icesheet_config  ! number of ice sheet config files

#ifndef HAVE_2003ARGS
    integer, external :: iargc
#endif
    character(len=100) :: argument
    integer, dimension(100) :: argumentIdx
    
    ! defaults
    commandline_resultsname = 'results'

    ! get number of arguments and file names
    numargs = NARGS

    ! reconstruct command line to store commandline_history
    call GETARG(0,commandline_history)

    do i = 1,numargs
       call GETARG(i,argument)
       commandline_history = trim(commandline_history)//" "//trim(argument)
    end do
    
    if (numargs > 0) then
       i = 0
       nfiles = 0
       ! loop over command line arguments
       do while (i < numargs)
          i = i + 1
          call GETARG(i,argument)
          ! check if it is an option
          if (argument(1:1) == '-') then
             select case (trim(argument))
             case ('-h')
                call glint_commandlineHelp()
                stop
             case ('-r')
                i = i + 1
                if (i > numargs) then
                   write(*,*) 'Error, expect name of output file to follow -o option'
                   call glint_commandlineHelp()
                   stop
                end if
                call GETARG(i,commandline_resultsname)
             case default
                write(*,*) 'Unknown option ',trim(argument)
                call glint_commandlineHelp()
                stop
             end select
          else
             ! it's not an option
             nfiles = nfiles+1
             argumentIdx(nfiles) = i
          end if
       end do

       if (numargs == 1) then      ! one ice sheet config file, no climate config file
          if (.not.allocated(commandline_configname)) &
               allocate(commandline_configname(1))
          call GETARG(i,commandline_configname(1))
       elseif (numargs > 1) then  ! one or more ice sheet config files, plus a climate config file
          if (.not.allocated(commandline_configname)) &
               allocate(commandline_configname(numargs-1))
          do i = 1, numargs-1
             call GETARG(i,commandline_configname(i))
          enddo
          ! assume the climate config file is listed last
          call GETARG(numargs,commandline_climatename)
       else
          write(*,*) 'Need at least one argument'
          call glint_commandlineHelp()
          stop
       end if
    else
       write(*,*) 'Enter name of climate configuration file'
       read(*,'(a)') commandline_climatename
       !TODO - Modify to allow > 1 ice sheet config files??
       write(*,*) 'Enter name of ice sheet configuration file to be read'
       allocate(commandline_configname(1))
       read(*,'(a)') commandline_configname
    end if

    ! If running multiple instances, concatenate the ice sheet config names into a single string
    commandline_configname_scalar = trim(commandline_configname(1))
    num_icesheet_config = size(commandline_configname)
    if (num_icesheet_config > 1) then
       do i = 2, num_icesheet_config
          commandline_configname_scalar = trim(commandline_configname_scalar)//' '//trim(commandline_configname(i))
       enddo
    endif

  end subroutine glint_GetCommandline

  !> print out command line
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
  !! Modified by William Lipscomb, March 2021
  subroutine glint_PrintCommandline()

    implicit none
    integer :: i

    write(*,*) 'Entire commandline'
    write(*,*) trim(commandline_history)
    write(*,*)
    write(*,*) 'commandline_climatename:   ', trim(commandline_climatename)
    write(*,*) 'commandline_configname(s): ', trim(commandline_configname_scalar)
    write(*,*) 'commandline_resultsname:   ', trim(commandline_resultsname)

  end subroutine glint_PrintCommandline

  !> print help message
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
  subroutine glint_commandlineHelp()
    implicit none
    character(len=500) :: pname

    call GETARG(0,pname)

    write(*,*) 'Usage: ',trim(pname),' [options] climname cfgname'
    write(*,*) 'where [options] are'
    write(*,*) '  -h:          this message'
    write(*,*) '  -r <fname>:  the name of the results file (default: results)'
  end subroutine glint_commandlineHelp

end module glint_commandline

