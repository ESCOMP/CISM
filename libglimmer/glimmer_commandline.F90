!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_commandline.F90 - part of the Community Ice Sheet Model (CISM)  
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

!> parsing common command line arguments
module glimmer_commandline

  use glimmer_global, only: fname_length

  use cism_parallel, only: main_task

  implicit none

  character(len=5000)         :: commandline_history     !< complete command line
  character(len=fname_length) :: commandline_configname  !< name of the configuration file
  character(len=fname_length) :: commandline_resultsname !< name of results file

contains

  !> get the command line and parse it
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
  subroutine glimmer_GetCommandline()

    implicit none

    integer numargs,nfiles
    integer :: i
#ifndef HAVE_2003ARGS
    integer, external :: iargc
#endif
    character(len=100) :: argument
    integer, dimension(100) :: argumentIdx
    
    ! defaults
    commandline_resultsname = 'results'

    if (main_task) then
       ! get number of arguments and file names
       numargs = NARGS
       ! reconstruct command line to store commandline_history
       call GETARG(0,commandline_history)
       do i=1,numargs
          call GETARG(i,argument)
          commandline_history = trim(commandline_history)//" "//trim(argument)
       end do
    
       if (numargs > 0) then
          i=0
          nfiles = 0
          ! loop over command line arguments
          do while (i < numargs)
             i = i + 1
             call GETARG(i,argument)
             ! check if it is an option
             if (argument(1:1) == '-') then
                select case (trim(argument))
                case ('-h')
                   call glimmer_commandlineHelp()
                   stop
                case ('-r')
                   i = i+1
                   if (i > numargs) then
                      write(*,*) 'Error, expect name of output file to follow -o option'
                      call glimmer_commandlineHelp()
                      stop
                   end if
                   call GETARG(i,commandline_resultsname)
                case default
                   write(*,*) 'Unkown option ',trim(argument)
                   call glimmer_commandlineHelp()
                   stop
                end select
             else
                ! it's not an option
                nfiles = nfiles+1
                argumentIdx(nfiles) = i
             end if
          end do
          if (nfiles > 0) then
             call GETARG(argumentIdx(1),commandline_configname)
          else
             write(*,*) 'Need at least one argument'
             call glimmer_commandlineHelp()
             stop
          end if
       else
          write(*,*) 'Enter name of GLIDE configuration file to be read'
          read(*,'(a)') commandline_configname
          !       commandline_configname = 'hump.config'
       end if   ! numargs > 0
    end if      ! main_task

  end subroutine glimmer_GetCommandline

  !> print out command line
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
  subroutine glimmer_PrintCommandline()

    implicit none

    if (main_task) then
       write(*,*) 'Entire commandline'
       write(*,*) trim(commandline_history)
       write(*,*)
       write(*,*) 'commandline_configname:  ', trim(commandline_configname)
       write(*,*) 'commandline_resultsname: ', trim(commandline_resultsname)
    endif
  end subroutine glimmer_PrintCommandline

  !> print help message
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
  subroutine glimmer_commandlineHelp()

    implicit none
    character(len=500) :: pname

    call GETARG(0,pname)

    if (main_task) then
       write(*,*) 'Usage: ',trim(pname),' [options] cfgname'
       write(*,*) 'where [options] are'
       write(*,*) '  -h:          this message'
       write(*,*) '  -r <fname>:  the name of the results file (default: results)'
    endif
  end subroutine glimmer_commandlineHelp
end module glimmer_commandline
