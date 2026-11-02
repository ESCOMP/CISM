!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_stop.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glide_stop

  use glimmer_model_registry, only : deregister_model

  implicit none

  !> module containing finalisation of glide
  !> this subroutine had to be split out from glide.f90 to avoid circular dependencies

contains

  subroutine glide_finalise(model,forcewrite_arg)

    !> finalise model instance

    use glimmer_ncio
    use glimmer_log
    use glide_types
    use glide_io
    use cism_profile
    implicit none
    type(glide_global_type) :: model                 !> model instance
    logical, optional, intent(in) :: forcewrite_arg  !> if true, then force a write to output files
    character(len=100) :: message

    logical :: forcewrite         !> if true, then force a write to output files

    forcewrite = .false.

    ! force write to output files if specified by the optional input argument
    if (present(forcewrite_arg)) then
       forcewrite = forcewrite_arg
    end if

    ! force write to output files if set by a model option
    if (model%options%forcewrite_final) then
       forcewrite = .true.
    endif

    call glide_io_writeall(model, model, forcewrite)

    call closeall_in(model%funits)
    call closeall_out(model%funits)

    call glide_deallocarr(model)
    call deregister_model(model%model_id)

    ! write some statistics
    call write_log('Some Stats')
    write(message,*) 'Maximum temperature iterations: ',model%temper%niter
    call write_log(message)

    ! close profile
#if (defined PROFILING || defined CCSMCOUPLED || defined CESMTIMERS)
    call profile_close(model%profile)
#endif

  end subroutine glide_finalise

end module glide_stop
