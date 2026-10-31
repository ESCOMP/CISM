!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_nc_custom.F90 - part of the Community Ice Sheet Model (CISM)  
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

#define NCO outfile%nc

module glide_nc_custom

  !module for filling in dimension variables

  use glimmer_global, only: dp

  implicit none

contains

  subroutine glide_nc_fillall(model, outfiles)

    !> fill dimension variables of all files
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none

    type(glide_global_type) :: model
    type(glimmer_nc_output),pointer,optional :: outfiles

    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       ! Note: Files written with one_file_per_write do not exist yet; they are filled when created.
       if (.not.oc%append .and. .not.oc%one_file_per_write) then
          call glide_nc_filldvars(oc, model)
       endif
       oc=>oc%next
    end do

  end subroutine glide_nc_fillall


  subroutine glide_nc_filldvars(outfile,    model)

    use glide_types
    use glimmer_ncdf
    use cism_parallel, only: parallel_inq_varid, parallel_put_var, parallel_enddef
    implicit none

    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model

    integer i,status,varid

    ! check if we are still in define mode and if so leave it
    if (NCO%define_mode) then
       status = parallel_enddef(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%define_mode = .FALSE.
    end if

    ! horizontal dimensions
    ! (x1,y1) is the unstaggered scalar grid
    ! (x0,y0) is the staggered velocity grid

    if (associated(model%funits%in_first)) then

       ! WHL, 8/5/19:
       ! The original code called subroutine distributed_put_var to gather model%general%x1, etc.
       !  from local tasks into a global array that was written to the output file.
       ! This does not work, in general, when computing on active blocks only, because the local versions
       !  of model%general%x1 may not span the global domain.
       ! The revised code calls parallel_put_var to write (x0,y0) and (x1,y1) to the output file
       !  based on the global arrays (x0_global,y0_global) and (x1_global,y1_global).
       ! This assumes that the global arrays were read in or computed at startup.

       status = parallel_inq_varid(NCO%id,'x1',varid)
       status = parallel_put_var(NCO%id,varid,model%general%x1_global)
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'y1',varid)
       status = parallel_put_var(NCO%id,varid,model%general%y1_global)
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'x0',varid)
       status = parallel_put_var(NCO%id,varid,model%general%x0_global)
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'y0',varid)
       status = parallel_put_var(NCO%id,varid,model%general%y0_global)
       call nc_errorhandle(__FILE__,__LINE__,status)

    else   ! no input file; compute the coordinates from the grid spacing

       ! Note: The original code (from early Glimmer) wrote one scalar per call without a 'start' index,
       !       so each call overwrote the first element, leaving the rest of the array as fill values.

       status = parallel_inq_varid(NCO%id,'x0',varid)
       status = parallel_put_var(NCO%id, varid, &
            (/ ((real(i,dp)-0.5d0)*model%numerics%dew, i = 1, model%parallel%global_ewn-1) /))
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'y0',varid)
       status = parallel_put_var(NCO%id, varid, &
            (/ ((real(i,dp)-0.5d0)*model%numerics%dns, i = 1, model%parallel%global_nsn-1) /))
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'x1',varid)
       status = parallel_put_var(NCO%id, varid, &
            (/ ((real(i,dp)-1.0d0)*model%numerics%dew, i = 1, model%parallel%global_ewn) /))
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'y1',varid)
       status = parallel_put_var(NCO%id, varid, &
            (/ ((real(i,dp)-1.0d0)*model%numerics%dns, i = 1, model%parallel%global_nsn) /))
       call nc_errorhandle(__FILE__,__LINE__,status)

    end if   ! associated(model%funits%in_first)

    ! layer interfaces

    status = parallel_inq_varid(NCO%id,'level',varid)
    status = parallel_put_var(NCO%id,varid,model%numerics%sigma)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! layer midpoints

    status = parallel_inq_varid(NCO%id,'staglevel',varid)
    status = parallel_put_var(NCO%id,varid,model%numerics%stagsigma)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! layer midpoints, plus upper and lower surfaces
    ! (e.g., temperature field in HO dycore)

    status = parallel_inq_varid(NCO%id,'stagwbndlevel',varid)
    status = parallel_put_var(NCO%id,varid,model%numerics%stagwbndsigma)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! lithosphere vertical coordinate

    if (model%options%gthf == GTHF_COMPUTE) then
       status = parallel_inq_varid(NCO%id,'lithoz',varid)
       status= parallel_put_var(NCO%id,varid,model%lithot%deltaz)
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

  end subroutine glide_nc_filldvars

end module glide_nc_custom
