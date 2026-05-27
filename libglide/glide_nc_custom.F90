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
       if (.not.oc%append) then
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

    else if(.not. associated(model%funits%in_first)) then

       ! WHL: The following code goes back to early Glimmer.
       !      I think it was used to define grid output variables when there was no input file.
       !      It has not been used recently; not sure it still works.

       status = parallel_inq_varid(NCO%id,'x0',varid)
       do i=1, model%general%ewn-1
          status=nf90_put_var(NCO%id,varid,(real(i,dp)-0.5d0)*model%numerics%dew)
       enddo
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'y0',varid)
       do i=1, model%general%nsn-1
          status=nf90_put_var(NCO%id,varid,(real(i,dp)-0.5d0)*model%numerics%dns)
       enddo
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'x1',varid)
       do i=1, model%general%ewn-1
          status=nf90_put_var(NCO%id,varid,(real(i,dp)-1.0d0)*model%numerics%dew)
       enddo
       call nc_errorhandle(__FILE__,__LINE__,status)

       status = parallel_inq_varid(NCO%id,'y1',varid)
       do i=1, model%general%nsn-1
          status=nf90_put_var(NCO%id,varid,(real(i,dp)-1.0d0)*model%numerics%dns)
       enddo
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

    !TODO - Uncomment to add an ocean level dimension
    ! ocean level dimension
!    status = parallel_inq_varid(NCO%id,'zocn',varid)
!    status= parallel_put_var(NCO%id,varid,model%ocean_data%zocn)
!    call nc_errorhandle(__FILE__,__LINE__,status)

    !TODO - Uncomment to add an atm level dimension
    ! atm level dimension
!    status = parallel_inq_varid(NCO%id,'zatm',varid)
!    status= parallel_put_var(NCO%id,varid,model%climate%zatm)
!    call nc_errorhandle(__FILE__,__LINE__,status)

    ! glacier dimension

    if (model%options%enable_glaciers) then
       status = parallel_inq_varid(NCO%id,'glacierid',varid)
       status= parallel_put_var(NCO%id,varid,model%glacier%glacierid)
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    ! axis coordinate (used for CalvingMIP output)

    if (model%options%which_ho_calvingmip_domain /= HO_CALVINGMIP_DOMAIN_NONE) then
       status = parallel_inq_varid(NCO%id,'axis',varid)
       status= parallel_put_var(NCO%id,varid,model%calving%axis)
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    ! clean up
    deallocate(x0_global, y0_global)
    deallocate(x1_global, y1_global)

  end subroutine glide_nc_filldvars

end module glide_nc_custom
