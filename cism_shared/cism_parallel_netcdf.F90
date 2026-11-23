!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_parallel_netcdf.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_parallel_netcdf

  ! Part of the cism_parallel module:
  !  netCDF wrappers: the main task does the I/O and broadcasts the result
  ! Other modules use these through cism_parallel.

  use netcdf
  use glimmer_global, only : dp, sp
  use cism_parallel_core
  use cism_parallel_broadcast

  implicit none

  interface parallel_def_var
     module procedure parallel_def_var_dimids
     module procedure parallel_def_var_nodimids
  end interface

  interface parallel_get_att
     module procedure parallel_get_att_character
     module procedure parallel_get_att_real4
     module procedure parallel_get_att_real4_1d
     module procedure parallel_get_att_real8
     module procedure parallel_get_att_real8_1d
  end interface

  interface parallel_get_var
     module procedure parallel_get_var_integer
     module procedure parallel_get_var_real4
     module procedure parallel_get_var_real8
     module procedure parallel_get_var_integer_1d
     module procedure parallel_get_var_real4_1d
     module procedure parallel_get_var_real8_1d
     module procedure parallel_get_var_integer_2d
     module procedure parallel_get_var_real8_2d
  end interface

  interface parallel_put_att
     module procedure parallel_put_att_character
     module procedure parallel_put_att_integer
     module procedure parallel_put_att_real4
     module procedure parallel_put_att_real4_1d
     module procedure parallel_put_att_real8
     module procedure parallel_put_att_real8_1d
  end interface

  interface parallel_put_var
     module procedure parallel_put_var_integer
     module procedure parallel_put_var_integer_1d
     module procedure parallel_put_var_real4
     module procedure parallel_put_var_real8
     module procedure parallel_put_var_real8_1d
  end interface

contains

!=======================================================================

  function parallel_close(ncid)

    implicit none
    integer :: ncid

    integer :: parallel_close

    ! begin
    if (main_task) parallel_close = nf90_close(ncid)
    call broadcast(parallel_close)

  end function parallel_close

!=======================================================================

  function parallel_create(path, cmode, ncid)

    implicit none
    integer :: cmode, ncid, parallel_create

    character(len=*) :: path

    ! begin
    if (main_task) parallel_create = nf90_create(path,cmode,ncid)
    call broadcast(parallel_create)
    call broadcast(ncid)

  end function parallel_create

!=======================================================================

  function parallel_def_dim(ncid, name, len, dimid)

    use netcdf
    implicit none
    integer :: dimid, len, ncid, parallel_def_dim
    type(parallel_type) :: parallel

    character(len=*) :: name

    ! begin
    if (main_task) parallel_def_dim = nf90_def_dim(ncid,name,len,dimid)
    call broadcast(parallel_def_dim)
    call broadcast(dimid)

  end function parallel_def_dim

!=======================================================================

  ! functions belonging to the parallel_def_var interface

  function parallel_def_var_dimids(ncid, name, xtype, dimids, varid)

    implicit none
    integer :: ncid, parallel_def_var_dimids, varid, xtype
    integer,dimension(:) :: dimids
    character(len=*) :: name

    ! begin
    if (main_task) parallel_def_var_dimids = &
         nf90_def_var(ncid,name,xtype,dimids,varid)
    call broadcast(parallel_def_var_dimids)
    call broadcast(varid)

  end function parallel_def_var_dimids


  function parallel_def_var_nodimids(ncid, name, xtype, varid)

    implicit none
    integer :: ncid, parallel_def_var_nodimids, varid, xtype
    character(len=*) :: name

    ! begin
    if (main_task) parallel_def_var_nodimids = &
         nf90_def_var(ncid,name,xtype,varid)
    call broadcast(parallel_def_var_nodimids)
    call broadcast(varid)

  end function parallel_def_var_nodimids

!=======================================================================

  function parallel_enddef(ncid)

    implicit none
    integer :: ncid,parallel_enddef

    ! begin
    if (main_task) parallel_enddef = nf90_enddef(ncid)
    call broadcast(parallel_enddef)

  end function parallel_enddef

!=======================================================================

  ! functions belonging to the parallel_get_att interface

  function parallel_get_att_character(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_get_att_character,varid
    character(len=*) :: name,values

    ! begin
    if (main_task) parallel_get_att_character = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_character)
    call broadcast(values)

  end function parallel_get_att_character


  function parallel_get_att_real4(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_get_att_real4,varid
    character(len=*) :: name
    real(sp) :: values

    ! begin
    if (main_task) parallel_get_att_real4 = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real4)
    call broadcast(values)

  end function parallel_get_att_real4


  function parallel_get_att_real4_1d(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_get_att_real4_1d, varid
    character(len=*) :: name
    real(sp),dimension(:) :: values

    ! begin
    if (main_task) parallel_get_att_real4_1d = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real4_1d)
    call broadcast(values)

  end function parallel_get_att_real4_1d


  function parallel_get_att_real8(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_get_att_real8, varid
    character(len=*) :: name
    real(dp) :: values

    ! begin
    if (main_task) parallel_get_att_real8 = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real8)
    call broadcast(values)

  end function parallel_get_att_real8


  function parallel_get_att_real8_1d(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_get_att_real8_1d, varid
    character(len=*) :: name
    real(dp),dimension(:) :: values

    ! begin
    if (main_task) parallel_get_att_real8_1d = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real8_1d)
    call broadcast(values)

  end function parallel_get_att_real8_1d

!=======================================================================

  ! functions belonging to the parallel_get_var interface

  !WHL - Added some parallel_get_var functions in analogy to parallel_put_var functions.
  !      The scalar and 1D functions broadcast values from main_task to other tasks.
  !      The 2D functions do not broadcast; they only bring a global array to main_task.
  !WHL Feb. 2025: Previously, these functions lacked an input 'start' argument.
  !      As a result, they always read the first time slice in each array.
  !      The 'start' argument is needed to get the correct time slice for time-varying fields.
  !      It is optional, since some calls do not need it.

  function parallel_get_var_integer(ncid, varid, values, start)

    implicit none
    integer :: ncid, parallel_get_var_integer, varid
    integer :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_get_var_integer = nf90_get_var(ncid,varid,values,start)
       else
          parallel_get_var_integer = nf90_get_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_get_var_integer)
    call broadcast(values)

  end function parallel_get_var_integer


  function parallel_get_var_real4(ncid, varid, values, start)

    implicit none
    integer :: ncid, parallel_get_var_real4, varid
    real(sp) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_get_var_real4 = nf90_get_var(ncid,varid,values,start)
       else
          parallel_get_var_real4 = nf90_get_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_get_var_real4)
    call broadcast(values)

  end function parallel_get_var_real4


  function parallel_get_var_real8(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_get_var_real8,varid
    real(dp) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_get_var_real8 = nf90_get_var(ncid,varid,values,start)
       else
          parallel_get_var_real8 = nf90_get_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_get_var_real8)
    call broadcast(values)

  end function parallel_get_var_real8


  function parallel_get_var_integer_1d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_get_var_integer_1d,varid
    integer,dimension(:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_get_var_integer_1d = nf90_get_var(ncid,varid,values,start)
       else
          parallel_get_var_integer_1d = nf90_get_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_get_var_integer_1d)
    call broadcast(values)

  end function parallel_get_var_integer_1d

  function parallel_get_var_real4_1d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_get_var_real4_1d,varid
    real(sp),dimension(:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_get_var_real4_1d = nf90_get_var(ncid,varid,values,start)
       else
          parallel_get_var_real4_1d = nf90_get_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_get_var_real4_1d)
    call broadcast(values)

  end function parallel_get_var_real4_1d

  function parallel_get_var_real8_1d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_get_var_real8_1d,varid
    real(dp),dimension(:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_get_var_real8_1d = nf90_get_var(ncid,varid,values,start)
       else
          parallel_get_var_real8_1d = nf90_get_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_get_var_real8_1d)
    call broadcast(values)

  end function parallel_get_var_real8_1d


  function parallel_get_var_integer_2d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_get_var_integer_2d,varid
    integer,dimension(:,:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_get_var_integer_2d = nf90_get_var(ncid,varid,values,start)
       else
          parallel_get_var_integer_2d = nf90_get_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_get_var_integer_2d)
!!    call broadcast(values)  ! no broadcast subroutine for 2D arrays

  end function parallel_get_var_integer_2d


  function parallel_get_var_real8_2d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_get_var_real8_2d,varid
    real(dp),dimension(:,:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_get_var_real8_2d = nf90_get_var(ncid,varid,values,start)
       else
          parallel_get_var_real8_2d = nf90_get_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_get_var_real8_2d)
!!    call broadcast(values)  ! no broadcast subroutine for 2D arrays

  end function parallel_get_var_real8_2d

!=======================================================================

  function parallel_inq_attname(ncid, varid, attnum, name)

    implicit none
    integer :: attnum, ncid, parallel_inq_attname, varid
    character(len=*) :: name

    ! begin
    if (main_task) parallel_inq_attname = &
         nf90_inq_attname(ncid,varid,attnum,name)
    call broadcast(parallel_inq_attname)
    call broadcast(name)

  end function parallel_inq_attname

!=======================================================================

  function parallel_inq_dimid(ncid, name, dimid)

    implicit none
    integer :: dimid, ncid, parallel_inq_dimid
    character(len=*) :: name

    ! begin
    if (main_task) parallel_inq_dimid = nf90_inq_dimid(ncid,name,dimid)
    call broadcast(parallel_inq_dimid)
    call broadcast(dimid)

  end function parallel_inq_dimid

!=======================================================================

  function parallel_inq_varid(ncid, name, varid)

    implicit none
    integer :: ncid,parallel_inq_varid,varid
    character(len=*) :: name

    ! begin
    if (main_task) parallel_inq_varid = nf90_inq_varid(ncid,name,varid)
    call broadcast(parallel_inq_varid)
    call broadcast(varid)

  end function parallel_inq_varid

!=======================================================================

  function parallel_inquire(ncid, nvariables)

    implicit none
    integer :: ncid,parallel_inquire,nvariables

    ! begin
    if (main_task) parallel_inquire = nf90_inquire(ncid,nvariables=nvariables)
    call broadcast(parallel_inquire)
    call broadcast(nvariables)

  end function parallel_inquire

!=======================================================================

  function parallel_inquire_dimension(ncid, dimid, name, len)

    implicit none
    integer :: dimid, ncid, parallel_inquire_dimension
    integer,optional :: len
    character(len=*),optional :: name
    
    integer :: l
    
    ! begin
    if (present(name)) then
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,name,len=l)
       call broadcast(name)
    else
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,len=l)
    end if
    call broadcast(parallel_inquire_dimension)
    if (present(len)) then
       call broadcast(l)
       len = l
    end if

  end function parallel_inquire_dimension

!=======================================================================

  function parallel_inquire_variable(ncid, varid, &
                                     name, ndims, dimids, natts)

    implicit none
    integer :: ncid, parallel_inquire_variable, varid
    integer,optional :: ndims, natts
    character(len=*),optional :: name
    integer,dimension(:),optional :: dimids

    integer :: nd,na

    ! begin
    if (present(name)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,name=name)
       call broadcast(parallel_inquire_variable)
       call broadcast(name)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (present(dimids)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,dimids=dimids)
       call broadcast(parallel_inquire_variable)
       call broadcast(dimids)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (main_task) parallel_inquire_variable = &
         nf90_inquire_variable(ncid,varid,ndims=nd,natts=na)
    call broadcast(parallel_inquire_variable)
    if (present(ndims)) then
       call broadcast(nd)
       ndims = nd
    end if
    if (present(natts)) then
       call broadcast(na)
       natts = na
    end if

  end function parallel_inquire_variable

!=======================================================================

  function parallel_open(path, mode, ncid)

    implicit none
    integer :: mode, ncid, parallel_open
    character(len=*) :: path

    ! begin
    if (main_task) parallel_open = nf90_open(path,mode,ncid)
    call broadcast(parallel_open)

  end function parallel_open

!=======================================================================

  ! functions belonging to the parallel_put_att interface

  function parallel_put_att_character(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_put_att_character, varid
    character(len=*) :: name, values

    ! begin
    if (main_task) parallel_put_att_character = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_character)

  end function parallel_put_att_character


  function parallel_put_att_integer(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_put_att_integer, varid
    character(len=*) :: name
    integer :: values

    ! begin
    if (main_task) parallel_put_att_integer = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_integer)

  end function parallel_put_att_integer


  function parallel_put_att_real4(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_put_att_real4,varid
    character(len=*) :: name
    real(sp) :: values

    ! begin
    if (main_task) parallel_put_att_real4 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4)

  end function parallel_put_att_real4


  function parallel_put_att_real4_1d(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_put_att_real4_1d,varid
    character(len=*) :: name
    real(sp),dimension(:) :: values

    ! begin
    if (main_task) parallel_put_att_real4_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4_1d)

  end function parallel_put_att_real4_1d


  function parallel_put_att_real8(ncid, varid, name, values)

    implicit none
    integer :: ncid, parallel_put_att_real8, varid
    character(len=*) :: name
    real(dp) :: values

    ! begin
    if (main_task) parallel_put_att_real8 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8)

  end function parallel_put_att_real8


  function parallel_put_att_real8_1d(ncid, varid, name, values)

    implicit none
    integer :: ncid,parallel_put_att_real8_1d,varid
    character(len=*) :: name
    real(dp),dimension(:) :: values

    ! begin
    if (main_task) parallel_put_att_real8_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8_1d)

  end function parallel_put_att_real8_1d

!=======================================================================

  ! functions belonging to the parallel_put_var interface

  function parallel_put_var_integer(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_integer,varid
    integer :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_integer = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_integer = nf90_put_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_put_var_integer)

  end function parallel_put_var_integer


  function parallel_put_var_integer_1d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_integer_1d,varid
    integer,dimension(:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_integer_1d = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_integer_1d = nf90_put_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_put_var_integer_1d)

  end function parallel_put_var_integer_1d


  function parallel_put_var_real4(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_real4,varid
    real(sp) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real4 = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real4 = nf90_put_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_put_var_real4)

  end function parallel_put_var_real4


  function parallel_put_var_real8(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_real8,varid
    real(dp) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real8 = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real8 = nf90_put_var(ncid,varid,values)
       endif
    endif
    call broadcast(parallel_put_var_real8)

  end function parallel_put_var_real8


  function parallel_put_var_real8_1d(ncid, varid, values, start)

    implicit none
    integer :: ncid,parallel_put_var_real8_1d,varid
    real(dp),dimension(:) :: values
    integer,dimension(:),optional :: start

    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values)
       end if
    end if
    call broadcast(parallel_put_var_real8_1d)

  end function parallel_put_var_real8_1d

!=======================================================================

  function parallel_redef(ncid)

    implicit none
    integer :: ncid,parallel_redef

    ! begin
    if (main_task) parallel_redef = nf90_redef(ncid)
    call broadcast(parallel_redef)

  end function parallel_redef

!=======================================================================

  function parallel_sync(ncid)

    implicit none
    integer :: ncid,parallel_sync

    ! begin
    if (main_task) parallel_sync = nf90_sync(ncid)
    call broadcast(parallel_sync)

  end function parallel_sync

!=======================================================================

end module cism_parallel_netcdf

!=======================================================================
