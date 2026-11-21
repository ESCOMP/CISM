!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_parallel_broadcast.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_parallel_broadcast

  ! Part of the cism_parallel module:
  !  broadcasts from the main task to all tasks (interface broadcast)
  ! Other modules use these through cism_parallel.

  use glimmer_global, only : dp, sp
  use cism_parallel_core

  implicit none

  interface broadcast
     module procedure broadcast_character
     module procedure broadcast_integer
     module procedure broadcast_integer_1d
     module procedure broadcast_integer_2d
     module procedure broadcast_logical
     module procedure broadcast_logical_1d
     module procedure broadcast_real4
     module procedure broadcast_real4_1d
     module procedure broadcast_real8     
     module procedure broadcast_real8_1d
     module procedure broadcast_real8_2d
  end interface

contains

!=======================================================================

  ! subroutines belonging to the broadcast interface

  subroutine broadcast_character(c, proc)

    use mpi_mod
    implicit none
    character(len=*) :: c
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror, n
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    n = len(c)
    call mpi_bcast(c,n,mpi_character,source,comm,ierror)

  end subroutine broadcast_character


  subroutine broadcast_integer(i, proc)

    use mpi_mod
    implicit none
    integer :: i
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(i,1,mpi_integer,source,comm,ierror)

  end subroutine broadcast_integer


  subroutine broadcast_integer_1d(a, proc)

    use mpi_mod
    implicit none
    integer,dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(a,size(a),mpi_integer,source,comm,ierror)

  end subroutine broadcast_integer_1d


  subroutine broadcast_logical(l, proc)

    use mpi_mod
    implicit none
    logical :: l
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(l,1,mpi_logical,source,comm,ierror)

  end subroutine broadcast_logical


  subroutine broadcast_logical_1d(l, proc)

    use mpi_mod
    implicit none
    logical,dimension(:) :: l
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(l,size(l),mpi_logical,source,comm,ierror)

  end subroutine broadcast_logical_1d


  subroutine broadcast_real4(r, proc)

    use mpi_mod
    implicit none
    real(sp) :: r
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(r,1,mpi_real4,source,comm,ierror)

  end subroutine broadcast_real4


  subroutine broadcast_real4_1d(a, proc)

    use mpi_mod
    implicit none
    real(sp),dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(a,size(a),mpi_real4,source,comm,ierror)

  end subroutine broadcast_real4_1d


  subroutine broadcast_real8(r, proc)

    use mpi_mod
    implicit none
    real(dp) :: r
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(r,1,mpi_real8,source,comm,ierror)

  end subroutine broadcast_real8


  subroutine broadcast_real8_1d(a, proc)

    use mpi_mod
    implicit none
    real(dp),dimension(:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    call mpi_bcast(a,size(a),mpi_real8,source,comm,ierror)

  end subroutine broadcast_real8_1d


  subroutine broadcast_integer_2d(a, proc)

    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    !    call mpi_bcast(a,size(a),mpi_real8,source,comm,ierror)
    call mpi_bcast(a,size(a),mpi_integer,source,comm,ierror)

  end subroutine broadcast_integer_2d


  subroutine broadcast_real8_2d(a, proc)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:) :: a
    integer, intent(in), optional :: proc  ! optional argument indicating which processor to broadcast from

    integer :: ierror
    integer :: source ! local variable indicating which processor to broadcast from

    ! begin
    if (present(proc)) then
       source = proc
    else
       source = main_rank
    endif
    !    call mpi_bcast(a,size(a),mpi_real8,source,comm,ierror)
    call mpi_bcast(a,size(a),mpi_real8,source,comm,ierror)

  end subroutine broadcast_real8_2d

!=======================================================================

end module cism_parallel_broadcast

!=======================================================================
