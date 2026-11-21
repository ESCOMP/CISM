!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_parallel_halo.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_parallel_halo

  ! Part of the cism_parallel module:
  !  halo updates on the scalar and staggered grids, and utilities that use them
  ! Other modules use these through cism_parallel.

  use glimmer_global, only : dp, sp, i8
  use glimmer_paramets, only: iulog
  use cism_parallel_core

  implicit none

  ! Debug and Verification Level
  integer, parameter :: DEBUG_LEVEL = 1
	! If > 0, then debug code executed.  Added for parallel_halo_verify()

  interface parallel_convert_nonhaloed_to_haloed
     module procedure parallel_convert_nonhaloed_to_haloed_real4_2d
     module procedure parallel_convert_nonhaloed_to_haloed_real8_2d
  end interface parallel_convert_nonhaloed_to_haloed
  
  interface parallel_halo
     module procedure parallel_halo_integer_2d
     module procedure parallel_halo_integer_3d
     module procedure parallel_halo_logical_2d
     module procedure parallel_halo_real4_2d
     module procedure parallel_halo_real8_2d
     module procedure parallel_halo_real8_3d
     module procedure parallel_halo_real8_4d
     module procedure parallel_halo_integer8_4d
  end interface

  interface parallel_halo_extrapolate
     module procedure parallel_halo_extrapolate_real8_1d
     module procedure parallel_halo_extrapolate_integer_2d
     module procedure parallel_halo_extrapolate_real8_2d
  end interface

  interface parallel_halo_tracers
     module procedure parallel_halo_tracers_real8_3d
     module procedure parallel_halo_tracers_real8_4d
  end interface

  interface parallel_halo_verify
     module procedure parallel_halo_verify_integer_2d
     module procedure parallel_halo_verify_real8_2d
     module procedure parallel_halo_verify_real8_3d
  end interface

  interface staggered_parallel_halo
     module procedure staggered_parallel_halo_integer_2d
     module procedure staggered_parallel_halo_integer_3d
     module procedure staggered_parallel_halo_real8_2d
     module procedure staggered_parallel_halo_real8_3d
     module procedure staggered_parallel_halo_real8_4d
  end interface

  interface staggered_parallel_halo_extrapolate
     module procedure staggered_parallel_halo_extrapolate_integer_2d
     module procedure staggered_parallel_halo_extrapolate_real8_2d
  end interface

contains

!=======================================================================

  ! subroutines belonging to the parallel_convert_nonhaloed_to_haloed interface

  subroutine parallel_convert_nonhaloed_to_haloed_real4_2d(input_no_halo, output_with_halo, parallel)

    ! Given an input array without halo cells, return an output array with halo cells
    real(sp),dimension(:,:), intent(in)  :: input_no_halo
    real(sp),dimension(:,:), intent(out) :: output_with_halo
    type(parallel_type) :: parallel
    
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn       &
         )

    if (size(input_no_halo,1) /= own_ewn .or. size(input_no_halo,2) /= own_nsn) then
       write(iulog,*) "Unexpected size for input_no_halo: ", &
            size(input_no_halo,1), size(input_no_halo,2)
       write(iulog,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_with_halo,1) /= local_ewn .or. size(output_with_halo,2) /= local_nsn) then
       write(iulog,*) "Unexpected size for output_with_halo: ", &
            size(output_with_halo,1), size(output_with_halo,2)
       write(iulog,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = &
         input_no_halo(1:own_ewn, 1:own_nsn)

    call parallel_halo(output_with_halo, parallel)
    
    end associate

  end subroutine parallel_convert_nonhaloed_to_haloed_real4_2d


  subroutine parallel_convert_nonhaloed_to_haloed_real8_2d(input_no_halo, output_with_halo, parallel)

    ! Given an input array without halo cells, return an output array with halo cells
    real(dp),dimension(:,:), intent(in)  :: input_no_halo
    real(dp),dimension(:,:), intent(out) :: output_with_halo
    type(parallel_type) :: parallel
    
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn       &
         )

    if (size(input_no_halo,1) /= own_ewn .or. size(input_no_halo,2) /= own_nsn) then
       write(iulog,*) "Unexpected size for input_no_halo: ", &
            size(input_no_halo,1), size(input_no_halo,2)
       write(iulog,*) "Expected size is: ", own_ewn, own_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    if (size(output_with_halo,1) /= local_ewn .or. size(output_with_halo,2) /= local_nsn) then
       write(iulog,*) "Unexpected size for output_with_halo: ", &
            size(output_with_halo,1), size(output_with_halo,2)
       write(iulog,*) "Expected size is: ", local_ewn, local_nsn
       call parallel_stop(__FILE__, __LINE__)
    end if

    output_with_halo(1+lhalo:local_ewn-uhalo, 1+lhalo:local_nsn-uhalo) = &
         input_no_halo(1:own_ewn, 1:own_nsn)

    call parallel_halo(output_with_halo, parallel)
    
    end associate

  end subroutine parallel_convert_nonhaloed_to_haloed_real8_2d

!=======================================================================

  subroutine parallel_global_edge_mask(global_edge_mask, parallel)

    ! Create a mask = 1 in locally owned cells at the edge of the global domain,
    ! = 0 elsewhere

    integer, dimension(:,:), intent(out) :: global_edge_mask
    type(parallel_type) :: parallel

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south)

    ! Check array dimensions

    ! unknown grid
    if (size(global_edge_mask,1)/=local_ewn .or. size(global_edge_mask,2)/=local_nsn) then
       write(iulog,*) "Unknown Grid: Size a=(", size(global_edge_mask,1), ",", size(global_edge_mask,2), &
            ") and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! Identify cells at the edge of the global domain

    global_edge_mask = 0

    if (this_rank >= east) then  ! at east edge of global domain
       global_edge_mask(local_ewn-uhalo,:) = 1
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       global_edge_mask(lhalo+1,:) = 1
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       global_edge_mask(:,local_nsn-uhalo) = 1
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       global_edge_mask(:,lhalo+1) = 1
    endif

    call parallel_halo(global_edge_mask, parallel)

    end associate

  end subroutine parallel_global_edge_mask

!=======================================================================

  ! subroutines belonging to the parallel_halo interface

  subroutine parallel_halo_integer_2d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel
    
    integer :: erequest, nrequest, srequest, wrequest, ierror
    integer,dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    integer,dimension(parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
       write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
      a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:) = 0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:) = 0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:) = 0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo) = 0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:) = 0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:) = 0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:) = 0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1) = 0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1) = 0
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1) = 0
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:) = 0
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:) = 0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_integer_2d


  subroutine parallel_halo_integer_3d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    integer,dimension(size(a,1), lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(size(a,1), uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(size(a,1), parallel%local_ewn, lhalo) :: nsend,srecv
    integer,dimension(size(a,1), parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,2)==local_ewn-1.and.size(a,3)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn.or.size(a,3)/=local_nsn) then
         write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:,:)

    nsend(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,local_nsn-uhalo+1:) = nrecv(:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,local_ewn-uhalo+1:,:) = 0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:lhalo,:) = 0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,local_nsn-uhalo+1:) = 0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:lhalo) = 0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,local_ewn-uhalo:,:) = 0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:lhalo+1,:) = 0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,local_nsn-uhalo:) = 0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:lhalo+1) = 0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:,:lhalo+1,:lhalo+1) = 0
       if (southeast_corner) a(:,local_ewn-lhalo:,:lhalo+1) = 0
       if (northeast_corner) a(:,local_ewn-lhalo:,local_nsn-lhalo:) = 0
       if (northwest_corner) a(:,:lhalo+1,local_nsn-lhalo:) = 0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_integer_3d


  subroutine parallel_halo_logical_2d(a, parallel)

    use mpi_mod
    implicit none
    logical,dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,nrequest,srequest,wrequest
    logical,dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    logical,dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    logical,dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    logical,dimension(parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
       write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_logical,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_logical,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_logical,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_logical,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_logical,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_logical,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_logical,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_logical,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                        ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:) = .false.
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:) = .false.
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:) = .false.
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo) = .false.
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:) = .false.
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:) = .false.
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:) = .false.
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1) = .false.
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1) = .false.
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1) = .false.
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:) = .false.
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:) = .false.

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_logical_2d


  subroutine parallel_halo_real4_2d(a, parallel)

    use mpi_mod
    implicit none
    real(sp),dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,nrequest,srequest,wrequest
    real(sp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(sp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(sp),dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    real(sp),dimension(parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
       write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real4,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real4,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real4,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real4,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real4,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real4,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real4,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real4,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                        ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:) = 0.
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:) = 0.
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:) = 0.
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo) = 0.
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:) = 0.
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:) = 0.
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:) = 0.
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1) = 0.
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1) = 0.
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1) = 0.
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:) = 0.
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:) = 0.

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_real4_2d


  subroutine parallel_halo_real8_2d(a, parallel, &
                                    periodic_offset_ew, periodic_offset_ns, zero_global_boundary_no_ice_bc)

    ! Added optional arguments for periodic offsets, to support ismip-hom test cases
    ! Also added an optional argument related to the no_ice BCs

    use mpi_mod
    implicit none
    real(dp),dimension(:,:) :: a
    type(parallel_type) :: parallel

    real(dp), intent(in), optional :: &
         periodic_offset_ew,     &! offset halo values by this amount
                                  ! if positive, the offset is positive for W halo, negative for E halo
         periodic_offset_ns       ! offset halo values by this amount
                                  ! if positive, the offset is positive for S halo, negative for N halo
    
    logical, intent(in), optional :: &
         zero_global_boundary_no_ice_bc  ! if true, then zero out values in grid cells adjacent
                                         ! to the global boundary when using no_ice BCs

    integer :: erequest,ierror,nrequest,srequest,wrequest
    real(dp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(parallel%local_ewn, uhalo) :: nrecv,ssend

    logical :: zero_global_boundary_no_ice  ! local version of zero_global_boundary_no_ice_bc

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    if (present(zero_global_boundary_no_ice_bc)) then
       zero_global_boundary_no_ice = zero_global_boundary_no_ice_bc
    else
       zero_global_boundary_no_ice = .true.
    endif

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) then
       write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ") and local_ewn and local_nsn = ", &
            local_ewn, ",", local_nsn
       call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    if (present(periodic_offset_ew)) then
       if (periodic_offset_ew /= 0.d0) then
          if (this_rank <= west) then   ! this proc lies at the west edge of the global domain
!             write(iulog,*) 'Offset at west edge: this_rank, west =', this_rank, west
             a(:lhalo,1+lhalo:local_nsn-uhalo) =   &
                a(:lhalo,1+lhalo:local_nsn-uhalo) + periodic_offset_ew
          endif
          if (this_rank >= east) then   ! this proc lies at the east edge of the global domain
!             write(iulog,*) 'Offset at east edge: this_rank, east =', this_rank, east
             a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) =    &
                a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) - periodic_offset_ew
          endif
       endif
    endif

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = nrecv(:,:)

    if (present(periodic_offset_ns)) then
       if (periodic_offset_ns /= 0.d0) then
          if (this_rank <= south) then  ! this proc lies at the south edge of the global domain
!             write(iulog,*) 'Offset at south edge: this_rank, south =', this_rank, south
             a(:,:lhalo) = a(:,:lhalo) + periodic_offset_ns
          endif
          if (this_rank >= north) then  ! this proc lies at the north edge of the global domain
!             write(iulog,*) 'Offset at north edge: this_rank, north =', this_rank, north
             a(:,local_nsn-uhalo+1:) = a(:,local_nsn-uhalo+1:) - periodic_offset_ns
          endif
       endif
    endif

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo) = 0.d0
       endif

    elseif (no_ice_bc) then

       if (zero_global_boundary_no_ice) then

          ! Set values to zero in cells adjacent to the global boundary;
          ! includes halo cells and one row of locally owned cells.

          if (this_rank >= east) then  ! at east edge of global domain
             a(local_ewn-uhalo:,:) = 0.d0
          endif

          if (this_rank <= west) then  ! at west edge of global domain
             a(:lhalo+1,:) = 0.d0
          endif

          if (this_rank >= north) then  ! at north edge of global domain
             a(:,local_nsn-uhalo:) = 0.d0
          endif

          if (this_rank <= south) then  ! at south edge of global domain
             a(:,:lhalo+1) = 0.d0
          endif

          ! Some interior blocks have a single cell at a corner of the global boundary.
          ! Set values in corner cells to zero, along with adjacent halo cells.
          if (southwest_corner) a(:lhalo+1,:lhalo+1) = 0.d0
          if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1) = 0.d0
          if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:) = 0.d0
          if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:) = 0.d0

       else  ! set values to zero in halo cells but not in locally owned cells

          if (this_rank >= east) then  ! at east edge of global domain
             a(local_ewn-uhalo+1:,:) = 0.d0
          endif

          if (this_rank <= west) then  ! at west edge of global domain
             a(:lhalo,:) = 0.d0
          endif

          if (this_rank >= north) then  ! at north edge of global domain
             a(:,local_nsn-uhalo+1:) = 0.d0
          endif

          if (this_rank <= south) then  ! at south edge of global domain
             a(:,:lhalo) = 0.d0
          endif

       endif   ! zero_global_boundary_no_ice

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_real8_2d


  subroutine parallel_halo_real8_3d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(dp),dimension(size(a,1), lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(size(a,1), uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(size(a,1), parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1), parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,2)==local_ewn-1.and.size(a,3)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn.or.size(a,3)/=local_nsn) then
         write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:,:)

    nsend(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,local_nsn-uhalo+1:) = nrecv(:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,local_ewn-uhalo+1:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:lhalo,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,local_nsn-uhalo+1:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:lhalo) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,local_ewn-uhalo:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:lhalo+1,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,local_nsn-uhalo:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:lhalo+1) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:,:lhalo+1,:lhalo+1) = 0.d0
       if (southeast_corner) a(:,local_ewn-lhalo:,:lhalo+1) = 0.d0
       if (northeast_corner) a(:,local_ewn-lhalo:,local_nsn-lhalo:) = 0.d0
       if (northwest_corner) a(:,:lhalo+1,local_nsn-lhalo:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_real8_3d


  subroutine parallel_halo_real8_4d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(dp),dimension(size(a,1), size(a,2), lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(size(a,1), size(a,2), uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(size(a,1), size(a,2), parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1), size(a,2), parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,3)==local_ewn-1.and.size(a,4)==local_nsn-1) return

    ! unknown grid
    if (size(a,3)/=local_ewn.or.size(a,4)/=local_nsn) then
         write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ",", size(a,4), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:,:) = &
         a(:,:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:,:,:)

    nsend(:,:,:,:) = a(:,:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:,:) = a(:,:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,:,:lhalo) = srecv(:,:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,:,local_nsn-uhalo+1:) = nrecv(:,:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,:,local_ewn-uhalo+1:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:,:lhalo,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,:,local_nsn-uhalo+1:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:,:lhalo) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,:,local_ewn-uhalo:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:,:lhalo+1,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,:,local_nsn-uhalo:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:,:lhalo+1) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:,:,:lhalo+1,:lhalo+1) = 0.d0
       if (southeast_corner) a(:,:,local_ewn-lhalo:,:lhalo+1) = 0.d0
       if (northeast_corner) a(:,:,local_ewn-lhalo:,local_nsn-lhalo:) = 0.d0
       if (northwest_corner) a(:,:,:lhalo+1,local_nsn-lhalo:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_real8_4d


  subroutine parallel_halo_integer8_4d(a, parallel)

    use mpi_mod
    implicit none
    integer(i8),dimension(:,:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    integer(i8),dimension(size(a,1), size(a,2), lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    integer(i8),dimension(size(a,1), size(a,2), uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    integer(i8),dimension(size(a,1), size(a,2), parallel%local_ewn, lhalo) :: nsend,srecv
    integer(i8),dimension(size(a,1), size(a,2), parallel%local_ewn, uhalo) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,3)==local_ewn-1.and.size(a,4)==local_nsn-1) return

    ! unknown grid
    if (size(a,3)/=local_ewn.or.size(a,4)/=local_nsn) then
         write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ",", size(a,4), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:,:) = &
         a(:,:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = erecv(:,:,:,:)

    nsend(:,:,:,:) = a(:,:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:,:) = a(:,:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,:,:lhalo) = srecv(:,:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,:,local_nsn-uhalo+1:) = nrecv(:,:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,:,local_ewn-uhalo+1:,:) = 0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:,:lhalo,:) = 0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,:,local_nsn-uhalo+1:) = 0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:,:lhalo) = 0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(:,:,local_ewn-uhalo:,:) = 0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:,:,:lhalo+1,:) = 0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,:,:,local_nsn-uhalo:) = 0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:,:,:lhalo+1) = 0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:,:,:lhalo+1,:lhalo+1) = 0
       if (southeast_corner) a(:,:,local_ewn-lhalo:,:lhalo+1) = 0
       if (northeast_corner) a(:,:,local_ewn-lhalo:,local_nsn-lhalo:) = 0
       if (northwest_corner) a(:,:,:lhalo+1,local_nsn-lhalo:) = 0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_integer8_4d

!=======================================================================

  ! subroutines for 1D halo updates

  subroutine parallel_halo_extrapolate_real8_1d(a, interval_in)

    !Note: Extrapolate a 1D real8 variable into halo cells to the east and west.
    !      Currently used only to compute halo values for grid cell coordinates.

    use mpi_mod
    implicit none
    real(dp),dimension(:) :: a
    type(parallel_type) :: parallel
    real(dp),intent(in), optional :: &
         interval_in  ! uniform difference between adjacent values, e.g. grid cell size dew or dns

    integer :: i
    real(dp) :: interval  ! local version of interval_in

    if (present(interval_in)) then
       interval = interval_in
    else
       interval = 0.0d0
    endif

    do i = 1, lhalo
       a(i) = a(lhalo+1) - interval*(lhalo+1-i)
    enddo

    do i = size(a)-uhalo+1, size(a)
       a(i) = a(size(a)-uhalo) + interval*(uhalo+i-size(a))
    enddo

  end subroutine parallel_halo_extrapolate_real8_1d

!=======================================================================

  ! subroutines belonging to the parallel_halo_extrapolate interface

  subroutine parallel_halo_extrapolate_integer_2d(a, parallel)

    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: i, j

    ! begin
    associate(  &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    ! Extrapolate the field into halo cells along the global boundary.

    ! First update the halos so that we are sure the interior halos are correct
    call parallel_halo(a, parallel)

! Useful for debugging small domains (the YYYY is just a tag for grepping the output,
!  particularly if you prepend the processor number, e.g. "0YYYY")
!  do j = 1, size(a,2)
!     write(iulog, "(i3, 'YYYY BEFORE row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    if (this_rank >= east) then  ! at east edge of global domain
       ! extrapolate eastward
       do i = size(a,1)-uhalo+1, size(a,1)
          a(i, :) = a(size(a,1)-uhalo, :)
       enddo
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! extrapolate westward
       do i = 1, lhalo
          a(i, :) = a(lhalo+1, :)
       enddo
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       ! extrapolate northward
       do j = size(a,2)-uhalo+1, size(a,2)
          a(:, j) = a(:, size(a,2)-uhalo)
       enddo
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! extrapolate southward
       do j = 1, lhalo
          a(:, j) = a(:, lhalo+1)
       enddo
    endif

! Useful for debugging small domains
!  do j = 1, size(a,2)
!     write(iulog, "(i3, 'YYYY AFTER  row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    end associate

  end subroutine parallel_halo_extrapolate_integer_2d


  subroutine parallel_halo_extrapolate_real8_2d(a, parallel)

    implicit none
    real(dp),dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: i, j

    ! begin
    associate(  &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    ! Extrapolate the field into halo cells along the global boundary.

    ! First update the halos so that we are sure the interior halos are correct
    call parallel_halo(a, parallel)

! Useful for debugging small domains (the YYYY is just a tag for grepping the output,
!  particularly if you prepend the processor number, e.g. "0YYYY")
!  do j = 1, size(a,2)
!     write(iulog, "(i3, 'YYYY BEFORE row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    if (this_rank >= east) then  ! at east edge of global domain
       ! extrapolate eastward
       do i = size(a,1)-uhalo+1, size(a,1)
          a(i, :) = a(size(a,1)-uhalo, :)
       enddo
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! extrapolate westward
       do i = 1, lhalo
          a(i, :) = a(lhalo+1, :)
       enddo
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       ! extrapolate northward
       do j = size(a,2)-uhalo+1, size(a,2)
          a(:, j) = a(:, size(a,2)-uhalo)
       enddo
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! extrapolate southward
       do j = 1, lhalo
          a(:, j) = a(:, lhalo+1)
       enddo
    endif

! Useful for debugging small domains
!  do j = 1, size(a,2)
!     write(iulog, "(i3, 'YYYY AFTER  row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    end associate

  end subroutine parallel_halo_extrapolate_real8_2d

!=======================================================================

  ! subroutines belonging to the parallel_halo_tracers interface

  subroutine parallel_halo_tracers_real8_3d(a, parallel)

    ! Custom halo routine for tracer arrays with dimension(nx,ny,ntracers)
    ! Will work for any 3D array with (nx,ny) in the first two slots

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: ierror, erequest, nrequest, srequest, wrequest
    real(dp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo, size(a,3)) :: esend,wrecv
    real(dp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo, size(a,3)) :: erecv,wsend
    real(dp),dimension(parallel%local_ewn, lhalo, size(a,3)) :: nsend,srecv
    real(dp),dimension(parallel%local_ewn, uhalo, size(a,3)) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn .or. size(a,2)/=local_nsn) then
         write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo,:)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo,:)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo,:) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo,:) = erecv(:,:,:)

    nsend(:,:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo,:)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,:)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo,:) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:,:) = nrecv(:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:,:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo,:) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:,:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1,:) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1,:) = 0.d0
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1,:) = 0.d0
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:,:) = 0.d0
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:,:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_tracers_real8_3d

  subroutine parallel_halo_tracers_real8_4d(a, parallel)

    ! Custom halo routine for tracer arrays with dimension(nx,ny,ntracers,nz)
    ! Will work for any 4D array with (nx,ny) in the first two slots

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:,:) :: a
    type(parallel_type) :: parallel

    integer :: ierror, erequest, nrequest, srequest, wrequest
    real(dp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo, size(a,3), size(a,4)) :: esend,wrecv
    real(dp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo, size(a,3), size(a,4)) :: erecv,wsend
    real(dp),dimension(parallel%local_ewn, lhalo, size(a,3), size(a,4)) :: nsend,srecv
    real(dp),dimension(parallel%local_ewn, uhalo, size(a,3), size(a,4)) :: nrecv,ssend

    ! begin
    associate(  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         southwest_corner  => parallel%southwest_corner,   &
         southeast_corner  => parallel%southeast_corner,   &
         northeast_corner  => parallel%northeast_corner,   &
         northwest_corner  => parallel%northwest_corner    &
         )

    ! staggered grid
    if (size(a,1)==local_ewn-1 .and. size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn .or. size(a,2)/=local_nsn) then
         write(iulog,*) "Unknown Grid: Size a=(", size(a,1), ",", size(a,2), ",", size(a,3), ",", size(a,4), ") &
                 &and local_ewn and local_nsn = ", local_ewn, ",", local_nsn
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo,:,:)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo,:,:)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo,:,:) = wrecv(:,:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo,:,:) = erecv(:,:,:,:)

    nsend(:,:,:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo,:,:)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,:,:)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:lhalo,:,:) = srecv(:,:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:,:,:) = nrecv(:,:,:,:)

    if (outflow_bc) then   ! set values in global halo to zero
                           ! interior halo cells should not be affected

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo+1:,:,:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo,:,:,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo+1:,:,:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo,:,:) = 0.d0
       endif

    elseif (no_ice_bc) then

       ! Set values to zero in cells adjacent to the global boundary;
       ! includes halo cells and one row of locally owned cells

       if (this_rank >= east) then  ! at east edge of global domain
          a(local_ewn-uhalo:,:,:,:) = 0.d0
       endif

       if (this_rank <= west) then  ! at west edge of global domain
          a(:lhalo+1,:,:,:) = 0.d0
       endif

       if (this_rank >= north) then  ! at north edge of global domain
          a(:,local_nsn-uhalo:,:,:) = 0.d0
       endif

       if (this_rank <= south) then  ! at south edge of global domain
          a(:,:lhalo+1,:,:) = 0.d0
       endif

       ! Some interior blocks have a single cell at a corner of the global boundary.
       ! Set values in corner cells to zero, along with adjacent halo cells.
       if (southwest_corner) a(:lhalo+1,:lhalo+1,:,:) = 0.d0
       if (southeast_corner) a(local_ewn-lhalo:,:lhalo+1,:,:) = 0.d0
       if (northeast_corner) a(local_ewn-lhalo:,local_nsn-lhalo:,:,:) = 0.d0
       if (northwest_corner) a(:lhalo+1,local_nsn-lhalo:,:,:) = 0.d0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine parallel_halo_tracers_real8_4d

!=======================================================================

  ! subroutines belonging to the parallel_halo_verify interface

  function parallel_halo_verify_integer_2d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel

    integer :: ierror, erequest, nrequest, srequest, wrequest
    integer,dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    integer,dimension(parallel%local_ewn, uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_integer_2d

    ! begin

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    ! ANY True if any value is true (LOGICAL)
    notverify_flag = ANY(a(:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:))

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:lhalo) /= srecv(:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,local_nsn-uhalo+1:) /= nrecv(:,:))

    ! if notverify_flag is TRUE, then there was some difference detected
    if (notverify_flag) then
         write(iulog,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_integer_2d = .NOT. notverify_flag  ! return if verified (True) or not verified (False)

    end associate

  end function parallel_halo_verify_integer_2d


  function parallel_halo_verify_real8_2d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:) :: a
    type(parallel_type) :: parallel
    
    integer :: ierror, erequest, nrequest, srequest, wrequest
    real(dp),dimension(lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(parallel%local_ewn, uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_real8_2d

    ! begin
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = &
         a(local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    notverify_flag = ANY(a(:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:))

    nsend(:,:) = a(:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:lhalo) /= srecv(:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,local_nsn-uhalo+1:) /= nrecv(:,:))

    if (notverify_flag) then
         write(iulog,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_real8_2d = .NOT. notverify_flag

    end associate

  end function parallel_halo_verify_real8_2d


  function parallel_halo_verify_real8_3d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:) :: a
    type(parallel_type) :: parallel
    
    integer :: ierror, erequest, nrequest, srequest, wrequest
    real(dp),dimension(size(a,1), lhalo, parallel%local_nsn-lhalo-uhalo) :: esend,wrecv
    real(dp),dimension(size(a,1), uhalo, parallel%local_nsn-lhalo-uhalo) :: erecv,wsend
    real(dp),dimension(size(a,1), parallel%local_ewn, lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1), parallel%local_ewn, uhalo) :: nrecv,ssend
    logical :: notverify_flag
    logical :: parallel_halo_verify_real8_3d

    ! begin
    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south         &
         )

    if (DEBUG_LEVEL <= 0) return

    ! staggered grid
    if (size(a,2)==local_ewn-1.and.size(a,3)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn.or.size(a,3)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = &
         a(:,local_ewn-uhalo-lhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    notverify_flag = ANY(a(:,:lhalo,1+lhalo:local_nsn-uhalo) /= wrecv(:,:,:))
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. &
      ANY(a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) /= erecv(:,:,:))

    nsend(:,:,:) = a(:,:,local_nsn-uhalo-lhalo+1:local_nsn-uhalo)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:,:lhalo) /= srecv(:,:,:))
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    notverify_flag = notverify_flag .OR. ANY(a(:,:,local_nsn-uhalo+1:) /= nrecv(:,:,:))

    if (notverify_flag) then
         write(iulog,*) "Halo Verify FAILED on processor ", this_rank
         ! call parallel_stop(__FILE__,__LINE__)
    endif

    parallel_halo_verify_real8_3d = .NOT. notverify_flag

    end associate

  end function parallel_halo_verify_real8_3d

!=======================================================================

  subroutine staggered_no_penetration_mask(umask, vmask, parallel)

    implicit none
    integer,dimension(:,:) :: umask, vmask  ! mask set to 1 wherever the outflow velocity should be zero
    type(parallel_type) :: parallel

    associate(  &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         east      => parallel%east,       &
         west      => parallel%west,       &
         north     => parallel%north,      &
         south     => parallel%south       &
         )

    ! initialize the no-penetration masks to 0
    umask(:,:) = 0
    vmask(:,:) = 0

    if (this_rank >= east) then  ! at east edge of global domain
       ! set u velocity mask = 1 at the east global boundary and vertices eastward
       umask(local_ewn-uhalo:,:) = 1
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! set u velocity mask = 1 at the west global boundary and vertices westward
       umask(:lhalo,:) = 1
    endif
    
    if (this_rank >= north) then  ! at north edge of global domain
       ! set v velocity mask = 1 at the north global boundary and vertices northward
       vmask(:,local_nsn-uhalo:) = 1
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! set v velocity mask = 1 at the south global boundary and vertices southward
       vmask(:,:lhalo) = 1
    endif

    call staggered_parallel_halo(umask, parallel)
    call staggered_parallel_halo(vmask, parallel)

    end associate

  end subroutine staggered_no_penetration_mask

!=======================================================================

  ! subroutines belonging to the staggered_parallel_halo interface

  !-----------------------------------------------------------------
  ! Comments on the staggered_parallel_halo subroutines:
  !
  ! The following subroutines implement a staggered grid halo update for integer and real(dp) arrays
  !  of various sizes.
  ! As the grid is staggered, the array 'a' is one smaller in both dimensions than an unstaggered array.
  !
  ! The grid is laid out from the SW, and the lower left corner is assigned to this_rank = 0.
  ! Its eastern neighbor is task_id = 1, proceeding rowwise and starting from the western edge.
  ! The South-most processes own one additional row of stagggered variables on the southern edge
  ! and have one less 'southern' halo row than other processes. Likewise, the West-most processes own one 
  ! additional column of staggered variables on the western edge and have one less 'western' halo column. 
  !
  ! The default BCs are periodic, in which case the southernmost and westernmost rows of the 
  ! global domain are filled with data from the northernmost and easternmost rows.
  !
  ! WHL: Comments on the implementation of outflow BC:
  !
  ! With outflow BCs, the southernmost and westernmost rows contain correct data (e.g., velocity
  ! at global boundaries) that should be passed to adjacent processors and should not be overwritten.
  ! This is implemented by defining variables called staggered_ilo, staggered_ihi, staggered_jlo
  ! and staggered_jhi, which denote the limits of locally owned staggered data on each processor.
  !
  ! For processors that include the western global boundary, staggered_ilo = staggered_lhalo; 
  !  for other processors (or for all processors when periodic_bc = T), staggered_ilo = staggered_lhalo+1. 
  ! Similarly, for processors that include the southern global boundary, staggered_jlo = staggered_lhalo; 
  !  for other processors (or for all processors when periodic_bc = T), staggered_jlo = staggered_lhalo+1.
  ! For all processors, staggered_ihi = (nx-1)-staggered_uhalo and staggered_jhi = (ny-1)-staggered_uhalo.
  ! Using these limits, we can correctly pass all and only the locally owned staggered data
  !  from each processor to its neighbors.
  ! Data are passed between the northern and southern global rows, and between the eastern and western rows,
  !  only when periodic_bc = T.
  !
  ! In addition to correctly passing staggered data along global boundaries, outflow BC imply
  !  that all staggered data are zero at points beyond the global boundary. The global boundary
  !  is delimited by:
  !    i = staggered_lhalo on the west boundary
  !    i = nx - 1 - staggered_uhalo on the east boundary
  !    j = staggered_lhalo on the south boundary
  !    j = ny - 1 - staggered_uhalo on the north boundary
  ! 
  ! TODO: Think about whether code simplifications are possible below using staggered_ilo/ihi/jlo/jhi.
  !
  !-----------------------------------------------------------------

  subroutine staggered_parallel_halo_integer_2d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a 2D integer field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    integer,dimension(staggered_lhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    integer,dimension(staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    integer,dimension(size(a,1),staggered_lhalo) :: nsend,srecv
    integer,dimension(size(a,1),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &       
         north     => parallel%north,          &    
         south     => parallel%south,          &      
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(iulog,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      wsend(:, 1:staggered_jhi-staggered_jlo+1) = &
           a(1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      esend(:, 1:staggered_jhi-staggered_jlo+1) = &
           a(size(a,1)-staggered_uhalo-staggered_lhalo+1:size(a,1)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(size(a,1)-staggered_uhalo+1:size(a,1), staggered_jlo:staggered_jhi) = &
          erecv(:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:) = &
        a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:) = &
        a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,size(a,2)-staggered_uhalo+1:size(a,2)) = nrecv(:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,1:staggered_lhalo) = srecv(:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary
    ! For no_ice BC, zero the field along and beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(1:staggered_lhalo-1, :) = 0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(size(a,1)-staggered_uhalo+1:size(a,1), :) = 0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, 1:staggered_lhalo-1) = 0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, size(a,2)-staggered_uhalo+1:size(a,2)) = 0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(1:staggered_lhalo, :) = 0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(size(a,1)-staggered_uhalo:size(a,1), :) = 0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, 1:staggered_lhalo) = 0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, size(a,2)-staggered_uhalo:size(a,2)) = 0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:staggered_lhalo, :staggered_lhalo) = 0
       if (southeast_corner) a(size(a,1)-staggered_uhalo:size(a,1), :staggered_lhalo) = 0
       if (northeast_corner) a(size(a,1)-staggered_uhalo:size(a,1), &
                               size(a,2)-staggered_uhalo:size(a,2)) = 0
       if (northwest_corner) a(:staggered_lhalo, size(a,2)-staggered_uhalo:size(a,2)) = 0

    endif   ! outflow or no_ice_bc

    end associate

  end subroutine staggered_parallel_halo_integer_2d


  subroutine staggered_parallel_halo_integer_3d(a, parallel)

    use mpi_mod
    implicit none
    integer,dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a 3D integer field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    integer,dimension(size(a,1), staggered_lhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    integer,dimension(size(a,1), staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    integer,dimension(size(a,1),size(a,2),staggered_lhalo) :: nsend,srecv
    integer,dimension(size(a,1),size(a,2),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &
         north     => parallel%north,          &
         south     => parallel%south,          &
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,2)/=local_ewn-1.or.size(a,3)/=local_nsn-1) then
         write(iulog,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      wsend(:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
       esend(:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:, size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(:, size(a,2)-staggered_uhalo+1:size(a,2), staggered_jlo:staggered_jhi) = &
          erecv(:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(:, 1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:,:) = &
        a(:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:,:) = &
        a(:,:,size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,:,size(a,3)-staggered_uhalo+1:size(a,3)) = nrecv(:,:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,:,1:staggered_lhalo) = srecv(:,:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then
          a(:, 1:staggered_lhalo-1, :) = 0
       endif

       if (this_rank >= east) then
          a(:, size(a,2)-staggered_uhalo+1:size(a,2), :) = 0
       endif

       if (this_rank <= south) then
          a(:, :, 1:staggered_lhalo-1) = 0
       endif

       if (this_rank >= north) then
          a(:, :, size(a,3)-staggered_uhalo+1:size(a,3)) = 0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(:, 1:staggered_lhalo, :) = 0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(:, size(a,2)-staggered_uhalo:size(a,2), :) = 0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, :, 1:staggered_lhalo) = 0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, :, size(a,3)-staggered_uhalo:size(a,3)) = 0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:, :staggered_lhalo, :staggered_lhalo) = 0
       if (southeast_corner) a(:, size(a,2)-staggered_uhalo:size(a,2), :staggered_lhalo) = 0
       if (northeast_corner) a(:, size(a,2)-staggered_uhalo:size(a,2), &
                               size(a,3)-staggered_uhalo:size(a,3)) = 0
       if (northwest_corner) a(:, :staggered_lhalo, size(a,3)-staggered_uhalo:size(a,3)) = 0

    endif   ! outflow or no_ice bc

    end associate

  end subroutine staggered_parallel_halo_integer_3d


  subroutine staggered_parallel_halo_real8_2d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a real(dp) 2D field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    real(dp),dimension(staggered_lhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    real(dp),dimension(staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    real(dp),dimension(size(a,1),staggered_lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &
         north     => parallel%north,          &
         south     => parallel%south,          &
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(iulog,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      wsend(:, 1:staggered_jhi-staggered_jlo+1) = &
           a(1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      esend(:, 1:staggered_jhi-staggered_jlo+1) = &
           a(size(a,1)-staggered_uhalo-staggered_lhalo+1:size(a,1)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(size(a,1)-staggered_uhalo+1:size(a,1), staggered_jlo:staggered_jhi) = &
          erecv(:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:) = &
        a(:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:) = &
        a(:,size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,size(a,2)-staggered_uhalo+1:size(a,2)) = nrecv(:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,1:staggered_lhalo) = srecv(:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(1:staggered_lhalo-1, :) = 0.0d0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(size(a,1)-staggered_uhalo+1:size(a,1), :) = 0.0d0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, 1:staggered_lhalo-1) = 0.0d0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, size(a,2)-staggered_uhalo+1:size(a,2)) = 0.0d0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(1:staggered_lhalo, :) = 0.0d0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(size(a,1)-staggered_uhalo:size(a,1), :) = 0.0d0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, 1:staggered_lhalo) = 0.0d0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, size(a,2)-staggered_uhalo:size(a,2)) = 0.0d0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:staggered_lhalo, :staggered_lhalo) = 0.0d0
       if (southeast_corner) a(size(a,1)-staggered_uhalo:size(a,1), :staggered_lhalo) = 0.0d0
       if (northeast_corner) a(size(a,1)-staggered_uhalo:size(a,1), &
                               size(a,2)-staggered_uhalo:size(a,2)) = 0.0d0
       if (northwest_corner) a(:staggered_lhalo, size(a,2)-staggered_uhalo:size(a,2)) = 0.0d0

    endif   ! outflow or no_ice_bc

    end associate

  end subroutine staggered_parallel_halo_real8_2d


  subroutine staggered_parallel_halo_real8_3d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a real(dp) 3D field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    real(dp),dimension(size(a,1), staggered_lhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    real(dp),dimension(size(a,1), staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    real(dp),dimension(size(a,1),size(a,2),staggered_lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1),size(a,2),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &
         north     => parallel%north,          &
         south     => parallel%south,          &
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,2)/=local_ewn-1 .or. size(a,3)/=local_nsn-1) then
         write(iulog,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    ! Send and receive east-west messages

    if (this_rank > west .or. periodic_bc) then
      wsend(:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
       esend(:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:, size(a,2)-staggered_uhalo-staggered_lhalo+1:size(a,2)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(:, size(a,2)-staggered_uhalo+1:size(a,2), staggered_jlo:staggered_jhi) = &
          erecv(:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(:, 1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    ! Send and receive north-south messages

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:,:) = &
        a(:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:,:) = &
        a(:,:,size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,:,size(a,3)-staggered_uhalo+1:size(a,3)) = nrecv(:,:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,:,1:staggered_lhalo) = srecv(:,:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then
          a(:, 1:staggered_lhalo-1, :) = 0.0d0
       endif

       if (this_rank >= east) then
          a(:, size(a,2)-staggered_uhalo+1:size(a,2), :) = 0.0d0
       endif

       if (this_rank <= south) then
          a(:, :, 1:staggered_lhalo-1) = 0.0d0
       endif

       if (this_rank >= north) then
          a(:, :, size(a,3)-staggered_uhalo+1:size(a,3)) = 0.0d0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(:, 1:staggered_lhalo, :) = 0.0d0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(:, size(a,2)-staggered_uhalo:size(a,2), :) = 0.0d0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:, :, 1:staggered_lhalo) = 0.0d0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:, :, size(a,3)-staggered_uhalo:size(a,3)) = 0.0d0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:, :staggered_lhalo, :staggered_lhalo) = 0.0d0
       if (southeast_corner) a(:, size(a,2)-staggered_uhalo:size(a,2), :staggered_lhalo) = 0.0d0
       if (northeast_corner) a(:, size(a,2)-staggered_uhalo:size(a,2), &
                               size(a,3)-staggered_uhalo:size(a,3)) = 0.0d0
       if (northwest_corner) a(:, :staggered_lhalo, size(a,3)-staggered_uhalo:size(a,3)) = 0.0d0

    endif   ! outflow or no_ice_bc

    end associate

  end subroutine staggered_parallel_halo_real8_3d


  subroutine staggered_parallel_halo_real8_4d(a, parallel)

    use mpi_mod
    implicit none
    real(dp),dimension(:,:,:,:) :: a
    type(parallel_type) :: parallel

    ! Implements a staggered grid halo update for a real(dp) 4D field

    integer :: ierror,nrequest,srequest,erequest,wrequest

    real(dp),dimension(size(a,1), size(a,2), staggered_lhalo,&
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: esend,wrecv
    real(dp),dimension(size(a,1), size(a,2), staggered_uhalo, &
         parallel%staggered_jhi-parallel%staggered_jlo+1) :: erecv,wsend
    real(dp),dimension(size(a,1),size(a,2),size(a,3),staggered_lhalo) :: nsend,srecv
    real(dp),dimension(size(a,1),size(a,2),size(a,3),staggered_uhalo) :: nrecv,ssend

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         local_ewn => parallel%local_ewn,      &
         local_nsn => parallel%local_nsn,      &
         east      => parallel%east,           &
         west      => parallel%west,           &
         north     => parallel%north,          &
         south     => parallel%south,          &
         staggered_jlo    => parallel%staggered_jlo,     &
         staggered_jhi    => parallel%staggered_jhi,     &
         southeast_corner => parallel%southeast_corner,  &
         southwest_corner => parallel%southwest_corner,  &
         northeast_corner => parallel%northeast_corner,  &
         northwest_corner => parallel%northwest_corner   &
         )

    ! Confirm staggered array
    if (size(a,3)/=local_ewn-1 .or. size(a,4)/=local_nsn-1) then
         write(iulog,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Prepost expected receives

    if (this_rank < east  .or. periodic_bc) then
      call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,comm,erequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,comm,wrequest,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,comm,nrequest,ierror)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,comm,srequest,ierror)
    endif

    if (this_rank > west .or. periodic_bc) then
      wsend(:,:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:,:, 1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1, &
             staggered_jlo:staggered_jhi)
      call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
       esend(:,:,:, 1:staggered_jhi-staggered_jlo+1) = &
           a(:,:, size(a,3)-staggered_uhalo-staggered_lhalo+1:size(a,3)-staggered_uhalo, &
             staggered_jlo:staggered_jhi)
      call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    endif

    if (this_rank < east .or. periodic_bc) then
      call mpi_wait(erequest,mpi_status_ignore,ierror)
      a(:,:, size(a,3)-staggered_uhalo+1:size(a,3), staggered_jlo:staggered_jhi) = &
          erecv(:,:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > west .or. periodic_bc) then
      call mpi_wait(wrequest,mpi_status_ignore,ierror)
      a(:,:, 1:staggered_lhalo, staggered_jlo:staggered_jhi) = &
          wrecv(:,:,:, 1:staggered_jhi-staggered_jlo+1)
    endif

    if (this_rank > south .or. periodic_bc) then
      ssend(:,:,:,:) = &
        a(:,:,:,1+staggered_lhalo:1+staggered_lhalo+staggered_uhalo-1)
      call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      nsend(:,:,:,:) = &
        a(:,:,:,size(a,4)-staggered_uhalo-staggered_lhalo+1:size(a,4)-staggered_uhalo)
      call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    endif

    if (this_rank < north .or. periodic_bc) then
      call mpi_wait(nrequest,mpi_status_ignore,ierror)
      a(:,:,:,size(a,4)-staggered_uhalo+1:size(a,4)) = nrecv(:,:,:,:)
    endif

    if (this_rank > south .or. periodic_bc) then
      call mpi_wait(srequest,mpi_status_ignore,ierror)
      a(:,:,:,1:staggered_lhalo) = srecv(:,:,:,:)
    endif

    ! For outflow BC, zero the field beyond the global boundary

    if (outflow_bc) then

       if (this_rank <= west) then
          a(:,:, 1:staggered_lhalo-1, :) = 0.0d0
       endif

       if (this_rank >= east) then
          a(:,:, size(a,3)-staggered_uhalo+1:size(a,3), :) = 0.0d0
       endif

       if (this_rank <= south) then
          a(:,:, :, 1:staggered_lhalo-1) = 0.0d0
       endif

       if (this_rank >= north) then
          a(:,:, :, size(a,4)-staggered_uhalo+1:size(a,4)) = 0.0d0
       endif

    elseif (no_ice_bc) then

       if (this_rank <= west) then  ! west edge of global domain
          a(:,:, 1:staggered_lhalo, :) = 0.0d0
       endif

       if (this_rank >= east) then  ! east edge of global domain
          a(:,:, size(a,3)-staggered_uhalo:size(a,3), :) = 0.0d0
       endif

       if (this_rank <= south) then  ! south edge of global domain
          a(:,:, :, 1:staggered_lhalo) = 0.0d0
       endif

       if (this_rank >= north) then  ! north edge of global domain
          a(:,:, :, size(a,4)-staggered_uhalo:size(a,4)) = 0.0d0
       endif

       ! Some interior blocks have a single vertex at a corner of the global boundary.
       ! Set values in corner vertices to zero, along with adjacent halo cells.
       if (southwest_corner) a(:,:, :staggered_lhalo, :staggered_lhalo) = 0.0d0
       if (southeast_corner) a(:,:, size(a,3)-staggered_uhalo:size(a,3), :staggered_lhalo) = 0.0d0
       if (northeast_corner) a(:,:, size(a,3)-staggered_uhalo:size(a,3), &
                               size(a,4)-staggered_uhalo:size(a,4)) = 0.0d0
       if (northwest_corner) a(:,:, :staggered_lhalo, size(a,4)-staggered_uhalo:size(a,4)) = 0.0d0

    endif   ! outflow or no_ice_bc

    end associate

  end subroutine staggered_parallel_halo_real8_4d

!=======================================================================

  ! subroutines belonging to the staggered_parallel_halo_extrapolate interface

  subroutine staggered_parallel_halo_extrapolate_integer_2d(a, parallel)

    implicit none
    integer,dimension(:,:) :: a
    type(parallel_type) :: parallel
    integer :: i, j

    ! begin

    associate(  &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         east   => parallel%east,    &
         west   => parallel%west,    &
         north  => parallel%north,   &
         south  => parallel%south    &
       )

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(iulog,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Extrapolate the staggered field into halo cells along the global boundary.
    ! Currently this is used only for kinbcmask.
    ! Note: The extrapolation region includes locally owned cells along
    !       the north and east boundaries of the global domain.

    ! First update the halos so that we are sure the interior halos are correct
    call staggered_parallel_halo(a, parallel)

    ! MJH Note: Modified code to now copy entire east and west columns rather than
    !  just the owned cells in those columns.  This avoids having the halos have
    !  potentially wrong information (i.e., a few cells in the corner don't get extrapolated into)

    if (this_rank >= east) then  ! at east edge of global domain
       ! extrapolate eastward
       do i = size(a,1)-staggered_uhalo, size(a,1)
          a(i, :) = a(size(a,1)-staggered_uhalo-1, :)
       enddo
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! extrapolate westward
       do i = 1, staggered_lhalo
          a(i, :) = a(staggered_lhalo+1, :)
       enddo
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       ! extrapolate northward
       do j = size(a,2)-staggered_uhalo, size(a,2)
          a(:, j) = a(:, size(a,2)-staggered_uhalo-1)
       enddo
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! extrapolate southward
       do j = 1, staggered_lhalo
          a(:, j) = a(:, staggered_lhalo+1)
       enddo
    endif

    end associate

  end subroutine staggered_parallel_halo_extrapolate_integer_2d


  subroutine staggered_parallel_halo_extrapolate_real8_2d(a, parallel)

    implicit none
    real(dp), dimension(:,:) :: a
    type(parallel_type) :: parallel
    integer :: i, j

    ! begin

    associate(  &
         local_ewn => parallel%local_ewn,  &
         local_nsn => parallel%local_nsn,  &
         east  => parallel%east,   &
         west  => parallel%west,   &
         north => parallel%north,  &
         south => parallel%south   &
       )

    ! Confirm staggered array
    if (size(a,1)/=local_ewn-1 .or. size(a,2)/=local_nsn-1) then
         write(iulog,*) "staggered_parallel_halo() requires staggered arrays."
         call parallel_stop(__FILE__,__LINE__)
    endif

    ! Extrapolate the staggered field into halo cells along the global boundary.
    ! Currently this is used only for kinbcmask.
    ! Note: The extrapolation region includes locally owned cells along
    !       the north and east boundaries of the global domain.

    ! First update the halos so that we are sure the interior halos are correct
    call staggered_parallel_halo(a, parallel)

    ! MJH Note: Modified code to now copy entire east and west columns rather than
    !  just the owned cells in those columns.  This avoids having the halos have
    !  potentially wrong information (i.e., a few cells in the corner don't get extrapolated into)

! Useful for debugging small domains (the YYYY is just a tag for grepping the output, particularly if you prepend the processor number, e.g. "0YYYY")
!  do j = 1, size(a,2)
!     write(iulog, "(i3, 'YYYY BEFORE row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo

    if (this_rank >= east) then  ! at east edge of global domain
       ! extrapolate eastward
       do i = size(a,1)-staggered_uhalo, size(a,1)
          a(i, :) = a(size(a,1)-staggered_uhalo-1, :)
       enddo
    endif

    if (this_rank <= west) then  ! at west edge of global domain
       ! extrapolate westward
       do i = 1, staggered_lhalo
          a(i, :) = a(staggered_lhalo+1, :)
       enddo
    endif

    if (this_rank >= north) then  ! at north edge of global domain
       ! extrapolate northward
       do j = size(a,2)-staggered_uhalo, size(a,2)
          a(:, j) = a(:, size(a,2)-staggered_uhalo-1)
       enddo
    endif

    if (this_rank <= south) then  ! at south edge of global domain
       ! extrapolate southward
       do j = 1, staggered_lhalo
          a(:, j) = a(:, staggered_lhalo+1)
       enddo
    endif

! Useful for debugging small domains
!  do j = 1, size(a,2)
!     write(iulog, "(i3, 'YYYY AFTER  row ', i3, 1000e9.2)")  this_rank, j, a(:,j)
!  enddo
    end associate

  end subroutine staggered_parallel_halo_extrapolate_real8_2d

!=======================================================================

end module cism_parallel_halo

!=======================================================================
