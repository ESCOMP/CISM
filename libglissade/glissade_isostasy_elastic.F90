!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_isostasy_elastic.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glissade_isostasy_elastic

  ! Code for an elastic lithosphere; typically combined with a relaxing asthenosphere
  !  Notes:
  !  * This module is based on the isostasy_elastic module from the old Glide code.
  !    It was copied to libglissade, renamed and restructured in July 2026.
  !  * The elastic lithosphere calculation is done on a single task.
  !    When running in parallel runs, data are gathered onto the main task for the computation,
  !    then scattered back to local processors.
  !    This procedure does not scale well, although is is manageable on a 4 km mesh as long as
  !    the update is done at a frequency of once every few decades or less.

  use glimmer_global, only : dp
  use glimmer_paramets, only: iulog
  use glide_types, only: isos_elastic
  use glimmer_log

  implicit none

  private
  public :: glissade_init_elastic, glissade_calc_elastic

  logical, parameter :: verbose_elastic = .true.  ! if true, print diagnostic messages

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine glissade_init_elastic(rbel, deltax)

    !> initialize elastic lithosphere calculations

    use glimmer_physcon, only : pi
    implicit none

    type(isos_elastic) :: rbel          !> structure holding elastic lithosphere data
    real(dp), intent(in) :: deltax      !> grid spacing

    ! local variables
    real(dp) :: a     ! radius of disk
    real(dp) :: r     ! distance from centre
    integer :: i,j

    real(dp), parameter :: r_lr = 6.0d0   ! influence of disk load at (0,0) is felt within a radius of r_lr*rbel_r

    ! calculate a so that a circle of radius a is equivalent to a square with size deltax
    a = deltax/sqrt(pi)

    ! initialise w
    call init_rbel(rbel, a)

    ! calculate size of operator
    rbel%wsize = int(r_lr*rbel%lr/deltax)

    ! allocate memory for operator
    allocate(rbel%w(0:rbel%wsize,0:rbel%wsize))

    ! calculating points within disk
    rbel%w(0,0) = rbel_iw(rbel,0.d0)
    r = deltax/rbel%lr
    rbel%w(0,1) = rbel_iw(rbel,r)
    rbel%w(1,0) = rbel%w(0,1)

    ! calculating points outside disk
    do j=0,rbel%wsize
       do i=2,rbel%wsize
          r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
          rbel%w(i,j) = rbel_ow(rbel,r)
       end do
    end do

    do j=2,rbel%wsize
       do i=0,1
          r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
          rbel%w(i,j) = rbel_ow(rbel,r)
       end do
    end do

    i=1
    j=1
    r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
    rbel%w(i,j) = rbel_ow(rbel,r)

#ifdef DEB_REBOUND
    open(1,file='w.dat',status='UNKNOWN')
    do j=0,rbel%wsize
       do i=0,rbel%wsize
          write(1,*) i,j,rbel%w(i,j)
       end do
    end do
    close(1)
#endif

  end subroutine glissade_init_elastic

!-------------------------------------------------------------------------

  subroutine glissade_calc_elastic(&
       rbel,                 &
       load_factors,         &
       load,                 &
       parallel)

    !> Calculate surface loading effect using elastic lithosphere approximation.
    !> Functionally equivalent to subroutine calc_elastic from Glimmer's original isostasy model.
    !> The main difference is that this subroutine uses a global gather and scatter to compute
    !>  the load for simulations on more than one task.

    use cism_parallel, only: this_rank, main_task, parallel_type, broadcast, &
         gather_var, scatter_var, parallel_halo, parallel_globalindex  !TODO - Remove scatter_var?

    implicit none

    ! input-output arguments
    type(isos_elastic) :: rbel                             !> structure holding elastic litho data
    real(dp), dimension(:,:), intent(in)  :: load_factors  !> load mass divided by mantle density
    real(dp), dimension(:,:), intent(out) :: load          !> loading effect due to load_factors

    type(parallel_type), intent(in) :: parallel            !> info for parallel communication

    ! local variables

    integer :: ewn, nsn                   ! grid dimensions on the local task; includes halo cells
    integer :: global_ewn, global_nsn     ! global grid dimensions

    integer :: i, j, n, m
    integer :: ig, jg                     ! global indices

    real(dp), dimension(:,:), allocatable :: &
         load_global,             & !> global version of the output 'load' array
         load_factors_global        !> global version of the input 'load_factors' array

    real(dp) :: local_sum_load, global_sum_load   !> diagnostic sums

    character(len=100) :: message

!    logical, parameter :: new_load_sum = .false.
    logical, parameter :: new_load_sum = .true.

    ! initialize

    ewn = size(load,1)
    nsn = size(load,2)
    global_ewn = parallel%global_ewn
    global_nsn = parallel%global_nsn

    load = 0.0d0

    if (verbose_elastic .and. main_task) then
       write(iulog,*) 'In glissade_calc_elastic'
       write(iulog,*) 'local ewn/nsn =', ewn, nsn
       write(iulog,*) 'global_ewn/nsn =', global_ewn, global_nsn
    endif

    ! Gather the local load_factors arrays onto the main task
    ! Note: global arrays are allocated in the subroutine
    call gather_var(load_factors, load_factors_global, parallel)


    if (new_load_sum) then

       if (verbose_elastic) then
          if (sum(load_factors_global) > 0.0d0) then
             write(iulog,*) 'my_task, sum(load_factors_global) =', &
                  this_rank, sum(load_factors_global)
          endif
          if (main_task) write(iulog,*) 'Allocate load_factors_global'
       endif

       ! allocate load_factors_global on tasks other than main
       if (.not.main_task) then
          if (allocated(load_factors_global)) deallocate(load_factors_global)
          allocate(load_factors_global(global_ewn,global_nsn))
       endif

       if (verbose_elastic .and. main_task) then
          write(iulog,*) 'Broadcasting ...'
       endif

       ! broadcast load_factors_global from main_task to all processors
       call broadcast(load_factors_global)

       if (verbose_elastic .and. main_task) then
          write(iulog,*) 'Broadcast done'
       endif

       if (sum(load_factors_global) == 0.0d0) then
          write(message,*) 'Error, calc_elastic, sum(load_factors_global) = 0, my_task =', this_rank
          call write_log(message)
       endif

       if (verbose_elastic .and. main_task) then
          write(iulog,*) 'Compute load locally on each task'
       endif

       do j = 1, nsn
          do i = 1, ewn
             call parallel_globalindex(i, j, ig, jg, parallel)

             ! Compute load terms by summing over cells in the radius of influence
             do n = max(1,jg-rbel%wsize), min(global_nsn,jg+rbel%wsize)
                do m = max(1,ig-rbel%wsize), min(global_ewn,ig+rbel%wsize)
                   load(i,j) = load(i,j) + load_factors_global(m,n) * rbel%w(abs(m-ig),abs(n-jg))
                end do
             end do

          enddo
       enddo

    else  ! do the sum on main_task and then scatter the solution

       ! allocate load_global
       allocate(load_global(global_ewn,global_nsn))
       load_global = 0.0d0

       if (main_task) then

          if (verbose_elastic) then
             write(iulog,*) 'Compute load on main_task'
          endif

          do j = 1, global_nsn

             if (verbose_elastic .and. main_task) then
                if (mod(j,100) == 0) write(iulog,*) 'j =', j   ! to see how fast the calculation is going
             endif

             do i = 1, global_ewn

                ! Compute load terms by summing over cells in the radius of influence
                do n = max(1,j-rbel%wsize), min(global_nsn,j+rbel%wsize)
                   do m = max(1,i-rbel%wsize), min(global_ewn,i+rbel%wsize)
                      load_global(i,j) = load_global(i,j) + load_factors_global(m,n) * rbel%w(abs(m-i),abs(n-j))
                   end do
                end do

             end do  ! i
          end do  ! j
       endif  ! main_task

       ! Scatter the load values back to local arrays
       ! Note: load_global is deallocated in the subroutine
       call scatter_var(load, load_global, parallel)

       ! scatter_var does not update the halo, so do an update here
       call parallel_halo(load, parallel)

    endif

    ! Deallocate global arrays
    deallocate(load_factors_global)

  end subroutine glissade_calc_elastic

!-------------------------------------------------------------------------

  subroutine init_rbel(rbel, a)

    !> initialize elastic lithosphere calculations

    use glimmer_physcon, only: rhom, grav
    use glissade_isostasy_kelvin, only: set_kelvin, dker, dkei, dber, dbei
    implicit none

    type(isos_elastic) :: rbel        !> structure holding elastic litho data
    real(dp), intent(in) :: a         !> radius of disk

    real(dp) :: dummy_a

    call set_kelvin(1.d-10,40)

    rbel%lr = (rbel%d/(rhom*grav))**0.25d0
    rbel%a  = a

    dummy_a = rbel%a/rbel%lr
    
    rbel%c1  =  dummy_a * dker(dummy_a)
    rbel%c2  = -dummy_a * dkei(dummy_a)
    rbel%cd3 =  dummy_a * dber(dummy_a)
    rbel%cd4 = -dummy_a * dbei(dummy_a)
    
  end subroutine init_rbel

!-------------------------------------------------------------------------

  function rbel_ow(rbel,r)

    use glissade_isostasy_kelvin, only: ker, kei

    !> calculate deflection outside disk

    implicit none
    real(dp) :: rbel_ow
    real(dp), intent(in) :: r          !> radius, r should be scaled with lr
    type(isos_elastic) :: rbel     !> structure holding elastic litho data
    
    rbel_ow = rbel%cd3*ker(r) + rbel%cd4*kei(r)

  end function rbel_ow

!-------------------------------------------------------------------------

  function rbel_iw(rbel,r)

    use glissade_isostasy_kelvin, only: ber, bei

    !> calculate deflection inside disk
    implicit none
    real(dp) :: rbel_iw
    real(dp), intent(in) :: r          !> radius, r should be scaled with lr
    type(isos_elastic) :: rbel         !> structure holding elastic litho data
    
    rbel_iw = 1.d0 + rbel%c1*ber(r) + rbel%c2*bei(r)

  end function rbel_iw

!-------------------------------------------------------------------------

end module glissade_isostasy_elastic

!-------------------------------------------------------------------------
