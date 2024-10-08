!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_noism.F90 - part of the Community Ice Sheet Model (CISM)  
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
!
! This module contains routines for setting ice velocity to zero 
!
! It is called with whichdycore = DYCORE_GLISSADE, which_ho_approx = HO_APPROX_NOISM. 
!
! Author: Heiko Goelzer
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_velo_noism

    use glimmer_global, only: dp

    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo, &
         parallel_halo, staggered_parallel_halo

    implicit none

    private
    public :: glissade_velo_noism_solve

  contains

!****************************************************************************

  subroutine glissade_velo_noism_solve(model,                &
                                     nx,     ny,     nz)

    !TODO - Remove nx, ny, nz from argument list?
    !       Would then have to allocate some local arrays.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    integer, intent(in) ::   &
       nx, ny,               &  ! number of grid cells in each direction
       nz                       ! number of vertical levels where velocity is computed
                                ! (same as model%general%upn)

    !----------------------------------------------------------------
    ! Local variables and pointers set to components of model derived type 
    !----------------------------------------------------------------

    real(dp), dimension(:,:,:), pointer ::  &
       uvel, vvel               ! velocity components (m/yr)

    type(parallel_type) :: &
         parallel               ! info for parallel communication

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: itest, jtest   ! coordinates of diagnostic point
    integer :: rtest          ! task number for processor containing diagnostic point

    integer :: i, j, k

    !--------------------------------------------------------
    ! Assign local pointers and variables to derived type components
    !--------------------------------------------------------

    parallel = model%parallel

    uvel     => model%velocity%uvel(:,:,:)
    vvel     => model%velocity%vvel(:,:,:)

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    do k = nz, 1, -1

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             uvel(k,i,j) = 0.d0
             vvel(k,i,j) = 0.d0
          enddo
       enddo
       
    enddo           ! k

    call staggered_parallel_halo(uvel, parallel)
    call staggered_parallel_halo(vvel, parallel)


  end subroutine glissade_velo_noism_solve

!*********************************************************************


!*********************************************************************

  end module glissade_velo_noism

!*********************************************************************
