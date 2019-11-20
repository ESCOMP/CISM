!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_basal_water.F90 - part of the Community Ice Sheet Model (CISM)  
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

!TODO - Support Jesse's water-routing code (or something similar) in parallel?
!       Currently supported only for serial Glide runs, in module glide_bwater.F90

module glissade_basal_water

   use glimmer_global, only: dp
   use glide_types

   implicit none

   private
   public :: glissade_basal_water_init, glissade_calcbwat

contains

  subroutine glissade_basal_water_init(model)

    use glimmer_paramets, only: thk0

    type(glide_global_type) :: model

    select case (model%options%which_ho_bwat)

    ! HO_BWAT_NONE:         basal water depth = 0
    ! HO_BWAT_CONSTANT:     basal water depth = prescribed constant
    ! HO_BWAT_LOCAL_TILL:   local basal till model with prescribed drainage rate
    ! HO_BWAT_SHAKTI:       SHAKTI basal hydrology model

    case(HO_BWAT_CONSTANT)

       ! Set a constant water thickness where ice is present
       where (model%geometry%thck > model%numerics%thklim)
          model%temper%bwat(:,:) = model%basal_physics%const_bwat / thk0
       elsewhere
          model%temper%bwat(:,:) = 0.0d0
       endwhere

    case default

       ! currently nothing to do for other options

    end select

  end subroutine glissade_basal_water_init


  subroutine glissade_calcbwat(which_ho_bwat,    &
                               basal_physics,    &
                               dt,               &
                               thck,             &
                               thklim,           &
                               bmlt,             &
                               bwat)

    ! Driver for updating basal water
    ! Note: This subroutine assumes SI units.
    ! Currently, only a few simple options are supported.

    use glimmer_physcon, only: rhow, scyr
    use glide_types

    !WHL - debug
    use parallel

    integer, intent(in) :: &
         which_ho_bwat     !> basal water options

    type(glide_basal_physics), intent(in)   :: basal_physics      ! basal physics object

    real(dp), intent(in) :: &
         dt,             & !> time step (s) 
         thklim            !> threshold for dynamically active ice (m)

    real(dp), dimension(:,:), intent(in) ::  &
         thck,           & !> ice thickness (m)
         bmlt              !> basal melt rate (m/s of ice)

    real(dp), dimension(:,:), intent(inout) ::  &
         bwat              !> basal water depth (m of water)

    ! local variables

    integer :: nx, ny    ! horizontal grid dimensions
    integer :: i, j

    real(dp) ::  &
         dbwat_dt        ! rate of change of bwat (m/s of water)

    select case (which_ho_bwat)

    ! HO_BWAT_NONE:         basal water depth = 0
    ! HO_BWAT_CONSTANT:     basal water depth = prescribed constant
    ! HO_BWAT_LOCAL_TILL:   local basal till model with prescribed drainage rate
    ! HO_BWAT_SHAKTI:       SHAKTI basal hydrology model

    case(HO_BWAT_NONE)

       bwat(:,:) = 0.0d0

    case(HO_BWAT_CONSTANT)

       ! Use a constant water thickness where ice is present, to force Tbed = Tpmp
       where (thck > thklim)
          bwat(:,:) = basal_physics%const_bwat
       elsewhere
          bwat(:,:) = 0.0d0
       endwhere

     case(HO_BWAT_LOCAL_TILL)

        nx = size(bwat,1)
        ny = size(bwat,2)

        do j = 1, ny
           do i = 1, nx

              ! compute new bwat, given source (bmlt) and sink (drainage)
              ! Note: bmlt > 0 for ice melting. Freeze-on will reduce bwat.
              dbwat_dt = bmlt(i,j)*rhoi/rhow - basal_physics%c_drainage/scyr  ! convert c_drainage from m/yr to m/s
              bwat(i,j) = bwat(i,j) + dbwat_dt*dt

              ! limit to the range [0, bwat_till_max]
              bwat(i,j) = min(bwat(i,j), basal_physics%bwat_till_max)
              bwat(i,j) = max(bwat(i,j), 0.0d0)

           enddo
        enddo

    end select

  end subroutine glissade_calcbwat
 
!---------------------------------------------------------------------------------
! ***************************************************
! SHAKTI subglacial hydrology model
! 2019 - ANS
! Based on the equations found in Sommers et al. (2018), GMD
!
! Note: SI units are used throughout this subroutine
! ***************************************************
!---------------------------------------------------------------------------------
 
  subroutine glissade_shakti(model)
!                             which_ho_bwat,            &
!                             basal_physics,            &
!                             dt,                       &
!                             ewn,            nsn,      &
!                             thck,           topg,     &
!                             flwa,           flwn,     &
!                             uvel,           vvel,     &
!                             bmlt,                     &
!                             head,                     &
!                             gap_height,               &
!                             effective_pressure,       &
!                             meltwater_input,          &
!                             moulin_input,             &
!                             reynolds,                 &
!                             head_gradient_mask_east,  &
!                             head_gradient_mask_north, &
!                             ice_mask_hydro,           &
!                             model,                    &
!                             bwat)

    use glimmer_physcon, only: rhoi, rhow, grav, c_t, c_w, nu_water 
    use parallel
    use glissade_velo_higher_pcg, only: pcg_solver_standard_2d_scalar
    !----------------------------------------------------
    ! SHAKTI - Input/output arguments
    ! ---------------------------------------------------
!    integer, intent(in) :: &
!         ewn, nsn          & ! grid dimensions
!         which_ho_bwat     & ! which basal option SHAKTI=3
!         dt                & ! time step (units?)
!         flwa              & ! flow law Parameter (Pa^-3 s^-1, I think) - ***check the name and units ***
!         flwn              & ! flow law exponent (unitless) -***check the name***

!    real(dp), dimension(:,:), intent(in) :: &
!         thck,            &! ice thickness (m)
!         topg,            &! bed topography (m)
!         meltwater_input  &! distributed meltwater input (m/s)
!         moulin_input     &! moulin point meltwater input (m^3/s)
!         bump_height      &! typical subglacial bump height (m)
!         bump_spacing     &! typical subglacial bump spacing (m)
!         uvel             &! sliding velocity (m/s?) ***only need basal velocity, also use vvel???***
!         vvel             &! sliding velocity (m/s?)
!         geothermal       &! geothermal heat flux (W/m2) - ***check the name of this***

!    real(dp), dimension(:,:), intent(in) :: &
!         ice_mask_hydro                &! ice mask (only use SHAKTI under ice)
!         head_gradient_mask_east       &! mask to set head gradient b.c. on east edge
!         head_gradient_mask_north      &! mask to set head gradient b.c. on north edge

!    real(dp), dimension(:,:), intent(inout) :: &
!         head             &! hydraulic head (m)
!         gap_height       &! subglacial gap height (m)
!         reynolds         &! Reynolds number (unitless)
!         bmlt             &! basal melt rate (units?) - will be altered by SHAKTI
!         bwat             &! basal water - do we need this?

!    real(dp), dimension(:,:), intent(out) :: &
!         effective pressure         &! effective pressure (Pa)

    type(glide_global_type), intent(inout) :: model       ! model instance

    integer, dimension(-1:1,-1:1) :: indxA                ! indices for cell and its neighbors, from 1 to 5
    real(dp), dimension(:,:,:), allocatable :: Ah         ! test matrix
    real(dp), dimension(:,:), allocatable ::   bh         ! rhs
    real(dp), dimension(:,:), allocatable ::   xh         ! solution
    logical, dimension(:,:), allocatable :: active_cell   ! true for active cells with a solution

         
    ! ----------------------------------------------
    ! SHAKTI - Local variables
    ! ----------------------------------------------
    ! *** Right now, from here down is based on Bill's glissade_test_matrix code ***
    ! *** Need to declare appropriate local SHAKTI variables ***

    integer :: nx, ny               ! grid dimensions
    integer :: itest, jtest, rtest  ! coordinates of diagnostic point
    real(dp) :: err                 ! solution error
    integer :: niters               ! number of iterations in solver
    integer :: i, j, m
    integer :: iglobal, jglobal
    integer :: K                    ! constant, uniform transmissivity (***for testing only!***)
    integer :: dx, dy               ! grid spacing

    nx = model%general%ewn
    ny = model%general%nsn
 
    K = 1.0d0
 
    dx = 1.0d0
    dy = 1.0d0

    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local
    rtest = model%numerics%rdiag_local

    if (this_rank == rtest) then
       print*, ' '
       print*, 'In glissade_test_matrix: itest, jtest, rtest =', itest, jtest, rtest
    endif 

    ! Set the indexing for the sparse matrix
    ! Indexing is as follows, with diagonal entries given an index of 3:
    !
    !        0  5  0
    !        2  3  4
    !        0  1  0

    indxA(:,:) = 0
    indxA( 0,-1) = 1  ! i =  0, j = -1
    indxA(-1, 0) = 2  ! i = -1, j =  0
    indxA( 0, 0) = 3  ! i =  0, j =  0
    indxA( 1, 0) = 4  ! i =  1, j =  0
    indxA( 0, 1) = 5  ! i =  0, j =  1

    ! Initialize
    allocate(active_cell(nx,ny))
    allocate(Ah(nx,ny,5))
    allocate(bh(nx,ny))
    allocate(xh(nx,ny))

    active_cell(:,:) = .false.
    Ah(:,:,:) = 0.0d0
    bh(:,:) = 0.0d0
    xh(:,:) = 0.0d0

    ! Set up a simple block-diagonal matrix, with 4 on the main diagonal and -1 on the off-diagonals.
    ! Set rhs = 1 for cells on the global boundary, with rhs = 0 elsewhere.
    ! This gives a solution of x = 1 in all cells on the domain.

    if (tasks == 1) then

       ! loop over locally owned cells
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             active_cell(i,j) = .true.
             Ah(i,j,3) = 4.0*K/dx**2d0 
             if (i == nhalo+1) then          ! westernmost cell with nonzero x; subdiag term = 0
                Ah(i,j,3) = 1.0d0            !     Dirichlet b.c. at outflow, h=zb 
                bh(i,j) = 0.0d0              !     rhs should be set to bed elevation (zero just for simple test) 
             elseif (i == nx-nhalo) then     ! easternmost cell with nonzero x; supdiag term = 0
                Ah(i,j,2) = -1.0*K/dx**2d0            !    Zero-flux Neumann b.c.
                bh(i,j) = 0.0d0  
             else                            ! interior cell
                Ah(i,j,2) = -1.0*K/dx**2d0
                Ah(i,j,4) = -1.0*K/dx**2d0
             endif
             if (j == nhalo+1) then          ! southernmost cell with nonzero x; subdiag term = 0
                Ah(i,j,5) = -1.0*K/dy**2d0            !    Zero-flux Neumann b.c.
                bh(i,j) = 0.0d0 
             elseif (j == ny-nhalo) then     ! northernmost cell with nonzero x; supdiag term = 0
                Ah(i,j,1) = -1.0*K/dy**2d0 
                bh(i,j) = 0.0d0 
             else                            ! interior cell
                Ah(i,j,1) = -1.0*K/dy**2d0
                Ah(i,j,5) = -1.0*K/dy**2d0
             endif
          enddo   ! i
       enddo   ! j

    else   ! tasks > 1
       ! loop over locally owned cells (iglobal = 1 to global_ewn, jglobal = 1 to global_nsn)
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo

             call parallel_globalindex(i, j, iglobal, jglobal)

             active_cell(i,j) = .true.
             Ah(i,j,3) =  4.0d0
             if (iglobal == 1) then                ! westernmost cell with nonzero x; subdiag term = 0
                Ah(i,j,4) = -1.0d0
                bh(i,j) = bh(i,j) + 1.0d0
             elseif (iglobal == global_ewn) then   ! easternmost cell with nonzero x; supdiag term = 0
                Ah(i,j,2) = -1.d0
                bh(i,j) = bh(i,j) + 1.0d0
             else                                  ! interior cell
                Ah(i,j,2) = -1.0d0
                Ah(i,j,4) = -1.0d0
             endif  ! iglobal
             if (jglobal == 1) then                ! southernmost cell with nonzero x; subdiag term = 0
                Ah(i,j,5) = -1.0d0
                bh(i,j) = bh(i,j) + 1.0d0
             elseif (jglobal == global_nsn) then   ! northernmost cell with nonzero x; supdiag term = 0
                Ah(i,j,1) = -1.0d0
                bh(i,j) = bh(i,j) + 1.0d0
             else                                  ! interior cell
                Ah(i,j,1) = -1.0d0
                Ah(i,j,5) = -1.0d0
             endif   ! jglobal
          enddo   ! i
       enddo   ! j

    endif   ! tasks

    do m = 1, 5
       call parallel_halo(Ah(:,:,m))
    enddo
    call parallel_halo(bh)
    call parallel_halo(active_cell)

    ! print out some matrix elements
    if (this_rank == rtest) then
       j = jtest
       print*, ' '
       print*, 'Matrix elements, itest, jtest =', itest, jtest
       do i = nhalo+1, nx-nhalo
          write(6,'(i4,l4,6f10.3)') i, active_cell(i,j), Ah(i,j,1:5), bh(i,j)
       enddo
    endif

    call pcg_solver_standard_2d_scalar(&
         nx,        ny,            &
         nhalo,                    &
         indxA,     size(Ah,3),    &
         active_cell,              &
         Ah,        bh,            &
         xh,                       &
         model%options%which_ho_precond,     &
         model%options%linear_solve_ncheck,  &
         err,       niters,        &
         itest, jtest, rtest)

    ! print out the solution on the local processor
    if (this_rank == rtest) then
       print*, ' '
       print*, 'Solution:'
       do j = ny, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx
             write(6,'(f5.1)',advance='no') xh(i,j)
          enddo
       print*, ' '
       enddo
       print*, 'err, niters:', err, niters
    endif

    ! clean up
    deallocate(active_cell)
    deallocate(Ah, bh, xh)

  end subroutine glissade_shakti

end module glissade_basal_water
