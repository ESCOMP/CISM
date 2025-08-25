!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_test.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module holds some test subroutines for the Glissade dynamical core
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glissade_test

  use glimmer_global, only: dp
  use glimmer_log
  use glimmer_utils, only: point_diag
  use glide_types

  implicit none

  private
  public :: glissade_test_halo, glissade_test_transport

contains

!=======================================================================

  subroutine glissade_test_halo(model)

    use cism_parallel, only: main_task, this_rank, uhalo, lhalo, staggered_lhalo, staggered_uhalo, &
         parallel_type, parallel_halo, staggered_parallel_halo, parallel_globalID_scalar, &
         parallel_globalindex, parallel_halo_tracers

    ! various tests of parallel halo updates

    ! print statements are formatted for a 30x30 global array of scalars
    ! (34x34 with nhalo = 2), as for the standard dome problem

    type(glide_global_type), intent(inout) :: model      ! model instance

    type(parallel_type) :: parallel   ! info for parallel communication

    integer, dimension (:,:), allocatable    ::  pgID       ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDr4     ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDr8     ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDr8_3d  ! unique global ID for parallel runs  

    integer, dimension (:,:), allocatable    ::  pgIDstagi  ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDstagr  ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDstagr3 ! unique global ID for parallel runs  

    logical,  dimension(:,:), allocatable   :: logvar
    integer,  dimension(:,:), allocatable   :: intvar
    real,     dimension(:,:), allocatable   :: r4var
    real(dp), dimension(:,:), allocatable   :: r8var
    real(dp), dimension(:,:,:), allocatable :: r8var_3d

    integer,  dimension(:,:), allocatable   :: intvar_2d_stag
    integer,  dimension(:,:,:), allocatable   :: intvar_3d_stag
    real(dp), dimension(:,:), allocatable :: r8var_2d_stag
    real(dp), dimension(:,:,:), allocatable :: r8var_3d_stag
    real(dp), dimension(:,:,:,:), allocatable :: r8var_4d_stag

    real(dp), dimension (:,:,:,:), allocatable ::  tracers

    integer :: i, j, k, ig, jg
    integer :: nx, ny, nz

    integer :: itest, jtest, rtest

    write(6,*) ' '
    write(6,*) 'In glissade_test_halo, this_rank =', this_rank

    nx = model%general%ewn
    ny = model%general%nsn
    nz = model%general%upn
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local
    rtest = model%numerics%rdiag_local
    parallel = model%parallel

    allocate(logvar(nx,ny))
    allocate(intvar(nx,ny))
    allocate(r4var(nx,ny))
    allocate(r8var(nx,ny))
    allocate(r8var_3d(nz,nx,ny))

    allocate(intvar_2d_stag(nx-1,ny-1))
    allocate(intvar_3d_stag(nz,nx-1,ny-1))
    allocate(r8var_2d_stag(nx-1,ny-1))
    allocate(r8var_3d_stag(nz,nx-1,ny-1))
    allocate(r8var_4d_stag(2,nz,nx-1,ny-1))

    allocate(pgID(nx,ny))
    allocate(pgIDr4(nx,ny))
    allocate(pgIDr8(nx,ny))
    allocate(pgIDr8_3d(nz,nx,ny))
    allocate(pgIDstagi(nx-1,ny-1))
    allocate(pgIDstagr(nx-1,ny-1))
    allocate(pgIDstagr3(nz,nx-1,ny-1))

    if (main_task) then
       write(6,*) ' '
       write(6,*) 'nx, ny, nz =', nx, ny, nz
       write(6,*) 'uhalo, lhalo =', uhalo, lhalo
       write(6,*) 'global_ewn, global_nsn =', parallel%global_ewn, parallel%global_nsn
       write(6,*) ' '
    endif

    write(6,*) 'this_rank, global_row/col offset =', &
         this_rank, parallel%global_row_offset, parallel%global_col_offset

    ! Test some standard parallel_halo routines for scalars: logical_2d, integer_2d, real4_2d, real8_2d, real8_3d

    ! logical 2D field

    logvar(:,:) = .false.

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       logvar(i,j) = .true.
    enddo
    enddo

    call point_diag(logvar, 'Before halo update, logvar', itest, jtest, rtest, 7, 7)
    call parallel_halo(logvar, parallel)
    call point_diag(logvar, 'After halo update, logvar', itest, jtest, rtest, 7, 7)

    ! staggered 2D int field

    intvar_2d_stag(:,:) = 0.0d0

    do j = staggered_lhalo, ny-staggered_uhalo-1
       do i = staggered_lhalo, nx-staggered_uhalo-1
          call parallel_globalindex(i, j, ig, jg, parallel)
          intvar_2d_stag(i,j) = 100*ig + jg
       enddo
    enddo

    call point_diag(intvar_2d_stag, 'Before halo update, intvar_2d_stag', itest, jtest, rtest, 7, 7)
    call staggered_parallel_halo(intvar_2d_stag, parallel)
    call point_diag(intvar_2d_stag, 'After halo update, intvar_2d_stag', itest, jtest, rtest, 7, 7)

    ! staggered 3D integer field

    intvar_3d_stag(:,:,:) = 0.0d0

    do j = staggered_lhalo, ny-staggered_uhalo-1
       do i = staggered_lhalo, nx-staggered_uhalo-1
          call parallel_globalindex(i, j, ig, jg, parallel)
          intvar_3d_stag(:,i,j) = 100*ig + jg
       enddo
    enddo

    call point_diag(intvar_3d_stag(1,:,:), 'Before halo update, intvar_3d_stag, k = 1', &
         itest, jtest, rtest, 7, 7)
    call staggered_parallel_halo(intvar_3d_stag, parallel)
    call point_diag(intvar_3d_stag(1,:,:), 'After halo update, intvar_3d_stag, k = 1', &
         itest, jtest, rtest, 7, 7)

    ! staggered 2D real field

    r8var_2d_stag(:,:) = 0.0d0

    do j = staggered_lhalo, ny-staggered_uhalo-1
       do i = staggered_lhalo, nx-staggered_uhalo-1
          call parallel_globalindex(i, j, ig, jg, parallel)
          r8var_2d_stag(i,j) = 100.d0*real(ig,dp) + real(jg,dp)
       enddo
    enddo

    call point_diag(r8var_2d_stag, 'Before halo update, r8var_2d_stag', itest, jtest, rtest, 7, 7)
    call staggered_parallel_halo(r8var_2d_stag, parallel)
    call point_diag(r8var_2d_stag, 'After halo update, r8var_2d_stag', itest, jtest, rtest, 7, 7)

    ! staggered 3D real field

    r8var_3d_stag(:,:,:) = 0.0d0

    do j = staggered_lhalo, ny-staggered_uhalo-1
       do i = staggered_lhalo, nx-staggered_uhalo-1
          call parallel_globalindex(i, j, ig, jg, parallel)
          r8var_3d_stag(:,i,j) = 100.d0*real(ig,dp) + real(jg,dp)
       enddo
    enddo

    call point_diag(r8var_3d_stag(1,:,:), 'Before halo update, r8var_3d_stag, k = 1', &
         itest, jtest, rtest, 7, 7)
    call staggered_parallel_halo(r8var_3d_stag, parallel)
    call point_diag(r8var_3d_stag(1,:,:), 'After halo update, r8var_3d_stag, k = 1', &
         itest, jtest, rtest, 7, 7)

    ! staggered 4D real field

    r8var_4d_stag(:,:,:,:) = 0.0d0

    do j = staggered_lhalo, ny-staggered_uhalo-1
       do i = staggered_lhalo, nx-staggered_uhalo-1
          call parallel_globalindex(i, j, ig, jg, parallel)
          r8var_4d_stag(:,:,i,j) = 100.d0*real(ig,dp) + real(jg,dp)
       enddo
    enddo

    call point_diag(r8var_4d_stag(1,1,:,:), 'Before halo update, r4var_3d_stag, k = 1, l = 1', &
         itest, jtest, rtest, 7, 7)
    call staggered_parallel_halo(r8var_4d_stag, parallel)
    call point_diag(r8var_4d_stag(1,1,:,:), 'After halo update, r4var_3d_stag, k = 1, l = 1', &
         itest, jtest, rtest, 7, 7)

    ! The next part of the code concerns parallel global IDs

    ! Compute parallel global ID for each grid cell

    pgID(:,:) = 0   ! integer

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgID(i,j) = parallel_globalID_scalar(i,j,nz,parallel)    ! function in parallel_mpi.F90
    enddo
    enddo

    call point_diag(pgID, 'Before halo update, pgID', itest, jtest, rtest, 7, 7)
    call parallel_halo(pgID, parallel)
    call point_diag(pgID, 'After halo update, pgID', itest, jtest, rtest, 7, 7)

    ! real 2D
    
    pgIDr4(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDr4(i,j) = real(parallel_globalID_scalar(i,j,nz,parallel))
    enddo
    enddo

    call point_diag(pgIDr4, 'Before halo update, pgIDr4', itest, jtest, rtest, 7, 7)
    call parallel_halo(pgIDr4, parallel)
    call point_diag(pgIDr4, 'After halo update, pgIDr4', itest, jtest, rtest, 7, 7)

    ! double 2D
    
    pgIDr8(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDr8(i,j) = real(parallel_globalID_scalar(i,j,nz,parallel), dp)
    enddo
    enddo

    call point_diag(pgIDr8, 'Before halo update, pgIDr8', itest, jtest, rtest, 7, 7)
    call parallel_halo(pgIDr8, parallel)
    call point_diag(pgIDr8, 'After halo update, pgIDr8', itest, jtest, rtest, 7, 7)

    ! double 3D

    pgIDr8_3d(:,:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          ! function in parallel_mpi.F90
          pgIDr8_3d(k,i,j) = real(parallel_globalID_scalar(i,j,nz,parallel),dp) + real(k,dp)
       enddo
    enddo
    enddo

    call point_diag(pgIDr8_3d(1,:,:), 'Before halo update, pgIDr8_3d, k = 1', itest, jtest, rtest, 7, 7)
    call parallel_halo(pgIDr8_3d, parallel)
    call point_diag(pgIDr8_3d(1,:,:), 'After halo update, pgIDr8_3d, k = 1', itest, jtest, rtest, 7, 7)

    ! Repeat for staggered variables

    ! First for an integer 2D field

    pgIDstagi(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDstagi(i,j) = parallel_globalID_scalar(i,j,nz,parallel)    ! function in parallel_mpi.F90
    enddo
    enddo

    call point_diag(pgIDstagi, 'Before halo update, pgIDstagi', itest, jtest, rtest, 7, 7)
    call staggered_parallel_halo(pgIDstagi, parallel)
    call point_diag(pgIDstagi, 'After halo update, pgIDstagi', itest, jtest, rtest, 7, 7)

    ! Then for a real 2D field

    pgIDstagr(:,:) = 0.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDstagr(i,j) = real(parallel_globalID_scalar(i,j,nz,parallel),dp)    ! function in parallel_mpi.F90
    enddo
    enddo

    call point_diag(pgIDstagr, 'Before halo update, pgIDstagr', itest, jtest, rtest, 7, 7)
    call staggered_parallel_halo(pgIDstagr, parallel)
    call point_diag(pgIDstagr, 'After halo update, pgIDstagr', itest, jtest, rtest, 7, 7)

    ! Then for a real 3D field

    pgIDstagr3(:,:,:) = 0.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          ! function in parallel_mpi.F90
          pgIDstagr3(k,i,j) = real(parallel_globalID_scalar(i,j,nz,parallel),dp) + real(k,dp)
       enddo
    enddo
    enddo

    call point_diag(pgIDstagr3(1,:,:), 'Before halo update, pgIDstagr3, k = 1', itest, jtest, rtest, 7, 7)
    call staggered_parallel_halo(pgIDstagr3, parallel)
    call point_diag(pgIDstagr3(1,:,:), 'After halo update, pgIDstagr3, k = 1', itest, jtest, rtest, 7, 7)

    deallocate(logvar)
    deallocate(intvar)
    deallocate(r4var)
    deallocate(r8var)
    deallocate(r8var_3d)
    deallocate(pgID)
    deallocate(pgIDr4)
    deallocate(pgIDr8)
    deallocate(pgIDr8_3d)
    deallocate(pgIDstagi)
    deallocate(pgIDstagr)
    deallocate(pgIDstagr3)

    ! Test tracer update routine

    allocate(tracers(nx,ny,2,2))  ! 2 tracers, 2 layers
    tracers(:,:,:,:) = 0.d0

    do k = 1, 2
       do j = 1+lhalo, ny-uhalo
       do i = 1+lhalo, nx-uhalo
          ! function in parallel_mpi.F90
          tracers(i,j,:,k) = real(parallel_globalID_scalar(i,j,k,parallel),dp) + real(k,dp)
       enddo
       enddo
    enddo

    call point_diag(tracers(:,:,1,1), 'Before halo update, tracers, k = 1, l = 1)', itest, jtest, rtest, 7, 7)
    call parallel_halo_tracers(tracers, parallel)
    call point_diag(tracers(:,:,1,1), 'After halo update, tracers, k = 1, l = 1)', itest, jtest, rtest, 7, 7)

    deallocate(tracers)

    stop

  end subroutine glissade_test_halo

!=======================================================================

  subroutine glissade_test_transport(model)

    use glissade_transport, only: glissade_transport_driver, &
         glissade_transport_setup_tracers, glissade_transport_finish_tracers
    use glimmer_physcon, only: pi, scyr

    !-------------------------------------------------------------------
    ! Test transport of a cylinder or block of ice once around the domain and
    ! back to the starting point, assuming uniform motion in a straight line.
    !
    ! Instructions for use:
    ! (1) At the top of this module, set test_transport = .true. and choose a
    !     value for thk_init.
    ! (2) Set the config file to run with the glissade dycore for one timestep, 
    !     with CF output written every timestep. 
    !     Note: Whatever the initial ice geometry in the config file
    !     (e.g., the dome test case), the ice extent will be preserved, but
    !     the thickness will be set everywhere to thk_init, giving a steep
    !     front at the ice margin that is challenging for transport schemes.
    ! (3) Comment out the call to glissade_diagnostic_variable_solve in
    !     cism_driver or simple_glide.  (It probable won't hurt to
    !     have it turned on, but will just use extra cycles.)
    !
    ! During the run, the following will happen:
    ! (1) During glissade_initialise, the ice thickness will be set to
    !     thk_init (everywhere there is ice).
    ! (2) During the first call to glissade_tstep, this subroutine 
    !     (glissade_test_transport) will be called.  The ice will be transported 
    !     around to its initial starting point (assuming periodic BCs).  
    ! (3) Then the model will return from glissade_tstep before doing any 
    !     other computations. Output will be written to netCDF and the code 
    !     will complete.
    !
    ! There should be two time slices in the output netCDF file.
    ! Compare the ice geometry at these two slices to see how diffusive 
    !  and/or dispersive the transport scheme is. A perfect scheme
    !  would leave the geometry unchanged.  
    !-------------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model      ! model instance

    type(parallel_type) :: parallel   ! info for parallel communication

    integer :: ntracers   ! number of tracers to be transported

    real(dp), dimension(:,:,:), allocatable :: uvel, vvel   ! uniform velocity field (m/yr)

    integer :: i, j, k, n
    integer :: nx, ny, nz
    real(dp) :: dx, dy

    real(dp), parameter :: umag = 100.      ! uniform speed (m/yr)

    ! Set angle of motion
    !WHL - Tested all of these angles (eastward, northward, and northeastward)
    real(dp), parameter :: theta = 0.d0      ! eastward
!    real(dp), parameter :: theta = pi/2.d0   ! northward
!    real(dp), parameter :: theta = pi/4.d0   ! northeastward

    real(dp), parameter :: thk = 500.d0

    real(dp) :: dt          ! time step in yr

    integer :: ntstep       ! run for this number of timesteps
    real(dp) :: theta_c    ! angle dividing paths through EW walls from paths thru NS walls
    real(dp) :: lenx       ! length of shortest path thru EW walls
    real(dp) :: leny       ! length of shortest path thru NS walls
    real(dp) :: len_path   ! length of path back to starting point
    real(dp) :: adv_cfl    ! advective CFL number

    integer :: itest, jtest, rtest

    logical :: do_upwind_transport  ! if true, do upwind transport

    ! Initialize

    dx = model%numerics%dew
    dy = model%numerics%dns
    nx = model%general%ewn
    ny = model%general%nsn
    nz = model%general%upn
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local
    rtest = model%numerics%rdiag_local
    parallel = model%parallel

    allocate(uvel(nz,nx-1,ny-1))
    allocate(vvel(nz,nx-1,ny-1))
    ! Find the length of the path around the domain and back to the starting point

    lenx = parallel%global_ewn * dx
    leny = parallel%global_nsn * dy
    theta_c = atan(leny/lenx)   ! 0 <= theta_c <= pi/2

    if ( (theta >= -theta_c   .and. theta <= theta_c) .or.   &
         (theta >= pi-theta_c .and. theta <= pi+theta_c) ) then
       ! path will cross east/west wall
       len_path = lenx / abs(cos(theta))
    else
       ! path will cross north/south wall
       len_path = leny / abs(sin(theta))
    endif

    ! Choose the number of time steps such that the ice will travel
    ! less than one grid cell per time step

    ntstep = nint(len_path/dx) + 10   ! 10 is arbitrary

    ! Choose the time step such that the ice will go around the domain
    ! exactly once in the chosen number of time steps.

    dt = len_path / (umag*ntstep)

    ! CFL check, just to be sure

    adv_cfl = max (dt*umag*cos(theta)/dx, dt*umag*sin(theta)/dy)
    
    if (adv_cfl >= 1.d0) then
       write(6,*) 'dt is too big for advective CFL; increase ntstep to', ntstep * adv_cfl
       stop
    endif

    ! Print some diagnostics

    if (main_task) then
       write(6,*) ' '
       write(6,*) 'In glissade_test_transport'
       write(6,*) 'nx, ny, nz =', nx, ny, nz
       write(6,*) 'len_path =', len_path
       write(6,*) 'umag (m/yr) =', umag
       write(6,*) 'dt (yr) =', dt
       write(6,*) 'ntstep =', ntstep
       write(6,*) 'theta (deg) =', theta * 180.d0/pi
    endif

    call point_diag(model%geometry%thck, 'Initial thck', itest, jtest, rtest, 7, 7)

    ! Set uniform ice speed everywhere

    do j = 1, ny-1
    do i = 1, nx-1
       do k = 1, nz
          uvel(k,i,j) = umag * cos(theta)
          vvel(k,i,j) = umag * sin(theta)
       enddo
    enddo
    enddo

    ! Determine which transport scheme

    if (model%options%whichevol == EVOL_UPWIND) then
       do_upwind_transport = .true.
    else
       do_upwind_transport = .false.
    endif

    ! Transport the ice around the domain

    do n = 1, ntstep

       ! call transport scheme
       ! Note: glissade_transport_driver expects dt in seconds, uvel/vvel in m/s

       call glissade_transport_setup_tracers(model)

       call glissade_transport_driver(dt*scyr,                                              &
                                      dx,                        dy,                        &
                                      nx,                        ny,                        &
                                      nz-1,                      model%numerics%sigma,      &
                                      parallel,                                             &
                                      model%numerics%idiag_local,                           &
                                      model%numerics%jdiag_local,                           &
                                      model%numerics%rdiag_local,                           &
                                      uvel(:,:,:)/scyr,          vvel(:,:,:)/scyr,          &
                                      model%geometry%thck(:,:),                             &
                                      model%geometry%ntracers,                              &
                                      model%geometry%tracers(:,:,:,:),                      &
                                      model%geometry%tracers_usrf(:,:,:),                   &
                                      model%geometry%tracers_lsrf(:,:,:),                   &
                                      model%options%which_ho_vertical_remap,                &
                                      upwind_transport_in = do_upwind_transport)

       call glissade_transport_finish_tracers(model)

       call point_diag(model%geometry%thck, 'New thck', itest, jtest, rtest, 7, 7)

    enddo  ! ntstep

    if (main_task) write(6,*) 'Done in glissade_test_transport'

    deallocate(uvel)
    deallocate(vvel)

  end subroutine glissade_test_transport

!=======================================================================

  end module glissade_test

!=======================================================================
