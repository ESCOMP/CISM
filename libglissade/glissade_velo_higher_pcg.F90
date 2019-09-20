!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_higher_pcg.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains subroutines called from glissade_velo_higher.F90 and used 
! to solve the problem Ax = b using the preconditioned conjugate gradient method.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

  module glissade_velo_higher_pcg

    use glimmer_global, only: dp
    use glide_types   ! for preconditioning options
    use glimmer_log
    use profile
    use parallel
    
    implicit none
    
    private
    public :: pcg_solver_standard_3d,  pcg_solver_standard_2d,  &
              pcg_solver_chrongear_3d, pcg_solver_chrongear_2d, &
              matvec_multiply_structured_3d
    
    logical :: verbose_pcg

    interface global_sum_staggered
       module procedure global_sum_staggered_3d_real8
       module procedure global_sum_staggered_3d_real8_nvar       
       module procedure global_sum_staggered_2d_real8
       module procedure global_sum_staggered_2d_real8_nvar       
    end interface

    ! linear solver settings
    !TODO - Pass in these solver settings as arguments?
    integer, parameter ::    &
       maxiters = 200        ! max number of linear iterations before quitting
!!       maxiters = 1000        ! max number of linear iterations before quitting

    real(dp), parameter ::   &
!!       tolerance = 1.d-11    ! tolerance for linear solver (old value; more stringent than necessary)
       tolerance = 1.d-08    ! tolerance for linear solver

    integer, parameter :: &
       solve_ncheck = 5      ! check for convergence every solve_ncheck iterations

    logical, parameter :: verbose_tridiag = .false.
!!    logical, parameter :: verbose_tridiag = .true.

    !TODO - Add verbose_pcg here.

  contains

!****************************************************************************

  subroutine pcg_solver_standard_3d(nx,        ny,            &
                                    nz,        nhalo,         &
                                    indxA,     active_vertex, &
                                    Auu,       Auv,           &
                                    Avu,       Avv,           &
                                    bu,        bv,            &
                                    xu,        xv,            &
                                    precond,   err,           &
                                    niters,                   &
                                    itest_in,  jtest_in,  rtest_in,   &
                                    verbose) 

    !---------------------------------------------------------------
    !  This subroutine uses a standard preconditioned conjugate-gradient algorithm
    !  to solve the equation $Ax=b$.
    !  Convergence is checked every {\em solve_ncheck} steps.
    !
    !  It is based on the barotropic solver in the POP ocean model 
    !  (author Phil Jones, LANL).  Input and output arrays are located
    !  on a structured (i,j,k) grid as defined in the glissade_velo_higher
    !  module.  The global matrix is sparse, but its nonzero elements
    !  are stored in four dense matrices called Auu, Avv, Auv, and Avu.
    !  Each matrix has 3x3x3 = 27 potential nonzero elements per
    !  node (i,j,k).
    !
    !  The current preconditioning options are
    !  (0) no preconditioning
    !  (1) diagonal preconditioning
    !  (2) preconditioning using a physics-based SIA solver
    ! 
    !  For the dome test case with higher-order dynamics, option (2) is best. 
    !
    !  Here is a schematic of the method implemented below for solving Ax = b:
    !
    !  halo_update(x0)
    !  r0 = b - A*x0
    !  d0 = 0
    !  eta0 = 1
    !
    !  while (not converged)
    !     solve Mz = r for z
    !     eta1 = (r,z)
    !     beta = eta1/eta0
    !     d = z + beta*d
    !     halo_update(d)
    !     eta0 = eta1
    !     q = Ad
    !     eta2 = (d,q)
    !     alpha = eta1/eta2
    !     x = x + alpha*d
    !     r = r - alpha*q (or occasionally, r = b - Ax)
    !     Check for convergence: err = sqrt(r,r)/sqrt(b,b) < tolerance
    !  end while
    !
    !  where x = solution (initial value = x0)
    !        d = conjugate direction vector (initial value = d0)
    !        r = residual vector (initial value = r0)
    !        M = preconditioning matrix
    !    (r,z) = dot product of vectors r and z
    !            and similarly for (d,q)
    !       
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
                                ! velocity grid has dimensions (nx-1,ny-1)
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers (for scalars)

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F
 
    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       xu, xv             ! u and v components of solution (i.e., uvel and vvel)

    integer, intent(in)  ::   &
       precond           ! = 0 for no preconditioning
                         ! = 1 for diagonal preconditioning (best option for SSA-dominated flow)
                         ! = 2 for preconditioning with SIA solver (works well for SIA-dominated flow)

    real(dp), intent(out) ::  &
       err                               ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                            ! iterations needed to solution

    integer, intent(in), optional :: &
       itest_in, jtest_in, rtest_in      ! point for debugging diagnostics

    logical, intent(in), optional :: &
       verbose                           ! if true, print diagnostic output
    
    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer ::  i, j, k      ! grid indices
    integer ::  iA, jA, kA   ! grid offsets ranging from -1 to 1
    integer ::  m            ! matrix element index
    integer ::  n            ! iteration counter

    real(dp) ::           &
       eta0, eta1, eta2,  &! scalar inner product results
       alpha,             &! eta1/eta2 = term in expression for new residual and solution
       beta                ! eta1/eta0 = term in expression for new direction vector

    ! vectors (each of these is split into u and v components)
    real(dp), dimension(nz,nx-1,ny-1) ::  &
       Adiagu, Adiagv,    &! diagonal terms of matrices Auu and Avv
       ru, rv,            &! residual vector (b-Ax)
       du, dv,            &! conjugate direction vector
       qu, qv,            &! A*d
       zu, zv,            &! solution of Mz = r (also used as a temporary vector)
       work0u, work0v      ! cg intermediate results

    real(dp) ::  &
       L2_resid,          &! L2 norm of residual vector Ax-b
       L2_rhs              ! L2 norm of rhs vector b
                           ! solver converges when L2_resid/L2_rhs < tolerance

    real(dp), dimension(-1:1,nz,nx-1,ny-1) ::  &
       Muu, Mvv            ! simplified SIA matrices for preconditioning

    integer :: itest, jtest, rtest

    if (present(itest_in)) then
       itest = itest_in
    else
       itest = nx/2
    endif

    if (present(itest_in)) then
       jtest = jtest_in
    else
       jtest = ny/2
    endif

    if (present(itest_in)) then
       rtest = rtest_in
    else
       rtest = 0
    endif

    if (present(verbose)) then
       verbose_pcg = verbose
    else
       verbose_pcg = .false.   ! for debugging
    endif

    if (verbose_pcg .and. main_task) then
       print*, 'Using native PCG solver (standard)'
       print*, 'tolerance, maxiters, precond =', tolerance, maxiters, precond
    endif

    ! Set up matrices for preconditioning

    call t_startf("pcg_precond_init")
    call setup_preconditioner_3d(nx,         ny,       &
                                 nz,                   &
                                 precond,    indxA,    &
                                 Auu,        Avv,      &
                                 Adiagu,     Adiagv,   &
                                 Muu,        Mvv)
    call t_stopf("pcg_precond_init")

    ! Compute initial residual and initialize the direction vector d
    ! Note: The matrix A must be complete for all rows corresponding to locally 
    !        owned vertices, and x must have the correct values in
    !        halo vertices bordering the locally owned vertices.
    !       Then y = Ax will be correct for locally owned vertices.

    ! Halo update for x (initial guess for velocity solution)

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(xu)
    call staggered_parallel_halo(xv)
    call t_stopf("pcg_halo_init")

    ! Compute A*x (use z as a temp vector for A*x)

    call t_startf("pcg_matmult_init")
    call matvec_multiply_structured_3d(nx,        ny,            &
                                       nz,        nhalo,         &
                                       indxA,     active_vertex, &
                                       Auu,       Auv,           &
                                       Avu,       Avv,           &
                                       xu,        xv,            &
                                       zu,        zv)
    call t_stopf("pcg_matmult_init")

    ! Compute the initial residual r(0) = b - Ax(0)
    ! This will be correct for locally owned vertices.

    call t_startf("pcg_vecupdate_init")
    ru(:,:,:) = bu(:,:,:) - zu(:,:,:)
    rv(:,:,:) = bv(:,:,:) - zv(:,:,:)
    call t_stopf("pcg_vecupdate_init")

    ! Initialize scalars and vectors

    niters = maxiters 
    eta0 = 1.d0

    du(:,:,:) = 0.d0
    dv(:,:,:) = 0.d0

    zu(:,:,:) = 0.d0
    zv(:,:,:) = 0.d0

    ! Compute the L2 norm of the RHS vectors
    ! (Goal is to obtain L2_resid/L2_rhs < tolerance)

    call t_startf("pcg_dotprod")
    work0u(:,:,:) = bu(:,:,:)*bu(:,:,:)    ! terms of dot product (b, b)
    work0v(:,:,:) = bv(:,:,:)*bv(:,:,:)
    call t_stopf("pcg_dotprod")

    ! find global sum of the squared L2 norm

    call t_startf("pcg_glbsum_init")
    call global_sum_staggered(nx,     ny,     &
                              nz,     nhalo,  &
                              L2_rhs,         &
                              work0u, work0v)
    call t_stopf("pcg_glbsum_init")

    ! take square root

    L2_rhs = sqrt(L2_rhs)       ! L2 norm of RHS

    ! iterate to solution

    iter_loop: do n = 1, maxiters

       call t_startf("pcg_precond")

       ! Compute PC(r) = solution z of Mz = r

       if (precond == 0) then      ! no preconditioning

           zu(:,:,:) = ru(:,:,:)         ! PC(r) = r     
           zv(:,:,:) = rv(:,:,:)         ! PC(r) = r    

       elseif (precond == 1 ) then  ! diagonal preconditioning

          do j = 1, ny-1
          do i = 1, nx-1
          do k = 1, nz
             if (Adiagu(k,i,j) /= 0.d0) then
                zu(k,i,j) = ru(k,i,j) / Adiagu(k,i,j)   ! PC(r), where PC is formed from diagonal elements of A
             else                                        
                zu(k,i,j) = 0.d0
             endif
             if (Adiagv(k,i,j) /= 0.d0) then
                zv(k,i,j) = rv(k,i,j) / Adiagv(k,i,j)  
             else                                        
                zv(k,i,j) = 0.d0
             endif
          enddo    ! k
          enddo    ! i
          enddo    ! j

       elseif (precond == 2) then   ! local vertical shallow-ice solver for preconditioning

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Muu,  ru,   zu)      ! solve Muu*zu = ru for zu 

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Mvv,  rv,   zv)      ! solve Mvv*zv = rv for zv

       endif    ! precond

       call t_stopf("pcg_precond")

       ! Compute the dot product eta1 = (r, PC(r))

       call t_startf("pcg_dotprod")
       work0u(:,:,:) = ru(:,:,:)*zu(:,:,:)    ! terms of dot product (r, PC(r))
       work0v(:,:,:) = rv(:,:,:)*zv(:,:,:)    
       call t_stopf("pcg_dotprod")

       call t_startf("pcg_glbsum_iter")
       call global_sum_staggered(nx,     ny,     &
                                 nz,     nhalo,  &
                                 eta1,           &
                                 work0u, work0v)
       call t_stopf("pcg_glbsum_iter")

       !WHL - If the SIA solver has failed due to singular matrices,
       !      then eta1 will be NaN.
 
       if (eta1 /= eta1) then  ! eta1 is NaN
          call write_log('PCG solver has failed, alpha = NaN', GM_FATAL)
       endif

       ! Update the conjugate direction vector d

       beta = eta1/eta0

       call t_startf("pcg_vecupdate")
       du(:,:,:) = zu(:,:,:) + beta*du(:,:,:)       ! d_(i+1) = PC(r_(i+1)) + beta_(i+1)*d_i
       dv(:,:,:) = zv(:,:,:) + beta*dv(:,:,:)       !
                                                    !                    (r_(i+1), PC(r_(i+1)))
                                                    ! where beta_(i+1) = --------------------  
                                                    !                        (r_i, PC(r_i)) 
                                                    ! Initially eta0 = 1  
                                                    ! For n >=2, eta0 = old eta1
       call t_stopf("pcg_vecupdate")

       ! Halo update for d

       call t_startf("pcg_halo_iter")
       call staggered_parallel_halo(du)
       call staggered_parallel_halo(dv)
       call t_stopf("pcg_halo_iter")
  
       ! Compute q = A*d
       ! This is the one matvec multiply required for each iteration

       call t_startf("pcg_matmult_iter")
       call matvec_multiply_structured_3d(nx,        ny,            &
                                          nz,        nhalo,         &
                                          indxA,     active_vertex, &
                                          Auu,       Auv,           &
                                          Avu,       Avv,           &
                                          du,        dv,            &
                                          qu,        qv)
       call t_stopf("pcg_matmult_iter")

       ! Copy old eta1 = (r, PC(r)) to eta0

       eta0 = eta1               ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i)) 

       ! Compute the dot product eta2 = (d, A*d)

       call t_startf("pcg_dotprod")
       work0u(:,:,:) = du(:,:,:) * qu(:,:,:)       ! terms of dot product (d, Ad)
       work0v(:,:,:) = dv(:,:,:) * qv(:,:,:)
       call t_stopf("pcg_dotprod")

       call t_startf("pcg_glbsum_iter")
       call global_sum_staggered(nx,     ny,     &
                                 nz,     nhalo,  &
                                 eta2,           &
                                 work0u, work0v)
       call t_stopf("pcg_glbsum_iter")

       ! Compute alpha
                              !          (r, PC(r))
       alpha = eta1/eta2      ! alpha = ----------
                              !          (d, A*d)
       
       !WHL - If eta2 = 0 (e.g., because all matrix entries are zero), then alpha = NaN
 
       if (alpha /= alpha) then  ! alpha is NaN
!!          write(6,*) 'eta1, eta2, alpha:', eta1, eta2, alpha
          call write_log('PCG solver has failed, alpha = NaN', GM_FATAL)
       endif

       ! Compute the new solution and residual

       call t_startf("pcg_vecupdate")
       xu(:,:,:) = xu(:,:,:) + alpha * du(:,:,:)    ! new solution, x_(i+1) = x_i + alpha*d
       xv(:,:,:) = xv(:,:,:) + alpha * dv(:,:,:)

       ru(:,:,:) = ru(:,:,:) - alpha * qu(:,:,:)    ! new residual, r_(i+1) = r_i - alpha*(Ad)
       rv(:,:,:) = rv(:,:,:) - alpha * qv(:,:,:)
       call t_stopf("pcg_vecupdate")

       ! Check for convergence every solve_ncheck iterations
       ! For convergence check, use r = b - Ax

       if (mod(n,solve_ncheck) == 0) then

          ! Halo update for x

          call t_startf("pcg_halo_resid")
          call staggered_parallel_halo(xu)
          call staggered_parallel_halo(xv)
          call t_stopf("pcg_halo_resid")

          ! Compute A*x (use z as a temp vector for A*x)
           
          call t_startf("pcg_matmult_resid")
          call matvec_multiply_structured_3d(nx,        ny,            &
                                             nz,        nhalo,         &
                                             indxA,     active_vertex, &
                                             Auu,       Auv,           &
                                             Avu,       Avv,           &
                                             xu,        xv,            &
                                             zu,        zv)
          call t_stopf("pcg_matmult_resid")

          ! Compute residual r = b - Ax

          call t_startf("pcg_vecupdate")
          ru(:,:,:) = bu(:,:,:) - zu(:,:,:)
          rv(:,:,:) = bv(:,:,:) - zv(:,:,:)
          call t_stopf("pcg_vecupdate")

          ! Compute squared L2 norm of (r, r)

          call t_startf("pcg_dotprod")
          work0u(:,:,:) = ru(:,:,:)*ru(:,:,:)   ! terms of dot product (r, r)
          work0v(:,:,:) = rv(:,:,:)*rv(:,:,:)
          call t_stopf("pcg_dotprod")

          call t_startf("pcg_glbsum_resid")
          call global_sum_staggered(nx,     ny,       &
                                    nz,     nhalo,    &
                                    L2_resid,         &
                                    work0u, work0v)
          call t_stopf("pcg_glbsum_resid")

          ! take square root
          L2_resid = sqrt(L2_resid)       ! L2 norm of residual

          ! compute normalized error
          err = L2_resid/L2_rhs

          if (verbose_pcg .and. main_task) then
!             print*, ' '
!             print*, 'iter, L2_resid, error =', n, L2_resid, err
          endif

          if (err < tolerance) then
             niters = n
             exit iter_loop
          endif            

       endif    ! solve_ncheck

    enddo iter_loop

!WHL - Without good preconditioning, convergence can be slow, but the solution after maxiters might be good enough.
 
    if (niters == maxiters) then
       if (verbose_pcg .and. main_task) then
          print*, 'Glissade PCG solver not converged'
          print*, 'niters, err, tolerance:', niters, err, tolerance
       endif
    endif

  end subroutine pcg_solver_standard_3d

!****************************************************************************

  subroutine pcg_solver_standard_2d(nx,        ny,            &
                                    nhalo,                    &
                                    indxA,     active_vertex, &
                                    Auu,       Auv,           &
                                    Avu,       Avv,           &
                                    bu,        bv,            &
                                    xu,        xv,            &
                                    precond,   err,           &
                                    niters,                   &
                                    itest_in,  jtest_in,  rtest_in,   &
                                    verbose) 

    !---------------------------------------------------------------
    !  This subroutine uses a standard preconditioned conjugate-gradient algorithm
    !  to solve the equation $Ax=b$.
    !  Convergence is checked every {\em solve_ncheck} steps.
    !
    !  It is similar to subroutine pcg_solver_standard_3d, but modified
    !  to solve for x and y at a single horizontal level, as in the
    !  shallow-shelf approximation.  See the comments in that subroutine
    !  (above) for more details on data structure and solver methods.
    !
    !  Input and output arrays are located on a structured (i,j) grid 
    !  as defined in the glissade_velo_higher module.  The global matrix 
    !  is sparse, but its nonzero element are stored in four dense matrices 
    !  called Auu, Avv, Auv, and Avu. Each matrix has 3x3 = 9 potential 
    !  nonzero elements per node (i,j).
    !
    !  The current preconditioning options are
    !  (0) no preconditioning
    !  (1) diagonal preconditioning
    !
    !  The SIA-based preconditioning optional is not available for a 2D solve.
    ! 
    !  TODO: Add a tridiagonal preconditioning option to this subroutine,
    !        as for subroutine pcg_solver_chrongear_2d.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
                              ! velocity grid has dimensions (nx-1,ny-1)
       nhalo                  ! number of halo layers (for scalars)

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA                  ! maps relative (x,y) coordinates to an index between 1 and 9

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for vertices (i,j) where velocity is computed, else F
 
    real(dp), dimension(9,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 9 (node and its nearest neighbors in x and y direction)
                              ! other dimensions = (x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    real(dp), dimension(nx-1,ny-1), intent(inout) ::   &
       xu, xv             ! u and v components of solution (i.e., uvel and vvel)

    integer, intent(in)  ::   &
       precond           ! = 0 for no preconditioning
                         ! = 1 for diagonal preconditioning (best option for SSA-dominated flow)

    real(dp), intent(out) ::  &
       err                               ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                            ! iterations needed to solution

    integer, intent(in), optional :: &
       itest_in, jtest_in, rtest_in      ! point for debugging diagnostics

    logical, intent(in), optional :: &
       verbose                           ! if true, print diagnostic output
    
    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer ::  i, j         ! grid indices
    integer ::  n            ! iteration counter

    real(dp) ::           &
       eta0, eta1, eta2,  &! scalar inner product results
       alpha,             &! eta1/eta2 = term in expression for new residual and solution
       beta                ! eta1/eta0 = term in expression for new direction vector

    ! vectors (each of these is split into u and v components)
    real(dp), dimension(nx-1,ny-1) ::  &
       Adiagu, Adiagv,    &! diagonal terms of matrices Auu and Avv
       ru, rv,            &! residual vector (b-Ax)
       du, dv,            &! conjugate direction vector
       qu, qv,            &! A*d
       zu, zv,            &! solution of Mz = r (also used as a temporary vector)
       work0u, work0v      ! cg intermediate results

    real(dp) ::  &
       L2_resid,          &! L2 norm of residual vector Ax-b
       L2_rhs              ! L2 norm of rhs vector b
                           ! solver converges when L2_resid/L2_rhs < tolerance

    integer :: itest, jtest, rtest

    if (present(itest_in)) then
       itest = itest_in
    else
       itest = nx/2
    endif

    if (present(itest_in)) then
       jtest = jtest_in
    else
       jtest = ny/2
    endif

    if (present(itest_in)) then
       rtest = rtest_in
    else
       rtest = 0
    endif

    if (present(verbose)) then
       verbose_pcg = verbose
    else
       verbose_pcg = .false.   ! for debugging
    endif

    if (verbose_pcg .and. main_task) then
       print*, 'Using native PCG solver (standard)'
       print*, 'tolerance, maxiters, precond =', tolerance, maxiters, precond
    endif

    ! Set up matrices for preconditioning

    !TODO - Add tridiagonal option

    call t_startf("pcg_precond_init")
    call setup_preconditioner_diag_2d(nx,         ny,       &
                                      indxA,                &
                                      Auu,        Avv,      &
                                      Adiagu,     Adiagv)
    call t_stopf("pcg_precond_init")

    ! Compute initial residual and initialize the direction vector d
    ! Note: The matrix A must be complete for all rows corresponding to locally 
    !        owned vertices, and x must have the correct values in
    !        halo vertices bordering the locally owned vertices.
    !       Then y = Ax will be correct for locally owned vertices.

    ! Halo update for x (initial guess for velocity solution)

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(xu)
    call staggered_parallel_halo(xv)
    call t_stopf("pcg_halo_init")

    ! Compute A*x (use z as a temp vector for A*x)

    call t_startf("pcg_matmult_init")
    call matvec_multiply_structured_2d(nx,        ny,            &
                                       nhalo,                    &
                                       indxA,     active_vertex, &
                                       Auu,       Auv,           &
                                       Avu,       Avv,           &
                                       xu,        xv,            &
                                       zu,        zv)
    call t_stopf("pcg_matmult_init")

    ! Compute the initial residual r(0) = b - Ax(0)
    ! This will be correct for locally owned vertices.

    call t_startf("pcg_vecupdate_init")
    ru(:,:) = bu(:,:) - zu(:,:)
    rv(:,:) = bv(:,:) - zv(:,:)
    call t_stopf("pcg_vecupdate_init")

    ! Initialize scalars and vectors

    niters = maxiters 
    eta0 = 1.d0

    du(:,:) = 0.d0
    dv(:,:) = 0.d0

    zu(:,:) = 0.d0
    zv(:,:) = 0.d0

    ! Compute the L2 norm of the RHS vectors
    ! (Goal is to obtain L2_resid/L2_rhs < tolerance)

    call t_startf("pcg_dotprod")
    work0u(:,:) = bu(:,:)*bu(:,:)    ! terms of dot product (b, b)
    work0v(:,:) = bv(:,:)*bv(:,:)
    call t_stopf("pcg_dotprod")

    ! find global sum of the squared L2 norm

    call t_startf("pcg_glbsum_init")
    call global_sum_staggered(nx,     ny,      &
                              nhalo,  L2_rhs,  &
                              work0u, work0v)
    call t_stopf("pcg_glbsum_init")

    ! take square root

    L2_rhs = sqrt(L2_rhs)       ! L2 norm of RHS

    ! iterate to solution

    iter_loop: do n = 1, maxiters

       call t_startf("pcg_precond")

       ! Compute PC(r) = solution z of Mz = r

       if (precond == 0) then      ! no preconditioning

           zu(:,:) = ru(:,:)         ! PC(r) = r     
           zv(:,:) = rv(:,:)         ! PC(r) = r    

       elseif (precond == 1) then  ! diagonal preconditioning

          do j = 1, ny-1
          do i = 1, nx-1
             if (Adiagu(i,j) /= 0.d0) then
                zu(i,j) = ru(i,j) / Adiagu(i,j)   ! PC(r), where PC is formed from diagonal elements of A
             else                                        
                zu(i,j) = 0.d0
             endif
             if (Adiagv(i,j) /= 0.d0) then
                zv(i,j) = rv(i,j) / Adiagv(i,j)  
             else                                        
                zv(i,j) = 0.d0
             endif
          enddo    ! i
          enddo    ! j

       endif    ! precond

       call t_stopf("pcg_precond")

       ! Compute the dot product eta1 = (r, PC(r))

       call t_startf("pcg_dotprod")
       work0u(:,:) = ru(:,:)*zu(:,:)    ! terms of dot product (r, PC(r))
       work0v(:,:) = rv(:,:)*zv(:,:)    
       call t_stopf("pcg_dotprod")

       call t_startf("pcg_glbsum_iter")
       call global_sum_staggered(nx,     ny,     &
                                 nhalo,  eta1,   &
                                 work0u, work0v)
       call t_stopf("pcg_glbsum_iter")

       !WHL - If the SIA solver has failed due to singular matrices,
       !      then eta1 will be NaN.
 
       if (eta1 /= eta1) then  ! eta1 is NaN
          call write_log('PCG solver has failed, eta1 = NaN', GM_FATAL)
       endif

       ! Update the conjugate direction vector d

       beta = eta1/eta0

       call t_startf("pcg_vecupdate")
       du(:,:) = zu(:,:) + beta*du(:,:)       ! d_(i+1) = PC(r_(i+1)) + beta_(i+1)*d_i
       dv(:,:) = zv(:,:) + beta*dv(:,:)       !
                                              !                    (r_(i+1), PC(r_(i+1)))
                                              ! where beta_(i+1) = --------------------  
                                              !                        (r_i, PC(r_i)) 
                                              ! Initially eta0 = 1  
                                              ! For n >=2, eta0 = old eta1
       call t_stopf("pcg_vecupdate")

       ! Halo update for d

       call t_startf("pcg_halo_iter")
       call staggered_parallel_halo(du)
       call staggered_parallel_halo(dv)
       call t_stopf("pcg_halo_iter")
  
       ! Compute q = A*d
       ! This is the one matvec multiply required for each iteration

       call t_startf("pcg_matmult_iter")
       call matvec_multiply_structured_2d(nx,        ny,            &
                                          nhalo,                    &
                                          indxA,     active_vertex, &
                                          Auu,       Auv,           &
                                          Avu,       Avv,           &
                                          du,        dv,            &
                                          qu,        qv)
       call t_stopf("pcg_matmult_iter")

       ! Copy old eta1 = (r, PC(r)) to eta0

       eta0 = eta1               ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i)) 

       ! Compute the dot product eta2 = (d, A*d)

       call t_startf("pcg_dotprod")
       work0u(:,:) = du(:,:) * qu(:,:)       ! terms of dot product (d, Ad)
       work0v(:,:) = dv(:,:) * qv(:,:)
       call t_stopf("pcg_dotprod")

       call t_startf("pcg_glbsum_iter")
       call global_sum_staggered(nx,     ny,     &
                                 nhalo,  eta2,   &
                                 work0u, work0v)
       call t_stopf("pcg_glbsum_iter")

       ! Compute alpha
                              !          (r, PC(r))
       alpha = eta1/eta2      ! alpha = ----------
                              !          (d, A*d)
       
       !WHL - If eta2 = 0 (e.g., because all matrix entries are zero), then alpha = NaN
 
       if (alpha /= alpha) then  ! alpha is NaN
!!          write(6,*) 'eta1, eta2, alpha:', eta1, eta2, alpha
          call write_log('PCG solver has failed, alpha = NaN', GM_FATAL)
       endif

       ! Compute the new solution and residual

       call t_startf("pcg_vecupdate")
       xu(:,:) = xu(:,:) + alpha * du(:,:)    ! new solution, x_(i+1) = x_i + alpha*d
       xv(:,:) = xv(:,:) + alpha * dv(:,:)

       ru(:,:) = ru(:,:) - alpha * qu(:,:)    ! new residual, r_(i+1) = r_i - alpha*(Ad)
       rv(:,:) = rv(:,:) - alpha * qv(:,:)
       call t_stopf("pcg_vecupdate")

       ! Check for convergence every solve_ncheck iterations
       ! For convergence check, use r = b - Ax

       if (mod(n,solve_ncheck) == 0) then

          ! Halo update for x

          call t_startf("pcg_halo_resid")
          call staggered_parallel_halo(xu)
          call staggered_parallel_halo(xv)
          call t_stopf("pcg_halo_resid")

          ! Compute A*x (use z as a temp vector for A*x)
           
          call t_startf("pcg_matmult_resid")
          call matvec_multiply_structured_2d(nx,        ny,            &
                                             nhalo,                    &
                                             indxA,     active_vertex, &
                                             Auu,       Auv,           &
                                             Avu,       Avv,           &
                                             xu,        xv,            &
                                             zu,        zv)
          call t_stopf("pcg_matmult_resid")

          ! Compute residual r = b - Ax

          call t_startf("pcg_vecupdate")
          ru(:,:) = bu(:,:) - zu(:,:)
          rv(:,:) = bv(:,:) - zv(:,:)
          call t_stopf("pcg_vecupdate")

          ! Compute squared L2 norm of (r, r)

          call t_startf("pcg_dotprod")
          work0u(:,:) = ru(:,:)*ru(:,:)   ! terms of dot product (r, r)
          work0v(:,:) = rv(:,:)*rv(:,:)
          call t_stopf("pcg_dotprod")

          call t_startf("pcg_glbsum_resid")
          call global_sum_staggered(nx,     ny,        &
                                    nhalo,  L2_resid,  &
                                    work0u, work0v)
          call t_stopf("pcg_glbsum_resid")

          ! take square root
          L2_resid = sqrt(L2_resid)       ! L2 norm of residual

          ! compute normalized error
          err = L2_resid/L2_rhs

          if (err < tolerance) then
             niters = n
             exit iter_loop
          endif            

       endif    ! solve_ncheck

    enddo iter_loop

!WHL - Without good preconditioning, convergence can be slow, but the solution after maxiters might be good enough.
 
    if (niters == maxiters) then
       if (verbose_pcg .and. main_task) then
          print*, 'Glissade PCG solver not converged'
          print*, 'niters, err, tolerance:', niters, err, tolerance
       endif
    endif

  end subroutine pcg_solver_standard_2d

!****************************************************************************

  subroutine pcg_solver_chrongear_3d(nx,        ny,            &
                                     nz,        nhalo,         &
                                     indxA,     active_vertex, &
                                     Auu,       Auv,           &
                                     Avu,       Avv,           &
                                     bu,        bv,            &
                                     xu,        xv,            &
                                     precond,   err,           &
                                     niters,                   &
                                     itest_in,  jtest_in,  rtest_in,   &
                                     verbose) 

    !---------------------------------------------------------------
    !  This subroutine uses a Chronopoulos-Gear preconditioned conjugate-gradient
    !  algorithm to solve the equation $Ax=b$.
    !
    !  It is based on the Chronopoulos-Gear PCG solver in the POP ocean model 
    !  (author Frank Bryan, NCAR). It is a rearranged conjugate gradient solver 
    !  that reduces the number of global reductions per iteration from two to one 
    !  (not counting the convergence check).  Convergence is checked every 
    !  {\em solve_ncheck} steps.
    !
    !     References are:
    !
    !     Chronopoulos, A.T., A Class of Parallel Iterative Methods Implemented on Multiprocessors,
    !        Ph.D. thesis, Technical Report UIUCDCS-R-86-1267, Department of Computer Science,
    !        University of Illinois, Urbana, Illinois, pp. 1-116, 1986.
    !
    !     Chronopoulos, A.T., and C.W. Gear. s-step iterative methods
    !        for symmetric linear systems. J. Comput. Appl. Math., 25(2),
    !        153-168, 1989.
    !
    !     Dongarra, J. and V. Eijkhout. LAPACK Working Note 159.
    !        Finite-choice algorithm optimization in conjugate gradients.
    !        Tech. Rep. ut-cs-03-502. Computer Science Department.
    !        University of Tennessee, Knoxville. 2003.
    !
    !     D Azevedo, E.F., V.L. Eijkhout, and C.H. Romine. LAPACK Working
    !        Note 56. Conjugate gradient algorithms with reduced
    !        synchronization overhead on distributed memory multiprocessors.
    !        Tech. Rep. CS-93-185.  Computer Science Department.
    !        University of Tennessee, Knoxville. 1993.
    !---------------------------------------------------------------
    !
    !  The input and output arrays are located on a structured (i,j,k) grid 
    !  as defined in the glissade_velo_higher module.  
    !  The global matrix is sparse, but its nonzero elements are stored in 
    !  four dense matrices called Auu, Avv, Auv, and Avu.
    !  Each matrix has 3x3x3 = 27 potential nonzero elements per node (i,j,k).
    !
    !  The current preconditioning options are
    !  (0) no preconditioning
    !  (1) diagonal preconditioning
    !  (2) preconditioning using a physics-based SIA solver
    ! 
    !  For the dome test case with higher-order dynamics, option (2) is best. 
    !
    !  Here is a schematic of the method implemented below for solving Ax = b:
    !
    !  Set up preconditioner M
    !  work0 = (b,b)
    !  bb = global_sum(work0)
    !
    !  First pass of algorithm:
    !  halo_update(x)
    !  r = b - A*x
    !  halo_update(r)
    !  solve Mz = r for z
    !  work(1) = (r,z)
    !  d = z
    !  q = A*d
    !  work(2) = (d,q)
    !  halo_update(q)
    !  rho_old = global_sum(work(1))
    !  sigma = global_sum(work(2))
    !  alpha = rho_old/sigma
    !  x = x + alpha*d
    !  r = r - alpha*q
    !
    !  Iterative loop:
    !  while (not converged)
    !     solve Mz = r for z
    !     Az = A*z
    !     work(1) = (r,z)
    !     work(2) = (Az,z)
    !     halo_update(Az)
    !     rho = global_sum(work(1))
    !     delta = global_sum(work(2))
    !     beta = rho/rho_old
    !     sigma = delta - beta^2 * sigma
    !     alpha = rho/sigma
    !     rho_old = rho
    !     d = z + beta*d
    !     q = Az + beta*q
    !     x = x + alpha*d
    !     r = r - alpha*q
    !     if (time to check convergence) then
    !        r = b - A*x
    !        work0 = (r,r)
    !        halo_update(r)
    !        rr = global_sum(work0)
    !        if (sqrt(r,r)/sqrt(b,b) < tolerance) exit
    !     endif
    !  end while
    !
    !  where x = solution vector
    !        d = conjugate direction vector
    !        r = residual vector
    !        M = preconditioning matrix
    !    (r,z) = dot product of vectors r and z
    !            and similarly for (Az,z), etc.
    !       
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
                                ! velocity grid has dimensions (nx-1,ny-1)
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers (for scalars)

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F
 
    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       xu, xv             ! u and v components of solution (i.e., uvel and vvel)

    integer, intent(in)  ::   &
       precond           ! = 0 for no preconditioning
                         ! = 1 for diagonal preconditioning (best option for SSA-dominated flow)
                         ! = 2 for preconditioning with SIA solver (works well for SIA-dominated flow)

    real(dp), intent(out) ::  &
       err                               ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                            ! iterations needed to solution

    integer, intent(in), optional :: &
       itest_in, jtest_in, rtest_in      ! point for debugging diagnostics

    logical, intent(in), optional :: &
       verbose                           ! if true, print diagnostic output
    
    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer ::  i, j, k      ! grid indices
    integer ::  n            ! iteration counter

    real(dp) ::           &
       alpha,             &! rho/sigma = term in expression for new residual and solution
       beta,              &! eta1/eta0 = term in expression for new direction vector
       rho,               &! global sum of (r,z)
       rho_old,           &! old value of rho
       delta,             &! global sum of (s,z)
       sigma               ! delta - beta^2 * sigma

    real(dp), dimension(2) ::   &
       gsum                ! result of global sum for dot products
   
    ! vectors (each of these is split into u and v components)
    real(dp), dimension(nz,nx-1,ny-1) ::  &
       Adiagu, Adiagv,    &! diagonal terms of matrices Auu and Avv
       ru, rv,            &! residual vector (b-Ax)
       du, dv,            &! conjugate direction vector
       zu, zv,            &! solution of Mz = r
       qu, qv,            &! vector used to adjust residual, r -> r - alpha*q
       Azu, Azv,          &! result of matvec multiply A*z
       worku, workv        ! intermediate results

    real(dp), dimension(nz,nx-1,ny-1,2) ::  &
       work2u, work2v      ! intermediate results

    real(dp) ::  &
       rr,                &! dot product (r,r)
       bb,                &! dot product (b,b)
       L2_resid,          &! L2 norm of residual = sqrt(r,r)
       L2_rhs              ! L2 norm of rhs vector = sqrt(b,b)
                           ! solver is converged when L2_resid/L2_rhs < tolerance

    real(dp), dimension(-1:1,nz,nx-1,ny-1) ::  &
       Muu, Mvv            ! simplified SIA matrices for preconditioning

    integer :: itest, jtest, rtest

    if (present(itest_in)) then
       itest = itest_in
    else
       itest = nx/2
    endif

    if (present(itest_in)) then
       jtest = jtest_in
    else
       jtest = ny/2
    endif

    if (present(itest_in)) then
       rtest = rtest_in
    else
       rtest = 0
    endif

    if (present(verbose)) then
       verbose_pcg = verbose
    else
       verbose_pcg = .false.   ! for debugging
    endif

    if (verbose_pcg .and. main_task) then
       print*, 'Using native PCG solver (Chronopoulos-Gear)'
       print*, 'tolerance, maxiters, precond =', tolerance, maxiters, precond
    endif

    !---- Set up matrices for preconditioning

    call t_startf("pcg_precond_init")
    call setup_preconditioner_3d(nx,         ny,       &
                                 nz,                   &
                                 precond,    indxA,    &
                                 Auu,        Avv,      &
                                 Adiagu,     Adiagv,   &
                                 Muu,        Mvv)
    call t_stopf("pcg_precond_init")

    !---- Initialize scalars and vectors

    niters = maxiters 
    ru(:,:,:) = 0.d0
    rv(:,:,:) = 0.d0
    du(:,:,:) = 0.d0
    dv(:,:,:) = 0.d0
    zu(:,:,:) = 0.d0
    zv(:,:,:) = 0.d0
    qu(:,:,:) = 0.d0
    qv(:,:,:) = 0.d0
    Azu(:,:,:) = 0.d0
    Azv(:,:,:) = 0.d0
    worku(:,:,:) = 0.d0
    workv(:,:,:) = 0.d0
    work2u(:,:,:,:) = 0.d0
    work2v(:,:,:,:) = 0.d0

    !---- Compute the L2 norm of the RHS vectors
    !---- (Goal is to obtain L2_resid/L2_rhs < tolerance)

    call t_startf("pcg_dotprod")
    worku(:,:,:) = bu(:,:,:)*bu(:,:,:)    ! terms of dot product (b, b)
    workv(:,:,:) = bv(:,:,:)*bv(:,:,:)
    call t_stopf("pcg_dotprod")

    ! find global sum of the squared L2 norm

    call t_startf("pcg_glbsum_init")
    call global_sum_staggered(nx,     ny,     &
                              nz,     nhalo,  &
                              bb,             &
                              worku,  workv)
    call t_stopf("pcg_glbsum_init")

    ! take square root

    L2_rhs = sqrt(bb)       ! L2 norm of RHS

    !---------------------------------------------------------------
    ! First pass of algorithm
    !---------------------------------------------------------------

    ! Note: The matrix A must be complete for all rows corresponding to locally 
    !        owned nodes, and x must have the correct values in
    !        halo nodes bordering the locally owned nodes.
    !       Then z = Ax will be correct for locally owned nodes.

    !---- Halo update for x (initial guess for velocity solution)

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(xu)
    call staggered_parallel_halo(xv)
    call t_stopf("pcg_halo_init")

    !---- Compute A*x   (use z as a temp vector for A*x)

    call t_startf("pcg_matmult_init")
    call matvec_multiply_structured_3d(nx,        ny,            &
                                       nz,        nhalo,         &
                                       indxA,     active_vertex, &
                                       Auu,       Auv,           &
                                       Avu,       Avv,           &
                                       xu,        xv,            &
                                       zu,        zv)
    call t_stopf("pcg_matmult_init")

    !---- Compute the initial residual r = b - A*x
    !---- This is correct for locally owned nodes.

    call t_startf("pcg_vecupdate")
    ru(:,:,:) = bu(:,:,:) - zu(:,:,:)
    rv(:,:,:) = bv(:,:,:) - zv(:,:,:)
    call t_stopf("pcg_vecupdate")

    !---- Halo update for residual

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(ru)
    call staggered_parallel_halo(rv)
    call t_stopf("pcg_halo_init")

    !---- Compute (PC)r = solution z of Mz = r
    !---- Since r was just updated in halo, z is correct in halo

    ! From here on, call timers with 'iter' suffix because this can be considered the first iteration
    call t_startf("pcg_precond_iter")

    if (precond == 0) then      ! no preconditioning

       zu(:,:,:) = ru(:,:,:)         ! PC(r) = r     
       zv(:,:,:) = rv(:,:,:)         ! PC(r) = r    

    elseif (precond == 1 ) then  ! diagonal preconditioning

       do j = 1, ny-1
       do i = 1, nx-1
       do k = 1, nz
          if (Adiagu(k,i,j) /= 0.d0) then
             zu(k,i,j) = ru(k,i,j) / Adiagu(k,i,j)   ! PC(r), where PC is formed from diagonal elements of A
          else                                        
             zu(k,i,j) = 0.d0
          endif
          if (Adiagv(k,i,j) /= 0.d0) then
             zv(k,i,j) = rv(k,i,j) / Adiagv(k,i,j)  
          else                                        
             zv(k,i,j) = 0.d0
          endif
       enddo    ! k
       enddo    ! i
       enddo    ! j

    elseif (precond == 2) then   ! local vertical shallow-ice solver for preconditioning

       call easy_sia_solver(nx,   ny,   nz,        &
                            active_vertex,         &
                            Muu,  ru,   zu)      ! solve Muu*zu = ru for zu 

       call easy_sia_solver(nx,   ny,   nz,        &
                            active_vertex,         &
                            Mvv,  rv,   zv)      ! solve Mvv*zv = rv for zv

    endif    ! precond

    call t_stopf("pcg_precond_iter")

    !---- Compute intermediate result for dot product (r,z)

    call t_startf("pcg_dotprod")
    work2u(:,:,:,1) = ru(:,:,:) * zu(:,:,:)
    work2v(:,:,:,1) = rv(:,:,:) * zv(:,:,:)
    call t_stopf("pcg_dotprod")

    !---- Compute the conjugate direction vector d
    !---- Since z is correct in halo, so is d

    du(:,:,:) = zu(:,:,:)
    dv(:,:,:) = zv(:,:,:)

    !---- Compute q = A*d
    !---- q is correct for locally owned nodes

    call t_startf("pcg_matmult_iter")
    call matvec_multiply_structured_3d(nx,        ny,            &
                                       nz,        nhalo,         &
                                       indxA,     active_vertex, &
                                       Auu,       Auv,           &
                                       Avu,       Avv,           &
                                       du,        dv,            &
                                       qu,        qv)
    call t_stopf("pcg_matmult_iter")

    !---- Compute intermediate result for dot product (d,q) = (d,Ad)

    call t_startf("pcg_dotprod")
    work2u(:,:,:,2) = du(:,:,:) * qu(:,:,:)
    work2v(:,:,:,2) = dv(:,:,:) * qv(:,:,:)
    call t_stopf("pcg_dotprod")

    !---- Find global sums of (r,z) and (d,q)

    call t_startf("pcg_glbsum_iter")
    call global_sum_staggered(nx,     ny,     &
                              nz,     nhalo,  &
                              gsum,           &
                              work2u, work2v)
    call t_stopf("pcg_glbsum_iter")

    !---- Halo update for q

    call t_startf("pcg_halo_iter")
    call staggered_parallel_halo(qu)
    call staggered_parallel_halo(qv)
    call t_stopf("pcg_halo_iter")

    rho_old = gsum(1)      ! (r,z) = (r, (PC)r)
    sigma = gsum(2)        ! (d,q) = (d, Ad)

    alpha = rho_old/sigma

    if (alpha /= alpha) then  ! alpha is NaN
!!       write(6,*) 'rho_old, sigma, alpha:', rho_old, sigma, alpha
       call write_log('Chron_Gear PCG solver has failed, alpha = NaN', GM_FATAL)
    endif

    !---- Update solution and residual
    !---- These are correct in halo

    call t_startf("pcg_vecupdate")
    xu(:,:,:) = xu(:,:,:) + alpha*du(:,:,:)
    xv(:,:,:) = xv(:,:,:) + alpha*dv(:,:,:)

    ru(:,:,:) = ru(:,:,:) - alpha*qu(:,:,:)     ! q = A*d
    rv(:,:,:) = rv(:,:,:) - alpha*qv(:,:,:)
    call t_stopf("pcg_vecupdate")

    !---------------------------------------------------------------
    ! Iterate to solution
    !---------------------------------------------------------------

    iter_loop: do n = 2, maxiters  ! first iteration done above

       !---- Compute PC(r) = solution z of Mz = r
       !---- z is correct in halo

       call t_startf("pcg_precond_iter")

       if (precond == 0) then      ! no preconditioning

           zu(:,:,:) = ru(:,:,:)         ! PC(r) = r
           zv(:,:,:) = rv(:,:,:)         ! PC(r) = r    

       elseif (precond == 1 ) then  ! diagonal preconditioning

          do j = 1, ny-1
          do i = 1, nx-1
          do k = 1, nz
             if (Adiagu(k,i,j) /= 0.d0) then
                zu(k,i,j) = ru(k,i,j) / Adiagu(k,i,j)   ! PC(r), where PC is formed from diagonal elements of A
             else
                zu(k,i,j) = 0.d0
             endif
             if (Adiagv(k,i,j) /= 0.d0) then
                zv(k,i,j) = rv(k,i,j) / Adiagv(k,i,j)
             else
                zv(k,i,j) = 0.d0
             endif
          enddo    ! k
          enddo    ! i
          enddo    ! j

       elseif (precond == 2) then   ! local vertical shallow-ice solver for preconditioning

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Muu,  ru,   zu)      ! solve Muu*zu = ru for zu 

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Mvv,  rv,   zv)      ! solve Mvv*zv = rv for zv

       endif    ! precond

       call t_stopf("pcg_precond_iter")

       !---- Compute Az = A*z
       !---- This is the one matvec multiply required per iteration
       !---- Az is correct for local owned nodes and needs a halo update (below)

       call t_startf("pcg_matmult_iter")
       call matvec_multiply_structured_3d(nx,        ny,            &
                                          nz,        nhalo,         &
                                          indxA,     active_vertex, &
                                          Auu,       Auv,           &
                                          Avu,       Avv,           &
                                          zu,        zv,            &
                                          Azu,       Azv)
       call t_stopf("pcg_matmult_iter")

       !---- Compute intermediate results for the dot products (r,z) and (Az,z)

       call t_startf("pcg_dotprod")
       work2u(:,:,:,1) = ru(:,:,:)*zu(:,:,:)     ! terms of dot product (r,z)
       work2v(:,:,:,1) = rv(:,:,:)*zv(:,:,:)    

       work2u(:,:,:,2) = Azu(:,:,:)*zu(:,:,:)    ! terms of dot product (A*z,z)
       work2v(:,:,:,2) = Azv(:,:,:)*zv(:,:,:)    
       call t_stopf("pcg_dotprod")

       ! Take the global sums of (r,z) and (Az,z)
       ! Two sums are combined here for efficiency;
       ! this is the one MPI global reduction per iteration.

       call t_startf("pcg_glbsum_iter")
       call global_sum_staggered(nx,     ny,     &
                                 nz,     nhalo,  &
                                 gsum,           &
                                 work2u, work2v)
       call t_stopf("pcg_glbsum_iter")

       !---- Halo update for Az
       !---- This is the one halo update required per iteration

       call t_startf("pcg_halo_iter")
       call staggered_parallel_halo(Azu)
       call staggered_parallel_halo(Azv)
       call t_stopf("pcg_halo_iter")

       !---- Compute some scalars

       rho = gsum(1)        ! (r,z)
       delta = gsum(2)      ! (Az,z)

       beta = rho / rho_old
       sigma = delta - beta*rho/alpha
       alpha = rho / sigma
       rho_old = rho        ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i))

       if (alpha /= alpha) then  ! alpha is NaN
!!          write(6,*) 'rho, sigma, alpha:', rho, sigma, alpha
          call write_log('Chron-Gear PCG solver has failed, alpha = NaN', GM_FATAL)
       endif

       !---- Update d and q
       !---- These are correct in halo

       call t_startf("pcg_vecupdate")

       du(:,:,:) = zu(:,:,:) + beta*du(:,:,:)       ! d_(i+1) = PC(r_(i+1)) + beta_(i+1)*d_i
       dv(:,:,:) = zv(:,:,:) + beta*dv(:,:,:)       !
                                                    !                    (r_(i+1), PC(r_(i+1)))
                                                    ! where beta_(i+1) = --------------------  
                                                    !                        (r_i, PC(r_i)) 
       qu(:,:,:) = Azu(:,:,:) + beta*qu(:,:,:)
       qv(:,:,:) = Azv(:,:,:) + beta*qv(:,:,:)

       !---- Update solution and residual
       !---- These are correct in halo

       xu(:,:,:) = xu(:,:,:) + alpha*du(:,:,:)
       xv(:,:,:) = xv(:,:,:) + alpha*dv(:,:,:)

       ru(:,:,:) = ru(:,:,:) - alpha*qu(:,:,:)
       rv(:,:,:) = rv(:,:,:) - alpha*qv(:,:,:)

       call t_stopf("pcg_vecupdate")

       !---------------------------------------------------------------
       ! Convergence check every solve_ncheck iterations
       !---------------------------------------------------------------

       if (mod(n,solve_ncheck) == 0) then    ! use r = b - Ax

          !---- Compute z = A*x  (use z as a temp vector for A*x)
           
          call t_startf("pcg_matmult_resid")
          call matvec_multiply_structured_3d(nx,        ny,            &
                                             nz,        nhalo,         &
                                             indxA,     active_vertex, &
                                             Auu,       Auv,           &
                                             Avu,       Avv,           &
                                             xu,        xv,            &
                                             zu,        zv)
          call t_stopf("pcg_matmult_resid")

          !---- Compute residual r = b - A*x

          call t_startf("pcg_vecupdate")
          ru(:,:,:) = bu(:,:,:) - zu(:,:,:)
          rv(:,:,:) = bv(:,:,:) - zv(:,:,:)
          call t_stopf("pcg_vecupdate")

          !---- Compute dot product (r, r)

          call t_startf("pcg_dotprod")
          worku(:,:,:) = ru(:,:,:)*ru(:,:,:)
          workv(:,:,:) = rv(:,:,:)*rv(:,:,:)
          call t_stopf("pcg_dotprod")

          call t_startf("pcg_glbsum_resid")
          call global_sum_staggered(nx,     ny,       &
                                    nz,     nhalo,    &
                                    rr,               &
                                    worku, workv)
          call t_stopf("pcg_glbsum_resid")

          L2_resid = sqrt(rr)          ! L2 norm of residual
          err = L2_resid/L2_rhs        ! normalized error

          if (verbose_pcg .and. main_task) then
!             print*, ' '
!             print*, 'iter, L2_resid, error =', n, L2_resid, err
          endif

          if (err < tolerance) then
             niters = n
             exit iter_loop
          endif            

          !---- Update residual in halo for next iteration

          call t_startf("pcg_halo_resid")
          call staggered_parallel_halo(ru)
          call staggered_parallel_halo(rv)
          call t_stopf("pcg_halo_resid")

       endif    ! solve_ncheck

    enddo iter_loop

    !WHL - Without good preconditioning, convergence can be slow, but the solution after maxiters might be good enough.
 
    if (niters == maxiters) then
       if (verbose_pcg .and. main_task) then
          print*, 'Glissade PCG solver not converged'
          print*, 'niters, err, tolerance:', niters, err, tolerance
       endif
    endif

  end subroutine pcg_solver_chrongear_3d

!****************************************************************************

  subroutine pcg_solver_chrongear_2d(nx,        ny,            &
                                     nhalo,                    &
                                     indxA_2d,  active_vertex, &
                                     Auu,       Auv,           &
                                     Avu,       Avv,           &
                                     bu,        bv,            &
                                     xu,        xv,            &
                                     precond,   err,           &
                                     niters,                   &
                                     itest_in,  jtest_in,  rtest_in,   &
                                     verbose) 

    !---------------------------------------------------------------
    !  This subroutine uses a Chronopoulos-Gear preconditioned conjugate-gradient
    !  algorithm to solve the equation $Ax=b$. (See references in subroutine above.)
    !
    !  It is similar to subroutine pcg_solver_chrongear_3d, but modified
    !  to solve for x and y at a single horizontal level, as in the
    !  shallow-shelf approximation.  See the comments in that subroutine
    !  (above) for more details on data structure and solver methods.
    !
    !  Input and output arrays are located on a structured (i,j) grid 
    !  as defined in the glissade_velo_higher module.  The global matrix 
    !  is sparse, but its nonzero element are stored in four dense matrices 
    !  called Auu, Avv, Auv, and Avu. Each matrix has 3x3 = 9 potential 
    !  nonzero elements per node (i,j).
    !
    !  The current preconditioning options for the 2D solver are
    !  (0) no preconditioning
    !  (1) diagonal preconditioning
    !  (3) tridiagonal preconditioning  !TODO - Add tridiagonal option for standard 2D solver
    !
    !  The SIA-based preconditioning option is not available for a 2D solve.
    !
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
                                ! velocity grid has dimensions (nx-1,ny-1)
       nhalo                    ! number of halo layers (for scalars)

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F
 
    real(dp), dimension(9,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 9 (node and its nearest neighbors in x and y direction)
                              ! other dimensions = (x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    real(dp), dimension(nx-1,ny-1), intent(inout) ::   &
       xu, xv             ! u and v components of solution (i.e., uvel and vvel)

    integer, intent(in)  ::   &
       precond           ! = 0 for no preconditioning
                         ! = 1 for diagonal preconditioning (best option for SSA-dominated flow)

    real(dp), intent(out) ::  &
       err                               ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                            ! iterations needed to solution

    integer, intent(in), optional :: &
       itest_in, jtest_in, rtest_in      ! point for debugging diagnostics

    logical, intent(in), optional :: &
       verbose                           ! if true, print diagnostic output
    
    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer ::  i, j         ! grid indices
    integer ::  m            ! matrix element index
    integer ::  n            ! iteration counter

    real(dp) ::           &
       alpha,             &! rho/sigma = term in expression for new residual and solution
       beta,              &! eta1/eta0 = term in expression for new direction vector
       rho,               &! global sum of (r,z)
       rho_old,           &! old value of rho
       delta,             &! global sum of (s,z)
       sigma               ! delta - beta^2 * sigma

    real(dp), dimension(2) ::   &
       gsum                ! result of global sum for dot products
   
    ! diagonal matrix elements
    real(dp), dimension(nx-1,ny-1) ::  &
       Adiagu, Adiagv      ! diagonal terms of matrices Auu and Avv

    ! tridiagonal matrix elements
    real(dp), dimension(:,:), allocatable :: &
         Asubdiag_u, Adiag_u, Asupdiag_u,  &  ! matrix entries from Auu for tridiagonal preconditioning
         Asubdiag_v, Adiag_v, Asupdiag_v      ! matrix entries from Avv for tridiagonal preconditioning

    real(dp), dimension(:,:), allocatable :: &
         omega_u, omega_v,         & ! work arrays for tridiagonal solve
         denom_u, denom_v,         &
         xuh_u,   xuh_v,           &
         xlh_u,   xlh_v,           &
         b_u,     b_v,             & ! rhs for tridiagonal solve
         x_u,     x_v                ! solution for tridiagonal solve

    ! vectors (each of these is split into u and v components)
    real(dp), dimension(nx-1,ny-1) ::  &
       ru, rv,            &! residual vector (b-Ax)
       du, dv,            &! conjugate direction vector
       zu, zv,            &! solution of Mz = r
       qu, qv,            &! vector used to adjust residual, r -> r - alpha*q
       Azu, Azv,          &! result of matvec multiply A*z
       worku, workv        ! intermediate results

    real(dp), dimension(nx-1,ny-1,2) ::  &
       work2u, work2v      ! intermediate results

    real(dp) ::  &
       rr,                &! dot product (r,r)
       bb,                &! dot product (b,b)
       L2_resid,          &! L2 norm of residual = sqrt(r,r)
       L2_rhs              ! L2 norm of rhs vector = sqrt(b,b)
                           ! solver is converged when L2_resid/L2_rhs < tolerance

    integer :: ilocal, jlocal   ! number of locally owned vertices in each direction

    integer :: itest, jtest, rtest
    integer :: ii, jj

    !WHL - debug
    real(dp) :: usum, usum_global, vsum, vsum_global
    logical, parameter :: &
!!         use_serial_tridiag_solver = .true.  ! if true, use the serial solver when tasks = 1
         use_serial_tridiag_solver = .false.   ! if false, use the parallel solver for any number of tasks, including tasks = 1

    if (present(itest_in)) then
       itest = itest_in
    else
       itest = nx/2
    endif

    if (present(itest_in)) then
       jtest = jtest_in
    else
       jtest = ny/2
    endif

    if (present(itest_in)) then
       rtest = rtest_in
    else
       rtest = 0
    endif

    if (present(verbose)) then
       verbose_pcg = verbose
    else
       verbose_pcg = .false.   ! for debugging
    endif

    !WHL - debug
!!    verbose_pcg = .true.

    if (verbose_pcg .and. main_task) then
       print*, 'Using native PCG solver (Chronopoulos-Gear)'
       print*, 'tolerance, maxiters, precond =', tolerance, maxiters, precond
    endif

    ! Compute array sizes for locally owned vertices
    ilocal = staggered_ihi - staggered_ilo + 1
    jlocal = staggered_jhi - staggered_jlo + 1

    !---- Set up matrices for preconditioning

    call t_startf("pcg_precond_init")

    if (precond == HO_PRECOND_DIAG) then

       call setup_preconditioner_diag_2d(nx,       ny,         &
                                         indxA_2d,             &
                                         Auu,      Avv,        &
                                         Adiagu,   Adiagv)

       !WHL - debug
       if (verbose_pcg .and. this_rank == rtest) then
          i = itest
          j = jtest
          print*, 'i, j, r =', i, j, this_rank
          print*, 'Au diag =', Adiagu(i,j)
          print*, 'Av diag =', Adiagv(i,j)
       endif

    elseif (precond == HO_PRECOND_TRIDIAG) then

       ! Allocate tridiagonal matrices
       ! Note: (i,j) indices are switced for the A_v matrices to reduce striding.

       allocate(Adiag_u   (ilocal,jlocal))
       allocate(Asubdiag_u(ilocal,jlocal))
       allocate(Asupdiag_u(ilocal,jlocal))

       allocate(Adiag_v   (jlocal,ilocal))
       allocate(Asubdiag_v(jlocal,ilocal))
       allocate(Asupdiag_v(jlocal,ilocal))

       allocate(omega_u(ilocal,jlocal))
       allocate(denom_u(ilocal,jlocal))
       allocate(xuh_u(ilocal,jlocal))
       allocate(xlh_u(ilocal,jlocal))
       allocate(b_u(ilocal,jlocal))
       allocate(x_u(ilocal,jlocal))

       allocate(omega_v(jlocal,ilocal))
       allocate(denom_v(jlocal,ilocal))
       allocate(xuh_v(jlocal,ilocal))
       allocate(xlh_v(jlocal,ilocal))
       allocate(b_v(jlocal,ilocal))
       allocate(x_v(jlocal,ilocal))

       ! Compute the entries of the tridiagonal matrices

       ! Extract tridiagonal matrix entries from Auu
       do j = 1, jlocal
          jj = j + staggered_jlo - 1
          do i = 1, ilocal
             ii = i + staggered_ilo - 1
             Asubdiag_u(i,j) = Auu(indxA_2d(-1,0),ii,jj)   ! subdiagonal elements
             Adiag_u   (i,j) = Auu(indxA_2d( 0,0),ii,jj)   ! diagonal elements
             Asupdiag_u(i,j) = Auu(indxA_2d( 1,0),ii,jj)   ! superdiagonal elements
          enddo
       enddo

       ! compute work arrays for the u solve in each matrix row
       call setup_preconditioner_tridiag_2d(ilocal,       jlocal,      &
!!                                            itest, jtest, rtest,       &
                                            itest - staggered_ilo + 1, &  ! itest referenced to (ilocal,jlocal) coordinates
                                            jtest - staggered_jlo + 1, &  ! jtest referenced to (ilocal,jlocal) coordinates
                                            rtest,                     &
                                            Adiag_u,                   &
                                            Asubdiag_u,   Asupdiag_u,  &
                                            omega_u,      denom_u,     &
                                            xuh_u,        xlh_u)

       ! Extract tridiagonal matrix entries from Avv
       do i = 1, ilocal
          ii = i + staggered_ilo - 1
          do j = 1, jlocal
             jj = j + staggered_jlo - 1
             Asubdiag_v(j,i) = Avv(indxA_2d(0,-1),ii,jj)   ! subdiagonal elements
             Adiag_v   (j,i) = Avv(indxA_2d(0, 0),ii,jj)   ! diagonal elements
             Asupdiag_v(j,i) = Avv(indxA_2d(0, 1),ii,jj)   ! superdiagonal elements
          enddo
       enddo

       ! compute work arrays for the v solve in each matrix column
       ! Note: The *_v arrays have dimensions (jlocal,ilocal) to reduce strides

       call setup_preconditioner_tridiag_2d(jlocal,       ilocal,      &
!!                                            itest, jtest, rtest,       &
                                            jtest - staggered_jlo + 1, &  ! jtest referenced to (jlocal,ilocal) coordinates
                                            itest - staggered_ilo + 1, &  ! itest referenced to (jlocal,ilocal) coordinates
                                            rtest,                     &
                                            Adiag_v,                   &
                                            Asubdiag_v,   Asupdiag_v,  &
                                            omega_v,      denom_v,     &
                                            xuh_v,        xlh_v)

    else    ! no preconditioner

       if (verbose_pcg .and. main_task) then
          print*, 'Using no preconditioner'
       endif

    endif   ! precond

    call t_stopf("pcg_precond_init")

    !---- Initialize scalars and vectors

    niters = maxiters
    ru(:,:) = 0.d0
    rv(:,:) = 0.d0
    du(:,:) = 0.d0
    dv(:,:) = 0.d0
    zu(:,:) = 0.d0
    zv(:,:) = 0.d0
    qu(:,:) = 0.d0
    qv(:,:) = 0.d0
    Azu(:,:) = 0.d0
    Azv(:,:) = 0.d0
    worku(:,:) = 0.d0
    workv(:,:) = 0.d0
    work2u(:,:,:) = 0.d0
    work2v(:,:,:) = 0.d0

    !---- Compute the L2 norm of the RHS vectors
    !---- (Goal is to obtain L2_resid/L2_rhs < tolerance)

    call t_startf("pcg_dotprod")
    worku(:,:) = bu(:,:)*bu(:,:)    ! terms of dot product (b, b)
    workv(:,:) = bv(:,:)*bv(:,:)
    call t_stopf("pcg_dotprod")

    ! find global sum of the squared L2 norm

    call t_startf("pcg_glbsum_init")
    call global_sum_staggered(nx,     ny,     &
                              nhalo,  bb,     &
                              worku,  workv)
    call t_stopf("pcg_glbsum_init")

    ! take square root

    L2_rhs = sqrt(bb)       ! L2 norm of RHS

    !---------------------------------------------------------------
    ! First pass of algorithm
    !---------------------------------------------------------------

    ! Note: The matrix A must be complete for all rows corresponding to locally 
    !        owned nodes, and x must have the correct values in
    !        halo nodes bordering the locally owned nodes.
    !       Then z = Ax will be correct for locally owned nodes.

    !---- Halo update for x (initial guess for velocity solution)

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(xu)
    call staggered_parallel_halo(xv)
    call t_stopf("pcg_halo_init")

    !---- Compute A*x   (use z as a temp vector for A*x)

    call t_startf("pcg_matmult_init")
    call matvec_multiply_structured_2d(nx,        ny,            &
                                       nhalo,                    &
                                       indxA_2d,  active_vertex, &
                                       Auu,       Auv,           &
                                       Avu,       Avv,           &
                                       xu,        xv,            &
                                       zu,        zv)
    call t_stopf("pcg_matmult_init")

    !---- Compute the initial residual r = b - A*x
    !---- This is correct for locally owned nodes.

    call t_startf("pcg_vecupdate")
    ru(:,:) = bu(:,:) - zu(:,:)
    rv(:,:) = bv(:,:) - zv(:,:)
    call t_stopf("pcg_vecupdate")

    !---- Halo update for residual

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(ru)
    call staggered_parallel_halo(rv)
    call t_stopf("pcg_halo_init")

    !---- Compute (PC)r = solution z of Mz = r
    !---- Since r was just updated in halo, z is correct in halo

    ! From here on, call timers with 'iter' suffix because this can be considered the first iteration
    call t_startf("pcg_precond_iter")

    if (precond == 0) then      ! no preconditioning

       zu(:,:) = ru(:,:)         ! PC(r) = r     
       zv(:,:) = rv(:,:)         ! PC(r) = r    

    elseif (precond == 1 ) then  ! diagonal preconditioning

       ! Solve Mz = r, which M is a diagonal matrix
       do j = 1, ny-1
       do i = 1, nx-1
          if (Adiagu(i,j) /= 0.d0) then
             zu(i,j) = ru(i,j) / Adiagu(i,j)   ! PC(r), where PC is formed from diagonal elements of A
          else                                        
             zu(i,j) = 0.d0
          endif
          if (Adiagv(i,j) /= 0.d0) then
             zv(i,j) = rv(i,j) / Adiagv(i,j)  
          else                                        
             zv(i,j) = 0.d0
          endif
       enddo    ! i
       enddo    ! j

       !WHL - debug
       if (verbose_pcg .and. this_rank == rtest) then
          i = itest
!          print*, ' '
!          print*, 'zv solve with diagonal precond, this_rank, i =', this_rank, i
!          print*, 'j, active, Adiagv, rv, zv, xv:'
!          do j = staggered_jhi, staggered_jlo, -1
!             write(6,'(i4, l4, 2f12.3, e12.3, f12.3)') j, active_vertex(i,j), Adiagv(i,j), rv(i,j), zv(i,j), xv(i,j)
!          enddo
       endif

    elseif (precond == 3 ) then  ! tridiagonal preconditioning

       ! Solve M*z = r, where M is a tridiagonal matrix

       if (tasks == 1 .and. use_serial_tridiag_solver) then

          call tridiag_solver_serial_2d(nx,           ny,         &
                                        itest, jtest, rtest,      &
                                        ilocal,       jlocal,     &
                                        Adiag_u,      Adiag_v,    &  ! entries of preconditioning matrix
                                        Asubdiag_u,   Asubdiag_v, &
                                        Asupdiag_u,   Asupdiag_v, &
                                        ru,           rv,         &  ! right hand side
                                        zu,           zv)            ! solution

       else

          ! convert ru(nx-1,ny-1) to b_u(ilocal,jlocal)

          do j = 1, jlocal
             jj = j + staggered_jlo - 1
             do i = 1, ilocal
                ii = i + staggered_ilo - 1
                b_u(i,j) = ru(ii,jj)
             enddo
          enddo

          call tridiag_solver_parallel_2d(ilocal,       jlocal,      &
                                          'row',                     &  ! tridiagonal solve for each row
!!                                          itest, jtest, rtest,      &
                                          itest - staggered_ilo + 1, &  ! itest referenced to (ilocal,jlocal) coordinates
                                          jtest - staggered_jlo + 1, &  ! jtest referenced to (ilocal,jlocal) coordinates
                                          rtest,                     &
                                          Adiag_u,                   &
                                          Asubdiag_u,   Asupdiag_u,  &
                                          omega_u,      denom_u,     &
                                          xuh_u,        xlh_u,       &
                                          b_u,          x_u)

          ! convert x_u(ilocal,jlocal) to zu(nx-1,ny-1)
          zu(:,:) = 0.0d0
          do j = 1, jlocal
             jj = j + staggered_jlo - 1
             do i = 1, ilocal
                ii = i + staggered_ilo - 1
                zu(ii,jj) = x_u(i,j)
             enddo
          enddo

          ! convert rv(nx-1,ny-1) to b_v(jlocal,ilocal)

          do i = 1, ilocal
             ii = i + staggered_ilo - 1
             do j = 1, jlocal
                jj = j + staggered_jlo - 1
                b_v(j,i) = rv(ii,jj)
             enddo
          enddo

          call tridiag_solver_parallel_2d(jlocal,       ilocal,      &
                                          'col',                     &  ! tridiagonal solve for each row
!!                                          itest, jtest, rtest,      &
                                          jtest - staggered_jlo + 1, &  ! jtest referenced to (jlocal,ilocal) coordinates
                                          itest - staggered_ilo + 1, &  ! itest referenced to (jlocal,ilocal) coordinates
                                          rtest,                     &
                                          Adiag_v,                   &
                                          Asubdiag_v,   Asupdiag_v,  &
                                          omega_v,      denom_v,     &
                                          xuh_v,        xlh_v,       &
                                          b_v,          x_v)

          ! convert x_v(jlocal,ilocal) to zv(nx-1,ny-1)

          zv(:,:) = 0.0d0
          do i = 1, ilocal
             ii = i + staggered_ilo - 1
             do j = 1, jlocal
                jj = j + staggered_jlo - 1
                zv(ii,jj) = x_v(j,i)
             enddo
          enddo

       endif  ! tasks = 1 & serial tridiag

       !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
       !TODO: See whether tridiag_solver_parallel_2d could be modified to provide zu and zv in halo cells?
       call staggered_parallel_halo(zu)
       call staggered_parallel_halo(zv)

    endif    ! precond

    call t_stopf("pcg_precond_iter")

    !---- Compute intermediate result for dot product (r,z)

    call t_startf("pcg_dotprod")
    work2u(:,:,1) = ru(:,:) * zu(:,:)
    work2v(:,:,1) = rv(:,:) * zv(:,:)
    call t_stopf("pcg_dotprod")

    !---- Compute the conjugate direction vector d
    !---- Since z is correct in halo, so is d

    du(:,:) = zu(:,:)
    dv(:,:) = zv(:,:)

    !---- Compute q = A*d
    !---- q is correct for locally owned nodes, provided d extends one layer into the halo

    call t_startf("pcg_matmult_iter")
    call matvec_multiply_structured_2d(nx,        ny,            &
                                       nhalo,                    &
                                       indxA_2d,  active_vertex, &
                                       Auu,       Auv,           &
                                       Avu,       Avv,           &
                                       du,        dv,            &
                                       qu,        qv)
    call t_stopf("pcg_matmult_iter")

    !WHL - debug
    usum = sum(qu(staggered_ilo:staggered_ihi,staggered_jlo:staggered_jhi))
    usum_global = parallel_reduce_sum(usum)
    vsum = sum(qv(staggered_ilo:staggered_ihi,staggered_jlo:staggered_jhi))
    vsum_global = parallel_reduce_sum(vsum)

    if (verbose_pcg .and. main_task) then
      print*, 'Prep: sum(qu), sum(qv) =', usum_global, vsum_global
    endif

    !---- Compute intermediate result for dot product (d,q) = (d,Ad)

    call t_startf("pcg_dotprod")
    work2u(:,:,2) = du(:,:) * qu(:,:)
    work2v(:,:,2) = dv(:,:) * qv(:,:)
    call t_stopf("pcg_dotprod")

    !---- Find global sums of (r,z) and (d,q)

    call t_startf("pcg_glbsum_iter")
    call global_sum_staggered(nx,     ny,     &
                              nhalo,  gsum,   &
                              work2u, work2v)
    call t_stopf("pcg_glbsum_iter")

    if (verbose_pcg .and. main_task) then
       print*, 'Prep: gsum(1), gsum(2) =', gsum(1), gsum(2)
    endif

    !---- Halo update for q

    call t_startf("pcg_halo_iter")
    call staggered_parallel_halo(qu)
    call staggered_parallel_halo(qv)
    call t_stopf("pcg_halo_iter")

    rho_old = gsum(1)      ! (r,z) = (r, (PC)r)
    sigma = gsum(2)        ! (d,q) = (d, Ad)

    alpha = rho_old/sigma

    if (alpha /= alpha) then  ! alpha is NaN
!!       write(6,*) 'rho_old, sigma, alpha:', rho_old, sigma, alpha
       call write_log('Chron_Gear PCG solver has failed, alpha = NaN', GM_FATAL)
    endif

    !---- Update solution and residual
    !---- These are correct in halo

    call t_startf("pcg_vecupdate")
    xu(:,:) = xu(:,:) + alpha*du(:,:)
    xv(:,:) = xv(:,:) + alpha*dv(:,:)

    ru(:,:) = ru(:,:) - alpha*qu(:,:)     ! q = A*d
    rv(:,:) = rv(:,:) - alpha*qv(:,:)
    call t_stopf("pcg_vecupdate")

    !---------------------------------------------------------------
    ! Iterate to solution
    !---------------------------------------------------------------

    iter_loop: do n = 2, maxiters  ! first iteration done above

       !---- Compute PC(r) = solution z of Mz = r
       !---- z is correct in halo

       call t_startf("pcg_precond_iter")

       if (precond == 0) then      ! no preconditioning

           zu(:,:) = ru(:,:)         ! PC(r) = r
           zv(:,:) = rv(:,:)         ! PC(r) = r    

       elseif (precond == 1 ) then  ! diagonal preconditioning

          do j = 1, ny-1
          do i = 1, nx-1
             if (Adiagu(i,j) /= 0.d0) then
                zu(i,j) = ru(i,j) / Adiagu(i,j)   ! PC(r), where PC is formed from diagonal elements of A
             else
                zu(i,j) = 0.d0
             endif
             if (Adiagv(i,j) /= 0.d0) then
                zv(i,j) = rv(i,j) / Adiagv(i,j)
             else
                zv(i,j) = 0.d0
             endif
          enddo    ! i
          enddo    ! j

       elseif (precond == 3) then   ! tridiagonal preconditioning

          if (tasks == 1 .and. use_serial_tridiag_solver) then

             ! Solve M*z = r, where M is a tridiagonal matrix
             call tridiag_solver_serial_2d(nx,           ny,         &
                                           itest, jtest, rtest,      &
                                           ilocal,       jlocal,     &
                                           Adiag_u,      Adiag_v,    &  ! entries of preconditioning matrix
                                           Asubdiag_u,   Asubdiag_v, &
                                           Asupdiag_u,   Asupdiag_v, &
                                           ru,           rv,         &  ! right hand side
                                           zu,           zv)            ! solution

          else   ! use parallel solver

             ! convert ru(nx-1,ny-1) to b_u(ilocal,jlocal)

             do j = 1, jlocal
                jj = j + staggered_jlo - 1
                do i = 1, ilocal
                   ii = i + staggered_ilo - 1
                   b_u(i,j) = ru(ii,jj)
                enddo
             enddo

             call tridiag_solver_parallel_2d(ilocal,       jlocal,      &
                                             'row',                     &  ! tridiagonal solve for each row
!!                                             itest, jtest, rtest,      &
                                             itest - staggered_ilo + 1, &  ! itest referenced to (ilocal,jlocal) coordinates
                                             jtest - staggered_jlo + 1, &  ! jtest referenced to (ilocal,jlocal) coordinates
                                             rtest,                     &
                                             Adiag_u,                   &
                                             Asubdiag_u,   Asupdiag_u,  &
                                             omega_u,      denom_u,     &
                                             xuh_u,        xlh_u,       &
                                             b_u,          x_u)

             ! convert x_u(ilocal,jlocal) to zu(nx-1,ny-1)
             zu(:,:) = 0.0d0
             do j = 1, jlocal
                jj = j + staggered_jlo - 1
                do i = 1, ilocal
                   ii = i + staggered_ilo - 1
                   zu(ii,jj) = x_u(i,j)
                enddo
             enddo

             ! convert rv(nx-1,ny-1) to b_v(jlocal,ilocal)

             do i = 1, ilocal
                ii = i + staggered_ilo - 1
                do j = 1, jlocal
                   jj = j + staggered_jlo - 1
                   b_v(j,i) = rv(ii,jj)
                enddo
             enddo

             call tridiag_solver_parallel_2d(jlocal,       ilocal,      &
                                             'col',                     &  ! tridiagonal solve for each row
!!                                             itest, jtest, rtest,      &
                                             jtest - staggered_jlo + 1, &  ! jtest referenced to (jlocal,ilocal) coordinates
                                             itest - staggered_ilo + 1, &  ! itest referenced to (jlocal,ilocal) coordinates
                                             rtest,                     &
                                             Adiag_v,                   &
                                             Asubdiag_v,   Asupdiag_v,  &
                                             omega_v,      denom_v,     &
                                             xuh_v,        xlh_v,       &
                                             b_v,          x_v)

             ! convert x_v(jlocal,ilocal) to zv(nx-1,ny-1)

             zv(:,:) = 0.0d0
             do i = 1, ilocal
                ii = i + staggered_ilo - 1
                do j = 1, jlocal
                   jj = j + staggered_jlo - 1
                   zv(ii,jj) = x_v(j,i)
                enddo
             enddo

             !Note: Need zu and zv in a row of halo cells so that A*z is correct in all locally owned cells
             call staggered_parallel_halo(zu)
             call staggered_parallel_halo(zv)

          endif   ! tasks = 1

       endif    ! precond

       call t_stopf("pcg_precond_iter")

       !---- Compute Az = A*z
       !---- This is the one matvec multiply required per iteration
       !---- Az is correct for locally owned nodes and needs a halo update (below)

       call t_startf("pcg_matmult_iter")
       call matvec_multiply_structured_2d(nx,        ny,            &
                                          nhalo,                    &
                                          indxA_2d,  active_vertex, &
                                          Auu,       Auv,           &
                                          Avu,       Avv,           &
                                          zu,        zv,            &
                                          Azu,       Azv)
       call t_stopf("pcg_matmult_iter")

       !---- Compute intermediate results for the dot products (r,z) and (Az,z)

       call t_startf("pcg_dotprod")
       work2u(:,:,1) = ru(:,:)*zu(:,:)     ! terms of dot product (r,z)
       work2v(:,:,1) = rv(:,:)*zv(:,:)    

       work2u(:,:,2) = Azu(:,:)*zu(:,:)    ! terms of dot product (A*z,z)
       work2v(:,:,2) = Azv(:,:)*zv(:,:)    
       call t_stopf("pcg_dotprod")

       ! Take the global sums of (r,z) and (Az,z)
       ! Two sums are combined here for efficiency;
       ! this is the one MPI global reduction per iteration.

       call t_startf("pcg_glbsum_iter")
       call global_sum_staggered(nx,     ny,     &
                                 nhalo,  gsum,   &
                                 work2u, work2v)
       call t_stopf("pcg_glbsum_iter")

       !---- Halo update for Az

       call t_startf("pcg_halo_iter")
       call staggered_parallel_halo(Azu)
       call staggered_parallel_halo(Azv)
       call t_stopf("pcg_halo_iter")

       if (verbose_pcg .and. main_task) then
          print*, 'iter, gsum(1), gsum(2) =', n, gsum(1), gsum(2)
       endif

       !---- Compute some scalars

       rho = gsum(1)        ! (r,z)
       delta = gsum(2)      ! (Az,z)

       beta = rho / rho_old
       sigma = delta - beta*rho/alpha
       alpha = rho / sigma
       rho_old = rho        ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i))

       if (alpha /= alpha) then  ! alpha is NaN
!!          write(6,*) 'rho, sigma, alpha:', rho, sigma, alpha
          call write_log('Chron_Gear PCG solver has failed, alpha = NaN', GM_FATAL)
       endif

       !---- Update d and q
       !---- These are correct in halo

       call t_startf("pcg_vecupdate")

       du(:,:) = zu(:,:) + beta*du(:,:)       ! d_(i+1) = PC(r_(i+1)) + beta_(i+1)*d_i
       dv(:,:) = zv(:,:) + beta*dv(:,:)       !
                                              !                    (r_(i+1), PC(r_(i+1)))
                                              ! where beta_(i+1) = --------------------  
                                              !                        (r_i, PC(r_i)) 
       qu(:,:) = Azu(:,:) + beta*qu(:,:)
       qv(:,:) = Azv(:,:) + beta*qv(:,:)

       !---- Update solution and residual
       !---- These are correct in halo

       xu(:,:) = xu(:,:) + alpha*du(:,:)
       xv(:,:) = xv(:,:) + alpha*dv(:,:)

       ru(:,:) = ru(:,:) - alpha*qu(:,:)
       rv(:,:) = rv(:,:) - alpha*qv(:,:)

       call t_stopf("pcg_vecupdate")

       !---------------------------------------------------------------
       ! Convergence check every solve_ncheck iterations
       !---------------------------------------------------------------

       if (mod(n,solve_ncheck) == 0) then    ! use r = b - Ax

          !---- Compute z = A*x  (use z as a temp vector for A*x)
           
          call t_startf("pcg_matmult_resid")
          call matvec_multiply_structured_2d(nx,        ny,            &
                                             nhalo,                    &
                                             indxA_2d,  active_vertex, &
                                             Auu,       Auv,           &
                                             Avu,       Avv,           &
                                             xu,        xv,            &
                                             zu,        zv)
          call t_stopf("pcg_matmult_resid")

          !---- Compute residual r = b - A*x

          call t_startf("pcg_vecupdate")
          ru(:,:) = bu(:,:) - zu(:,:)
          rv(:,:) = bv(:,:) - zv(:,:)
          call t_stopf("pcg_vecupdate")

          !---- Compute dot product (r, r)

          call t_startf("pcg_dotprod")
          worku(:,:) = ru(:,:)*ru(:,:)
          workv(:,:) = rv(:,:)*rv(:,:)
          call t_stopf("pcg_dotprod")

          call t_startf("pcg_glbsum_resid")
          call global_sum_staggered(nx,     ny,       &
                                    nhalo,  rr,       &
                                    worku,  workv)
          call t_stopf("pcg_glbsum_resid")

          L2_resid = sqrt(rr)          ! L2 norm of residual
          err = L2_resid/L2_rhs        ! normalized error

          if (verbose_pcg .and. main_task) then
!             print*, ' '
!             print*, 'iter, L2_resid, error =', n, L2_resid, err
          endif

          if (err < tolerance) then
             niters = n
             exit iter_loop
          endif            

          !---- Update residual in halo for next iteration

          call t_startf("pcg_halo_resid")
          call staggered_parallel_halo(ru)
          call staggered_parallel_halo(rv)
          call t_stopf("pcg_halo_resid")

       endif    ! solve_ncheck

    enddo iter_loop

    ! Clean up
    if (allocated(Adiag_u)) deallocate(Adiag_u, Asubdiag_u, Asupdiag_u)
    if (allocated(Adiag_v)) deallocate(Adiag_v, Asubdiag_v, Asupdiag_v)
    if (allocated(omega_u)) deallocate(omega_u, denom_u, xuh_u, xlh_u, b_u, x_u)
    if (allocated(omega_v)) deallocate(omega_v, denom_v, xuh_v, xlh_v, b_v, x_v)

    !WHL - Without good preconditioning, convergence can be slow, but the solution after maxiters might be good enough.
 
    if (niters == maxiters) then
       if (verbose_pcg .and. main_task) then
          print*, 'Glissade PCG solver not converged'
          print*, 'niters, err, tolerance:', niters, err, tolerance
       endif
    endif

  end subroutine pcg_solver_chrongear_2d

!****************************************************************************

  subroutine setup_preconditioner_3d(nx,         ny,       &
                                     nz,                   &
                                     precond,    indxA,    &
                                     Auu,        Avv,      &
                                     Adiagu,     Adiagv,   &
                                     Muu,        Mvv)

    ! Set up preconditioning matrices using one of several options

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
                              ! velocity grid has dimensions (nx-1,ny-1)
       nz                     ! number of vertical levels where velocity is computed

    integer, intent(in)  ::   &
       precond                ! = 0 for no preconditioning
                              ! = 1 for diagonal preconditioning
                              ! = 2 for preconditioning with SIA solver (works well for SIA-dominated flow)

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                  ! maps relative (x,y,z) coordinates to an index between 1 and 27

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Avv               ! two out of the four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::   &
       Adiagu, Adiagv         ! matrices for diagonal preconditioning 

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(out) ::   &
       Muu, Mvv               ! preconditioning matrices based on shallow-ice approximation

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, k, m

    ! Initialize

    Adiagu(:,:,:) = 0.d0
    Adiagv(:,:,:) = 0.d0
    Muu (:,:,:,:) = 0.d0    
    Mvv (:,:,:,:) = 0.d0    

    if (precond == HO_PRECOND_DIAG) then    ! form diagonal matrix for preconditioning

       if (verbose_pcg .and. main_task) then
          print*, 'Using diagonal matrix for preconditioning'
       endif  ! verbose_pcg

       m = indxA(0,0,0)
       Adiagu(:,:,:) = Auu(m,:,:,:)
       Adiagv(:,:,:) = Avv(m,:,:,:)

    elseif (precond == HO_PRECOND_SIA) then  ! form SIA matrices Muu and Mvv with vertical coupling only

       if (verbose_pcg .and. main_task) then
          print*, 'Using shallow-ice preconditioner'
       endif  ! verbose_pcg

       do j = 1, ny-1
       do i = 1, nx-1
       do k = 1, nz
           ! Remove horizontal coupling by using only the iA=0, jA=0 term in each layer.

            !WHL - Summing over the terms in each layer does not work for simple shelf problems
            !      because the matrix can be singular.
!           Muu(-1,k,i,j) = sum(Auu(1:9,k,i,j))
!           Mvv(-1,k,i,j) = sum(Avv(1:9,k,i,j))

!           Muu( 0,k,i,j) = sum(Auu(10:18,k,i,j))
!           Mvv( 0,k,i,j) = sum(Avv(10:18,k,i,j))

!           Muu( 1,k,i,j) = sum(Auu(19:27,k,i,j))
!           Mvv( 1,k,i,j) = sum(Avv(19:27,k,i,j))

           ! WHL: Taking the (0,0) term in each layer does not give singular matrices for
           !       the confined-shelf and circular-shelf problems.
           !      The solution converges for these problems even though the preconditioner 
           !       is not expected to be very good.

           m = indxA(0,0,-1)
           Muu(-1,k,i,j) = Auu(m,k,i,j)
           Mvv(-1,k,i,j) = Avv(m,k,i,j)

           m = indxA(0,0,0)
           Muu( 0,k,i,j) = Auu(m,k,i,j)
           Mvv( 0,k,i,j) = Avv(m,k,i,j)

           m = indxA(0,0,1)
           Muu( 1,k,i,j) = Auu(m,k,i,j)
           Mvv( 1,k,i,j) = Avv(m,k,i,j)
       enddo
       enddo
       enddo

    else   ! no preconditioning

       if (verbose_pcg .and. main_task) then
          print*, 'Using no preconditioner'
       endif

    endif      ! precond

  end subroutine setup_preconditioner_3d

!****************************************************************************

  subroutine setup_preconditioner_diag_2d(nx,         ny,      &
                                          indxA_2d,            &
                                          Auu,        Avv,     &
                                          Adiagu,     Adiagv)    ! entries of diagonal preconditioner

    ! Set up diagonal preconditioning matrices for 2D SSA-type solve

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
         nx, ny                 ! horizontal grid dimensions (for scalars)
                                ! velocity grid has dimensions (nx-1,ny-1)

    integer, dimension(-1:1,-1:1), intent(in) :: &
         indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    real(dp), dimension(9,nx-1,ny-1), intent(in) ::   &
         Auu, Avv               ! two out of the four components of assembled matrix
                                ! 1st dimension = 9 (node and its nearest neighbors in x and y direction)
                                ! other dimensions = (z,x,y) indices
                                !
                                !    Auu  | Auv
                                !    _____|____
                                !    Avu  | Avv
                                !         |

    real(dp), dimension(nx-1,ny-1), intent(out) ::   &
         Adiagu, Adiagv         ! matrix entries for diagonal preconditioning

    integer :: i, j, m

    ! Form diagonal matrix for preconditioning

    if (verbose_pcg .and. main_task) then
       print*, 'Using diagonal matrix for preconditioning'
    endif  ! verbose_pcg

    m = indxA_2d(0,0)
    Adiagu(:,:) = Auu(m,:,:)
    Adiagv(:,:) = Avv(m,:,:)

  end subroutine setup_preconditioner_diag_2d

!****************************************************************************

  subroutine setup_preconditioner_tridiag_2d(ilocal,       jlocal,    &
                                             itest, jtest, rtest,     &
                                             Adiag,                   &
                                             Asubdiag,     Asupdiag,  &
                                             omega,        denom,     &
                                             xuh,          xlh)

    ! Set up tridiagonal preconditioning matrices for 2D SSA-type solve

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer :: &
         ilocal, jlocal         ! size of input/output arrays; number of locally owned vertices in each direction

    integer, intent(in) :: itest, jtest, rtest   ! coordinates of diagnostic point

    real(dp), dimension(ilocal,jlocal), intent(inout) :: &
         Asubdiag, Adiag, Asupdiag   ! matrix entries from Au for tridiagonal preconditioning
                                     ! output Adiag = reciprocal of diagonal entry

    real(dp), dimension(ilocal,jlocal), intent(out) :: &
         omega,                    & ! work arrays for tridiagonal solve
         denom,                    &
         xuh,                      &
         xlh

    ! local variables

    real(dp), dimension(ilocal) :: omega_uh
    real(dp) :: denom_uh

    integer :: i, j

    ! Form tridiagonal matrix for preconditioning

    if (verbose_pcg .and. main_task) then
       print*, 'Using tridiagonal matrix for preconditioning'
    endif  ! verbose_pcg

    !WHL - debug
    if (verbose_tridiag .and. main_task) then
       print*, 'In setup_preconditioner_tridiag_2d: itest, jtest, rtest =', itest, jtest, rtest
    endif

    ! Notes on array dimensions ilocal and jlocal:
    ! The outer loops (over either rows or columns) are from 1 to jlocal.
    ! The inner loops have size(ilocal), i.e., the size of the tridiagonal system in each row or column.
    !
    !       If solving the Auu matrix problem in rows, then ilocal = staggered_ihi - staggered_ilo + 1
    !                                                       jlocal = staggered_jhi - staggered_jlo + 1
    !
    !       If solving the Avv matrix problem in columns, then the i and j indices are reversed to reduce strides:
    !                                                       ilocal = staggered_jhi - staggered_jlo + 1
    !                                                       jlocal = staggered_ihi - staggered_ilo + 1

    ! Modify entries as needed so that the tridiagonal matrix is nonsingular.
    ! For inactive vertices with zero diagonal entries, set diag = 1, subdiag = supdiag = 0
    ! (so the solution is x = rhs)
    where (Adiag == 0.0d0)
       Adiag    = 1.0d0
       Asubdiag = 0.0d0
       Asupdiag = 0.0d0
    endwhere

    ! initialize work arrays
    omega = 0.0d0
    denom = 0.0d0
    xuh = 0.0d0
    xlh = 0.0d0
    omega_uh = 0.0d0

    ! Loop over locally owned rows (for the Auu matrix problem) or columns (for the Avv problem)

    do j = 1, jlocal

       ! Forward elimination
       ! Note: denom -> 1/denom to speed up future computations
       omega(1,j) = Asupdiag(1,j) / Adiag(1,j)
       do i = 2, ilocal
          denom(i,j) = Adiag(i,j) - Asubdiag(i,j)*omega(i-1,j)
          if (denom(i,j) == 0.0d0) then
             call write_log('ERROR: divzero in setup_preconditioner_tridiag_2d', GM_FATAL)
          else
             denom(i,j) = 1.0d0 / denom(i,j)
          endif
          omega(i,j) = Asupdiag(i,j) * denom(i,j)
       enddo

       ! Back substitution
       ! Note: omega_uh and denom_uh are different from omega and denom, which are intent(out)

       xlh(ilocal,j) = -omega(ilocal,j)
       omega_uh(ilocal) = Asubdiag(ilocal,j) / Adiag(ilocal,j)

       do i = ilocal-1, 1, -1
          xlh(i,j) = -omega(i,j)*xlh(i+1,j)
          denom_uh = Adiag(i,j) - Asupdiag(i,j)*omega_uh(i+1)
          if (denom_uh == 0.0d0) then
             call write_log('ERROR: divzero in setup_preconditioner_tridiag_2d', GM_FATAL)
          endif
          omega_uh(i) = Asubdiag(i,j) / denom_uh
       enddo

       ! Forward substitution
       xuh(1,j) = -omega_uh(1)
       do i = 2, ilocal
          xuh(i,j) = -omega_uh(i)*xuh(i-1,j)
       enddo

    enddo   ! j

    ! Take reciprocal of Adiag, to replace division with multiplication below.
    ! Note: Any zero values have been set to 1 above.
    Adiag = 1.0d0/Adiag

  end subroutine setup_preconditioner_tridiag_2d

!****************************************************************************

  subroutine tridiag_solver_serial_2d(nx,           ny,         &
                                      itest, jtest, rtest,      &
                                      ilocal,       jlocal,     &
                                      Adiag_u,      Adiag_v,    &
                                      Asubdiag_u,   Asubdiag_v, &
                                      Asupdiag_u,   Asupdiag_v, &
                                      bu,           bv,         &
                                      xu,           xv)

    use glimmer_utils, only: tridiag

    integer, intent(in) :: &
         nx, ny              ! horizontal grid dimensions (for scalars)

    integer, intent(in) :: itest, jtest, rtest   ! diagnostic only

    integer, intent(in) :: &
         ilocal, jlocal      ! size of vectors passed to subroutine tridiag

    real(dp), dimension(ilocal,jlocal), intent(in), optional :: &
         Asubdiag_u, Adiag_u, Asupdiag_u  ! matrix entries from Auu for tridiagonal preconditioning

    real(dp), dimension(jlocal,ilocal), intent(in), optional :: &
         Asubdiag_v, Adiag_v, Asupdiag_v  ! matrix entries from Avv for tridiagonal preconditioning
                                          ! Note: 1st index is j, 2nd index is i

    !TODO - Change index order for bv and xv?
    real(dp), intent(in), dimension(nx-1,ny-1) :: &
         bu, bv                 ! right-hand side vectors

    real(dp), intent(out), dimension(nx-1,ny-1) :: &
         xu, xv                 ! solution vectors

    ! local variables

!    real(dp), dimension(ilocal) :: &
!         yu                             ! rhs for Au tridiagonal solve
                                         ! ilocal = staggered_ihi - staggered_ilo + 1
!    real(dp), dimension(jlocal) :: &
!         yv                             ! rhs for Av tridiagonal solve
                                         ! jlocal = staggered_jhi - staggered_jlo + 1

    integer :: i, j, ii, jj

    !-------------------------------------------------------------------------------------
    ! Solve a tridiagonal system of the form A*x = b.
    ! The solution x is split into u and v components:
    !
    !   |Au   0 | |xu|   |bu|
    !   |       | |  | = |  |
    !   | 0   Av| |xv|   |bv|
    !
    !   In the u matrix, the subdiag entry for vertex(i,j) refers to vertex(i-1,j),
    !    and the superdiag entry refers to vertex(i+1,j).
    !   There is no connection between the last vertex on row (:,j-1) and the first vertex on row (:,j),
    !    or between the last vertex on row (:,j) and the first vertex on row (:,j+1).
    !
    !   Similarly, in the v matrix, the subdiag entry for vertex(i,j) refers to vertex(i,j-1),
    !    and the superdiag entry refers to vertex(i,j+1).
    !   There is no connection between the last vertex on row (i-1,:) and the first vertex on row (i,:),
    !    or between the last vertex on row (i,:) and the first vertex on row (i+1,:).
    !
    !   So we can solve a distinct tridiagonal problem Au(:,j)*xu(:,j) = bu(:,j) for each row j,
    !    and a distinct problem Av(i,:)*xv(i,:) = bv(i,:) for each column i.
    !-------------------------------------------------------------------------------------

    !WHL - debug
    if (verbose_tridiag .and. this_rank == rtest) then
       print*, 'In tridiag_solver_2d_serial: itest, jtest, rtest =', itest, jtest, rtest
    endif

    !-------------------------------------------------------------------------------------
    ! Solve a tridiagonal problem Au*xu = bu for each row j, from the bottom to the top of the global domain.
    !-------------------------------------------------------------------------------------

    ! Initialize the solution vector
    xu(:,:) = 0.0d0

    ! Loop over locally owned rows
    do j = 1, jlocal
       jj = j + staggered_jlo - 1   ! j = 1 corresponds to jj = staggered_jlo

       !Note: Code commented out since it seems OK not to set yv = 0 anywhere.
       ! Initialize the rhs vector to be passed to subroutine tridiag.
!!       yu(:) = bu(staggered_ilo:staggered_ihi, jj)

       ! Modify entries as needed so that the tridiagonal matrix is nonsingular.
       ! The matrix entries were set above, so it only remains to set yu = 0 here.
!!       where (Adiag_u(:,j) == 1.0d0) yu = 0.0d0

       ! Solve for x in row j
       ! Note: Could speed up by inlining the tridiag calculation and precomputing the 'aa' array.
       !       But probably not worth the effort for a serial problem.
       call tridiag(Asubdiag_u(:,j), Adiag_u(:,j), Asupdiag_u(:,j), &
                    xu(staggered_ilo:staggered_ihi,jj), &
                    bu(staggered_ilo:staggered_ihi,jj))

       !WHL - debug
       if (verbose_tridiag .and. this_rank == rtest .and. jj == itest) then
          print*, ' '
          print*, 'xu tridiag solve, this_rank, j =', this_rank, j
          print*, 'i, subdiag_u, diag_u, supdiag_u, rhs_u, x_u:'
          do i = ilocal, 1, -1
             ii = i + staggered_ilo - 1
             write(6,'(i4, 5f12.3)') &
                  i, Asubdiag_u(i,j), Adiag_u(i,j), Asupdiag_u(i,j), bu(ii,jj), xu(ii,jj)
          enddo
       endif

    enddo  ! j

    !-------------------------------------------------------------------------------------
    ! Solve a tridiagonal problem Av*xv = bv for each column i, from the left to the right of the global domain
    !-------------------------------------------------------------------------------------

    ! Initialize the solution vector
    xv(:,:) = 0.0d0

    ! Loop over locally owned columns
    do i = 1, ilocal

       ii = i + staggered_ilo - 1   ! i = 1 corresponds to ii = staggered_ilo

       !Note: Code commented out since it seems OK not to set yv = 0 anywhere.
       ! Initialize the rhs vector to be passed to subroutine tridiag.
!!       yv(:) = bv(ii, staggered_jlo:staggered_jhi)

       ! Modify entries as needed so that the tridiagonal matrix is nonsingular.
       ! The matrix entries were set above, so it only remains to set yv = 0 here.
!!       where (Adiag_v(:,i) == 1.0d0) yv = 0.0d0

       ! Solve for x in column i
       call tridiag(Asubdiag_v(:,i), Adiag_v(:,i), Asupdiag_v(:,i), &
                    xv(ii,staggered_jlo:staggered_jhi), &
                    bv(ii,staggered_jlo:staggered_jhi))

       !WHL - debug
       if (verbose_tridiag .and. this_rank == rtest .and. ii == itest) then
          print*, ' '
          print*, 'xv tridiag solve, this_rank, i =', this_rank, i
          print*, 'j, subdiag_v, diag_v, supdiag_v, rhs_v, x:'
          do j = jlocal, 1, -1
             jj = j + staggered_jlo - 1
             write(6,'(i4, 5f12.3)') &
                  j, Asubdiag_v(j,i), Adiag_v(j,i), Asupdiag_v(j,i), bv(ii,jj), xv(ii,jj)
          enddo
       endif

    enddo  ! i

    ! Diagnostics
    if (verbose_tridiag .and. this_rank == rtest) then
       ii = itest
       jj = jtest
       i = itest - staggered_ilo + 1
       j = jtest - staggered_jlo + 1
       print*, ' '
       print*, 'Done in tridiag_solver_serial_2d, i, j, r =', ii, jj, this_rank
       print*, 'Au(1:3):', Asubdiag_u(i,j), Adiag_u(i,j), Asupdiag_u(i,j)
       print*, 'bu:', bu(ii,jj)
       print*, 'xu:', xu(ii,jj)
       print*, 'Av(1:3):', Asubdiag_v(j,i), Adiag_v(j,i), Asupdiag_v(j,i)
       print*, 'bv:', bv(ii,jj)
       print*, 'xv:', xv(ii,jj)
    endif

  end subroutine tridiag_solver_serial_2d

!****************************************************************************

  subroutine tridiag_solver_parallel_2d(ilocal, jlocal,           &
                                        tridiag_solver_flag,      &
                                        itest, jtest, rtest,      &
                                        Adiag,                    &
                                        Asubdiag,     Asupdiag,   &
                                        omega,        denom,      &
                                        xuh,          xlh,        &
                                        bu,           xu)


    use glimmer_utils, only: tridiag

    integer, intent(in) :: &
         ilocal, jlocal            ! size of input/output arrays; number of locally owned vertices in each direction

    character(len=3), intent(in) ::  &
         tridiag_solver_flag       ! either 'row' or 'col' depending on whether solving over rows or columns

    integer, intent(in) :: itest, jtest, rtest   ! coordinates of diagnostic point

    real(dp), dimension(ilocal,jlocal), intent(in) :: &
         Asubdiag, Adiag, Asupdiag    ! matrix entries from A for tridiagonal preconditioning
                                      !       Adiag = reciprocal of diagonal entry
                                      !       Asupdiag is not used below; included for diagnostics only

    real(dp), dimension(ilocal,jlocal), intent(in) :: &
         omega,                     & ! work arrays for Au tridiagonal solve
         denom,                     & !
         xuh,                       & ! upper homogenous solution
         xlh,                       & ! lower homogenous solution
         bu                           ! right-hand side vector

    real(dp), dimension(ilocal,jlocal), intent(out) :: &
         xu                           ! solution vector

    ! local variables

    logical :: main_task_rc           ! either main_task_row or main_task_col
    integer :: tasks_rc               ! either tasks_row or tasks_col
    integer :: i, j, n, p, m
    integer :: ibase

    real(dp), dimension(ilocal,jlocal) :: &
         xr,                  & ! particular solution for Au tridiagonal solve
         gamma                  ! work array for Au tridiagonal solve

    real(dp), dimension(8,jlocal) :: &
         outdata                ! data computed locally for global tridiagonal problem

    real(dp), dimension(2,jlocal) :: &
         local_coeffs           ! coefficients for tridiagonal solution, scattered to local tasks

    real(dp), dimension(:,:), allocatable :: &
         global_outdata,      & ! outdata, gathered to main_task
         global_coeffs          ! coefficients for tridiagonal solution, computed on main_task

    real(dp), dimension(:), allocatable :: &
         subdiag, diag, supdiag, rhs, coeffs  ! matrix entries for the triadiagonal problem solved on main_task

    !-------------------------------------------------------------------------------------
    ! Solve a tridiagonal system of the form A*x = b.
    ! The larger tridiagonal problem is split into u and v components:
    !
    !   |Au   0 | |xu|   |bu|
    !   |       | |  | = |  |
    !   | 0   Av| |xv|   |bv|
    !
    ! Depending on the input arguments, this subroutine solves either Au*xu = bu, or Av*xv = bv.
    !
    !   In the u matrix, the subdiag entry for vertex(i,j) refers to vertex(i-1,j),
    !    and the superdiag entry refers to vertex(i+1,j).
    !   There is no connection between the last vertex on row (:,j-1) and the first vertex on row (:,j),
    !    or between the last vertex on row (:,j) and the first vertex on row (:,j+1).
    !
    !   Similarly, in the v matrix, the subdiag entry for vertex(i,j) refers to vertex(i,j-1),
    !    and the superdiag entry refers to vertex(i,j+1).
    !   There is no connection between the last vertex on row (i-1,:) and the first vertex on row (i,:),
    !    or between the last vertex on row (i,:) and the first vertex on row (i+1,:).
    !
    !   So we can solve a distinct tridiagonal problem Au(:,j)*xu(:,j) = bu(:,j) for each row j,
    !    and a distinct problem Av(i,:)*xv(i,:) = bv(i,:) for each column i.
    !
    ! The parallel algorithm is based on this paper:
    !   N. Mattor, T. Williams, and D. W. Hewett, 1995: Algorithm for solving tridiagonal
    !   matrix problems in parallel, Parallel Computing, 21, 1769-1782.
    !-------------------------------------------------------------------------------------

    !WHL - debug
    if (verbose_tridiag .and. main_task) then
       print*, ' '
       print*, 'In tridiag_solver_parallel_2d: itest, jtest, rtest =', itest, jtest, rtest
       if (tridiag_solver_flag == 'row') then
          print*, 'Solving a tridiag problem for each matrix row'
       elseif (tridiag_solver_flag == 'col') then
          print*, 'Solving a tridiag problem for each matrix column'
       else
          call write_log('ERROR: Invalid value for tridiag_solver flag; should be "row" or "col"', GM_FATAL)
       endif
    endif

    ! Notes on array dimensions ilocal and jlocal:
    ! The outer loops (over either rows or columns) are from 1 to jlocal.
    ! The inner loops have size(ilocal), i.e., the size of the tridiagonal system in each row or column.
    !
    ! Note: If solving the Auu matrix problem in rows, then ilocal = staggered_ihi - staggered_ilo + 1
    !                                                       jlocal = staggered_jhi - staggered_jlo + 1
    !
    !       If solving the Avv matrix problem in columns, then the i and j indices are reversed to reduce strides:
    !                                                       ilocal = staggered_jhi - staggered_jlo + 1
    !                                                       jlocal = staggered_ihi - staggered_ilo + 1

    if (tridiag_solver_flag == 'row') then
       tasks_rc = tasks_row
       main_task_rc = main_task_row
    elseif (tridiag_solver_flag == 'col') then
       tasks_rc = tasks_col
       main_task_rc = main_task_col
    endif

    if (verbose_tridiag .and. main_task) then
       print*, 'tasks_rc, main_task_rc =', tasks_rc, main_task_rc
       call flush(6)
    endif

    !-------------------------------------------------------------------------------------
    ! Solve a tridiagonal problem Au*xu = bu for each row, from the bottom to the top of the global domain.
    ! (Alternatively, solve Av*xv = bv for each column, from the left to the right of the global domain.)
    !-------------------------------------------------------------------------------------

    ! Initialize the solution arrays
    xu(:,:) = 0.0d0
    xr(:,:) = 0.0d0

    ! Loop over locally owned rows (or columns)
    ! Note: The input arrays Asubdiag, Adiag and Asupdiag have indices (i,j).
    !       The local arrays omega, gamma, y, xr, xlh and xuh also have indices (i,j).
    do j = 1, jlocal

       ! Forward elimination
       ! Note: denom -> 1/denom in setup subroutine
       !       Adiag -> 1/Adiag in setup subroutine

       gamma(1,j) = bu(1,j) * Adiag(1,j)   ! actually bu/Adiag
       do i = 2, ilocal
          gamma(i,j) = (bu(i,j) - Asubdiag(i,j)*gamma(i-1,j)) * denom(i,j)   ! actually (bu-Asubdiag*gamma)/denom
       enddo

       ! Back substitution
       xr(ilocal,j) = gamma(ilocal,j)
       do i = ilocal-1, 1, -1
          xr(i,j) = gamma(i,j) - omega(i,j)*xr(i+1,j)
       enddo

       ! Write contributions of this task to the outdata array

       outdata(1,j) = -1.0d0
       outdata(2,j) = xuh(1,j)
       outdata(3,j) = xlh(1,j)
       outdata(4,j) = -xr(1,j)
       outdata(5,j) = xuh(ilocal,j)
       outdata(6,j) = xlh(ilocal,j)
       outdata(7,j) = -1.0d0
       outdata(8,j) = -xr(ilocal,j)

       !WHL - debug
       if (verbose_tridiag .and. this_rank == rtest .and. j == jtest) then
          print*, ' '
          print*, 'jtest =', jtest
          print*, 'After forward substitution, i, omega, gamma, bu, xuh, xlh, xr:'
          do i = ilocal, 1, -1
             write(6,'(i4, 6e12.3)') i, omega(i,j), gamma(i,j), bu(i,j), xuh(i,j), xlh(i,j), xr(i,j)
          enddo
          print*, ' '
          print*, 'outdata(1:4)'
          write(6,'(4e12.3)') outdata(1:4,j)
          print*, 'outdata(5:8)'
          write(6,'(4e12.3)') outdata(5:8,j)
       endif

    enddo  ! j

    if (verbose_tridiag .and. main_task) then
       print*, 'tridiag solve on main task'
       call flush(6)
    endif

    call t_startf("pcg_tridiag_main_task")

    if (tasks_rc > 1) then

       ! Use the row- or column-based communicator to gather outdata_u on main_task_row.
       ! The global array is allocated in the subroutine.

       if (tridiag_solver_flag == 'row') then
          call distributed_gather_var_row(outdata, global_outdata)
       elseif (tridiag_solver_flag == 'col') then
          call distributed_gather_var_col(outdata, global_outdata)
       endif

       ! On the main task of each row or column, put outdata in tridiagonal form and solve.

       if (main_task_rc) then

          allocate(global_coeffs(2*tasks_rc,jlocal))
          global_coeffs(:,:) = 0.0d0

          allocate(diag(2*tasks_rc-2))
          allocate(subdiag(2*tasks_rc-2))
          allocate(supdiag(2*tasks_rc-2))
          allocate(rhs(2*tasks_rc-2))
          allocate(coeffs(2*tasks_rc-2))

          do j = 1, jlocal
             do n = 1, tasks_rc

                ibase = 8*(n-1)  ! base index of global_outdata array
                p = 2*n-2        ! matrix row

                if (n==1) then   ! westernmost or southernmost task
                   subdiag(p+1) = global_outdata(ibase+5,j)
                   diag(p+1)    = global_outdata(ibase+6,j)
                   supdiag(p+1) = global_outdata(ibase+7,j)
                   rhs(p+1)     = global_outdata(ibase+8,j)
                elseif (n==tasks_rc) then   ! easternmost or northernmost task
                   subdiag(p) = global_outdata(ibase+1,j)
                   diag(p)    = global_outdata(ibase+2,j)
                   supdiag(p) = global_outdata(ibase+3,j)
                   rhs(p)     = global_outdata(ibase+4,j)
                else
                   subdiag(p)   = global_outdata(ibase+1,j)
                   diag(p)      = global_outdata(ibase+2,j)
                   supdiag(p)   = global_outdata(ibase+3,j)
                   rhs(p)       = global_outdata(ibase+4,j)
                   subdiag(p+1) = global_outdata(ibase+5,j)
                   diag(p+1)    = global_outdata(ibase+6,j)
                   supdiag(p+1) = global_outdata(ibase+7,j)
                   rhs(p+1)     = global_outdata(ibase+8,j)
                endif

             enddo   ! tasks_rc

             if (verbose_tridiag .and. this_rank==rtest .and. j==jtest) then
                print*, 'Fill global_outdata array, this_rank, j =', this_rank, j
                do m = 1, 8*tasks_rc
                   write(6,'(e12.3)',advance='no') global_outdata(m,j)
                enddo
             endif

             ! Check for zeroes along the main diagonal.
             ! Where diag = 0, set diag = 1 and subdiag = supdiag = rhs = 0, giving coeffs = 0.
             where (diag == 0.0d0)
                diag = 1.0d0
                subdiag = 0.0d0
                supdiag = 0.0d0
                rhs = 0.0d0
             endwhere

             call tridiag(subdiag, diag, supdiag, coeffs, rhs)

             global_coeffs(2:2*tasks_rc-1,j) = coeffs(:)

             if (verbose_tridiag .and. j==jtest) then
                print*, ' '
                print*, 'Solved global tridiag problem'
                print*, 'subdiag, diag, supdiag, rhs, coeffs:'
                do m = 1, 2*tasks_rc-2
                   write(6,'(5e12.3)') subdiag(m), diag(m), supdiag(m), rhs(m), coeffs(m)
                enddo
             endif

          enddo   ! j

          deallocate(diag, subdiag, supdiag, rhs, coeffs)

       endif   ! main_task_rc

       ! Scatter the coefficients back to the local tasks.
       ! Each task receives 2 coefficients for each value of j:
       !  local_coeffs_u(1,:) = uh_coeff, and local_coeffs_u(2,:) = lh_coeff.
       ! Note: uh_coeff = 0 on the westernmost task and lh_coeff = 0 on the easternmost task.

       local_coeffs(:,:) = 0.d0

       if (tridiag_solver_flag == 'row') then
          call distributed_scatter_var_row(local_coeffs, global_coeffs)
       elseif (tridiag_solver_flag == 'col') then
          call distributed_scatter_var_col(local_coeffs, global_coeffs)
       endif

       ! Use the coefficients to combine xr, xuh and xlh into the full solution xu.
       ! Note that xu has dimensions (nx-1,ny-1).

       do j = 1, jlocal
          do i = 1, ilocal
             xu(i,j) = xr(i,j) + local_coeffs(1,j)*xuh(i,j) + local_coeffs(2,j)*xlh(i,j)
          enddo
       enddo   ! j

       if (allocated(global_coeffs)) deallocate(global_coeffs)

    else    ! tasks_row = 1

       xu(:,:) = xr(:,:)

    endif   ! tasks_rc > 1

    call t_stopf("pcg_tridiag_main_task")

    !WHL - debug
    if (verbose_tridiag .and. this_rank == rtest) then
       print*, 'ilocal, jlocal =', ilocal, jlocal
       print*, 'size(xu) =', size(xu,1), size(xu,2)
       do j = 1, jlocal
          if (j == jtest) then
             print*, ' '
             print*, 'xu tridiag solve, this_rank, j =', this_rank, j
             print*, 'i, xu:'
             do i = ilocal, 1, -1
!!                write(6,'(i4, e12.3)') i, xu(i,j)
                print*, i, xu(i,j)
             enddo
          endif
       enddo   ! j
    endif   ! verbose_tridiag

  end subroutine tridiag_solver_parallel_2d

!****************************************************************************

  subroutine global_sum_staggered_3d_real8(nx,     ny,         &
                                           nz,     nhalo,      &
                                           global_sum,         &
                                           work1,  work2)

     ! Sum one or two local arrays on the staggered grid, then take the global sum.

     integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz,                 &  ! number of vertical layers at which velocity is computed
       nhalo                  ! number of halo layers (for scalars)

     real(dp), intent(out) :: global_sum   ! global sum

     real(dp), intent(in), dimension(nz,nx-1,ny-1) :: work1            ! local array
     real(dp), intent(in), dimension(nz,nx-1,ny-1), optional :: work2  ! local array

     integer :: i, j, k
     real(dp) :: local_sum

     local_sum = 0.d0

     ! sum over locally owned velocity points
     
     if (present(work2)) then
        do j = staggered_jlo, staggered_jhi
           do i = staggered_ilo, staggered_ihi
              do k = 1, nz
                 local_sum = local_sum + work1(k,i,j) + work2(k,i,j)
              enddo
           enddo
        enddo
     else
        do j = staggered_jlo, staggered_jhi
           do i = staggered_ilo, staggered_ihi
              do k = 1, nz
                 local_sum = local_sum + work1(k,i,j)    
              enddo
           enddo
        enddo
     endif

     ! take the global sum

     global_sum = parallel_reduce_sum(local_sum)

    end subroutine global_sum_staggered_3d_real8

!****************************************************************************

  subroutine global_sum_staggered_3d_real8_nvar(nx,     ny,           &
                                                nz,     nhalo,        &
                                                global_sum,           &
                                                work1,  work2)

     ! Sum one or two local arrays on the staggered grid, then take the global sum.

     integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz,                 &  ! number of vertical layers at which velocity is computed
       nhalo                  ! number of halo layers (for scalars)

     real(dp), intent(out), dimension(:) :: global_sum   ! global sum

     real(dp), intent(in), dimension(nz,nx-1,ny-1,size(global_sum)) :: work1            ! local array
     real(dp), intent(in), dimension(nz,nx-1,ny-1,size(global_sum)), optional :: work2  ! local array

     integer :: i, j, k, n, nvar
     real(dp), dimension(size(global_sum)) :: local_sum

     nvar = size(global_sum)

     local_sum(:) = 0.d0

     do n = 1, nvar

        ! sum over locally owned velocity points

        if (present(work2)) then
           do j = staggered_jlo, staggered_jhi
              do i = staggered_ilo, staggered_ihi
                 do k = 1, nz
                    local_sum(n) = local_sum(n) + work1(k,i,j,n) + work2(k,i,j,n)
                 enddo
              enddo
           enddo
        else
           do j = staggered_jlo, staggered_jhi
              do i = staggered_ilo, staggered_ihi
                 do k = 1, nz
                    local_sum(n) = local_sum(n) + work1(k,i,j,n)    
                 enddo
              enddo
           enddo
        endif

     enddo   ! nvar

     ! take the global sum

     global_sum(:) = parallel_reduce_sum(local_sum(:))

    end subroutine global_sum_staggered_3d_real8_nvar

!****************************************************************************

  subroutine global_sum_staggered_2d_real8(nx,     ny,            &
                                           nhalo,  global_sum,    &
                                           work1,  work2)

     ! Sum one or two local arrays on the staggered grid, then take the global sum.

     integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nhalo                  ! number of halo layers (for scalars)

     real(dp), intent(out) :: global_sum   ! global sum

     real(dp), intent(in), dimension(nx-1,ny-1) :: work1            ! local array
     real(dp), intent(in), dimension(nx-1,ny-1), optional :: work2  ! local array

     integer :: i, j
     real(dp) :: local_sum

     local_sum = 0.d0

     ! sum over locally owned velocity points
     
     if (present(work2)) then
        do j = staggered_jlo, staggered_jhi
           do i = staggered_ilo, staggered_ihi
              local_sum = local_sum + work1(i,j) + work2(i,j)
           enddo
        enddo
     else
        do j = staggered_jlo, staggered_jhi
           do i = staggered_ilo, staggered_ihi
              local_sum = local_sum + work1(i,j) 
           enddo
        enddo
     endif

     ! take the global sum

     global_sum = parallel_reduce_sum(local_sum)

    end subroutine global_sum_staggered_2d_real8

!****************************************************************************

  subroutine global_sum_staggered_2d_real8_nvar(nx,     ny,             &
                                                nhalo,  global_sum,     &
                                                work1,  work2)

     ! Sum one or two local arrays on the staggered grid, then take the global sum.

     integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nhalo                  ! number of halo layers (for scalars)

     real(dp), intent(out), dimension(:) :: global_sum   ! global sum

     real(dp), intent(in), dimension(nx-1,ny-1,size(global_sum)) :: work1            ! local array
     real(dp), intent(in), dimension(nx-1,ny-1,size(global_sum)), optional :: work2  ! local array

     integer :: i, j, n, nvar
     real(dp), dimension(size(global_sum)) :: local_sum

     nvar = size(global_sum)

     local_sum(:) = 0.d0

     do n = 1, nvar

        ! sum over locally owned velocity points

        if (present(work2)) then
           do j = staggered_jlo, staggered_jhi
              do i = staggered_ilo, staggered_ihi
                 local_sum(n) = local_sum(n) + work1(i,j,n) + work2(i,j,n)
              enddo
           enddo
        else
           do j = staggered_jlo, staggered_jhi
              do i = staggered_ilo, staggered_ihi
                 local_sum(n) = local_sum(n) + work1(i,j,n)    
              enddo
           enddo
        endif

     enddo   ! nvar

     ! take the global sum

     global_sum(:) = parallel_reduce_sum(local_sum(:))

    end subroutine global_sum_staggered_2d_real8_nvar

!****************************************************************************

  subroutine matvec_multiply_structured_3d(nx,        ny,            &
                                           nz,        nhalo,         &
                                           indxA,     active_vertex, &
                                           Auu,       Auv,           &
                                           Avu,       Avv,           &
                                           xu,        xv,            &
                                           yu,        yv)

    !---------------------------------------------------------------
    ! Compute the matrix-vector product $y = Ax$.
    !
    ! The A matrices should have complete matrix elements for all
    !  rows corresponding to locally owned vertices.
    ! The terms of x should be correct for all locally owned vertices
    !  and also for all halo vertices adjacent to locally owned vertices.
    ! The resulting y will then be correct for locally owned vertices.
    !
    ! TODO: Are the matvec_multiply routines as efficient as possible?
    !       E.g., could use the active_vertex array to set up indirect addressing and avoid an 'if'
    !       Could replace the three short iA/jA/kA loops with long multadds
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz,                 &  ! number of vertical layers at which velocity is computed
       nhalo                  ! number of halo layers (for scalars)

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27
    
    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       xu, xv             ! current guess for solution


    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::  &
       yu, yv             ! y = Ax

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, k, m
    integer :: iA, jA, kA
    
    ! Initialize the result vector.

    yu(:,:,:) = 0.d0
    yv(:,:,:) = 0.d0

    ! Compute y = Ax

    ! Loop over locally owned vertices
    ! Note: For periodic BC, the southern and western rows of the global domain are halo cells.
    !          For all processors, staggered_ilo = staggered_jlo = staggered_lhalo + 1.
    !       For outflow BC, the southern and western rows of the global domain are locally owned.
    !          For processors owning these rows, staggered_ilo = staggered_jlo = staggered_lhalo.

    do j = staggered_jlo, staggered_jhi
    do i = staggered_ilo, staggered_ihi

       if (active_vertex(i,j)) then

          do k = 1, nz

             do kA = -1,1
             do jA = -1,1
             do iA = -1,1

                if ( (k+kA >= 1 .and. k+kA <= nz)         &
                                .and.                     &
                     (i+iA >= 1 .and. i+iA <= nx-1)       &
                                .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   m = indxA(iA,jA,kA)

                   yu(k,i,j) = yu(k,i,j)   &
                             + Auu(m,k,i,j)*xu(k+kA,i+iA,j+jA)  &
                             + Auv(m,k,i,j)*xv(k+kA,i+iA,j+jA)

                   yv(k,i,j) = yv(k,i,j)   &
                             + Avu(m,k,i,j)*xu(k+kA,i+iA,j+jA)  &
                             + Avv(m,k,i,j)*xv(k+kA,i+iA,j+jA)
 
                endif   ! k+kA, i+iA, j+jA in bounds

             enddo   ! kA
             enddo   ! iA
             enddo   ! jA

          enddo   ! k

       endif   ! active_vertex

    enddo   ! i
    enddo   ! j

  end subroutine matvec_multiply_structured_3d
 
!****************************************************************************

  subroutine matvec_multiply_structured_2d(nx,        ny,            &
                                           nhalo,                    &
                                           indxA_2d,  active_vertex, &
                                           Auu,       Auv,           &
                                           Avu,       Avv,           &
                                           xu,        xv,            &
                                           yu,        yv)

    !---------------------------------------------------------------
    ! Compute the matrix-vector product $y = Ax$.
    !
    ! This subroutine is similar to subroutine matrix_vector_structured,
    ! but modified to solve for x and y at a single horizontal level, 
    ! as in the shallow-shelf approximation.  
    !
    ! The A matrices should have complete matrix elements for all
    !  rows corresponding to locally owned vertices.
    ! The terms of x should be correct for all locally owned vertices
    !  and also for all halo vertices adjacent to locally owned vertices.
    ! The resulting y will then be correct for locally owned vertices.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nhalo                  ! number of halo layers (for scalars)

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9
    
    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(9,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xu, xv             ! current guess for solution


    real(dp), dimension(nx-1,ny-1), intent(out) ::  &
       yu, yv             ! y = Ax

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, m
    integer :: iA, jA
    
    ! Initialize the result vector.

    yu(:,:) = 0.d0
    yv(:,:) = 0.d0

    ! Compute y = Ax
    ! Loop over locally owned vertices
    ! Note: For periodic BC, the southern and western rows of the global domain are halo cells.
    !          For all processors, staggered_ilo = staggered_jlo = staggered_lhalo + 1.
    !       For outflow BC, the southern and western rows of the global domain are locally owned.
    !          For processors owning these rows, staggered_ilo = staggered_jlo = staggered_lhalo.

    do j = staggered_jlo, staggered_jhi
    do i = staggered_ilo, staggered_ihi

       if (active_vertex(i,j)) then

          do jA = -1,1
             do iA = -1,1

                if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                                .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   m = indxA_2d(iA,jA)

                   yu(i,j) = yu(i,j)   &
                           + Auu(m,i,j)*xu(i+iA,j+jA)  &
                           + Auv(m,i,j)*xv(i+iA,j+jA)

                   yv(i,j) = yv(i,j)   &
                           + Avu(m,i,j)*xu(i+iA,j+jA)  &
                           + Avv(m,i,j)*xv(i+iA,j+jA)
 
                endif   ! i+iA, j+jA in bounds

             enddo   ! iA
          enddo      ! jA

       endif   ! active_vertex

    enddo   ! i
    enddo   ! j

  end subroutine matvec_multiply_structured_2d
 
!****************************************************************************

  subroutine easy_sia_solver(nx,   ny,   nz,    &
                             active_vertex,     & 
                             A,    b,    x)

    !---------------------------------------------------------------
    ! Solve the problem Ax = b where A is a local shallow-ice matrix,
    !  with coupling in the vertical but not the horizontal.
    ! We simply solve a tridiagonal matrix for each column.
    !---------------------------------------------------------------
   
    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz                     ! number of vertical levels
           
    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(in) ::   &
       A                      ! matrix with vertical coupling only
                              ! 1st dimension = node and its upper and lower neighbors

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       b                      ! right-hand side

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       x                      ! solution

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    real(dp), dimension(nz) ::  &
       sbdiag,         &      ! subdiagonal matrix entries
       diag,           &      ! diagonal matrix entries
       spdiag,         &      ! superdiagonal matrix entries
       rhs,            &      ! right-hand side
       soln                   ! tridiagonal solution

    integer :: i, j, k

    do j = 1, ny-1
    do i = 1, nx-1

       if (active_vertex(i,j)) then

          ! initialize rhs and solution 

          rhs(:)  = b(:,i,j)
          soln(:) = x(:,i,j)

          ! top layer

          k = 1
          sbdiag(k) = 0.d0
          diag(k)   = A(0,k,i,j)
          spdiag(k) = A(1,k,i,j)

          ! intermediate layers

          do k = 2, nz-1
             sbdiag(k) = A(-1,k,i,j)
             diag(k)   = A( 0,k,i,j)
             spdiag(k) = A( 1,k,i,j)
          enddo

          ! bottom layer

          k = nz 
          sbdiag(k) = A(-1,k,i,j)
          diag(k)   = A( 0,k,i,j)
          spdiag(k) = 0.d0
         
          ! solve

          call tridiag_solver(nz,    sbdiag,   &
                              diag,  spdiag,   &
                              rhs,   soln)

          x(:,i,j) = soln(:)

       endif  ! active_vertex

    enddo     ! i
    enddo     ! j

  end subroutine easy_sia_solver
  
!****************************************************************************

  subroutine tridiag_solver(order,    sbdiag,   &
                            diag,     spdiag,   &
                            rhs,      soln)

    !---------------------------------------------------------------
    ! Solve a 1D tridiagonal matrix problem.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       order            ! matrix dimension

    real(dp), dimension(order), intent(in) :: &
       sbdiag,        & ! sub-diagonal matrix elements
       diag,          & ! diagonal matrix elements
       spdiag,        & ! super-diagonal matrix elements
       rhs              ! right-hand side

    real(dp), dimension(order), intent(inout) :: &
       soln             ! solution vector

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: &
       k               ! row counter

    real(dp) :: &
       beta            ! temporary matrix variable

    real(dp), dimension(order) :: &
       gamma           ! temporary matrix variable

    ! Solve

    beta = diag(1)
    soln(1) = rhs(1) / beta

    do k = 2, order
       gamma(k) = spdiag(k-1) / beta
       beta = diag(k) - sbdiag(k)*gamma(k)
       soln(k) = (rhs(k) - sbdiag(k)*soln(k-1)) / beta
    enddo

    do k = order-1, 1, -1
       soln(k) = soln(k) - gamma(k+1)*soln(k+1)
    enddo

  end subroutine tridiag_solver

!****************************************************************************

end module glissade_velo_higher_pcg

!****************************************************************************
