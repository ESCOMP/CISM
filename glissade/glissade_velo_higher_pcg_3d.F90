!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_higher_pcg_3d.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains the 3D solvers called from glissade_velo_higher.F90 to solve
! the problem Ax = b using the preconditioned conjugate gradient method, together with
! their preconditioners and matrix-vector multiply. The 2D solvers are in
! glissade_velo_higher_pcg_2d.
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

  module glissade_velo_higher_pcg_3d

    use glimmer_global, only: dp
    use glimmer_paramets, only: iulog
    use glissade_types   ! for preconditioning options
    use glimmer_log
    use cism_profile, only: t_startf, t_stopf
    use cism_parallel, only: this_rank, main_task, &
         parallel_type, staggered_parallel_halo, parallel_global_sum_stagger
    use glissade_velo_higher_pcg_2d, only: &
         setup_preconditioner_tridiag_local_2d, setup_preconditioner_tridiag_global_2d, &
         tridiag_solver_local_2d, tridiag_solver_global_2d

    implicit none
    
    private
    public :: pcg_solver_standard_3d, pcg_solver_chrongear_3d, &
              matvec_multiply_structured_3d

    logical :: verbose_pcg = .false.
    logical :: verbose_tridiag = .false.

  contains

!****************************************************************************

  subroutine pcg_solver_standard_3d(nx,        ny,            &
                                    nz,        parallel,      &
                                    indxA_2d,  indxA_3d,      &
                                    active_vertex,            &
                                    Auu,       Auv,           &
                                    Avu,       Avv,           &
                                    bu,        bv,            &
                                    xu,        xv,            &
                                    precond,   linear_solve_ncheck,  &
                                    tolerance, maxiters,      &
                                    err,       niters,        &
                                    itest, jtest, rtest)

    !---------------------------------------------------------------
    !  This subroutine uses a standard preconditioned conjugate-gradient algorithm
    !  to solve the equation $Ax=b$.
    !  Convergence is checked every {\em linear_solve_ncheck} steps.
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
       nz                       ! number of vertical levels where velocity is computed

    type(parallel_type), intent(in) :: &
       parallel               ! info for parallel communication

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA_3d               ! maps relative (x,y,z) coordinates to an index between 1 and 27

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

    integer, intent(in)  :: &
       linear_solve_ncheck          ! number of iterations between convergence checks in the linear solver

    integer, intent(in) ::    &
       maxiters                     ! max number of linear iterations before quitting

    real(dp), intent(in) ::   &
       tolerance                    ! tolerance for linear solver

    real(dp), intent(out) ::  &
       err                          ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                       ! iterations needed to solution

    integer, intent(in) :: &
       itest, jtest, rtest          ! point for debugging diagnostics

    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer :: i, j, k        ! grid indices
    integer :: iA, jA, kA     ! grid offsets ranging from -1 to 1
    integer :: m              ! matrix element index
    integer :: ilocal, jlocal ! number of locally owned vertices in each direction
    integer :: iter           ! iteration counter

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

    ! arrays for tridiagonal preconditioning
    ! Note: 2D diagonal entries are Adiag_u and Adiag_v; distinct from 3D Adiagu and Adiagv above
    real(dp), dimension(:,:), allocatable :: &
         Asubdiag_u, Adiag_u, Asupdiag_u,  &  ! matrix entries from Auu for tridiagonal preconditioning
         Asubdiag_v, Adiag_v, Asupdiag_v      ! matrix entries from Avv for tridiagonal preconditioning

    real(dp), dimension(:,:), allocatable :: &
         omega_u, omega_v,         & ! work arrays for tridiagonal solve
	 denom_u, denom_v,         &
         xuh_u,   xuh_v,           &
         xlh_u,   xlh_v

    ! Note: These two matrices are global in the EW and NS dimensions, respectively.
    !       Each holds 8 pieces of information for each task on each row or column.
    !       Since only 2 of these 8 pieces of information change from one iteration to the next,
    !        it is more efficient to gather the remaining information once and pass the arrays
    !        with intent(inout), than to declare the arrays in subroutine tridiag_solver_global_2d
    !        and gather all the information every time the subroutine is called.
    ! TODO: Revisit this. Is the efficiency gain large enough to justify the extra complexity?

    real(dp), dimension(:,:), allocatable :: &
         gather_data_row,   &   ! arrays for gathering data from every task on a row or column
         gather_data_col

    integer :: &
         tasks_row,         &   ! number of tasks per row and column for tridiagonal solve
         tasks_col

    logical :: first_time       ! true on the first subroutine call (iter = 1), false thereafter

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    !WHL - debug
    integer :: iu_max, ju_max, iv_max, jv_max
    real(dp) :: ru_max, rv_max

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    tasks_row = parallel%tasks_row
    tasks_col = parallel%tasks_col

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Using native PCG solver (standard)'
       write(iulog,*) 'tolerance, maxiters, precond =', tolerance, maxiters, precond
    endif

    ! Compute array sizes for locally owned vertices
    ilocal = staggered_ihi - staggered_ilo + 1
    jlocal = staggered_jhi - staggered_jlo + 1

    ! Set up matrices for preconditioning

    call t_startf("pcg_precond_init")

    if (precond == PRECONDITIONER_NONE) then    ! no preconditioner

       if (verbose_pcg .and. this_rank == rtest) then
          write(iulog,*) 'Using no preconditioner'
       endif

    elseif (precond == PRECONDITIONER_DIAGONAL) then

       call setup_preconditioner_diag_3d(nx,       ny,         &
                                         nz,       indxA_3d,   &
                                         Auu,      Avv,        &
                                         Adiagu,   Adiagv)

       !WHL - debug
       if (verbose_pcg .and. this_rank == rtest) then
          i = itest
          j = jtest
          write(iulog,*) 'i, j, r =', i, j, this_rank
          write(iulog,*) 'Auu diag =', Adiagu(:,i,j)
          write(iulog,*) 'Avv diag =', Adiagv(:,i,j)
       endif

    elseif (precond == PRECONDITIONER_SIA) then

       call setup_preconditioner_sia_3d(nx,        ny,       &
                                        nz,        indxA_3d, &
                                        Auu,       Avv,      &
                                        Muu,       Mvv)

       if (verbose_pcg .and. this_rank == rtest) then
          j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'i, k, Muu_sia, Mvv_sia:'
          do i = staggered_ihi, staggered_ilo, -1
             write(iulog,*) ' '
             do k = 1, nz
                write(iulog,'(2i4, 6e13.5)') i, k, Muu(:,k,i,j), Mvv(:,k,i,j)
             enddo
          enddo  ! i
       endif

    elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then

       ! Allocate tridiagonal preconditioning matrices
       allocate(Adiag_u   (nx-1,ny-1))
       allocate(Asubdiag_u(nx-1,ny-1))
       allocate(Asupdiag_u(nx-1,ny-1))
       allocate(omega_u   (nx-1,ny-1))
       allocate(denom_u   (nx-1,ny-1))

       allocate(Adiag_v   (nx-1,ny-1))
       allocate(Asubdiag_v(nx-1,ny-1))
       allocate(Asupdiag_v(nx-1,ny-1))
       allocate(omega_v   (nx-1,ny-1))
       allocate(denom_v   (nx-1,ny-1))

       ! Compute arrays for tridiagonal preconditioning

       call setup_preconditioner_tridiag_local_3d(&
            nx,            ny,            &
            nz,            parallel,      &
            active_vertex,                &
            indxA_2d,      indxA_3d,      &
            itest, jtest,  rtest,         &
            Auu,           Avv,           &
            Muu,           Mvv,           &
            Adiag_u,       Adiag_v,       &
            Asubdiag_u,    Asubdiag_v,    &
            Asupdiag_u,    Asupdiag_v,    &
            omega_u,       omega_v,       &
            denom_u,       denom_v)

    elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then

       ! Allocate tridiagonal preconditioning matrices
       ! Note: (i,j) indices are switched for the A_v matrices to reduce striding.
       allocate(Adiag_u   (ilocal,jlocal))
       allocate(Asubdiag_u(ilocal,jlocal))
       allocate(Asupdiag_u(ilocal,jlocal))
       allocate(omega_u(ilocal,jlocal))
       allocate(denom_u(ilocal,jlocal))
       allocate(xuh_u(ilocal,jlocal))
       allocate(xlh_u(ilocal,jlocal))

       allocate(Adiag_v   (jlocal,ilocal))
       allocate(Asubdiag_v(jlocal,ilocal))
       allocate(Asupdiag_v(jlocal,ilocal))
       allocate(omega_v(jlocal,ilocal))
       allocate(denom_v(jlocal,ilocal))
       allocate(xuh_v(jlocal,ilocal))
       allocate(xlh_v(jlocal,ilocal))

       ! These two matrices are for gathering data from all tasks on a given row or column.
       allocate(gather_data_row(8*tasks_row,jlocal))
       allocate(gather_data_col(8*tasks_col,ilocal))
       gather_data_row = 0.0d0
       gather_data_col = 0.0d0

       call setup_preconditioner_tridiag_global_3d(&
            nx,            ny,            &
            nz,            parallel,      &
            active_vertex,                &
            indxA_2d,      indxA_3d,      &
            ilocal,        jlocal,        &
            itest, jtest,  rtest,         &
            Auu,           Avv,           &
            Muu,           Mvv,           &
            Adiag_u,       Adiag_v,       &
            Asubdiag_u,    Asubdiag_v,    &
            Asupdiag_u,    Asupdiag_v,    &
            omega_u,       omega_v,       &
            denom_u,       denom_v,       &
            xuh_u,         xuh_v,         &
            xlh_u,         xlh_v)

    endif   ! precond

    call t_stopf("pcg_precond_init")

    ! Compute initial residual and initialize the direction vector d
    ! Note: The matrix A must be complete for all rows corresponding to locally 
    !        owned vertices, and x must have the correct values in
    !        halo vertices bordering the locally owned vertices.
    !       Then y = Ax will be correct for locally owned vertices.

    ! Halo update for x (initial guess for velocity solution)

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(xu, parallel)
    call staggered_parallel_halo(xv, parallel)
    call t_stopf("pcg_halo_init")

    ! Compute A*x (use z as a temp vector for A*x)

    call t_startf("pcg_matmult_init")
    call matvec_multiply_structured_3d(nx,            ny,            &
                                       nz,            parallel,      &
                                       indxA_3d,      active_vertex, &
                                       Auu,           Auv,           &
                                       Avu,           Avv,           &
                                       xu,            xv,            &
                                       zu,            zv)
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
    L2_rhs = parallel_global_sum_stagger(work0u, parallel, arr2=work0v)
    call t_stopf("pcg_glbsum_init")

    ! take square root

    L2_rhs = sqrt(L2_rhs)       ! L2 norm of RHS

    ! iterate to solution

    iter_loop: do iter = 1, maxiters

       call t_startf("pcg_precond")

       ! Compute PC(r) = solution z of Mz = r

       if (precond == PRECONDITIONER_NONE) then      ! no preconditioning

           zu(:,:,:) = ru(:,:,:)         ! PC(r) = r     
           zv(:,:,:) = rv(:,:,:)         ! PC(r) = r    

       elseif (precond == PRECONDITIONER_DIAGONAL ) then  ! diagonal preconditioning

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

       elseif (precond == PRECONDITIONER_SIA) then   ! local vertical shallow-ice solver for preconditioning

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Muu,  ru,   zu)      ! solve Muu*zu = ru for zu 

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Mvv,  rv,   zv)      ! solve Mvv*zv = rv for zv

       elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then

          ! Use a local tridiagonal solver to find an approximate solution of A*z = r

          call tridiag_solver_local_3d(&
               nx,            ny,            &
               nz,            parallel,      &
               active_vertex,                &
               itest, jtest,  rtest,         &
               Adiag_u,       Adiag_v,       &  ! entries of 2D preconditioning matrix
               Asubdiag_u,    Asubdiag_v,    &
               Asupdiag_u,    Asupdiag_v,    &
               omega_u,       omega_v,       &
               denom_u,       denom_v,       &
               Muu,           Mvv,           &  ! entries of SIA matrix
               ru,            rv,            &  ! 3D residual
               zu,            zv)               ! approximate solution of Az = r

       elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then

          ! Use a global tridiagonal solver to find an approximate solution of A*z = r

          if (iter == 1) then
             first_time = .true.
          else
             first_time = .false.
          endif

          call tridiag_solver_global_3d(&
               nx,              ny,              &
               nz,              parallel,        &
               active_vertex,                    &
               ilocal,          jlocal,          &
               tasks_row,       tasks_col,       &
               itest,  jtest,   rtest,           &
               Adiag_u,         Adiag_v,         &  ! entries of 2D preconditioning matrix
               Asubdiag_u,      Asubdiag_v,      &
               Asupdiag_u,      Asupdiag_v,      &
               omega_u,         omega_v,         &
               denom_u,         denom_v,         &
               xuh_u,           xuh_v,           &
               xlh_u,           xlh_v,           &
               Muu,             Mvv,             &  ! entries of SIA matrix
               gather_data_row, gather_data_col, &
               first_time,                       &
               ru,              rv,              &  ! 3D residual
               zu,              zv)                 ! approximate solution of Az = r

       endif    ! precond

       call t_stopf("pcg_precond")

       ! Compute the dot product eta1 = (r, PC(r))

       call t_startf("pcg_dotprod")
       work0u(:,:,:) = ru(:,:,:)*zu(:,:,:)    ! terms of dot product (r, PC(r))
       work0v(:,:,:) = rv(:,:,:)*zv(:,:,:)    
       call t_stopf("pcg_dotprod")

       call t_startf("pcg_glbsum_iter")
       eta1 = parallel_global_sum_stagger(work0u, parallel, arr2=work0v)
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
       call staggered_parallel_halo(du, parallel)
       call staggered_parallel_halo(dv, parallel)
       call t_stopf("pcg_halo_iter")
  
       ! Compute q = A*d
       ! This is the one matvec multiply required for each iteration

       call t_startf("pcg_matmult_iter")
       call matvec_multiply_structured_3d(nx,            ny,            &
                                          nz,            parallel,      &
                                          indxA_3d,      active_vertex, &
                                          Auu,           Auv,           &
                                          Avu,           Avv,           &
                                          du,            dv,            &
                                          qu,            qv)
       call t_stopf("pcg_matmult_iter")

       ! Copy old eta1 = (r, PC(r)) to eta0

       eta0 = eta1               ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i)) 

       ! Compute the dot product eta2 = (d, A*d)

       call t_startf("pcg_dotprod")
       work0u(:,:,:) = du(:,:,:) * qu(:,:,:)       ! terms of dot product (d, Ad)
       work0v(:,:,:) = dv(:,:,:) * qv(:,:,:)
       call t_stopf("pcg_dotprod")

       call t_startf("pcg_glbsum_iter")
       eta2 = parallel_global_sum_stagger(work0u, parallel, arr2=work0v)
       call t_stopf("pcg_glbsum_iter")

       ! Compute alpha
                              !          (r, PC(r))
       alpha = eta1/eta2      ! alpha = ----------
                              !          (d, A*d)
       
       !WHL - If eta2 = 0 (e.g., because all matrix entries are zero), then alpha = NaN
 
       if (alpha /= alpha) then  ! alpha is NaN
!!          write(iulog,*) 'eta1, eta2, alpha:', eta1, eta2, alpha
          call write_log('PCG solver has failed, alpha = NaN', GM_FATAL)
       endif

       ! Compute the new solution and residual

       call t_startf("pcg_vecupdate")
       xu(:,:,:) = xu(:,:,:) + alpha * du(:,:,:)    ! new solution, x_(i+1) = x_i + alpha*d
       xv(:,:,:) = xv(:,:,:) + alpha * dv(:,:,:)

       ru(:,:,:) = ru(:,:,:) - alpha * qu(:,:,:)    ! new residual, r_(i+1) = r_i - alpha*(Ad)
       rv(:,:,:) = rv(:,:,:) - alpha * qv(:,:,:)
       call t_stopf("pcg_vecupdate")

       if (verbose_pcg .and. this_rank == rtest) then
          j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'iter =', iter
          write(iulog,*) 'i, k, xu, xv, ru, rv:'
          do i = itest-3, itest+3
             write(iulog,*) ' '
             do k = 1, nz
                write(iulog,'(2i4, 4f16.10)') i, k, xu(k,i,j), xv(k,i,j), ru(k,i,j), rv(k,i,j)
             enddo
          enddo  ! i
       endif

       ! Check for convergence every linear_solve_ncheck iterations.
       ! Also check at iter = 5, to reduce iterations when the nonlinear solver is close to convergence.
       ! TODO: Check at iter = linear_solve_ncheck/2 instead of 5?  This would be answer-changing.
       !
       ! For convergence check, use r = b - Ax

       if (mod(iter, linear_solve_ncheck) == 0 .or. iter == 5) then

          if (verbose_pcg .and. main_task) then
             write(iulog,*) ' '
             write(iulog,*) 'Check convergence, iter =', iter
          endif

          ! Halo update for x

          call t_startf("pcg_halo_resid")
          call staggered_parallel_halo(xu, parallel)
          call staggered_parallel_halo(xv, parallel)
          call t_stopf("pcg_halo_resid")

          ! Compute A*x (use z as a temp vector for A*x)
           
          call t_startf("pcg_matmult_resid")
          call matvec_multiply_structured_3d(nx,            ny,            &
                                             nz,            parallel,      &
                                             indxA_3d,      active_vertex, &
                                             Auu,           Auv,           &
                                             Avu,           Avv,           &
                                             xu,            xv,            &
                                             zu,            zv)
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
          L2_resid = parallel_global_sum_stagger(work0u, parallel, arr2=work0v)
          call t_stopf("pcg_glbsum_resid")

          ! take square root
          L2_resid = sqrt(L2_resid)       ! L2 norm of residual

          ! compute normalized error
          err = L2_resid/L2_rhs

          if (verbose_pcg .and. main_task) then
             write(iulog,*) 'iter, L2_resid, L2_rhs, error =', iter, L2_resid, L2_rhs, err
          endif

          !WHL - debug
          if (verbose_pcg .and. this_rank == rtest) then
             ru_max = 0.d0
             rv_max = 0.d0
             iu_max = 0
             ju_max = 0
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (sum(abs(ru(:,i,j))) > ru_max) then
                      ru_max = sum(abs(ru(:,i,j)))
                      iu_max = i
                      ju_max = j
                   endif
                   if (sum(abs(rv(:,i,j))) > rv_max) then
                      rv_max = sum(abs(rv(:,i,j)))
                      iv_max = i
                      jv_max = j
                   endif
                enddo
             enddo
             write(iulog,*) 'r, i, j, ru_max:', this_rank, iu_max, ju_max, ru_max
             write(iulog,*) 'r, i, j, rv_max:', this_rank, iv_max, jv_max, rv_max
          endif

          ! If converged, then exit the loop.
          ! Note: Without good preconditioning, convergence can be slow,
          !       but the solution after maxiters_chrongear might be good enough.

          if (err < tolerance) then
             niters = iter
             if (verbose_pcg .and. main_task) then
                write(iulog,*) 'Glissade PCG solver has converged, iter =', niters
                write(iulog,*) ' '
             endif
             exit iter_loop
          elseif (iter == maxiters) then
             if (verbose_pcg .and. main_task) then
                write(iulog,*) 'Glissade PCG solver did not converge'
                write(iulog,*) 'iter, err, tolerance:', iter, err, tolerance
                write(iulog,*) ' '
             endif
          endif

       endif    ! linear_solve_ncheck

    enddo iter_loop

    ! Clean up
    if (allocated(Adiag_u))    deallocate(Adiag_u, Adiag_v)
    if (allocated(Asubdiag_u)) deallocate(Asubdiag_u, Asubdiag_v)
    if (allocated(Asupdiag_u)) deallocate(Asupdiag_u, Asupdiag_v)
    if (allocated(omega_u))    deallocate(omega_u, omega_v)
    if (allocated(denom_u))    deallocate(denom_u, denom_v)
    if (allocated(xuh_u))      deallocate(xuh_u, xuh_v)
    if (allocated(xlh_u))      deallocate(xlh_u, xlh_v)
    if (allocated(gather_data_row)) deallocate(gather_data_row)
    if (allocated(gather_data_col)) deallocate(gather_data_col)

  end subroutine pcg_solver_standard_3d

!****************************************************************************
    
  subroutine pcg_solver_chrongear_3d(nx,        ny,            &
                                     nz,        parallel,      &
                                     indxA_2d,  indxA_3d,      &
                                     active_vertex,            &
                                     Auu,       Auv,           &
                                     Avu,       Avv,           &
                                     bu,        bv,            &
                                     xu,        xv,            &
                                     precond,   linear_solve_ncheck,  &
                                     tolerance, maxiters,      &
                                     err,       niters,        &
                                     itest, jtest, rtest)

    !---------------------------------------------------------------
    !  This subroutine uses a Chronopoulos-Gear preconditioned conjugate-gradient
    !  algorithm to solve the equation $Ax=b$.
    !
    !  It is based on the Chronopoulos-Gear PCG solver in the POP ocean model 
    !  (author Frank Bryan, NCAR). It is a rearranged conjugate gradient solver 
    !  that reduces the number of global reductions per iteration from two to one 
    !  (not counting the convergence check).  Convergence is checked every 
    !  {\em linear_solve_ncheck} steps.
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
       nx, ny,                       &  ! horizontal grid dimensions (for scalars)
                                        ! velocity grid has dimensions (nx-1,ny-1)
       nz                               ! number of vertical levels where velocity is computed

    type(parallel_type), intent(in) :: &
       parallel               ! info for parallel communication

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA_3d               ! maps relative (x,y,z) coordinates to an index between 1 and 27

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

    integer, intent(in)  :: &
       linear_solve_ncheck         ! number of iterations between convergence checks in the linear solver

    integer, intent(in) ::    &
       maxiters                    ! max number of linear iterations before quitting

    real(dp), intent(in) ::   &
       tolerance                   ! tolerance for linear solver

    real(dp), intent(out) ::  &
       err                         ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                      ! iterations needed to solution

    integer, intent(in) :: &
       itest, jtest, rtest         ! point for debugging diagnostics

    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer :: i, j, k, m     ! grid indices
    integer :: ii, jj
    integer :: ilocal, jlocal ! number of locally owned vertices in each direction
    integer :: iter           ! iteration counter
    integer :: maxiters_chrongear  ! max number of linear iterations before quitting

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

    ! arrays for tridiagonal preconditioning
    ! Note: 2D diagonal entries are Adiag_u and Adiag_v; distinct from 3D Adiagu and Adiagv above

    real(dp), dimension(:,:), allocatable :: &
         Asubdiag_u, Adiag_u, Asupdiag_u,  &  ! matrix entries from Auu for tridiagonal preconditioning
         Asubdiag_v, Adiag_v, Asupdiag_v      ! matrix entries from Avv for tridiagonal preconditioning

    real(dp), dimension(:,:), allocatable :: &
         omega_u, omega_v,         & ! work arrays for tridiagonal solve
         denom_u, denom_v,         &
         xuh_u,   xuh_v,           &
         xlh_u,   xlh_v

    ! Note: These two matrices are global in the EW and NS dimensions, respectively.
    !       Each holds 8 pieces of information for each task on each row or column.
    !       Since only 2 of these 8 pieces of information change from one iteration to the next,
    !        it is more efficient to gather the remaining information once and pass the arrays
    !        with intent(inout), than to declare the arrays in subroutine tridiag_solver_global_2d
    !        and gather all the information every time the subroutine is called.
    ! TODO: Revisit this. Is the efficiency gain large enough to justify the extra complexity?

    real(dp), dimension(:,:), allocatable :: &
         gather_data_row,   &   ! arrays for gathering data from every task on a row or column
         gather_data_col

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    integer :: &
         tasks_row,         &   ! number of tasks per row and column for tridiagonal solve
         tasks_col

    !WHL - debug
    integer :: iu_max, ju_max, iv_max, jv_max
    real(dp) :: ru_max, rv_max

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    tasks_row = parallel%tasks_row
    tasks_col = parallel%tasks_col

    ! Note: maxiters_tridiag commented out here, because the BP tridiagonal solver
    !       tends not to converge as well as the 2D version.
    ! TODO: Make maxiters a config option.

    ! Set the maximum number of linear iterations.
    ! Typically allow up to 200 iterations with diagonal preconditioning, but only 100
    !  with tridiagonal, which usually converges faster.

    !TODO - Test whether maxiters_tridiag (currently = 100) is sufficient for convergence with 3D solver
!!    if (precond == PRECONDITIONER_TRIDIAG_LOCAL .or. precond == PRECONDITIONER_TRIDIAG_GLOBAL) then
!!       maxiters_chrongear = maxiters_tridiag
!!    else
       maxiters_chrongear = maxiters
!!    endif

    if (verbose_pcg .and. this_rank == rtest) then
       write(iulog,*) 'Using native PCG solver (Chronopoulos-Gear)'
       write(iulog,*) 'tolerance, maxiters, precond =', tolerance, maxiters_chrongear, precond
    endif

    ! Compute array sizes for locally owned vertices
    ilocal = staggered_ihi - staggered_ilo + 1
    jlocal = staggered_jhi - staggered_jlo + 1

    !---- Set up matrices for preconditioning

    call t_startf("pcg_precond_init")

    if (precond == PRECONDITIONER_NONE) then    ! no preconditioner

       if (verbose_pcg .and. this_rank == rtest) then
          write(iulog,*) 'Using no preconditioner'
       endif

    elseif (precond == PRECONDITIONER_DIAGONAL) then

       call setup_preconditioner_diag_3d(nx,       ny,         &
                                         nz,       indxA_3d,   &
                                         Auu,      Avv,        &
                                         Adiagu,   Adiagv)

       !WHL - debug
       if (verbose_pcg .and. this_rank == rtest) then
          i = itest
          j = jtest
          write(iulog,*) 'i, j, r =', i, j, this_rank
          write(iulog,*) 'Auu diag =', Adiagu(:,i,j)
          write(iulog,*) 'Avv diag =', Adiagv(:,i,j)
       endif

    elseif (precond == PRECONDITIONER_SIA) then

       call setup_preconditioner_sia_3d(nx,        ny,       &
                                        nz,        indxA_3d, &
                                        Auu,       Avv,      &
                                        Muu,       Mvv)

       if (verbose_pcg .and. this_rank == rtest) then
          j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'i, k, Muu_sia, Mvv_sia:'
          do i = staggered_ihi, staggered_ilo, -1
             write(iulog,*) ' '
             do k = 1, nz
                write(iulog,'(2i4, 6e13.5)') i, k, Muu(:,k,i,j), Mvv(:,k,i,j)
             enddo
          enddo  ! i
       endif

    elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then

       ! Allocate tridiagonal preconditioning matrices
       allocate(Adiag_u   (nx-1,ny-1))
       allocate(Asubdiag_u(nx-1,ny-1))
       allocate(Asupdiag_u(nx-1,ny-1))
       allocate(omega_u   (nx-1,ny-1))
       allocate(denom_u   (nx-1,ny-1))

       allocate(Adiag_v   (nx-1,ny-1))
       allocate(Asubdiag_v(nx-1,ny-1))
       allocate(Asupdiag_v(nx-1,ny-1))
       allocate(omega_v   (nx-1,ny-1))
       allocate(denom_v   (nx-1,ny-1))

       ! Compute arrays for tridiagonal preconditioning

       call setup_preconditioner_tridiag_local_3d(&
            nx,            ny,            &
            nz,            parallel,      &
            active_vertex,                &
            indxA_2d,      indxA_3d,      &
            itest, jtest,  rtest,         &
            Auu,           Avv,           &
            Muu,           Mvv,           &
            Adiag_u,       Adiag_v,       &
            Asubdiag_u,    Asubdiag_v,    &
            Asupdiag_u,    Asupdiag_v,    &
            omega_u,       omega_v,       &
            denom_u,       denom_v)
       
    elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then

       ! Allocate tridiagonal preconditioning matrices
       ! Note: (i,j) indices are switched for the A_v matrices to reduce striding.
       allocate(Adiag_u   (ilocal,jlocal))
       allocate(Asubdiag_u(ilocal,jlocal))
       allocate(Asupdiag_u(ilocal,jlocal))
       allocate(omega_u(ilocal,jlocal))
       allocate(denom_u(ilocal,jlocal))
       allocate(xuh_u(ilocal,jlocal))
       allocate(xlh_u(ilocal,jlocal))

       allocate(Adiag_v   (jlocal,ilocal))
       allocate(Asubdiag_v(jlocal,ilocal))
       allocate(Asupdiag_v(jlocal,ilocal))
       allocate(omega_v(jlocal,ilocal))
       allocate(denom_v(jlocal,ilocal))
       allocate(xuh_v(jlocal,ilocal))
       allocate(xlh_v(jlocal,ilocal))

       ! These two matrices are for gathering data from all tasks on a given row or column.
       allocate(gather_data_row(8*tasks_row,jlocal))
       allocate(gather_data_col(8*tasks_col,ilocal))
       gather_data_row = 0.0d0
       gather_data_col = 0.0d0

       ! Compute arrays for tridiagonal preconditioning

       call setup_preconditioner_tridiag_global_3d(&
            nx,            ny,            &
            nz,            parallel,      &
            active_vertex,                &
            indxA_2d,      indxA_3d,      &
            ilocal,        jlocal,        &
            itest, jtest,  rtest,         &
            Auu,           Avv,           &
            Muu,           Mvv,           &
            Adiag_u,       Adiag_v,       &
            Asubdiag_u,    Asubdiag_v,    &
            Asupdiag_u,    Asupdiag_v,    &
            omega_u,       omega_v,       &
            denom_u,       denom_v,       &
            xuh_u,         xuh_v,         &
            xlh_u,         xlh_v)

    endif   ! precond

    call t_stopf("pcg_precond_init")

    !---- Initialize scalars and vectors

    niters = maxiters_chrongear
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
    bb = parallel_global_sum_stagger(worku, parallel, arr2=workv)
    call t_stopf("pcg_glbsum_init")

    ! take square root

    L2_rhs = sqrt(bb)       ! L2 norm of RHS

    !---------------------------------------------------------------
    ! First pass of algorithm
    !---------------------------------------------------------------

    iter = 1

    ! Note: The matrix A must be complete for all rows corresponding to locally 
    !        owned nodes, and x must have the correct values in
    !        halo nodes bordering the locally owned nodes.
    !       Then z = Ax will be correct for locally owned nodes.

    !---- Halo update for x (initial guess for velocity solution)

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(xu, parallel)
    call staggered_parallel_halo(xv, parallel)
    call t_stopf("pcg_halo_init")

    !---- Compute A*x   (use z as a temp vector for A*x)

    call t_startf("pcg_matmult_init")
    call matvec_multiply_structured_3d(nx,            ny,            &
                                       nz,            parallel,      &
                                       indxA_3d,      active_vertex, &
                                       Auu,           Auv,           &
                                       Avu,           Avv,           &
                                       xu,            xv,            &
                                       zu,            zv)
    call t_stopf("pcg_matmult_init")

    !---- Compute the initial residual r = b - A*x
    !---- This is correct for locally owned nodes.

    call t_startf("pcg_vecupdate")
    ru(:,:,:) = bu(:,:,:) - zu(:,:,:)
    rv(:,:,:) = bv(:,:,:) - zv(:,:,:)
    call t_stopf("pcg_vecupdate")

    !---- Halo update for residual

    call t_startf("pcg_halo_init")
    call staggered_parallel_halo(ru, parallel)
    call staggered_parallel_halo(rv, parallel)
    call t_stopf("pcg_halo_init")

    !---- Compute (PC)r = solution z of Mz = r
    !---- Since r was just updated in halo, z is correct in halo

    ! From here on, call timers with 'iter' suffix because this can be considered the first iteration
    call t_startf("pcg_precond_iter")

    if (precond == PRECONDITIONER_NONE) then   ! no preconditioning

       zu(:,:,:) = ru(:,:,:)         ! PC(r) = r     
       zv(:,:,:) = rv(:,:,:)         ! PC(r) = r    

    elseif (precond == PRECONDITIONER_DIAGONAL ) then  ! diagonal preconditioning

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

    elseif (precond == PRECONDITIONER_SIA) then   ! local vertical shallow-ice solver for preconditioning

       call easy_sia_solver(nx,   ny,   nz,        &
                            active_vertex,         &
                            Muu,  ru,   zu)      ! solve Muu*zu = ru for zu 

       call easy_sia_solver(nx,   ny,   nz,        &
                            active_vertex,         &
                            Mvv,  rv,   zv)      ! solve Mvv*zv = rv for zv

    elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then

       ! Use a local tridiagonal solver to find an approximate solution of A*z = r

       call tridiag_solver_local_3d(&
            nx,            ny,            &
            nz,            parallel,      &
            active_vertex,                &
            itest, jtest,  rtest,         &
            Adiag_u,       Adiag_v,       &  ! entries of 2D preconditioning matrix
            Asubdiag_u,    Asubdiag_v,    &
            Asupdiag_u,    Asupdiag_v,    &
            omega_u,       omega_v,       &
            denom_u,       denom_v,       &
            Muu,           Mvv,           &  ! entries of SIA matrix
            ru,            rv,            &  ! 3D residual
            zu,            zv)               ! approximate solution of Az = r

    elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then

       ! Use a global tridiagonal solver to find an approximate solution of A*z = r

       call tridiag_solver_global_3d(&
            nx,              ny,              &
            nz,              parallel,        &
            active_vertex,                    &
            ilocal,          jlocal,          &
            tasks_row,       tasks_col,       &
            itest,  jtest,   rtest,           &
            Adiag_u,         Adiag_v,         &  ! entries of 2D preconditioning matrix
            Asubdiag_u,      Asubdiag_v,      &
            Asupdiag_u,      Asupdiag_v,      &
            omega_u,         omega_v,         &
            denom_u,         denom_v,         &
            xuh_u,           xuh_v,           &
            xlh_u,           xlh_v,           &
            Muu,             Mvv,             &  ! entries of SIA matrix
            gather_data_row, gather_data_col, &
            .true.,                           &  ! first_time = T (first iteration)
            ru,              rv,              &  ! 3D residual
            zu,              zv)                 ! approximate solution of Az = r

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
    call matvec_multiply_structured_3d(nx,            ny,            &
                                       nz,            parallel,      &
                                       indxA_3d,      active_vertex, &
                                       Auu,           Auv,           &
                                       Avu,           Avv,           &
                                       du,            dv,            &
                                       qu,            qv)
    call t_stopf("pcg_matmult_iter")

    !---- Compute intermediate result for dot product (d,q) = (d,Ad)

    call t_startf("pcg_dotprod")
    work2u(:,:,:,2) = du(:,:,:) * qu(:,:,:)
    work2v(:,:,:,2) = dv(:,:,:) * qv(:,:,:)
    call t_stopf("pcg_dotprod")

    !---- Find global sums of (r,z) and (d,q)

    call t_startf("pcg_glbsum_iter")
    gsum = parallel_global_sum_stagger(work2u, 2, parallel, arr2=work2v)
    call t_stopf("pcg_glbsum_iter")

    !---- Halo update for q

    call t_startf("pcg_halo_iter")
    call staggered_parallel_halo(qu, parallel)
    call staggered_parallel_halo(qv, parallel)
    call t_stopf("pcg_halo_iter")

    rho_old = gsum(1)      ! (r,z) = (r, (PC)r)
    sigma = gsum(2)        ! (d,q) = (d, Ad)

    alpha = rho_old/sigma

    if (alpha /= alpha) then  ! alpha is NaN
!!       write(iulog,*) 'rho_old, sigma, alpha:', rho_old, sigma, alpha
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

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'iter =', iter
       write(iulog,*) 'i, k, xu, xv, ru, rv:'
       do i = itest-3, itest+3
          write(iulog,*) ' '
          do k = 1, nz
             write(iulog,'(2i4, 4f16.10)') i, k, xu(k,i,j), xv(k,i,j), ru(k,i,j), rv(k,i,j)
          enddo
       enddo  ! i
    endif

    !---------------------------------------------------------------
    ! Iterate to solution
    !---------------------------------------------------------------

    iter_loop: do iter = 2, maxiters_chrongear  ! first iteration done above

       !---- Compute PC(r) = solution z of Mz = r
       !---- z is correct in halo

       call t_startf("pcg_precond_iter")

       if (precond == PRECONDITIONER_NONE) then   ! no preconditioning

           zu(:,:,:) = ru(:,:,:)         ! PC(r) = r
           zv(:,:,:) = rv(:,:,:)         ! PC(r) = r    

       elseif (precond == PRECONDITIONER_DIAGONAL ) then  ! diagonal preconditioning

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

       elseif (precond == PRECONDITIONER_SIA) then   ! local vertical shallow-ice solver for preconditioning

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Muu,  ru,   zu)      ! solve Muu*zu = ru for zu 

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Mvv,  rv,   zv)      ! solve Mvv*zv = rv for zv

       elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then

          ! Use a local tridiagonal solver to find an approximate solution of A*z = r

          call tridiag_solver_local_3d(&
               nx,            ny,            &
               nz,            parallel,      &
               active_vertex,                &
               itest, jtest,  rtest,         &
               Adiag_u,       Adiag_v,       &  ! entries of 2D preconditioning matrix
               Asubdiag_u,    Asubdiag_v,    &
               Asupdiag_u,    Asupdiag_v,    &
               omega_u,       omega_v,       &
               denom_u,       denom_v,       &
               Muu,           Mvv,           &  ! entries of SIA matrix
               ru,            rv,            &  ! 3D residual
               zu,            zv)               ! approximate solution of Az = r

       elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then   ! tridiagonal preconditioning with global solve

          ! Use a global tridiagonal solver to find an approximate solution of A*z = r

          call tridiag_solver_global_3d(&
               nx,              ny,              &
               nz,              parallel,        &
               active_vertex,                    &
               ilocal,          jlocal,          &
               tasks_row,       tasks_col,       &
               itest,  jtest,   rtest,           &
               Adiag_u,         Adiag_v,         &  ! entries of 2D preconditioning matrix
               Asubdiag_u,      Asubdiag_v,      &
               Asupdiag_u,      Asupdiag_v,      &
               omega_u,         omega_v,         &
               denom_u,         denom_v,         &
               xuh_u,           xuh_v,           &
               xlh_u,           xlh_v,           &
               Muu,             Mvv,             &  ! entries of SIA matrix
               gather_data_row, gather_data_col, &
               .false.,                          &  ! first_time = F (iteration 2+)
               ru,              rv,              &  ! 3D residual
               zu,              zv)                 ! approximate solution of Az = r

       endif    ! precond

       call t_stopf("pcg_precond_iter")

       !---- Compute Az = A*z
       !---- This is the one matvec multiply required per iteration
       !---- Az is correct for local owned nodes and needs a halo update (below)

       call t_startf("pcg_matmult_iter")
       call matvec_multiply_structured_3d(nx,            ny,            &
                                          nz,            parallel,      &
                                          indxA_3d,      active_vertex, &
                                          Auu,           Auv,           &
                                          Avu,           Avv,           &
                                          zu,            zv,            &
                                          Azu,           Azv)
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
       gsum = parallel_global_sum_stagger(work2u, 2, parallel, arr2=work2v)
       call t_stopf("pcg_glbsum_iter")

       !---- Halo update for Az
       !---- This is the one halo update required per iteration

       call t_startf("pcg_halo_iter")
       call staggered_parallel_halo(Azu, parallel)
       call staggered_parallel_halo(Azv, parallel)
       call t_stopf("pcg_halo_iter")

       !---- Compute some scalars

       rho = gsum(1)        ! (r,z)
       delta = gsum(2)      ! (Az,z)

       beta = rho / rho_old
       sigma = delta - beta*rho/alpha
       alpha = rho / sigma
       rho_old = rho        ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i))

       if (alpha /= alpha) then  ! alpha is NaN
!!          write(iulog,*) 'rho, sigma, alpha:', rho, sigma, alpha
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

       if (verbose_pcg .and. this_rank == rtest) then
          j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'iter =', iter
          write(iulog,*) 'i, k, xu, xv, ru, rv:'
          do i = itest-3, itest+3
             write(iulog,*) ' '
             do k = 1, nz
                write(iulog,'(2i4, 4f16.10)') i, k, xu(k,i,j), xv(k,i,j), ru(k,i,j), rv(k,i,j)
             enddo
          enddo  ! i
       endif

       ! Check for convergence every linear_solve_ncheck iterations.
       ! Also check at iter = 5, to reduce iterations when the nonlinear solver is close to convergence.
       ! TODO: Check at iter = linear_solve_ncheck/2 instead of 5?  This would be answer-changing.
       !
       ! For convergence check, use r = b - Ax

       if (mod(iter, linear_solve_ncheck) == 0 .or. iter == 5) then

          if (verbose_pcg .and. main_task) then
             write(iulog,*) ' '
             write(iulog,*) 'Check convergence, iter =', iter
          endif

          !---- Compute z = A*x  (use z as a temp vector for A*x)
           
          call t_startf("pcg_matmult_resid")
          call matvec_multiply_structured_3d(nx,            ny,            &
                                             nz,            parallel,      &
                                             indxA_3d,      active_vertex, &
                                             Auu,           Auv,           &
                                             Avu,           Avv,           &
                                             xu,            xv,            &
                                             zu,            zv)
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
          rr = parallel_global_sum_stagger(worku, parallel, arr2=workv)
          call t_stopf("pcg_glbsum_resid")

          ! take square root
          L2_resid = sqrt(rr)

          ! compute normalized error
          err = L2_resid/L2_rhs

          if (verbose_pcg .and. main_task) then
             write(iulog,*) 'iter, L2_resid, L2_rhs, error =', iter, L2_resid, L2_rhs, err
          endif

          !WHL - debug
          if (verbose_pcg .and. this_rank == rtest) then
             ru_max = 0.d0
             rv_max = 0.d0
             iu_max = 0
             ju_max = 0
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (sum(abs(ru(:,i,j))) > ru_max) then
                      ru_max = sum(abs(ru(:,i,j)))
                      iu_max = i
                      ju_max = j
                   endif
                   if (sum(abs(rv(:,i,j))) > rv_max) then
                      rv_max = sum(abs(rv(:,i,j)))
                      iv_max = i
                      jv_max = j
                   endif
                enddo
             enddo
             write(iulog,*) 'r, i, j, ru_max:', this_rank, iu_max, ju_max, ru_max
             write(iulog,*) 'r, i, j, rv_max:', this_rank, iv_max, jv_max, rv_max
          endif

          ! If converged, then exit the loop.
          ! Note: Without good preconditioning, convergence can be slow,
          !       but the solution after maxiters_chrongear might be good enough.

          if (err < tolerance) then
             niters = iter
             if (verbose_pcg .and. main_task) then
                write(iulog,*) 'Glissade PCG solver has converged, iter =', niters
                write(iulog,*) ' '
             endif
             exit iter_loop
          elseif (iter == maxiters_chrongear) then
             if (verbose_pcg .and. main_task) then
                write(iulog,*) 'Glissade PCG solver did not converge'
                write(iulog,*) 'iter, err, tolerance:', iter, err, tolerance
                write(iulog,*) ' '
             endif
          endif

          !---- Update residual in halo for next iteration

          call t_startf("pcg_halo_resid")
          call staggered_parallel_halo(ru, parallel)
          call staggered_parallel_halo(rv, parallel)
          call t_stopf("pcg_halo_resid")

       endif    ! linear_solve_ncheck

    enddo iter_loop

    ! Clean up
    if (allocated(Adiag_u))    deallocate(Adiag_u, Adiag_v)
    if (allocated(Asubdiag_u)) deallocate(Asubdiag_u, Asubdiag_v)
    if (allocated(Asupdiag_u)) deallocate(Asupdiag_u, Asupdiag_v)
    if (allocated(omega_u))    deallocate(omega_u, omega_v)
    if (allocated(denom_u))    deallocate(denom_u, denom_v)
    if (allocated(xuh_u))      deallocate(xuh_u, xuh_v)
    if (allocated(xlh_u))      deallocate(xlh_u, xlh_v)
    if (allocated(gather_data_row)) deallocate(gather_data_row)
    if (allocated(gather_data_col)) deallocate(gather_data_col)

  end subroutine pcg_solver_chrongear_3d

!****************************************************************************

  subroutine setup_preconditioner_diag_3d(nx,         ny,         &
                                          nz,         indxA_3d,   &
                                          Auu,        Avv,        &
                                          Adiagu,     Adiagv)

    ! Set up diagonal preconditioning matrices for 3D solve
    ! Note: This is likely to be a poor preconditioner for 3D problems.

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
         nx, ny,             &  ! horizontal grid dimensions (for scalars)
                                ! velocity grid has dimensions (nx-1,ny-1)
         nz                     ! number of vertical levels where velocity is computed

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
         indxA_3d               ! maps relative (x,y,z) coordinates to an index between 1 and 27

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
         Adiagu, Adiagv         ! matrix entries for diagonal preconditioning

    integer :: i, j, m

    ! Form diagonal matrix for preconditioning

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Setting up diagonal preconditioner'
    endif  ! verbose_pcg

    m = indxA_3d(0,0,0)
    Adiagu(:,:,:) = Auu(m,:,:,:)
    Adiagv(:,:,:) = Avv(m,:,:,:)

  end subroutine setup_preconditioner_diag_3d

!****************************************************************************

  subroutine setup_preconditioner_sia_3d(nx,        ny,       &
                                         nz,        indxA_3d, &
                                         Auu,       Avv,      &
                                         Muu,       Mvv)

    ! Set up preconditioning matrices for the SIA-based solver

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
                              ! velocity grid has dimensions (nx-1,ny-1)
       nz                     ! number of vertical levels where velocity is computed

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA_3d               ! maps relative (x,y,z) coordinates to an index between 1 and 27

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Avv               ! two out of the four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(out) ::   &
       Muu, Mvv               ! preconditioning matrices based on shallow-ice approximation

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, k, m

    ! Initialize
    Muu (:,:,:,:) = 0.d0    
    Mvv (:,:,:,:) = 0.d0    

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Setting up shallow-ice preconditioner'
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

           m = indxA_3d(0,0,-1)
           Muu(-1,k,i,j) = Auu(m,k,i,j)
           Mvv(-1,k,i,j) = Avv(m,k,i,j)

           m = indxA_3d(0,0,0)
           Muu( 0,k,i,j) = Auu(m,k,i,j)
           Mvv( 0,k,i,j) = Avv(m,k,i,j)

           m = indxA_3d(0,0,1)
           Muu( 1,k,i,j) = Auu(m,k,i,j)
           Mvv( 1,k,i,j) = Avv(m,k,i,j)
        enddo
      enddo
    enddo

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Done in shallow-ice preconditioner'
    endif  ! verbose_pcg

  end subroutine setup_preconditioner_sia_3d

!****************************************************************************

  subroutine setup_preconditioner_tridiag_local_3d(&
       nx,            ny,            &
       nz,            parallel,      &
       active_vertex,                &
       indxA_2d,      indxA_3d,      &
       itest, jtest,  rtest,         &
       Auu,           Avv,           &
       Muu,           Mvv,           &
       Adiag_u,       Adiag_v,       &
       Asubdiag_u,    Asubdiag_v,    &
       Asupdiag_u,    Asupdiag_v,    &
       omega_u,       omega_v,       &
       denom_u,       denom_v)

    ! Set up arrays for the local tridiagonal preconditioner.
    ! This preconditioner uses a combination of horizontal tridiagonal solutionss and
    !  vertical SIA solutions, so we compute arrays needed for both solutions.

    integer, intent(in) :: &
         nx, ny,             &  ! horizontal grid dimensions (for scalars);
                                ! velocity grid has dimensions (nx-1,ny-1)
         nz                     ! number of vertical levels where velocity is computed

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    logical, dimension(nx-1,ny-1), intent(in) ::   &
         active_vertex          ! T for columns (i,j) where velocity is computed, else F

    integer, dimension(-1:1,-1:1), intent(in) :: &
         indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
         indxA_3d               ! maps relative (x,y,z) coordinates to an index between 1 and 27

    integer, intent(in) :: itest, jtest, rtest   ! coordinates of diagnostic point

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
         Auu, Avv               ! two out of the four components of assembled matrix
                                ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                                ! other dimensions = (z,x,y) indices
                                !
                                !    Auu  | Auv
                                !    _____|____
                                !    Avu  | Avv
                                !         |

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(out) ::   &
         Muu, Mvv               ! preconditioning matrices based on shallow-ice approximation

    real(dp), dimension(nx-1,ny-1), intent(out) :: &
         Asubdiag_u, Adiag_u, Asupdiag_u,  &  ! matrix entries from Auu for tridiagonal preconditioning
         Asubdiag_v, Adiag_v, Asupdiag_v,  &  ! matrix entries from Avv for tridiagonal preconditioning
         omega_u, omega_v,                 &  ! work arrays for tridiagonal solve
         denom_u, denom_v

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, k, m

    real(dp), dimension(:,:,:), allocatable :: &
         Auu_vsum, Avv_vsum          ! arrays to hold vertical sums of matrix elements

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Setting up local tridiagonal preconditioner'
    endif  ! verbose_pcg

    ! Construct arrays consisting of vertical sums of matrix elements.
    ! The resulting arrays are approximations of the 2D matrices for an SSA solution.
    ! Elements m, m+9 and m+18 (for m = 1 to 9) are vertically stacked,
    !  corresponding to (k-1:k+1,ii,jj) for a given (ii,jj) adjacent to cell (i,j).
    ! TODO: Switch index order for 3D arrays in glissade_velo_higher.F90, putting m last.

    allocate(Auu_vsum(nx-1,ny-1,9))
    allocate(Avv_vsum(nx-1,ny-1,9))

    do j = 1, ny-1
       do i = 1, nx-1
          do m = 1, 9
             Auu_vsum(i,j,m) = 0.0d0
             Avv_vsum(i,j,m) = 0.0d0
             do k = 1, nz
                Auu_vsum(i,j,m) = Auu_vsum(i,j,m) + Auu(m,k,i,j) + Auu(m+9,k,i,j) + Auu(m+18,k,i,j)
                Avv_vsum(i,j,m) = Avv_vsum(i,j,m) + Avv(m,k,i,j) + Avv(m+9,k,i,j) + Avv(m+18,k,i,j)
             enddo
          enddo
       enddo
    enddo

    ! Adjust the matrix for Dirichlet points.
    ! At these points, Auu has a diagonal element of 1.0 at each vertical level.
    ! We want the vertically summed matrix also to have a diagonal element of 1.0, instead of 1.0*nz.
    do j = 1, ny-1
       do i = 1, nx-1
          if (.not.active_vertex(i,j)) then
             Auu_vsum(i,j,:) = 0.0d0
             Avv_vsum(i,j,:) = 0.0d0
             Auu_vsum(i,j,5) = 1.0d0
             Avv_vsum(i,j,5) = 1.0d0
          endif
       enddo
    enddo

    !WHL - debug
!    if (verbose_tridiag .and. this_rank==rtest) then
    if (verbose_pcg .and. this_rank==rtest) then
       i = itest
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'r, i, j =', this_rank, i, j
       write(iulog,*) 'Auu(nz, 1: 9) =', Auu(1:9,nz,i,j)
       write(iulog,*) 'Auu(nz,10:18) =', Auu(10:18,nz,i,j)
       write(iulog,*) 'Auu(nz,19:27) =', Auu(19:27,nz,i,j)
       write(iulog,*) ' '
       write(iulog,*) 'm, Auu_vsum(i,j,m), Avv_vsum(i,j,m):'
       do m = 1,9
          write(iulog,*) m, Auu_vsum(i,j,m), Avv_vsum(i,j,m)
       enddo
    endif

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'Before global tridiag PC setup, r, j =', rtest, jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, Auu_vsum(1), Auu_sum(5), Auu_sum(9):'
       do i = itest-3, itest+3
          write(iulog,'(i4, 3f15.7)') i, Auu_vsum(i,j,1), Auu_vsum(i,j,5), Auu_vsum(i,j,9)
       enddo
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_u, Asubdiag_u, Asupdiag_u, active'
       do i = itest-3, itest+3
          write(iulog,'(i4, 3f15.7, l4)') &
               i, Adiag_u(i,j), Asubdiag_u(i,j), Asupdiag_u(i,j), active_vertex(i,j)
       enddo
    endif

    call setup_preconditioner_tridiag_local_2d(&
         nx,            ny,             &
         parallel,      indxA_2d,       &
         itest, jtest,  rtest,          &
         Auu_vsum,      Avv_vsum,       &
         Adiag_u,       Adiag_v,        &
         Asubdiag_u,    Asubdiag_v,     &
         Asupdiag_u,    Asupdiag_v,     &
         omega_u,       omega_v,        &
         denom_u,       denom_v)

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'After local tridiag PC setup, r, j =', rtest, jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_u, Asubdiag_u, Asupdiag_u, omega_u, denom_u:'
       do i = itest-3, itest+3
          write(iulog,'(i4, 5f15.7)') &
               i, Adiag_u(i,j), Asubdiag_u(i,j), Asupdiag_u(i,j), omega_u(i,j), denom_u(i,j)
       enddo
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_v, Asubdiag_v, Asupdiag_v, omega_v, denom_v:'
       do i = itest-3, itest+3
          write(iulog,'(i4, 5f15.7)') &
               i, Adiag_v(i,j), Asubdiag_v(i,j), Asupdiag_v(i,j), omega_v(i,j), denom_v(i,j)
       enddo
    endif

    ! Set up a shallow-ice preconditioner to introduce vertical shear in the column.

    call setup_preconditioner_sia_3d(&
         nx,        ny,       &
         nz,        indxA_3d, &
         Auu,       Avv,      &
         Muu,       Mvv)

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) 'SIA PC with tridiag:'
       write(iulog,*) ' '
       write(iulog,*) 'i, k, Muu_sia, Mvv_sia:'
       do i = itest-3, itest+3
          write(iulog,*) ' '
          do k = 1, nz
             write(iulog,'(2i4, 6e13.5)') i, k, Muu(:,k,i,j), Mvv(:,k,i,j)
          enddo
       enddo
    endif

  end subroutine setup_preconditioner_tridiag_local_3d

!****************************************************************************

  subroutine setup_preconditioner_tridiag_global_3d(&
       nx,            ny,            &
       nz,            parallel,      &
       active_vertex,                &
       indxA_2d,      indxA_3d,      &
       ilocal,        jlocal,        &
       itest, jtest,  rtest,         &
       Auu,           Avv,           &
       Muu,           Mvv,           &
       Adiag_u,       Adiag_v,       &
       Asubdiag_u,    Asubdiag_v,    &
       Asupdiag_u,    Asupdiag_v,    &
       omega_u,       omega_v,       &
       denom_u,       denom_v,       &
       xuh_u,         xuh_v,         &
       xlh_u,         xlh_v)

    ! Set up arrays for the local tridiagonal preconditioner.
    ! This preconditioner uses a combination of horizontal tridiagonal solutionss and
    !  vertical SIA solutions, so we compute arrays needed for both solutions.

    integer, intent(in) :: &
         nx, ny,                &  ! horizontal grid dimensions (for scalars)
                                   ! velocity grid has dimensions (nx-1,ny-1)
         nz                        ! number of vertical levels where velocity is computed

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    logical, dimension(nx-1,ny-1), intent(in) ::   &
         active_vertex          ! T for columns (i,j) where velocity is computed, else F

    integer, dimension(-1:1,-1:1), intent(in) :: &
         indxA_2d             ! maps relative (x,y) coordinates to an index between 1 and 9

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
         indxA_3d             ! maps relative (x,y,z) coordinates to an index between 1 and 27

    integer, intent(in) :: &
         ilocal, jlocal       ! size of input/output arrays; number of locally owned vertices in each direction

    integer, intent(in) :: itest, jtest, rtest   ! coordinates of diagnostic point

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
         Auu, Avv             ! two out of the four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(out) ::   &
         Muu, Mvv               ! preconditioning matrices based on shallow-ice approximation


    real(dp), dimension(ilocal, jlocal), intent(out) :: &
         Asubdiag_u, Adiag_u, Asupdiag_u,  &  ! matrix entries from Auu for tridiagonal preconditioning
         omega_u, denom_u,                 &  ! work arrays for tridiagonal solve
         xuh_u,   xlh_u

    ! Note: For the v arrays, the i and j indices are switched to reduce strides
    real(dp), dimension(jlocal, ilocal), intent(out) :: &
         Asubdiag_v, Adiag_v, Asupdiag_v,  &  ! matrix entries from Avv for tridiagonal preconditioning
         omega_v, denom_v,                 &  ! work arrays for tridiagonal solve
         xuh_v,   xlh_v

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, k, m
    integer :: ii, jj

    real(dp), dimension(:,:,:), allocatable :: &
         Auu_vsum, Avv_vsum          ! arrays to hold vertical sums of matrix elements

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Setting up global tridiagonal preconditioner'
    endif  ! verbose_pcg

    ! Construct arrays consisting of vertical sums of matrix elements.
    ! The resulting arrays are approximations of the 2D matrices for an SSA solution.
    ! Elements m, m+9 and m+18 (for m = 1 to 9) are vertically stacked,
    !  corresponding to (k-1:k+1,ii,jj) for a given (ii,jj) adjacent to cell (i,j).
    ! TODO: Switch index order for 3D arrays in glissade_velo_higher.F90, putting m last.
    allocate(Auu_vsum(nx-1,ny-1,9))
    allocate(Avv_vsum(nx-1,ny-1,9))

    do j = 1, ny-1
       do i = 1, nx-1
          do m = 1, 9
             Auu_vsum(i,j,m) = 0.0d0
             Avv_vsum(i,j,m) = 0.0d0
             do k = 1, nz
                Auu_vsum(i,j,m) = Auu_vsum(i,j,m) + Auu(m,k,i,j) + Auu(m+9,k,i,j) + Auu(m+18,k,i,j)
                Avv_vsum(i,j,m) = Avv_vsum(i,j,m) + Avv(m,k,i,j) + Avv(m+9,k,i,j) + Avv(m+18,k,i,j)
             enddo
          enddo
       enddo
    enddo

    ! Adjust the matrix for Dirichlet points.
    ! At these points, Auu has a diagonal element of 1.0 at each vertical level.
    ! We want the vertically summed matrix also to have a diagonal element of 1.0, instead of 1.0*nz.
    do j = 1, ny-1
       do i = 1, nx-1
          if (.not.active_vertex(i,j)) then
             Auu_vsum(i,j,5) = 1.0d0
             Avv_vsum(i,j,5) = 1.0d0
          endif
       enddo
    enddo

    !WHL - debug
    do m = 1, 9
       call staggered_parallel_halo(Auu_vsum(:,:,m), parallel)
       call staggered_parallel_halo(Avv_vsum(:,:,m), parallel)
    enddo

    !WHL - debug
!!    if (verbose_tridiag .and. this_rank==rtest) then
    if (verbose_pcg .and. this_rank==rtest) then
       i = itest
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'r, i, j =', this_rank, i, j
       write(iulog,*) 'Auu(nz, 1: 9) =', Auu(1:9,nz,i,j)
       write(iulog,*) 'Auu(nz,10:18) =', Auu(10:18,nz,i,j)
       write(iulog,*) 'Auu(nz,19:27) =', Auu(19:27,nz,i,j)
       write(iulog,*) ' '
       write(iulog,*) 'm, Auu_vsum(i,j,m), Avv_vsum(i,j,m):'
       do m = 1,9
          write(iulog,*) m, Auu_vsum(i,j,m), Avv_vsum(i,j,m)
       enddo
    endif

    ! Extract tridiagonal matrix entries from Auu_vsum
    do j = 1, jlocal
       jj = j + staggered_jlo - 1
       do i = 1, ilocal
          ii = i + staggered_ilo - 1
          Asubdiag_u(i,j) = Auu_vsum(ii,jj,indxA_2d(-1,0))   ! subdiagonal elements
          Adiag_u   (i,j) = Auu_vsum(ii,jj,indxA_2d( 0,0))   ! diagonal elements
          Asupdiag_u(i,j) = Auu_vsum(ii,jj,indxA_2d( 1,0))   ! superdiagonal elements
       enddo
    enddo

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'Before global tridiag PC setup, r, j =', rtest, jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, Auu_vsum(1), Auu_sum(5), Auu_sum(9):'
       do i = itest-5, itest+5
          write(iulog,'(i4, 3f15.7)') i, Auu_vsum(i,j,1), Auu_vsum(i,j,5), Auu_vsum(i,j,9)
       enddo
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_u, Asubdiag_u, Asupdiag_u, active'
       do i = itest-3, itest+3
          write(iulog,'(i4, 3f15.7, l4)') &
               i, Adiag_u(i,j), Asubdiag_u(i,j), Asupdiag_u(i,j), active_vertex(i,j)
       enddo
    endif

    ! compute work arrays for the u solve in each matrix row
    call setup_preconditioner_tridiag_global_2d(&
         ilocal,       jlocal,      &
!!          itest, jtest, rtest,       &
         itest - staggered_ilo + 1, &  ! itest referenced to (ilocal,jlocal) coordinates
         jtest - staggered_jlo + 1, &  ! jtest referenced to (ilocal,jlocal) coordinates
         rtest,                     &
         Adiag_u,                   &
         Asubdiag_u,   Asupdiag_u,  &
         omega_u,      denom_u,     &
         xuh_u,        xlh_u)

    ! Extract tridiagonal matrix entries from Avv_vsum
    do i = 1, ilocal
       ii = i + staggered_ilo - 1
       do j = 1, jlocal
          jj = j + staggered_jlo - 1
          Asubdiag_v(j,i) = Avv_vsum(ii,jj,indxA_2d(0,-1))   ! subdiagonal elements
          Adiag_v   (j,i) = Avv_vsum(ii,jj,indxA_2d(0, 0))   ! diagonal elements
          Asupdiag_v(j,i) = Avv_vsum(ii,jj,indxA_2d(0, 1))   ! superdiagonal elements
       enddo
    enddo

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'Before global tridiag PC setup, r, j =', rtest, jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_v, Asubdiag_v, Asupdiag_v, active'
       do i = itest-3, itest+3
          write(iulog,'(i4, 3f15.7, l4)') &
               i, Adiag_v(i,j), Asubdiag_v(i,j), Asupdiag_v(i,j), active_vertex(i,j)
       enddo
    endif

    ! compute work arrays for the v solve in each matrix column
    ! Note: The *_v arrays have dimensions (jlocal,ilocal) to reduce strides

    call setup_preconditioner_tridiag_global_2d(&
         jlocal,       ilocal,      &
         !!          itest, jtest, rtest,       &
         jtest - staggered_jlo + 1, &  ! jtest referenced to (jlocal,ilocal) coordinates
         itest - staggered_ilo + 1, &  ! itest referenced to (jlocal,ilocal) coordinates
         rtest,                     &
         Adiag_v,                   &
         Asubdiag_v,   Asupdiag_v,  &
         omega_v,      denom_v,     &
         xuh_v,        xlh_v)

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'After global tridiag PC setup, r, j =', rtest, jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_u, Asubdiag_u, Asupdiag_u, omega_u, denom_u, xuh_u, xlh_u, active_vertex:'
       do i = itest-3, itest+3
          write(iulog,'(i4, 7f15.7, l4)') &
               i, Adiag_u(i,j), Asubdiag_u(i,j), Asupdiag_u(i,j), &
               omega_u(i,j), denom_u(i,j), xuh_u(i,j), xlh_u(i,j), active_vertex(i,j)
       enddo
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_v, Asubdiag_v, Asupdiag_v, omega_v, denom_v, xuh_v, xlh_v::'
       do i = itest-3, itest+3
          write(iulog,'(i4, 7f15.7)') &
               i, Adiag_v(i,j), Asubdiag_v(i,j), Asupdiag_v(i,j), &
               omega_v(i,j), denom_v(i,j), xuh_v(i,j), xlh_v(i,j)
       enddo
    endif

    ! Set up a shallow-ice preconditioner to introduce vertical shear in the column.

    call setup_preconditioner_sia_3d(&
         nx,        ny,       &
         nz,        indxA_3d, &
         Auu,       Avv,      &
         Muu,       Mvv)

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) 'SIA PC with tridiag:'
       write(iulog,*) ' '
       write(iulog,*) 'i, k, Muu_sia, Mvv_sia:'
       do i = itest-3, itest+3
          write(iulog,*) ' '
          do k = 1, nz
             write(iulog,'(2i4, 6e13.5)') i, k, Muu(:,k,i,j), Mvv(:,k,i,j)
          enddo
       enddo
    endif

  end subroutine setup_preconditioner_tridiag_global_3d

!****************************************************************************

  subroutine tridiag_solver_local_3d(&
       nx,            ny,            &
       nz,            parallel,      &
       active_vertex,                &
       itest, jtest,  rtest,         &
       Adiag_u,       Adiag_v,       &  ! tridiagonal entries of 3D matrix
       Asubdiag_u,    Asubdiag_v,    &
       Asupdiag_u,    Asupdiag_v,    &
       omega_u,       omega_v,       &
       denom_u,       denom_v,       &
       Muu,           Mvv,           &  ! entries of SIA preconditioning matrix
       ru,            rv,            &
       zu,            zv)

    ! Find the approximate solution of A*z = r, where A and r are the 3D matrix and rhs.
    ! This is done in several steps:
    ! (1) Solve a local (1-processor) tridiagonal problem, A_2d * z_2d = r_2d, where A_2d and r_2d are
    !     vertically summed version of the 3D matrix and rhs.
    ! (2) Solve a vertical SIA-based problem, M_sia*z_sia = r_sia for each velocity component in each local column.
    ! (3) Assume an approximate solution z(k) = z_2d + z_sia(k) at each level k in each column.
    !
    ! Without the SIA adjustment (2), the preconditioning cannot generate solutions with vertical shear.

    integer, intent(in) :: &
         nx, ny,                &  ! horizontal grid dimensions (for scalars);
                                   ! velocity grid has dimensions (nx-1,ny-1)
         nz                        ! number of vertical levels where velocity is computed

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    logical, dimension(nx-1,ny-1), intent(in) ::   &
         active_vertex          ! T for columns (i,j) where velocity is computed, else F

    integer, intent(in) :: &
         itest, jtest, rtest    ! coordinates of diagnostic point

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
         Asubdiag_u, Adiag_u, Asupdiag_u,  &  ! matrix entries from Auu for tridiagonal preconditioning
         Asubdiag_v, Adiag_v, Asupdiag_v,  &  ! matrix entries from Avv for tridiagonal preconditioning
         omega_u, omega_v,                 &  ! work arrays for tridiagonal solve
         denom_u, denom_v

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(in) ::   &
         Muu, Mvv               ! preconditioning matrices based on shallow-ice approximation

    real(dp), dimension(nz,nx-1,ny-1), intent(in) :: &
         ru,         rv         ! 3d residual; rhs of A*z = r

    real(dp), dimension(nz,nx-1,ny-1), intent(out) :: &
         zu,         zv         ! 3d solution of A*z = r

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j

    real(dp), dimension(nx-1,ny-1) ::  &
         ru_vsum, rv_vsum,     & ! vertical sum of residual
         zu_vsum, zv_vsum        ! solution of vertically summed problem

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Applying local tridiagonal preconditioner'
    endif

    do j = 1, ny-1
       do i = 1, nx-1
          ru_vsum(i,j) = sum(ru(:,i,j))
          rv_vsum(i,j) = sum(rv(:,i,j))
       enddo
    enddo

    if (verbose_tridiag .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'jtest =', jtest
       write(iulog,*) 'i, ru_vsum, rv_vsum:'
       do i = itest-3, itest+3
          write(iulog,'(i4, 2f15.10)') i, ru_vsum(i,j), rv_vsum(i,j)
       enddo
    endif

    ! Solve A_2d*z_2d = r_2d

    call tridiag_solver_local_2d(nx,            ny,            &
                                 parallel,                     &
                                 itest, jtest,  rtest,         &
                                 Adiag_u,       Adiag_v,       &  ! entries of 2D preconditioning matrix
                                 Asubdiag_u,    Asubdiag_v,    &
                                 Asupdiag_u,    Asupdiag_v,    &
                                 omega_u,       omega_v,       &
                                 denom_u,       denom_v,       &
                                 ru_vsum,       rv_vsum,       &  ! right hand side, vertically integrated
                                 zu_vsum,       zv_vsum)          ! solution to vertically integrated problem

    !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
    !TODO: See whether tridiag solvers could be modified to provide zu and zv in halo cells?
    call staggered_parallel_halo(zu_vsum, parallel)
    call staggered_parallel_halo(zv_vsum, parallel)

    ! Solve a tridiagonal problem for zu and zv in each problem, to account for vertical shear

    call easy_sia_solver(nx,   ny,   nz,        &
                         active_vertex,         &
                         Muu,  ru,   zu)      ! solve Muu*zu = ru for zu

    call easy_sia_solver(nx,   ny,   nz,        &
                         active_vertex,         &
                         Mvv,  rv,   zv)      ! solve Mvv*zv = rv for zv

    !WHL - debug
    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, zu_ssa, zu_sia(1), zu_sia(nz), zu_sum(1):'
       do i = itest-3, itest+3
          write(iulog,'(i4, 4e16.8)') i, zu_vsum(i,j), zu(1,i,j), zu(nz,i,j), zu_vsum(i,j) + zu(1,i,j)
       enddo  ! i
       write(iulog,*) ' '
       write(iulog,*) 'i, zv_ssa, zv_sia(1), zv_sia(nz), zv_sum(1):'
       do i = itest-3, itest+3
          write(iulog,'(i4, 4e16.8)') i, zv_vsum(i,j), zv(1,i,j), zv(nz,i,j), zv_vsum(i,j) + zv(1,i,j)
       enddo
    endif

    ! Combine the 2d solution (zu_vsum and zv_vsum) with the SIA solution
    do j = 1, ny-1
       do i = 1, nx-1
          zu(:,i,j) = zu(:,i,j) + zu_vsum(i,j)
          zv(:,i,j) = zv(:,i,j) + zv_vsum(i,j)
       enddo
    enddo

  end subroutine tridiag_solver_local_3d

!****************************************************************************

  subroutine tridiag_solver_global_3d(&
       nx,              ny,              &
       nz,              parallel,        &
       active_vertex,                    &
       ilocal,          jlocal,          &
       tasks_row,       tasks_col,       &
       itest,  jtest,   rtest,           &
       Adiag_u,         Adiag_v,         &
       Asubdiag_u,      Asubdiag_v,      &
       Asupdiag_u,      Asupdiag_v,      &
       omega_u,         omega_v,         &
       denom_u,         denom_v,         &
       xuh_u,           xuh_v,           &
       xlh_u,           xlh_v,           &
       Muu,             Mvv,             &
       gather_data_row, gather_data_col, &
       first_time,                       &
       ru,              rv,              &
       zu,              zv)

    ! Find the approximate solution of A*z = r, where A and r are the 3D matrix and rhs.
    ! This is done in several steps:
    ! (1) Solve a global tridiagonal problem, A_2d * z_2d = r_2d, where A_2d and r_2d are
    !     vertically summed version of the 3D matrix and rhs.
    ! (2) Solve a vertical SIA-based problem, M_sia*z_sia = r_sia for each velocity component in each local column.
    ! (3) Assume an approximate solution z(k) = z_2d + z_sia(k) at each level k in each column.
    !
    ! Without the SIA adjustment (2), the preconditioning cannot generate solutions with vertical shear.

    integer, intent(in) :: &
         nx, ny,                       &  ! horizontal grid dimensions (for scalars);
                                          ! velocity grid has dimensions (nx-1,ny-1)
         nz                               ! number of vertical levels where velocity is computed

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    logical, dimension(nx-1,ny-1), intent(in) ::   &
         active_vertex          ! T for columns (i,j) where velocity is computed, else F

    integer, intent(in) :: &
         ilocal, jlocal,      & ! size of input/output arrays; number of locally owned vertices in each direction
         tasks_row, tasks_col   ! number of tasks per row and column for tridiagonal solve

    integer, intent(in) :: &
         itest, jtest, rtest    ! coordinates of diagnostic point

    real(dp), dimension(ilocal,jlocal), intent(in) :: &
         Asubdiag_u, Adiag_u, Asupdiag_u,  &  ! matrix entries from Auu for tridiagonal preconditioning
         omega_u,    denom_u,              &  ! work arrays for tridiagonal solve
         xuh_u,      xlh_u

    ! Note: For the v arrays, the i and j indicess are switched to reduce strides
    real(dp), dimension(jlocal,ilocal), intent(in) :: &
         Asubdiag_v, Adiag_v, Asupdiag_v,  &  ! matrix entries from Avv for tridiagonal preconditioning
         omega_v,    denom_v,              &  ! work arrays for tridiagonal solve
         xuh_v,      xlh_v

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(in) ::   &
         Muu, Mvv               ! preconditioning matrices based on shallow-ice approximation

    real(dp), dimension(8*tasks_row,jlocal), intent(inout) ::  &
         gather_data_row        ! array for gathering data from every task on a row

    !WHL - Earlier code version had 8*tasks_row here; pretty sure that was a bug
    real(dp), dimension(8*tasks_col,ilocal), intent(inout) ::  &
         gather_data_col        ! array for gathering data from every task on a column

    logical, intent(in) :: &
         first_time             ! flag, = T the first iteration of the solver, else = F

    real(dp), dimension(nz,nx-1,ny-1), intent(in) :: &
         ru,         rv         ! 3d residual; rhs of A*z = r

    real(dp), dimension(nz,nx-1,ny-1), intent(out) :: &
         zu,         zv         ! 3d solution of A*z = r

    ! local variables

    integer :: i, j, k, ii, jj

    real(dp), dimension(nx-1,ny-1) ::  &
         ru_vsum, rv_vsum,    & ! vertical sum of residual
         zu_vsum, zv_vsum       ! solution of vertically summed problem

    real(dp), dimension(ilocal, jlocal) ::  &
         b_u,     x_u           ! rhs and solution for tridiagonal solve, u component

    real(dp), dimension(jlocal, ilocal) ::  &
         b_v,     x_v           ! rhs and solution for tridiagonal solve, v component

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Applying global tridiagonal preconditioner'
    endif

    ! Compute vertical sums of ru and rv
    do j = 1, ny-1
       do i = 1, nx-1
          ru_vsum(i,j) = sum(ru(:,i,j))
          rv_vsum(i,j) = sum(rv(:,i,j))
       enddo
    enddo

    if (verbose_tridiag .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'jtest =', jtest
       write(iulog,*) 'i, ru_vsum, rv_vsum:'
       do i = itest-3, itest+3
          write(iulog,'(i4, 2f15.10)') i, ru_vsum(i,j), rv_vsum(i,j)
       enddo
    endif

    ! convert ru_vsum(nx-1,ny-1) to b_u(ilocal,jlocal)

    do j = 1, jlocal
       jj = j + staggered_jlo - 1
       do i = 1, ilocal
          ii = i + staggered_ilo - 1
          b_u(i,j) = ru_vsum(ii,jj)
       enddo
    enddo


    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'Before global tridiag PC u solve, r, j =', rtest, jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_u, Asubdiag_u, Asupdiag_u, b_u:'
       do i = itest-3, itest+3
          write(iulog,'(i4, 4e16.8)') i, Adiag_u(i,j), Asubdiag_u(i,j), Asupdiag_u(i,j), b_u(i,j)
       enddo
    endif

    call tridiag_solver_global_2d(ilocal,       jlocal,      &
                                  parallel,     tasks_row,   &
                                  'row',                     &  ! tridiagonal solve for each row
                                  itest - staggered_ilo + 1, &  ! itest referenced to (ilocal,jlocal) coordinates
                                  jtest - staggered_jlo + 1, &  ! jtest referenced to (ilocal,jlocal) coordinates
                                  rtest,                     &
                                  Adiag_u,                   &
                                  Asubdiag_u,   Asupdiag_u,  &
                                  omega_u,      denom_u,     &
                                  xuh_u,        xlh_u,       &
                                  b_u,          x_u,         &
                                  first_time,   gather_data_row)

    ! convert x_u(ilocal,jlocal) to zu_vsum(nx-1,ny-1)
    zu_vsum(:,:) = 0.0d0
    do j = 1, jlocal
       jj = j + staggered_jlo - 1
       do i = 1, ilocal
          ii = i + staggered_ilo - 1
          zu_vsum(ii,jj) = x_u(i,j)
       enddo
    enddo

    ! convert rv_vsum(nx-1,ny-1) to b_v(jlocal,ilocal)
    do i = 1, ilocal
       ii = i + staggered_ilo - 1
       do j = 1, jlocal
          jj = j + staggered_jlo - 1
          b_v(j,i) = rv_vsum(ii,jj)
       enddo
    enddo

    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'Before global tridiag PC v solve, r, j =', rtest, jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, Adiag_v, Asubdiag_v, Asupdiag_v, b_v:'
       do i = itest-3, itest+3
          write(iulog,'(i4, 4e16.8)') i, Adiag_v(i,j), Asubdiag_v(i,j), Asupdiag_v(i,j), b_v(i,j)
       enddo
    endif

    call tridiag_solver_global_2d(jlocal,       ilocal,      &
                                  parallel,     tasks_col,   &
                                  'col',                     &  ! tridiagonal solve for each column
                                  jtest - staggered_jlo + 1, &  ! jtest referenced to (jlocal,ilocal) coordinates
                                  itest - staggered_ilo + 1, &  ! itest referenced to (jlocal,ilocal) coordinates
                                  rtest,                     &
                                  Adiag_v,                   &
                                  Asubdiag_v,   Asupdiag_v,  &
                                  omega_v,      denom_v,     &
                                  xuh_v,        xlh_v,       &
                                  b_v,          x_v,         &
                                  first_time,   gather_data_col)

    ! convert x_v(jlocal,ilocal) to zv_vsum(nx-1,ny-1)
    zv_vsum(:,:) = 0.0d0
    do i = 1, ilocal
       ii = i + staggered_ilo - 1
       do j = 1, jlocal
          jj = j + staggered_jlo - 1
          zv_vsum(ii,jj) = x_v(j,i)
       enddo
    enddo

    !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
    !TODO: See whether tridiag_solver_local_2d could be modified to provide zu and zv in halo cells?
    call staggered_parallel_halo(zu_vsum, parallel)
    call staggered_parallel_halo(zv_vsum, parallel)

    ! Solve a tridiagonal problem for zu and zv in each problem, to account for vertical shear

    call easy_sia_solver(nx,   ny,   nz,        &
                         active_vertex,         &
                         Muu,  ru,   zu)      ! solve Muu*zu = ru for zu

    call easy_sia_solver(nx,   ny,   nz,        &
                         active_vertex,         &
                         Mvv,  rv,   zv)      ! solve Mvv*zv = rv for zv

    !WHL - debug
    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'i, zu_ssa, zu_sia(1), zu_sia(nz), zu_sum(1):'
       do i = itest-3, itest+3
          write(iulog,'(i4, 4e16.8)') i, zu_vsum(i,j), zu(1,i,j), zu(nz,i,j), zu_vsum(i,j) + zu(1,i,j)
       enddo  ! i
       write(iulog,*) ' '
       write(iulog,*) 'i, zv_ssa, zv_sia(1), zv_sia(nz), zv_sum(1):'
       do i = itest-3, itest+3
          write(iulog,'(i4, 4e16.8)') i, zv_vsum(i,j), zv(1,i,j), zv(nz,i,j), zv_vsum(i,j) + zv(1,i,j)
       enddo
    endif

    ! Combine the SSA-type solution (zu_vsum and zv_vsum) with the SIA solution
    do j = 1, ny-1
       do i = 1, nx-1
          zu(:,i,j) = zu(:,i,j) + zu_vsum(i,j)
          zv(:,i,j) = zv(:,i,j) + zv_vsum(i,j)
       enddo
    enddo

    !WHL - debug
    if (verbose_pcg .and. this_rank == rtest) write(iulog,*) 'Done in tridiag_solver_global_3d'

  end subroutine tridiag_solver_global_3d

!****************************************************************************

  subroutine matvec_multiply_structured_3d(nx,        ny,              &
                                           nz,        parallel,        &
                                           indxA_3d,  active_vertex,   &
                                           Auu,       Auv,             &
                                           Avu,       Avv,             &
                                           xu,        xv,              &
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
         nz                     ! number of vertical layers at which velocity is computed

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
         indxA_3d               ! maps relative (x,y,z) coordinates to an index between 1 and 27
    
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

    integer :: &
       staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
       staggered_jlo, staggered_jhi

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

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

                   m = indxA_3d(iA,jA,kA)

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

end module glissade_velo_higher_pcg_3d

!****************************************************************************
