!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_higher_pcg_2d.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains the 2D solvers called from glissade_velo_higher.F90 to solve
! the problem Ax = b using the preconditioned conjugate gradient method, together with
! their preconditioners and matrix-vector multiply. The tridiagonal preconditioner
! setups and solvers are public, since the 3D solvers in glissade_velo_higher_pcg_3d
! also use them.
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

  module glissade_velo_higher_pcg_2d

    use glimmer_global, only: dp
    use glimmer_paramets, only: iulog
    use glissade_types   ! for preconditioning options
    use glimmer_log
    use cism_profile, only: t_startf, t_stopf
    use cism_parallel, only: this_rank, main_task, &
         parallel_type, staggered_parallel_halo, parallel_global_sum_stagger

    implicit none
    
    private
    public :: pcg_solver_standard_2d, pcg_solver_chrongear_2d, &
              setup_preconditioner_tridiag_local_2d, setup_preconditioner_tridiag_global_2d, &
              tridiag_solver_local_2d, tridiag_solver_global_2d

    logical :: verbose_pcg = .false.
    logical :: verbose_tridiag = .false.

  contains

!****************************************************************************

  subroutine pcg_solver_standard_2d(nx,        ny,            &
                                    parallel,                 &
                                    indxA_2d,  active_vertex, &
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
       nx, ny                 ! horizontal grid dimensions (for scalars)
                              ! velocity grid has dimensions (nx-1,ny-1)

    type(parallel_type), intent(in) :: &
         parallel             ! info for parallel communication

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for vertices (i,j) where velocity is computed, else F
 
    real(dp), dimension(nx-1,ny-1,9), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 3rd dimension = 9 (node and its nearest neighbors in x and y direction)
                              ! 1st and 2nd dimensions = (x,y) indices
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

    integer ::  i, j, ii, jj ! grid indices
    integer ::  iter         ! iteration counter

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

    ! tridiagonal matrix elements
    real(dp), dimension(:,:), allocatable :: &
         Asubdiag_u, Adiag_u, Asupdiag_u,  &  ! matrix entries from Auu for tridiagonal preconditioning
         Asubdiag_v, Adiag_v, Asupdiag_v      ! matrix entries from Avv for tridiagonal preconditioning
    
    real(dp), dimension(:,:), allocatable :: &
         omega_u, omega_v,         & ! work arrays for tridiagonal solve
         denom_u, denom_v,         &
         xuh_u,   xuh_v,           &
         xlh_u,   xlh_v

    real(dp), dimension(:,:), allocatable :: &
         b_u, b_v, x_u, x_v

    real(dp), dimension(:,:), allocatable :: &
         gather_data_row,   &   ! arrays for gathering data from every task on a row or column
         gather_data_col

    integer :: ilocal, jlocal   ! number of locally owned vertices in each direction   

    integer :: &
         tasks_row,         &   ! number of tasks per row and column for tridiagonal solve
         tasks_col

    logical :: first_time       ! true on the first subroutine call (iter = 1), false thereafter

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    !TODO - Apply the following for tridiag PCs?
!    integer, parameter :: &
!         maxiters_tridiag = 100 ! max number of linear iterations for tridiagonal preconditioning,
!                                ! which generally leads to faster convergence than diagonal preconditioning 

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

       if (verbose_pcg .and. main_task) then
          write(iulog,*) 'Using no preconditioner'
       endif

    elseif (precond == PRECONDITIONER_DIAGONAL) then

       if (verbose_pcg .and. main_task) then
          write(iulog,*) 'Using diagonal matrix for preconditioning'
       endif  ! verbose_pcg

       call setup_preconditioner_diag_2d(nx,         ny,       &
                                         indxA_2d,             &
                                         Auu,        Avv,      &
                                         Adiagu,     Adiagv)

    elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then

       !WHL - debug
       if (verbose_tridiag .and. this_rank==rtest) then
          i = itest
          j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'r, i, j =', this_rank, i, j
          write(iulog,*) 'Auu =', Auu(i,j,:)
          write(iulog,*) 'Avv =', Avv(i,j,:)
       endif

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

       call setup_preconditioner_tridiag_local_2d(&
            nx,            ny,             &
            parallel,      indxA_2d,       &
            itest,         jtest,  rtest,  &
            Auu,           Avv,            &
            Adiag_u,       Adiag_v,        &
            Asubdiag_u,    Asubdiag_v,     &
            Asupdiag_u,    Asupdiag_v,     &
            omega_u,       omega_v,        &
            denom_u,       denom_v)

    elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then

       ! Allocate tridiagonal matrices
       ! Note: (i,j) indices are switced for the A_v matrices to reduce striding.

       allocate(Adiag_u   (ilocal,jlocal))
       allocate(Asubdiag_u(ilocal,jlocal))
       allocate(Asupdiag_u(ilocal,jlocal))
       allocate(omega_u(ilocal,jlocal))
       allocate(denom_u(ilocal,jlocal))
       allocate(xuh_u(ilocal,jlocal))
       allocate(xlh_u(ilocal,jlocal))
       allocate(b_u(ilocal,jlocal))
       allocate(x_u(ilocal,jlocal))

       allocate(Adiag_v   (jlocal,ilocal))
       allocate(Asubdiag_v(jlocal,ilocal))
       allocate(Asupdiag_v(jlocal,ilocal))
       allocate(omega_v(jlocal,ilocal))
       allocate(denom_v(jlocal,ilocal))
       allocate(xuh_v(jlocal,ilocal))
       allocate(xlh_v(jlocal,ilocal))
       allocate(b_v(jlocal,ilocal))
       allocate(x_v(jlocal,ilocal))

       ! These two matrices are for gathering data from all tasks on a given row or column.
       allocate(gather_data_row(8*tasks_row,jlocal))
       allocate(gather_data_col(8*tasks_col,ilocal))
       gather_data_row = 0.0d0
       gather_data_col = 0.0d0

       ! Compute the entries of the tridiagonal matrices

       ! Extract tridiagonal matrix entries from Auu
       do j = 1, jlocal
          jj = j + staggered_jlo - 1
          do i = 1, ilocal
             ii = i + staggered_ilo - 1
             Asubdiag_u(i,j) = Auu(ii,jj,indxA_2d(-1,0))   ! subdiagonal elements
             Adiag_u   (i,j) = Auu(ii,jj,indxA_2d( 0,0))   ! diagonal elements
             Asupdiag_u(i,j) = Auu(ii,jj,indxA_2d( 1,0))   ! superdiagonal elements
          enddo
       enddo

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

       ! Extract tridiagonal matrix entries from Avv

       do i = 1, ilocal
          ii = i + staggered_ilo - 1
          do j = 1, jlocal
             jj = j + staggered_jlo - 1
             Asubdiag_v(j,i) = Avv(ii,jj,indxA_2d(0,-1))   ! subdiagonal elements
             Adiag_v   (j,i) = Avv(ii,jj,indxA_2d(0, 0))   ! diagonal elements
             Asupdiag_v(j,i) = Avv(ii,jj,indxA_2d(0, 1))   ! superdiagonal elements
          enddo
       enddo

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

    endif   ! precond

    if (verbose_pcg .and. main_task) write(iulog,*) 'Done in PC setup'

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
    call matvec_multiply_structured_2d(nx,            ny,            &
                                       parallel,                     &
                                       indxA_2d,      active_vertex, &
                                       Auu,           Auv,           &
                                       Avu,           Avv,           &
                                       xu,            xv,            &
                                       zu,            zv)
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
    L2_rhs = parallel_global_sum_stagger(work0u, parallel, arr2=work0v)
    call t_stopf("pcg_glbsum_init")

    ! take square root

    L2_rhs = sqrt(L2_rhs)       ! L2 norm of RHS

    ! iterate to solution

    iter_loop: do iter = 1, maxiters

       if (verbose_pcg .and. main_task) then
!          write(iulog,*) 'iter =', iter
       endif

       call t_startf("pcg_precond")

       ! Compute PC(r) = solution z of Mz = r

       if (precond == PRECONDITIONER_NONE) then      ! no preconditioning

           zu(:,:) = ru(:,:)         ! PC(r) = r     
           zv(:,:) = rv(:,:)         ! PC(r) = r    

       elseif (precond == PRECONDITIONER_DIAGONAL) then  ! diagonal preconditioning

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

       elseif(precond == PRECONDITIONER_TRIDIAG_LOCAL) then  ! local

          if (verbose_tridiag .and. this_rank == rtest) then
             i = itest
             j = jtest
             write(iulog,*) 'Residual:'
             write(iulog,*) 'r, i, j, ru:', this_rank, i, j, ru(i,j)
             write(iulog,*) 'r, i, j, rv:', this_rank, i, j, rv(i,j)
             write(iulog,*) ' '
             write(iulog,*) 'jtest =', jtest
             write(iulog,*) 'i, ru, rv:'
             do i = itest-3, itest+3
                write(iulog,'(i4, 2f15.10)') i, ru(i,j), rv(i,j)
             enddo
          endif

          if (verbose_pcg .and. main_task) then
             write(iulog,*) 'call tridiag_solver_local_2d'
          endif

          ! Solve M*z = r, where M is a local tridiagonal matrix (one matrix per task)
          !TODO - Test a local solver that can compute zu and zv in the halo
          !       (to avoid the halo update below)

          call tridiag_solver_local_2d(&
               nx,            ny,            &
               parallel,                     &
               itest, jtest,  rtest,         &
               Adiag_u,       Adiag_v,       &  ! entries of preconditioning matrix
               Asubdiag_u,    Asubdiag_v,    &
               Asupdiag_u,    Asupdiag_v,    &
               omega_u,       omega_v,       &
               denom_u,       denom_v,       &
               ru,            rv,            &  ! right hand side
               zu,            zv)               ! solution
       
          !WHL - debug
          if (verbose_pcg .and. this_rank == rtest) then
             j = jtest
             write(iulog,*) ' '
             write(iulog,*) 'tridiag solve: i, ru, rv, zu, zv:'
             do i = itest-3, itest+3
                write(iulog,'(i4, 4f16.10)') i, ru(i,j), rv(i,j), zu(i,j), zv(i,j)
             enddo
          endif

          !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
          !TODO: See whether tridiag solvers could be modified to provide zu and zv in halo cells?
          call staggered_parallel_halo(zu, parallel)
          call staggered_parallel_halo(zv, parallel)

       elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then   ! tridiagonal preconditioning with global solve

          ! Use a global tridiagonal solver to find an approximate solution of A*z = r
          if (iter == 1) then
             first_time = .true.
          else
             first_time = .false.
          endif

          ! convert ru(nx-1,ny-1) to b_u(ilocal,jlocal)
          do j = 1, jlocal
             jj = j + staggered_jlo - 1
             do i = 1, ilocal
                ii = i + staggered_ilo - 1
                b_u(i,j) = ru(ii,jj)
             enddo
          enddo

          ! Solve M*z = r, where M is a global tridiagonal matrix

          call tridiag_solver_global_2d(&
               ilocal,       jlocal,      &
               parallel,     tasks_row,   &
               'row',                     &  ! tridiagonal solve for each row
!!                                          itest, jtest, rtest,      &
               itest - staggered_ilo + 1, &  ! itest referenced to (ilocal,jlocal) coordinates                                  
               jtest - staggered_jlo + 1, &  ! jtest referenced to (ilocal,jlocal) coordinates
               rtest,                     &
               Adiag_u,                   &
               Asubdiag_u,   Asupdiag_u,  &
               omega_u,      denom_u,     &
               xuh_u,        xlh_u,       &
               b_u,          x_u,         &
               first_time,                &
               gather_data_row)

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

          call tridiag_solver_global_2d(&
               jlocal,       ilocal,      &
               parallel,     tasks_col,   &
               'col',                     &  ! tridiagonal solve for each column
               !!                                          itest, jtest, rtest,      &
               jtest - staggered_jlo + 1, &  ! jtest referenced to (jlocal,ilocal) coordinates
               itest - staggered_ilo + 1, &  ! itest referenced to (jlocal,ilocal) coordinates
               rtest,                     &
               Adiag_v,                   &
               Asubdiag_v,   Asupdiag_v,  &
               omega_v,      denom_v,     &
               xuh_v,        xlh_v,       &
               b_v,          x_v,         &
               first_time,                &
               gather_data_col)

          ! convert x_v(jlocal,ilocal) to zv(nx-1,ny-1)
          zv(:,:) = 0.0d0
          do i = 1, ilocal
             ii = i + staggered_ilo - 1
             do j = 1, jlocal
                jj = j + staggered_jlo - 1
                zv(ii,jj) = x_v(j,i)
             enddo
          enddo

          !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
          !TODO: See whether tridiag_solver_local_2d could be modified to provide zu and zv in halo cells?
          call staggered_parallel_halo(zu, parallel)
          call staggered_parallel_halo(zv, parallel)
          
       endif    ! precond

       call t_stopf("pcg_precond")

       ! Compute the dot product eta1 = (r, PC(r))

       call t_startf("pcg_dotprod")
       work0u(:,:) = ru(:,:)*zu(:,:)    ! terms of dot product (r, PC(r))
       work0v(:,:) = rv(:,:)*zv(:,:)    
       call t_stopf("pcg_dotprod")

       call t_startf("pcg_glbsum_iter")
       eta1 = parallel_global_sum_stagger(work0u, parallel, arr2=work0v)
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
       call staggered_parallel_halo(du, parallel)
       call staggered_parallel_halo(dv, parallel)
       call t_stopf("pcg_halo_iter")
  
       ! Compute q = A*d
       ! This is the one matvec multiply required for each iteration

       call t_startf("pcg_matmult_iter")
       call matvec_multiply_structured_2d(nx,            ny,            &
                                          parallel,                     &
                                          indxA_2d,      active_vertex, &
                                          Auu,           Auv,           &
                                          Avu,           Avv,           &
                                          du,            dv,            &
                                          qu,            qv)
       call t_stopf("pcg_matmult_iter")

       ! Copy old eta1 = (r, PC(r)) to eta0

       eta0 = eta1               ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i)) 

       ! Compute the dot product eta2 = (d, A*d)

       call t_startf("pcg_dotprod")
       work0u(:,:) = du(:,:) * qu(:,:)       ! terms of dot product (d, Ad)
       work0v(:,:) = dv(:,:) * qv(:,:)
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
       xu(:,:) = xu(:,:) + alpha * du(:,:)    ! new solution, x_(i+1) = x_i + alpha*d
       xv(:,:) = xv(:,:) + alpha * dv(:,:)

       ru(:,:) = ru(:,:) - alpha * qu(:,:)    ! new residual, r_(i+1) = r_i - alpha*(Ad)
       rv(:,:) = rv(:,:) - alpha * qv(:,:)
       call t_stopf("pcg_vecupdate")

       ! Check for convergence every linear_solve_ncheck iterations.
       ! Also check at iter = 5, to reduce iterations when the nonlinear solver is close to convergence.
       ! TODO: Check at iter = linear_solve_ncheck/2 instead of 5?  This would be answer-changing.
       !
       ! For convergence check, use r = b - Ax

       if (mod(iter, linear_solve_ncheck) == 0 .or. iter == 5) then

          ! Halo update for x

          call t_startf("pcg_halo_resid")
          call staggered_parallel_halo(xu, parallel)
          call staggered_parallel_halo(xv, parallel)
          call t_stopf("pcg_halo_resid")

          ! Compute A*x (use z as a temp vector for A*x)
           
          call t_startf("pcg_matmult_resid")
          call matvec_multiply_structured_2d(nx,            ny,            &
                                             parallel,                     &
                                             indxA_2d,      active_vertex, &
                                             Auu,           Auv,           &
                                             Avu,           Avv,           &
                                             xu,            xv,            &
                                             zu,            zv)
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
          L2_resid = parallel_global_sum_stagger(work0u, parallel, arr2=work0v)
          call t_stopf("pcg_glbsum_resid")

          ! take square root
          L2_resid = sqrt(L2_resid)       ! L2 norm of residual

          ! compute normalized error
          err = L2_resid/L2_rhs

          !WHL - debug
          if (verbose_pcg .and. this_rank == rtest) then
             ru_max = 0.d0
             rv_max = 0.d0
             iu_max = 0
             ju_max = 0
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (abs(ru(i,j)) > ru_max) then
                      ru_max = ru(i,j)
                      iu_max = i
                      ju_max = j
                   endif
                   if (abs(rv(i,j)) > rv_max) then
                      rv_max = rv(i,j)
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
    if (allocated(b_u))        deallocate(b_u, b_v)
    if (allocated(x_u))        deallocate(x_u, x_v)
    if (allocated(gather_data_row)) deallocate(gather_data_row)
    if (allocated(gather_data_col)) deallocate(gather_data_col)

  end subroutine pcg_solver_standard_2d

!****************************************************************************

  subroutine pcg_solver_chrongear_2d(nx,        ny,            &
                                     parallel,                 &
                                     indxA_2d,  active_vertex, &
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
    !  The current preconditioning options for the solver are
    !  (0) no preconditioning
    !  (1) diagonal preconditioning
    !  (3) local tridiagonal preconditioning
    !  (4) global tridiagonal preconditioning
    !  The SIA-based preconditioning option is not available for a 2D solve.
    !
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                 ! horizontal grid dimensions (for scalars)
                              ! velocity grid has dimensions (nx-1,ny-1)

    type(parallel_type), intent(in) :: &
       parallel               ! info for parallel communication

    integer, dimension(-1:1,-1:1), intent(in) :: &
       indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F
 
    real(dp), dimension(nx-1,ny-1,9), intent(in) ::   &
       Auu, Auv, &            ! four components of assembled matrix
       Avu, Avv               ! 3rd dimension = 9 (node and its nearest neighbors in x and y direction)
                              ! 1st and 2nd dimensions = (x,y) indices
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
       itest, jtest, rtest      ! point for debugging diagnostics

    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer ::  i, j, ii, jj ! grid indices
    integer ::  m            ! matrix element index
    integer ::  iter         ! iteration counter

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
         xlh_u,   xlh_v

    real(dp), dimension(:,:), allocatable :: &
         b_u, b_v, x_u, x_v

    real(dp), dimension(:,:), allocatable :: &
         gather_data_row,   &   ! arrays for gathering data from every task on a row or column
         gather_data_col

    integer :: ilocal, jlocal   ! number of locally owned vertices in each direction

    integer :: maxiters_chrongear  ! max number of linear iterations before quitting

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    integer :: &
         tasks_row,         &   ! number of tasks per row and column for tridiagonal solve
         tasks_col

    integer, parameter :: &
       maxiters_tridiag = 100  ! max number of linear iterations for tridiagonal preconditioning,
                               ! which generally leads to faster convergence than diagonal preconditioning

    !WHL - debug
    real(dp) :: usum, usum_global, vsum, vsum_global

    !WHL - debug
    integer :: iu_max, ju_max, iv_max, jv_max
    real(dp) :: ru_max, rv_max
    real(dp) :: sum_temp

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    tasks_row = parallel%tasks_row
    tasks_col = parallel%tasks_col

    ! Set the maximum number of linear iterations.
    ! Typically allow up to 200 iterations with diagonal preconditioning, but only 100
    !  with tridiagonal, which usually converges faster.

    if (precond == PRECONDITIONER_TRIDIAG_LOCAL .or. precond == PRECONDITIONER_TRIDIAG_GLOBAL) then
       maxiters_chrongear = maxiters_tridiag
    else
       maxiters_chrongear = maxiters
    endif

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Using native PCG solver (Chronopoulos-Gear)'
       write(iulog,*) 'tolerance, maxiters, precond =', tolerance, maxiters_chrongear, precond
    endif

    ! Compute array sizes for locally owned vertices
    ilocal = staggered_ihi - staggered_ilo + 1
    jlocal = staggered_jhi - staggered_jlo + 1

    !---- Set up matrices for preconditioning

    call t_startf("pcg_precond_init")

    if (precond == PRECONDITIONER_NONE) then    ! no preconditioner

       if (verbose_pcg .and. main_task) then
          write(iulog,*) 'Using no preconditioner'
       endif

    elseif (precond == PRECONDITIONER_DIAGONAL) then

       call setup_preconditioner_diag_2d(nx,       ny,         &
                                         indxA_2d,             &
                                         Auu,      Avv,        &
                                         Adiagu,   Adiagv)

       !WHL - debug
       if (verbose_pcg .and. this_rank == rtest) then
          i = itest
          j = jtest
          write(iulog,*) 'i, j, r =', i, j, this_rank
          write(iulog,*) 'Au diag =', Adiagu(i,j)
          write(iulog,*) 'Av diag =', Adiagv(i,j)
       endif

    elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then

       !WHL - debug
       if (verbose_tridiag .and. this_rank==rtest) then
          i = itest
          j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'r, i, j =', this_rank, i, j
          write(iulog,*) 'Auu =', Auu(i,j,:)
          write(iulog,*) 'Avv =', Avv(i,j,:)
       endif

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

       call setup_preconditioner_tridiag_local_2d(&
            nx,            ny,             &
            parallel,      indxA_2d,       &
            itest,         jtest,  rtest,  &
            Auu,           Avv,            &
            Adiag_u,       Adiag_v,        &
            Asubdiag_u,    Asubdiag_v,     &
            Asupdiag_u,    Asupdiag_v,     &
            omega_u,       omega_v,        &
            denom_u,       denom_v)

    elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then

       ! Allocate tridiagonal preconditioning matrices
       ! Note: (i,j) indices are switced for the A_v matrices to reduce striding.

       allocate(Adiag_u   (ilocal,jlocal))
       allocate(Asubdiag_u(ilocal,jlocal))
       allocate(Asupdiag_u(ilocal,jlocal))
       allocate(omega_u(ilocal,jlocal))
       allocate(denom_u(ilocal,jlocal))
       allocate(xuh_u(ilocal,jlocal))
       allocate(xlh_u(ilocal,jlocal))
       allocate(b_u(ilocal,jlocal))
       allocate(x_u(ilocal,jlocal))

       allocate(Adiag_v   (jlocal,ilocal))
       allocate(Asubdiag_v(jlocal,ilocal))
       allocate(Asupdiag_v(jlocal,ilocal))
       allocate(omega_v(jlocal,ilocal))
       allocate(denom_v(jlocal,ilocal))
       allocate(xuh_v(jlocal,ilocal))
       allocate(xlh_v(jlocal,ilocal))
       allocate(b_v(jlocal,ilocal))
       allocate(x_v(jlocal,ilocal))

       ! These two matrices are for gathering data from all tasks on a given row or column.
       allocate(gather_data_row(8*tasks_row,jlocal))
       allocate(gather_data_col(8*tasks_col,ilocal))
       gather_data_row = 0.0d0
       gather_data_col = 0.0d0

       ! Compute the entries of the tridiagonal matrices

       ! Extract tridiagonal matrix entries from Auu
       do j = 1, jlocal
          jj = j + staggered_jlo - 1
          do i = 1, ilocal
             ii = i + staggered_ilo - 1
             Asubdiag_u(i,j) = Auu(ii,jj,indxA_2d(-1,0))   ! subdiagonal elements
             Adiag_u   (i,j) = Auu(ii,jj,indxA_2d( 0,0))   ! diagonal elements
             Asupdiag_u(i,j) = Auu(ii,jj,indxA_2d( 1,0))   ! superdiagonal elements
          enddo
       enddo

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

       ! Extract tridiagonal matrix entries from Avv
       do i = 1, ilocal
          ii = i + staggered_ilo - 1
          do j = 1, jlocal
             jj = j + staggered_jlo - 1
             Asubdiag_v(j,i) = Avv(ii,jj,indxA_2d(0,-1))   ! subdiagonal elements
             Adiag_v   (j,i) = Avv(ii,jj,indxA_2d(0, 0))   ! diagonal elements
             Asupdiag_v(j,i) = Avv(ii,jj,indxA_2d(0, 1))   ! superdiagonal elements
          enddo
       enddo

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

    endif   ! precond

    !WHL - debug
    if (verbose_pcg .and. main_task) write(iulog,*) 'Done in PC setup'

    call t_stopf("pcg_precond_init")

    !---- Initialize scalars and vectors

    niters = maxiters_chrongear
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
    bb = parallel_global_sum_stagger(worku, parallel, arr2=workv)
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
    call staggered_parallel_halo(xu, parallel)
    call staggered_parallel_halo(xv, parallel)
    call t_stopf("pcg_halo_init")

    !---- Compute A*x   (use z as a temp vector for A*x)

    call t_startf("pcg_matmult_init")
    call matvec_multiply_structured_2d(nx,            ny,            &
                                       parallel,                     &
                                       indxA_2d,      active_vertex, &
                                       Auu,           Auv,           &
                                       Avu,           Avv,           &
                                       xu,            xv,            &
                                       zu,            zv)
    call t_stopf("pcg_matmult_init")

    !---- Compute the initial residual r = b - A*x
    !---- This is correct for locally owned nodes.

    call t_startf("pcg_vecupdate")
    ru(:,:) = bu(:,:) - zu(:,:)
    rv(:,:) = bv(:,:) - zv(:,:)
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

    if (precond == PRECONDITIONER_NONE) then      ! no preconditioning

       zu(:,:) = ru(:,:)         ! PC(r) = r     
       zv(:,:) = rv(:,:)         ! PC(r) = r    

    elseif (precond == PRECONDITIONER_DIAGONAL) then  ! diagonal preconditioning

       ! Solve Mz = r, where M is a diagonal matrix
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
!          i = itest
!          write(iulog,*) ' '
!          write(iulog,*) 'zv solve with diagonal precond, this_rank, i =', this_rank, i
!          write(iulog,*) 'j, active, Adiagv, rv, zv, xv:'
!          do j = staggered_jhi, staggered_jlo, -1
!             write(iulog,'(i4, l4, 2f12.3, e12.3, f12.3)') j, active_vertex(i,j), Adiagv(i,j), rv(i,j), zv(i,j), xv(i,j)
!          enddo
       endif

    elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then  ! local

       if (verbose_tridiag .and. this_rank == rtest) then
          i = itest
          j = jtest
          write(iulog,*) 'Residual:'
          write(iulog,*) 'r, i, j, ru:', this_rank, i, j, ru(i,j)
          write(iulog,*) 'r, i, j, rv:', this_rank, i, j, rv(i,j)
          write(iulog,*) ' '
          write(iulog,*) 'jtest =', jtest
          write(iulog,*) 'i, ru, rv:'
          do i = itest-3, itest+3
             write(iulog,'(i4, 2f15.10)') i, ru(i,j), rv(i,j)
          enddo
       endif

       if (verbose_pcg .and. main_task) then
          write(iulog,*) 'call tridiag_solver_local_2d'
       endif

       ! Solve M*z = r, where M is a local tridiagonal matrix (one matrix per task)

       !TODO - Test a local solver that can compute zu and zv in the halo
       !       (to avoid the halo update below)

       call tridiag_solver_local_2d(nx,            ny,            &
                                    parallel,                     &
                                    itest, jtest,  rtest,         &
                                    Adiag_u,       Adiag_v,       &  ! entries of preconditioning matrix
                                    Asubdiag_u,    Asubdiag_v,    &
                                    Asupdiag_u,    Asupdiag_v,    &
                                    omega_u,       omega_v,       &
                                    denom_u,       denom_v,       &
                                    ru,            rv,            &  ! right hand side
                                    zu,            zv)               ! solution

       !WHL - debug
       if (verbose_pcg .and. this_rank == rtest) then
          j = jtest
          write(iulog,*) ' '
          write(iulog,*) 'tridiag solve: i, ru, rv, zu, zv:'
          do i = itest-3, itest+3
             write(iulog,'(i4, 4f16.10)') i, ru(i,j), rv(i,j), zu(i,j), zv(i,j)
          enddo
       endif

       !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
       !TODO: See whether tridiag solvers could be modified to provide zu and zv in halo cells?
       call staggered_parallel_halo(zu, parallel)
       call staggered_parallel_halo(zv, parallel)

    elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then   ! tridiagonal preconditioning with global solve

       ! convert ru(nx-1,ny-1) to b_u(ilocal,jlocal)
       do j = 1, jlocal
          jj = j + staggered_jlo - 1
          do i = 1, ilocal
             ii = i + staggered_ilo - 1
             b_u(i,j) = ru(ii,jj)
          enddo
       enddo

       ! Solve M*z = r, where M is a global tridiagonal matrix

       call tridiag_solver_global_2d(ilocal,       jlocal,      &
                                     parallel,     tasks_row,   &
                                     'row',                     &  ! tridiagonal solve for each row
!!                                          itest, jtest, rtest,      &
                                     itest - staggered_ilo + 1, &  ! itest referenced to (ilocal,jlocal) coordinates
                                     jtest - staggered_jlo + 1, &  ! jtest referenced to (ilocal,jlocal) coordinates
                                     rtest,                     &
                                     Adiag_u,                   &
                                     Asubdiag_u,   Asupdiag_u,  &
                                     omega_u,      denom_u,     &
                                     xuh_u,        xlh_u,       &
                                     b_u,          x_u,         &
                                     .true.,                    &  ! first_time
                                     gather_data_row)

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

       call tridiag_solver_global_2d(jlocal,       ilocal,      &
                                     parallel,     tasks_col,   &
                                     'col',                     &  ! tridiagonal solve for each column
!!                                          itest, jtest, rtest,      &
                                     jtest - staggered_jlo + 1, &  ! jtest referenced to (jlocal,ilocal) coordinates
                                     itest - staggered_ilo + 1, &  ! itest referenced to (jlocal,ilocal) coordinates
                                     rtest,                     &
                                     Adiag_v,                   &
                                     Asubdiag_v,   Asupdiag_v,  &
                                     omega_v,      denom_v,     &
                                     xuh_v,        xlh_v,       &
                                     b_v,          x_v,         &
                                     .true.,                    &  ! first_time
                                     gather_data_col)

       ! convert x_v(jlocal,ilocal) to zv(nx-1,ny-1)

       zv(:,:) = 0.0d0
       do i = 1, ilocal
          ii = i + staggered_ilo - 1
          do j = 1, jlocal
             jj = j + staggered_jlo - 1
             zv(ii,jj) = x_v(j,i)
          enddo
       enddo

       !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
       !TODO: See whether tridiag_solver_local_2d could be modified to provide zu and zv in halo cells?
       call staggered_parallel_halo(zu, parallel)
       call staggered_parallel_halo(zv, parallel)

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
    call matvec_multiply_structured_2d(nx,            ny,            &
                                       parallel,                     &
                                       indxA_2d,      active_vertex, &
                                       Auu,           Auv,           &
                                       Avu,           Avv,           &
                                       du,            dv,            &
                                       qu,            qv)
    call t_stopf("pcg_matmult_iter")

    !---- Compute intermediate result for dot product (d,q) = (d,Ad)

    call t_startf("pcg_dotprod")
    work2u(:,:,2) = du(:,:) * qu(:,:)
    work2v(:,:,2) = dv(:,:) * qv(:,:)
    call t_stopf("pcg_dotprod")

    !---- Find global sums of (r,z) and (d,q)

    call t_startf("pcg_glbsum_iter")
    gsum = parallel_global_sum_stagger(work2u, 2, parallel, arr2=work2v)  ! nflds = 2
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
    xu(:,:) = xu(:,:) + alpha*du(:,:)
    xv(:,:) = xv(:,:) + alpha*dv(:,:)

    ru(:,:) = ru(:,:) - alpha*qu(:,:)     ! q = A*d
    rv(:,:) = rv(:,:) - alpha*qv(:,:)
    call t_stopf("pcg_vecupdate")

    !WHL - debug
    if (verbose_pcg .and. this_rank == rtest) then
       j = jtest
!!       write(iulog,*) ' '
!!       write(iulog,*) 'iter = 1: i, xu, xv, ru, rv:'
       do i = itest-3, itest+3
!!          write(iulog,'(i4, 4f16.10)') i, xu(i,j), xv(i,j), ru(i,j), rv(i,j)
       enddo  ! i
    endif


    !---------------------------------------------------------------
    ! Iterate to solution
    !---------------------------------------------------------------

    iter_loop: do iter = 2, maxiters_chrongear  ! first iteration done above

       if (verbose_pcg .and. main_task) then
!          write(iulog,*) 'iter =', iter
       endif

       !---- Compute PC(r) = solution z of Mz = r
       !---- z is correct in halo

       call t_startf("pcg_precond_iter")

       if (precond == PRECONDITIONER_NONE) then      ! no preconditioning

           zu(:,:) = ru(:,:)         ! PC(r) = r
           zv(:,:) = rv(:,:)         ! PC(r) = r    

       elseif (precond == PRECONDITIONER_DIAGONAL) then  ! diagonal preconditioning

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

       elseif (precond == PRECONDITIONER_TRIDIAG_LOCAL) then   ! tridiagonal preconditioning with local solve

          ! Solve M*z = r, where M is a local tridiagonal matrix (one matrix per task)

          !TODO - Test a local solver that can compute zu and zv in the halo
          !       (to avoid the halo update below)

          !WHL - debug
          if (verbose_tridiag .and. this_rank == rtest) then
             i = itest
             j = jtest
!             write(iulog,*) 'Residual:'
!             write(iulog,*) 'r, i, j, ru:', this_rank, i, j, ru(i,j)
!             write(iulog,*) 'r, i, j, rv:', this_rank, i, j, rv(i,j)
!             write(iulog,*) ' '
!             write(iulog,*) 'jtest =', jtest
!             write(iulog,*) 'i, ru, rv:'
!             do i = staggered_ihi, staggered_ilo, -1
!                write(iulog,'(i4, 2f15.10)') i, ru(i,j), rv(i,j)
!             enddo
          endif


          !WHL - debug
          if (verbose_pcg .and. this_rank == rtest) then
             write(iulog,*) '   call tridiag_solver_local_2d'
          endif

          call tridiag_solver_local_2d(nx,            ny,            &
                                       parallel,                     &
                                       itest, jtest,  rtest,         &
                                       Adiag_u,       Adiag_v,       &  ! entries of preconditioning matrix
                                       Asubdiag_u,    Asubdiag_v,    &
                                       Asupdiag_u,    Asupdiag_v,    &
                                       omega_u,       omega_v,       &
                                       denom_u,       denom_v,       &
                                       ru,            rv,            &  ! right hand side
                                       zu,            zv)               ! solution

          !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
          !TODO: See whether tridiag solvers could be modified to provide zu and zv in halo cells?
          call staggered_parallel_halo(zu, parallel)
          call staggered_parallel_halo(zv, parallel)

          if (verbose_pcg .and. this_rank == rtest) then
             j = jtest
             write(iulog,*) ' '
             write(iulog,*) 'tridiag solve: i, ru, rv, zu, zv:'
             do i = itest-3, itest+3
                write(iulog,'(i4, 4f16.10)') i, ru(i,j), rv(i,j), zu(i,j), zv(i,j)
             enddo
          endif

       elseif (precond == PRECONDITIONER_TRIDIAG_GLOBAL) then   ! tridiagonal preconditioning with global solve

          !WHL - debug
          if (verbose_tridiag .and. this_rank == rtest) then
             j = jtest
             write(iulog,*) ' '
             write(iulog,*) 'jtest =', jtest
             write(iulog,*) 'i, ru, rv:'
             do i = itest-3, itest+3
                write(iulog,'(i4, 2f15.10)') i, ru(i,j), rv(i,j)
             enddo
          endif

          ! convert ru(nx-1,ny-1) to b_u(ilocal,jlocal)

          do j = 1, jlocal
             jj = j + staggered_jlo - 1
             do i = 1, ilocal
                ii = i + staggered_ilo - 1
                b_u(i,j) = ru(ii,jj)
             enddo
          enddo

          !WHL - debug
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
!!                                             itest, jtest, rtest,      &
                                        itest - staggered_ilo + 1, &  ! itest referenced to (ilocal,jlocal) coordinates
                                        jtest - staggered_jlo + 1, &  ! jtest referenced to (ilocal,jlocal) coordinates
                                        rtest,                     &
                                        Adiag_u,                   &
                                        Asubdiag_u,   Asupdiag_u,  &
                                        omega_u,      denom_u,     &
                                        xuh_u,        xlh_u,       &
                                        b_u,          x_u,         &
                                        .false.,                   &  ! first_time
                                        gather_data_row)

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

          call tridiag_solver_global_2d(jlocal,       ilocal,      &
                                        parallel,     tasks_col,   &
                                        'col',                     &  ! tridiagonal solve for each column
!!                                             itest, jtest, rtest,      &
                                        jtest - staggered_jlo + 1, &  ! jtest referenced to (jlocal,ilocal) coordinates
                                        itest - staggered_ilo + 1, &  ! itest referenced to (jlocal,ilocal) coordinates
                                        rtest,                     &
                                        Adiag_v,                   &
                                        Asubdiag_v,   Asupdiag_v,  &
                                        omega_v,      denom_v,     &
                                        xuh_v,        xlh_v,       &
                                        b_v,          x_v,         &
                                        .false.,                   &  ! first_time
                                        gather_data_col)

          ! convert x_v(jlocal,ilocal) to zv(nx-1,ny-1)

          zv(:,:) = 0.0d0
          do i = 1, ilocal
             ii = i + staggered_ilo - 1
             do j = 1, jlocal
                jj = j + staggered_jlo - 1
                zv(ii,jj) = x_v(j,i)
             enddo
          enddo

          !Note: Need zu and zv in a row of halo cells so that q = A*d is correct in all locally owned cells
          !TODO: See whether tridiag solvers could be modified to provide zu and zv in halo cells?
          call staggered_parallel_halo(zu, parallel)
          call staggered_parallel_halo(zv, parallel)

       endif    ! precond

       call t_stopf("pcg_precond_iter")

       !---- Compute Az = A*z
       !---- This is the one matvec multiply required per iteration
       !---- Az is correct for locally owned nodes and needs a halo update (below)

       call t_startf("pcg_matmult_iter")
       call matvec_multiply_structured_2d(nx,            ny,            &
                                          parallel,                     &
                                          indxA_2d,      active_vertex, &
                                          Auu,           Auv,           &
                                          Avu,           Avv,           &
                                          zu,            zv,            &
                                          Azu,           Azv)
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
       gsum = parallel_global_sum_stagger(work2u, 2, parallel, arr2=work2v)  ! nflds = 2
       call t_stopf("pcg_glbsum_iter")

       !---- Halo update for Az

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

       ! Check for convergence every linear_solve_ncheck iterations.
       ! Also check at iter = 5, to reduce iterations when the nonlinear solver is close to convergence.
       ! TODO: Check at iter = linear_solve_ncheck/2 instead of 5?  This would be answer-changing.
       !
       ! For convergence check, use r = b - Ax

       if (mod(iter, linear_solve_ncheck) == 0 .or. iter == 5) then

          if (verbose_pcg .and. main_task) then
             write(iulog,*) ' '
             write(iulog,*) '   check convergence, iter =', iter
          endif

          !---- Compute z = A*x  (use z as a temp vector for A*x)
           
          !WHL - debug - don't think this is needed, but try just in case
!!          call staggered_parallel_halo(xu, parallel)
!!          call staggered_parallel_halo(xv, parallel)

          call t_startf("pcg_matmult_resid")
          call matvec_multiply_structured_2d(nx,            ny,            &
                                             parallel,                     &
                                             indxA_2d,      active_vertex, &
                                             Auu,           Auv,           &
                                             Avu,           Avv,           &
                                             xu,            xv,            &
                                             zu,            zv)
          call t_stopf("pcg_matmult_resid")

          !---- Compute residual r = b - A*x

          !WHL - debug - don't think this is needed, but try just in case
!!          call staggered_parallel_halo(bu, parallel)
!!          call staggered_parallel_halo(bv, parallel)
          !WHL - debug - don't think this is needed, but try just in case
!!          call staggered_parallel_halo(zu, parallel)
!!          call staggered_parallel_halo(zv, parallel)

          call t_startf("pcg_vecupdate")
          ru(:,:) = bu(:,:) - zu(:,:)
          rv(:,:) = bv(:,:) - zv(:,:)
          call t_stopf("pcg_vecupdate")

          !WHL - debug - don't think this is needed, but try just in case
!!          call staggered_parallel_halo(ru, parallel)
!!          call staggered_parallel_halo(rv, parallel)

          !---- Compute dot product (r, r)

          call t_startf("pcg_dotprod")
          worku(:,:) = ru(:,:)*ru(:,:)
          workv(:,:) = rv(:,:)*rv(:,:)
          call t_stopf("pcg_dotprod")

          call t_startf("pcg_glbsum_resid")
          rr = parallel_global_sum_stagger(worku, parallel, arr2=workv)
          call t_stopf("pcg_glbsum_resid")

          L2_resid = sqrt(rr)          ! L2 norm of residual
          err = L2_resid/L2_rhs        ! normalized error

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
                   if (abs(ru(i,j)) > ru_max) then
                      ru_max = ru(i,j)
                      iu_max = i
                      ju_max = j
                   endif
                   if (abs(rv(i,j)) > rv_max) then
                      rv_max = rv(i,j)
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
    if (allocated(b_u))        deallocate(b_u, b_v)
    if (allocated(x_u))        deallocate(x_u, x_v)
    if (allocated(gather_data_row)) deallocate(gather_data_row)
    if (allocated(gather_data_col)) deallocate(gather_data_col)

  end subroutine pcg_solver_chrongear_2d

!****************************************************************************

  subroutine setup_preconditioner_diag_2d(&
       nx,         ny,      &
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

    real(dp), dimension(nx-1,ny-1,9), intent(in) ::   &
         Auu, Avv               ! two out of the four components of assembled matrix
                                ! 3rd dimension = 9 (node and its nearest neighbors in x and y direction)
                                ! 1st and 2nd dimensions = (x,y) indices
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
       write(iulog,*) 'Using diagonal matrix for preconditioning'
    endif  ! verbose_pcg

    m = indxA_2d(0,0)
    Adiagu(:,:) = Auu(:,:,m)
    Adiagv(:,:) = Avv(:,:,m)

  end subroutine setup_preconditioner_diag_2d

!****************************************************************************

  subroutine setup_preconditioner_tridiag_local_2d(nx,            ny,            &
                                                   parallel,      indxA_2d,      &
                                                   itest, jtest,  rtest,         &
                                                   Auu,           Avv,           &
                                                   Adiag_u,       Adiag_v,       &
                                                   Asubdiag_u,    Asubdiag_v,    &
                                                   Asupdiag_u,    Asupdiag_v,    &
                                                   omega_u,       omega_v,       &
                                                   denom_u,       denom_v)

    ! Set up some arrays that are used repeatedly for tridiagonal preconditioning

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
         nx, ny                 ! horizontal grid dimensions (for scalars)
                                ! velocity grid has dimensions (nx-1,ny-1)

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    integer, intent(in) :: itest, jtest, rtest   ! coordinates of diagnostic point

    integer, dimension(-1:1,-1:1), intent(in) :: &
         indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    real(dp), dimension(nx-1,ny-1,9), intent(in) ::   &
         Auu, Avv               ! two out of the four components of assembled matrix
                                ! 3rd dimension = 9 (node and its nearest neighbors in x and y direction)
                                ! 1st and 2nd dimensions = (x,y) indices
                                !
                                !    Auu  | Auv
                                !    _____|____
                                !    Avu  | Avv
                                !         |

    real(dp), dimension(nx-1,ny-1), intent(out) ::   &
         Adiag_u,    Adiag_v,      &   ! matrix entries for diagonal preconditioning
         Asubdiag_u, Asubdiag_v,   &
         Asupdiag_u, Asupdiag_v,   &
         omega_u,    omega_v,      &   ! work arrays for diagonal preconditioning
         denom_u,    denom_v

    integer :: i, j

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    if (verbose_pcg .and. main_task) then
       write(iulog,*) 'Using local tridiagonal matrix for preconditioning'
    endif  ! verbose_pcg

    ! Compute arrays that will be used repeatedly for solving Au*xu = bu over rows of the grid
    ! Note: These arrays have nonzero values for locally owned vertices,
    !       i.e. in the range (staggered_ilo:staggered_ihi, staggered_jlo:staggered_jhi)
    ! TODO: Consider extending to halo vertices so that halo updates can be avoided later.

    ! Extract tridiagonal elements from the global matrix Auu
    Asubdiag_u(:,:) = Auu(:,:,indxA_2d(-1,0))   ! subdiagonal elements
    Adiag_u   (:,:) = Auu(:,:,indxA_2d( 0,0))   ! diagonal elements
    Asupdiag_u(:,:) = Auu(:,:,indxA_2d( 1,0))   ! superdiagonal elements

    ! Modify entries as needed so that the tridiagonal matrix is nonsingular.
    ! For inactive vertices with zero diagonal entries, set diag = 1, subdiag = supdiag = 0
    ! (so the solution is x = rhs)
    where (Adiag_u == 0.0d0)
       Adiag_u    = 1.0d0
       Asubdiag_u = 0.0d0
       Asupdiag_u = 0.0d0
    endwhere

    ! initialize work arrays
    omega_u = 0.0d0
    denom_u = 0.0d0

    ! Forward elimination
    ! Note: denom -> 1/denom to speed up future computations

    do j = staggered_jlo, staggered_jhi
       i = staggered_ilo
       omega_u(i,j) = Asupdiag_u(i,j) / Adiag_u(i,j)
       do i = staggered_ilo+1, staggered_ihi
          denom_u(i,j) = Adiag_u(i,j) - Asubdiag_u(i,j)*omega_u(i-1,j)
          if (denom_u(i,j) == 0.0d0) then
             call write_log('ERROR: divzero in setup_preconditioner_tridiag_local_2d', GM_FATAL)
          else
             denom_u(i,j) = 1.0d0 / denom_u(i,j)
          endif
          omega_u(i,j) = Asupdiag_u(i,j) * denom_u(i,j)
       enddo   ! i
    enddo   ! j

    ! Take reciprocal of Adiag, to replace division with multiplication in tridiag_solver
    ! Note: Any zero values have been set to 1 above.
    Adiag_u = 1.0d0 / Adiag_u

    ! Compute arrays that will be used repeatedly for solving Av*xv = bv over columns of the grid

    ! Extract tridiagonal elements from the global matrix Auu
    Asubdiag_v(:,:) = Avv(:,:,indxA_2d( 0,-1))   ! subdiagonal elements
    Adiag_v   (:,:) = Avv(:,:,indxA_2d( 0, 0))   ! diagonal elements
    Asupdiag_v(:,:) = Avv(:,:,indxA_2d( 0, 1))   ! superdiagonal elements

    ! Modify entries as needed so that the tridiagonal matrix is nonsingular.
    ! For inactive vertices with zero diagonal entries, set diag = 1, subdiag = supdiag = 0
    ! (so the solution is x = rhs)
    where (Adiag_v == 0.0d0)
       Adiag_v    = 1.0d0
       Asubdiag_v = 0.0d0
       Asupdiag_v = 0.0d0
    endwhere

    ! initialize work arrays
    omega_v = 0.0d0
    denom_v = 0.0d0

    ! Forward elimination
    ! Note: denom -> 1/denom to speed up future computations

    do i = staggered_ilo, staggered_ihi
       j = staggered_jlo
       omega_v(i,j) = Asupdiag_v(i,j) / Adiag_v(i,j)
       do j = staggered_jlo+1, staggered_jhi
          denom_v(i,j) = Adiag_v(i,j) - Asubdiag_v(i,j)*omega_v(i,j-1)
          if (denom_v(i,j) == 0.0d0) then
             call write_log('ERROR: divzero in setup_preconditioner_tridiag_2d', GM_FATAL)
          else
             denom_v(i,j) = 1.0d0 / denom_v(i,j)
          endif
          omega_v(i,j) = Asupdiag_v(i,j) * denom_v(i,j)
       enddo   ! j
    enddo   ! i

    ! Take reciprocal of Adiag, to replace division with multiplication in tridiag_solver
    ! Note: Any zero values have been set to 1 above.
    Adiag_v = 1.0d0 / Adiag_v

  end subroutine setup_preconditioner_tridiag_local_2d

!****************************************************************************

  subroutine setup_preconditioner_tridiag_global_2d(ilocal,       jlocal,    &
                                                    itest, jtest, rtest,     &
                                                    Adiag,                   &
                                                    Asubdiag,     Asupdiag,  &
                                                    omega,        denom,     &
                                                    xuh,          xlh)

    ! Set up tridiagonal preconditioning matrices for 2D SSA-type solve

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
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

    !WHL - debug
    if (verbose_tridiag .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'In setup_preconditioner_tridiag_global_2d: itest, jtest, rtest =', itest, jtest, rtest
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

       !WHL - debug
       if (verbose_tridiag .and. this_rank == rtest .and. j == 10) then
          write(iulog,*) ' '
          write(iulog,*) 'j =', j
          write(iulog,*) 'i, Adiag, Asupdiag, Asubdiag, omega, denom, xlh, xuh:'
          do i = 1, ilocal
             write(iulog,*) i, Adiag(i,j), Asubdiag(i,j), Asupdiag(i,j), omega(i,j), denom(i,j), xlh(i,j), xuh(i,j)
          enddo
       endif


    enddo   ! j

    ! Take reciprocal of Adiag, to replace division with multiplication below.
    ! Note: Any zero values have been set to 1 above.
    Adiag = 1.0d0/Adiag

  end subroutine setup_preconditioner_tridiag_global_2d

!****************************************************************************

  subroutine tridiag_solver_local_2d(nx,            ny,            &
                                     parallel,                     &
                                     itest, jtest,  rtest,         &
                                     Adiag_u,       Adiag_v,       &
                                     Asubdiag_u,    Asubdiag_v,    &
                                     Asupdiag_u,    Asupdiag_v,    &
                                     omega_u,       omega_v,       &
                                     denom_u,       denom_v,       &
                                     bu,            bv,            &
                                     xu,            xv)

    use glimmer_utils, only: tridiag

    integer, intent(in) :: &
         nx, ny                    ! horizontal grid dimensions (for scalars)

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    integer, intent(in) :: itest, jtest, rtest   ! diagnostic only

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
         Adiag_u,    Adiag_v,        &  ! matrix entries for tridiagonal preconditioning
         Asubdiag_u, Asubdiag_v,     &
         Asupdiag_u, Asupdiag_v,     &
         omega_u,    omega_v,        &  ! work arrays for tridiagonal preconditioning
         denom_u,    denom_v

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
         bu, bv                 ! right-hand side vectors

    real(dp), dimension(nx-1,ny-1), intent(out) :: &
         xu, xv                 ! solution vectors

    ! local variables

    real(dp), dimension(nx-1,ny-1) :: &
         gamma_u, gamma_v        ! work arrays

    integer :: i, j

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

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
    !
    !-------------------------------------------------------------------------------------

    !WHL - debug
    if (verbose_tridiag .and. this_rank == rtest) then
       write(iulog,*) ' '
       write(iulog,*) 'In tridiag_solver_local_2d: itest, jtest, rtest =', itest, jtest, rtest
    endif

    !-------------------------------------------------------------------------------------
    ! Solve a tridiagonal problem Au*xu = bu for each row j, from the bottom to the top of the global domain.
    ! The matrix is solved using the standard Thomas algorithm, with terms omega and denom
    !  precomputed in the setup subroutine (since they are independent of the rhs).
    !-------------------------------------------------------------------------------------

    ! Initialize
    xu = 0.0d0
    gamma_u = 0.0d0

    ! Loop over locally owned rows
    !TODO - Expand the loop range by 1 into halo layers?

    do j = staggered_jlo, staggered_jhi

       ! Forward elimination
       ! Note: denom -> 1/denom in setup subroutine
       !       Adiag -> 1/Adiag in setup subroutine

       i = staggered_ilo
       gamma_u(i,j) = bu(i,j) * Adiag_u(i,j)   ! actually bu/Adiag
       do i = staggered_ilo+1, staggered_ihi
          gamma_u(i,j) = (bu(i,j) - Asubdiag_u(i,j)*gamma_u(i-1,j)) * denom_u(i,j)   ! actually (bu-Asubdiag_u*gamma_u)/denom_u
       enddo

       i = staggered_ihi
       xu(i,j) = gamma_u(i,j)
       do i = staggered_ihi-1, 1, -1
          xu(i,j) = gamma_u(i,j) - omega_u(i,j)*xu(i+1,j)
       enddo

       !WHL - debug
       if (verbose_tridiag .and. this_rank == rtest .and. j == jtest) then
          write(iulog,*) ' '
          write(iulog,*) 'jtest =', jtest
          write(iulog,*) 'After forward substitution, i, omega, gamma, denom, bu, xu:'
          do i = staggered_ihi, staggered_ilo, -1
             write(iulog,'(i4, 5e15.7)') i, omega_u(i,j), gamma_u(i,j), denom_u(i,j), bu(i,j), xu(i,j)
          enddo
       endif

    enddo  ! j

    !WHL - debug
    if (verbose_tridiag .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       write(iulog,*) 'xu tridiag solve, this_rank, j =', this_rank, j
       write(iulog,*) 'i, xu:'
       do i = staggered_ihi, staggered_ilo, -1
          write(iulog,*) i, xu(i,j)
       enddo
    endif   ! verbose_tridiag

    !-------------------------------------------------------------------------------------
    ! Solve a tridiagonal problem Av*xv = bv for each column i, from the left to the right of the global domain
    !-------------------------------------------------------------------------------------

    ! Initialize
    xv = 0.0d0
    gamma_v = 0.0d0

    ! Loop over locally owned columns
    do i = staggered_ilo, staggered_ihi

       ! Forward elimination
       ! Note: denom -> 1/denom in setup subroutine
       !       Adiag -> 1/Adiag in setup subroutine

       j = staggered_jlo
       gamma_v(i,j) = bv(i,j) * Adiag_v(i,j)   ! actually bv/Adiag_v
       do j = staggered_jlo+1, staggered_jhi
          gamma_v(i,j) = (bv(i,j) - Asubdiag_v(i,j)*gamma_v(i,j-1)) * denom_v(i,j)   ! actually (bv-Asubdiag_v*gamma_v)/denom_v
       enddo

       ! Back substitution
       j = staggered_jhi
       xv(i,j) = gamma_v(i,j)
       do j = staggered_jhi-1, 1, -1
          xv(i,j) = gamma_v(i,j) - omega_v(i,j)*xv(i,j+1)
       enddo

       if (verbose_tridiag .and. this_rank == rtest .and. i == itest) then
          write(iulog,*) ' '
          write(iulog,*) 'itest =', itest
          write(iulog,*) 'After forward substitution, j, omega, gamma, denom, bv, xv:'
          do j = staggered_jhi, staggered_jlo, -1
             write(iulog,'(i4, 5e15.7)') j, omega_v(i,j), gamma_v(i,j), denom_v(i,j), bv(i,j), xv(i,j)
          enddo
       endif

    enddo  ! i

    !WHL - debug
    if (verbose_tridiag .and. this_rank == rtest) then
       i = itest
       write(iulog,*) ' '
       write(iulog,*) 'xv tridiag solve, this_rank, i =', this_rank, i
       write(iulog,*) 'j, xv:'
       do j = staggered_jhi, staggered_jlo, -1
          write(iulog,*) j, xv(i,j)
       enddo
    endif   ! verbose_tridiag

  end subroutine tridiag_solver_local_2d

!****************************************************************************

  subroutine tridiag_solver_global_2d(ilocal,       jlocal,     &
                                      parallel,     tasks_rc,   &
                                      tridiag_solver_flag,      &
                                      itest, jtest, rtest,      &
                                      Adiag,                    &
                                      Asubdiag,     Asupdiag,   &
                                      omega,        denom,      &
                                      xuh,          xlh,        &
                                      bu,           xu,         &
                                      first_time,   gather_data)

    use glimmer_utils, only: tridiag
    use cism_parallel, only: gather_var_row, gather_var_col, &
         gather_all_var_row, gather_all_var_col, scatter_var_row, scatter_var_col

    integer, intent(in) :: &
         ilocal, jlocal            ! size of input/output arrays; number of locally owned vertices in each direction

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    character(len=3), intent(in) ::  &
         tridiag_solver_flag       ! either 'row' or 'col' depending on whether solving over rows or columns

    integer, intent(in) :: tasks_rc   ! number of tasks in this row or column

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

    logical, intent(in) :: &
         first_time                   ! if true, then gather_data is computed in full
                                      ! if false, then only some elements of gather_data need to be recomputed

    real(dp), dimension(8*tasks_rc,jlocal), intent(inout) :: &
         gather_data                  ! array for gathering data from all tasks on a row or column

    ! local variables

    logical :: main_task_rc           ! either main_task_row or main_task_col
    integer :: this_rank_rc           ! either this_rank_row or this_rank_col
    integer :: comm_rc                ! either comm_row or comm_col
    integer :: i, j, n, p, m
    integer :: i1, i2
    integer :: ibase, ibase_gather

    real(dp), dimension(ilocal,jlocal) :: &
         xr,                  & ! particular solution for Au tridiagonal solve
         gamma                  ! work array for Au tridiagonal solve

    real(dp), dimension(2,jlocal) :: &
         local_coeffs           ! coefficients for tridiagonal solution, scattered to local tasks

    real(dp), dimension(:,:), allocatable :: &
         outdata,             & ! data computed locally for global tridiagonal problem
         gather_data2,        & ! temporary array for gathering data from multiple tasks on a row or column
         global_coeffs          ! coefficients for tridiagonal solution, computed on main_task

    real(dp), dimension(:), allocatable :: &
         subdiag, diag, supdiag, rhs, coeffs  ! matrix entries for the triadiagonal problem solved on main_task

    ! Note: I have tested both options (gather_all = T and F) on Cheyenne, and found not much difference in computing cost.
!    logical, parameter :: gather_all = .false.  ! if false, then call gather routines, solve on main task, and scatter
    logical, parameter :: gather_all = .true.   ! if true, then call gather_all routines and solve on all tasks

    ! debug parameter; set tridiag_solver_scope = 0 to use this subroutine as a local tridiagonal solver
    integer, parameter :: &
!!         tridiag_solver_scope = 0   ! 0 = local, 1 = global
         tridiag_solver_scope = 1   ! 0 = local, 1 = global

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
       write(iulog,*) ' '
       write(iulog,*) 'In tridiag_solver_global_2d: itest, jtest, rtest =', itest, jtest, rtest
       if (tridiag_solver_flag == 'row') then
          write(iulog,*) 'Solving a tridiag problem for each matrix row'
       elseif (tridiag_solver_flag == 'col') then
          write(iulog,*) 'Solving a tridiag problem for each matrix column'
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
    !
    ! Setting up the input arrays in this form requires some overhead, but reduces striding for the Avv problem,
    ! and also makes it possible to solve the same subroutine twice (once for Auu and once for Avv), instead of
    ! calling a subroutine twice as long with much duplicated logic.

    if (tridiag_solver_flag == 'row') then
       comm_rc = parallel%comm_row
       main_task_rc = parallel%main_task_row
       this_rank_rc = parallel%this_rank_row
    elseif (tridiag_solver_flag == 'col') then
       comm_rc = parallel%comm_col
       main_task_rc = parallel%main_task_col
       this_rank_rc = parallel%this_rank_col
    endif

    if (verbose_tridiag .and. main_task) then
       write(iulog,*) 'tasks_rc, main_task_rc =', tasks_rc, main_task_rc
       call flush(6)
    endif

    !-------------------------------------------------------------------------------------
    ! Solve a tridiagonal problem Au*xu = bu for each row, from the bottom to the top of the global domain.
    ! (Alternatively, solve Av*xv = bv for each column, from the left to the right of the global domain.)
    !-------------------------------------------------------------------------------------

    ! Initialize the solution arrays
    xu(:,:) = 0.0d0
    xr(:,:) = 0.0d0

    if (first_time) then
       allocate(outdata(8,jlocal))  ! compute the global array in full
    else
       allocate(outdata(2,jlocal))  ! compute only array elements that depend on the rhs and need to be updated
    endif
    outdata = 0.0d0

    ! Loop over locally owned rows or columns, gathering data for the global problem.
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

       if (first_time) then   ! fill all 8 values
          outdata(1,j) = -1.0d0
          outdata(2,j) = xuh(1,j)
          outdata(3,j) = xlh(1,j)
          outdata(4,j) = -xr(1,j)
          outdata(5,j) = xuh(ilocal,j)
          outdata(6,j) = xlh(ilocal,j)
          outdata(7,j) = -1.0d0
          outdata(8,j) = -xr(ilocal,j)
       else   ! fill just the two rhs values that need to be updated
          outdata(1,j) = -xr(1,j)       ! element 4 above
          outdata(2,j) = -xr(ilocal,j)  ! element 8 above
       endif

       !WHL - debug
       if (verbose_tridiag .and. this_rank == rtest .and. j == jtest) then
          write(iulog,*) ' '
          if (tridiag_solver_flag == 'row') then
             write(iulog,*) 'jtest =', jtest
             write(iulog,*) 'After forward substitution, i, omega, gamma, denom, bu, xr, xuh, xlh:'
             do i = ilocal, 1, -1
                write(iulog,'(i4, 7e15.7)') i, omega(i,j), gamma(i,j), denom(i,j), bu(i,j), xr(i,j), xuh(i,j), xlh(i,j)
             enddo
          elseif (tridiag_solver_flag == 'col') then
             write(iulog,*) 'itest =', jtest   ! for columns, what we call jtest in the subroutine is really itest
             write(iulog,*) 'After forward substitution, j, omega, gamma, denom, bv, xr, xuh, xlh:'
             do i = jlocal, 1, -1
                write(iulog,'(i4, 7e15.7)') i, omega(i,j), gamma(i,j), denom(i,j), bu(i,j), xr(i,j), xuh(i,j), xlh(i,j)
             enddo
          endif
          write(iulog,*) ' '
          write(iulog,*) 'outdata(1:8)'
          write(iulog,'(8e12.3)') outdata(1:8,j)
       endif

    enddo  ! j

    call t_startf("pcg_tridiag_main_task")

    if (tasks_rc > 1 .and. tridiag_solver_scope == 1) then  ! scope = 1 implies a global tridiagonal solver

       if (gather_all) then   ! gather global data to all tasks in each row or column, and solve on all tasks

          ! Use the row- or column-based communicator to gather outdata on all tasks in this row or column.
          ! The global array is allocated in the subroutine.

          if (tridiag_solver_flag == 'row') then
             call t_startf("pcg_tridiag_gather_row")
             call gather_all_var_row(outdata, gather_data2, parallel)
             call t_stopf ("pcg_tridiag_gather_row")
          elseif (tridiag_solver_flag == 'col') then
             call t_startf("pcg_tridiag_gather_col")
             call gather_all_var_col(outdata, gather_data2, parallel)
             call t_stopf ("pcg_tridiag_gather_col")
          endif

       else   ! gather data to main_rank_rc where it will be solved

          ! Use the row- or column-based communicator to gather outdata on main_task_row or main_task_col.
          ! The global array is allocated in the subroutine.

          if (tridiag_solver_flag == 'row') then
             call t_startf("pcg_tridiag_gather_row")
             call gather_var_row(outdata, gather_data2, parallel)
             call t_stopf ("pcg_tridiag_gather_row")
          elseif (tridiag_solver_flag == 'col') then
             call t_startf("pcg_tridiag_gather_col")
             call gather_var_col(outdata, gather_data2, parallel)
             call t_stopf ("pcg_tridiag_gather_col")
          endif

       endif   ! gather_all

       ! Use gather_data to form the elements of a tridiagonal matrix, and solve the global problem.
       ! If gather_all = F, we solve on the main task of the row or column, and then scatter to other tasks.
       ! For gather_all = T, we solve the same problem on all tasks of a given row or column.

       if (main_task_rc .or. gather_all) then

          if (first_time) then

             ! copy in full from gather_data2
             gather_data = gather_data2

          else

             ! old value of gather_data was passed into subroutine;
             ! replace the values associated with the new rhs
             do n = 1, tasks_rc
                ibase = 8*(n-1)  ! base index of gather_data array
                ibase_gather = 2*(n-1)  ! base index of gather_data2 array
                do j = 1, jlocal
                   gather_data(ibase+4,j) = gather_data2(ibase_gather+1,j)
                   gather_data(ibase+8,j) = gather_data2(ibase_gather+2,j)
                enddo
             enddo

          endif

          allocate(global_coeffs(2*tasks_rc,jlocal))
          global_coeffs(:,:) = 0.0d0

          allocate(diag(2*tasks_rc-2))
          allocate(subdiag(2*tasks_rc-2))
          allocate(supdiag(2*tasks_rc-2))
          allocate(rhs(2*tasks_rc-2))
          allocate(coeffs(2*tasks_rc-2))

          do j = 1, jlocal
             do n = 1, tasks_rc

                ibase = 8*(n-1)  ! base index of gather_data array
                p = 2*n-2        ! matrix row

                if (n==1) then   ! westernmost or southernmost task
                   subdiag(p+1) = gather_data(ibase+5,j)
                   diag(p+1)    = gather_data(ibase+6,j)
                   supdiag(p+1) = gather_data(ibase+7,j)
                   rhs(p+1)     = gather_data(ibase+8,j)
                elseif (n==tasks_rc) then   ! easternmost or northernmost task
                   subdiag(p) = gather_data(ibase+1,j)
                   diag(p)    = gather_data(ibase+2,j)
                   supdiag(p) = gather_data(ibase+3,j)
                   rhs(p)     = gather_data(ibase+4,j)
                else
                   subdiag(p)   = gather_data(ibase+1,j)
                   diag(p)      = gather_data(ibase+2,j)
                   supdiag(p)   = gather_data(ibase+3,j)
                   rhs(p)       = gather_data(ibase+4,j)
                   subdiag(p+1) = gather_data(ibase+5,j)
                   diag(p+1)    = gather_data(ibase+6,j)
                   supdiag(p+1) = gather_data(ibase+7,j)
                   rhs(p+1)     = gather_data(ibase+8,j)
                endif

             enddo   ! n

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
                write(iulog,*) ' '
                write(iulog,*) 'Solved global tridiag problem'
                write(iulog,*) 'subdiag, diag, supdiag, rhs, coeffs:'
                do m = 1, 2*tasks_rc-2
                   write(iulog,'(5e12.3)') subdiag(m), diag(m), supdiag(m), rhs(m), coeffs(m)
                enddo
             endif

          enddo   ! j

          deallocate(diag, subdiag, supdiag, rhs, coeffs)

       endif   ! main_task_rc or gather_all

       if (gather_all) then

          ! Given the global array of coefficients, extract uh_coeff and lh_coeff for this task.
          i1 = 2*this_rank_rc + 1
          i2 = 2*this_rank_rc + 2
          local_coeffs(1:2,:) = global_coeffs(i1:i2,:)

       else   ! gather_all = F

          ! Scatter the coefficients to the local tasks.
          ! Each task receives 2 coefficients for each value of j:
          !  local_coeffs_u(1,:) = uh_coeff, and local_coeffs_u(2,:) = lh_coeff.
          ! Note: uh_coeff = 0 on the westernmost task and lh_coeff = 0 on the easternmost task.

          local_coeffs(:,:) = 0.d0

          if (tridiag_solver_flag == 'row') then
             call t_startf("pcg_tridiag_scatter_row")
             call scatter_var_row(local_coeffs, global_coeffs, parallel)
             call t_stopf ("pcg_tridiag_scatter_row")
          elseif (tridiag_solver_flag == 'col') then
             call t_startf("pcg_tridiag_scatter_col")
             call scatter_var_col(local_coeffs, global_coeffs, parallel)
             call t_stopf ("pcg_tridiag_scatter_col")
          endif

       endif   ! gather_all

       ! Use the coefficients to combine xr, xuh and xlh into the full solution xu.
       ! Note that xu has dimensions (nx-1,ny-1).

       do j = 1, jlocal
          do i = 1, ilocal
             xu(i,j) = xr(i,j) + local_coeffs(1,j)*xuh(i,j) + local_coeffs(2,j)*xlh(i,j)
          enddo
       enddo   ! j

    else    ! tasks_rc = 1 or tridiag_solver_scope = 0 (local solve)

       xu(:,:) = xr(:,:)

    endif   ! tasks_rc > 1

    call t_stopf("pcg_tridiag_main_task")

    ! clean up
    if (allocated(global_coeffs)) deallocate(global_coeffs)
    if (allocated(gather_data2)) deallocate(gather_data2)
    deallocate(outdata)

    !WHL - debug
    if (verbose_tridiag .and. this_rank == rtest) then
       j = jtest
       write(iulog,*) ' '
       if (tridiag_solver_flag == 'row') then
          write(iulog,*) 'xu tridiag solve, this_rank, j =', this_rank, j
          write(iulog,*) 'i, xu:'
       elseif (tridiag_solver_flag == 'col') then
          write(iulog,*) 'xv tridiag solve, this_rank, i =', this_rank, j
          write(iulog,*) 'j, xv:'
       endif
       do i = ilocal, 1, -1
          write(iulog,*) i, xu(i,j)
       enddo
    endif   ! verbose_tridiag

  end subroutine tridiag_solver_global_2d
 
!****************************************************************************

  subroutine matvec_multiply_structured_2d(nx,            ny,            &
                                           parallel,                     &
                                           indxA_2d,      active_vertex, &
                                           Auu,           Auv,           &
                                           Avu,           Avv,           &
                                           xu,            xv,            &
                                           yu,            yv)

    !---------------------------------------------------------------
    ! Compute the matrix-vector product $y = Ax$.
    !
    ! This subroutine is similar to subroutine matvec_multiply_structured_3d,
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
         nx, ny                 ! horizontal grid dimensions (for scalars)

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    integer, dimension(-1:1,-1:1), intent(in) :: &
         indxA_2d               ! maps relative (x,y) coordinates to an index between 1 and 9

    logical, dimension(nx-1,ny-1), intent(in) ::   &
         active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(nx-1,ny-1,9), intent(in) ::   &
         Auu, Auv, Avu, Avv     ! four components of assembled matrix
                                ! first two dimensions = (x,y) indices
                                ! 3rd dimension = 9 (node and its nearest neighbors in x and y directions)
                                !
                                !    Auu  | Auv
                                !    _____|____
                                !    Avu  | Avv
                                !           |

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
         xu, xv             ! current guess for solution


    real(dp), dimension(nx-1,ny-1), intent(out) ::  &
         yu, yv             ! y = Ax

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, m
    integer :: iA, jA

    integer :: &
       staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
       staggered_jlo, staggered_jhi

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    ! Initialize the result vector.

    yu(:,:) = 0.d0
    yv(:,:) = 0.d0

    ! Compute y = Ax
    ! Loop over locally owned vertices
    ! Note: For periodic BC, the southern and western rows of the global domain are halo cells.
    !          For all processors, staggered_ilo = staggered_jlo = staggered_lhalo + 1.
    !       For outflow BC, the southern and western rows of the global domain are locally owned.
    !          For processors owning these rows, staggered_ilo = staggered_jlo = staggered_lhalo.

    do jA = -1,1
       do iA = -1,1

          ! Get the third index of the Auu and Avv matrices
          m = indxA_2d(iA,jA)

          !dir$ ivdep
          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi

                if (active_vertex(i,j)) then

                   if ( (i+iA >= 1 .and. i+iA <= nx-1)       &
                                   .and.                     &
                        (j+jA >= 1 .and. j+jA <= ny-1) ) then

                      yu(i,j) = yu(i,j)   &
                              + Auu(i,j,m)*xu(i+iA,j+jA)  &
                              + Auv(i,j,m)*xv(i+iA,j+jA)

                      yv(i,j) = yv(i,j)   &
                              + Avu(i,j,m)*xu(i+iA,j+jA)  &
                              + Avv(i,j,m)*xv(i+iA,j+jA)

                   endif   ! i+iA, j+jA in bounds

                endif   ! active_vertex

             enddo   ! i
          enddo   ! j

       enddo   ! iA
    enddo   ! jA

  end subroutine matvec_multiply_structured_2d

!****************************************************************************

end module glissade_velo_higher_pcg_2d

!****************************************************************************
